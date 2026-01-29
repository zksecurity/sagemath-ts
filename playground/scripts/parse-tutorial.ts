/**
 * Parses a tutorial TypeScript file and extracts structured lesson content.
 *
 * Tutorial files have this structure:
 * 1. Top-level JSDoc comment → intro paragraphs
 * 2. Import statement → extracted for transformation
 * 3. Example sections delimited by `// ===...Example N:...===`
 * 4. Within sections: JSDoc comments → paragraphs, code → code blocks
 */

export type LessonBlock =
  | { type: 'paragraph'; text: string }
  | { type: 'code'; code: string; prefix: string };

export type ParsedLesson = {
  id: string;
  title: string;
  blocks: LessonBlock[];
};

const EXAMPLE_DELIMITER = /^\/\/ =+\n\/\/ Example \d+:/m;
const EXAMPLE_HEADER = /^\/\/ =+\n\/\/ (Example \d+: .+)\n\/\/ =+$/m;

/**
 * Extract the top-level JSDoc comment from a file.
 */
function extractTopJsDoc(content: string): { doc: string; rest: string } | null {
  const match = content.match(/^\/\*\*\n([\s\S]*?)\n \*\/\n/);
  if (!match) return null;

  const rawDoc = match[1]!;
  // Remove leading " * " from each line
  const doc = rawDoc
    .split('\n')
    .map((line) => line.replace(/^ \* ?/, ''))
    .join('\n')
    .trim();

  return { doc, rest: content.slice(match[0].length) };
}

/**
 * Extract imports from the content.
 * Handles both single-line and multi-line import statements.
 */
function extractImports(content: string): { imports: string[]; rest: string } {
  const imports: string[] = [];
  let rest = content;

  // Match import statements from 'sagemath-ts' (single-line or multi-line)
  const importPattern = /import\s*\{([^}]+)\}\s*from\s*['"]sagemath-ts['"];?\n?/gs;
  let match;

  while ((match = importPattern.exec(content)) !== null) {
    // Extract names from the import
    const names = match[1]!
      .split(',')
      .map((s) => s.trim())
      .filter((s) => s.length > 0);
    imports.push(...names);

    // Remove the import statement from rest
    rest = rest.replace(match[0], '');
  }

  // Also remove any other import statements (non-sagemath imports)
  rest = rest.replace(/import\s+.*from\s+['"][^'"]+['"];?\n?/g, '');

  return { imports, rest: rest.trim() };
}

/**
 * Split content into sections by example delimiters.
 */
function splitIntoSections(content: string): string[] {
  // Split on the example header pattern
  const parts = content.split(/\n(?=\/\/ =+\n\/\/ Example \d+:)/);
  return parts.filter((p) => p.trim());
}

/**
 * Extract JSDoc comments from a section.
 */
function extractJsDocComments(code: string): { paragraphs: string[]; code: string } {
  const paragraphs: string[] = [];
  let remaining = code;

  // Find all JSDoc comments
  const jsDocPattern = /\/\*\*\n([\s\S]*?)\n \*\//g;
  let match;

  while ((match = jsDocPattern.exec(code)) !== null) {
    const rawDoc = match[1]!;
    const doc = rawDoc
      .split('\n')
      .map((line) => line.replace(/^ \* ?/, ''))
      .join('\n')
      .trim();
    if (doc) {
      paragraphs.push(doc);
    }
    // Remove this JSDoc from remaining code
    remaining = remaining.replace(match[0], '');
  }

  return { paragraphs, code: remaining };
}

// Bun transpiler for stripping TypeScript types
const transpiler = new Bun.Transpiler({ loader: 'ts' });

/**
 * Strip TypeScript type annotations using Bun's transpiler.
 */
function stripTypeScript(code: string): string {
  try {
    return transpiler.transformSync(code);
  } catch {
    // If transpilation fails, return original code
    return code;
  }
}

/**
 * Transform code for browser execution:
 * - console.log → print
 * - Remove example header console.log
 * - Strip TypeScript type annotations
 */
function transformCode(code: string): string {
  let transformed = code
    // Remove the example header console.log (e.g., console.log("=== Example 1: ... ===\n");)
    .replace(/console\.log\s*\(\s*["']=== Example \d+:.*===.*["']\s*\);\s*\n?/g, '')
    // Transform remaining console.log to print
    .replace(/console\.log\(/g, 'print(');

  // Strip TypeScript types
  transformed = stripTypeScript(transformed);

  return transformed
    // Clean up multiple blank lines
    .replace(/\n{3,}/g, '\n\n')
    .trim();
}

/**
 * Clean up code by removing empty lines at start/end and normalizing whitespace.
 */
function cleanCode(code: string): string {
  return code
    .split('\n')
    .map((line) => line.trimEnd())
    .join('\n')
    .trim();
}

/**
 * Parse a section into paragraphs and code.
 */
function parseSection(section: string): { title: string; paragraphs: string[]; code: string } {
  // Extract title from header
  const headerMatch = section.match(/^\/\/ =+\n\/\/ (Example \d+: .+)\n\/\/ =+\n/);
  const title = headerMatch ? headerMatch[1]! : 'Example';

  // Remove the header
  let content = headerMatch ? section.slice(headerMatch[0].length) : section;

  // Extract JSDoc comments as paragraphs
  const { paragraphs, code } = extractJsDocComments(content);

  return {
    title,
    paragraphs,
    code: cleanCode(transformCode(code)),
  };
}

/**
 * Convert intro JSDoc to paragraphs, splitting on ## headers.
 */
function introToParagraphs(doc: string): string[] {
  // Split by ## headers, keeping the content
  const sections = doc.split(/\n(?=## )/);
  const paragraphs: string[] = [];

  for (const section of sections) {
    // Remove # Title header (first section might have it)
    let content = section.replace(/^# .+\n+/, '').trim();
    if (!content) continue;

    // For ## sections, include the header
    if (content.startsWith('## ')) {
      // Convert ## Header to bold text
      content = content.replace(/^## (.+)\n/, '**$1**: ');
    }

    // Clean up markdown artifacts
    content = content
      .replace(/\n\s*-\s+/g, '\n• ') // Convert - bullets to •
      .replace(/\n\s*\d+\.\s+/g, (m) => m.replace(/\d+\./, (n) => n)) // Keep numbered lists
      .trim();

    if (content) {
      paragraphs.push(content);
    }
  }

  return paragraphs;
}

/**
 * Parse a tutorial file into a lesson structure.
 */
export function parseTutorial(content: string, id: string): ParsedLesson {
  const blocks: LessonBlock[] = [];

  // 1. Extract top JSDoc as intro
  const topDoc = extractTopJsDoc(content);
  if (!topDoc) {
    throw new Error(`Tutorial must start with a JSDoc comment`);
  }

  // Convert intro to paragraphs
  const introParagraphs = introToParagraphs(topDoc.doc);
  for (const text of introParagraphs) {
    blocks.push({ type: 'paragraph', text });
  }

  // 2. Extract imports
  const { imports, rest } = extractImports(topDoc.rest);
  const sageImport = imports.length > 0 ? `const { ${imports.join(', ')} } = Sage;\n\n` : '';

  // 3. Split into example sections
  const sections = splitIntoSections(rest);

  // 4. Parse each section and build cumulative prefix
  let cumulativeCode = '';

  for (const section of sections) {
    const { title, paragraphs, code } = parseSection(section);

    // Skip exercise sections for now (they have different semantics)
    if (title.toLowerCase().includes('exercise')) {
      continue;
    }

    // Add paragraphs from this section
    for (const text of paragraphs) {
      blocks.push({ type: 'paragraph', text });
    }

    // Add code block with cumulative prefix
    if (code) {
      blocks.push({
        type: 'code',
        code,
        prefix: sageImport + cumulativeCode,
      });

      // Add this code to cumulative prefix for next blocks
      cumulativeCode += code + '\n\n';
    }
  }

  // Extract title from first paragraph or use id
  const titleMatch = topDoc.doc.match(/^# (.+)$/m);
  const title = titleMatch ? titleMatch[1]! : id;

  return { id, title, blocks };
}

/**
 * Parse a tutorial file from disk.
 */
export async function parseTutorialFile(filePath: string): Promise<ParsedLesson> {
  const content = await Bun.file(filePath).text();

  // Generate id from filename
  const id = filePath
    .split('/')
    .pop()!
    .replace(/\.ts$/, '');

  return parseTutorial(content, id);
}
