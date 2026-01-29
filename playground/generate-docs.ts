import { readdirSync, readFileSync, statSync, writeFileSync } from 'fs';
import { join, resolve } from 'path';
import * as ts from 'typescript';

const rootDir = new URL('..', import.meta.url).pathname;
const sourceRoot = join(rootDir, 'packages', 'sagemath-ts', 'src');
const outputPath = join(rootDir, 'playground', 'docs-data.json');

type DocItem = { name: string; kind: string; doc?: string; from: string };

type DocSection = { title: string; items: DocItem[] };

function walk(dir: string): string[] {
  const entries = readdirSync(dir);
  const files: string[] = [];
  for (const entry of entries) {
    const full = join(dir, entry);
    const stats = statSync(full);
    if (stats.isDirectory()) {
      files.push(...walk(full));
    } else if (stats.isFile() && full.endsWith('.ts') && !full.endsWith('.test.ts')) {
      files.push(full);
    }
  }
  return files;
}

function getDoc(node: ts.Node): string | undefined {
  const docs = ts.getJSDocCommentsAndTags(node);
  const comment = docs
    .map((doc) => (ts.isJSDoc(doc) && doc.comment ? String(doc.comment) : undefined))
    .filter(Boolean)
    .join('\n');
  return comment || undefined;
}

function exportKind(node: ts.Node): string {
  if (ts.isFunctionDeclaration(node)) return 'function';
  if (ts.isClassDeclaration(node)) return 'class';
  if (ts.isInterfaceDeclaration(node)) return 'interface';
  if (ts.isTypeAliasDeclaration(node)) return 'type';
  if (ts.isEnumDeclaration(node)) return 'enum';
  if (ts.isVariableStatement(node)) return 'const';
  return 'export';
}

function hasExportModifier(node: ts.Node): boolean {
  return Boolean(node.modifiers?.some((mod) => mod.kind === ts.SyntaxKind.ExportKeyword));
}

function collectExports(filePath: string): DocItem[] {
  const sourceText = readFileSync(filePath, 'utf8');
  const source = ts.createSourceFile(filePath, sourceText, ts.ScriptTarget.ES2022, true);
  const relative = filePath.replace(sourceRoot + '/', '');
  const items: DocItem[] = [];

  source.forEachChild((node) => {
    if (!hasExportModifier(node)) return;

    if (ts.isFunctionDeclaration(node) && node.name) {
      items.push({
        name: node.name.text,
        kind: exportKind(node),
        doc: getDoc(node),
        from: relative,
      });
      return;
    }

    if ((ts.isClassDeclaration(node) || ts.isInterfaceDeclaration(node) || ts.isEnumDeclaration(node)) && node.name) {
      items.push({
        name: node.name.text,
        kind: exportKind(node),
        doc: getDoc(node),
        from: relative,
      });
      return;
    }

    if (ts.isTypeAliasDeclaration(node)) {
      items.push({
        name: node.name.text,
        kind: exportKind(node),
        doc: getDoc(node),
        from: relative,
      });
      return;
    }

    if (ts.isVariableStatement(node)) {
      for (const decl of node.declarationList.declarations) {
        if (ts.isIdentifier(decl.name)) {
          items.push({
            name: decl.name.text,
            kind: exportKind(node),
            doc: getDoc(node),
            from: relative,
          });
        }
      }
    }
  });

  return items;
}

const files = walk(sourceRoot);
const sectionsMap = new Map<string, DocSection>();

for (const filePath of files) {
  const items = collectExports(filePath);
  if (items.length === 0) continue;
  const title = filePath
    .replace(sourceRoot + '/', '')
    .replace(/\.(t|j)s$/, '');
  sectionsMap.set(title, { title, items });
}

const sections = Array.from(sectionsMap.values()).sort((a, b) => a.title.localeCompare(b.title));

const payload = {
  generatedAt: new Date().toISOString(),
  source: 'packages/sagemath-ts/src',
  sections,
};

writeFileSync(outputPath, JSON.stringify(payload, null, 2));
console.log(`Docs data written to ${outputPath}`);
