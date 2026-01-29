/**
 * Generate playground lessons from tutorial files.
 *
 * This script:
 * 1. Auto-discovers all .ts files in tutorial/ directory
 * 2. Parses each tutorial using parse-tutorial.ts
 * 3. Generates a lessons.generated.ts file for the playground
 */

import { readdirSync, statSync } from 'fs';
import { join, relative } from 'path';
import { parseTutorialFile, type ParsedLesson } from './parse-tutorial.ts';

const rootDir = new URL('../..', import.meta.url).pathname;
const tutorialDir = join(rootDir, 'tutorial');
const outputPath = join(rootDir, 'playground', 'lessons', 'lessons.generated.ts');

type TutorialFile = {
  partDir: string; // e.g., 'part1-foundations'
  chapterDir: string; // e.g., '01-integers-and-divisibility'
  file: string; // e.g., 'gcd.ts'
  fullPath: string;
  order: number; // derived from chapter number
};

/**
 * Discover all tutorial files in the tutorial directory.
 */
function discoverTutorials(): TutorialFile[] {
  const tutorials: TutorialFile[] = [];

  const parts = readdirSync(tutorialDir)
    .filter((entry) => entry.startsWith('part') && statSync(join(tutorialDir, entry)).isDirectory())
    .sort();

  for (const partDir of parts) {
    const partPath = join(tutorialDir, partDir);
    const chapters = readdirSync(partPath)
      .filter((entry) => {
        const fullPath = join(partPath, entry);
        return statSync(fullPath).isDirectory() && /^\d+-/.test(entry);
      })
      .sort();

    for (const chapterDir of chapters) {
      const chapterPath = join(partPath, chapterDir);
      const files = readdirSync(chapterPath)
        .filter((f) => f.endsWith('.ts') && !f.startsWith('index'))
        .sort();

      // Extract chapter number for ordering (e.g., "01-integers" -> 1)
      const chapterMatch = chapterDir.match(/^(\d+)-/);
      const chapterNum = chapterMatch ? parseInt(chapterMatch[1]!, 10) : 0;

      for (const file of files) {
        tutorials.push({
          partDir,
          chapterDir,
          file,
          fullPath: join(chapterPath, file),
          order: chapterNum * 100 + files.indexOf(file),
        });
      }
    }
  }

  return tutorials.sort((a, b) => a.order - b.order);
}

/**
 * Generate a lesson ID from tutorial path.
 */
function generateLessonId(tutorial: TutorialFile): string {
  // Use filename without extension as ID
  return tutorial.file.replace(/\.ts$/, '');
}

/**
 * Build navigation mapping (chapter path -> first lesson ID).
 */
function buildNavMapping(tutorials: TutorialFile[]): Record<string, string> {
  const mapping: Record<string, string> = {};
  for (const tutorial of tutorials) {
    const chapterPath = `${tutorial.partDir}/${tutorial.chapterDir}`;
    if (!mapping[chapterPath]) {
      mapping[chapterPath] = generateLessonId(tutorial);
    }
  }
  return mapping;
}

/**
 * Parse all discovered tutorials and return lessons.
 */
async function parseAllTutorials(tutorials: TutorialFile[]): Promise<ParsedLesson[]> {
  const lessons: ParsedLesson[] = [];

  for (const tutorial of tutorials) {
    try {
      const lesson = await parseTutorialFile(tutorial.fullPath);
      lesson.id = generateLessonId(tutorial);
      lessons.push(lesson);

      const relativePath = relative(tutorialDir, tutorial.fullPath);
      console.log(`✓ Parsed: ${relativePath} → ${lesson.id}`);
    } catch (err) {
      const relativePath = relative(tutorialDir, tutorial.fullPath);
      console.error(`✗ Failed to parse ${relativePath}:`, err instanceof Error ? err.message : err);
    }
  }

  return lessons;
}

/**
 * Generate TypeScript source for lessons.
 */
function generateLessonsSource(lessons: ParsedLesson[], navMapping: Record<string, string>): string {
  // Serialize lessons to TypeScript
  const lessonsJson = JSON.stringify(lessons, null, 2);

  return `/**
 * AUTO-GENERATED FILE - DO NOT EDIT
 * Generated from tutorial files by playground/scripts/generate-lessons.ts
 * Run: bun playground/scripts/generate-lessons.ts
 */

export type LessonBlock =
  | { type: 'paragraph'; text: string }
  | { type: 'code'; code: string; prefix: string };

export type Lesson = {
  id: string;
  title: string;
  blocks: LessonBlock[];
};

/**
 * Maps chapter paths to lesson IDs for navigation.
 * Chapters with mappings will be enabled in the sidebar.
 */
export const lessonNavMapping: Record<string, string> = ${JSON.stringify(navMapping, null, 2)};

/**
 * Generated lessons from tutorial files.
 */
export const lessons: Lesson[] = ${lessonsJson};
`;
}

// Main execution
const tutorials = discoverTutorials();
console.log(`Found ${tutorials.length} tutorial files\n`);

const lessons = await parseAllTutorials(tutorials);
const navMapping = buildNavMapping(tutorials);
const source = generateLessonsSource(lessons, navMapping);

await Bun.write(outputPath, source);
console.log(`\n✓ Generated ${lessons.length} lessons → ${relative(rootDir, outputPath)}`);
