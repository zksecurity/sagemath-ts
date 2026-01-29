/**
 * Tutorial Test Suite
 *
 * This test file ensures all tutorial files are runnable.
 * Run with: bun test tutorial/tutorial.test.ts
 */

import { describe, it } from 'bun:test';
import { readdirSync, statSync } from 'fs';
import { join, relative, dirname } from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);
const TUTORIAL_DIR = __dirname;
const ROOT_DIR = join(TUTORIAL_DIR, '..');

interface TutorialFile {
  part: string;
  chapter: string;
  file: string;
  path: string;
}

function discoverTutorials(): TutorialFile[] {
  const tutorials: TutorialFile[] = [];

  const parts = readdirSync(TUTORIAL_DIR).filter(f => {
    const fullPath = join(TUTORIAL_DIR, f);
    return statSync(fullPath).isDirectory() && f.startsWith('part');
  }).sort();

  for (const part of parts) {
    const partPath = join(TUTORIAL_DIR, part);
    const chapters = readdirSync(partPath).filter(f => {
      const fullPath = join(partPath, f);
      return statSync(fullPath).isDirectory();
    }).sort();

    for (const chapter of chapters) {
      const chapterPath = join(partPath, chapter);
      const files = readdirSync(chapterPath).filter(f =>
        f.endsWith('.ts') && !f.endsWith('.test.ts') && !f.endsWith('.d.ts')
      ).sort();

      for (const file of files) {
        tutorials.push({
          part,
          chapter,
          file,
          path: join(chapterPath, file),
        });
      }
    }
  }

  return tutorials;
}

const tutorials = discoverTutorials();

// Group tutorials by part and chapter for better test organization
const tutorialsByPart = tutorials.reduce((acc, t) => {
  if (!acc[t.part]) acc[t.part] = {};
  if (!acc[t.part][t.chapter]) acc[t.part][t.chapter] = [];
  acc[t.part][t.chapter].push(t);
  return acc;
}, {} as Record<string, Record<string, TutorialFile[]>>);

describe('Cryptography Tutorials', () => {
  for (const [part, chapters] of Object.entries(tutorialsByPart)) {
    describe(part, () => {
      for (const [chapter, files] of Object.entries(chapters)) {
        describe(chapter, () => {
          for (const tutorial of files) {
            it(`should run ${tutorial.file}`, async () => {
              const relativePath = relative(ROOT_DIR, tutorial.path);
              const proc = Bun.spawn([process.execPath, relativePath], {
                cwd: ROOT_DIR,
                stdout: 'pipe',
                stderr: 'pipe',
              });

              const timeout = setTimeout(() => proc.kill(), 30_000);

              const [stdout, stderr, exitCode] = await Promise.all([
                new Response(proc.stdout).text(),
                new Response(proc.stderr).text(),
                proc.exited,
              ]);

              clearTimeout(timeout);

              if (exitCode !== 0) {
                throw new Error(
                  `Tutorial failed to run:\n` +
                  `  File: ${relativePath}\n` +
                  `  Error: ${stderr || stdout || `exit code ${exitCode}`}`
                );
              }
            });
          }
        });
      }
    });
  }
});
