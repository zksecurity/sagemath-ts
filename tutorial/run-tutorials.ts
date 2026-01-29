/**
 * Tutorial Runner
 *
 * This script discovers and runs all tutorial files in the tutorial directory.
 * Each tutorial file should be a standalone TypeScript file that demonstrates
 * cryptographic concepts using the sagemath-ts library.
 *
 * Usage:
 *   bun tutorial/run-tutorials.ts [part] [chapter]
 *
 * Examples:
 *   bun tutorial/run-tutorials.ts                    # Run all tutorials
 *   bun tutorial/run-tutorials.ts part1-foundations  # Run Part 1 only
 *   bun tutorial/run-tutorials.ts part1-foundations 01-integers-and-divisibility  # Run specific chapter
 */

import { readdirSync, statSync, existsSync } from 'fs';
import { join, relative, dirname } from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);
const TUTORIAL_DIR = __dirname;

interface TutorialResult {
  path: string;
  success: boolean;
  output: string;
  error?: string;
  duration: number;
}

async function runTutorialFile(filePath: string): Promise<TutorialResult> {
  const startTime = Date.now();

  return new Promise((resolve) => {
    const proc = Bun.spawn([process.execPath, filePath], {
      cwd: join(TUTORIAL_DIR, '..'),
      stdout: 'pipe',
      stderr: 'pipe',
    });

    Promise.all([
      new Response(proc.stdout).text(),
      new Response(proc.stderr).text(),
      proc.exited,
    ])
      .then(([stdout, stderr, code]) => {
        const duration = Date.now() - startTime;
        resolve({
          path: relative(TUTORIAL_DIR, filePath),
          success: code === 0,
          output: stdout,
          error: stderr || undefined,
          duration,
        });
      })
      .catch((err) => {
        const duration = Date.now() - startTime;
        resolve({
          path: relative(TUTORIAL_DIR, filePath),
          success: false,
          output: '',
          error: err instanceof Error ? err.message : String(err),
          duration,
        });
      });
  });
}

function discoverTutorials(dir: string, partFilter?: string, chapterFilter?: string): string[] {
  const tutorials: string[] = [];

  const parts = readdirSync(dir).filter(f => {
    const fullPath = join(dir, f);
    return statSync(fullPath).isDirectory() && f.startsWith('part');
  }).sort();

  for (const part of parts) {
    if (partFilter && !part.includes(partFilter)) continue;

    const partPath = join(dir, part);
    const chapters = readdirSync(partPath).filter(f => {
      const fullPath = join(partPath, f);
      return statSync(fullPath).isDirectory();
    }).sort();

    for (const chapter of chapters) {
      if (chapterFilter && !chapter.includes(chapterFilter)) continue;

      const chapterPath = join(partPath, chapter);
      const files = readdirSync(chapterPath).filter(f =>
        f.endsWith('.ts') && !f.endsWith('.test.ts') && !f.endsWith('.d.ts')
      ).sort();

      for (const file of files) {
        tutorials.push(join(chapterPath, file));
      }
    }
  }

  return tutorials;
}

async function main() {
  const [partFilter, chapterFilter] = process.argv.slice(2);

  console.log('🎓 Cryptography Tutorial Runner\n');
  console.log('=' .repeat(60));

  if (partFilter) {
    console.log(`Filter: part=${partFilter}${chapterFilter ? `, chapter=${chapterFilter}` : ''}`);
  }
  console.log('');

  const tutorials = discoverTutorials(TUTORIAL_DIR, partFilter, chapterFilter);

  if (tutorials.length === 0) {
    console.log('No tutorials found.');
    process.exit(0);
  }

  console.log(`Found ${tutorials.length} tutorial file(s)\n`);

  const results: TutorialResult[] = [];
  let currentPart = '';
  let currentChapter = '';

  for (const tutorial of tutorials) {
    const parts = relative(TUTORIAL_DIR, tutorial).split('/');
    const part = parts[0];
    const chapter = parts[1];

    if (part !== currentPart) {
      currentPart = part;
      console.log(`\n📚 ${part}`);
      console.log('-'.repeat(40));
    }

    if (chapter !== currentChapter) {
      currentChapter = chapter;
      console.log(`\n  📖 ${chapter}`);
    }

    process.stdout.write(`    ▸ ${parts[2]}... `);

    const result = await runTutorialFile(tutorial);
    results.push(result);

    if (result.success) {
      console.log(`✅ (${result.duration}ms)`);
    } else {
      console.log(`❌ (${result.duration}ms)`);
      if (result.error) {
        console.log(`      Error: ${result.error.split('\n')[0]}`);
      }
    }
  }

  // Summary
  console.log('\n' + '='.repeat(60));
  console.log('📊 Summary\n');

  const passed = results.filter(r => r.success).length;
  const failed = results.filter(r => !r.success).length;
  const totalDuration = results.reduce((sum, r) => sum + r.duration, 0);

  console.log(`  Total:    ${results.length} tutorials`);
  console.log(`  Passed:   ${passed} ✅`);
  console.log(`  Failed:   ${failed} ❌`);
  console.log(`  Duration: ${(totalDuration / 1000).toFixed(2)}s`);

  if (failed > 0) {
    console.log('\n❌ Failed tutorials:');
    for (const result of results.filter(r => !r.success)) {
      console.log(`  - ${result.path}`);
      if (result.error) {
        const errorLines = result.error.split('\n').slice(0, 3);
        for (const line of errorLines) {
          console.log(`      ${line}`);
        }
      }
    }
    process.exit(1);
  }

  console.log('\n✅ All tutorials passed!');
  process.exit(0);
}

main().catch(console.error);
