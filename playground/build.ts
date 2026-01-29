import { mkdirSync } from 'fs';
import { join } from 'path';

const rootDir = new URL('..', import.meta.url).pathname;
const outDir = join(rootDir, 'playground', 'dist');

mkdirSync(outDir, { recursive: true });

// Generate lessons from tutorials, then build HTML + docs data
await import('./scripts/generate-lessons.ts');
await import('./build-lessons.ts');
await import('./generate-docs.ts');

const entrypoints = [
  join(rootDir, 'playground', 'sage-browser.ts'),
  join(rootDir, 'playground', 'main.js'),
  join(rootDir, 'playground', 'course.js'),
  join(rootDir, 'playground', 'docs.js'),
  join(rootDir, 'playground', 'landing.js'),
];

const result = await Bun.build({
  entrypoints,
  outdir: outDir,
  format: 'esm',
  target: 'browser',
  splitting: false,
  minify: false,
});

if (!result.success) {
  for (const message of result.logs) {
    console.error(message);
  }
  throw new Error('Playground bundle failed');
}

const outputFiles = result.outputs.map((out) => out.path).join(', ');
console.log(`Playground bundle written to ${outputFiles}`);
