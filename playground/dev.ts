/**
 * Development server with hot reload.
 * Watches tutorial/ and playground/ for changes and rebuilds automatically.
 */

import { watch } from 'fs';
import { join } from 'path';

const rootDir = new URL('..', import.meta.url).pathname;
const playgroundDir = join(rootDir, 'playground');
const tutorialDir = join(rootDir, 'tutorial');
const packagesDir = join(rootDir, 'packages');

let buildInProgress = false;
let pendingBuild = false;

async function build() {
  if (buildInProgress) {
    pendingBuild = true;
    return;
  }

  buildInProgress = true;
  console.log('\n🔨 Building...');
  const start = performance.now();

  try {
    const proc = Bun.spawn(['bun', 'playground/build.ts'], {
      cwd: rootDir,
      stdout: 'pipe',
      stderr: 'pipe',
    });

    const [stdout, stderr] = await Promise.all([
      new Response(proc.stdout).text(),
      new Response(proc.stderr).text(),
    ]);

    const exitCode = await proc.exited;

    if (exitCode !== 0) {
      console.error('❌ Build failed:');
      if (stderr) console.error(stderr);
      if (stdout) console.log(stdout);
    } else {
      const elapsed = (performance.now() - start).toFixed(0);
      console.log(`✅ Built in ${elapsed}ms`);
    }
  } catch (err) {
    console.error('❌ Build error:', err);
  } finally {
    buildInProgress = false;
    if (pendingBuild) {
      pendingBuild = false;
      await build();
    }
  }
}

// Debounce file changes
let debounceTimer: ReturnType<typeof setTimeout> | null = null;

function scheduleRebuild(path: string) {
  if (debounceTimer) {
    clearTimeout(debounceTimer);
  }
  console.log(`📝 Changed: ${path}`);
  debounceTimer = setTimeout(() => {
    debounceTimer = null;
    build();
  }, 100);
}

// Watch directories
function watchDir(dir: string, name: string) {
  try {
    watch(dir, { recursive: true }, (event, filename) => {
      if (!filename) return;
      // Ignore generated files and dist
      if (filename.includes('lessons.generated') || filename.includes('/dist/')) return;
      // Only watch relevant files
      if (filename.endsWith('.ts') || filename.endsWith('.js') || filename.endsWith('.css')) {
        scheduleRebuild(`${name}/${filename}`);
      }
    });
    console.log(`👀 Watching ${name}/`);
  } catch (err) {
    console.error(`Failed to watch ${dir}:`, err);
  }
}

// Initial build
await build();

// Start watching
watchDir(tutorialDir, 'tutorial');
watchDir(playgroundDir, 'playground');
watchDir(packagesDir, 'packages');

// Start server
const serverModule = await import('./server.ts');

console.log('\n🔥 Hot reload enabled. Edit files to rebuild automatically.\n');
