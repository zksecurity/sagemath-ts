const rootUrl = new URL('.', import.meta.url);
const rootPath = rootUrl.pathname;
const bundleUrl = new URL('./dist/sage-browser.js', rootUrl);
const mainBundleUrl = new URL('./dist/main.js', rootUrl);
const courseBundleUrl = new URL('./dist/course.js', rootUrl);
const docsBundleUrl = new URL('./dist/docs.js', rootUrl);
const landingBundleUrl = new URL('./dist/landing.js', rootUrl);

async function ensureBundle() {
  const bundleFile = Bun.file(bundleUrl);
  const mainFile = Bun.file(mainBundleUrl);
  const courseFile = Bun.file(courseBundleUrl);
  const docsFile = Bun.file(docsBundleUrl);
  const landingFile = Bun.file(landingBundleUrl);
  if (
    (await bundleFile.exists()) &&
    (await mainFile.exists()) &&
    (await courseFile.exists()) &&
    (await docsFile.exists()) &&
    (await landingFile.exists())
  ) {
    return;
  }

  console.log('Building browser bundle...');
  const build = Bun.spawnSync([process.execPath, new URL('./build.ts', rootUrl).pathname], {
    cwd: rootPath,
    stdout: 'inherit',
    stderr: 'inherit',
  });
  if (!build.success) {
    throw new Error('Failed to build playground bundle');
  }
}

await ensureBundle();

function sanitizePath(pathname: string): string | null {
  const cleaned = pathname.replace(/\0/g, '');
  if (cleaned.includes('..')) return null;
  return cleaned;
}

function contentTypeFor(pathname: string): string | undefined {
  if (pathname.endsWith('.css')) return 'text/css';
  if (pathname.endsWith('.js')) return 'text/javascript';
  if (pathname.endsWith('.html')) return 'text/html';
  if (pathname.endsWith('.json')) return 'application/json';
  return undefined;
}

const preferredPort = Number(Bun.env.PORT ?? 5173);

function createServer(port: number) {
  return Bun.serve({
    port,
    async fetch(request) {
      const url = new URL(request.url);
      const safePath = sanitizePath(url.pathname);
      if (!safePath) {
        return new Response('Bad path', { status: 400 });
      }

      const relative = safePath === '/' ? '/index.html' : safePath;
      const fileUrl = new URL(`.${relative}`, rootUrl);
      const file = Bun.file(fileUrl);

      if (!(await file.exists())) {
        return new Response('Not found', { status: 404 });
      }

      const type = contentTypeFor(relative);
      return new Response(file, {
        headers: type ? { 'Content-Type': type } : undefined,
      });
    },
  });
}

let server: ReturnType<typeof Bun.serve> | null = null;
let port = preferredPort;
for (let attempt = 0; attempt < 5; attempt += 1) {
  try {
    server = createServer(port);
    break;
  } catch (err) {
    console.error(`Port ${port} failed:`, err instanceof Error ? err.message : err);
    port += 1;
  }
}

if (!server) {
  const base = 30000 + Math.floor(Math.random() * 1000);
  for (let attempt = 0; attempt < 10; attempt += 1) {
    try {
      server = createServer(base + attempt);
      break;
    } catch (err) {
      console.error(`Port ${base + attempt} failed:`, err instanceof Error ? err.message : err);
      continue;
    }
  }
  if (!server) {
    throw new Error('Failed to start playground server on ports 5173-5177');
  }
}

console.log(`Playground running at http://localhost:${server.port}`);
