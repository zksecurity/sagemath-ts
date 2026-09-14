import { expect, test } from 'bun:test';
import { mkdtempSync, readFileSync, rmSync, statSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { tmpdir } from 'node:os';

const cli = join(import.meta.dir, '../scripts/quiet.ts');
async function quiet(args: string[]) {
  const child = Bun.spawn([process.execPath, cli, ...args], { stdout: 'pipe', stderr: 'pipe' });
  const [stdout, stderr, exit] = await Promise.all([
    new Response(child.stdout).text(), new Response(child.stderr).text(), child.exited,
  ]);
  return { stdout, stderr, exit };
}
function logPath(output: string) {
  const path = /Log: (.+)\n/.exec(output)?.[1];
  if (!path) throw Error('Wrapper did not report a log');
  return path;
}
function removeLog(path: string) { rmSync(join(path, '..'), { recursive: true, force: true }); }

test('multi-megabyte single lines and repeated successes stay bounded; early failures and replay survive', async () => {
  const result = await quiet(['--', process.execPath, '-e', `
    console.log('Property seed: 12345; runs per generator: 25');
    console.error('error: early failure');
    console.error('Expected: 3'); console.error('Received: 4');
    console.log('9'.repeat(4 * 1024 * 1024));
    for (let i = 0; i < 20000; i++) console.log('(pass) irrelevant success ' + i);
    console.log('20000 pass'); console.log('1 fail'); process.exitCode = 7;
  `]);
  const path = logPath(result.stdout);
  try {
    expect(result.exit).toBe(7);
    expect(result.stdout.length).toBeLessThanOrEqual(6000);
    expect(result.stdout).toContain('early failure');
    expect(result.stdout).toContain('Expected: 3');
    expect(result.stdout).toContain('Received: 4');
    expect(result.stdout).toContain('Property seed: 12345');
    expect(result.stdout).toContain('20000 pass');
    expect(result.stdout).toContain('1 fail');
    expect(statSync(path).size).toBeGreaterThan(4 * 1024 * 1024);
    expect(readFileSync(path, 'utf8')).toContain('9'.repeat(4 * 1024 * 1024));
    const inspected = await quiet(['show', path, '--match', 'early failure']);
    expect(inspected.stdout).toContain('early failure');
    expect(inspected.stdout.length).toBeLessThanOrEqual(6000);
  } finally { removeLog(path); }
}, 15000);

test('argument boundaries are preserved without shell interpretation, and raw UTF-8/control bytes are retained', async () => {
  const argument = 'literal `echo BAD` $(echo BAD) ; spaced';
  const result = await quiet(['--', process.execPath, '-e', `
    process.stdout.write('\u001b[31m🙂 café\u001b[0m\\n');
    console.log(process.argv[1]);
  `, argument]);
  const path = logPath(result.stdout);
  try {
    expect(result.exit).toBe(0);
    expect(result.stdout).toContain(argument);
    expect(result.stdout).toContain('🙂 café');
    expect(result.stdout).not.toContain('\u001b');
    expect(readFileSync(path, 'utf8')).toBe('\u001b[31m🙂 café\u001b[0m\n' + argument + '\n');
  } finally { removeLog(path); }
});

test('log inspection finds text beyond a huge line prefix, handles CRLF, and enforces record/character limits', async () => {
  const directory = mkdtempSync(join(tmpdir(), 'quiet-inspect-test-'));
  const path = join(directory, 'output.log');
  try {
    writeFileSync(path, 'start\r\n' + 'a'.repeat(65526) + 'NEEDLE' + 'b'.repeat(200000) + '\r\nend\n');
    const result = await quiet(['show', path, '--match', 'NEEDLE', '--max-chars', '1000']);
    expect(result.exit).toBe(0);
    expect(result.stdout.length).toBeLessThanOrEqual(1000);
    expect(result.stdout).toContain('L2:');
    expect(result.stdout).toContain('chars]');
    expect(result.stdout).toContain('NEEDLE');
    expect(result.stdout).toContain('Scanned 3 records; selected 1');
    const page = await quiet(['show', path, '--from', '3', '--lines', '1']);
    expect(page.stdout).toContain('L3: end');
    expect(page.stdout).not.toContain('L2:');
  } finally { rmSync(directory, { recursive: true, force: true }); }
});

test('tight output budget and a missing executable still report failure', async () => {
  const result = await quiet(['--max-chars', '1000', '--', 'sagemath-command-that-does-not-exist']);
  const path = logPath(result.stdout);
  try {
    expect(result.exit).toBe(127);
    expect(result.stdout.length).toBeLessThanOrEqual(1000);
    expect(result.stdout).toContain('Exit: 127');
  } finally { removeLog(path); }
  expect((await quiet(['--max-chars', '0', '--', 'anything'])).exit).toBe(2);
});

test('timeout terminates descendants and returns 124 rather than swallowing failure', async () => {
  const directory = mkdtempSync(join(tmpdir(), 'quiet-child-test-'));
  const marker = join(directory, 'child.txt');
  const result = await quiet(['--timeout', '1', '--', process.execPath, '-e', `
    const {spawn} = require('node:child_process');
    const child = spawn(process.execPath, ['-e', 'setInterval(() => {}, 1000)'], {stdio:'inherit'});
    require('node:fs').writeFileSync(${JSON.stringify(marker)}, String(child.pid));
    setInterval(() => {}, 1000);
  `]);
  const path = logPath(result.stdout);
  try {
    expect(result.exit).toBe(124);
    expect(result.stdout).toContain('(timeout)');
    if (process.platform !== 'win32') {
      const pid = Number(readFileSync(marker, 'utf8'));
      expect(() => process.kill(pid, 0)).toThrow();
    }
  } finally {
    removeLog(path); rmSync(directory, { recursive: true, force: true });
  }
}, 10000);

test('cancellation also kills descendants that ignore SIGTERM and detach their stdio', async () => {
  if (process.platform === 'win32') return;
  const directory = mkdtempSync(join(tmpdir(), 'quiet-ignored-child-'));
  const marker = join(directory, 'child.txt');
  const grandchild = `process.on('SIGTERM',()=>{}); require('node:fs').writeFileSync(${JSON.stringify(marker)}, String(process.pid)); setInterval(()=>{},1000);`;
  const result = await quiet(['--timeout', '1', '--', process.execPath, '-e', `
    require('node:child_process').spawn(process.execPath, ['-e', ${JSON.stringify(grandchild)}], {stdio:'ignore'});
    setInterval(()=>{},1000);
  `]);
  const path = logPath(result.stdout);
  try {
    expect(result.exit).toBe(124);
    const pid = Number(readFileSync(marker, 'utf8'));
    expect(() => process.kill(pid, 0)).toThrow();
  } finally {
    removeLog(path); rmSync(directory, { recursive: true, force: true });
  }
}, 10000);
