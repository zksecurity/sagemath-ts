import { expect, test } from 'bun:test';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

const cli = fileURLToPath(new URL('./cli.ts', import.meta.url));
function run(input: string, args: string[] = []) {
  return spawnSync(process.execPath, [cli, '--input', '-', ...args], { input, encoding: 'utf8' });
}

test('CLI handles local synthetic JSON records and safe reports', () => {
  const result = run(JSON.stringify([{ type: 'rsa', n: '101', e: '3', label: 'DO_NOT_ECHO' }]));
  expect(result.status).toBe(0);
  expect(result.stderr).toBe('');
  expect(JSON.parse(result.stdout).records[0].checks[0].status).toBe('matched');
  expect(result.stdout).not.toContain('DO_NOT_ECHO');
  expect(result.stdout).not.toContain('101');
});

test('CLI does not echo invalid options, paths or malformed source content', () => {
  for (const [input, args] of [
    ['SECRET_BAD_JSON', []],
    ['{}', []],
    ['[]', ['--unknown', 'SECRET_OPTION']],
    ['[]', ['--max-keys', '-1']],
    ['[]', ['--max-keys', 'NaN']],
    ['[]', ['--input', '/SECRET_PATH']],
  ] as Array<[string, string[]]>) {
    const result = run(input, args);
    expect(result.status).toBe(2);
    expect(result.stdout).toBe('');
    expect(result.stderr.trim()).toBe('weak-keys: invalid input or options; use --help');
  }
});

test('CLI rejects oversized input before JSON parsing', () => {
  const result = run(' '.repeat(1048577));
  expect(result.status).toBe(2);
  expect(result.stdout).toBe('');
});

test('CLI documents its bounded scope without running a scan', () => {
  const result = spawnSync(process.execPath, [cli, '--help'], { encoding: 'utf8' });
  expect(result.status).toBe(0);
  expect(result.stdout).toContain('No network');
});
