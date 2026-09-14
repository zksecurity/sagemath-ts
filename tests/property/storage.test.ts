import { expect, test } from 'bun:test';
import { mkdtempSync, readFileSync, writeFileSync, rmSync, existsSync, statSync } from 'node:fs';
import { join } from 'node:path';
import { tmpdir } from 'node:os';
import { gzipSync, gunzipSync } from 'node:zlib';
import { readText, isCaseFile, caseArea } from './storage.js';
import { writeJSONArray, writeJSONArrayGzip } from './json.js';
import { normalizeFile } from './normalize-cases.js';
import { isCanonicalText } from './case-format.js';

test('compressed transcripts preserve exact existing serialization across chunks', async () => {
  const dir = mkdtempSync(join(tmpdir(), 'sage-gzip-'));
  try {
    for (const count of [0, 1, 64, 65, 129]) {
      const rows = Array.from({ length: count }, (_, seed) => ({
        seed,
        result: JSON.stringify(['🙂', '\ud800', 'é', '\n', seed]),
        error: null,
      }));
      const plain = join(dir, 'rows.json'),
        compressed = `${plain}.gz`;
      writeJSONArray(plain, rows);
      await writeJSONArrayGzip(compressed, rows);
      expect(gunzipSync(readFileSync(compressed))).toEqual(readFileSync(plain));
      expect(JSON.parse(readText(compressed))).toEqual(rows);
    }
    writeFileSync(join(dir, 'bad.gz'), 'not gzip');
    expect(() => readText(join(dir, 'bad.gz'))).toThrow();
  } finally {
    rmSync(dir, { recursive: true });
  }
});

test('case discovery retains plain and compressed area names', () => {
  for (const name of ['arith.cases.json', 'arith.cases.json.gz']) {
    expect(isCaseFile(name)).toBe(true);
    expect(caseArea(name)).toBe('arith');
  }
  for (const name of ['arith.json.gz', 'arith.cases.json.tmp', 'arith.cases.json.gz.tmp'])
    expect(isCaseFile(name)).toBe(false);
});

test('normalization preserves gzip and compresses large new suites without losing cases', async () => {
  const dir = mkdtempSync(join(tmpdir(), 'sage-case-gzip-'));
  try {
    const suite = { module: 'example', cases: [{ function: 'f', rows: [[1, 2]] }] };
    const path = join(dir, 'example.cases.json.gz');
    writeFileSync(path, gzipSync(JSON.stringify(suite)));
    expect((await normalizeFile(path)).cases).toBe(1);
    expect(isCanonicalText(readText(path))).toBe(true);
    expect(JSON.parse(readText(path))).toEqual(suite);

    const large = join(dir, 'large.cases.json');
    const rows = Array.from({ length: 1100 }, (_, seed) => [seed, Array(600).fill(0)]);
    writeFileSync(large, JSON.stringify({ module: 'large', cases: [{ function: 'f', rows }] }));
    expect(statSync(large).size).toBeGreaterThan(1_000_000);
    expect((await normalizeFile(large)).cases).toBe(rows.length);
    expect(existsSync(large)).toBe(false);
    expect(isCanonicalText(readText(`${large}.gz`))).toBe(true);
    expect(JSON.parse(readText(`${large}.gz`)).cases[0].rows).toEqual(rows);
  } finally {
    rmSync(dir, { recursive: true });
  }
});

test('property input definitions stay small and contain no compressed bulk corpora', () => {
  const root = join(import.meta.dir, '../..');
  let bytes = 0;
  for (const file of new Bun.Glob('tests/property/cases/*').scanSync({ cwd: root, onlyFiles: true })) {
    expect(file.endsWith('.cases.json'), file).toBe(true);
    const size = statSync(join(root, file)).size;
    expect(size, `${file}: use generators, not stored sweeps`).toBeLessThan(1_000_000);
    bytes += size;
  }
  expect(bytes).toBeLessThan(1_000_000);
  expect(statSync(join(root, 'tests/property/native-suites.json')).size).toBeLessThan(400_000);
  expect([...new Bun.Glob('packages/*/src/**/*.fixtures.json').scanSync({ cwd: root })]).toHaveLength(0);
  // These five small, hand-written resource/boundary controls are not bulk snapshots.
  const controls = new Set(['lll_resources.native.json', 'ideal_intersection.native.json',
    'trace_norm.bounds.native.json', 'valuation.hangs.native.json', 'valuation.resource.native.json']);
  for (const file of new Bun.Glob('packages/*/src/**/*.native.json{,.gz}').scanSync({ cwd: root, onlyFiles: true })) {
    expect(controls.has(file.split('/').at(-1)!), `${file}: use live native comparisons`).toBe(true);
    expect(statSync(join(root, file)).size).toBeLessThan(25_000);
  }
});

test('failed gzip serialization preserves the last complete transcript', async () => {
  const dir = mkdtempSync(join(tmpdir(), 'sage-atomic-gzip-'));
  try {
    const path = join(dir, 'rows.json.gz');
    await writeJSONArrayGzip(path, [{ result: 'previous native result' }]);
    const before = readFileSync(path);
    const circular: { self?: unknown } = {};
    circular.self = circular;
    await expect(writeJSONArrayGzip(path, [circular])).rejects.toThrow();
    expect(readFileSync(path)).toEqual(before);
    expect(existsSync(`${path}.tmp`)).toBe(false);
  } finally { rmSync(dir, { recursive: true }); }
});
