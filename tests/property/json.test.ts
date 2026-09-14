import { expect, test } from 'bun:test';
import { mkdtempSync, readFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { writeJSONArray } from './json.js';

test('chunked transcripts preserve Unicode, nested result JSON and chunk boundaries', () => {
  const dir = mkdtempSync(join(tmpdir(), 'sage-transcript-'));
  try {
    for (const count of [0, 1, 63, 64, 65, 129]) {
      const rows = Array.from({ length: count }, (_, seed) => ({
        seed,
        result: JSON.stringify(['🙂', '\ud800', '\n', 'é', seed]),
        error: null,
      }));
      const path = join(dir, 'output.json');
      writeJSONArray(path, rows);
      const text = readFileSync(path, 'utf8');
      expect(JSON.parse(text)).toEqual(rows);
      expect([...text].every((c) => c.charCodeAt(0) < 128)).toBe(true);
    }
  } finally {
    rmSync(dir, { recursive: true });
  }
});
