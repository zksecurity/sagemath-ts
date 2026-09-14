import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./lzz_pX_minpoly.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import { MinPolySeq } from './lzz_pX1.js';
export function ntl_word_minpoly(p: bigint, m: bigint, a: bigint[]): string {
  return JSON.stringify(MinPolySeq(a, Number(m), p), (_, v) =>
    typeof v === 'bigint' ? String(v) : v
  );
}

const arg = (s: string): bigint | bigint[] =>
  s.startsWith('[')
    ? s.slice(1, -1).trim()
      ? s
          .slice(1, -1)
          .split(',')
          .map((v) => BigInt(v.trim()))
      : []
    : BigInt(s);
for (const [i, row] of fixtures.entries())
  test('native NTL word minimum sequence ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = (ntl_word_minpoly as Function)(...row.args.map(arg));
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    expect({ result, error, errorType }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
  }, 10000);
test('NTL sequence reconstruction preserves inputs', () => {
  const a = [0n, 1n, 1n, 2n, 3n, 5n],
    h = MinPolySeq(a, 3, 101n);
  expect(h).toEqual([100n, 100n, 1n]);
  h[0] = 99n;
  expect(a).toEqual([0n, 1n, 1n, 2n, 3n, 5n]);
});
test('NTL sequence reconstruction ignores unused suffix and returns fresh arrays', () => {
  const a = [0n, 1n, 1n, 2n, 3n, 5n],
    h = MinPolySeq(a, 3, 101n),
    g = MinPolySeq([...a, 999n, -1n], 3, 101n);
  expect(g).toEqual(h);
  expect(g).not.toBe(h);
  g[0] = 99n;
  expect(h).toEqual([100n, 100n, 1n]);
});
