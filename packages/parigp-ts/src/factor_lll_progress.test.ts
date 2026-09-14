import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';
const fixtures = await loadLiveNative(import.meta.url, "./factor_lll_progress.native.json");
import { LLL_check_progress } from './QX_factor.js';
import { itor, shiftr } from './qfb.js';
export function pari_lll_progress(
  n0: bigint,
  final: bigint,
  m: bigint,
  n: bigint,
  p: bigint,
  shift: bigint,
  bound: bigint,
  flat: bigint[]
): string {
  const A = Array.from({ length: Number(n) }, (_, j) =>
    Array.from({ length: Number(m) }, (_, i) => flat[i * Number(n) + j]!)
  );
  return JSON.stringify(
    LLL_check_progress(shiftr(itor(bound, Number(p)), Number(shift)), Number(n0), A, !!final),
    (_, v) => (typeof v === 'bigint' ? String(v) : v)
  );
}

const invoke = pari_lll_progress;

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
  test('bundled PARI factor dependency ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = (invoke as any)(...row.args.map(arg));
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

test('factor progress preserves its source lattice', () => {
  const A = [
      [1n, 0n, 0n],
      [0n, 1n, 0n],
      [0n, 0n, 2n],
    ],
    saved = structuredClone(A);
  LLL_check_progress(itor(3n, 64), 2, A, true);
  expect(A).toEqual(saved);
});
