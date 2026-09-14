import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./hnf_permutation.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import { ZM_hnfperm, ZM_hnf_knapsack, hnfperm } from './hnf_snf.js';
export function pari_hnfperm(op: bigint, m: bigint, n: bigint, flat: bigint[]): string {
  const A = Array.from({ length: Number(n) }, (_, j) =>
    Array.from({ length: Number(m) }, (_, i) => flat[i * Number(n) + j]!)
  );
  const out =
    op === 4n
      ? ZM_hnf_knapsack(A)
      : op === 5n
        ? hnfperm(A)
        : ZM_hnfperm(A, !!(op & 1n), !!(op & 2n));
  return JSON.stringify(out, (_, v) =>
    typeof v === 'bigint' || typeof v === 'number' ? String(v) : v
  );
}

const invoke = pari_hnfperm;

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

test('permuted HNF and knapsack preserve input columns', () => {
  const A = [
      [2n, 1n, 0n],
      [3n, -1n, 1n],
    ],
    saved = structuredClone(A);
  ZM_hnfperm(A, true, true);
  expect(A).toEqual(saved);
  ZM_hnf_knapsack(A);
  expect(A).toEqual(saved);
});
