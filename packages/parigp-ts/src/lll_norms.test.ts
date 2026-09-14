import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./lll_norms.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import { ZM_lll_norms } from './lll.js';
export function pari_lll_norms(
  op: bigint,
  flag: bigint,
  mode: bigint,
  m: bigint,
  n: bigint,
  p: bigint,
  shift: bigint,
  flat: bigint[]
): string {
  const A = Array.from({ length: Number(n) }, (_, j) =>
    Array.from({ length: Number(m) }, (_, i) => flat[i * Number(n) + j]!)
  );
  const r = ZM_lll_norms(A, op === 75n ? 0.75 : op === 999n ? 0.999 : 0.99, Number(flag));
  const out = (x: any): unknown =>
    typeof x === 'bigint'
      ? String(x)
      : Array.isArray(x)
        ? x.map(out)
        : x && typeof x === 'object'
          ? [x.s, String(x.e), String(x.m), x.p]
          : x;
  return JSON.stringify(out(r));
}

export function pari_lll_norms_resource(...args: Parameters<typeof pari_lll_norms>): string {
  try {
    pari_lll_norms(...args);
  } catch (e) {
    if (
      e instanceof RangeError &&
      e.message === 'fplll_dpe: Gram-Schmidt coefficient requires an unrepresentable shift'
    )
      return 'resource_failure';
    throw e;
  }
  return 'unexpected_success';
}

const functions = { pari_lll_norms, pari_lll_norms_resource };

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
  test('native adaptive LLL norms ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = (functions[row.function as keyof typeof functions] as any)(...row.args.map(arg));
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    expect({ result, error, errorType }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
  }, 5000);

test('LLL norms preserve caller input', () => {
  const B = [
      [4n, 1n, 0n],
      [3n, 2n, 0n],
      [1n, 1n, 1n],
    ],
    saved = structuredClone(B);
  ZM_lll_norms(B);
  expect(B).toEqual(saved);
});
