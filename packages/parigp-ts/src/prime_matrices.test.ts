import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const cases = (await loadLiveNative(import.meta.url, "./prime_matrices.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import { FpM_inv, FpM_image, FpM_suppl } from './alglin1.js';
import { ZM_hnfmodprime } from './hnf_snf.js';
import { F2m_gauss_pivot, F2m_gauss } from './F2v.js';
export function pari_prime_matrix(
  op: bigint,
  mm: bigint,
  nn: bigint,
  p: bigint,
  flat: bigint[]
): string {
  const m = Number(mm),
    n = Number(nn),
    A = [
      [],
      ...Array.from({ length: n }, (_, j) => [
        0n,
        ...Array.from({ length: m }, (_, i) => flat[i * n + j]!),
      ]),
    ];
  let out: unknown;
  if (op === 0n)
    out =
      FpM_inv(A, p)
        ?.slice(1)
        .map((c) => c.slice(1)) ?? null;
  else if (op === 1n)
    out = FpM_image(A, p)
      .slice(1)
      .map((c) => c.slice(1));
  else if (op === 2n)
    out = FpM_suppl(A, p)
      .slice(1)
      .map((c) => c.slice(1));
  else if (op === 3n)
    out = ZM_hnfmodprime(
      A.slice(1).map((c) => c.slice(1)),
      p
    );
  else {
    const a = A.slice(1).map((c) =>
      c.slice(1).reduce((v, x, i) => v | ((x & 1n) << BigInt(i)), 0n)
    );
    if (op === 4n) {
      const [d, r] = F2m_gauss_pivot(a, m);
      out = [d?.map(BigInt) ?? 0n, BigInt(r)];
    } else {
      const x = F2m_gauss(
        a,
        Array.from({ length: m }, (_, i) => 1n << BigInt(i)),
        m
      );
      out = x?.map((c) => Array.from({ length: n }, (_, i) => (c >> BigInt(i)) & 1n)) ?? null;
    }
  }
  return JSON.stringify(out, (_, x) => (typeof x === 'bigint' ? String(x) : x));
}
const functions = { pari_prime_matrix };
const arg = (s: string): bigint | bigint[] =>
  s.startsWith('[')
    ? s.slice(1, -1).trim()
      ? s
          .slice(1, -1)
          .split(',')
          .map((v) => BigInt(v.trim()))
      : []
    : BigInt(s);
for (const [i, row] of cases.filter((row) => Object.hasOwn(functions, row.function)).entries())
  test('bundled prime matrix ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = (functions as any)[row.function](...row.args.map(arg));
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    expect({ result, error, errorType }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
  });
