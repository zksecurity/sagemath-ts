import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./factor_bounds.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import { rationalPair } from './_rational_polynomial.js';
import { Mignotte_bound, Beauzamy_bound, factor_bound, root_bound } from './QX_factor.js';
import { vecbinomial } from './bibli2.js';
import { ZX_Z_eval } from './ZX.js';
import { ceil_safe } from './gen3.js';
import { powruhalf } from './trans1.js';
import { cmpir, cmpri } from './kernel/none/level1.js';
import { itor, shiftr, type MpReal, type MpComplex } from './qfb.js';
export function pari_factor_bounds(
  op: bigint,
  p: bigint,
  shift: bigint,
  n: bigint,
  a: bigint | bigint[]
): string {
  let result: unknown;
  if (op === 0n)
    result = ceil_safe(
      n === 0n
        ? shiftr(itor(a as bigint, Number(p)), Number(shift))
        : (a as bigint | [bigint, bigint])
    );
  else if (op === 1n) result = powruhalf(shiftr(itor(a as bigint, Number(p)), Number(shift)), n);
  else if (op === 6n) result = vecbinomial(Number(n));
  else if (op === 10n) {
    const [a0, b0] = a as bigint[];
    result = rationalPair(a0!, b0!);
  } else if (op === 7n || op === 8n) {
    const [x, y] = a as bigint[],
      r = shiftr(itor(y!, Number(p)), Number(shift));
    result = op === 7n ? cmpir(x!, r) : cmpri(r, x!);
  } else if (op === 9n) {
    const [y, ...coefficients] = a as bigint[];
    result = ZX_Z_eval(coefficients, y!);
  } else
    result =
      op === 2n
        ? Mignotte_bound(a as bigint[])
        : op === 3n
          ? Beauzamy_bound(a as bigint[])
          : op === 4n
            ? factor_bound(a as bigint[])
            : root_bound(a as bigint[]);
  const encode = (x: unknown): unknown =>
    typeof x === 'bigint' || typeof x === 'number'
      ? String(x)
      : Array.isArray(x)
        ? x.map(encode)
        : x && typeof x === 'object'
          ? 're' in x
            ? [encode((x as MpComplex).re), encode((x as MpComplex).im)]
            : [(x as MpReal).s, String((x as MpReal).e), String((x as MpReal).m), (x as MpReal).p]
          : x;
  return JSON.stringify(encode(result));
}

const invoke = pari_factor_bounds;

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

test('integer-real comparisons return native integer zero, not negative zero', () => {
  const z = itor(0n, 64),
    one = itor(1n, 64);
  for (const c of [cmpir(0n, z), cmpri(z, 0n), cmpir(1n, one), cmpri(one, 1n)]) {
    expect(c).toBe(0);
    expect(Object.is(c, -0)).toBe(false);
  }
});
test('polynomial factor bounds and sparse evaluation preserve coefficients', () => {
  const A = [-2n, 0n, 0n, 1n],
    saved = A.slice();
  Mignotte_bound(A);
  Beauzamy_bound(A);
  factor_bound(A);
  root_bound(A);
  ZX_Z_eval(A, -17n);
  expect(A).toEqual(saved);
});
