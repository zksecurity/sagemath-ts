import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const cases = (await loadLiveNative(import.meta.url, "./rational_polynomials.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import { RgX_rem, RgXQ_mul, RgXQ_norm, RgXQ_trace } from './RgX.js';
import { QX_mul, QX_ZX_rem, ZX_rem } from './ZX.js';
import { QXQ_mul, QXQ_norm } from './polarit3.js';
import { gnorm, gtrace } from './alglin2.js';
import { type RationalPair, type RationalPolynomialData } from './_rational_polynomial.js';
export function pari_rational_polynomial(
  op: bigint,
  aa: bigint[],
  ad: bigint,
  bb: bigint[],
  bd: bigint,
  tt: bigint[],
  td: bigint
): string {
  const a: RationalPolynomialData = [aa, ad],
    b: RationalPolynomialData = [bb, bd],
    T: RationalPolynomialData = [tt, td];
  let z: unknown;
  if (op === 0n) z = RgX_rem(a, b);
  else if (op === 1n) z = QX_mul(a, b);
  else if (op === 2n) z = QX_ZX_rem(a, b[0]);
  else if (op === 3n) z = QXQ_mul(a, b, T[0]);
  else if (op === 4n) z = RgXQ_mul(a, b, T);
  else if (op === 10n) z = [ZX_rem(a[0], b[0]), 1n];
  else {
    const q: RationalPair =
      op === 5n
        ? QXQ_norm(a, T[0])
        : op === 6n
          ? RgXQ_norm(a, T)
          : op === 7n
            ? RgXQ_trace(a, T)
            : op === 8n
              ? gnorm({ value: RgX_rem(a, T), modulus: T })
              : gtrace({ value: RgX_rem(a, T), modulus: T });
    z = q[1] === 1n ? String(q[0]) : q[0] + '/' + q[1];
  }
  return JSON.stringify(z, (_, x) => (typeof x === 'bigint' ? String(x) : x));
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
for (const [i, row] of cases.entries())
  test('original trace/norm rational_polynomials ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = (pari_rational_polynomial as any)(...row.args.map(arg));
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
