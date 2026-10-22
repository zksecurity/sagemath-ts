import { expect, test } from 'bun:test';
import { mpfr_init2, mpfr_pow, mpfr_set_d, mpfr_set_z } from './index.js';
import type { mpfr_t } from './types.js';

function pow(x: bigint | number, y: number): [mpfr_t, number] {
  const a = mpfr_init2(53);
  const b = mpfr_init2(53);
  const z = mpfr_init2(53);
  if (typeof x === 'bigint') mpfr_set_z(a, x);
  else mpfr_set_d(a, x);
  mpfr_set_d(b, y);
  const ternary = mpfr_pow(z, a, b);
  return [z, ternary];
}

/** The exact value `mantissa * 2^(exponent - precision)` as `[numerator, log2 denominator]`. */
function value(z: mpfr_t): [bigint, number] {
  let m = BigInt(z.sign) * z.mantissa;
  let e = z.exponent - z.precision;
  while (m !== 0n && (m & 1n) === 0n) {
    m >>= 1n;
    e++;
  }
  return e >= 0 ? [m << BigInt(e), 0] : [m, -e];
}

// Expected values are Sage's ZZ(x)**RR(y), i.e. MPFR 4.2.1 at 53 bits.
test('mpfr_pow rounds inexact powers correctly', () => {
  expect(value(pow(3n, 0.5)[0])).toEqual([3900231685776981n, 51]);
  expect(value(pow(1n << 100n, 0.375)[0])).toEqual([6369051672525773n, 15]);
  expect(value(pow(2n, -0.5)[0])).toEqual([6369051672525773n, 53]);
  expect(value(pow(7n, -3)[0])).toEqual([6722574370885405n, 61]);
  expect(pow(3n, 0.5)[1]).toBe(-1);
});

test('mpfr_pow returns exact powers exactly', () => {
  expect(pow(9n, 0.5)).toEqual([expect.objectContaining({ kind: 'finite' }), 0]);
  expect(value(pow(9n, 0.5)[0])).toEqual([3n, 0]);
  expect(value(pow(6561n, 0.25)[0])).toEqual([9n, 0]);
  expect(value(pow(1n << 100n, 0.25)[0])).toEqual([1n << 25n, 0]);
  expect(value(pow(-2n, 3)[0])).toEqual([-8n, 0]);
});

test('mpfr_pow special values follow C99', () => {
  expect(pow(Number.NaN, 0)[0].kind).toBe('finite');
  expect(value(pow(1n, Number.NaN)[0])).toEqual([1n, 0]);
  expect(pow(-2n, 0.5)[0].kind).toBe('nan');
  expect(pow(0, -1)[0]).toMatchObject({ kind: 'inf', sign: 1 });
  expect(pow(-0, -3)[0]).toMatchObject({ kind: 'inf', sign: -1 });
  expect(pow(0.5, Number.POSITIVE_INFINITY)[0]).toMatchObject({ kind: 'zero', sign: 1 });
  expect(pow(1.5, Number.POSITIVE_INFINITY)[0]).toMatchObject({ kind: 'inf', sign: 1 });
  expect(value(pow(-1n, Number.NEGATIVE_INFINITY)[0])).toEqual([1n, 0]);
  expect(pow(Number.NEGATIVE_INFINITY, 3)[0]).toMatchObject({ kind: 'inf', sign: -1 });
});
