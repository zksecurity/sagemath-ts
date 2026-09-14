import { expect, spyOn, test } from 'bun:test';
import * as P from './ff.js';
import * as arithmetic from './_polynomial_division.js';
import { PariError } from './errors.js';

test('exported multiplication and inversion reduce signed operands correctly', () => {
  expect(P.Fp_inv(-1n, 17n)).toBe(16n);
  expect(P.Fp_mul(-1n, 1n, 17n)).toBe(16n);
  expect(P.Fp_addmul(0n, -1n, 1n, 17n)).toBe(16n);
  expect(P.Fp_red(-1n, -17n)).toBe(16n);
});

test('modii and remii preserve their different zero-divisor phases', () => {
  expect(P.Fp_red(0n, 0n)).toBe(0n);
  expect(P.Fp_mul(0n, 1n, 0n)).toBe(0n);
  expect(P.Fp_add(1n, -1n, 0n)).toBe(0n);
  expect(P.Fp_addmul(-1n, 1n, 1n, 0n)).toBe(0n);
  expect(() => P.Fp_sqr(0n, 0n)).toThrow(new PariError('impossible inverse in dvmdii: 0.'));
  expect(() => P.Fp_inv(1n, 0n)).toThrow(new PariError('impossible inverse in Fp_inv: Mod(1, 0).'));
});

test('word-modulus division returns zero before inversion but zero divisors use Fp_inv', () => {
  const inverse = spyOn(arithmetic, 'inverseCoefficient');
  try {
    expect(P.Fp_div(0n, 17n, 17n)).toBe(0n);
    expect(inverse).toHaveBeenCalledTimes(0);
    expect(() => P.Fp_div(0n, 0n, 17n)).toThrow(
      new PariError('impossible inverse in Fp_inv: Mod(17, 17).')
    );
    expect(inverse).toHaveBeenCalledTimes(1);
  } finally {
    inverse.mockRestore();
  }
});

test('large-modulus word division uses exact division and the divisor inverse', () => {
  const p = (1n << 64n) + 13n,
    inverse = spyOn(arithmetic, 'inverseCoefficient');
  try {
    expect(P.Fp_div(1n, 2n, p)).toBe((p + 1n) / 2n);
    expect(inverse.mock.calls[0]).toEqual([1n, 2n, true]);
    expect(P.Fp_div(1n, 2n, -p)).toBe((1n - p) / 2n);
    expect(() => P.Fp_div(0n, 2n, p + 1n)).toThrow(
      new PariError('impossible inverse in Fl_inv: Mod(0, 2).')
    );
    expect(P.Fp_div(0n, 1n, 0n)).toBe(0n);
  } finally {
    inverse.mockRestore();
  }
});

test('centering and doubling retain native single-comparison arithmetic', () => {
  expect(P.Fp_center(-10n, 17n)).toBe(-27n);
  expect(P.Fp_double(35n, 17n)).toBe(53n);
  expect(P.Fp_double(1n, 0n)).toBe(2n);
  expect(P.Fp_halve(-3n, 2n)).toBe(0n);
});
