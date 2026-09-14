import { expect, spyOn, test } from 'bun:test';
import { Fp_pow } from './ff.js';
import * as arithmetic from './_polynomial_division.js';
import * as powering from './bb_group.js';
import { PariError } from './errors.js';

// Permanent native counterparts: ff_pari_fp_power in the shared property area.
test('zero exponent uses divisibility before any modular operation', () => {
  expect(Fp_pow(0n, 0n, 0n)).toBe(0n);
  expect(Fp_pow(17n, 0n, 17n)).toBe(0n);
  expect(Fp_pow(1n, 0n, 0n)).toBe(1n);
  expect(Fp_pow(0n, 2n, 0n)).toBe(0n);
  expect(() => Fp_pow(1n, 2n, 0n)).toThrow(new PariError('impossible inverse in dvmdii: 0.'));
});

test('negative exponent width selects the native inverse and payload', () => {
  const inverse = spyOn(arithmetic, 'inverseCoefficient');
  try {
    expect(() => Fp_pow(0n, -1n, -17n)).toThrow(
      new PariError('impossible inverse in Fl_inv: Mod(0, 17).')
    );
    expect(inverse.mock.calls[0]).toEqual([0n, 17n, true]);
    expect(() => Fp_pow(0n, -(1n << 64n), -17n)).toThrow(
      new PariError('impossible inverse in Fp_inv: Mod(17, -17).')
    );
    expect(inverse.mock.calls[1]).toEqual([0n, -17n, false]);
    expect(Fp_pow(-1n, -1n, 17n)).toBe(16n);
  } finally {
    inverse.mockRestore();
  }
});

test('large exponents use PARI windows and fused base-two powering', () => {
  const window = spyOn(powering, 'gen_pow_i'),
    fold = spyOn(powering, 'gen_pow_fold');
  try {
    expect(Fp_pow(3n, (1n << 64n) + 1n, -17n)).toBe(20n);
    expect(window).toHaveBeenCalledTimes(1);
    expect(Fp_pow(-1n, (1n << 64n) + 1n, -17n)).toBe(-50n);
    expect(Fp_pow(2n, (1n << 64n) + 1n, 17n)).toBe(2n);
    expect(fold).toHaveBeenCalledTimes(1);
  } finally {
    window.mockRestore();
    fold.mockRestore();
  }
});

test('negative Montgomery doubling reports its native stack-exhaustion boundary', () => {
  const k = (1n << 64n) + 1n;
  expect(() => Fp_pow(2n, k, -((1n << 64n) - 59n))).toThrow(
    new RangeError('Fp_pow: negative modulus makes Montgomery doubling nonterminating')
  );
  expect(Fp_pow(2n, k, -((1n << 65n) - 59n))).toBe(71239208179867822193n);
});
