import { expect, spyOn, test } from 'bun:test';
import * as matrices from './FpV.js';
import * as native from './FpX.js';
import * as composition from './_polynomial_composition.js';
import { FpXQ_autpow, FpXQ_autpowers } from './galconj.js';
import { FpX_mul, FpX_Fp_mul } from './ffinit.js';
import { Flxq_powers, Flx_Flxq_eval } from './Flx.js';
import { brent_kung_optpow } from './RgX.js';

test('PARI signed polynomial coefficients reduce canonically', () => {
  expect(FpX_mul([1n, -1n], [1n, 1n], 3n)).toEqual([1n, 0n, 2n]);
  expect(FpX_Fp_mul([-1n, 1n], 1n, 3n)).toEqual([2n, 1n]);
  expect(native.FpXQ_powers([-1n], 2, [1n, 0n, 1n], 3n)[1]).toEqual([-1n]);
  expect(native.FpXQ_powers([-1n], 3, [1n, 0n, 1n], 3n)[1]).toEqual([2n]);
});

test('PARI modular evaluation preserves original truncation and safe word contracts', () => {
  const p = (1n << 64n) + 13n;
  expect(native.FpX_FpXQ_eval([0n, 1n], [1n, 2n, 3n], [1n, 0n, 1n], p)).toEqual([1n, 2n]);
  expect(native.FpX_FpXQ_eval([1n], [], [1n], p)).toEqual([]);
  expect(() => Flx_Flxq_eval([0n, 1n], [1n, 2n, 3n], [1n, 0n, 1n], 5n)).toThrow(
    'power polynomial exceeds modulus degree'
  );
  expect(brent_kung_optpow(10, 2, 1)).toBe(5);
  expect(Flxq_powers([1n, 1n], 3, [1n, 0n, 1n], 5n)).toEqual([[1n], [1n, 1n], [0n, 2n], [3n, 2n]]);
  expect(
    native.FpXQ_auttrace(
      [
        [0n, 1n],
        [1n, 1n],
      ],
      3n,
      [1n, 0n, 1n],
      5n
    )
  ).toEqual([
    [0n, 1n],
    [3n, 3n],
  ]);
});

test('PARI automorphism callers reuse tables and perform logarithmic matrix compositions', () => {
  const p = (1n << 64n) + 13n,
    x = [0n, 1n],
    T = [1n, 1n, 1n];
  const product = spyOn(matrices, 'FpM_mul');
  try {
    expect(FpXQ_autpow(x, 2 ** 20 + 1, T, p)).toEqual(x);
    expect(product).toHaveBeenCalledTimes(21);
  } finally {
    product.mockRestore();
  }
  // The prepared quotient context enters the shared table kernel directly.
  const powers = spyOn(composition, 'quotientPowers');
  try {
    expect(
      FpXQ_autpowers(x, 9, T, p)
        .slice(1)
        .every((v) => v.join(',') === '0,1')
    ).toBe(true);
    expect(powers).toHaveBeenCalledTimes(1);
  } finally {
    powers.mockRestore();
  }
});
