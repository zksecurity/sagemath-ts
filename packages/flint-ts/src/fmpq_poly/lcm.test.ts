import { expect, test } from 'bun:test';
import { _fmpq_poly_lcm } from './lcm.js';
import { _fmpz_poly_lcm } from '../fmpz_poly/lcm.js';
import { fmpq_poly_get_numerator, fmpq_poly_get_denominator } from './get_numerator_denominator.js';

test('native LCM kernels retain integer content and monic rational output', () => {
  expect(_fmpz_poly_lcm([-2n, 2n], [3n, 3n])).toEqual([-6n, 0n, 6n]);
  expect(_fmpq_poly_lcm([-2n, 2n], [3n, 3n])).toEqual([[-1n, 0n, 1n], 1n]);
  expect(_fmpq_poly_lcm([1n, 2n], [1n, 2n])).toEqual([[1n, 2n], 2n]);
  expect(_fmpz_poly_lcm([], [1n])).toEqual([]);
});

test('native rational polynomial accessors expose canonical integer storage', () => {
  const a = [-1n, 0n, 3n];
  expect(fmpq_poly_get_numerator(a, 2n)).toEqual(a);
  expect(fmpq_poly_get_numerator(a, 2n)).not.toBe(a);
  expect(fmpq_poly_get_denominator(a, 2n)).toBe(2n);
});
