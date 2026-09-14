import { expect, test } from 'bun:test';
import { _fmpz_poly_pow } from './pow.js';
import { _fmpq_poly_pow } from '../fmpq_poly/pow.js';

test('FLINT power valuation, exact binomial coefficients and raw rational denominator', () => {
  // Original fmpz_poly_pow/fmpq_poly_pow results; broader shared cases compare native calls.
  expect(_fmpz_poly_pow([0n, 2n, -1n], 3n)).toEqual([0n, 0n, 0n, 8n, -12n, 6n, -1n]);
  expect(_fmpz_poly_pow([], 0n)).toEqual([1n]);
  expect(_fmpz_poly_pow([-1n], (1n << 64n) - 1n)).toEqual([-1n]);
  expect(_fmpq_poly_pow([2n], 2n, 3n)).toEqual([[8n], 8n]);
  expect(_fmpq_poly_pow([], 6n, 3n)).toEqual([[], 1n]);
});
