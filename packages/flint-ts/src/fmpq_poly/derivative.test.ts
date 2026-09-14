import { expect, test } from 'bun:test';
import { _fmpz_poly_derivative } from '../fmpz_poly/derivative.js';
import { _fmpq_poly_derivative } from './derivative.js';

test('native derivatives normalize outputs and cancel the rational denominator', () => {
  const a = [1n, 0n, 3n, 0n];
  expect(_fmpz_poly_derivative(a)).toEqual([0n, 6n]);
  expect(_fmpq_poly_derivative(a, 6n)).toEqual([[0n, 1n], 1n]);
  expect(_fmpq_poly_derivative([7n], 11n)).toEqual([[], 1n]);
  expect(a).toEqual([1n, 0n, 3n, 0n]);
});
