import { expect, test } from 'bun:test';
import { _fmpq_poly_resultant } from './resultant.js';

// Direct native comparisons: polynomial_ops.cases.json / flint_fmpq_poly_resultant.
test('rational resultant clears content and denominators without changing its inputs', () => {
  const a = [2n, 4n, 0n],
    b = [-6n, -9n];
  expect(_fmpq_poly_resultant(a, 2n, b, 3n)).toEqual([-1n, 1n]);
  expect(_fmpq_poly_resultant(b, 3n, a, 2n)).toEqual([1n, 1n]);
  expect(_fmpq_poly_resultant([1n, 1n], 2n, [1n, -1n], 3n)).toEqual([1n, 3n]);
  expect(a).toEqual([2n, 4n, 0n]);
  expect(b).toEqual([-6n, -9n]);
});

test('rational resultant preserves constant, zero and shared-factor conventions', () => {
  expect(_fmpq_poly_resultant([], 1n, [2n], 3n)).toEqual([0n, 1n]);
  expect(_fmpq_poly_resultant([2n], 3n, [4n], 5n)).toEqual([1n, 1n]);
  expect(_fmpq_poly_resultant([2n], 3n, [1n, 0n, 1n], 7n)).toEqual([4n, 9n]);
  expect(_fmpq_poly_resultant([1n, 2n, 1n], 3n, [1n, 1n], 5n)).toEqual([0n, 1n]);
});
