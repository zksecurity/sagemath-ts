import { expect, test } from 'bun:test';
import { _fmpq_poly_gcd } from './gcd.js';

// Sage comparisons for these kernels are in polynomial_ops.cases.json.
test('rational GCD removes numerator content and returns a positive common denominator', () => {
  const a = [-6n, -9n, 0n],
    b = [4n, 6n];
  expect(_fmpq_poly_gcd(a, b)).toEqual([[2n, 3n], 3n]);
  expect(a).toEqual([-6n, -9n, 0n]);
  expect(b).toEqual([4n, 6n]);
  expect(_fmpq_poly_gcd([], a)).toEqual([[2n, 3n], 3n]);
  expect(_fmpq_poly_gcd([], [])).toEqual([[], 1n]);
});
