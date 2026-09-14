import { expect, test } from 'bun:test';
import { _fmpz_poly_xgcd } from './xgcd.js';

// Expected coefficients are compared independently with native FLINT in polynomial_ops.
test('integer Bezout coefficients preserve native length ordering and input arrays', () => {
  const a = [2n, 1n], b = [3n, 0n, 0n, 1n];
  expect(_fmpz_poly_xgcd(a, b)).toEqual([5n, [4n, -2n, 1n], [-1n]]);
  expect(_fmpz_poly_xgcd(b, a)).toEqual([5n, [-1n], [4n, -2n, 1n]]);
  expect(a).toEqual([2n, 1n]);
  expect(b).toEqual([3n, 0n, 0n, 1n]);
});
test('integer native XGCD reports a vanishing resultant and rejects undefined constants', () => {
  expect(_fmpz_poly_xgcd([2n, 2n], [4n, 6n, 2n])).toEqual([0n, [], []]);
  expect(_fmpz_poly_xgcd([], [1n, 1n])).toEqual([0n, [], []]);
  expect(() => _fmpz_poly_xgcd([1n], [2n])).toThrow(RangeError);
});
