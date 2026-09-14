import { expect, test } from 'bun:test';
import { _nmod_poly_xgcd } from './xgcd.js';

// These coefficients and cutoff branches are compared directly with native FLINT.
test('modular XGCD handles a constant remainder above the half-GCD cutoff', () => {
  const b = [2n, 3n, ...Array<bigint>(97).fill(0n), 1n];
  const a = [...b];
  a[0] = 3n;
  const before = [...a];
  expect(_nmod_poly_xgcd(a, b, 7n)).toEqual([[1n], [1n], [6n]]);
  expect(_nmod_poly_xgcd(b, a, 7n)).toEqual([[1n], [6n], [1n]]);
  expect(a).toEqual(before);
});
test('native modular XGCD zero and nonunit boundaries', () => {
  expect(_nmod_poly_xgcd([], [], 7n)).toEqual([[], [], []]);
  expect(_nmod_poly_xgcd([2n, 2n], [], 7n)).toEqual([[1n, 1n], [4n], []]);
  expect(() => _nmod_poly_xgcd([1n], [2n], 14n)).toThrow('coefficient is not invertible');
});
