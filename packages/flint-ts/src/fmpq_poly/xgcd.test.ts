import { expect, test } from 'bun:test';
import { _fmpq_poly_xgcd } from './xgcd.js';

test('rational native XGCD removes common factors and restores rational scaling', () => {
  expect(_fmpq_poly_xgcd([2n, 2n], 6n, [4n, 6n, 2n], 10n)).toEqual([
    [[1n, 1n], 1n], [[3n], 1n], [[], 1n],
  ]);
  expect(_fmpq_poly_xgcd([2n, 2n], 6n, [], 1n)).toEqual([
    [[1n, 1n], 1n], [[3n], 1n], [[], 1n],
  ]);
});
test('rational native zero outputs are independently allocated and denominators validated', () => {
  const [g, s, t] = _fmpq_poly_xgcd([], 1n, [], 1n);
  expect(g).toEqual([[], 1n]);
  expect(g[0]).not.toBe(s[0]);
  expect(s[0]).not.toBe(t[0]);
  expect(() => _fmpq_poly_xgcd([1n], 0n, [1n], 1n)).toThrow(RangeError);
});
