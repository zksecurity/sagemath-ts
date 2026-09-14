import { expect, test } from 'bun:test';
import { _fmpz_poly_divrem } from './divrem.js';

test('FLINT divrem preserves signed floor quotients and the small-magnitude skip', () => {
  expect(_fmpz_poly_divrem([-7n, -5n], [3n])).toEqual([
    [-3n, -2n],
    [2n, 1n],
  ]);
  expect(_fmpz_poly_divrem([7n, 5n], [-3n])).toEqual([
    [-3n, -2n],
    [-2n, -1n],
  ]);
  expect(_fmpz_poly_divrem([-1n], [2n])).toEqual([[], [-1n]]);
  expect(_fmpz_poly_divrem([0n, 0n, 1n], [-1n, 2n])).toEqual([[], [0n, 0n, 1n]]);
});
test('FLINT exact flag checks leading quotients without demanding a zero remainder', () => {
  const a = Object.freeze([1n, 2n]),
    b = Object.freeze([0n, 2n]);
  expect(_fmpz_poly_divrem(a, b, true)).toEqual([[1n], [1n]]);
  expect(_fmpz_poly_divrem([1n, 3n], b, true)).toBeNull();
  expect(a).toEqual([1n, 2n]);
  expect(b).toEqual([0n, 2n]);
});
