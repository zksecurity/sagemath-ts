import { expect, test } from 'bun:test';
import { _fmpz_poly_resultant } from './resultant.js';

// Native FLINT comparisons cover these formulas and their dispatched kernels.
test('integer resultant preserves orientation, content and zero conventions', () => {
  expect(_fmpz_poly_resultant([2n, 4n], [-6n, -9n])).toBe(-6n);
  expect(_fmpz_poly_resultant([-6n, -9n], [2n, 4n])).toBe(6n);
  expect(_fmpz_poly_resultant([-2n], [1n, 0n, 0n, 1n])).toBe(-8n);
  expect(_fmpz_poly_resultant([], [1n])).toBe(0n);
});
test('modular resultant reconstruction skips primes dividing leading coefficients', () => {
  const p = (1n << 63n) + 29n,
    q = (1n << 63n) + 99n;
  const a = [1n, ...Array<bigint>(143).fill(0n), p];
  const b = [2n, ...Array<bigint>(143).fill(0n), q];
  expect(_fmpz_poly_resultant(a, b)).toBe((2n * p - q) ** 144n);
  expect(a[0]).toBe(1n);
  expect(b[0]).toBe(2n);
});
