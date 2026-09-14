import { expect, test } from 'bun:test';
import { fmpz_get_d } from './get.js';

test('bundled FLINT truncates at the mantissa boundary instead of rounding', () => {
  // Native fmpz_get_d, permanently compared in polynomial_ops.
  expect(fmpz_get_d(9007199254740995n)).toBe(9007199254740994);
  expect(fmpz_get_d(-9007199254740995n)).toBe(-9007199254740994);
  expect(fmpz_get_d((1n << 1024n) - 1n)).toBe(Number.MAX_VALUE);
  expect(fmpz_get_d(1n << 1024n)).toBe(Infinity);
});
