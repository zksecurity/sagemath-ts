import { expect, test } from 'bun:test';
import { Z_pvalrem } from './gen2.js';

test('native valuation/unit threshold and signs', () => {
  expect(Z_pvalrem(-11n * 3n ** 127n, 3n)).toEqual([127, -11n]);
  expect(Z_pvalrem(11n << 256n, 2n)).toEqual([256, 11n]);
  expect(Z_pvalrem(7n, 3n)).toEqual([0, 7n]);
  expect(Z_pvalrem(-5n * 18446744073709551629n ** 17n, 18446744073709551629n)).toEqual([17, -5n]);
});

test('valuation primitive enforces the native nonzero, nonunit precondition', () => {
  expect(() => Z_pvalrem(0n, 2n)).toThrow('Z_pvalrem requires nonzero n and p > 1');
  expect(() => Z_pvalrem(3n, 1n)).toThrow('Z_pvalrem requires nonzero n and p > 1');
});
