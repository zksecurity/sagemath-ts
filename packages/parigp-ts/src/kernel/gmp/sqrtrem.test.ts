import { expect, test } from 'bun:test';
import { sqrtr, sqrtr_abs, sqrti, itor, real_0_bit } from '../../qfb.js';
import { sqrti as buchSqrti } from '../../buch.js';
import { sqrtremi } from './mp.js';

test('native even-exponent root retains its last-guard correction', () => {
  expect(sqrtr_abs({ s: 1, e: 0, p: 64, m: (1n << 63n) + 1n })).toEqual({
    s: 1, e: 0, p: 64, m: (1n << 63n) + 1n,
  });
});

test('negative real roots have an exact zero real component', () => {
  expect(sqrtr(itor(-4n, 64))).toEqual({ re: 0n, im: { s: 1, e: 1, p: 64, m: 1n << 63n } });
});

test('zero root accuracy keeps the full exponent before signed halving', () => {
  expect(sqrtr(real_0_bit(-(2 ** 32) - 1))).toEqual({ s: 0, e: -2147483649, m: 0n, p: 0 });
});

test('native integer kernels read the magnitude and retain exact remainders', () => {
  expect(sqrti(-9n)).toBe(3n);
  expect(buchSqrti(-9n)).toBe(3n);
  expect(sqrtremi((1n << 64n) + 1n)).toEqual([1n << 32n, 1n]);
  expect(sqrtremi((1n << 4096n) - 1n)).toEqual([(1n << 2048n) - 1n, (1n << 2049n) - 2n]);
});
