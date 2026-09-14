import { expect, test } from 'bun:test';
import { gmp_primesieve, mpz_fac_ui, mpz_oddfac_1 } from './gmp_factorial.js';

test('GMP factorial native preconditions reject invalid unsigned inputs and flags', () => {
  for (const n of [-1n, 1n << 64n]) {
    expect(() => mpz_fac_ui(n)).toThrow(RangeError);
    expect(() => mpz_oddfac_1(n)).toThrow(RangeError);
    expect(() => gmp_primesieve(n)).toThrow(RangeError);
  }
  expect(() => gmp_primesieve(4n)).toThrow(RangeError);
  expect(() => gmp_primesieve((1n << 64n) - 1n)).toThrow('Invalid typed array length');
  expect(() => mpz_oddfac_1(251n, 1)).toThrow(RangeError);
  expect(() => mpz_oddfac_1(252n, 2 as 1)).toThrow(RangeError);
});
