import { expect, test } from 'bun:test';
import { mpz_remove } from './gmp.js';

test('GMP factor removal preserves signs, composite bases and large exponents', () => {
  expect(mpz_remove(0n, 3n)).toEqual([0n, 0n]);
  expect(mpz_remove(-11n * 3n ** 4096n, 3n)).toEqual([4096n, -11n]);
  expect(mpz_remove(2n ** 25n * 3n ** 17n, 6n)).toEqual([17n, 256n]);
  expect(mpz_remove(-5n << 262144n, 2n)).toEqual([262144n, -5n]);
  expect(() => mpz_remove(1n, 1n)).toThrow(
    'You can only compute the valuation with respect to a integer larger than 1.'
  );
});
