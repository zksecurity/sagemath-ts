import { expect, test } from 'bun:test';
import { is_357_power, is_pth_power, Z_isanypower, forprime } from './ifactor.js';

test('negative perfect powers retain only odd exponents and a negative base', () => {
  expect(Z_isanypower(-64n)).toEqual([3, -4n]);
  expect(Z_isanypower(-729n)).toEqual([3, -9n]);
  expect(Z_isanypower(-27n)).toEqual([3, -3n]);
  expect(Z_isanypower(-16n)).toEqual([0, -16n]);
  expect(Z_isanypower(-1n)).toEqual([0, -1n]);
  expect(Z_isanypower(-37n)).toEqual([0, -37n]);
});
test('word and multiword perfect-power searches have native priorities', () => {
  expect(is_357_power(3n ** 21n, 7)).toEqual([3, 3n ** 7n, 5]);
  expect(is_357_power(103n ** 21n, 7)).toEqual([7, 103n ** 3n, 5]);
});
test('native residue and valuation filters update the returned power mask', () => {
  expect(is_357_power(32n, 7)).toEqual([5, 2n, 2]);
  expect(is_357_power(128n, 7)).toEqual([7, 2n, 4]);
  expect(is_357_power(243n, 3)).toEqual([5, 3n, 2]);
});
test('an exhausted exponent search returns native zero and consumes the cutoff prime', () => {
  const primes = forprime(11, 100);
  let last = 0;
  const iterator = {
    next: () => {
      const value = primes.next();
      last = value.done ? 0 : value.value;
      return value.done ? null : value.value;
    },
  };
  expect(is_pth_power(103n, iterator, 1)).toEqual([0, 103n]);
  expect(last).toBe(11);
  expect(is_pth_power(103n ** 13n, iterator, 1)).toEqual([13, 103n]);
  expect(last).toBe(13);
});
