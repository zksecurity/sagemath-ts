import { expect, test } from 'bun:test';
import { subgrouplist } from './subgroup.js';
import { absZ_factor_limit_strict } from './ifactor1.js';
import { absZ_factor_limit_strict_default } from './ifactor.js';
import { WordPrimeIterator } from './language/forprime.js';

test('unbounded and scalar-bound subgroup lists retain native order', () => {
  expect(subgrouplist([4n])).toEqual([[[4n]], [[2n]], [[1n]]]);
  expect(subgrouplist([4n], 2n)).toEqual([[[2n]], [[1n]]]);
  expect(subgrouplist([4n], [2n])).toEqual([[[2n]]]);
});

test('strict word and multiword cofactors retain different native power detection', () => {
  expect(absZ_factor_limit_strict(3n ** 11n, 1n)).toEqual([[], [3n ** 11n, 1n]]);
  expect(absZ_factor_limit_strict(3n ** 39n, 1n)).toEqual([[], [3n ** 13n, 3n]]);
  expect(absZ_factor_limit_strict(3n ** 41n, 1n)).toEqual([[], [3n, 41n]]);
});

test('cached primes can identify a pure power beyond the requested factor limit', () => {
  expect(absZ_factor_limit_strict(509n ** 3n, 257n)).toEqual([[[509n, 3n]], null]);
  expect(absZ_factor_limit_strict(509n ** 3n, -257n)).toEqual([[[509n, 3n]], null]);
});

test('default factorization shares the strict engine without eager prime-table cycles', () => {
  expect(absZ_factor_limit_strict_default(81n)).toEqual([[[3n, 4n]], null]);
  const primes = new WordPrimeIterator(251n);
  expect([primes.next(), primes.next(), primes.next()]).toEqual([251n, 257n, 263n]);
});

test('native bounds and conversion errors remain exact', () => {
  for (const [run, message] of [
    [() => subgrouplist([2n], []), 'incorrect type in subgroup (t_VEC).'],
    [() => subgrouplist([2n], [1n, 2n]), 'incorrect type in subgroup (t_VEC).'],
    [() => subgrouplist([2n], 0n), 'domain error in subgroup: index bound <= 0'],
    [() => absZ_factor_limit_strict(0n, 1n << 64n), 'overflow in t_INT-->ulong assignment.'],
  ] as const) {
    try {
      run();
      throw new Error('expected native rejection');
    } catch (e) {
      expect((e as Error).name).toBe('PariError');
      expect((e as Error).message).toBe(message);
    }
  }
});
