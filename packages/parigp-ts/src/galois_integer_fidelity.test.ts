import { expect, spyOn, test } from 'bun:test';
import { ugcd, ulcm, factoru_small, Forprime, logint } from './galconj.js';
import * as factors from './ifactor.js';
import { PariError } from './errors.js';

// Native counterparts, including prior timeouts: ff_pari_galois_integer.
test('word helpers preserve magnitude conversion and native lcm overflow', () => {
  expect(ugcd(-6, 17)).toBe(1);
  expect(ulcm(-6, 17)).toBe(102);
  expect(ulcm(3, 2 ** 63)).toBe(2 ** 63);
  expect(ulcm(0, 0)).toBe(0);
  expect(() => ugcd(0, 2 ** 64)).toThrow(new PariError('overflow in t_INT-->ulong assignment.'));
  expect(() => ulcm(0, 2 ** 64)).toThrow(new PariError('overflow in t_INT-->ulong assignment.'));
});

test('small factorization delegates to the backend, including large primes and zero', () => {
  const factor = spyOn(factors, 'factoru');
  try {
    expect(factoru_small(9007199254740881)).toEqual([[9007199254740881, 1]]);
    expect(factor.mock.calls[0]?.[0]).toBe(9007199254740881n);
    expect(factoru_small(0)).toEqual([[0, 1]]);
    expect(factoru_small(-17)).toEqual([[17, 1]]);
  } finally {
    factor.mockRestore();
  }
});

test('integer logarithms validate the base first and use exact power boundaries', () => {
  expect(() => logint(1n, 1n)).toThrow(new PariError('domain error in logint: b <= 1'));
  expect(() => logint(0n, 0n)).toThrow(new PariError('domain error in logint: b <= 1'));
  expect(() => logint(0n, 2n)).toThrow(new PariError('domain error in logint: x <= 0'));
  expect(logint(3n ** 8192n - 1n, 3n)).toBe(8191);
  expect(logint(3n ** 8192n, 3n)).toBe(8192);
  expect(logint(3n ** 8192n + 1n, 3n)).toBe(8192);
});

test('prime iterator retains integer progress above the exact-number boundary', () => {
  const p = new Forprime(2 ** 53);
  expect(Array.from({ length: 4 }, () => BigInt(p.next()))).toEqual([
    9007199254740996n,
    9007199254741032n,
    9007199254741048n,
    9007199254741068n,
  ]);
  const end = new Forprime(Number((1n << 64n) - 2048n));
  const values = Array.from({ length: 100 }, () => end.next());
  expect(values[0]).toBeGreaterThan(0);
  expect(values[99]).toBe(0);
  expect(end.next()).toBe(0);
});

test('prime iterator reproduces the initialized table transition', () => {
  expect(new Forprime(499978).next()).toBe(499979);
  expect(new Forprime(499979).next()).toBe(500009);
  expect(new Forprime(-17).next()).toBe(17);
  expect(() => new Forprime(2 ** 64)).toThrow(
    new PariError('overflow in t_INT-->ulong assignment.')
  );
});
