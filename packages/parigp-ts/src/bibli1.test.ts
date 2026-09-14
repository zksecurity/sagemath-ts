import { expect, test } from 'bun:test';
import { algdep } from './bibli1.js';

test('PARI algdep double conversion, raw polynomial and degree dispatch', () => {
  expect(algdep(Math.sqrt(2), 2n)).toEqual([-2n, 0n, 1n]);
  expect(algdep(17 / 9, 1n)).toEqual([-17n, 9n]);
  expect(algdep(0, -2n)).toEqual([0n, 1n]);
  expect(algdep(1, 0n)).toEqual([1n]);
  expect(() => algdep(1, -1n)).toThrow('domain error in algdep: degree < 0');
  expect(() => algdep(NaN, 2n)).toThrow('overflow in dbltor [NaN or Infinity]');
  expect(() => algdep(1, 1n << 63n)).toThrow('Python int too large to convert to C long');
});
