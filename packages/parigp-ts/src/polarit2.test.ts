import { expect, test } from 'bun:test';
import { dbltor, rtodbl } from './kernel/none/mp_indep.js';
import { resultant } from './polarit2.js';

// These bit-level cases are also compared directly to libpari in polynomial_ops.
test('native PARI double conversion preserves exponent cutoffs and signed-zero loss', () => {
  expect(Object.is(rtodbl(dbltor(-0)), 0)).toBe(true);
  expect(rtodbl(dbltor(Number.MIN_VALUE))).toBe(0);
  expect(rtodbl(dbltor(2 ** -1022))).toBe(2 ** -1022);
  expect(rtodbl(dbltor(1 / 3))).toBe(1 / 3);
  expect(() => rtodbl(dbltor(Number.MAX_VALUE))).toThrow('overflow');
  expect(() => dbltor(Number.NaN)).toThrow('overflow');
});

test('native real resultant handles constants in both positions and converts before zero checks', () => {
  expect(resultant([1, 2, 3, 4], [0.5])).toBe(0.125);
  expect(resultant([0.5], [1, 2, 3, 4])).toBe(0.125);
  expect(resultant([2], [3])).toBe(1);
  expect(resultant([], [3])).toBe(0);
  expect(resultant([1, 1], [Number.MIN_VALUE])).toBe(0);
  expect(() => resultant([], [Number.NaN])).toThrow('overflow');
});

test('native real resultant retains inexact leading zeros and their polynomial degree', () => {
  expect(resultant([1, 0], [1, 1])).toBe(-1);
  expect(resultant([1, 0, 0], [1, 1])).toBe(1);
  expect(resultant([0, 0, 1, 0], [1, 1])).toBe(-1);
  expect(resultant([0, 0], [1, 1])).toBe(0);
  expect(resultant([1, 0, 1], [-2, 0, 0, 1])).toBe(5);
});
