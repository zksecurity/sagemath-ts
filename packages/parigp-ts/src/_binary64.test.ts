import { expect, test } from 'bun:test';
import { fma, frexp, ldexp } from './_binary64.js';

test('native fma rounds only once across cancellation and underflow', () => {
  expect(fma(1 + 2 ** -52, 1 - 2 ** -52, -1)).toBe(-(2 ** -104));
  expect(fma(2 ** -1074, 0.5, 2 ** -1074)).toBe(2 ** -1073);
  expect(fma(1e308, 1e308, -Infinity)).toBe(-Infinity);
  expect(fma(-1e308, 1e308, Infinity)).toBe(Infinity);
  expect(fma(Infinity, 0, 1)).toBeNaN();
  expect(Object.is(fma(-0, 1, -0), -0)).toBe(true);
});

test('native ldexp rounds subnormal ties once and preserves signed zero', () => {
  expect(ldexp(1 + 2 ** -52, -1075)).toBe(2 ** -1074);
  expect(ldexp(1, -1075)).toBe(0);
  expect(Object.is(ldexp(-1, -1075), -0)).toBe(true);
  expect(ldexp(Number.MAX_VALUE, 1)).toBe(Infinity);
  expect(ldexp(2 ** -1074, 1074)).toBe(1);
});

test('native frexp includes subnormals and exceptional values', () => {
  expect(frexp(2 ** -1074)).toEqual([0.5, -1073]);
  expect(frexp(-3)).toEqual([-0.75, 2]);
  expect(frexp(Infinity)).toEqual([Infinity, 0]);
  expect(frexp(-0)).toEqual([-0, 0]);
});
