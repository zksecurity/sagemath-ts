import { expect, test } from 'bun:test';
import { sympol_eval, fixedfieldsympol, fixedfieldorbits } from './galconj.js';
import { Fp_powu } from './arith1.js';

// Native galconj.c:883/898 reduces each Newton sum, then combines over Z.
test('weighted Newton sums preserve sign and values above the modulus', () => {
  const O = [[], [0n, 2n], [0n, 3n]];
  expect(sympol_eval({ v: [0, 3, -2], w: [0, 1, 2] }, O, 5n)).toEqual([0n, -2n, 1n]);
  expect(sympol_eval({ v: [0, 4], w: [0, 1] }, O, 5n)).toEqual([0n, 8n, 12n]);
});

test('zero weights return scalar zero but cancellation remains a column', () => {
  const O = [[], [0n, 2n]];
  expect(sympol_eval({ v: [0], w: [0] }, O, 5n)).toBe(0n);
  expect(sympol_eval({ v: [0, 0], w: [0, 1] }, O, 5n)).toBe(0n);
  expect(sympol_eval({ v: [0, 1, -1], w: [0, 1, 1] }, O, 5n)).toEqual([0n, 0n]);
});

test('Newton exponents use native unsigned casts and modular powering', () => {
  expect(sympol_eval({ v: [0, 1], w: [0, -1] }, [[], [0n, 2n]], 5n)).toEqual([0n, 3n]);
  expect(sympol_eval({ v: [0, 1], w: [0, 2 ** 30] }, [[], [0n, 2n]], 5n)).toEqual([0n, 1n]);
});

test('large-modulus first powers retain the signed remainder', () => {
  const p = (1n << 65n) + 13n;
  expect(sympol_eval({ v: [0, 1], w: [0, 1] }, [[], [0n, -2n]], p)).toEqual([0n, -2n]);
});

test('unsigned power preserves the native small-exponent shortcuts', () => {
  const p = (1n << 65n) + 13n;
  expect(Fp_powu(0n, 0n, 5n)).toBe(1n);
  expect(Fp_powu(p + 2n, 1n, p)).toBe(p + 2n);
  expect(Fp_powu(p + 2n, 2n, p)).toBe(4n);
  expect(Fp_powu(2n, -1n, 5n)).toBe(3n);
});

test('fixed-field orbit selection and separating weights match native', () => {
  const O = fixedfieldorbits([[], [0, 1, 2], [0, 3, 4]], [0n, 1n, 4n, 2n, 3n]);
  expect(O).toEqual([[], [0n, 1n, 4n], [0n, 2n, 3n]]);
  expect(fixedfieldsympol(O, 7n)).toEqual({ v: [0, 1], w: [0, 2] });
});
