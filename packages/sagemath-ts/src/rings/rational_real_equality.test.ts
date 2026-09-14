import { expect, test } from 'bun:test';
import { Rational } from './rational.js';

test('rational equality accepts fractional host numbers through RDF', () => {
  expect(new Rational(1n, 2n).eq(0.5)).toBe(true);
  expect(new Rational(1n, 10n).eq(0.1)).toBe(true);
  expect(new Rational(1n, 10n).eq(0.2)).toBe(false);
});

test('rational equality applies native real rounding to large integral numbers', () => {
  const q = new Rational(2n ** 100n + 1n);
  expect(q.eq(2 ** 100)).toBe(true);
  expect(q.eq(2n ** 100n)).toBe(false);
});

test('rational equality follows real underflow and signed zero', () => {
  expect(new Rational(1n, 2n ** 2000n).eq(0)).toBe(true);
  expect(new Rational(-1n, 2n ** 2000n).eq(-0)).toBe(true);
  expect(Rational.zero().eq(-0)).toBe(true);
});

test('rational equality follows real overflow and infinity signs', () => {
  expect(new Rational(2n ** 2000n).eq(Infinity)).toBe(true);
  expect(new Rational(-(2n ** 2000n)).eq(-Infinity)).toBe(true);
  expect(new Rational(2n ** 2000n).eq(-Infinity)).toBe(false);
});

test('rational equality with NaN is false', () => {
  expect(Rational.zero().eq(NaN)).toBe(false);
  expect(new Rational(2n ** 2000n).eq(NaN)).toBe(false);
});
