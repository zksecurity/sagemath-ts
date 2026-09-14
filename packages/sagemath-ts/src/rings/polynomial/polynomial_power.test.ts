import { expect, test } from 'bun:test';
import { GF } from '../finite_rings/finite_field_constructor.js';
import { QQ } from '../rational_field.js';
import { Polynomial } from './polynomial_element.js';
import { PolynomialRing } from './polynomial_ring.js';

test('native polynomial negative powers and backend zero conventions', () => {
  const Q = new PolynomialRing(QQ, 'x');
  const f = Q.__call__([1n, 1n]);
  expect(f.pow(-2n).toString()).toBe('1/(x^2 + 2*x + 1)');
  expect(() => Q.zero().pow(-1n)).toThrow('negative exponent in power of zero');
  const F = new PolynomialRing(GF(7n), 'x');
  expect(F.zero().pow(-1n).toString()).toBe('0');
  expect(F.zero().pow(-1n)).toBeInstanceOf(Polynomial);
  expect(f.pow(1n)).not.toBe(f);
});

test('power overloads require narrowing for negative and dynamic exponents', () => {
  const f = new PolynomialRing(QQ, 'x').gen();
  const positive: typeof f = f.pow(2n);
  expect(positive.degree()).toBe(2);
  const negative = f.pow(-1n);
  // @ts-expect-error A negative result can belong to the fraction field.
  const negativePolynomial: typeof f = negative;
  const e: bigint = -1n;
  const dynamic = f.pow(e);
  // @ts-expect-error Dynamic exponents retain the fraction alternative.
  const dynamicPolynomial: typeof f = dynamic;
  const mixed = f.pow(-1n as -1n | 1n);
  // @ts-expect-error A literal union containing a negative value is not nonnegative.
  const mixedPolynomial: typeof f = mixed;
  expect([negativePolynomial, dynamicPolynomial, mixedPolynomial].map(String)).toEqual([
    '1/x',
    '1/x',
    '1/x',
  ]);
});
