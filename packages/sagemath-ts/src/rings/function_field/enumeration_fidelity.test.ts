import { expect, test } from 'bun:test';
import type { PrimeField, PrimeFieldElement } from '../finite_rings/finite_field_extension.js';
import type { RationalFunctionField_global } from './function_field_rational.js';
import { GF } from '../finite_rings/index.js';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
import { FunctionField } from './index.js';

test('a place prefix never requests the materialized constant field list', () => {
  const k = GF(2305843009213693951n) as PrimeField;
  const previous = Object.getOwnPropertyDescriptor(k, 'list');
  Object.defineProperty(k, 'list', {
    configurable: true,
    value() {
      throw new Error('eager field list');
    },
  });
  try {
    const K = FunctionField(k, 'x') as RationalFunctionField_global<PrimeFieldElement>;
    expect(String(K.get_place(1))).toBe('Place (x)');
    const places = K._places_finite(1);
    expect(String(places.next().value)).toBe('Place (x)');
    expect(String(places.next().value)).toBe('Place (x + 1)');
  } finally {
    if (previous) Object.defineProperty(k, 'list', previous);
    else delete (k as unknown as { list?: unknown }).list;
  }
});

test('polynomial enumeration preserves bounds and independent coefficient rows', () => {
  const R = new PolynomialRing(GF(3n) as PrimeField, 'x');
  const g = R.polynomials({ max_degree: 1 });
  const first = g.next().value!;
  expect(String(first)).toBe('0');
  expect(String(g.next().value)).toBe('1');
  expect(String(first)).toBe('0');
  expect([...R.monics({ max_degree: 1 })].map(String)).toEqual(['1', 'x', 'x + 1', 'x + 2']);
  expect([...R.polynomials({ max_degree: -1 })].map(String)).toEqual(['0']);
  expect(() => R.polynomials()).toThrow('you should pass exactly one of of_degree and max_degree');
});

test('multivariate display and projective constructor preserve native errors', async () => {
  const { MPolynomialRing } = await import('../polynomial/multi_polynomial_ring.js');
  const { HyperellipticCurve } = await import('../../schemes/hyperelliptic_curves/constructor.js');
  const k = GF(5n) as PrimeField;
  const R = new MPolynomialRing(k, ['x', 'y']);
  expect(String(R.monomial([1, 0], k.__call__(3n)))).toBe('-2*x');
  const T = new PolynomialRing(k, 'x');
  expect(() => HyperellipticCurve(T.one())).toThrow('(x1^2 - x2^2)/x2^2 cannot be converted to a polynomial');
});
