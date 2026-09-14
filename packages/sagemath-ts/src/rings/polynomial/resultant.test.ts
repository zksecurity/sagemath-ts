import { expect, test } from 'bun:test';
import { PrimeField } from '../finite_rings/finite_field_extension.js';
import { Zmod } from '../finite_rings/integer_mod_ring.js';
import { Integer } from '../integer_ring.js';
import { QQ } from '../rational_field.js';
import { Rational } from '../rational.js';
import { RDF } from '../real_double.js';
import { PolynomialRing } from './polynomial_ring.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';

const ZZ = {
  zero: () => new Integer(0n),
  one: () => new Integer(1n),
  __call__: (x: unknown) => new Integer(x as bigint),
  is_field: () => false,
  toString: () => 'Integer Ring',
} as unknown as CoefficientRing<RingElement>;

// The permanent mixed-parent/scalar matrix in polynomials.cases.json calls Sage.
test('resultant promotes rational scalars and accepts proof only on the integer backend', () => {
  const R = new PolynomialRing(ZZ, 'x');
  const f = R.__call__([1n, 0n, 1n]);
  expect(String(f.resultant(new Rational(1n, 2n)))).toBe('1/4');
  expect(String(f.resultant(2n, { proof: false }))).toBe('4');
  expect(() => f.resultant(new Rational(1n, 2n), { proof: false })).toThrow(
    "unexpected keyword argument 'proof'"
  );
});

test('word resultant retains native constant aliasing and composite zero errors', () => {
  const R = new PolynomialRing(new PrimeField(7n), 'x');
  const c = R.__call__([3n]);
  expect(String(c.resultant(c))).toBe('0');
  expect(String(c.resultant(R.__call__([3n])))).toBe('1');
  const C = new PolynomialRing(Zmod(14n), 'x');
  expect(() => C.zero().resultant(C.one())).toThrow(
    'The Sylvester matrix is not defined for zero polynomials'
  );
});

test('real scalar resultants use PARI conversion and preserve the resulting RDF parent', () => {
  const R = new PolynomialRing(QQ, 'x');
  const f = R.__call__([1n, 0n, 1n]);
  const result = f.resultant(0.5);
  expect(String(result)).toBe('0.25');
  expect((result as unknown as { parent: unknown }).parent).toBe(RDF);
  expect(String(R.__call__([1n, 1n]).resultant(Number.MIN_VALUE))).toBe('0.0');
  expect(String(R.one().resultant(Number.NaN))).toBe('1.0');
  expect(() => f.resultant(Number.NaN)).toThrow('array must not contain infs or NaNs');
});

test('Sage RDF resultants normalize trailing zeros before PARI conversion', () => {
  const R = new PolynomialRing(RDF, 'x');
  expect(String(R.__call__([1, 0]).resultant(R.__call__([1, 1])))).toBe('1.0');
  expect(String(R.__call__([1, 0, 1]).resultant(R.__call__([-2, 0, 0, 1])))).toBe('5.0');
});
