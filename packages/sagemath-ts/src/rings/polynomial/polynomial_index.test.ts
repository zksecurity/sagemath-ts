import { expect, test } from 'bun:test';
import { Integer } from '../integer_ring.js';
import { Rational } from '../rational.js';
import { QQ } from '../rational_field.js';
import { PrimeField, FiniteFieldExtension } from '../finite_rings/finite_field_extension.js';
import { PolynomialRing } from './polynomial_ring.js';

// Every behavior below is also executed against Sage in polynomials.cases.json.
test('polynomial indices distinguish exact indexing from Cython truncation', () => {
  const f = new PolynomialRing(QQ).__call__([1n, 2n, 3n]);
  expect(String(f.getCoeff(new Integer(1n)))).toBe('2');
  expect(() => f.getCoeff(new Rational(3n, 2n))).toThrow(
    'unable to convert rational 3/2 to an integer'
  );
  expect(f.truncate(new Rational(3n, 2n)).coeffs.map(String)).toEqual(['1']);
  expect(f.reverse(new Rational(3n, 2n)).coeffs.map(String)).toEqual(['2', '1']);
  expect(f.reverse(-1n).isZero()).toBe(true);
  expect(() => f.getCoeff(2n ** 63n)).toThrow('Sage Integer too large to convert to C long');
});
test('generic shifts retain Sage zero shortcuts and asymmetric fractional handling', () => {
  const R = new PolynomialRing(QQ),
    f = R.__call__([1n, 2n, 3n]);
  const zero = R.zero();
  expect(zero.shift('bad')).toBe(zero);
  expect(f.shift(new Rational(-3n, 2n))?.coeffs.map(String)).toEqual(['2', '3']);
  expect(() => f.shift(new Rational(2n))).toThrow('unsupported operand parent(s) for *');
  expect(f.shift(NaN)).toBeNull();
  expect(() => f.shift(2n ** 80n)).toThrow('index-sized integer');
});
test('native and generic backends preserve truncation and shift identity', () => {
  const F = new PrimeField(7n),
    f = new PolynomialRing(F).__call__([1n, 2n, 3n]);
  expect(f.shift(0n)).toBe(f);
  expect(f.truncate(3n)).toBe(f);
  expect(() => f.shift(2n ** 31n)).toThrow('value too large to convert to int');
  const E = new FiniteFieldExtension(7n, 2);
  const g = new PolynomialRing(E).__call__([1n, 2n, 3n]);
  expect(g.shift(0n)).not.toBe(g);
  const h = new PolynomialRing(new PolynomialRing(QQ, 'y'), 'x').__call__([1n, 2n, 3n]);
  expect(h.truncate(-1n).coeffs.map(String)).toEqual(['1', '2']);
  expect(h.truncate(5n)).not.toBe(h);
});
test('reverse preserves finite scalar arithmetic before converting its degree', () => {
  const F = new PrimeField(7n);
  const f = new PolynomialRing(QQ).__call__([1n, 2n, 3n]);
  expect(f.reverse(F.__call__(-1n)).isZero()).toBe(true);
  const g = new PolynomialRing(new PrimeField(2n)).__call__([1n, 0n, 1n, 1n, 1n, 1n, 1n, 1n, 1n]);
  // len(g)=9 becomes 2 modulo 7; comparing with degree+1=4 pads two zeros.
  expect(g.reverse(F.__call__(3n)).degree()).toBe(10);
  const a = new FiniteFieldExtension(7n, 2).gen();
  expect(() => f.truncate(a)).toThrow('element is not in the prime field');
  expect(() => f.getCoeff(a)).toThrow('cannot be interpreted as an integer');
});
