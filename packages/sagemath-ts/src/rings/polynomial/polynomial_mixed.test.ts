import { expect, test } from 'bun:test';
import { QQ } from '../rational_field.js';
import { Rational } from '../rational.js';
import { PrimeField, FiniteFieldExtension } from '../finite_rings/finite_field_extension.js';
import { Zmod } from '../finite_rings/integer_mod_ring.js';
import { PolynomialRing } from './polynomial_ring.js';

test('mixed polynomials promote their parent before arithmetic and coefficient reduction', () => {
  const F = new PrimeField(7n),
    E = new FiniteFieldExtension(7n, 2);
  const f = new PolynomialRing(F).__call__([1n, 2n]);
  const h = new PolynomialRing(E).__call__([3n, 1n]);
  for (const g of [f.add(h), f.sub(h), f.mul(h)]) expect(g.parent.base_ring === E).toBe(true);
  expect(f.add(h).coeffs.map(String)).toEqual(['4', '3']);
  const Z = new PolynomialRing(Zmod(14n));
  expect(Z.__call__([0n, 7n]).eq(new PolynomialRing(F).zero())).toBe(true);
  expect(Z.__call__([0n, 7n]).add(new PolynomialRing(F).zero()).isZero()).toBe(true);
});
test('common quotient coefficients determine the result parent', () => {
  const f = new PolynomialRing(Zmod(6n)).__call__([1n, 2n]);
  const h = new PolynomialRing(Zmod(10n)).__call__([3n, 1n]);
  const g = f.mul(h);
  expect(String(g.parent)).toBe(
    'Univariate Polynomial Ring in x over Ring of integers modulo 2 (using GF2X)'
  );
  expect(g.coeffs.map(String)).toEqual(['1', '1']);
});
test('incompatible characteristic and variable names are rejected before zero shortcuts', () => {
  const Q = new PolynomialRing(QQ),
    F = new PolynomialRing(new PrimeField(7n));
  expect(() => Q.zero().add(F.zero())).toThrow('unsupported operand parent(s) for +');
  expect(Q.zero().eq(F.zero())).toBe(false);
  const y = new PolynomialRing(QQ, 'y').zero();
  expect(() => Q.zero().mul(y)).toThrow('unsupported operand parent(s) for *');
  expect(Q.zero().eq(y)).toBe(false);
});
test('nested coefficients use Sage sign, unit and parentheses formatting', () => {
  const T = new PolynomialRing(QQ, 't'),
    X = new PolynomialRing(T, 'x');
  const c = T.__call__([new Rational(1n), new Rational(-1n)]);
  expect(String(X.__call__([T.one(), c]))).toBe('(-t + 1)*x + 1');
  expect(String(T.__call__([1n, -2n, -1n]))).toBe('-t^2 - 2*t + 1');
  expect(String(X.__call__([T.zero(), T.gen().scalar_mul(new Rational(2n))]))).toBe('2*t*x');
});

// Permanent Sage comparisons for these examples live in polynomials.cases.json.
test('polynomial equality coerces integer, rational and finite scalars', () => {
  const Q = new PolynomialRing(QQ);
  expect(Q.one().eq(1n)).toBe(true);
  expect(Q.__call__([new Rational(1n, 2n)]).eq(new Rational(1n, 2n))).toBe(true);
  expect(Q.one().eq('1')).toBe(false);
  const F = new PrimeField(7n);
  expect(new PolynomialRing(Zmod(14n)).__call__([1n, 7n]).eq(F.one())).toBe(true);
  expect(Q.one().eq(F.one())).toBe(false);
});
test('polynomial comparisons with nonintegral numbers follow real-double coercion', () => {
  const Q = new PolynomialRing(QQ);
  expect(Q.__call__([new Rational(1n, 3n)]).eq(1 / 3)).toBe(true);
  expect(Q.__call__([1n << 1100n]).eq(Infinity)).toBe(true);
  expect(
    Q.__call__([new Rational(2n, 1n << 1075n), new Rational(1n, 1n << 1075n)]).eq(Number.MIN_VALUE)
  ).toBe(true);
  expect(Q.zero().eq(NaN)).toBe(false);
});
test('coefficient polynomials embed as constants in the outer polynomial ring', () => {
  const T = new PolynomialRing(QQ, 't'),
    X = new PolynomialRing(T, 'x');
  const t = T.gen(),
    f = X.__call__([t]);
  expect(f.eq(t)).toBe(true);
  expect(f.add(t).toString()).toBe('2*t');
  expect(t.mul(X.gen()).toString()).toBe('t*x');
  expect(f.quo_rem(t).map(String)).toEqual(['1', '0']);
  const z = X.zero();
  expect(z.mul(t)).toBe(z);
  expect(t.mul(z)).toBe(z);
});
