import { expect, test } from 'bun:test';
import { QQ } from '../rational_field.js';
import { Rational } from '../rational.js';
import { PrimeField, FiniteFieldExtension } from '../finite_rings/finite_field_extension.js';
import { Zmod } from '../finite_rings/integer_mod_ring.js';
import { PolynomialRing } from './polynomial_ring.js';

// Every case below is also compared with the original Sage implementation.
test('GCD zero shortcuts retain nonmonic native finite operands across equivalent parents', () => {
  const F = new PrimeField(7n),
    R = new PolynomialRing(F),
    S = new PolynomialRing(F);
  const f = R.__call__([2n, 2n]);
  expect(S.zero().gcd(f)).toBe(f);
  expect(f.gcd(S.zero())).toBe(f);
  expect(f.gcd(S.__call__([2n, 2n])).coeffs.map(String)).toEqual(['1', '1']);
});
test('GCD rational normalization allocates and promotes before division', () => {
  const R = new PolynomialRing(QQ);
  const f = R.__call__([new Rational(2n, 3n), new Rational(1n)]);
  const g = f.gcd(R.zero());
  expect(g.eq(f)).toBe(true);
  expect(g).not.toBe(f);
});
test('composite GCD distinguishes equal nonunit inputs from native calculation failure', () => {
  const R = new PolynomialRing(Zmod(14n));
  const f = R.__call__([2n, 2n]);
  expect(() => f.gcd(f)).toThrow('leading coefficient must be invertible');
  expect(() => R.__call__([1n, 3n, 1n]).gcd(R.__call__([1n, 1n, 1n]))).toThrow(
    'FLINT gcd calculation failed'
  );
  expect(
    R.one()
      .gcd(R.__call__([2n]))
      .coeffs.map(String)
  ).toEqual(['1']);
});
test('extension GCD allocates for nontrivial division even when the divisor is monic', () => {
  const R = new PolynomialRing(new FiniteFieldExtension(7n, 2));
  const f = R.__call__([1n, 2n, 1n]),
    h = R.__call__([1n, 1n]);
  const g = f.gcd(h);
  expect(g.eq(h)).toBe(true);
  expect(g).not.toBe(h);
});
