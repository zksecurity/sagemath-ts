import { expect, test } from 'bun:test';
import { QQ } from '../rational_field.js';
import { Rational } from '../rational.js';
import { PrimeField } from '../finite_rings/finite_field_extension.js';
import { Zmod } from '../finite_rings/integer_mod_ring.js';
import { PolynomialRing } from './polynomial_ring.js';
import { NTLError } from '../../errors.js';

test('word-modular polynomial division validates the divisor before degree shortcuts', () => {
  const R = new PolynomialRing(Zmod(14n));
  expect(() => R.zero().quo_rem(R.__call__([1n, 2n]))).toThrow(
    'Leading coefficient of a must be invertible.'
  );
  const F = new PolynomialRing(new PrimeField(7n));
  const error = (() => {
    try {
      F.one().quo_rem(F.zero());
    } catch (e) {
      return e as Error;
    }
  })();
  expect(error?.name).toBe('ZeroDivisionError');
  expect(error?.message).toBe('');
});
test('division preserves generic zero identity and native NTL errors', () => {
  const R = new PolynomialRing(QQ),
    zero = R.zero();
  const [q, r] = zero.quo_rem(R.gen());
  expect(q).toBe(zero);
  expect(r).toBe(zero);
  const N = new PolynomialRing(new PrimeField((1n << 127n) - 1n));
  expect(() => N.one().quo_rem(N.zero())).toThrow(NTLError);
  expect(() => N.one().quo_rem(N.zero())).toThrow('ZZ_pX: division by zero');
});
test('nested polynomial coefficients support exact division and reject fractional coefficients', () => {
  const T = new PolynomialRing(QQ, 't'),
    X = new PolynomialRing(T, 'x');
  const c = T.gen().add(T.one());
  const [q, r] = X.__call__([c, c]).quo_rem(X.__call__([c]));
  expect(String(q)).toBe('x + 1');
  expect(r.isZero()).toBe(true);
  expect(() => X.one().quo_rem(X.__call__([c]))).toThrow('division non exact');
  const half = T.__call__(new Rational(1n, 2n));
  expect(String(X.one().quo_rem(X.__call__([half]))[0])).toBe('2');
});
