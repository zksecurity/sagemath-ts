import { expect, test } from 'bun:test';
import { QQ } from '../rational_field.js';
import { PrimeField } from '../finite_rings/finite_field_extension.js';
import { Zmod } from '../finite_rings/integer_mod_ring.js';
import { PolynomialRing } from './polynomial_ring.js';

// Permanent comparisons include native result coefficients and the full alias matrix.
test('finite XGCD preserves its zero shortcut and nonmonic operand', () => {
  const F = new PrimeField(7n),
    R = new PolynomialRing(F),
    S = new PolynomialRing(F);
  const f = S.__call__([2n, 2n]);
  const [g, s, t] = R.zero().xgcd(f);
  expect(g).toBe(f);
  expect(s.isZero()).toBe(true);
  expect(t.coeffs.map(String)).toEqual(['1']);
});
test('QQ XGCD returns three distinct zeros and NTL returns a unit first cofactor', () => {
  const R = new PolynomialRing(QQ),
    [g, s, t] = R.zero().xgcd(R.zero());
  expect([g.isZero(), s.isZero(), t.isZero()]).toEqual([true, true, true]);
  expect(g).not.toBe(s);
  expect(g).not.toBe(t);
  expect(s).not.toBe(t);
  const N = new PolynomialRing(new PrimeField((1n << 127n) - 1n));
  expect(N.zero().xgcd(N.zero()).map(String)).toEqual(['0', '1', '0']);
});
test('composite XGCD catches native inversion failure', () => {
  const R = new PolynomialRing(Zmod(14n));
  expect(() => R.one().xgcd(R.__call__([2n]))).toThrow(
    'non-invertible elements encountered during XGCD'
  );
});
