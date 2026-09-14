import { expect, test } from 'bun:test';
import { QQ } from '../rational_field.js';
import { Integer } from '../integer_ring.js';
import { Zmod } from '../finite_rings/integer_mod_ring.js';
import { FractionField_generic } from '../fraction_field.js';
import { FractionFieldElement } from '../fraction_field_element.js';
import { PolynomialRing } from './polynomial_ring.js';

test('Sage degree gaps scale the remainder by a negative power', () => {
  const R = new PolynomialRing(QQ, 'x');
  const [q, r] = R.__call__(1n).pseudo_quo_rem(R.__call__([2n, 0n, 0n, 2n]));
  expect([String(q), String(r)]).toEqual(['0', '1/4']);
  const Z = new PolynomialRing(Zmod(14n), 'x');
  expect(() => Z.__call__(1n).pseudo_quo_rem(Z.__call__([2n, 0n, 0n, 2n]))).toThrow(
    'inverse of Mod(4, 14) does not exist'
  );
});
test('original operand protocols distinguish scalars and constant polynomials', () => {
  const R = new PolynomialRing(QQ, 'x');
  const [scalar] = R.zero().pseudo_quo_rem(new Integer(2n));
  const [poly] = R.zero().pseudo_quo_rem(R.__call__(2n));
  expect(String(scalar.parent)).toBe(String(R));
  expect(poly).toBeInstanceOf(FractionFieldElement);
  expect(String(poly.parent)).toBe(`Fraction Field of ${R}`);
  expect(() => R.gen().pseudo_quo_rem(false)).toThrow("'bool' object has no attribute 'is_zero'");
  expect(() => R.gen().pseudo_quo_rem(0n)).toThrow('Pseudo-division by zero is not possible');
});
test('constant-denominator fraction data and arithmetic agree with Sage', () => {
  const R = new PolynomialRing(QQ, 'x'),
    K = new FractionField_generic(R);
  const a = K.__call__([1n, 2n], 3n),
    b = K.__call__(2n, 3n);
  expect([
    String(a.numerator()),
    String(a.denominator()),
    String(a.add(b)),
    String(a.sub(b)),
    String(a.mul(b)),
    String(a.div(b)),
    String(a.neg()),
  ]).toEqual([
    '2/3*x + 1/3',
    '1',
    '2/3*x + 1',
    '2/3*x - 1/3',
    '4/9*x + 2/9',
    'x + 1/2',
    '-2/3*x - 1/3',
  ]);
  expect(K.__call__(a)).toBe(a);
  expect(K.ring()).toBe(R);
  expect(K.is_field()).toBe(true);
  expect(K.zero().isZero()).toBe(true);
  expect(K.one().eq(1n)).toBe(true);
});
test('zero-ring scalars do not acquire a synthetic polynomial coercion', () => {
  const R = new PolynomialRing(Zmod(14n), 'x'),
    z = Zmod(1n).zero();
  expect(R.zero().eq(z)).toBe(false);
  expect(R.__call__(2n).eq(z)).toBe(false);
});
