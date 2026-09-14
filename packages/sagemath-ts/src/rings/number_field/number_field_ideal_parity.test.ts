import { describe, expect, test } from 'bun:test';
import { Integer } from '../integer_ring.js';
import { Rational } from '../rational.js';
import { NumberField, RationalPolynomial } from './number_field.js';
import type { AbsoluteOrder } from './order.js';

const field = (coefficients: bigint[]) =>
  new NumberField(RationalPolynomial.fromBigInts(coefficients), 'a');

describe('Sage ideal inversion and different regressions', () => {
  test.each([
    [[-5n, 0n, 1n], '5'],
    [[23n, 0n, 1n], '23'],
    [[-8n, -2n, -1n, 1n], '503'],
  ] as const)('uses the maximal-order trace dual for %s', (coefficients, norm) => {
    const K = field([...coefficients]);
    const D = K.different();
    expect(D.norm().toString()).toBe(norm);
    expect(K.different()).toBe(D);
    const O = K.maximal_order() as AbsoluteOrder;
    expect(O.different()).toBe(D);
    expect(O.codifferent().norm().toString()).toBe(`1/${norm}`);
    expect(D.mul(O.codifferent()).eq(K.ideal(1n))).toBe(true);
  });

  test('inverts ideals in a shifted quadratic presentation', () => {
    const K = field([-1n, -1n, 1n]);
    const I = K.ideal(3n, K.gen().add(1n));
    expect(I.inverse().mul(I).eq(K.ideal(1n))).toBe(true);
    expect(I.inverse().inverse().eq(I)).toBe(true);
  });

  test('inverts redundant fractional generators in degree three', () => {
    const K = field([-8n, -2n, -1n, 1n]);
    const b = K.gen().add(1n);
    const I = K.ideal(new Rational(2n, 3n), b.div(3n), b.mul(new Rational(2n, 3n)));
    expect(I.inverse().mul(I).eq(K.ideal(1n))).toBe(true);
    expect(I.pow(-2n).eq(I.mul(I).inverse())).toBe(true);
  });

  test('preserves exponent-one identity and accepts Integer wrappers', () => {
    const K = field([23n, 0n, 1n]);
    const I = K.ideal(2n, K.gen().add(1n));
    expect(I.pow(1n)).toBe(I);
    expect(I.pow(new Integer(1n))).toBe(I);
    expect(I.pow(new Integer(-2n)).eq(I.inverse().mul(I.inverse()))).toBe(true);
  });

  test('bounds generator growth during large nonprincipal powers', () => {
    const K = field([23n, 0n, 1n]);
    const I = K.ideal(2n, K.gen().add(1n));
    const result = I.pow(100n);
    expect(result.ngens()).toBeLessThanOrEqual(K.degree());
    expect(result.norm().toString()).toBe((I.norm().numerator ** 100n).toString());
    expect(result.mul(I.pow(-100n)).eq(K.ideal(1n))).toBe(true);
  });

  test('matches zero-ideal power identity and both Sage exponent error paths', () => {
    const K = field([-2n, 0n, 1n]);
    const zero = K.ideal(0n);
    expect(zero.pow(0n).eq(K.ideal(1n))).toBe(true);
    expect(zero.pow(1n)).toBe(zero);
    expect(zero.pow(new Integer(2n)).is_zero()).toBe(true);
    expect(() => zero.inverse()).toThrow("bad operand type for unary ~: 'NumberFieldIdeal'");
    expect(() => zero.pow(-1n)).toThrow(
      "unsupported operand type(s) for ** or pow(): 'NumberFieldIdeal' and 'int'"
    );
    expect(() => zero.pow(new Integer(-1n))).toThrow(
      "bad operand type for unary ~: 'NumberFieldIdeal'"
    );
  });
});
