import { Integer, ZZ } from './integer_ring.js';
import { Zmod } from './finite_rings/integer_mod_ring.js';
import { GF2 } from './finite_rings/gf2.js';
import { expect, test } from 'bun:test';
import { PolynomialRing } from './polynomial/polynomial_ring.js';
import type { CoefficientRing, RingElement } from './polynomial/polynomial_element.js';
import { QQ } from './rational_field.js';
import { Rational } from './rational.js';
import { PrimeField } from './finite_rings/finite_field_extension.js';
import {
  FractionField,
  FractionField_generic,
  FractionField_1poly_field,
} from './fraction_field.js';
import {
  FractionFieldElement,
  FractionFieldElement_1poly_field,
} from './fraction_field_element.js';
import { FpT, FpTElement } from './fraction_field_FpT.js';
import { ZeroDivisionError } from '../errors.js';

test('fraction factories choose cached native and generic parents', () => {
  const R = new PolynomialRing(QQ, 'x'),
    K = R.fraction_field();
  expect(FractionField(R)).toBe(K);
  expect(K).toBeInstanceOf(FractionField_1poly_field);
  const f = K.__call__([1n, 2n, 1n], [1n, 1n]);
  expect(f).toBeInstanceOf(FractionFieldElement_1poly_field);
  expect(String(f)).toBe('x + 1');
  expect(String(f.inv())).toBe('1/(x + 1)');
  expect(String(f.mul(f.inv()))).toBe('1');
  expect(K.__call__(f)).toBe(f);
  expect(K.zero()).toBe(K.zero());
  const generic = new FractionField_generic(R),
    raw = new FractionFieldElement(generic, R.__call__([1n, 2n, 1n]), R.__call__([1n, 1n]), {
      reduce: false,
    });
  expect(String(raw)).toBe('(x^2 + 2*x + 1)/(x + 1)');
  raw.reduce();
  expect(String(raw)).toBe('x + 1');
});
test('native finite fractions preserve unreduced representation and safe zero inversion', () => {
  const R = new PolynomialRing(new PrimeField(7n), 'x'),
    K = R.fraction_field();
  expect(K).toBeInstanceOf(FpT);
  const raw = new FpTElement(K, R.zero(), R.gen(), { reduce: false });
  expect(String(raw)).toBe('0/x');
  expect(raw.isZero()).toBe(false);
  expect(raw.eq(K.zero())).toBe(false);
  expect(raw.numerator()).not.toBe(raw.numerator());
  expect(raw.denominator()).not.toBe(raw.denominator());
  expect(raw.numer().eq(raw.numerator())).toBe(true);
  expect(raw.denom().eq(raw.denominator())).toBe(true);
  expect(() => raw.inv()).toThrow(ZeroDivisionError);
  expect(() => raw.pow(-2n)).toThrow(ZeroDivisionError);
  const f = K.__call__([1n, 1n], [6n, 1n]);
  expect(String(f.mul(f.inv()))).toBe('1');
  expect('reduce' in f).toBe(false);
});

test('fraction sections preserve generic identity and normalize native storage', () => {
  const R = new PolynomialRing(QQ, 'x'),
    K = R.fraction_field();
  const f = K.__call__([1n, 2n, 1n], [1n, 1n]);
  expect(R.__call__(f)).toBe(f.numerator());
  const raw = new FractionFieldElement(new FractionField_generic(R), R.gen(), R.gen(), {
    reduce: false,
  });
  expect(() => R.__call__(raw)).toThrow('fraction must have unit denominator');
  const T = new PolynomialRing(new PrimeField(7n), 'x'),
    F = T.fraction_field();
  const native = new FpTElement(F, T.gen(), T.gen(), { reduce: false });
  expect(String(native)).toBe('x/x');
  expect(String(T.__call__(native))).toBe('1');
  expect(String(native)).toBe('1');
});

test('nested fraction coefficients retain inverses and their canonical equality parent', () => {
  const U = new PolynomialRing(QQ, 't'),
    F = U.fraction_field();
  const R = new PolynomialRing(F, 'x'),
    t = F.__call__(U.gen());
  expect(R.__call__(t).eq(t)).toBe(true);
  expect(U.gen().eq(t)).toBe(true);
  expect(R.gen().eq(t)).toBe(false);
  expect(F.__call__(R.__call__(t))).toBe(t);
  expect(String(U.__call__(R.__call__(t)))).toBe('t');
  expect(() => F.__call__(R.gen())).toThrow('not a constant polynomial');
  expect(
    R.one()
      .pseudo_quo_rem(R.__call__([1n, 0n, U.gen()]))
      .map(String)
  ).toEqual(['0', '1/t']);
});

test('scalar sections differ from direct conversion hooks and preserve original errors', () => {
  const R = new PolynomialRing(QQ, 'x'),
    x = R.gen();
  expect(() => QQ.__call__(x)).toThrow('x is not a constant polynomial');
  expect(() => x._rational_()).toThrow('cannot convert nonconstant polynomial');
  expect(String(new PrimeField(7n).__call__(R.__call__(new Rational(3n, 2n))))).toBe('5');
  expect(() => Zmod(14n).__call__(R.__call__(new Rational(3n, 2n)))).toThrow(
    'inverse of Mod(2, 14) does not exist'
  );
  expect(GF2.__call__(R.__call__(3n))._integer_()).toBe(1n);
  expect(String(GF2.__call__(1n)._rational_())).toBe('1');
});

test('generic scalar conversion reduces in place and uses the ring unit inverse', () => {
  const integers = {
    zero: () => new Integer(0n),
    one: () => new Integer(1n),
    __call__: (x: unknown) => new Integer(x as bigint),
    is_field: () => false,
    toString: () => 'Integer Ring',
  };
  const R = new PolynomialRing(integers as unknown as CoefficientRing<RingElement>, 'x'),
    K = new FractionField_generic(R);
  const raw = new FractionFieldElement(K, R.zero(), R.gen(), { reduce: false });
  const zero = raw.numerator();
  expect(String(QQ.__call__(raw))).toBe('0');
  expect(raw.numerator()).not.toBe(zero);
  expect(String(raw.denominator())).toBe('1');
  const half = new FractionFieldElement(K, 1n, 2n);
  expect(String(QQ.__call__(half))).toBe('1/2');
  expect(() => ZZ.__call__(half)).toThrow('inverse does not exist');
  expect(() => half._conversion(Zmod(14n))).toThrow('element is not a unit');
});

test('native scalar sections normalize only on the native conversion path', () => {
  const R = new PolynomialRing(new PrimeField(7n), 'x'),
    K = R.fraction_field();
  const raw = new FpTElement(K, R.gen(), R.gen(), { reduce: false });
  expect(() => QQ.__call__(raw as unknown as bigint)).toThrow(
    'unable to convert x/x to a rational'
  );
  expect(String(raw)).toBe('x/x');
  expect(ZZ.__call__(raw as unknown as bigint)).toBe(1n);
  expect(String(raw)).toBe('1');
});
