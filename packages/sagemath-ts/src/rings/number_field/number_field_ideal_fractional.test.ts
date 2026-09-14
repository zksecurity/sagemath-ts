import { expect, test } from 'bun:test';
import { NumberField, QuadraticField, RationalPolynomial } from './number_field.js';
import { Rational } from '../rational.js';

test('fractional intersections do not truncate rational generators', () => {
  const K = QuadraticField.create(-2n);
  const half = K.ideal(new Rational(1n, 2n));
  expect(half.intersection(half).eq(half)).toBe(true);
  expect(half.intersection(K.ideal(new Rational(2n, 3n))).eq(K.ideal(2n))).toBe(true);
  expect(K.ideal(0n).intersection(half).is_zero()).toBe(true);
});

test('numerator and denominator are coprime integral ideals', () => {
  const K = QuadraticField.create(-1n, 'i');
  const I = K.ideal(K.gen().mul(4n).add(3n).div(5n));
  const N = I.numerator(),
    D = I.denominator();
  expect(N.norm().toString()).toBe('5');
  expect(D.norm().toString()).toBe('5');
  expect(N.is_integral() && D.is_integral()).toBe(true);
  expect(N.add(D).eq(K.ideal(1n))).toBe(true);
  expect(N.div(D).eq(I)).toBe(true);
  expect(I.numerator()).toBe(N);
  expect(I.denominator()).toBe(D);
});

test('integral numerators are cached fresh ideals, not the input object', () => {
  const K = QuadraticField.create(5n);
  const I = K.ideal(K.gen().add(1n).div(2n));
  expect(I.denominator().eq(K.ideal(1n))).toBe(true);
  expect(I.numerator().eq(I)).toBe(true);
  expect(I.numerator()).not.toBe(I);
  expect(I.numerator()).toBe(I.numerator());
});

test('fractional coprimality tests numerator and denominator supports', () => {
  const K = QuadraticField.create(-1n, 'i');
  const P = K.ideal(K.gen().add(2n)),
    Q = K.ideal(K.gen().neg().add(2n));
  expect(P.inverse().is_coprime(Q.pow(3n))).toBe(true);
  expect(P.inverse().is_coprime(P)).toBe(false);
  expect(P.inverse().is_coprime(P.inverse())).toBe(false);
  expect(P.is_coprime(K.ideal(0n))).toBe(false);
  expect(K.ideal(1n).is_coprime(K.ideal(0n))).toBe(true);
});

test('zero ideals lack fractional-only methods in Sage', () => {
  const I = QuadraticField.create(2n).ideal(0n);
  for (const method of ['denominator', 'numerator', 'divides', 'is_coprime'] as const) {
    expect(() =>
      method === 'divides' || method === 'is_coprime' ? I[method](I) : I[method]()
    ).toThrow(`'NumberFieldIdeal' object has no attribute '${method}'`);
  }
});

test('division preserves principal, PARI and monoid zero paths', () => {
  const K = QuadraticField.create(2n);
  const zero = K.ideal(0n),
    one = K.ideal(1n);
  expect(() => one.div(zero)).toThrow('number field element division by zero');
  expect(() => K.ideal(1n, 1n).div(zero)).toThrow('impossible inverse in ginv: 0');
  expect(() => zero.div(one)).toThrow('[;] has unsupported PARI type t_MAT');
  expect(() => zero.div(zero)).toThrow("bad operand type for unary ~: 'NumberFieldIdeal'");
  const L = new NumberField(RationalPolynomial.fromBigInts([-3n, 1n]), 'a');
  expect(L.ideal(0n).div(L.ideal(1n)).is_zero()).toBe(true);
});

test('zero products preserve the generator-count dispatch', () => {
  const K = QuadraticField.create(2n);
  expect(K.ideal(0n).mul(K.ideal(1n)).is_zero()).toBe(true);
  expect(() => K.ideal(0n).mul(K.ideal(1n, 1n))).toThrow('[;] has unsupported PARI type t_MAT');
});

test('empty and all-zero generator lists normalize to a fresh zero ideal', () => {
  const K = QuadraticField.create(5n);
  for (const I of [K.ideal(), K.ideal(0n), K.ideal(0n, 0n, 0n)]) {
    expect(I.ngens()).toBe(1);
    expect(I.gen(0).is_zero()).toBe(true);
    expect(I).not.toBe(K.ideal(0n));
  }
});
