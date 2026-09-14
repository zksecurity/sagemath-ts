import { expect, test } from 'bun:test';
import { NumberField, QuadraticField, RationalPolynomial } from './number_field.js';
import { Rational } from '../rational.js';

test('ideal integral bases use the maximal order beyond the equation order', () => {
  const K = QuadraticField.create(5n);
  const basis = K.ideal(1n).integral_basis();
  expect(basis.length).toBe(2);
  expect(basis.some((a) => a.list().some((c) => c.denominator === 2n))).toBe(true);
});

test('fractional and redundant cubic ideals have a full HNF basis', () => {
  const K = new NumberField(RationalPolynomial.fromBigInts([-8n, -2n, -1n, 1n]), 'a');
  const I = K.ideal(
    new Rational(2n, 3n),
    K.gen().add(1n).div(3n),
    K.gen().add(1n).mul(new Rational(2n, 3n))
  );
  const basis = I.integral_basis();
  expect(basis.length).toBe(3);
  expect(basis.every((a) => I.contains(a))).toBe(true);
  expect(I.integral_basis()).not.toBe(basis);
  expect(I.free_module()).toBe(I.free_module());
  expect(I.free_module().rank).toBe(3);
});

test('zero ideals have an empty basis and rank zero', () => {
  const I = QuadraticField.create(5n).ideal(0n);
  expect(I.integral_basis()).toEqual([]);
  expect(I.zk_basis()).toEqual([]);
  expect(I.free_module()).toEqual({ basis: [], rank: 0 });
  expect(I.free_module()).toBe(I.free_module());
});

test('two generators regenerate ideals even when the first supplied nonrational generator does not', () => {
  const K = QuadraticField.create(-5n);
  const I = K.ideal(9n, K.gen().mul(3n), K.gen().add(2n));
  const pair = I.gens_two();
  expect(pair[0].toString()).toBe('3');
  expect(K.ideal(...pair).eq(I)).toBe(true);
  expect(I.gens_two()).toBe(pair);
  expect(pair.every((a) => a.parent() === K)).toBe(true);
});

test('two-generator intersections with QQ retain fractional denominators', () => {
  const K = QuadraticField.create(-5n);
  const I = K.ideal(new Rational(3n, 2n));
  const [a, b] = I.gens_two();
  expect(a.toString()).toBe('3/2');
  expect(b.is_zero()).toBe(true);
  expect(a.parent()).toBe(K);
  expect(b.parent()).toBe(K);
  expect(K.ideal(a, b).eq(I)).toBe(true);
});

test('the second generator vanishes exactly for rational ideals', () => {
  const K = QuadraticField.create(-5n);
  for (const I of [K.ideal(0n), K.ideal(12n), K.ideal(12n, 24n)]) {
    expect(I.gens_two()[1].is_zero()).toBe(true);
  }
  expect(K.ideal(K.gen().add(2n)).gens_two()[1].is_zero()).toBe(false);
});
