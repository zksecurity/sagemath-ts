import { expect, test } from 'bun:test';
import { QQ } from '../rational_field.js';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
import { FunctionField, PlaceSet, FunctionFieldValuationRing } from './index.js';
import type { ConstantField, ConstantFieldElement } from './constant_field.js';

const K = FunctionField(QQ as unknown as ConstantField<ConstantFieldElement>, 'x');

test('places and valuation rings preserve native parent and equal-place caching', () => {
  const S = K.place_set(), I = K.maximal_order().ideal(K.gen());
  const P = S.__call__(I), Q = S.__call__(I), V = P.valuation_ring();
  expect(new PlaceSet(K)).toBe(S);
  expect(S.__call__(P)).toBe(P);
  expect(P).not.toBe(Q);
  expect(Q.valuation_ring()).toBe(V);
  expect(new FunctionFieldValuationRing(K, Q)).toBe(V);
  expect(P.residue_field()).toBe(V.residue_field());
  expect(P.residue_field('a')).toBe(V.residue_field('a'));
  expect(P.residue_field('a')).not.toBe(V.residue_field());
});

test('valuation coercion retains an existing element and public residue maps check their domain', () => {
  const P = K.maximal_order().ideal(K.gen()).place(), V = P.valuation_ring();
  const f = K.gen().add(K.one());
  expect(K.__call__(f)).toBe(f);
  expect(V.__call__(f)).toBe(f);
  expect(() => P.residue_field()[2](K.gen().inv())).toThrow(
    "1/x fails to convert into the map's domain Valuation ring at Place (x), but a `pushforward` method is not properly implemented"
  );
  expect(() => P._residue_field()[2](K.gen().inv())).toThrow('not in the valuation ring');
});

test('QQ residue reduction retains the native polynomial-numerator call', () => {
  const P = K.maximal_order().ideal(K.gen()).place();
  expect(String(P.residue_field()[2](K.__call__(2n).div(K.gen().mul(K.__call__(2n)).add(K.__call__(3n)))))).toBe('1/3');
});

test('QQ polynomial numerator clears coefficients into a stable ZZ parent', () => {
  const R = new PolynomialRing(QQ, 't'), f = R.__call__([QQ.__call__(1n).div(QQ.__call__(3n)), QQ.__call__(1n).div(QQ.__call__(2n))]);
  expect(f.denominator()).toBe(6n);
  expect(String(f.numerator())).toBe('3*t + 2');
  expect(String(f.numerator().parent.base_ring)).toBe('Integer Ring');
  expect(f.numerator().parent).toBe(f.numerator().parent);
  expect(String(f.lcm(R.__call__([1n, 1n])))).toBe('t^2 + 5/3*t + 2/3');
});
