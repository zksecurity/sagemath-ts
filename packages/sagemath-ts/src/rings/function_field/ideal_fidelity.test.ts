import { expect, test } from 'bun:test';
import { QQ } from '../rational_field.js';
import { FunctionField, FunctionFieldIdeal_rational } from './index.js';
import type { ConstantField, ConstantFieldElement } from './constant_field.js';

const K = FunctionField(QQ as unknown as ConstantField<ConstantFieldElement>, 'x');

test('ideal parent, monoid identity and exponent one preserve native objects', () => {
  const O = K.maximal_order(), I = new FunctionFieldIdeal_rational(O, K.gen().mul(K.__call__(2n)));
  expect(I.parent()).toBe(O.ideal_monoid());
  expect(I.pow(1n)).toBe(I);
  expect(String(I.pow(1n))).toContain('2*x');
  expect(I.pow(0n)).toBe(I.parent().one());
  expect(I.parent().one()).toBe(I.parent().one());
  expect(String(I.pow(3n).gens()[0])).toBe('x^3');
});

test('zero-ideal membership preserves native division errors for zero and nonzero elements', () => {
  for (const O of [K.maximal_order(), K.maximal_order_infinite()]) {
    const I = O.ideal(K.zero());
    for (const f of [K.zero(), K.one(), 0n, 1n])
      expect(() => I.contains(f)).toThrow('fraction field element division by zero');
    expect(O.ideal(K.one()).contains(1n)).toBe(true);
  }
});
