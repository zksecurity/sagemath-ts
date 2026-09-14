import { expect, test } from 'bun:test';
import { Integer } from '../integer_ring.js';
import { QQ } from '../rational_field.js';
import { PrimeField } from '../finite_rings/finite_field_extension.js';
import { Zmod } from '../finite_rings/integer_mod_ring.js';
import { PolynomialRing } from './polynomial_ring.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';
const integers: CoefficientRing<Integer & RingElement> = {
  zero: () => new Integer(0n) as Integer & RingElement,
  one: () => new Integer(1n) as Integer & RingElement,
  __call__: (x: unknown) => new Integer(x as bigint) as Integer & RingElement,
  is_field: () => false,
  toString: () => 'Integer Ring',
} as CoefficientRing<Integer & RingElement>;

test('monic promotes integer coefficients exactly to QQ, including a leading minus one', () => {
  const R = new PolynomialRing(integers, 'x');
  for (const [coeffs, expected] of [
    [
      [2n, 4n],
      ['1/2', '1'],
    ],
    [
      [3n, 2n, 6n],
      ['1/2', '1/3', '1'],
    ],
    [
      [1n, -1n],
      ['-1', '1'],
    ],
  ] as const) {
    const f = R.__call__(coeffs.slice());
    const g = f.monic();
    expect(Object.is(g.parent.base_ring, QQ)).toBe(true);
    expect(g.coeffs.map(String)).toEqual([...expected]);
    expect(f._monic().coeffs.map(String)).toEqual([...expected]);
  }
});
test('monic integer identity uses the coefficient ring one', () => {
  const R = new PolynomialRing(integers, 'x');
  const f = R.__call__([2n, 1n]);
  expect(f.is_monic()).toBe(true);
  expect(f.monic()).toBe(f);
});
test('zero normalization preserves the distinct Sage backend errors', () => {
  expect(() => new PolynomialRing(integers).zero().monic()).toThrow('rational division by zero');
  expect(() => new PolynomialRing(QQ).zero().monic()).toThrow('rational division by zero');
  expect(() => new PolynomialRing(new PrimeField(7n)).zero().monic()).toThrow(
    'leading coefficient must be invertible'
  );
  expect(() => new PolynomialRing(new PrimeField(2n)).zero().monic()).toThrow(
    'inverse of Mod(0, 2) does not exist'
  );
});
test('FLINT normalization allocates an output and rejects nonunit leading coefficients', () => {
  const R = new PolynomialRing(Zmod(14n));
  const f = R.__call__([2n, 1n]);
  expect(f.monic()).not.toBe(f);
  expect(f.monic().eq(f)).toBe(true);
  expect(() => R.__call__([1n, 2n]).monic()).toThrow('leading coefficient must be invertible');
});
