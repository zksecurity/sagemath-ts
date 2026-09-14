import { expect, test } from 'bun:test';
import { Rational } from '../rings/rational.js';
import { Integer } from '../rings/integer_ring.js';
import {
  NumberField,
  QuadraticField,
  RationalPolynomial,
} from '../rings/number_field/number_field.js';
import {
  bsgs,
  discrete_log,
  has_order,
  multiple,
  order_from_bounds,
  order_from_multiple,
  pohlig_hellman,
} from './generic.js';

test('number-field scalar action accepts exact integers and rationals', () => {
  const a = QuadraticField.create(2n, 'a').gen();
  expect(a.mul(3n).list().map(String)).toEqual(['0', '3']);
  expect(a.mul(new Integer(-3n)).list().map(String)).toEqual(['0', '-3']);
  expect(a.mul(new Rational(7n, 3n)).list().map(String)).toEqual(['0', '7/3']);
});

test('additive multiple uses the number-field scalar action', () => {
  const a = QuadraticField.create(5n, 'a').gen();
  for (const n of [-3n, 0n, 1n, 12n])
    expect(multiple(a, n, '+').list().map(String)).toEqual(['0', String(n)]);
});

test('BSGS accepts Sage identity predicates on number fields', () => {
  const a = QuadraticField.create(2n, 'a').gen();
  expect(bsgs(a, a.mul(31n), [0n, 64n], '+')).toBe(31n);
  expect(bsgs(a, a.pow(31n), [0n, 64n], '*')).toBe(31n);
});

test('number-field torsion works through both discrete-log entrypoints', () => {
  const a = QuadraticField.create(-1n, 'a').gen();
  expect(discrete_log(a.pow(3n), a, 4n, '*')).toBe(3n);
  expect(pohlig_hellman(a.pow(3n), a, 4n, undefined, '*')).toBe(3n);
});

test('number-field exact orders distinguish multiples', () => {
  const a = QuadraticField.create(-1n, 'a').gen();
  expect(has_order(a, 4n, '*')).toBe(true);
  expect(has_order(a, 12n, '*')).toBe(false);
  expect(order_from_multiple(a, 12n, undefined, '*')).toBe(4n);
  expect(order_from_bounds(a, [1n, 16n], undefined, '*')).toBe(4n);
});

test('cubic fields support scalar actions and additive BSGS', () => {
  const K = new NumberField(RationalPolynomial.fromBigInts([-2n, 0n, 0n, 1n]), 'a');
  expect(K.gen().mul(new Rational(-7n, 3n)).list().map(String)).toEqual(['0', '-7/3', '0']);
  expect(bsgs(K.gen(), K.gen().mul(64n), [0n, 64n], '+')).toBe(64n);
});
