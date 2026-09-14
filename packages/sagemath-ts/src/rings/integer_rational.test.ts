import { expect, test } from 'bun:test';
import {
  type GroupOps,
  bsgs,
  has_order,
  multiple,
  order_from_multiple,
  parseGroupOps,
} from '../groups/generic.js';
import { Integer } from './integer_ring.js';
import { Rational } from './rational.js';

test('integer-left arithmetic coerces rational operands into QQ', () => {
  const a = new Integer(2n),
    b = new Rational(1n, 3n);
  expect(String(a.add(b))).toBe('7/3');
  expect(String(a.sub(b))).toBe('5/3');
  expect(String(a.mul(b))).toBe('2/3');
  expect(a.mul(new Rational(3n))).toBeInstanceOf(Rational);
  const integerProduct: Integer = a.mul(3n);
  expect(integerProduct.value).toBe(6n);
});

test('integer-rational equality is symmetric', () => {
  const a = new Integer(2n),
    b = new Rational(2n);
  expect(a.eq(b)).toBe(true);
  expect(b.eq(a)).toBe(true);
  expect(a.eq(new Rational(7n, 3n))).toBe(false);
});

test('integer-rational ordering is exact beyond binary64 precision', () => {
  const n = 2n ** 256n,
    a = new Integer(n);
  const below = new Rational(3n * n - 1n, 3n),
    above = new Rational(3n * n + 1n, 3n);
  expect(a.gt(below)).toBe(true);
  expect(a.ge(above)).toBe(false);
  expect(a.lt(above)).toBe(true);
  expect(a.le(below)).toBe(false);
});

test('integer BSGS resolves native inversion and rational intermediate values', () => {
  expect(bsgs(new Integer(2n), new Integer(8n), [0n, 10n], '*')).toBe(3n);
  expect(bsgs(new Integer(2n), new Integer(2n ** 31n), [0n, 64n], '*')).toBe(31n);
});

test('integer group identities and exact finite orders work', () => {
  expect(has_order(new Integer(-1n), 2n, '*')).toBe(true);
  expect(has_order(new Integer(-1n), 4n, '*')).toBe(false);
  expect(order_from_multiple(new Integer(-1n), 12n, undefined, '*')).toBe(2n);
});

test('standard group parsing obtains identities for integer wrappers', () => {
  const ops: GroupOps<Integer | Rational> = parseGroupOps(
    '*',
    undefined,
    undefined,
    undefined,
    new Integer(2n)
  );
  expect(String(ops.identity)).toBe('1');
  expect(String(ops.inverse(new Integer(2n)))).toBe('1/2');
  expect(ops.isIdentity(new Rational(1n))).toBe(true);
});

test('standard group parsing accepts rational wrappers without stored parents', () => {
  const a = new Rational(2n, 3n),
    ops = parseGroupOps('+', undefined, undefined, undefined, a);
  expect(String(ops.identity)).toBe('0');
  expect(String(ops.power(a, -3n))).toBe('-2');
});

test('integer multiplicative multiple exposes rational reciprocals', () => {
  const r: Integer | Rational = multiple(new Integer(2n), -3n, '*');
  expect(r).toBeInstanceOf(Rational);
  expect(String(r)).toBe('1/8');
});
