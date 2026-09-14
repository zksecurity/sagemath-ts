import { expect, test } from 'bun:test';
import { Integer } from '../integer_ring.js';
import { PolynomialRing } from './polynomial_ring.js';
import { Polynomial, type CoefficientRing, type RingElement } from './polynomial_element.js';
const integers = {
  zero: () => new Integer(0n), one: () => new Integer(1n),
  __call__: (x: unknown) => new Integer(x as bigint),
  is_field: () => false, toString: () => 'Integer Ring',
} as CoefficientRing<Integer & RingElement>;
const R = new PolynomialRing(integers, 'x');

test('Sage integer XGCD returns Integer constants and retains polynomial zero shortcuts', () => {
  const result = R.__call__(6n).xgcd(R.__call__(9n));
  expect(result.map(String)).toEqual(['3', '-1', '1']);
  expect(result.every((g) => g instanceof Integer)).toBe(true);
  const f = R.__call__([2n, 2n]);
  const [g, s, t] = f.xgcd(R.zero());
  expect(g).toBe(f);
  expect(s instanceof Integer).toBe(true);
  expect(t instanceof Integer).toBe(true);
  expect([String(s), String(t)]).toEqual(['1', '0']);
});
test('Sage integer XGCD clears the rational Bezout triple for shared factors', () => {
  const result = R.__call__([2n, 2n]).xgcd(R.__call__([4n, 6n, 2n]));
  expect(result.map(String)).toEqual(['2*x + 2', '1', '0']);
  expect(result.every((g) => g instanceof Polynomial)).toBe(true);
});
