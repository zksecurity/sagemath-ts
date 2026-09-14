import { expect, test } from 'bun:test';
import { GF } from '../finite_rings/index.js';
import { QQ } from '../rational_field.js';
import { PolynomialRing } from './polynomial_ring.js';

test('polynomial dictionary constructors preserve backend coefficient rules', () => {
  const R = new PolynomialRing(QQ, 'x');
  expect(String(R.__call__(new Map([[2n, 3n], [0n, -2n]])))).toBe('3*x^2 - 2');
  expect(String(R.__call__({ 2: 3n, 0: -2n }))).toBe('3*x^2 - 2');
  expect(() => R.__call__(new Map([[0n, [1n, 2n]]]))).toThrow('unable to convert [1, 2] to a rational');
  const S = new PolynomialRing(GF(3n) as never, 'x');
  expect(String(S.__call__(new Map([[2n, [1n, 1n]]])))).toBe('x^3 + x^2');
  expect(() => S.__call__(new Map([[[2n], 1n]]))).toThrow("'tuple' object cannot be interpreted as an integer");
});
