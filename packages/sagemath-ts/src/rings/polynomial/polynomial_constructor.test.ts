/** Constructor regressions paired with tests/property/{python,typescript}/areas/polynomials. */
import { expect, test } from 'bun:test';
import { GFpn, PrimeField } from '../finite_rings/finite_field_extension.js';
import { Zmod } from '../finite_rings/integer_mod_ring.js';
import { QQ } from '../rational_field.js';
import { PolynomialRing } from './polynomial_ring.js';

test('coerce coefficients before normalization and retain the target field', () => {
  const F = new PrimeField(7n);
  const R = new PolynomialRing(F, 'x');
  const input = [Zmod(14n).__call__(9n), 3n, 7n];
  const f = R.__call__(input);
  expect(f.coeffs.map(String)).toEqual(['2', '3']);
  expect(f.coeffs.every((c) => c.parent === F)).toBe(true);
  expect(input.map(String)).toEqual(['9', '3', '7']);
  expect(R.__call__(f)).toBe(f);
});

test('QQ singleton lists and GF2X recursive coefficient lists use backend rules', () => {
  const R = new PolynomialRing(QQ, 'x');
  expect(R.__call__([[2n, 1n]]).toString()).toBe('x + 2');
  expect(() => R.__call__([[2n, 1n], 1n])).toThrow('unable to convert [2, 1] to a rational');
  const S = new PolynomialRing(new PrimeField(2n), 'x');
  expect(
    S.__call__([
      [1n, 1n],
      [1n, 1n],
    ]).toString()
  ).toBe('x^2 + 1');
});

test('embed the exact polynomial base parent as a constant', () => {
  const R = new PolynomialRing(QQ, 'y');
  const S = new PolynomialRing(R, 'x');
  const f = R.__call__([2n, 1n]);
  const g = S.__call__(f);
  expect(g.degree()).toBe(0);
  expect(g.coeffs[0]).toBe(f);
  expect(S.__call__()).toEqual(S.zero());
});

test('validate generator indices before returning x', () => {
  const R = new PolynomialRing(QQ, 'x');
  expect(R.gen(0n).toString()).toBe('x');
  expect(() => R.gen(1n)).toThrow('generator n not defined');
  expect(() => R.gen([])).toThrow("unhashable type: 'list'");
});

test('coefficient parents provide rational lifts and compatible extension conversions', () => {
  expect(QQ.__call__(new PrimeField(7n).__call__(-2n)).toString()).toBe('5');
  const F = GFpn(7n, 2);
  expect(F.__call__(Zmod(14n).__call__(9n)).toString()).toBe('2');
  expect(() => new PrimeField(7n).__call__(F.gen())).toThrow('is not in the image');
});
