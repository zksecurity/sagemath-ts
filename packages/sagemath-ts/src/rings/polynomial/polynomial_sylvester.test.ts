import { expect, test } from 'bun:test';
import { Integer } from '../integer_ring.js';
import { QQ } from '../rational_field.js';
import { Rational } from '../rational.js';
import { GF2 } from '../finite_rings/gf2.js';
import { GFpn } from '../finite_rings/finite_field_extension.js';
import { PolynomialRing } from './polynomial_ring.js';

// Each behavior also runs through the permanent original-Sage comparison area.
test('Sylvester scalar operands, zero errors and optional variable conversion', () => {
  const R = new PolynomialRing(QQ, 'x');
  const f = R.__call__([1n, 0n, 1n]);
  expect(f.sylvester_matrix(new Rational(2n)).map((row) => row.map(String))).toEqual([
    ['2', '0'],
    ['0', '2'],
  ]);
  expect(() => f.sylvester_matrix(R.zero())).toThrow(
    'The Sylvester matrix is not defined for zero polynomials'
  );
  expect(() => f.sylvester_matrix(f, 'x')).toThrow("'str' object has no attribute 'parent'");
  expect(() => R.one().sylvester_matrix(new Integer(2n))).toThrow('tuple index out of range');
  expect(f.sylvester_matrix(R.one(), false).map((row) => row.map(String))).toEqual([
    ['1', '0'],
    ['0', '1'],
  ]);
});
test('explicit finite polynomial conversion differs from a QQ coefficient list', () => {
  const R = new PolynomialRing(QQ, 'x');
  const a = GFpn(7n, 2, undefined, 'a').gen();
  expect(R.__call__(a).coeffs.map(String)).toEqual(['0', '1']);
  expect(() => R.__call__([a])).toThrow('unable to convert a to a rational');
  const B = new PolynomialRing(GF2, 'x');
  expect(B.__call__([new Integer(3n), new Rational(3n)]).coeffs.map(String)).toEqual(['1', '1']);
});
test('Sylvester preserves the native small PARI-field matrix backend error', () => {
  const R = new PolynomialRing(GFpn(7n, 2, undefined, 'a'), 'x');
  expect(() => R.one().sylvester_matrix(R.one())).toThrow(
    "'FiniteField_pari_ffelt_with_category' object has no attribute '_cache'"
  );
});
