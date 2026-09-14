import { expect, test } from 'bun:test';
import { PrimeField } from '../finite_rings/finite_field_extension.js';
import { QQ } from '../rational_field.js';
import { Rational } from '../rational.js';
import { RDF } from '../real_double.js';
import { PolynomialRing } from './polynomial_ring.js';

test('derivative supports counts and method aliases with zero-count identity', () => {
  const R = new PolynomialRing(QQ, 'x');
  const f = R.__call__([-1n, 0n, new Rational(1n, 2n), 0n, -1n]);
  expect(f.derivative(R.gen(), 2n).coeffs.map(String)).toEqual(['1', '0', '-12']);
  expect(f.diff).toBe(f.derivative);
  expect(f.differentiate).toBe(f.derivative);
  expect(f.diff([])).toBe(f);
  expect(f.differentiate(0n)).toBe(f);
  expect(f.gradient()[0]!.coeffs.map(String)).toEqual(['0', '1', '0', '-4']);
});

test('generic zero derivatives retain identity while QQ uses a fresh native result', () => {
  const R = new PolynomialRing(new PrimeField(7n), 'x'),
    f = R.zero();
  expect(f.derivative()).toBe(f);
  expect(f.gradient()[0]).toBe(f);
  const q = new PolynomialRing(QQ, 'x').zero();
  expect(q.derivative()).not.toBe(q);
  expect(() => q._derivative('x')).toThrow('cannot differentiate with respect to x');
});

test('nested derivative variables recurse into polynomial coefficients', () => {
  const R = new PolynomialRing(QQ, 't'),
    t = R.gen(),
    S = new PolynomialRing(R, 'x');
  const f = S.__call__([t.pow(2), t]);
  expect(f.derivative(t).coeffs.map(String)).toEqual(['2*t', '1']);
  expect(f.derivative(S.gen(), t).coeffs.map(String)).toEqual(['1']);
});

test('real derivative uses one scalar multiplication with native rounding', () => {
  const R = new PolynomialRing(RDF, 'x');
  const f = R.__call__([...Array(11).fill(0), 6392024.733746229]);
  expect(f.derivative().getCoeff(10).value).toBe(70312272.07120852);
});
