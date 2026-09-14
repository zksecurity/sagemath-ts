import { expect, test } from 'bun:test';
import { Integer } from './integer_ring.js';
import { Rational } from './rational.js';
import { RDF } from './real_double.js';

// Original Sage comparisons include IEEE bit patterns, formatting and identity.
test('RDF coefficient conversion uses Sage Integer overflow and rational conversion', () => {
  expect(RDF.__call__(new Integer(1n << 1200n)).value).toBe(Infinity);
  expect(String(RDF.__call__(-(1n << 1200n)))).toBe('-infinity');
  expect(RDF.__call__(new Rational((1n << 1200n) + 1n, 1n << 1200n)).value).toBe(1);
  expect(RDF.__call__((1n << 54n) - 1n).value).toBe(2 ** 54);
  const x = RDF.__call__(1 / 3);
  expect(RDF.__call__(x)).toBe(x);
  expect(RDF.zero()).toBe(RDF.zero());
});

test('RDF arithmetic retains IEEE exceptional values and Python double formatting', () => {
  expect(String(RDF.__call__(-0))).toBe('-0.0');
  expect(String(RDF.__call__(1e-5))).toBe('1e-05');
  expect(String(RDF.__call__(1e16))).toBe('1e+16');
  expect(String(RDF.one().div(RDF.zero()))).toBe('+infinity');
  expect(String(RDF.zero().div(RDF.zero()))).toBe('NaN');
  expect(RDF.__call__(Number.NaN).eq(RDF.__call__(Number.NaN))).toBe(false);
});
