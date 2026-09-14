import { expect, test } from 'bun:test';
import { multiple } from '../groups/generic.js';
import { Integer } from './integer_ring.js';
import { Rational } from './rational.js';
import { RDF } from './real_double.js';

test('RDF arithmetic coerces integer and rational scalars', () => {
  const a = RDF.__call__(1.5);
  expect(a.add(2n).value).toBe(3.5);
  expect(a.sub(new Integer(2n)).value).toBe(-0.5);
  expect(a.mul(new Rational(1n, 3n)).value).toBe(0.5);
  expect(a.div(2n).value).toBe(0.75);
  expect(a.eq(new Rational(3n, 2n))).toBe(true);
});

test('RDF arithmetic accepts Python-bool equivalents', () => {
  const a = RDF.__call__(2);
  expect(a.add(true).value).toBe(3);
  expect(a.sub(false).value).toBe(2);
  expect(a.mul(false).value).toBe(0);
  expect(a.div(true).value).toBe(2);
  expect(RDF.one().eq(true)).toBe(true);
});

test('RDF scalar division preserves IEEE exceptional results', () => {
  expect(RDF.one().div(0n).value).toBe(Infinity);
  expect(RDF.__call__(-1).div(false).value).toBe(-Infinity);
  expect(RDF.zero().div(new Rational(0n)).value).toBeNaN();
});

test('RDF reuses self only for native integer zero addition and one multiplication', () => {
  const a = RDF.__call__(-0);
  expect(a.add(0n)).toBe(a);
  expect(a.add(false)).toBe(a);
  expect(a.mul(1n)).toBe(a);
  expect(a.mul(true)).toBe(a);
  expect(a.add(new Integer(0n))).not.toBe(a);
  expect(a.mul(new Rational(1n))).not.toBe(a);
  expect(a.sub(0n)).not.toBe(a);
  expect(a.div(1n)).not.toBe(a);
});

test('generic multiple uses the additive identity for zero repetitions', () => {
  for (const x of [-0, -1, Infinity, NaN]) {
    expect(Object.is(multiple(RDF.__call__(x), 0n, '+').value, 0)).toBe(true);
  }
});

test('generic multiple preserves the native idempotence shortcut', () => {
  const a = RDF.__call__(-0);
  expect(multiple(a, 2n, '*')).toBe(a);
  expect(Object.is(multiple(a, 2n, '*').value, -0)).toBe(true);
});

test('generic negative multiple uses inversion before binary multiplication', () => {
  expect(multiple(RDF.zero(), -10000n, '*').value).toBe(Infinity);
  expect(Object.is(multiple(RDF.__call__(-Infinity), -2n, '*').value, -0)).toBe(true);
});

test('generic additive multiple follows the native binary rounding schedule', () => {
  const result = multiple(RDF.__call__(1e-200), 31n, '+');
  const view = new DataView(new ArrayBuffer(8));
  view.setFloat64(0, result.value);
  expect(view.getBigUint64(0)).toBe(0x16b7ba9d84a47feen);
});
