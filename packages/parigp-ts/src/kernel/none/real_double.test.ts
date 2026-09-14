import { expect, test } from 'bun:test';
import { dbltor, rtodbl } from './mp_indep.js';
import { dbltor as buchDbltor, rtodbl as buchRtodbl } from '../../buch.js';
import { PariError } from '../../errors.js';

test('binary64 conversion preserves native zero accuracy and whole-word working precision', () => {
  expect(buchDbltor(-0)).toEqual({ s: 0, e: -1023, m: 0n, p: 0 });
  expect(buchDbltor(0, 200)).toEqual({ s: 0, e: -1023, m: 0n, p: 256 });
  expect(buchDbltor(0, 1024)).toEqual({ s: 0, e: -1024, m: 0n, p: 1024 });
  expect(buchDbltor(Number.MIN_VALUE)).toEqual({ s: 1, e: -1074, m: 1n << 63n, p: 64 });
});

test('binary64 output retains native tie rounding and subnormal cutoff', () => {
  for (const convert of [rtodbl, buchRtodbl]) {
    expect(convert({ s: 1, e: 0, p: 64, m: (1n << 63n) + 1024n })).toBe(1 + 2 ** -52);
    expect(Object.is(convert({ s: -1, e: -1024, p: 64, m: 1n << 63n }), 0)).toBe(true);
    expect(Object.is(convert({ s: -1, e: -1023, p: 64, m: 1n << 63n }), -0)).toBe(true);
  }
});

test('binary64 overflow reports the native PARI exception', () => {
  for (const convert of [dbltor, buchDbltor]) {
    for (const d of [NaN, Infinity, -Infinity]) {
      expect(() => convert(d)).toThrow(new PariError('overflow in dbltor [NaN or Infinity]'));
      expect(() => convert(d)).toThrow(PariError);
    }
  }
  for (const convert of [rtodbl, buchRtodbl]) {
    expect(() => convert(dbltor(Number.MAX_VALUE))).toThrow(PariError);
    expect(() => convert({ s: 1, e: 1022, p: 64, m: (1n << 64n) - 1n }))
      .toThrow('overflow in t_REAL->double conversion');
  }
});
