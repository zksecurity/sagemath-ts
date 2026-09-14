import { expect, test } from 'bun:test';
import { cmprr } from './cmp.js';
import { cmprr as buchCmprr } from '../../buch.js';
import { itor, real_0_bit } from '../../qfb.js';

test('native comparison observes the absolute accuracy of zero', () => {
  for (const compare of [cmprr, buchCmprr]) {
    expect(compare(real_0_bit(0), itor(1n, 64))).toBe(0);
    expect(compare(itor(-1n, 64), real_0_bit(0))).toBe(0);
    expect(compare(real_0_bit(-1), itor(1n, 64))).toBe(-1);
    expect(compare(itor(-1n, 64), real_0_bit(-1))).toBe(-1);
  }
});

test('native comparison never allocates according to the exponent gap', () => {
  for (const compare of [cmprr, buchCmprr]) {
    const tiny = { ...itor(1n, 64), e: -(2 ** 52) };
    const huge = { ...itor(1n, 64), e: 2 ** 52 };
    expect(compare(tiny, huge)).toBe(-1);
    expect(compare(huge, tiny)).toBe(1);
    expect(compare({ ...tiny, s: -1 }, { ...huge, s: -1 })).toBe(1);
  }
});

test('native comparison includes low words after a common mantissa prefix', () => {
  for (const compare of [cmprr, buchCmprr]) {
    const a = itor(1n, 64), b = itor(1n, 192);
    expect(compare(a, b)).toBe(0);
    expect(compare(a, { ...b, m: b.m + 1n })).toBe(-1);
    expect(compare({ ...b, m: b.m + 1n }, a)).toBe(1);
  }
});
