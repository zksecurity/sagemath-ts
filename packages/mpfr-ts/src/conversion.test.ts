import { expect, test } from 'bun:test';
import { mpfr_init2, mpfr_set_str, mpfr_get_d, mpfr_get_str } from './index.js';

test('MPFR conversion retains digits beyond binary64 and rounds ties to even', () => {
  const x = mpfr_init2(128);
  expect(mpfr_set_str(x, '9007199254740993')).toBe(0);
  expect(mpfr_get_d(x)).toBe(9007199254740992);
  expect(mpfr_get_str(10, 16, x)).toEqual(['9007199254740993', 16]);
  expect(mpfr_set_str(x, '9007199254740995')).toBe(0);
  expect(mpfr_get_d(x)).toBe(9007199254740996);
});

test('MPFR mutable parsing preserves empty input but resets NaN sign', () => {
  const x = mpfr_init2(53);
  mpfr_set_str(x, '-1');
  expect(mpfr_set_str(x, '')).toBe(-1);
  expect(mpfr_get_d(x)).toBe(-1);
  expect(mpfr_set_str(x, '-nan')).toBe(0);
  expect(x.sign).toBe(1);
  expect(Number.isNaN(mpfr_get_d(x))).toBe(true);
  expect(mpfr_set_str(x, '1.25junk')).toBe(-1);
  expect(mpfr_get_d(x)).toBe(1.25);
});

test('unsupported native conversion domains fail explicitly', () => {
  const x = mpfr_init2(53);
  expect(() => mpfr_init2(0)).toThrow('invalid MPFR precision');
  expect(() => mpfr_init2(4097)).toThrow('precision above 4096 bits');
  expect(() => mpfr_set_str(x, '1', 2)).toThrow('non-decimal');
  expect(() => mpfr_set_str(x, '1', 10, 'RNDU')).toThrow();
  expect(() => mpfr_set_str(x, '1e4097')).toThrow('decimal exponent');
  expect(() => mpfr_set_str(x, '1'.repeat(4097))).toThrow('decimal mantissa');
  expect(() => mpfr_get_str(2, 0, x)).toThrow('non-decimal');
  expect(() => mpfr_get_str(10, -1, x)).toThrow('invalid MPFR digit count');
  expect(() => mpfr_get_str(10, 4097, x)).toThrow('output digits');
});
