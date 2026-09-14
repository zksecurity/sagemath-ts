import { expect, spyOn, test } from 'bun:test';
import { FpXQXn_mul, FpXQXn_sqr, FpXQX_dotproduct, random_FpXQX } from './FpXX.js';
import { FlxqXn_mul, random_FlxqX } from './FlxX.js';
import { FlxqX_dotproduct } from './Flx.js';
import { pari_init_rand, setrand, getrand } from './random.js';
import * as word from './Flx.js';

test('native generic truncation precedes reduction while word truncation follows it', () => {
  const T = [1n, 0n, 2n],
    a = [[1n, 1n]];
  expect(FpXQXn_mul(a, a, 0, T, 4n)).toEqual([]);
  expect(() => FlxqXn_mul(a, a, 0, T, 4n)).toThrow('impossible inverse in Fl_inv: Mod(2, 4).');
});

test('generic scalar truncated squares retain native integer coefficients', () => {
  const a = [18n, 19n],
    T = [1n, 0n, 1n];
  expect(FpXQXn_sqr(a, 2, T, 17n)).toEqual([324n, 684n]);
  expect(FpXQXn_mul(a, a, 2, T, 17n)).toEqual([1n, 4n]);
});

test('dot products reduce after cancellation, preserving the coefficient tag', () => {
  const T = [1n, 0n, 2n],
    a = [
      [1n, 1n],
      [1n, 1n],
    ],
    b = [
      [1n, 1n],
      [3n, 3n],
    ];
  expect(FpXQX_dotproduct(a, b, T, 4n)).toEqual([]);
  const remainder = spyOn(word, 'Flx_rem');
  try {
    expect(FlxqX_dotproduct(a, b, T, 4n)).toEqual([]);
    expect(remainder).toHaveBeenCalledTimes(1);
  } finally {
    remainder.mockRestore();
  }
  expect(FpXQX_dotproduct([1n, 1n], [1n, -1n], T, 4n)).toBe(0n);
});

test('extension random projections preserve coefficient order and unused parameters', () => {
  const saved = getrand();
  try {
    pari_init_rand();
    setrand(1n);
    expect(random_FpXQX(2, [1n, 0n, 1n], 17n)).toEqual([
      [13n, 11n],
      [11n, 3n],
    ]);
    const state = getrand();
    pari_init_rand();
    setrand(1n);
    expect(random_FlxqX(2, [1n, 0n, 1n], 17n)).toEqual([
      [13n, 11n],
      [11n, 3n],
    ]);
    expect(getrand()).toBe(state);
    expect(random_FlxqX(3, [1n], 1n)).toEqual([]);
    expect(random_FpXQX(0, [], 0n)).toEqual([]);
    expect(getrand()).toBe(state);
  } finally {
    setrand(saved);
  }
});
