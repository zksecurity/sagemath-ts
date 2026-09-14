import { expect, spyOn, test } from 'bun:test';
import * as word from './Flx.js';
import * as field from './FpX.js';
import * as native from './FpX_factor.js';
import { FpX_mul } from './ffinit.js';
import { FpX_is_totally_split, FpX_roots } from './galconj.js';
import { getrand, setrand } from './random.js';

test('native roots reduce before degree checks and preserve zero errors', () => {
  expect(() => FpX_roots([], 7n)).toThrow('zero polynomial in FpX_roots.');
  expect(() => native.Flx_roots([], 7n)).toThrow('zero polynomial in Flx_roots.');
  expect(FpX_roots([1n, 7n], 7n)).toEqual([]);
  expect(() => FpX_roots([0n, 7n], 7n)).toThrow('zero polynomial in FpX_roots.');
});

test('total splitting retains raw-degree and unsigned-word zero boundaries', () => {
  expect(FpX_is_totally_split([7n], 7n)).toBe(true);
  expect(FpX_is_totally_split([0n, 7n], 7n)).toBe(true);
  expect(FpX_is_totally_split([], 7n)).toBe(false);
  expect(() => FpX_is_totally_split([], (1n << 64n) + 13n)).toThrow(
    'impossible inverse in FpX_divrem: 0.'
  );
  expect(FpX_is_totally_split([1n, 2n, 1n], 7n)).toBe(false);
});

test('general word roots use native signed sorting while quadratic roots use unsigned sorting', () => {
  const p = (1n << 64n) - 59n;
  const q = FpX_mul([p - 1n, 1n], [1n, 1n], p);
  const f = FpX_mul(q, [p - 2n, 1n], p);
  expect(native.Flx_roots(q, p)).toEqual([1n, p - 1n]);
  expect(FpX_roots(f, p)).toEqual([p - 1n, 1n, 2n]);
});

test('word root counts use quadratic and binary power shortcuts', () => {
  const gcd = spyOn(word, 'Flx_gcd');
  const inverse = spyOn(word, 'Flx_invBarrett');
  try {
    expect(native.Flx_nbroots([6n, 0n, 1n], 7n)).toBe(2);
    expect(gcd).toHaveBeenCalledTimes(0);
    const f = [1n, ...new Array<bigint>(127).fill(0n), 1n];
    expect(native.Flx_nbroots(f, 2n)).toBe(1);
    expect(gcd).toHaveBeenCalledTimes(1);
    expect(inverse).toHaveBeenCalledTimes(0);
  } finally {
    gcd.mockRestore();
    inverse.mockRestore();
  }
});

test('word root splitting delegates to Flx and preserves global random state', () => {
  const saved = getrand();
  const gcd = spyOn(word, 'Flx_gcd');
  const fieldGcd = spyOn(field, 'FpX_gcd');
  try {
    setrand(42n);
    const seeded = getrand();
    expect(FpX_roots([16n, 0n, 0n, 1n], 17n)).toEqual([1n]);
    expect(gcd.mock.calls.length).toBeGreaterThan(0);
    expect(fieldGcd).toHaveBeenCalledTimes(0);
    expect(getrand()).toBe(seeded);
  } finally {
    gcd.mockRestore();
    fieldGcd.mockRestore();
    setrand(saved);
  }
});
