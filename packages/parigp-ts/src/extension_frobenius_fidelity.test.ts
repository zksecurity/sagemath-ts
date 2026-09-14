import { expect, spyOn, test } from 'bun:test';
import * as api from './FpXQX_factor.js';
import * as word from './Flx.js';
import * as generic from './FpX.js';
import * as automorphism from './_extension_automorphism.js';
import * as quotient from './_extension_quotient.js';
import { PariError } from './errors.js';

test('full Frobenius selects the native composition or direct-power branch', () => {
  const aut = spyOn(automorphism, 'extensionAutomorphism'),
    pow = spyOn(quotient, 'extensionQuotient');
  try {
    for (const mode of [0, 1, 2])
      for (const degree of [9, 16]) {
        aut.mockClear();
        pow.mockClear();
        const S = [1n, ...Array<bigint>(degree - 1).fill(0n), 1n];
        if (mode === 0) api.FpXQX_Frobenius(S, [1n, 0n, 1n], 3n);
        else if (mode === 1)
          api.FlxqX_Frobenius(
            S.map((c) => (c ? [c] : [])),
            [1n, 0n, 1n],
            3n
          );
        else api.F2xqX_Frobenius(S, 7n);
        const compose = mode !== 2 && degree === 9;
        expect(aut.mock.calls.length).toBe(compose ? 1 : 0);
        if (!compose)
          expect(pow.mock.calls.some((c) => c[1] === 5 && c[3] === (mode === 2 ? 4n : 9n))).toBe(
            true
          );
      }
  } finally {
    aut.mockRestore();
    pow.mockRestore();
  }
});

test('word characteristic-two Frobenius skips a fresh reciprocal and half Frobenius reuses one', () => {
  const inv = spyOn(word, 'Flx_invBarrett');
  const T = [1n, ...Array<bigint>(90).fill(0n), 1n];
  try {
    expect(word.Flx_Frobenius(T, 2n)).toEqual([0n, 0n, 1n]);
    expect(inv).toHaveBeenCalledTimes(0);
    expect(api.FlxqXQ_halfFrobenius([[1n]], [[1n], [], [1n]], T, 2n)).toEqual([[1n]]);
    expect(inv).toHaveBeenCalledTimes(1);
  } finally {
    inv.mockRestore();
  }
});

test('generic half Frobenius reuses prepared inner and outer caches', () => {
  const inv = spyOn(generic, 'FpX_invBarrett'),
    pow = spyOn(quotient, 'extensionQuotient');
  const T = [1n, ...Array<bigint>(35).fill(0n), 1n],
    S = [1n, ...Array<bigint>(9).fill(0n), 1n];
  try {
    expect(api.FpXQXQ_halfFrobenius([1n], S, T, (1n << 64n) + 13n)).toEqual([1n]);
    expect(inv).toHaveBeenCalledTimes(1);
    const initial = pow.mock.calls[0]!;
    expect(Array.isArray(initial[7])).toBe(false);
    expect(initial[8]).toBeDefined();
    expect(pow.mock.calls.every((c) => c[7] === initial[7] && c[8] === initial[8])).toBe(true);
  } finally {
    inv.mockRestore();
    pow.mockRestore();
  }
});

test('half Frobenius preserves word conversion at the unsigned boundary', () => {
  const T = [1n, ...Array<bigint>(35).fill(0n), 2n],
    S = [1n, 0n, 1n],
    p = (1n << 63n) + 2n;
  expect(() => api.FpXQXQ_halfFrobenius([1n], S, T, p)).toThrow(PariError);
  expect(() => api.FpXQXQ_halfFrobenius([1n], S, T, p)).toThrow(
    'impossible inverse in Fl_inv: Mod(2, 9223372036854775810).'
  );
  expect(() => api.FpXQX_Frobenius(S, T, p)).toThrow(
    'impossible inverse in Fp_inv: Mod(2, 9223372036854775810).'
  );
});

test('constant outer quotients and source arrays retain their native behavior', () => {
  const S = [[1n, 0n], [], [1n], []],
    T = [3n, 0n, 1n, 0n],
    saved = structuredClone([S, T]);
  expect(api.FpXQX_Frobenius(S, T, 17n)).toEqual([0n, 1n]);
  expect([S, T]).toEqual(saved);
  expect(api.FpXQX_Frobenius([1n], T, 17n)).toEqual([]);
  expect(api.FlxqX_Frobenius([[1n]], T, 17n)).toEqual([]);
  expect(api.F2xqX_Frobenius([1n], 7n)).toEqual([]);
});

test('generic extension moduli must retain positive degree after word conversion', () => {
  for (const run of [
    () => api.FpXQX_Frobenius([1n, 0n, 1n], [1n, 0n, 2n], 2n),
    () => api.FpXQXQ_halfFrobenius([1n], [1n, 0n, 1n], [1n, 0n, 2n], 2n),
  ]) {
    expect(run).toThrow(RangeError);
    expect(run).toThrow('extension modulus must have positive degree');
  }
});
