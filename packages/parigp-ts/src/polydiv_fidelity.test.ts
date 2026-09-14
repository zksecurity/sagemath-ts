import { expect, spyOn, test } from 'bun:test';
import * as kernels from './_polynomial_division.js';
import * as word from './Flx.js';
import * as field from './FpX.js';
import { FpX_divrem, FpX_rem } from './ffinit.js';

test('PARI division distinguishes remainder shortcuts, normalized zeros and inverse errors', () => {
  expect(FpX_rem([1n, 2n], [2n], 8n)).toEqual([]);
  expect(() => FpX_divrem([1n, 2n], [2n], 8n)).toThrow('impossible inverse in Fp_inv: Mod(2, 8).');
  expect(FpX_divrem([-1n], [1n, 1n], 17n)).toEqual([[], [16n]]);
  expect(() => FpX_divrem([], [0n], 17n)).toThrow('impossible inverse in FpX_divrem: 0.');
  expect(FpX_divrem([1n, -2n, -3n], [1n, 1n], (1n << 64n) + 13n)).toEqual([[1n, -3n], []]);
});

test('PARI division preserves independent word remainder and quotient cutoffs', () => {
  const inverse = spyOn(word, 'Flx_invBarrett');
  try {
    for (const [p, divcut, remcut] of [
      [17n, 161, 159],
      [3037000494n, 14, 89],
    ] as const) {
      const T = [1n, 1n, 1n];
      inverse.mockClear();
      word.Flx_divrem(new Array(divcut - 1).fill(1n), T, p);
      expect(inverse).not.toHaveBeenCalled();
      word.Flx_divrem(new Array(divcut).fill(1n), T, p);
      expect(inverse).toHaveBeenCalledTimes(1);
      inverse.mockClear();
      word.Flx_rem(new Array(remcut - 1).fill(1n), T, p);
      expect(inverse).not.toHaveBeenCalled();
      word.Flx_rem(new Array(remcut).fill(1n), T, p);
      expect(inverse).toHaveBeenCalledTimes(1);
    }
  } finally {
    inverse.mockRestore();
  }
});

test('PARI field division delegates word inputs and uses native inverse thresholds', () => {
  const wi = spyOn(word, 'Flx_divrem'),
    fi = spyOn(field, 'FpX_invBarrett'),
    inverse = spyOn(kernels, 'inverseBarrett');
  try {
    FpX_divrem([1n, 2n, 1n], [1n, 1n], 17n);
    expect(wi).toHaveBeenCalledTimes(1);
    const p = (1n << 64n) + 13n,
      T = [1n, 1n, 1n];
    FpX_divrem(new Array(112).fill(1n), T, p);
    expect(fi).not.toHaveBeenCalled();
    FpX_divrem(new Array(113).fill(1n), T, p);
    expect(fi).toHaveBeenCalledTimes(1);
    fi.mockClear();
    FpX_rem(new Array(110).fill(1n), T, p);
    expect(fi).not.toHaveBeenCalled();
    FpX_rem(new Array(111).fill(1n), T, p);
    expect(fi).toHaveBeenCalledTimes(1);
    for (const [modulus, cutoff, isWord] of [
      [17n, 200, true],
      [3037000494n, 22, true],
      [p, 110, false],
    ] as const) {
      const inv = isWord ? word.Flx_invBarrett : field.FpX_invBarrett;
      inverse.mockClear();
      inv(new Array(cutoff - 1).fill(1n), modulus);
      expect(inverse.mock.calls[0]![3]).toBe(true);
      inv(new Array(cutoff).fill(1n), modulus);
      expect(inverse.mock.calls[1]![3]).toBe(false);
    }
  } finally {
    wi.mockRestore();
    fi.mockRestore();
    inverse.mockRestore();
  }
});
