import { expect, spyOn, test } from 'bun:test';
import { FpXQ_autpow, FpXQ_autpowers } from './galconj.js';
import * as generic from './FpX.js';
import * as word from './Flx.js';

test('automorphism caches precede native trace/unit and zero-input shortcuts', () => {
  const T = [1n, ...new Array<bigint>(35).fill(0n), 2n];
  const error = 'impossible inverse in Fp_inv: Mod(2, 4).';
  expect(() => generic.FpXQ_auttrace([[1n], [1n]], 1n, T, 4n)).toThrow(error);
  expect(() => FpXQ_autpowers([], 0, T, 4n)).toThrow(error);
  expect(() => FpXQ_autpowers([], 3, T, 4n)).toThrow(error);
  expect(() => FpXQ_autpow([], 3, T, 4n)).toThrow(error);
  expect(FpXQ_autpow([], 0, T, 4n)).toEqual([0n, 1n]);
  expect(FpXQ_autpow([], 1, T, 4n)).toEqual([]);
});

test('raw automorphism copies canonicalize storage without reducing coefficients', () => {
  const T = [3n, 0n, 1n, 0n];
  expect(FpXQ_autpowers([18n, 0n], 1, T, 17n)).toEqual([[], [0n, 1n], [18n]]);
  expect(
    generic.FpXQ_auttrace(
      [
        [18n, 0n],
        [19n, 0n],
      ],
      1n,
      T,
      17n
    )
  ).toEqual([[18n], [19n]]);
});

test('word composition skips its lower cache cutoff when coefficient reduction yields zero', () => {
  const p = (1n << 64n) - 2n;
  const T = [3n, ...new Array<bigint>(29).fill(0n), 2n];
  expect(FpXQ_autpow([p, 0n], 2, T, p)).toEqual([]);
  expect(() => FpXQ_autpow([1n], 2, T, p)).toThrow(
    'impossible inverse in Fl_inv: Mod(2, 18446744073709551614).'
  );
});

test('automorphism composition retains the prepared generic reciprocal through word calls', () => {
  const T = [3n, ...new Array<bigint>(35).fill(0n), 1n];
  const cache = spyOn(generic, 'FpX_invBarrett');
  const wordCache = spyOn(word, 'Flx_invBarrett');
  try {
    expect(FpXQ_autpowers([0n, 1n], 9, T, 17n)).toEqual([
      [],
      ...Array.from({ length: 10 }, () => [0n, 1n]),
    ]);
    expect(cache.mock.calls.length).toBe(1);
    expect(wordCache.mock.calls.length).toBe(0);
  } finally {
    cache.mockRestore();
    wordCache.mockRestore();
  }
});
