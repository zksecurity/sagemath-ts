import { expect, spyOn, test } from 'bun:test';
import * as kernels from './_polynomial_gcd.js';
import * as word from './Flx.js';
import * as field from './FpX.js';
import { FpX_gcd } from './ffinit.js';
import { FpX_extgcd, FpXQ_inv, ZpX_liftfact, bezout_lift_fact } from './galconj.js';
const sample = (n: number, p: bigint, salt: bigint): bigint[] => {
  let state = salt;
  const mask = (1n << 128n) - 1n;
  return Array.from({ length: n }, () => {
    state = (state * 6364136223846793005n + 1442695040888963407n) & mask;
    return state % p;
  });
};
test('PARI gcd and Bezout preserve their native scale, and callers invert it explicitly', () => {
  expect(FpX_gcd([1n, 2n, 3n], [1n, 1n], 5n)).toEqual([2n]);
  expect(FpX_extgcd([1n, 2n, 3n], [1n, 1n], 5n)).toEqual([[2n], [1n], [1n, 2n]]);
  expect(FpXQ_inv([2n], [5n, 0n, 1n], 3n)).toEqual([2n]);
  expect(ZpX_liftfact([9n, 8n, 1n], [[], [2n, 1n], [0n, 1n]], 3n, 3)).toEqual([
    [],
    [17n, 1n],
    [18n, 1n],
  ]);
  expect(bezout_lift_fact([9n, 8n, 1n], [[], [2n, 1n], [0n, 1n]], 3n, 3)).toEqual([
    [],
    [18n, 1n],
    [10n, 26n],
  ]);
});
test('PARI half-GCD recurses at each native threshold instead of running full-size Euclid', () => {
  const base = spyOn(kernels, 'polynomialHalfGcdBasecase');
  try {
    for (const [p, cutoff, isWord] of [
      [17n, 120, true],
      [(1n << 64n) - 59n, 31, true],
      [(1n << 64n) + 13n, 45, false],
    ] as const) {
      const run = isWord ? word.Flx_halfgcd : field.FpX_halfgcd;
      base.mockClear();
      run(sample(cutoff - 1, p, 17n), sample(cutoff - 2, p, 43n), p);
      expect(base).toHaveBeenCalled();
      expect(base.mock.calls[0]![0].length).toBe(cutoff - 1);
      base.mockClear();
      run(sample(2 * cutoff, p, 17n), sample(2 * cutoff - 1, p, 43n), p);
      expect(base).toHaveBeenCalled();
      expect(base.mock.calls.every((call) => call[0].length < cutoff)).toBe(true);
    }
  } finally {
    base.mockRestore();
  }
});
test('PARI large gcds select half-GCD at their distinct native word and field cutoffs', () => {
  const half = spyOn(kernels, 'polynomialHalfGcd');
  try {
    for (const [p, cutoff] of [
      [17n, 426],
      [(1n << 64n) - 59n, 937],
      [(1n << 64n) + 13n, 194],
    ] as const) {
      half.mockClear();
      FpX_gcd(sample(cutoff, p, 17n), sample(cutoff - 1, p, 43n), p);
      expect(half).not.toHaveBeenCalled();
      FpX_gcd(sample(cutoff + 1, p, 17n), sample(cutoff, p, 43n), p);
      expect(half).toHaveBeenCalled();
    }
  } finally {
    half.mockRestore();
  }
});

test('PARI quotient powers preserve the native initial base and signed-word dispatch', async () => {
  const { FpXQ_pow } = await import('./ffinit.js');
  const { FpX_split_part, FpX_nbroots } = await import('./galconj.js');
  const ws = spyOn(word, 'Flx_sqr'),
    fs = spyOn(field, 'FpX_sqr');
  try {
    FpXQ_pow([-1n, 2n, -3n], 2n, [1n, 0n, 1n], 17n);
    expect(fs).not.toHaveBeenCalled();
    expect(ws.mock.calls[0]![0]).toEqual([16n, 2n, 14n]);
    const p = (1n << 63n) + 29n;
    FpXQ_pow([-1n, 2n, -3n], 2n, [1n, 0n, 1n], p);
    expect(fs).toHaveBeenCalledTimes(1);
    expect(fs.mock.calls[0]![0]).toEqual([-1n, 2n, -3n]);
    expect(FpXQ_pow([1n], -1n, [], p)).toEqual([1n]);
    expect(FpX_split_part([17n], 17n)).toEqual([17n]);
    expect(FpX_nbroots([17n], 17n)).toBe(0);
    expect(FpX_split_part([1n, 17n], 17n)).toEqual([1n, 17n]);
  } finally {
    ws.mockRestore();
    fs.mockRestore();
  }
});
