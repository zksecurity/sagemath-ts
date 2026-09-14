import { expect, spyOn, test } from 'bun:test';
import * as packing from './_polynomial_packing.js';
import * as word from './Flx.js';
import * as integer from './ZX.js';
import { FpX_mul, FpX_sqr, FpXQ_powers } from './FpX.js';

test('PARI packed polynomial products preserve signed carries and valuations', () => {
  const x = [0n, 0n, -1n, 1n, -1n];
  expect(integer.ZX_sqr(x)).toEqual([0n, 0n, 0n, 0n, 1n, -2n, 3n, -2n, 1n]);
  expect(integer.ZX_mul(x, [1n, -1n, 1n])).toEqual([0n, 0n, -1n, 2n, -3n, 2n, -1n]);
  expect(word.Flx_mul([0n, 2n, 2n], [2n, 2n], 4n)).toEqual([]);
  expect(FpX_sqr([-1n, 1n, -1n], 3n)).toEqual([1n, 1n, 0n, 1n, 1n]);
});

test('PARI word products and squares use their distinct native packing cutoffs', () => {
  const spy = spyOn(packing, 'packUnsigned');
  try {
    for (const [p, mulCut, squareCut] of [
      [17n, 30, 37],
      [3037000494n, 8, 14],
    ] as const) {
      spy.mockClear();
      word.Flx_mul(new Array(mulCut - 1).fill(1n), new Array(mulCut - 1).fill(1n), p);
      expect(spy).not.toHaveBeenCalled();
      word.Flx_mul(new Array(mulCut).fill(1n), new Array(mulCut).fill(1n), p);
      expect(spy).toHaveBeenCalledTimes(2);
      spy.mockClear();
      word.Flx_sqr(new Array(squareCut - 1).fill(1n), p);
      expect(spy).not.toHaveBeenCalled();
      word.Flx_sqr(new Array(squareCut).fill(1n), p);
      expect(spy).toHaveBeenCalledTimes(1);
    }
  } finally {
    spy.mockRestore();
  }
});

test('PARI field products delegate and quotient power tables call the square kernel', () => {
  const wm = spyOn(word, 'Flx_mul'),
    zm = spyOn(integer, 'ZX_mul'),
    square = spyOn(word, 'Flx_sqr');
  try {
    FpX_mul([1n, 1n], [1n, 1n], 17n);
    expect(wm).toHaveBeenCalledTimes(1);
    FpX_mul([1n, 1n], [1n, 1n], 1n << 64n);
    expect(zm).toHaveBeenCalledTimes(1);
    FpXQ_powers([1n, 1n], 4, [1n, 0n, 1n], 17n);
    expect(square).toHaveBeenCalledTimes(2);
  } finally {
    wm.mockRestore();
    zm.mockRestore();
    square.mockRestore();
  }
});

test('PARI products split oversized temporaries while every result coefficient remains representable', () => {
  const p = (1n << 64n) - 59n,
    n = 4096,
    a = new Array<bigint>(n).fill(p - 1n);
  const triangle = (i: number, length: number) => BigInt(i < length ? i + 1 : 2 * length - 1 - i);
  const product = word.Flx_mul(a, a, p),
    square = word.Flx_sqr(a, p);
  expect(product).toHaveLength(2 * n - 1);
  expect(square).toHaveLength(2 * n - 1);
  expect(product.every((c, i) => c === triangle(i, n))).toBe(true);
  expect(square.every((c, i) => c === triangle(i, n))).toBe(true);
  const b = new Array<bigint>(1024).fill(1n << 512n),
    integerSquare = integer.ZX_sqr(b);
  expect(integerSquare).toHaveLength(2047);
  expect(integerSquare.every((c, i) => c === triangle(i, 1024) << 1024n)).toBe(true);
});
