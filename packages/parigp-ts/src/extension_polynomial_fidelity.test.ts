import { expect, spyOn, test } from 'bun:test';
import * as integer from './ZX.js';
import * as word from './Flx.js';
import * as binary from './F2x.js';
import { FpXQX_mul, FpXQX_sqr, FpXQX_red } from './FpXX.js';
import { FlxqX_mul, FlxqX_sqr, FlxqX_normalize } from './FlxX.js';
import { FpXQX_normalize } from './polarit3.js';

test('generic extension squaring preserves native integer coefficient tags', () => {
  const square = spyOn(integer, 'ZX_sqr');
  try {
    expect(FpXQX_sqr([4n, 4n], [2n, 0n, 1n], 5n)).toEqual([16n, 32n, 16n]);
    expect(square).toHaveBeenCalledTimes(1);
    expect(FpXQX_sqr([[4n], [4n]], [2n, 0n, 1n], 5n)).toEqual([[1n], [2n], [1n]]);
    expect(square).toHaveBeenCalledTimes(2);
  } finally {
    square.mockRestore();
  }
});

test('extension products use native integer and word Kronecker products', () => {
  const intMul = spyOn(integer, 'ZX_mul'),
    wordMul = spyOn(word, 'Flx_mul');
  try {
    const a = [[1n, 1n], [1n]],
      T = [2n, 0n, 1n];
    expect(FpXQX_mul(a, a, T, 5n)).toEqual([[4n, 2n], [2n, 2n], [1n]]);
    expect(intMul.mock.calls.length).toBeGreaterThan(0);
    expect(FlxqX_mul(a, a, T, 5n)).toEqual([[4n, 2n], [2n, 2n], [1n]]);
    expect(wordMul.mock.calls.length).toBeGreaterThan(0);
    expect(FlxqX_sqr(a, T, 5n)).toEqual([[4n, 2n], [2n, 2n], [1n]]);
  } finally {
    intMul.mockRestore();
    wordMul.mockRestore();
  }
});

test('generic monic normalization preserves raw lower coefficients while word normalization reduces them', () => {
  const T = [2n, 0n, 1n],
    a = [[1n, 0n, 1n], [1n]];
  const inverse = spyOn(word, '_Flx_extgcd');
  try {
    expect(FpXQX_normalize(a, T, 5n)).toEqual([[1n, 0n, 1n], 1n]);
    expect(FlxqX_normalize(a, T, 5n)).toEqual([[4n], [1n]]);
    expect(inverse).toHaveBeenCalledTimes(1);
    expect(FpXQX_red([5n, [5n], 1n, [1n]], T, 5n)).toEqual([0n, [], 1n, [1n]]);
  } finally {
    inverse.mockRestore();
  }
});

test('binary extension squaring uses coefficient squares and inserts zero odd coefficients', () => {
  const square = spyOn(binary, 'F2xq_sqr');
  try {
    expect(binary.F2xqX_sqr([3n, 1n], 7n)).toEqual([2n, 0n, 1n]);
    expect(square).toHaveBeenCalledTimes(2);
    expect(binary.F2xqX_normalize([7n, 1n], 7n)).toEqual([0n, 1n]);
    expect(binary.F2xqX_mul([3n, 1n], [3n, 1n], 7n)).toEqual([2n, 0n, 1n]);
  } finally {
    square.mockRestore();
  }
});

test('native extension inverse errors keep the original unreduced polynomial', () => {
  expect(() => FpXQX_normalize([[-1n, 0n, -1n]], [1n, 0n, 1n], 5n)).toThrow(
    'impossible inverse in FpXQ_inv: -y^2 - 1.'
  );
  expect(() => FlxqX_normalize([[1n, 0n, 1n]], [1n, 0n, 1n], 5n)).toThrow(
    'impossible inverse in Flxq_inv: y^2 + 1.'
  );
  expect(() => binary.F2xqX_normalize([7n], 7n)).toThrow(
    'impossible inverse in F2xq_inv: y^2 + y + 1.'
  );
});
