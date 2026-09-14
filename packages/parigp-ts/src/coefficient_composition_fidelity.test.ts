import { expect, spyOn, test } from 'bun:test';
import * as generic from './FpX.js';
import * as binary from './F2x.js';
import { FpXY_FpXQ_evalx, FpXY_FpXQV_evalx } from './FpXX.js';
import { FlxY_Flxq_evalx } from './FlxX.js';
import { coefficientSubstitution } from './_coefficient_composition.js';

test('native coefficient substitution preserves scalar and polynomial tags', () => {
  expect(FpXY_FpXQ_evalx([18n, [18n, 0n]], [1n], [3n, 0n, 1n], 17n)).toEqual([18n, [1n]]);
  expect(FlxY_Flxq_evalx([[1n], [1n, 1n]], [2n], [3n, 0n, 1n], 17n)).toEqual([[1n], [3n]]);
  expect(binary.F2xY_F2xq_evalx([3n, 5n], 2n, 7n)).toEqual([3n, 2n]);
  expect(FpXY_FpXQV_evalx([-19n, 0n, 23n], [], [], 0n)).toEqual([-19n, 0n, 23n]);
});

test('raw coefficient tables preserve native per-polynomial reciprocal preparation', () => {
  const T = [3n, ...Array<bigint>(35).fill(0n), 1n];
  const inverse = generic.FpX_invBarrett(T, 17n);
  const cache = spyOn(generic, 'FpX_invBarrett');
  try {
    expect(FpXY_FpXQV_evalx([[1n], [1n]], [[1n]], T, 17n)).toEqual([[1n], [1n]]);
    expect(cache).toHaveBeenCalledTimes(2);
    cache.mockClear();
    expect(coefficientSubstitution(0, 1, 17n, T, [[1n], [1n]], [], [[1n]], inverse)).toEqual([
      [1n],
      [1n],
    ]);
    expect(cache).toHaveBeenCalledTimes(0);
    expect(FpXY_FpXQV_evalx([18n, 19n], [], T, 17n)).toEqual([18n, 19n]);
    expect(cache).toHaveBeenCalledTimes(0);
  } finally {
    cache.mockRestore();
  }
});

test('empty direct substitution still builds powers before skipping coefficients', () => {
  const T = [1n, ...Array<bigint>(35).fill(0n), 2n];
  expect(() => FpXY_FpXQ_evalx([], [1n], T, 4n)).toThrow(
    'impossible inverse in Fp_inv: Mod(2, 4).'
  );
  expect(FpXY_FpXQV_evalx([], [], T, 4n)).toEqual([]);
});

test('binary scalar composition uses a square-root power table and giant steps', () => {
  const multiply = spyOn(binary, 'F2xq_mul');
  const square = spyOn(binary, 'F2xq_sqr');
  try {
    expect(binary.F2x_F2xq_eval((1n << 4096n) - 1n, 2n, 7n)).toBe(1n);
    expect(multiply.mock.calls.length + square.mock.calls.length).toBeLessThanOrEqual(128);
    multiply.mockClear();
    square.mockClear();
    expect(binary.F2x_F2xqV_eval(0n, [], 0n)).toBe(0n);
    expect(multiply).toHaveBeenCalledTimes(0);
    expect(square).toHaveBeenCalledTimes(0);
  } finally {
    multiply.mockRestore();
    square.mockRestore();
  }
});
