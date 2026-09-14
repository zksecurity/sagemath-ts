import { expect, spyOn, test } from 'bun:test';
import * as integers from './ZV.js';
import * as wordPolynomial from './Flx.js';
import { FlxqM_mul } from './alglin1.js';
import { FpXQX_FpXQXQ_eval, FpXQX_FpXQXQV_eval } from './FpXX.js';
import { FlxqX_FlxqXQ_eval, FlxqX_FlxqXQV_eval } from './FlxX.js';
import { F2xqX_F2xqXQ_eval, F2xqX_F2xqXQV_eval } from './F2x.js';

test('native extension composition preserves coefficient representations and supplied tables', () => {
  const T = [3n, 0n, 1n],
    S = [1n, 0n, 1n];
  expect(FpXQX_FpXQXQ_eval([1n, 1n, 1n], [1n, 1n], S, T, 17n)).toEqual([2n, 3n]);
  expect(FpXQX_FpXQXQV_eval([1n, 1n, 1n], [[1n], [1n, 1n]], S, T, 17n)).toEqual([2n, 3n]);
  expect(FlxqX_FlxqXQ_eval([[1n], [1n], [1n]], [[1n], [1n]], [[1n], [], [1n]], T, 17n)).toEqual([
    [2n],
    [3n],
  ]);
  expect(
    FlxqX_FlxqXQV_eval([[1n], [1n], [1n]], [[[1n]], [[1n], [1n]]], [[1n], [], [1n]], T, 17n)
  ).toEqual([[2n], [3n]]);
  expect(F2xqX_F2xqXQ_eval([1n, 2n, 1n], [2n, 1n], [3n, 1n, 1n], 7n)).toEqual([2n, 3n]);
  expect(F2xqX_F2xqXQV_eval([1n, 2n, 1n], [[1n], [2n, 1n]], [3n, 1n, 1n], 7n)).toEqual([2n, 3n]);
});

test('word extension matrices delegate packed entries to native integer-matrix schedules', () => {
  const A = [[], [[], [1n], [2n]], [[], [3n], [4n]]];
  const spy = spyOn(integers, 'ZM_mul');
  try {
    expect(FlxqM_mul(A, A, [3n, 0n, 1n], 17n)).toEqual([[], [[], [7n], [10n]], [[], [15n], [5n]]]);
    expect(spy.mock.calls.length).toBe(1);
    expect(spy.mock.calls[0]![0]).toEqual([[], [0n, 1n, 2n], [0n, 3n, 4n]]);
  } finally {
    spy.mockRestore();
  }
});

test('word composition shares the prepared inner reciprocal with its matrix backend', () => {
  const T = [1n, ...new Array<bigint>(90).fill(0n), 1n],
    S = [[1n], [], [1n]];
  const spy = spyOn(wordPolynomial, 'Flx_invBarrett');
  try {
    expect(FlxqX_FlxqXQ_eval([[1n], [1n], [1n]], [[1n], [1n]], S, T, 17n)).toEqual([[2n], [3n]]);
    expect(spy.mock.calls.length).toBe(1);
    expect(spy.mock.calls[0]![0]).toEqual(T);
  } finally {
    spy.mockRestore();
  }
});
