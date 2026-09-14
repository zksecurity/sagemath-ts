import { expect, spyOn, test } from 'bun:test';
import * as R from './FpXQX_factor.js';
import { FpXX_deriv } from './FpXX.js';
import { FlxX_deriv } from './FlxX.js';
import { FpX_nbroots } from './galconj.js';
import { residue, inverseCoefficient } from './_polynomial_division.js';
import * as frobenius from './_extension_frobenius.js';
import { PariError } from './errors.js';

test('extension root counts preserve prime dispatch and unscaled split parts', () => {
  const T = [1n, 0n, 1n],
    S = [1n, 0n, 1n];
  expect(R.FpXQX_split_part([2n, 0n, 2n], T, 3n)).toEqual([2n, 0n, 2n]);
  expect(R.FpXQX_nbroots(S, T, 3n)).toBe(2);
  expect(R.FlxqX_nbroots([[1n], [], [1n]], T, 3n)).toBe(2);
  expect(R.F2xqX_nbroots([1n, 1n, 1n], 7n)).toBe(2);
  expect(R.FqX_nbroots(S, null, 3n)).toBe(0);
  expect(R.FqX_nbroots(S, T, 3n)).toBe(2);
});

test('degree shortcuts preserve native conversion and unused inner moduli', () => {
  const large = (1n << 64n) + 13n;
  expect(R.FpXQX_nbroots([], [], 17n)).toBe(-1);
  expect(R.FpXQX_nbroots([17n], [], 17n)).toBe(-1);
  expect(R.FpXQX_nbroots([large], [], large)).toBe(0);
  expect(R.FpXQX_nbroots([1n, 17n], [], 17n)).toBe(0);
  expect(R.FpXQX_nbroots([1n, large], [], large)).toBe(1);
  expect(R.FlxqX_nbroots([[], [1n]], [], 3n)).toBe(1);
  expect(R.F2xqX_nbroots([1n, 7n], 0n)).toBe(1);
});

test('derivatives preserve coefficient tags and detect repeated factors', () => {
  const f = [1n, [2n, 3n], 4n],
    saved = structuredClone(f);
  expect(FpXX_deriv(f, 17n)).toEqual([[2n, 3n], 8n]);
  expect(FlxX_deriv([[1n], [2n, 3n], [4n]], 17n)).toEqual([[2n, 3n], [8n]]);
  expect(f).toEqual(saved);
  expect(R.FpXQX_is_squarefree([1n, -2n, 1n], [1n, 0n, 1n], 3n)).toBe(false);
  expect(R.FlxqX_is_squarefree([[1n], [], [1n]], [1n, 0n, 1n], 3n)).toBe(true);
});

test('zero scalar reductions precede derivative division and preserve the reached error', () => {
  expect(residue(0n, 0n)).toBe(0n);
  expect(FpXX_deriv([1n], 0n)).toEqual([]);
  expect(() => FpXX_deriv([0n, 1n], 0n)).toThrow('impossible inverse in dvmdii: 0.');
  expect(() => FpXX_deriv([0n, 0n, [1n]], 0n)).toThrow('impossible inverse in umodui: 0.');
  expect(() => FpXX_deriv([0n, 0n, [0n, 1n]], 0n)).toThrow('impossible inverse in umodui: 0.');
});

test('signed inverse residues and zero-modulus failures match native invmod', () => {
  expect(inverseCoefficient(-1n, -17n, false)).toBe(16n);
  expect(inverseCoefficient(0n, -1n, false)).toBe(0n);
  expect(() => inverseCoefficient(2n, -6n, false)).toThrow(PariError);
  expect(() => inverseCoefficient(2n, -6n, false)).toThrow(
    'impossible inverse in Fp_inv: Mod(2, -6).'
  );
  expect(() => inverseCoefficient(1n, 0n, false)).toThrow(
    'impossible inverse in Fp_inv: Mod(1, 0).'
  );
  expect(FpX_nbroots([1n, 0n, 1n], -((1n << 64n) + 13n))).toBe(0);
});

test('nonlinear root counts use Frobenius while degree shortcuts avoid it', () => {
  const f = spyOn(frobenius, 'extensionFrobenius');
  try {
    expect(R.FpXQX_nbroots([0n, 1n], [], 3n)).toBe(1);
    expect(f).toHaveBeenCalledTimes(0);
    expect(R.FpXQX_nbroots([1n, 0n, 1n], [1n, 0n, 1n], 3n)).toBe(2);
    expect(f).toHaveBeenCalledTimes(1);
    expect(f.mock.calls[0]![0]).toBe(1);
  } finally {
    f.mockRestore();
  }
});
