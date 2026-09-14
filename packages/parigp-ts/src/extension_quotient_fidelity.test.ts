import { expect, spyOn, test } from 'bun:test';
import * as generic from './FpXX.js';
import * as word from './FlxX.js';
import * as binary from './F2x.js';
import * as wordPolynomial from './Flx.js';
import { extensionField } from './_extension_field.js';
import { inverseCoefficient } from './_polynomial_division.js';

test('quotient inverse and powers preserve native coefficient representations', () => {
  const T = [3n, 0n, 1n],
    S = [1n, 0n, 1n];
  expect(generic.FpXQXQ_inv([1n, 1n], S, T, 17n)).toEqual([9n, 8n]);
  expect(generic.FpXQXQ_invsafe(S, S, T, 17n)).toBeNull();
  expect(generic.FpXQXQ_powers([1n, 1n], 3, S, T, 17n)).toEqual([
    [1n],
    [1n, 1n],
    [0n, 2n],
    [15n, 2n],
  ]);
  expect(word.FlxqXQ_powu([[1n], [1n]], 3n, [[1n], [], [1n]], T, 17n)).toEqual([[15n], [2n]]);
  expect(binary.F2xqXQ_inv([2n, 1n], [3n, 1n, 1n], 7n)).toEqual([2n, 3n]);
  expect(binary.F2xqXQ_powers([2n, 1n], 3, [3n, 1n, 1n], 7n)).toEqual([
    [1n],
    [2n, 1n],
    [0n, 1n],
    [3n, 3n],
  ]);
});

test('word powering prepares one inner reciprocal and keeps zero/unit shortcuts ahead of it', () => {
  const T = [1n, ...new Array<bigint>(90).fill(0n), 1n];
  const S = [[1n], [], [1n]],
    a = [[0n, 1n], [1n]];
  const spy = spyOn(wordPolynomial, 'Flx_invBarrett');
  try {
    expect(word.FlxqXQ_pow(a, 0n, S, T, 17n)).toEqual([[1n]]);
    expect(word.FlxqXQ_pow(a, 1n, S, T, 17n)).toEqual(a);
    expect(spy.mock.calls.length).toBe(0);
    word.FlxqXQ_pow(a, 13n, S, T, 17n);
    expect(spy.mock.calls.length).toBe(1);
    expect(spy.mock.calls[0]![0]).toEqual(T);
  } finally {
    spy.mockRestore();
  }
});

test('generic word powering retains the native cached scalar-modulus conversion', () => {
  const T = [3n, 0n, 1n],
    S = new Array<bigint>(11).fill(1n);
  const cache = generic.FpXQX_get_red(S, T, 17n);
  expect(Array.isArray(cache)).toBe(false);
  expect(generic.FpXQXQ_pow([1n, 1n], 2n, cache, T, 17n)).toEqual([]);
  const tagged = S.map((c) => [c]),
    taggedCache = generic.FpXQX_get_red(tagged, T, 17n);
  expect(() => generic.FpXQXQ_pow([1n, 1n], 2n, taggedCache, T, 17n)).toThrow(
    'generic word powering cannot convert cached polynomial coefficients'
  );
  expect(generic.FpXQXQ_pow([[1n]], 1n, taggedCache, T, 17n)).toEqual([[1n]]);
});

test('native inverse display clips polynomial and modular-integer payloads', () => {
  const p = (1n << 6000n) - 1n;
  expect(() => inverseCoefficient(p, p, false)).toThrow(
    'impossible inverse in Fp_inv: \n  ***  (...) Huge t_INTMOD omitted; you can access it via dbg_err().'
  );
  const T = (1n << 300n) - 1n;
  expect(() => extensionField(2, T, 2n).inv(T)).toThrow(
    'impossible inverse in F2xq_inv: \n  ***  (...) Huge t_POL omitted; you can access it via dbg_err().'
  );
});
