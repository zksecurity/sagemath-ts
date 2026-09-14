import { expect, spyOn, test } from 'bun:test';
import * as word from './Flx.js';
import * as field from './FpX.js';
import * as factors from './FpX_factor.js';
import * as minimal from './_polynomial_minpoly.js';
import { _galconj_factor_squarefree_irreducibles } from './galconj.js';
import { getrand, setrand } from './random.js';

test('native normalization keeps forced leading ones and inverse payloads', () => {
  expect(field.FpX_normalize([1n, 2n], 1n)).toEqual([0n, 1n]);
  expect(field.FpX_normalize([0n, 0n], 7n)).toEqual([]);
  expect(field.FpX_normalize([8n, 1n], 7n)).toEqual([8n, 1n]);
  expect(() => field.FpX_normalize([1n, 6n], 9n)).toThrow(
    'impossible inverse in Fp_inv: Mod(3, 9).'
  );
  expect(() => word.Flx_normalize([1n, 6n], 9n)).toThrow(
    'impossible inverse in Fl_inv: Mod(6, 9).'
  );
  expect(() => factors.Flx_factor([], 7n)).toThrow('impossible inverse in Fl_inv: Mod(0, 7).');
  expect(factors.FpX_factor([], 7n)).toEqual([[[], 1]]);
});

test('word squarefree components use native word GCD and division', () => {
  const wordGcd = spyOn(word, 'Flx_gcd');
  const fieldGcd = spyOn(field, 'FpX_gcd');
  const wordDivision = spyOn(word, 'Flx_divrem');
  try {
    expect(factors.FpX_factor_squarefree([1n, 3n, 3n, 1n], 5n)).toEqual([[1n], [1n], [1n, 1n]]);
    expect(wordGcd.mock.calls.length).toBeGreaterThan(0);
    expect(wordDivision.mock.calls.length).toBeGreaterThan(0);
    expect(fieldGcd).toHaveBeenCalledTimes(0);
  } finally {
    wordGcd.mockRestore();
    fieldGcd.mockRestore();
    wordDivision.mockRestore();
  }
});

test('the squarefree Galois helper delegates to native full factorization', () => {
  const state = getrand();
  const factor = spyOn(factors, 'FpX_factor');
  try {
    setrand(1n);
    expect(_galconj_factor_squarefree_irreducibles([0n, 2n, 2n, 1n], 5n)).toEqual([
      [0n, 1n],
      [3n, 1n],
      [4n, 1n],
    ]);
    expect(factor).toHaveBeenCalledTimes(1);
  } finally {
    factor.mockRestore();
    setrand(state);
  }
});

test('Shoup chooses the minimal-polynomial splitter for two quadratic factors', () => {
  const state = getrand();
  const minpoly = spyOn(minimal, 'polynomialMinimalPolynomial');
  try {
    setrand(1n);
    expect(factors.FpX_factor([1n, 0n, 0n, 0n, 1n], 5n)).toEqual([
      [[2n, 0n, 1n], 1],
      [[3n, 0n, 1n], 1],
    ]);
    expect(minpoly.mock.calls.length).toBeGreaterThan(0);
    expect(minpoly.mock.calls.every((call) => call[4] !== undefined)).toBe(true);
    minpoly.mockClear();
    setrand(1n);
    expect(factors.FpX_factor([0n, 2n, 2n, 1n], 5n)).toHaveLength(3);
    expect(minpoly).toHaveBeenCalledTimes(0);
  } finally {
    minpoly.mockRestore();
    setrand(state);
  }
});
