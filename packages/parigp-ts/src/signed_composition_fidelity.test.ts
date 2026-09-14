import { expect, spyOn, test } from 'bun:test';
import * as P from './FpX.js';
import { FpXQ_autpow } from './galconj.js';
import { FpM_mul } from './FpV.js';
import * as word from './Flx.js';
import * as contexts from './_polynomial_quotient.js';
import * as binary from './F2v.js';
import { PariError } from './errors.js';

test('signed power tables retain the native count-dependent initial conversion', () => {
  const x = [-1n],
    T = [1n, 1n];
  expect(P.FpXQ_powers(x, 2, T, -2n)).toEqual([[1n], [-1n], [1n]]);
  expect(P.FpXQ_powers(x, 3, T, -2n)).toEqual([[1n], [1n], [1n], [1n]]);
  expect(x).toEqual([-1n]);
  expect(T).toEqual([1n, 1n]);
});

test('direct and table composition preserve signed modular matrix arithmetic', () => {
  expect(P.FpX_FpXQ_eval([1n, 2n, 3n], [1n], [1n, 1n], -17n)).toEqual([6n]);
  expect(P.FpX_FpXQV_eval([1n], [[1n]], [1n, 1n], -((1n << 64n) + 13n))).toEqual([1n]);
  expect(P.FpX_FpXQ_eval([], [], [], 0n)).toEqual([]);
});

test('matrix shortcuts and zero products precede zero-modulus reduction', () => {
  expect(FpM_mul([[]], [[]], 0n)).toEqual([[]]);
  expect(FpM_mul([[]], [[], [0n], [0n]], -1n)).toEqual([[], [0n], [0n]]);
  const a = [[], [0n, 1n], [0n, 1n]],
    b = [[], [0n, 1n, -1n]];
  expect(FpM_mul(a, b, 0n)).toEqual([[], [0n, 0n]]);
  expect(() => FpM_mul([[], [0n, 1n]], [[], [0n, 1n]], 0n)).toThrow(
    new PariError('impossible inverse in dvmdii: 0.')
  );
});

test('negative characteristic two uses the binary matrix backend', () => {
  const mul = spyOn(binary, 'F2m_mul');
  const a = [[], [0n, 1n, 3n], [0n, 2n, 4n]],
    saved = structuredClone(a);
  try {
    expect(FpM_mul(a, a, -2n)).toEqual([[], [0n, 1n, 1n], [0n, 0n, 0n]]);
    expect(mul).toHaveBeenCalledTimes(1);
    expect(a).toEqual(saved);
  } finally {
    mul.mockRestore();
  }
});

test('signed trace power tables convert and reuse the generic reciprocal', () => {
  const T = [1n, ...Array<bigint>(38).fill(0n), 1n];
  const pair: [bigint[], bigint[]] = [Array<bigint>(15).fill(1n), Array<bigint>(12).fill(2n)];
  const expected = P.FpXQ_auttrace(pair, 2n, T, 17n);
  const generic = spyOn(P, 'FpX_invBarrett'),
    flx = spyOn(word, 'Flx_invBarrett');
  try {
    expect(P.FpXQ_auttrace(pair, 2n, T, -17n)).toEqual(expected);
    expect(generic).toHaveBeenCalledTimes(1);
    expect(flx).toHaveBeenCalledTimes(0);
  } finally {
    generic.mockRestore();
    flx.mockRestore();
  }
});

test('automorphism direct composition follows the signed word dispatch', () => {
  const ctx = spyOn(contexts, 'polynomialQuotient');
  try {
    expect(FpXQ_autpow([1n, 1n], 2, [1n, 0n, 1n], -17n)).toEqual([2n, 1n]);
    expect(ctx.mock.calls.some((args) => args[1] === 17n && args[2] === true)).toBe(true);
  } finally {
    ctx.mockRestore();
  }
});
