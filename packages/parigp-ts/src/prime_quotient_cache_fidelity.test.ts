import { expect, spyOn, test } from 'bun:test';
import { FpXQ_pow } from './ffinit.js';
import { FpXQ_inv, FpXQ_powBig, FpX_Frobenius } from './galconj.js';
import { PariError } from './errors.js';
import * as generic from './FpX.js';
import * as word from './Flx.js';

function error(run: () => unknown): Error {
  try {
    run();
  } catch (e) {
    return e as Error;
  }
  throw new Error('expected native inverse error');
}

test('signed quotient power prepares reciprocals after its zero and unit shortcuts', () => {
  const T = [1n, ...Array<bigint>(90).fill(0n), 2n];
  expect(FpXQ_pow([1n], 0n, T, 4n)).toEqual([1n]);
  expect(FpXQ_pow([1n], 1n, T, 4n)).toEqual([1n]);
  for (const power of [FpXQ_pow, FpXQ_powBig]) {
    const e = error(() => power([1n], 2n, T, 4n));
    expect(e).toBeInstanceOf(PariError);
    expect(e.message).toBe('impossible inverse in Fl_inv: Mod(2, 4).');
  }
  expect(error(() => FpX_Frobenius(T, 4n)).message).toBe(
    'impossible inverse in Fl_inv: Mod(2, 4).'
  );
});

test('generic quotient powers retain the signed-word dispatch boundary', () => {
  const p = (1n << 63n) + 2n,
    T = [1n, ...Array<bigint>(35).fill(0n), 2n];
  const e = error(() => FpXQ_pow([1n], 2n, T, p));
  expect(e).toBeInstanceOf(PariError);
  expect(e.message).toBe('impossible inverse in Fp_inv: Mod(2, 9223372036854775810).');
});

test('inverse errors preserve the native backend, input and safe-inverse phase', () => {
  const T = [1n, 0n, 1n];
  expect(error(() => FpXQ_pow([], -1n, T, 17n)).message).toBe('impossible inverse in FpXQ_inv: 0.');
  expect(error(() => FpXQ_pow([], -2n, T, 17n)).message).toBe('impossible inverse in Flxq_inv: 0.');
  expect(error(() => FpXQ_inv([], T, 17n))).toBeInstanceOf(PariError);
  expect(error(() => FpXQ_inv([], T, 17n, undefined)).message).toBe(
    'impossible inverse in FpXQ_inv: 0.'
  );
  expect(error(() => FpXQ_inv([], [2n], 4n)).message).toBe('impossible inverse in FpXQ_inv: 0.');
  expect(error(() => FpXQ_pow([], -2n, [2n], 4n)).message).toBe(
    'impossible inverse in Fl_inv: Mod(2, 4).'
  );
});

test('the Hensel overload initializes through the full unsigned-word backend', () => {
  const T = [1n, 0n, 1n];
  for (const p of [17n, (1n << 64n) - 59n, (1n << 64n) + 13n]) {
    const backend = p < 1n << 64n ? 'Flxq_inv' : 'FpXQ_inv';
    expect(error(() => FpXQ_inv([], T, p, p)).message).toBe(`impossible inverse in ${backend}: 0.`);
  }
});

test('powering reuses one native reciprocal across the exponent window', () => {
  const g = spyOn(generic, 'FpX_invBarrett'),
    w = spyOn(word, 'Flx_invBarrett');
  try {
    expect(
      FpXQ_pow([1n], 513n, [1n, ...Array<bigint>(35).fill(0n), 1n], (1n << 64n) + 13n)
    ).toEqual([1n]);
    expect(g).toHaveBeenCalledTimes(1);
    expect(FpXQ_pow([1n], 513n, [1n, ...Array<bigint>(90).fill(0n), 1n], 17n)).toEqual([1n]);
    expect(w).toHaveBeenCalledTimes(1);
  } finally {
    g.mockRestore();
    w.mockRestore();
  }
});

test('native inverse display uses the original polynomial and its large-payload cutoff', () => {
  const a = [17n * 10n ** 2000n],
    T = [1n, 0n, 1n];
  expect(error(() => FpXQ_inv(a, T, 17n)).message).toBe(
    'impossible inverse in FpXQ_inv: \n  ***  (...) Huge t_POL omitted; you can access it via dbg_err().'
  );
  expect(error(() => FpXQ_pow(a, -2n, T, 17n)).message).toBe('impossible inverse in Flxq_inv: 0.');
});
