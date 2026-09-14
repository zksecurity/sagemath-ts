import { expect, test } from 'bun:test';
import { FpX_red, FpX_add, FpX_sub, FpX_neg, FpX_Fp_mul } from './ffinit.js';
import { PariError } from './errors.js';

test('prime polynomial linear operations follow native signed-modulus reduction', () => {
  expect(FpX_red([-1n], -17n)).toEqual([16n]);
  expect(FpX_add([-1n], [], -17n)).toEqual([16n]);
  expect(FpX_sub([], [1n], -17n)).toEqual([16n]);
  expect(FpX_neg([1n], -17n)).toEqual([16n]);
  expect(FpX_Fp_mul([-1n], 1n, -17n)).toEqual([16n]);
});

test('zero modulus preserves canonical empty, cancellation and raw-zero-scalar shortcuts', () => {
  expect(FpX_red([0n, 0n], 0n)).toEqual([]);
  expect(FpX_neg([0n], 0n)).toEqual([]);
  expect(FpX_add([1n], [-1n], 0n)).toEqual([]);
  expect(FpX_add([0n, 1n], [0n, -1n], 0n)).toEqual([]);
  expect(FpX_sub([0n, 1n], [0n, 1n], 0n)).toEqual([]);
  expect(FpX_Fp_mul([1n, 1n, 1n], 0n, 0n)).toEqual([]);
  expect(FpX_Fp_mul([0n, 0n], 1n, 0n)).toEqual([]);
});

test('reached zero-modulus reductions preserve the full PARI diagnostic', () => {
  for (const run of [
    () => FpX_red([1n], 0n),
    () => FpX_add([1n], [1n], 0n),
    () => FpX_sub([1n], [], 0n),
    () => FpX_neg([1n], 0n),
    () => FpX_Fp_mul([1n], 1n, 0n),
    () => FpX_add([0n, 1n], [0n], 0n),
  ]) {
    expect(run).toThrow(PariError);
    expect(run).toThrow('impossible inverse in dvmdii: 0.');
  }
});

test('canonicalization preserves input storage', () => {
  const a = [-1n, 18n, 0n],
    saved = a.slice();
  expect(FpX_red(a, -17n)).toEqual([16n, 1n]);
  expect(a).toEqual(saved);
});
