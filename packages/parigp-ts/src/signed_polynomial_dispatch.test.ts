import { expect, spyOn, test } from 'bun:test';
import * as P from './FpX.js';
import * as word from './Flx.js';
import * as integer from './ZX.js';
import { PariError } from './errors.js';

test('signed word division converts coefficients before quotient construction', () => {
  const a = [1n, 2n],
    b = [0n, 1n];
  expect(P.FpX_divrem(a, b, -2n)).toEqual([[], [1n]]);
  expect(a).toEqual([1n, 2n]);
  expect(b).toEqual([0n, 1n]);
  expect(P.FpX_extgcd([], [1n], -2n)).toEqual([[1n], [], [1n]]);
});

test('signed backend conversion preserves half-GCD storage and transformations', () => {
  expect(P.FpX_halfgcd_all([], [2n], -2n)).toEqual([
    [
      [[], [1n]],
      [[1n], []],
    ],
    [],
    [],
  ]);
  expect(P.FpX_halfgcd([2n], [], -2n)).toEqual([
    [[], [1n]],
    [[1n], []],
  ]);
});

test('native constant shortcuts retain signed generic inverse errors', () => {
  expect(P.FpX_rem([1n], [2n], -6n)).toEqual([]);
  expect(() => P.FpX_divrem([1n], [2n], -6n)).toThrow(
    new PariError('impossible inverse in Fp_inv: Mod(2, -6).')
  );
  expect(() => P.FpX_gcd([0n, 1n], [1n, 2n], -6n)).toThrow(
    new PariError('impossible inverse in Fl_inv: Mod(2, 6).')
  );
});

test('reciprocal division dispatch retains distinct native zero-divisor errors', () => {
  const a = Array<bigint>(110).fill(1n);
  expect(() => P.FpX_divrem(a, [], -2n)).toThrow(
    new PariError('impossible inverse in Flx_divrem: Vecsmall([0]).')
  );
  expect(() => P.FpX_rem(a, [], -2n)).toThrow(
    new PariError('impossible inverse in Fl_inv: Mod(0, 2).')
  );
});

test('multiplication and squaring use word magnitude without moving the cutoff', () => {
  const mul = spyOn(word, 'Flx_mul'),
    sqr = spyOn(word, 'Flx_sqr'),
    zx = spyOn(integer, 'ZX_mul');
  try {
    expect(P.FpX_mul([1n, 1n], [1n, 1n], -17n)).toEqual([1n, 2n, 1n]);
    expect(mul.mock.calls[0]![2]).toBe(17n);
    expect(P.FpX_sqr([1n, 1n], -17n)).toEqual([1n, 2n, 1n]);
    expect(sqr.mock.calls[0]![1]).toBe(17n);
    const count = mul.mock.calls.length;
    expect(P.FpX_mul([1n, 1n], [1n, 1n], -((1n << 64n) + 13n))).toEqual([1n, 2n, 1n]);
    expect(mul.mock.calls.length).toBe(count);
    expect(zx).toHaveBeenCalledTimes(1);
  } finally {
    mul.mockRestore();
    sqr.mockRestore();
    zx.mockRestore();
  }
});
