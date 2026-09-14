import { expect, spyOn, test } from 'bun:test';
import * as generic from './FpXX.js';
import * as word from './FlxX.js';
import * as binary from './F2x.js';

test('generic word division collapses scalar tags while its early remainder keeps them', () => {
  const T = [3n, 0n, 1n],
    a = [[1n]],
    b = [[1n], [1n]];
  expect(generic.FpXQX_divrem(a, b, T, 17n)).toEqual([[], [1n]]);
  expect(generic.FpXQX_rem(a, b, T, 17n)).toEqual([[1n]]);
  expect(generic.FpXQX_divrem(a, b, T, (1n << 64n) + 13n)).toEqual([[], [[1n]]]);
});

test('remainder-only constant division avoids inversion and preserves native zero errors', () => {
  const T = [3n, 0n, 1n];
  expect(word.FlxqX_rem([[1n]], [[2n]], T, 4n)).toEqual([]);
  expect(() => word.FlxqX_div([[1n]], [[2n]], T, 4n)).toThrow();
  expect(() => generic.FpXQX_divrem([], [], T, 17n)).toThrow(
    'impossible inverse in FlxqX_divrem: 0.'
  );
  expect(() => generic.FpXQX_divrem([], [], T, (1n << 64n) + 13n)).toThrow(
    'impossible inverse in FpX_divrem: 0.'
  );
});

test('native reduction objects are reused at the exact get_red boundary', () => {
  const T = [3n, 0n, 1n],
    S = Array.from({ length: 16 }, () => [1n]);
  const cached = word.FlxqX_get_red(S, T, 17n);
  expect(Array.isArray(cached)).toBe(false);
  expect(word.FlxqX_get_red(cached, T, 17n)).toBe(cached);
  expect(Array.isArray(word.FlxqX_get_red(S.slice(0, 15), T, 17n))).toBe(true);
  expect(word.FlxqX_divrem(S, cached, T, 17n)).toEqual([[[1n]], []]);
});

test('Newton reciprocal inversion delegates to native coefficient squaring', () => {
  const square = spyOn(binary, 'F2x_sqr');
  try {
    const S = new Array<bigint>(49).fill(1n);
    S[48] = 3n;
    binary.F2xqX_invBarrett(S, 7n);
    expect(square.mock.calls.length).toBeGreaterThan(0);
  } finally {
    square.mockRestore();
  }
});

test('binary remainder preserves the bundled PARI stale-buffer branch', () => {
  const a = [1n, ...new Array<bigint>(99).fill(0n), 1n],
    b = [0n, 0n, 1n];
  expect(binary.F2xqX_rem(a, b, 7n)).toEqual(a);
  expect(binary.F2xqX_divrem(a, b, 7n)[1]).toEqual([1n]);
  expect(binary.F2xqX_rem([1n, ...new Array<bigint>(98).fill(0n), 1n], b, 7n)).toEqual([1n]);
});

test('Barrett traversal copies a linear amount of large-array data', () => {
  const n = 4096,
    a = Array.from({ length: n + 1 }, (_, i) => BigInt(i % 4));
  a[n] = 1n;
  const slice = Array.prototype.slice;
  let copied = 0,
    result: [bigint[], bigint[]];
  Array.prototype.slice = function (start?: number, end?: number) {
    if (this.length >= n / 2) {
      const len = this.length,
        from =
          start === undefined ? 0 : start < 0 ? Math.max(len + start, 0) : Math.min(start, len);
      const to = end === undefined ? len : end < 0 ? Math.max(len + end, 0) : Math.min(end, len);
      copied += Math.max(0, to - from);
    }
    return slice.call(this, start, end);
  };
  try {
    result = binary.F2xqX_divrem(a, [1n, 1n, 1n], 7n);
  } finally {
    Array.prototype.slice = slice;
  }
  expect(result![0].length).toBe(n - 1);
  expect(copied).toBeLessThanOrEqual(20 * n);
});

test('binary Kronecker width admits coefficients at the modulus degree', () => {
  expect(binary.F2xqX_mul([7n], [1n], 7n)).toEqual([]);
  expect(binary.F2xqX_mul([7n, 1n], [7n, 1n], 7n)).toEqual([0n, 0n, 1n]);
  expect(() => binary.F2xqX_mul([1n << 64n], [1n], 7n)).toThrow(RangeError);
  expect(() => binary.F2xqX_mul([1n << 64n], [1n], 7n)).toThrow(
    'binary coefficient degree exceeds the native packing bound'
  );
});

test('word and binary basecase reciprocals preserve native raw initial terms', () => {
  const S = Array.from({ length: 20 }, () => [3n, 0n, 1n]);
  S[19] = [1n];
  const wordResult = word.FlxqX_invBarrett(S, [3n, 0n, 1n], 17n);
  expect(wordResult[1]).toEqual([14n, 0n, 16n]);
  expect(binary.F2xqX_invBarrett([...new Array<bigint>(19).fill(7n), 1n], 7n))
    .toEqual([1n, ...new Array<bigint>(17).fill(7n)]);
  expect(generic.FpXQX_invBarrett(S, [3n, 0n, 1n], (1n << 64n) + 13n))
    .toEqual([1n]);
});
