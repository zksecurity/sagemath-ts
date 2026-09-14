import { expect, test, spyOn } from 'bun:test';
import * as generic from './FpXX.js';
import * as word from './FlxX.js';
import * as binary from './F2x.js';
import * as division from './_extension_division.js';
import { extensionGcd } from './_extension_gcd.js';

test('extension gcd retains the unscaled native last remainder', () => {
  const T = [3n, 0n, 1n];
  expect(generic.FpXQX_gcd([2n, 2n], [4n, 4n], T, 17n)).toEqual([4n, 4n]);
  expect(word.FlxqX_gcd([[2n], [2n]], [[4n], [4n]], T, 17n)).toEqual([[4n], [4n]]);
  expect(generic.FpXQX_extgcd([1n, 0n, 1n], [1n, 1n], T, 17n)).toEqual([[2n], [1n], [1n, 16n]]);
});

test('native V-only output avoids a zero-first-input quotient error', () => {
  const T = [3n, 0n, 1n];
  expect(() => word.FlxqX_extgcd([], [[1n]], T, 17n)).toThrow(
    'impossible inverse in FlxqX_divrem: 0.'
  );
  expect(extensionGcd(1, 3, 17n, T, [], [[1n]])).toEqual([[[1n]], [[1n]]]);
  expect(extensionGcd(2, 3, 2n, 7n, [], [])).toEqual([[], []]);
});

test('raw binary half-GCD follows native coefficient packing', () => {
  expect(binary.F2xqX_halfgcd([1n], [8n, 1n], 7n)).toEqual([
    [[1n], []],
    [[1n, 1n], [1n]],
  ]);
  expect(binary.F2xqX_mul([32n, 1n], [1n], 7n)).toEqual([]);
  expect(binary.F2xqX_mul([32n], [1n], 7n)).toEqual([0n, 1n]);
});

test('generic word half-GCD delegates its Euclidean divisions to the word backend', () => {
  const spy = spyOn(division, 'extensionDivision');
  try {
    generic.FpXQX_halfgcd([[1n], [1n], [1n]], [[2n], [1n]], [3n, 0n, 1n], 17n);
    expect(spy.mock.calls.length).toBeGreaterThan(0);
    expect(spy.mock.calls.every((args) => args[0] === 1)).toBe(true);
  } finally {
    spy.mockRestore();
  }
});
