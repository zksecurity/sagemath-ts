import { expect, test } from 'bun:test';
import { _fmpz_poly_pseudo_divrem_basecase as basecase } from './pseudo_divrem_basecase.js';
import { _fmpz_poly_pseudo_divrem_divconquer as divide } from './pseudo_divrem_divconquer.js';
import { _fmpz_poly_pseudo_rem as rem } from './pseudo_rem.js';
import { _fmpq_poly_rem as rationalRem } from '../fmpq_poly/rem.js';

test('pseudo division tracks required leading-coefficient powers and raw zero padding', () => {
  expect(basecase([1n, 0n, 1n], [1n, 2n])).toEqual([[-1n, 2n], [5n], 2n]);
  expect(divide([1n, 0n, 1n], [1n, 2n])).toEqual([[-1n, 2n], [5n], 2n]);
  expect(rem([1n, 0n, 1n], [1n, 2n])).toEqual([[5n], 2n]);
  expect(basecase([0n, 0n, 1n], [0n, 1n])).toEqual([[0n, 1n], [0n], 0n]);
  expect(divide([0n, 0n, 1n], [0n, 1n])).toEqual([[0n, 1n], [], 0n]);
});

test('rational remainder clears and canonicalizes denominators', () => {
  expect(rationalRem([1n, 0n, 1n], 7n, [1n, 2n], 3n)).toEqual([[5n], 28n]);
  expect(rationalRem([2n, 4n], 6n, [1n, 0n, 1n], 3n)).toEqual([[1n, 2n], 3n]);
  expect(rationalRem([], 7n, [2n], 3n)).toEqual([[], 1n]);
  expect(rationalRem([1n, 2n], 7n, [2n], 3n)).toEqual([[], 1n]);
});

test('array adapters reject invalid denominator and divisor buffers', () => {
  expect(() => basecase([1n], [])).toThrow(RangeError);
  expect(() => basecase([1n], [1n, 0n])).toThrow(RangeError);
  expect(() => divide([1n], [0n])).toThrow(RangeError);
  expect(() => rationalRem([1n], 0n, [1n], 1n)).toThrow(RangeError);
  expect(() => rationalRem([1n], 1n, [1n], -1n)).toThrow(RangeError);
  expect(() => rationalRem([1n], 1n, [], 1n)).toThrow(RangeError);
});
