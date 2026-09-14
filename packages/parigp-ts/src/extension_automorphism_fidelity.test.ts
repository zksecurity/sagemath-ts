import { expect, spyOn, test } from 'bun:test';
import { FpXQXQ_autpow, FpXQXQ_auttrace, FpXQXQ_autsum } from './FpXX.js';
import { FlxqXQ_auttrace, FlxqXQ_autpow } from './FlxX.js';
import { F2xqXQ_auttrace } from './F2x.js';
import * as generic from './FpX.js';
import * as division from './_extension_division.js';

test('native extension trace adds while autsum multiplies its aggregate', () => {
  const T = [1n, 0n, 1n],
    S = [1n, 0n, 1n],
    phi = [0n, 1n],
    B = [0n, 1n];
  expect(FpXQXQ_auttrace([B, [2n, 1n]], 3n, S, T, 17n)).toEqual([B, [6n, 3n]]);
  expect(FpXQXQ_autsum([phi, B, [2n, 1n]], 3n, S, T, 17n)).toEqual([phi, B, [2n, 11n]]);
  expect(
    FlxqXQ_auttrace(
      [
        [[], [1n]],
        [[1n], [1n]],
      ],
      3n,
      [[1n], [], [1n]],
      T,
      17n
    )
  ).toEqual([
    [[], [1n]],
    [[3n], [3n]],
  ]);
  expect(F2xqXQ_auttrace([2n, [0n, 1n], [1n, 1n]], 3n, [1n, 0n, 1n], 7n)).toEqual([
    2n,
    [0n, 1n],
    [1n, 1n],
  ]);
});

test('signed native automorphism exponents preserve the unsigned C cast', () => {
  const B = [0n, 1n],
    S = [1n, 0n, 1n],
    T = [1n, 0n, 1n];
  expect(FpXQXQ_auttrace([B, [1n]], -1n, S, T, 17n)).toEqual([B, []]);
  expect(FpXQXQ_auttrace([B, [1n]], -(1n << 63n), S, T, 17n)).toEqual([B, [9n]]);
  expect(() => FpXQXQ_auttrace([B, [1n]], 0n, S, T, 17n)).toThrow(
    'automorphism exponent must be a nonzero signed word'
  );
});

test('unit-count extension automorphisms still prepare native reciprocals', () => {
  const T = [1n, 0n, 1n],
    S = [1n, ...Array<bigint>(9).fill(0n), 2n];
  expect(() =>
    FpXQXQ_autpow(
      [
        [0n, 1n],
        [0n, 1n],
      ],
      1n,
      S,
      T,
      4n
    )
  ).toThrow('impossible inverse in Fp_inv: Mod(2, 4).');
  const wordS = [[1n], ...Array.from({ length: 14 }, () => [] as bigint[]), [2n]];
  expect(() =>
    FlxqXQ_autpow(
      [
        [0n, 1n],
        [[], [1n]],
      ],
      1n,
      wordS,
      T,
      4n
    )
  ).toThrow('impossible inverse in Fl_inv: Mod(2, 4).');
  expect(
    FpXQXQ_autpow(
      [
        [0n, 1n, 0n],
        [18n, [19n, 0n], 0n],
      ],
      1n,
      [1n],
      T,
      17n
    )
  ).toEqual([
    [0n, 1n],
    [18n, [19n]],
  ]);
});

test('prepared inner and outer reciprocals survive successive automorphism compositions', () => {
  const T = [3n, ...Array<bigint>(35).fill(0n), 1n];
  const S = [1n, ...Array<bigint>(9).fill(0n), 1n];
  const B = [...Array<bigint>(9).fill(0n), 1n];
  const inner = spyOn(generic, 'FpX_invBarrett');
  const outer = spyOn(division, 'extensionGetRed');
  try {
    expect(FpXQXQ_autpow([[0n, 1n], B], 5n, S, T, 17n)).toEqual([[0n, 1n], B]);
    expect(inner).toHaveBeenCalledTimes(1);
    const cache = outer.mock.results[0]!.value;
    expect(Array.isArray(cache)).toBe(false);
    expect(outer.mock.calls.length).toBeGreaterThan(1);
    expect(outer.mock.calls.slice(1).every((args) => args[1] === cache)).toBe(true);
  } finally {
    inner.mockRestore();
    outer.mockRestore();
  }
});
