import { expect, spyOn, test } from 'bun:test';
import { FpXQXQ_minpoly } from './FpXX.js';
import { FlxqXQ_minpoly } from './FlxX.js';
import { getrand, setrand } from './random.js';
import * as division from './_extension_division.js';
import * as quotient from './_extension_quotient.js';

test('extension minimal polynomials preserve native coefficient tags', () => {
  const saved = getrand();
  try {
    setrand(1n);
    expect(FpXQXQ_minpoly([0n, 1n], [1n, 0n, 1n], [3n, 0n, 1n], 17n)).toEqual([[1n], [], 1n]);
    setrand(1n);
    expect(FlxqXQ_minpoly([[], [1n]], [[1n], [], [1n]], [3n, 0n, 1n], 17n)).toEqual([
      [1n],
      [],
      [1n],
    ]);
    setrand(1n);
    expect(FpXQXQ_minpoly([[1n, 1n]], [1n, 0n, 1n], [3n, 0n, 1n], 17n)).toEqual([[16n, 16n], 1n]);
  } finally {
    setrand(saved);
  }
});

test('minimal polynomial retains multiplicity for a repeated outer modulus', () => {
  const saved = getrand();
  try {
    setrand(1n);
    expect(FpXQXQ_minpoly([0n, 1n], [1n, 0n, 0n, 0n, 1n], [1n, 1n, 1n], 2n)).toEqual([
      [1n],
      [],
      [],
      [],
      1n,
    ]);
    setrand(1n);
    expect(FlxqXQ_minpoly([[], [1n]], [[1n], [], [], [], [1n]], [1n, 1n, 1n], 2n)).toEqual([
      [1n],
      [],
      [],
      [],
      [1n],
    ]);
  } finally {
    setrand(saved);
  }
});

test('native outer reciprocal errors precede the power table inner preparation', () => {
  const genericT = [1n, ...Array<bigint>(35).fill(0n), 2n];
  const wordT = [1n, ...Array<bigint>(90).fill(0n), 2n];
  const saved = getrand();
  for (const [degree, inverse] of [
    [9, 2],
    [10, 3],
  ]) {
    const S = [[1n], ...Array.from({ length: degree! - 1 }, () => [] as bigint[]), [3n]];
    expect(() => FpXQXQ_minpoly([], S, genericT, 6n)).toThrow(
      `impossible inverse in ${degree === 9 ? 'Fp_inv' : 'Fl_inv'}: Mod(${inverse}, 6).`
    );
  }
  for (const [degree, inverse] of [
    [14, 2],
    [15, 3],
  ]) {
    const S = [[1n], ...Array.from({ length: degree! - 1 }, () => [] as bigint[]), [3n]];
    expect(() => FlxqXQ_minpoly([], S, wordT, 6n)).toThrow(
      `impossible inverse in Fl_inv: Mod(${inverse}, 6).`
    );
  }
  expect(() => FpXQXQ_minpoly([], [1n, ...Array<bigint>(9).fill(0n), 3n], genericT, 6n)).toThrow(
    'impossible inverse in Fp_inv: Mod(3, 6).'
  );
  expect(getrand()).toBe(saved);
});

test('minimal polynomial reuses one outer cache and one initial power table', () => {
  const S = [1n, ...Array<bigint>(9).fill(0n), 1n],
    saved = getrand();
  const outer = spyOn(division, 'extensionGetRed'),
    powers = spyOn(quotient, 'extensionQuotient');
  try {
    setrand(1n);
    expect(FpXQXQ_minpoly([0n, 1n], S, [3n, 0n, 1n], 17n)).toEqual([
      [1n],
      ...Array.from({ length: 9 }, () => []),
      1n,
    ]);
    const cache = outer.mock.results[0]!.value;
    expect(Array.isArray(cache)).toBe(false);
    expect(outer.mock.calls.length).toBeGreaterThan(1);
    expect(outer.mock.calls.slice(1).every((args) => args[1] === cache)).toBe(true);
    expect(powers.mock.calls.filter((args) => args[1] === 6).length).toBe(1);
  } finally {
    outer.mockRestore();
    powers.mockRestore();
    setrand(saved);
  }
});
