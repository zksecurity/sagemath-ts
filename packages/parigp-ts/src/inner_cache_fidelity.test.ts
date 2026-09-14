import { expect, spyOn, test } from 'bun:test';
import * as generic from './FpX.js';
import * as word from './Flx.js';

test('native reciprocal preparation precedes power-table zero and unit shortcuts', () => {
  const genericT = [1n, ...new Array<bigint>(35).fill(0n), 2n];
  const wordT = [1n, ...new Array<bigint>(90).fill(0n), 2n];
  for (const n of [0, 1, 2]) {
    expect(() => generic.FpXQ_powers([1n], n, genericT, 4n)).toThrow(
      'impossible inverse in Fp_inv: Mod(2, 4).'
    );
    expect(() => word.Flxq_powers([1n], n, wordT, 4n)).toThrow(
      'impossible inverse in Fl_inv: Mod(2, 4).'
    );
  }
});

test('native constant evaluation prepares a reciprocal but canonical zero skips it', () => {
  const T = [1n, ...new Array<bigint>(90).fill(0n), 2n];
  expect(() => generic.FpX_FpXQV_eval([1n], [[1n]], T, 4n)).toThrow(
    'impossible inverse in Fp_inv: Mod(2, 4).'
  );
  expect(() => word.Flx_FlxqV_eval([1n], [[1n]], T, 4n)).toThrow(
    'impossible inverse in Fl_inv: Mod(2, 4).'
  );
  expect(generic.FpX_FpXQV_eval([0n, 0n], [], T, 4n)).toEqual([]);
  expect(word.Flx_FlxqV_eval([0n, 0n], [], T, 4n)).toEqual([]);
  expect(generic.FpX_FpXQ_eval([0n, 0n], [1n], T, 4n)).toEqual([]);
  expect(word.Flx_Flxq_eval([0n, 0n], [1n], T, 4n)).toEqual([]);
});

test('canonical storage retains unreduced first coefficients and uses word arithmetic', () => {
  expect(generic.FpXQ_powers([18n, 0n], 1, [3n, 0n, 1n], 17n)).toEqual([[1n], [18n]]);
  expect(word.Flxq_powers([1n, 0n], 1, [3n, 0n, 1n], 17n)).toEqual([[1n], [1n]]);
  const cache = spyOn(word, 'Flx_invBarrett'),
    square = spyOn(word, 'Flx_sqr');
  try {
    const T = [1n, ...new Array<bigint>(90).fill(0n), 1n];
    expect(word.Flxq_powers([1n, 1n], 3, T, 17n)).toEqual([
      [1n],
      [1n, 1n],
      [1n, 2n, 1n],
      [1n, 3n, 3n, 1n],
    ]);
    expect(cache.mock.calls.length).toBe(1);
    expect(square.mock.calls.length).toBeGreaterThan(0);
  } finally {
    cache.mockRestore();
    square.mockRestore();
  }
});
