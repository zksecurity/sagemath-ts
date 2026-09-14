import { expect, test } from 'bun:test';
import { FpXQ_powers, FpXQ_autpow, FpXQ_autpowers } from './galconj.js';

test('PARI quotient power vectors retain their initial unreduced entries', () => {
  const x = [1n, 1n, 1n];
  expect(FpXQ_powers(x, 2, [1n], 3n)).toEqual([[1n], x, []]);
  expect(FpXQ_autpowers(x, 1, [1n], 3n)).toEqual([[], [0n, 1n], x]);
  const v = FpXQ_powers(x, 1, [1n], 3n);
  v[1]![0] = 2n;
  expect(x).toEqual([1n, 1n, 1n]);
});

test('PARI automorphism powering preserves native binary composition order', () => {
  // Actual bundled PARI result. For a general substitution, quotient reduction
  // makes binary composition observably different from repeated substitution.
  expect(FpXQ_autpow([1n, 1n, 0n, 1n, 1n], 4, [1n, 0n, 0n, 1n, 0n, 1n], 2n)).toEqual([
    0n,
    1n,
    1n,
    1n,
  ]);
  expect(FpXQ_autpow([0n, 1n], Number.MAX_SAFE_INTEGER, [1n, 1n, 0n, 1n], 2n)).toEqual([0n, 1n]);
  for (const fn of [FpXQ_powers, FpXQ_autpow, FpXQ_autpowers]) {
    expect(() => fn([1n], -1, [1n, 0n, 1n], 3n)).toThrow('power count must be nonnegative');
  }
});
