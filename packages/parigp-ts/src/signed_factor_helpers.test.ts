import { expect, test } from 'bun:test';
import { is_kth_power } from './ifactor.js';
import { PariError } from './errors.js';
import { mpqsInternals } from './mpqs.js';

test('negative k-th powers preserve the native unsupported-root error', () => {
  for (const [n, k] of [
    [-8n, 3],
    [-27n, 1],
    [-1n, 11],
  ] as const) {
    let caught: unknown;
    try {
      is_kth_power(n, k);
    } catch (e) {
      caught = e;
    }
    expect(caught?.constructor).toBe(PariError);
    expect((caught as Error).message).toBe('sorry, sqrtnr for x < 0 is not yet implemented.');
  }
});
test('negative k-th powers still stop at the native modular filters', () => {
  expect(is_kth_power(-16n, 3)).toBeNull();
  expect(is_kth_power(-1n, 2)).toBeNull();
});
test('MPQS class-group relations use the magnitude of their discriminant for inversion', () => {
  const h = mpqsInternals.newHandle();
  h.N = -23n;
  h.size_of_FB = 6;
  mpqsInternals.mpqs_FB_ctor(h).p.set([0, 0, 2, 3, 5, 7, 11, 13]);
  expect(
    mpqsInternals.combine_large_primes(
      h,
      29,
      { Y: 79n, relp: [1048578, 3145731] },
      { Y: 21n, relp: [2097154] },
      1
    )
  ).toEqual({ Y: 0n, relp: [-3145726, 3145731] });
});

test('MPQS factor-base square roots select the native smaller root', () => {
  expect(mpqsInternals.Fl_sqrt(4, 5)).toBe(2);
  expect(mpqsInternals.Fl_sqrt(2, 7)).toBe(3);
  expect(mpqsInternals.Fl_sqrt(3, 13)).toBe(4);
});

test('large root exponents use the exact bit bound without allocating a giant power', () => {
  expect(is_kth_power(2n, 17886697)).toBeNull();
  expect(is_kth_power(2n, 17886698)).toBeNull();
});
