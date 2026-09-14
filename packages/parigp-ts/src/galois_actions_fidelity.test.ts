import { expect, test } from 'bun:test';
import { PariError } from './errors.js';
import { galois_group, galoisfixedfield, galoisinit, galoissubgroups } from './galconj.js';

function nativeError(fn: () => unknown, message: string) {
  let caught: unknown;
  try {
    fn();
  } catch (e) {
    caught = e;
  }
  expect(caught?.constructor).toBe(PariError);
  expect((caught as Error).message).toBe(message);
}

test('fixed-field flag validation precedes permutation shape validation', () => {
  const gal = galoisinit([1n, 0n, 1n])!;
  for (const flag of [-1, 3])
    nativeError(() => galoisfixedfield(gal, [], flag as 0), 'invalid flag in galoisfixedfield.');
});

test('fixed-field permutation length errors preserve the native type', () => {
  const gal = galoisinit([1n, 0n, 1n])!;
  for (const perm of [[0], [0, 1], [0, 1, 2, 1]])
    nativeError(
      () => galoisfixedfield(gal, perm, 0),
      'incorrect type in galoisfixedfield (t_VECSMALL).'
    );
});

test('malformed generator/order pairs follow each native group validator', () => {
  const gal = galoisinit([1n, 0n, 1n])!;
  const group = { gen: [[], [0, 2, 1]], ord: [0] };
  for (const flag of [0, 1, 2])
    nativeError(
      () => galoisfixedfield(gal, group, flag as 0),
      'incorrect type in galoisfixedfield (t_VEC).'
    );
  nativeError(() => galoissubgroups(group), 'incorrect type in checkgal (t_VEC).');
});

test('full-group flag-two factorization handles a vanishing trace bound', () => {
  for (const pol of [
    [1n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 1n],
    [1n, -1n, 0n, 1n, -1n, 1n, 0n, -1n, 1n],
  ]) {
    const gal = galoisinit(pol)!;
    const full = galois_group(gal);
    const expected = {
      P: pol[1] === 0n ? [0n, 1n] : [-1n, 1n],
      S: { num: pol[1] === 0n ? [] : [1n], den: 1n },
      factors: [pol.map((c) => ({ num: c ? [c] : [], den: 1n }))],
    };
    expect(galoisfixedfield(gal, full, 2)).toEqual(expected);
    expect(galoisfixedfield(gal, full.gen, 2)).toEqual(expected);
  }
});
