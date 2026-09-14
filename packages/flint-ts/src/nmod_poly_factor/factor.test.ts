import { expect, test } from 'bun:test';
import {
  flint_rand_init,
  flint_rand_get_seed,
  nmod_poly_factor,
  nmod_poly_factor_equal_deg_prob,
  nmod_poly_deflate,
  nmod_poly_inflate,
  nmod_poly_remove,
  _nmod_vec_randtest,
} from '../index.js';

test('bundled FLINT native factor insertion order and split random state', () => {
  // Direct bundled C results; the full comparative area also checks every output word.
  const f = [4n, 0n, 0n, 0n, 1n];
  expect(nmod_poly_factor(f, 5n)).toEqual([
    1n,
    [
      [[4n, 1n], 1],
      [[2n, 1n], 1],
      [[3n, 1n], 1],
      [[1n, 1n], 1],
    ],
  ]);
  const state = flint_rand_init();
  expect(nmod_poly_factor_equal_deg_prob(state, f, 1, 5n)).toEqual([4n, 1n]);
  expect(flint_rand_get_seed(state)).toEqual([12846242740235074558n, 16846056370759765703n]);
});

test('deflation preserves native storage length; inflation zero evaluates at one', () => {
  expect(nmod_poly_deflate([1n, 0n, 0n, 1n], 2, 5n)).toEqual([1n, 0n]);
  expect(nmod_poly_inflate([1n, 2n, 2n], 0, 5n)).toEqual([]);
  expect(nmod_poly_remove([1n, 2n, 1n], [1n, 1n], 5n)).toEqual([[1n], 2]);
});

test('unsafe FLINT factor inputs fail before consuming random state', () => {
  const state = flint_rand_init(),
    before = flint_rand_get_seed(state);
  expect(() => _nmod_vec_randtest(state, -1, 5n)).toThrow('length must be nonnegative');
  expect(() => nmod_poly_factor_equal_deg_prob(state, [1n, 1n], 1, 5n)).toThrow(
    'degree at least two'
  );
  expect(flint_rand_get_seed(state)).toEqual(before);
  expect(() => nmod_poly_deflate([1n], 0, 5n)).toThrow('deflation must be positive');
  expect(() => nmod_poly_inflate([1n], -1, 5n)).toThrow('inflation must be nonnegative');
  expect(() => nmod_poly_remove([1n], [1n], 5n)).toThrow('factor must have positive degree');
});
