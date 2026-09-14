import { expect, test } from 'bun:test';
import {
  mzd_init,
  mzd_ple,
  mzd_pluq,
  mzd_echelonize,
  mzd_echelonize_m4ri,
  mzd_echelonize_pluq,
  mzd_trsm_upper_left,
  mzd_trsm_lower_left,
} from './index.js';

test('Sage binary factorization doctest retains compact factors and transpositions', () => {
  const A = mzd_init(4, 4, [10n, 14n, 8n, 6n]);
  const E = mzd_ple(A),
    U = mzd_pluq(A);
  expect(E.rank).toBe(3);
  expect(E.matrix.rows).toEqual([9n, 3n, 4n, 7n]);
  expect(U.matrix.rows).toEqual([5n, 3n, 4n, 7n]);
  expect(E.P).toEqual([0, 1, 2, 3]);
  expect(E.Q).toEqual([1, 2, 3, 3]);
  for (const reduce of [mzd_echelonize, mzd_echelonize_m4ri, mzd_echelonize_pluq])
    expect(reduce(A).matrix.rows).toEqual([2n, 4n, 8n, 0n]);
  expect(A.rows).toEqual([10n, 14n, 8n, 6n]);
});

test('native triangular adapters solve unit triangular systems', () => {
  const B = mzd_init(2, 2, [1n, 2n]);
  expect(mzd_trsm_upper_left(mzd_init(2, 2, [3n, 2n]), B).rows).toEqual([3n, 2n]);
  expect(mzd_trsm_lower_left(mzd_init(2, 2, [1n, 3n]), B).rows).toEqual([1n, 3n]);
  expect(() => mzd_trsm_upper_left(mzd_init(1, 2), B)).toThrow('incompatible triangular');
});
