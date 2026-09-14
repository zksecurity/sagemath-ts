import type { mzd_t } from './mzd.js';
import { _mzd_echelonize_m4ri } from './brilliantrussian.js';
/** echelonform.c: explicit M4RI elimination. Returns an owned matrix and rank. */
export function mzd_echelonize_m4ri(a: mzd_t, full = true, k = 0): { matrix: mzd_t; rank: number } {
  return _mzd_echelonize_m4ri(a, full, k);
}

import { mzd_init, mzd_submatrix } from './mzd.js';
import { mzd_ple, mzd_pluq } from './ple.js';
import { swapColumns } from './ple_russian.js';
import { mzd_trsm_upper_left } from './triangular.js';
/** echelonform.c: factorization, triangular solve, and inverse permutation. */
export function mzd_echelonize_pluq(a: mzd_t, full = true): { matrix: mzd_t; rank: number } {
  if (!a.nrows || !a.ncols) return { matrix: mzd_init(a.nrows, a.ncols, a.rows), rank: 0 };
  const result = full ? mzd_pluq(a) : mzd_ple(a),
    r = result.rank,
    rows = Array.from(result.matrix.rows);
  if (full) {
    if (r !== a.ncols) {
      const U = mzd_submatrix(result.matrix, 0, 0, r, r),
        B = mzd_submatrix(result.matrix, 0, r, r, a.ncols);
      const solved = mzd_trsm_upper_left(U, B);
      for (let i = 0; i < r; i++) rows[i] = (1n << BigInt(i)) | (solved.rows[i]! << BigInt(r));
    } else for (let i = 0; i < r; i++) rows[i] = 1n << BigInt(i);
    for (let j = a.ncols - 1; j >= 0; j--) swapColumns(rows, j, result.Q[j]!, 0, r || a.nrows);
  } else
    for (let i = 0; i < r; i++)
      rows[i] = (rows[i]! & (-1n << BigInt(i + 1))) | (1n << BigInt(result.Q[i]!));
  for (let i = r; i < a.nrows; i++) rows[i] = 0n;
  return { matrix: mzd_init(a.nrows, a.ncols, rows), rank: r };
}
/** echelonform.c: M4RI panels with original density-based PLUQ crossover. */
export function mzd_echelonize(a: mzd_t, full = true): { matrix: mzd_t; rank: number } {
  return _mzd_echelonize_m4ri(a, full, 0, mzd_echelonize_pluq);
}
