import { type mzd_t, mzd_init, mzd_add, mzd_submatrix as slice } from './mzd.js';
import { mzd_mul } from './strassen.js';
import { trsmRussian } from './triangular_russian.js';
/** triangular.c: scalar base case, Russian middle range, block recursion. */
function solve(a: mzd_t, b: mzd_t, upper: boolean, cutoff: number): mzd_t {
  const m = b.nrows,
    n = b.ncols;
  if (a.nrows !== a.ncols || a.ncols !== m)
    throw new RangeError('incompatible triangular solve dimensions');
  if (!m || !n) return mzd_init(m, n, b.rows);
  if (m <= 64) {
    const rows = Array.from(b.rows);
    for (let z = 0; z < m; z++) {
      const i = upper ? m - z - 1 : z;
      for (let j = upper ? i + 1 : 0; j < (upper ? m : i); j++)
        if ((a.rows[i]! >> BigInt(j)) & 1n) rows[i] ^= rows[j]!;
    }
    return mzd_init(m, n, rows);
  }
  if (m <= 2048) return trsmRussian(a, b, upper);
  const split = Math.floor(Math.ceil(m / 64) / 2) * 64;
  const a0 = slice(a, 0, 0, split, split),
    a1 = slice(a, split, split, m, m);
  let b0 = slice(b, 0, 0, split, n),
    b1 = slice(b, split, 0, m, n);
  if (upper) {
    b1 = solve(a1, b1, true, cutoff);
    b0 = solve(a0, mzd_add(b0, mzd_mul(slice(a, 0, split, split, m), b1, cutoff)), true, cutoff);
  } else {
    b0 = solve(a0, b0, false, cutoff);
    b1 = solve(a1, mzd_add(b1, mzd_mul(slice(a, split, 0, m, split), b0, cutoff)), false, cutoff);
  }
  return mzd_init(m, n, [...b0.rows, ...b1.rows]);
}
export function mzd_trsm_upper_left(a: mzd_t, b: mzd_t, cutoff = 0): mzd_t {
  return solve(a, b, true, cutoff);
}
export function mzd_trsm_lower_left(a: mzd_t, b: mzd_t, cutoff = 0): mzd_t {
  return solve(a, b, false, cutoff);
}
