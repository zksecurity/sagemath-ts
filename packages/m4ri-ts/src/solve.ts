import { type mzd_t, mzd_init, mzd_submatrix as slice, mzd_add } from './mzd.js';
import { mzd_pluq } from './ple.js';
import { mzd_trsm_lower_left, mzd_trsm_upper_left } from './triangular.js';
import { mzd_mul } from './strassen.js';

/** solve.c:mzd_solve_left. Owned copies replace the two mutated native arguments.
 * Preserve the native consistency flag, including the skipped first padding row.
 * @see Deviation: Binary Linear Solve and Kernel Adapters
 */
export function mzd_solve_left(
  a: mzd_t,
  b: mzd_t,
  cutoff = 0,
  check = true
): {
  matrix: mzd_t;
  rhs: mzd_t;
  status: number;
} {
  if (a.ncols > b.nrows || b.nrows !== Math.max(a.nrows, a.ncols))
    throw new RangeError('M4RI solve requires max(A.nrows, A.ncols) right-hand rows');
  const nativeFault = () => {
    throw Object.assign(new Error('Segmentation fault'), { name: 'RuntimeError' });
  };
  if (!b.ncols && check && b.nrows > a.nrows + 1) nativeFault();
  if (check && b.nrows > a.nrows && b.rows.slice(a.nrows + 1).some((row) => row !== 0n))
    return {
      matrix: mzd_init(a.nrows, a.ncols, a.rows),
      rhs: mzd_init(b.nrows, b.ncols, b.rows),
      status: -1,
    };
  const factor = mzd_pluq(a, cutoff),
    rank = factor.rank,
    rows = Array.from(b.rows);
  for (let i = 0; i < a.nrows; i++) [rows[i], rows[factor.P[i]!]] = [rows[factor.P[i]!]!, rows[i]!];
  const LU = slice(factor.matrix, 0, 0, rank, rank);
  // Empty native B rows are indexed at width-1 by scalar substitutions and
  // zero tests. Russian tables also access that word even for diagonal LU.
  if (
    !b.ncols &&
    (rank > 64 ||
      LU.rows.some((row, i) => (row & ~(1n << BigInt(i))) !== 0n) ||
      (check && (rank < a.nrows || b.nrows > a.nrows)))
  )
    nativeFault();
  let upper = mzd_trsm_lower_left(LU, mzd_init(rank, b.ncols, rows.slice(0, rank)), cutoff);
  let status = 0;
  if (check) {
    const lower = mzd_add(
      mzd_init(a.nrows - rank, b.ncols, rows.slice(rank, a.nrows)),
      mzd_mul(slice(factor.matrix, rank, 0, a.nrows, rank), upper, cutoff)
    );
    for (let i = rank; i < a.nrows; i++) rows[i] = lower.rows[i - rank]!;
    for (let i = a.nrows; i < b.nrows; i++) rows[i] = 0n;
    if (lower.rows.some((row) => row !== 0n)) status = -1;
  }
  upper = mzd_trsm_upper_left(LU, upper, cutoff);
  for (let i = 0; i < rank; i++) rows[i] = upper.rows[i]!;
  if (!check) for (let i = rank; i < b.nrows; i++) rows[i] = 0n;
  for (let i = a.ncols - 1; i >= 0; i--)
    [rows[i], rows[factor.Q[i]!]] = [rows[factor.Q[i]!]!, rows[i]!];
  return { matrix: factor.matrix, rhs: mzd_init(b.nrows, b.ncols, rows), status };
}

/** solve.c:mzd_kernel_left_pluq. The kernel's columns solve A*X=0.
 * Return the factored input and kernel, or null when the native kernel is trivial.
 * @see Deviation: Binary Linear Solve and Kernel Adapters
 */
export function mzd_kernel_left_pluq(
  a: mzd_t,
  cutoff = 0
): { matrix: mzd_t; kernel: mzd_t | null } {
  const factor = mzd_pluq(a, cutoff),
    rank = factor.rank;
  if (rank === a.ncols) return { matrix: factor.matrix, kernel: null };
  const free = a.ncols - rank;
  const upper = mzd_trsm_upper_left(
    slice(factor.matrix, 0, 0, rank, rank),
    slice(factor.matrix, 0, rank, rank, a.ncols),
    cutoff
  );
  const rows = [...upper.rows, ...Array.from({ length: free }, (_, i) => 1n << BigInt(i))];
  for (let i = a.ncols - 1; i >= 0; i--)
    [rows[i], rows[factor.Q[i]!]] = [rows[factor.Q[i]!]!, rows[i]!];
  return { matrix: factor.matrix, kernel: mzd_init(a.ncols, free, rows) };
}
