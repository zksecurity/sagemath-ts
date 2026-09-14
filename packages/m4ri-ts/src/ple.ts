import { type mzd_t, mzd_init, mzd_add, mzd_submatrix as slice } from './mzd.js';
import { _mzd_ple_russian, swapColumns, type PLE } from './ple_russian.js';
import { mzd_trsm_lower_left } from './triangular.js';
import { mzd_mul } from './strassen.js';
/** ple.c: recursive PLE with a Russian base case and Schur complement. */
export function mzd_ple(a: mzd_t, cutoff = 0): PLE {
  const m = a.nrows,
    n = a.ncols;
  if (m > 0 && n === 0)
    throw Object.assign(new Error('Segmentation fault'), { name: 'RuntimeError' });
  let active = m;
  while (active > 0 && a.rows[active - 1] === 0n) active--;
  if (!active || !n)
    return {
      matrix: mzd_init(m, n, a.rows),
      rank: 0,
      P: Array.from({ length: m }, (_, i) => i),
      Q: Array.from({ length: n }, (_, i) => i),
    };
  if (n <= 64 || Math.ceil(n / 64) * m <= 524288) return _mzd_ple_russian(a);
  const split = Math.floor(Math.ceil(n / 64) / 2) * 64;
  const first = mzd_ple(slice(a, 0, 0, active, split), cutoff),
    r1 = first.rank;
  let right = Array.from(slice(a, 0, split, active, n).rows);
  if (r1)
    for (let i = 0; i < active; i++)
      [right[i], right[first.P[i]!]] = [right[first.P[i]!]!, right[i]!];
  let upper = slice(mzd_init(active, n - split, right), 0, 0, r1, n - split);
  let lower = slice(mzd_init(active, n - split, right), r1, 0, active, n - split);
  if (r1) {
    upper = mzd_trsm_lower_left(slice(first.matrix, 0, 0, r1, r1), upper, cutoff);
    lower = mzd_add(lower, mzd_mul(slice(first.matrix, r1, 0, active, r1), upper, cutoff));
  }
  const second = mzd_ple(lower, cutoff),
    r2 = second.rank;
  const left = Array.from(first.matrix.rows);
  const lowMask = (1n << BigInt(r1)) - 1n;
  for (let i = 0; i < active - r1; i++) {
    const j = second.P[i]!,
      delta = (left[r1 + i]! ^ left[r1 + j]!) & lowMask;
    left[r1 + i] ^= delta;
    left[r1 + j] ^= delta;
  }
  const rows = Array.from({ length: m }, (_, i) =>
    i >= active
      ? a.rows[i]!
      : left[i]! | ((i < r1 ? upper.rows[i]! : second.matrix.rows[i - r1]!) << BigInt(split))
  );
  const P = Array.from({ length: m }, (_, i) =>
    i < r1 ? first.P[i]! : i < active ? second.P[i - r1]! + r1 : i
  );
  const Q = Array.from({ length: n }, (_, i) =>
    i < split ? first.Q[i]! : second.Q[i - split]! + split
  );
  for (let j = 0; j < r2; j++) Q[r1 + j] = second.Q[j]! + split;
  if (r1 !== split) {
    for (let j = 0; j < r2; j++) swapColumns(rows, r1 + j, split + j, r1 + j, r1 + r2);
    const keep = (1n << BigInt(r1)) - 1n,
      mask = (1n << BigInt(r2)) - 1n;
    for (let i = r1 + r2; i < active; i++)
      rows[i] = (rows[i]! & keep) | (((rows[i]! >> BigInt(split)) & mask) << BigInt(r1));
  }
  return { matrix: mzd_init(m, n, rows), rank: r1 + r2, P, Q };
}
/** ple.c: PLUQ additionally permutes the triangular upper part. */
export function mzd_pluq(a: mzd_t, cutoff = 0): PLE {
  const result = mzd_ple(a, cutoff),
    rows = Array.from(result.matrix.rows);
  for (let j = 0; j < a.ncols; j++) swapColumns(rows, j, result.Q[j]!, 0, Math.min(j, result.rank));
  return { ...result, matrix: mzd_init(a.nrows, a.ncols, rows) };
}

/** ple.c: explicitly requested scalar PLE/PLUQ algorithms. */
function naive(a: mzd_t, pluq: boolean): PLE {
  const rows = Array.from(a.rows),
    m = a.nrows,
    n = a.ncols,
    P = Array.from({ length: m }, (_, i) => i),
    Q = Array.from({ length: n }, (_, i) => i);
  let rank = 0,
    column = 0;
  while (rank < m && column < n) {
    let found = -1,
      pivot = n;
    for (let j = column; j < n && found < 0; j++)
      for (let i = rank; i < m; i++)
        if ((rows[i]! >> BigInt(j)) & 1n) {
          found = i;
          pivot = j;
          break;
        }
    if (found < 0) break;
    P[rank] = found;
    Q[rank] = pivot;
    [rows[rank], rows[found]] = [rows[found]!, rows[rank]!];
    if (pluq) {
      swapColumns(rows, rank, pivot);
      pivot = rank;
    }
    for (let i = rank + 1; i < m; i++)
      if ((rows[i]! >> BigInt(pivot)) & 1n) rows[i] ^= rows[rank]! & (-1n << BigInt(pivot + 1));
    rank++;
    column = pivot + 1;
  }
  if (!pluq) for (let j = 0; j < rank; j++) if (Q[j]! > j) swapColumns(rows, Q[j]!, j, j, m);
  return { matrix: mzd_init(m, n, rows), rank, P, Q };
}
export function _mzd_ple_naive(a: mzd_t): PLE {
  return naive(a, false);
}
export function _mzd_pluq_naive(a: mzd_t): PLE {
  return naive(a, true);
}
