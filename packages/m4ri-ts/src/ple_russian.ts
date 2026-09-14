import { nativeReadBits } from './brilliantrussian.js';
import { type mzd_t, mzd_init } from './mzd.js';

export interface PLE {
  matrix: mzd_t;
  rank: number;
  P: number[];
  Q: number[];
}
const bits = (v: bigint, start: number, count: number) =>
  Number((v >> BigInt(start)) & ((1n << BigInt(count)) - 1n));
const interval = (start: number, end: number) =>
  (-1n << BigInt(start)) & ((1n << BigInt(end)) - 1n);
function swapRows(rows: bigint[], i: number, j: number, mask: bigint): void {
  const delta = (rows[i]! ^ rows[j]!) & mask;
  rows[i] ^= delta;
  rows[j] ^= delta;
}
export function swapColumns(
  rows: bigint[],
  i: number,
  j: number,
  start = 0,
  end = rows.length
): void {
  const mask = (1n << BigInt(i)) | (1n << BigInt(j));
  for (let r = start; r < end; r++)
    if (((rows[r]! >> BigInt(i)) ^ (rows[r]! >> BigInt(j))) & 1n) rows[r] ^= mask;
}
/** ple_russian.c:_mzd_ple_russian. Split panels preserve the transformation
 * matrix, then up to seven Gray-code multiplication/elimination tables update
 * the remaining rows. Native cache-line copies become packed bit intervals.
 */
export function _mzd_ple_russian(a: mzd_t, k = 0): PLE {
  const rows = Array.from(a.rows),
    m = a.nrows,
    n = a.ncols,
    width = Math.ceil(n / 64);
  const P = Array.from({ length: m }, (_, i) => i),
    Q = Array.from({ length: n }, (_, i) => i);
  if (!m || !n) return { matrix: mzd_init(m, n, rows), rank: 0, P, Q };
  if (k === 0)
    k = Math.max(
      2,
      Math.min(
        8,
        Math.trunc(Math.log2(4194304 / 8 / width / 7)),
        Math.round(0.75 * Math.floor(Math.log2(Math.min(m, n))))
      )
    );
  if (!Number.isInteger(k) || k < 1 || k > 16) throw new RangeError('invalid PLE table parameter');
  let r = 0,
    c = 0,
    kk = 7 * k;
  while (c < n && r < m) {
    kk = Math.min(kk, n - c);
    const split =
      Math.min(Math.max(Math.floor((c + kk) / 64) + 1, Math.floor(c / 64) + 8), width) * 64;
    const low = (1n << BigInt(split)) - 1n,
      high = -1n << BigInt(split);
    const pivots: number[] = [],
      done: number[] = [];
    let rank = 0;
    // _mzd_ple_submatrix: lazy elimination of the bounded panel.
    for (let pos = 0; pos < kk; pos++) {
      let found = -1;
      for (let i = r + rank; i < m; i++) {
        if (nativeReadBits(rows[i]!, c, pos + 1)) {
          for (let l = 0; l < rank; l++)
            if (done[l]! < i) {
              if (bits(rows[i]!, c + pivots[l]!, 1))
                rows[i] ^= rows[r + l]! & interval(c + pivots[l]! + 1, split);
              done[l] = i;
            }
          if (bits(rows[i]!, c + pos, 1)) {
            found = i;
            break;
          }
        }
      }
      if (found >= 0) {
        P[r + rank] = found;
        swapRows(rows, found, r + rank, low);
        Q[r + rank] = c + pos;
        pivots.push(pos);
        done.push(found);
        rank++;
      }
    }
    const doneRow = rank < kk ? m - 1 : Math.max(...done);
    for (let l = 0; l < rank && c + pivots[l]! < n - 1; l++)
      for (let i = done[l]! + 1; i <= doneRow; i++)
        if (bits(rows[i]!, c + pivots[l]!, 1))
          rows[i] ^= rows[r + l]! & interval(c + pivots[l]! + 1, split);
    // _mzd_ple_a10: replay swaps and triangular updates on the right block.
    if (split < width * 64) {
      for (let i = r; i < r + rank; i++) swapRows(rows, i, P[i]!, high);
      for (let i = 1; i < rank; i++)
        for (let j = 0; j < i; j++)
          if (bits(rows[r + i]!, c + pivots[j]!, 1)) rows[r + i] ^= rows[r + j]! & high;
    }
    const wordStart = Math.floor(c / 64) * 64;
    const U = Array.from(
      { length: rank },
      (_, i) => rows[r + i]! & ~interval(wordStart, c + pivots[i]!)
    );
    if (rank === 0) {
      c += kk;
      let foundRow = -1,
        foundCol = n;
      for (let i = r; i < m; i++) {
        const tail = rows[i]! >> BigInt(c);
        if (tail) {
          const j = c + (tail & -tail).toString(2).length - 1;
          if (j < foundCol) {
            foundRow = i;
            foundCol = j;
          }
        }
      }
      if (foundRow < 0) break;
      P[r] = foundRow;
      Q[r] = foundCol;
      swapRows(rows, r, foundRow, -1n);
      for (let i = r + 1; i < m; i++)
        if (bits(rows[i]!, foundCol, 1)) rows[i] ^= rows[r]! & (-1n << BigInt(foundCol + 1));
      c = foundCol + 1;
      r++;
      continue;
    }
    let ntables = 1;
    for (let count = 7; count >= 2; count--)
      if (kk >= (count - 1) * k && kk >= count) {
        ntables = count;
        break;
      }
    const sizes = Array.from(
      { length: ntables },
      (_, i) =>
        Math.floor(kk / ntables) + Number(i !== ntables - 1 && kk % ntables >= ntables - 1 - i)
    );
    let column = 0,
      pivot = 0;
    const tables = sizes.map((size) => {
      const start = column,
        first = pivot;
      while (pivot < rank && pivots[pivot]! < start + size) pivot++;
      const count = pivot - first,
        table = Array<bigint>(2 ** count).fill(0n),
        M = new Map<number, number>([[0, 0]]),
        E = new Map<number, number>([[0, 0]]);
      const writeStart = Math.floor((c + start) / 64) * 64;
      for (let i = 1; i < table.length; i++) {
        const gray = i ^ (i >> 1),
          previous = (i - 1) ^ ((i - 1) >> 1),
          bit = 31 - Math.clz32(gray ^ previous);
        table[i] = (table[i - 1]! ^ U[first + bit]!) & (-1n << BigInt(writeStart));
        let spread = 0;
        for (let j = 0; j < count; j++)
          if ((gray >> j) & 1) spread += 2 ** (pivots[first + j]! - start);
        M.set(spread, i);
        if (rank === kk) E.set(bits(table[i]!, c + start, size), i);
      }
      if (rank === kk)
        for (let i = 1; i < table.length; i++)
          table[i] ^= BigInt(i ^ (i >> 1)) << BigInt(c + start);
      column += size;
      return { size, start, table, M, E };
    });
    // Already scanned panel rows need only their high block updated.
    if (split < width * 64)
      for (let i = r + rank; i <= doneRow; i++) {
        let sum = 0n;
        for (const T of tables) sum ^= T.table[T.M.get(bits(rows[i]!, c + T.start, T.size))!]!;
        rows[i] ^= sum & high;
      }
    // Unscanned rows use the elimination lookup, preserving the lower factor.
    for (let i = doneRow + 1; i < m; i++) {
      let value = nativeReadBits(rows[i]!, c, kk),
        sum = 0n;
      for (const T of tables) {
        const key = Number((value >> BigInt(T.start % 64)) & ((1n << BigInt(T.size)) - 1n));
        const term = T.table[T.E.get(key)!]!;
        sum ^= term;
        value ^= nativeReadBits(term, c, Math.min(64, n - c));
      }
      rows[i] ^= sum;
    }
    c += kk;
    r += rank;
  }
  // Compress L using the original forward transpositions.
  for (let j = 0; j < r; j++) if (Q[j]! > j) swapColumns(rows, Q[j]!, j, j, r);
  for (let j = 0; j < r; j++) swapColumns(rows, j, Q[j]!, r, m);
  return { matrix: mzd_init(m, n, rows), rank: r, P, Q };
}

/** ple_russian.c: Russian PLUQ adds the upper-triangle column permutation. */
export function _mzd_pluq_russian(a: mzd_t, k = 0): PLE {
  const result = _mzd_ple_russian(a, k),
    rows = Array.from(result.matrix.rows);
  for (let j = 0; j < a.ncols; j++) swapColumns(rows, j, result.Q[j]!, 0, Math.min(j, a.nrows));
  return { ...result, matrix: mzd_init(a.nrows, a.ncols, rows) };
}
