import { type mzd_t, mzd_init, mzd_mul_naive, _mzd_density } from './mzd.js';
/** M4RI's Gray-code combination table and inverse indexing. */
export function mzd_make_table(
  a: mzd_t,
  r: number,
  k: number
): { rows: bigint[]; lookup: number[] } {
  if (!Number.isInteger(k) || k < 0 || k > 16 || !Number.isInteger(r) || r < 0 || r + k > a.nrows)
    throw new RangeError('invalid M4RI table rows');
  const rows = Array<bigint>(2 ** k).fill(0n),
    lookup = Array<number>(2 ** k).fill(0);
  for (let i = 1; i < rows.length; i++) {
    const gray = i ^ (i >> 1),
      previous = (i - 1) ^ ((i - 1) >> 1),
      bit = 31 - Math.clz32(gray ^ previous);
    rows[i] = rows[i - 1]! ^ a.rows[r + bit]!;
    lookup[gray] = i;
  }
  return { rows, lookup };
}
/** M4RI multiplication with eight Gray-code tables and original small-input dispatch.
 * A portable 4 MiB cache tuning model matches the pinned oracle build; scheduling
 * and cache detection do not affect GF(2) results.
 */
export function mzd_mul_m4rm(a: mzd_t, b: mzd_t, k = 0): mzd_t {
  if (a.ncols !== b.nrows) throw new RangeError('incompatible M4RI multiplication dimensions');
  if (!Number.isInteger(k)) throw new RangeError('M4RI table parameter must be an integer');
  if (b.ncols < 54 || a.nrows < 16) return mzd_mul_naive(a, b);
  if (a.ncols === 0) return mzd_init(a.nrows, b.ncols);
  if (k === 0) {
    const width = Math.ceil(b.ncols / 64),
      cache = 4194304;
    k = Math.trunc(Math.log2(cache / 64 / width));
    if (cache - 64 * 2 ** k * width > 64 * 2 ** (k + 1) * width - cache) k++;
    k = Math.min(k, Math.round(0.75 * Math.floor(Math.log2(Math.min(a.nrows, a.ncols, b.ncols)))));
  }
  k = Math.max(2, Math.min(8, k));
  const rows = Array<bigint>(a.nrows).fill(0n),
    kk = 8 * k,
    end = Math.floor(a.ncols / kk),
    mask = (1n << BigInt(k)) - 1n;
  for (let giant = 0; giant < a.nrows; giant += 2048) {
    for (let i = 0; i < end; i++) {
      const tables = Array.from({ length: 8 }, (_, z) => mzd_make_table(b, kk * i + k * z, k));
      for (let j = giant; j < Math.min(giant + 2048, a.nrows); j++) {
        let bits = a.rows[j]! >> BigInt(kk * i);
        for (let z = 0; z < 8; z++) {
          const table = tables[z]!;
          rows[j] ^= table.rows[table.lookup[Number(bits & mask)]!]!;
          bits >>= BigInt(k);
        }
      }
    }
  }
  for (let offset = kk * end; offset < a.ncols; offset += k) {
    const size = Math.min(k, a.ncols - offset),
      table = mzd_make_table(b, offset, size),
      lastMask = (1n << BigInt(size)) - 1n;
    for (let j = 0; j < a.nrows; j++)
      rows[j] ^= table.rows[table.lookup[Number((a.rows[j]! >> BigInt(offset)) & lastMask)]!]!;
  }
  return { nrows: a.nrows, ncols: b.ncols, rows };
}

/** Native panel reduction. Rows outside the panel are cleared lazily while
 * searching for pivots, then updated with up to six Gray-code tables.
 */
function gaussPanel(
  rows: bigint[],
  r: number,
  c: number,
  end: number,
  k: number,
  full: boolean
): number {
  let foundRows = 0;
  for (let j = c; j < c + k; j++) {
    let found = false;
    for (let i = r + foundRows; i < end; i++) {
      const original = nativeReadBits(rows[i]!, c, j - c + 1);
      for (let l = 0; l < j - c; l++) {
        const bit = full ? (original >> BigInt(l % 64)) & 1n : (rows[i]! >> BigInt(c + l)) & 1n;
        if (bit) rows[i] ^= rows[r + l]! & (-1n << BigInt(c + l));
      }
      if ((rows[i]! >> BigInt(j)) & 1n) {
        const pivot = r + foundRows;
        [rows[i], rows[pivot]] = [rows[pivot]!, rows[i]!];
        if (full)
          for (let l = r; l < pivot; l++)
            if ((rows[l]! >> BigInt(j)) & 1n) rows[l] ^= rows[pivot]! & (-1n << BigInt(j));
        foundRows++;
        found = true;
        break;
      }
    }
    if (!found) break;
  }
  return foundRows;
}
function topPanel(rows: bigint[], r: number, c: number, k: number): void {
  for (let j = 0; j < k; j++)
    for (let i = r; i < r + j; i++)
      if ((rows[i]! >> BigInt(c + j)) & 1n) rows[i] ^= rows[r + j]! & (-1n << BigInt(c + j));
}
function panelSizes(size: number, k: number): number[] {
  const n = Math.min(6, Math.ceil(size / k));
  if (n === 2) return [Math.floor(size / 2), Math.ceil(size / 2)];
  return Array.from(
    { length: n },
    (_, i) => Math.floor(size / n) + Number(size % n >= n - 1 - i && i !== n - 1)
  );
}
function clearPanel(
  rows: bigint[],
  ncols: number,
  r: number,
  c: number,
  size: number,
  k: number,
  start: number,
  end: number
): void {
  if (!size || start >= end) return;
  let offset = 0;
  const tables = panelSizes(size, k).map((width) => {
    const a = { nrows: rows.length, ncols, rows },
      table = mzd_make_table(a, r + offset, width);
    const shift = BigInt(c + offset),
      mask = (1n << BigInt(width)) - 1n;
    offset += width;
    return { table, shift, mask };
  });
  const tail = -1n << BigInt(c);
  for (let i = start; i < end; i++) {
    const original = nativeReadBits(rows[i]!, c, size) << BigInt(c);
    let sum = 0n;
    for (const { table, shift, mask } of tables)
      sum ^= table.rows[table.lookup[Number((original >> shift) & mask)]!]!;
    rows[i] ^= sum & tail;
  }
}
/** mzd_find_pivot: first nonzero column, earliest row at that column. */
function nextPivot(rows: bigint[], r: number, c: number): [number, number] | undefined {
  let best: bigint | undefined,
    row = 0;
  for (let i = r; i < rows.length; i++) {
    const v = rows[i]! >> BigInt(c);
    if (v) {
      const low = v & -v;
      if (best === undefined || low < best) {
        best = low;
        row = i;
        if (low === 1n) break;
      }
    }
  }
  return best === undefined ? undefined : [row, c + best.toString(2).length - 1];
}
/** brilliantrussian.c:_mzd_echelonize_m4ri without heuristic PLUQ crossover.
 * An owned result replaces native in-place mutation. Automatic tuning uses the
 * same portable 4 MiB cache model as multiplication.
 */
export function _mzd_echelonize_m4ri(
  a: mzd_t,
  full: boolean,
  k = 0,
  crossover?: (a: mzd_t, full: boolean) => { matrix: mzd_t; rank: number }
): { matrix: mzd_t; rank: number } {
  if (!Number.isInteger(k) || k < 0 || k > 16)
    throw new RangeError('invalid M4RI echelon table parameter');
  const rows = Array.from(a.rows),
    n = a.ncols;
  if (!a.nrows || !n) return { matrix: mzd_init(a.nrows, n, rows), rank: 0 };
  if (k === 0) {
    k = Math.min(
      7,
      Math.max(1, Math.trunc(0.75 * (1 + Math.floor(Math.log2(Math.min(a.nrows, n))))))
    );
    if (0.75 * 2 ** k * n > 4194304 / 2) k--;
  }
  let kk = 6 * k,
    r = 0,
    c = 0,
    lastCheck = 0;
  while (c < n) {
    if (crossover && r < rows.length && (c === 0 || c > lastCheck + 256)) {
      lastCheck = c;
      if (_mzd_density({ nrows: a.nrows, ncols: n, rows }, 32, r, c) >= 0.15) {
        const start = Math.floor(c / 64) * 64;
        const result = crossover(
          mzd_init(
            rows.length - r,
            n - start,
            rows.slice(r).map((v) => v >> BigInt(start))
          ),
          full
        );
        const lowerMask = (1n << BigInt(start)) - 1n;
        for (let i = r; i < rows.length; i++)
          rows[i] = (rows[i]! & lowerMask) | (result.matrix.rows[i - r]! << BigInt(start));
        if (full && r > 0) topEchelon(rows, n, r, c, r);
        r += result.rank;
        break;
      }
    }
    kk = Math.min(kk, n - c);
    const count = gaussPanel(rows, r, c, rows.length, kk, full);
    const saved = full ? [] : rows.slice(r, r + count);
    if (!full) topPanel(rows, r, c, count);
    if (count === kk) clearPanel(rows, n, r, c, count, k, r + count, rows.length);
    if (full) clearPanel(rows, n, r, c, count, k, 0, r);
    if (!full) for (let i = 0; i < count; i++) rows[r + i] = saved[i]!;
    r += count;
    c += count;
    if (kk !== count) {
      const pivot = nextPivot(rows, r, c);
      if (!pivot) break;
      const [i, j] = pivot;
      c = j;
      [rows[r], rows[i]] = [rows[i]!, rows[r]!];
    }
  }
  return { matrix: mzd_init(a.nrows, n, rows), rank: r };
}

/** Original top-echelon pass after heuristic PLUQ crossover. */
function topEchelon(rows: bigint[], n: number, r: number, c: number, max: number): void {
  let k = Math.min(
    7,
    Math.max(1, Math.trunc(0.75 * (1 + Math.floor(Math.log2(Math.min(max, n))))))
  );
  if (0.75 * 2 ** k * n > 4194304 / 2) k--;
  let kk = 6 * k;
  while (c < n) {
    kk = Math.min(kk, n - c);
    const count = gaussPanel(rows, r, c, Math.min(rows.length, r + kk), kk, true);
    clearPanel(rows, n, r, c, count, k, 0, Math.min(r, max));
    r += count;
    c += count;
    if (kk !== count) c++;
  }
}

/** mzd.h:mzd_read_bits. The pinned release disables assertions; oversized
 * panels retain its 64-bit shift behavior on the supported native oracle.
 */
export function nativeReadBits(row: bigint, c: number, n: number): bigint {
  const spot = c % 64,
    block = Math.floor(c / 64),
    spill = spot + n - 64,
    mask = (1n << 64n) - 1n;
  const current = (row >> BigInt(block * 64)) & mask;
  const value =
    spill <= 0
      ? current << BigInt(-spill & 63)
      : (((row >> BigInt((block + 1) * 64)) & mask) << BigInt((64 - spill) & 63)) |
        (current >> BigInt(spill & 63));
  return (value & mask) >> BigInt((64 - n) & 63);
}

/** brilliantrussian.c:mzd_inv_m4ri. Align the identity window to a native word,
 * fully reduce the augmented matrix, and copy that window. The native k argument
 * is ignored. Like the native kernel, this does not check nonsingularity.
 * @see Deviation: Native Binary Inverse Adapter
 */
export function mzd_inv_m4ri(a: mzd_t, _k = 0): mzd_t {
  if (a.nrows !== a.ncols) throw new RangeError('M4RI inverse requires a square matrix');
  const n = a.nrows,
    nr = 64 * Math.ceil(n / 64);
  const augmented = mzd_init(
    n,
    2 * nr,
    a.rows.map((row, i) => row | (1n << BigInt(nr + i)))
  );
  const reduced = _mzd_echelonize_m4ri(augmented, true, 0).matrix;
  return mzd_init(
    n,
    n,
    reduced.rows.map((row) => row >> BigInt(nr))
  );
}
