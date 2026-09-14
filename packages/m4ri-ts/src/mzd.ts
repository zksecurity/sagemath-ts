/** Packed row adapter for M4RI mzd_t. Bit j is column j; padding bits are zero.
 * Source: M4RI 20251207 mzd.c; native pointers/windows become immutable values.
 */
export interface mzd_t {
  readonly nrows: number;
  readonly ncols: number;
  readonly rows: readonly bigint[];
}
export function mzd_init(nrows: number, ncols: number, rows?: readonly bigint[]): mzd_t {
  if (!Number.isSafeInteger(nrows) || !Number.isSafeInteger(ncols) || nrows < 0 || ncols < 0)
    throw new RangeError('invalid M4RI matrix dimensions');
  if (rows !== undefined && rows.length !== nrows) throw new RangeError('incorrect M4RI row count');
  const mask = (1n << BigInt(ncols)) - 1n;
  return {
    nrows,
    ncols,
    rows: rows === undefined ? Array<bigint>(nrows).fill(0n) : rows.map((r) => r & mask),
  };
}
export function mzd_add(a: mzd_t, b: mzd_t): mzd_t {
  if (a.nrows !== b.nrows || a.ncols !== b.ncols)
    throw new RangeError('incompatible M4RI addition dimensions');
  return { nrows: a.nrows, ncols: a.ncols, rows: a.rows.map((r, i) => r ^ b.rows[i]!) };
}
export function mzd_init_window(a: mzd_t, r0: number, c0: number, r1: number, c1: number): mzd_t {
  if (
    ![r0, c0, r1, c1].every(Number.isSafeInteger) ||
    r0 < 0 ||
    r1 < r0 ||
    r0 > a.nrows ||
    c0 < 0 ||
    c1 < c0 ||
    c1 > a.ncols ||
    c0 % 64 !== 0
  )
    throw new RangeError('invalid M4RI window');
  const rows = a.rows.slice(r0, Math.min(r1, a.nrows));
  return mzd_init(
    rows.length,
    c1 - c0,
    rows.map((r) => r >> BigInt(c0))
  );
}
/** mzd.c:mzd_submatrix: copy aligned words or extract shifted bit windows.
 * The adapter owns its result and applies Sage's zero-dimensional shortcut.
 */
export function mzd_submatrix(a: mzd_t, r0: number, c0: number, r1: number, c1: number): mzd_t {
  if (
    ![r0, c0, r1, c1].every(Number.isSafeInteger) ||
    r0 < 0 ||
    c0 < 0 ||
    r1 < r0 ||
    c1 < c0 ||
    r1 > a.nrows ||
    c1 > a.ncols
  )
    throw new RangeError('invalid M4RI submatrix');
  return mzd_init(
    r1 - r0,
    c1 - c0,
    a.rows.slice(r0, r1).map((row) => row >> BigInt(c0))
  );
}
export function mzd_mul_naive(a: mzd_t, b: mzd_t): mzd_t {
  if (a.ncols !== b.nrows) throw new RangeError('incompatible M4RI multiplication dimensions');
  const rows = Array<bigint>(a.nrows).fill(0n);
  if (b.ncols < 54) {
    // mzd_mul_naive transposes narrow right operands and takes packed dot parities.
    const transposed = Array<bigint>(b.ncols).fill(0n);
    for (let i = 0; i < b.nrows; i++)
      for (let j = 0; j < b.ncols; j++)
        if ((b.rows[i]! >> BigInt(j)) & 1n) transposed[j] |= 1n << BigInt(i);
    let fold = 1;
    while (fold < a.ncols) fold *= 2;
    for (let i = 0; i < a.nrows; i++)
      for (let j = 0; j < b.ncols; j++) {
        let dot = a.rows[i]! & transposed[j]!;
        for (let shift = fold / 2; shift >= 1; shift /= 2) dot ^= dot >> BigInt(shift);
        rows[i] |= (dot & 1n) << BigInt(j);
      }
  } else {
    // _mzd_mul_va combines whole rows selected by each left row.
    for (let i = 0; i < a.nrows; i++)
      for (let j = 0; j < a.ncols; j++) if ((a.rows[i]! >> BigInt(j)) & 1n) rows[i] ^= b.rows[j]!;
  }
  return { nrows: a.nrows, ncols: b.ncols, rows };
}

/** mzd.c:m4ri_bitcount's six parallel mask-and-add steps. */
function m4ri_bitcount(word: bigint): bigint {
  let n = BigInt.asUintN(64, word);
  for (let c = 0; c < 6; c++) {
    const shift = 1n << BigInt(c),
      mask = ((1n << 64n) - 1n) / ((1n << shift) + 1n);
    n = (n & mask) + ((n >> shift) & mask);
  }
  return n;
}
/** mzd.c:_mzd_density/mzd_density, including the original full-word tail rule.
 * Native zero-width row dereferences become RuntimeError instead of terminating JS.
 * @see Deviation: Binary Matrix Density
 */
export function mzd_density(a: mzd_t, res = 0): number {
  return _mzd_density(a, res, 0, 0);
}
export function _mzd_density(a: mzd_t, res: number, r: number, c: number): number {
  if (!Number.isSafeInteger(res)) throw new RangeError('invalid M4RI density resolution');
  const width = Math.ceil(a.ncols / 64);
  if (a.nrows > 0 && width === 0)
    throw Object.assign(new Error('Segmentation fault'), { name: 'RuntimeError' });
  let count = 0n,
    total = 0n;
  if (width === 1) {
    for (const row of a.rows.slice(r)) count += m4ri_bitcount(row & (-1n << BigInt(c)));
    return Number(count) / (1.0 * a.ncols * a.nrows);
  }
  if (res === 0) res = Math.floor(width / 100);
  if (res < 1) res = 1;
  for (const row of a.rows.slice(r)) {
    count += c < 64 ? m4ri_bitcount(row & (-1n << BigInt(c))) : 0n;
    total += 64n;
    for (let j = Math.max(1, Math.floor(c / 64)); j < width - 1; j += res) {
      count += m4ri_bitcount(row >> BigInt(64 * j));
      total += 64n;
    }
    const tail = a.ncols % 64;
    count += m4ri_bitcount(
      (row >> BigInt(64 * Math.floor(a.ncols / 64))) & ((1n << BigInt(tail)) - 1n)
    );
    total += BigInt(tail);
  }
  return Number(count) / Number(total);
}

/** mzd.c:mzd_concat; copy the left words and append the right bit interval.
 * @see Deviation: Binary Basic Operations and Column Cache
 */
export function mzd_concat(a: mzd_t, b: mzd_t): mzd_t {
  if (a.nrows !== b.nrows) throw new RangeError('incompatible M4RI concatenation dimensions');
  return mzd_init(
    a.nrows,
    a.ncols + b.ncols,
    a.rows.map((row, i) => row | (b.rows[i]! << BigInt(a.ncols)))
  );
}

/** mzd.c:mzd_transpose. Recursive cache blocks and native six-stage bit swaps.
 * Owned packed rows replace native destination pointers and windows.
 * @see Deviation: Binary Basic Operations and Column Cache
 */
export function mzd_transpose(a: mzd_t): mzd_t {
  if (!a.nrows || !a.ncols) {
    if (a.nrows !== a.ncols) throw Object.assign(new Error('Aborted'), { name: 'RuntimeError' });
    return mzd_init(a.ncols, a.nrows);
  }
  const rows = Array<bigint>(a.ncols).fill(0n);
  const masks = [
    0x5555555555555555n,
    0x3333333333333333n,
    0x0f0f0f0f0f0f0f0fn,
    0x00ff00ff00ff00ffn,
    0x0000ffff0000ffffn,
    0x00000000ffffffffn,
  ];
  const block = (r: number, c: number, m: number, n: number): void => {
    const size = Math.max(m, n);
    if (size > 512) {
      const align = size <= 768 ? 64 : 512,
        split = Math.ceil(Math.floor(size / 2) / align) * align;
      if (m >= n) {
        block(r, c, split, n);
        block(r + split, c, m - split, n);
      } else {
        block(r, c, m, split);
        block(r, c + split, m, n - split);
      }
      return;
    }
    for (let i = 0; i < m; i += 64)
      for (let j = 0; j < n; j += 64) {
        const h = Math.min(64, m - i),
          w = Math.min(64, n - j),
          mask = (1n << BigInt(w)) - 1n;
        const t = Array.from({ length: 64 }, (_, k) =>
          k < h ? (a.rows[r + i + k]! >> BigInt(c + j)) & mask : 0n
        );
        for (let stage = 0; stage < 6; stage++) {
          const step = 2 ** stage,
            shift = BigInt(step),
            swapMask = masks[stage]!;
          for (let start = 0; start < 64; start += 2 * step)
            for (let k = start; k < start + step; k++) {
              const xor = ((t[k]! >> shift) ^ t[k + step]!) & swapMask;
              t[k] ^= xor << shift;
              t[k + step] ^= xor;
            }
        }
        for (let k = 0; k < w; k++) rows[c + j + k] |= t[k]! << BigInt(r + i);
      }
  };
  block(0, 0, a.nrows, a.ncols);
  return mzd_init(a.ncols, a.nrows, rows);
}
