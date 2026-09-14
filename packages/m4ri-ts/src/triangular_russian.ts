import { type mzd_t, mzd_init } from './mzd.js';
import { mzd_make_table } from './brilliantrussian.js';
/** triangular_russian.c: eight-table forward/back substitution. */
export function trsmRussian(a: mzd_t, b: mzd_t, upper: boolean): mzd_t {
  const rows = Array.from(b.rows),
    m = b.nrows,
    n = b.ncols;
  if (!m || !n) return mzd_init(m, n, rows);
  const k = Math.max(
    2,
    Math.min(
      8,
      Math.trunc(Math.log2(4194304 / 8 / Math.ceil(n / 64) / 8)),
      Math.round(0.75 * Math.floor(Math.log2(Math.min(m, n))))
    )
  );
  const panel = (start: number, count: number) => {
    for (let z = 0; z < count; z++) {
      const i = upper ? start + count - z - 1 : start + z;
      for (let v = 0; v < z; v++) {
        const j = upper ? i + v + 1 : start + v;
        if ((a.rows[i]! >> BigInt(j)) & 1n) rows[i] ^= rows[j]!;
      }
    }
  };
  const update = (start: number, count: number, tableSize: number) => {
    const tables = Array.from({ length: count / tableSize }, (_, z) => ({
      offset: start + z * tableSize,
      table: mzd_make_table({ nrows: m, ncols: n, rows }, start + z * tableSize, tableSize),
    }));
    const mask = (1n << BigInt(tableSize)) - 1n;
    for (let i = upper ? 0 : start + count; i < (upper ? start : m); i++) {
      let sum = 0n;
      for (const { offset, table } of tables)
        sum ^= table.rows[table.lookup[Number((a.rows[i]! >> BigInt(offset)) & mask)]!]!;
      rows[i] ^= sum;
    }
  };
  let done = 0;
  for (; done < m - 8 * k; done += 8 * k) {
    const start = upper ? m - done - 8 * k : done;
    panel(start, 8 * k);
    update(start, 8 * k, k);
  }
  while (done < m) {
    const size = Math.min(k, m - done),
      start = upper ? m - done - size : done;
    panel(start, size);
    update(start, size, size);
    done += size;
  }
  return mzd_init(m, n, rows);
}
