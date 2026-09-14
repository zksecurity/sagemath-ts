/** PARI F3v.c packed ternary kernel; one BigInt stores the native two-bit trits. */
export function F3m_ker(columns: readonly bigint[], rows: number): bigint[] {
  const x = columns.slice(),
    n = x.length,
    d = Array<number>(n).fill(-1);
  let unused = (1n << BigInt(rows)) - 1n;
  const mask = ((1n << BigInt(2 * rows)) - 1n) / 3n;
  const swap = (a: bigint) => ((a & mask) << 1n) | ((a >> 1n) & mask);
  const plus = (a: bigint, b: bigint) => {
    const c = a ^ b ^ swap(a & b);
    return c & ~swap(c);
  };
  const coeff = (a: bigint, j: number) => (a >> BigInt(2 * j)) & 3n;
  for (let k = 0; k < n; k++) {
    let j = 0;
    for (; j < rows; j++) if ((unused >> BigInt(j)) & 1n && coeff(x[k]!, j)) break;
    if (j === rows) continue;
    const pos = BigInt(2 * j),
      pivot = coeff(x[k]!, j);
    x[k] = x[k]! & ~(3n << pos);
    unused &= ~(1n << BigInt(j));
    d[k] = j;
    for (let i = k + 1; i < n; i++) {
      const u = coeff(x[i]!, j);
      if (u) x[i] = plus(x[i]!, u === pivot ? swap(x[k]!) : x[k]!);
    }
    x[k] = x[k]! | (2n << pos);
    if (pivot === 1n)
      for (let i = k + 1; i < n; i++) if (coeff(x[i]!, j)) x[i] = x[i]! ^ (3n << pos);
  }
  const out: bigint[] = [];
  for (let k = 0; k < n; k++)
    if (d[k] === -1) {
      let v = 1n << BigInt(2 * k);
      for (let i = 0; i < k; i++) if (d[i] !== -1) v |= coeff(x[k]!, d[i]!) << BigInt(2 * i);
      out.push(v);
    }
  return out;
}
