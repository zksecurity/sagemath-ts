import { random_zv } from './random.js';
/** PARI basemath/F2v.c column elimination and native kernel basis order.
 * Packed columns have the row count supplied separately.
 * @see Deviation: PARI packed binary-polynomial kernels
 */
export function F2m_ker_sp(columns: bigint[], rows: number, deplin = 0): bigint[] | bigint | null {
  if (!Number.isSafeInteger(rows) || rows < 0)
    throw new RangeError('row count must be nonnegative');
  const mask = (1n << BigInt(rows)) - 1n;
  if (columns.some((c) => c < 0n || (c & mask) !== c))
    throw new RangeError('column exceeds matrix row count');
  if (!columns.length) return columns;
  let unused = mask;
  const pivots: number[] = [];
  for (let k = 0; k < columns.length; k++) {
    const available = columns[k]! & unused;
    if (!available) {
      if (deplin) {
        let v = 1n << BigInt(k);
        for (let i = 0; i < k; i++)
          if ((columns[k]! >> BigInt(pivots[i]!)) & 1n) v |= 1n << BigInt(i);
        return v;
      }
      pivots.push(-1);
    } else {
      const lowest = available & -available,
        j = lowest.toString(2).length - 1;
      unused ^= lowest;
      pivots.push(j);
      columns[k]! ^= lowest;
      for (let i = k + 1; i < columns.length; i++)
        if (columns[i]! & lowest) columns[i]! ^= columns[k]!;
      columns[k]! |= lowest;
    }
  }
  if (deplin) return null;
  const out: bigint[] = [];
  for (let k = 0; k < columns.length; k++)
    if (pivots[k] === -1) {
      let v = 1n << BigInt(k);
      for (let i = 0; i < k; i++)
        if (pivots[i] !== -1 && (columns[k]! >> BigInt(pivots[i]!)) & 1n) v |= 1n << BigInt(i);
      out.push(v);
    }
  return out;
}
export function F2m_ker(columns: readonly bigint[], rows: number): bigint[] {
  return F2m_ker_sp(columns.slice(), rows) as bigint[];
}

/** PARI F2m_mul: XOR the left columns selected by each right column.
 * @see Deviation: PARI matrix product representation
 */
export function F2m_mul(x: readonly bigint[], y: readonly bigint[], rows: number): bigint[] {
  if (!Number.isSafeInteger(rows) || rows < 0)
    throw new RangeError('row count must be nonnegative');
  if (x.some((c) => c < 0n || c >> BigInt(rows)) || y.some((c) => c < 0n || c >> BigInt(x.length)))
    throw new RangeError('column exceeds matrix row count');
  return y.map((c) => {
    let result = 0n;
    for (let i = 0; i < x.length; i++) if ((c >> BigInt(i)) & 1n) result ^= x[i]!;
    return result;
  });
}

// F2v.c:671–1077, 64-bit block Lanczos; based on Jason Papadopoulos's
// GPLv2+ implementation as incorporated by PARI. Word arrays represent the
// 64 simultaneous vectors by one unsigned word per matrix coordinate.
const WORD_MASK = (1n << 64n) - 1n;
const identity = () => Array.from({ length: 64 }, (_, i) => 1n << BigInt(i));
function transpose(x: bigint[]): bigint[] {
  const z = Array<bigint>(64).fill(0n);
  for (let i = 0; i < x.length; i++)
    for (let j = 0; j < 64; j++) if (x[i]! & (1n << BigInt(j))) z[j]! |= 1n << BigInt(i);
  return z;
}
function mul(a: bigint[], b: bigint[]): bigint[] {
  return a.map((ai) => {
    let s = 0n;
    for (let j = 0; ai; j++, ai >>= 1n) if (ai & 1n) s ^= b[j]!;
    return s;
  });
}
function mulAdd(v: bigint[], x: bigint[], y: bigint[]): void {
  const c = Array<bigint>(2048).fill(0n);
  for (let j = 0; j < 8; j++)
    for (let i = 1; i < 256; i++) {
      const low = i & -i;
      c[j * 256 + i] = c[j * 256 + (i ^ low)]! ^ x[j * 8 + 31 - Math.clz32(low)]!;
    }
  for (let i = 0; i < y.length; i++) {
    let word = v[i]!;
    for (let j = 0; j < 8; j++, word >>= 8n) y[i]! ^= c[j * 256 + Number(word & 255n)]!;
  }
}
function transmul(x: bigint[], y: bigint[]): bigint[] {
  const c = Array<bigint>(2048).fill(0n),
    z = Array<bigint>(64).fill(0n);
  for (let i = 0; i < x.length; i++) {
    let xi = x[i]!;
    for (let j = 0; j < 8; j++, xi >>= 8n) c[j * 256 + Number(xi & 255n)]! ^= y[i]!;
  }
  for (let i = 0; i < 8; i++)
    for (let j = 0; j < 256; j++)
      if ((j >> i) & 1) for (let k = 0; k < 8; k++) z[k * 8 + i]! ^= c[k * 256 + j]!;
  return z;
}
function nonsingular(t: bigint[], last: number[]): [bigint[], number[]] | null {
  const M0 = t.slice(),
    M1 = identity(),
    s: number[] = [];
  const old = new Set(last),
    cols = Array.from({ length: 64 }, (_, i) => i)
      .filter((i) => !old.has(i))
      .concat(last.slice().reverse());
  for (let i = 0; i < 64; i++) {
    const ri = cols[i]!,
      mask = 1n << BigInt(ri);
    let j = i;
    for (; j < 64; j++) if (M0[cols[j]!]! & mask) break;
    const success = j < 64;
    if (!success) {
      for (j = i; j < 64; j++) if (M1[cols[j]!]! & mask) break;
      if (j === 64) return null;
    }
    const rj = cols[j]!;
    [M0[rj], M0[ri]] = [M0[ri]!, M0[rj]!];
    [M1[rj], M1[ri]] = [M1[ri]!, M1[rj]!];
    for (const row of cols)
      if (row !== ri && (success ? M0 : M1)[row]! & mask) {
        M0[row]! ^= M0[ri]!;
        M1[row]! ^= M1[ri]!;
      }
    if (success) s.push(ri);
    else M0[ri] = M1[ri] = 0n;
  }
  if (new Set([...s, ...last]).size !== 64) return null;
  return [M1, s];
}
function sparseTransmul(x: bigint[], A: readonly (readonly number[])[], rows: number): bigint[] {
  const b = Array<bigint>(rows).fill(0n);
  for (let i = 0; i < A.length; i++) for (const j of A[i]!) b[j - 1]! ^= x[i]!;
  return b;
}
function sparseMul(x: bigint[], A: readonly (readonly number[])[]): bigint[] {
  return A.map((c) => {
    let s = 0n;
    for (const j of c) s ^= x[j - 1]!;
    return s;
  });
}
function addid(f: bigint[]): void {
  for (let i = 0; i < 64; i++) f[i]! ^= 1n << BigInt(i);
}
function mask(f: bigint[], m: bigint): void {
  for (let i = 0; i < f.length; i++) f[i]! &= m;
}
function independent(x: bigint[], rows: number): number[] {
  x = x.slice();
  let unused = (1n << BigInt(rows)) - 1n;
  const p: number[] = [];
  for (let k = 0; k < x.length; k++) {
    const a = x[k]! & unused;
    if (!a) continue;
    const low = a & -a;
    unused ^= low;
    p.push(k);
    for (let i = k + 1; i < x.length; i++) if (x[i]! & low) x[i]! ^= x[k]!;
  }
  return p;
}
function block_lanczos(B: readonly (readonly number[])[], rows: number): bigint[] | null {
  const n = B.length,
    zero = () => Array<bigint>(64).fill(0n);
  let v1 = Array<bigint>(n).fill(0n),
    v2 = v1.slice(),
    vt1 = zero(),
    vt21 = zero(),
    winv1 = zero(),
    winv2 = zero();
  let s = Array.from({ length: 64 }, (_, i) => i),
    mask1 = WORD_MASK;
  let x = random_zv(n).map((v) => BigInt.asUintN(64, v));
  const w = sparseMul(sparseTransmul(x, B, rows), B);
  let v0 = w;
  for (;;) {
    const next = sparseMul(sparseTransmul(v0, B, rows), B),
      vt0 = transmul(v0, next);
    if (vt0.every((v) => !v)) break;
    const vt20 = transmul(next, next),
      inv = nonsingular(vt0, s);
    if (!inv) return null;
    const winv0 = inv[0];
    s = inv[1];
    let mask0 = 0n;
    for (const i of s) mask0 |= 1n << BigInt(i);
    let d = mul(
      winv0,
      vt20.map((v, i) => (v & mask0) ^ vt0[i]!)
    );
    addid(d);
    const e = mul(winv1, vt0);
    mask(e, mask0);
    let f = mul(vt1, winv1);
    addid(f);
    f = mul(winv2, f);
    f = mul(
      f,
      vt21.map((v, i) => ((v & mask1) ^ vt1[i]!) & mask0)
    );
    mask(next, mask0);
    mulAdd(v0, d, next);
    mulAdd(v1, e, next);
    mulAdd(v2, f, next);
    d = mul(winv0, transmul(v0, w));
    mulAdd(v0, d, x);
    v2 = v1;
    v1 = v0;
    v0 = next;
    winv2 = winv1;
    winv1 = winv0;
    vt1 = vt0;
    vt21 = vt20;
    mask1 = mask0;
  }
  const image = transpose(sparseTransmul(x, B, rows)).concat(
    transpose(sparseTransmul(v0, B, rows))
  );
  x = transpose(x).concat(transpose(v0));
  const p = independent(x, n);
  return F2m_mul(
    p.map((i) => x[i]!),
    F2m_ker(
      p.map((i) => image[i]!),
      rows
    ),
    n
  );
}
/** Native alglin1.c:4851 singleton-column elimination. Sparse row indices and
 * the returned column indices are one-based; outer arrays are zero-based.
 * @see Deviation: PARI sparse binary kernel representation
 */
export function F2Ms_colelim(M: readonly (readonly number[])[], rows: number): number[] {
  const active = Array<boolean>(M.length).fill(true),
    weight = Array<number>(rows).fill(0);
  for (const c of M) for (const j of c) weight[j - 1]!++;
  let changed: boolean;
  do {
    changed = false;
    for (let i = 0; i < M.length; i++)
      if (active[i] && M[i]!.some((j) => weight[j - 1] === 1)) {
        active[i] = false;
        changed = true;
        for (const j of M[i]!) weight[j - 1]!--;
      }
  } while (changed);
  return M.map((_, i) => i + 1).filter((i) => active[i - 1]);
}
/** Native F2v.c:1063 sparse kernel: dense through 640 rows, block Lanczos above.
 * Results are zero-based BigInt column bitsets. The sparse path consumes the
 * native global random stream and may return a proper subspace of the kernel.
 * @see Deviation: PARI sparse binary kernel representation
 */
export function F2Ms_ker(M: readonly (readonly number[])[], rows: number): bigint[] {
  if (rows <= 640)
    return F2m_ker(
      M.map((c) => c.reduce((v, j) => v | (1n << BigInt(j - 1)), 0n)),
      rows
    );
  const p = F2Ms_colelim(M, rows),
    B = p.map((i) => M[i - 1]!);
  let R: bigint[] | null;
  do {
    R = block_lanczos(B, rows);
  } while (!R);
  return R.map((v) => {
    let w = 0n;
    for (let i = 0; i < p.length; i++) if (v & (1n << BigInt(i))) w |= 1n << BigInt(p[i]! - 1);
    return w;
  });
}

/** Native F2m_gauss_pivot; packed columns, one-based pivot rows and nullity.
 * @see Deviation: PARI prime-decomposition matrix adapters
 */
export function F2m_gauss_pivot(columns: bigint[], rows: number): [number[] | null, number] {
  if (!columns.length) return [null, 0];
  const x = columns.slice(),
    d: number[] = [],
    mask = (1n << BigInt(rows)) - 1n;
  let unused = mask,
    r = 0;
  for (let k = 0; k < x.length; k++) {
    const v = x[k]! & unused;
    if (!v) {
      r++;
      d.push(0);
      continue;
    }
    const bit = v & -v,
      j = bit.toString(2).length;
    unused ^= bit;
    d.push(j);
    for (let i = k + 1; i < x.length; i++) if (x[i]! & bit) x[i]! ^= x[k]!;
  }
  return [d, r];
}
/** Native F2m_gauss_sp, with exact column-wise clearing and pivot preference.
 * @see Deviation: PARI prime-decomposition matrix adapters
 */
export function F2m_gauss(A: bigint[], B: bigint[], rows: number): bigint[] | null {
  if (!A.length) return [];
  const a = A.slice(),
    b = B.slice(),
    d = Array<number>(rows).fill(-1);
  for (let i = 0; i < a.length; i++) {
    let ai = a[i]!,
      k = i;
    if (k >= rows || d[k] !== -1 || !(ai & (1n << BigInt(k)))) {
      for (k = 0; k < rows; k++) if (d[k] === -1 && ai & (1n << BigInt(k))) break;
    }
    if (k === rows) return null;
    d[k] = i;
    const bit = 1n << BigInt(k);
    ai &= ~bit;
    for (let j = 0; j < a.length; j++) if (a[j]! & bit) a[j]! ^= ai;
    for (let j = 0; j < b.length; j++) if (b[j]! & bit) b[j]! ^= ai;
  }
  return b.map((c) =>
    d.reduce((v, j, i) => (j >= 0 && c & (1n << BigInt(i)) ? v | (1n << BigInt(j)) : v), 0n)
  );
}
