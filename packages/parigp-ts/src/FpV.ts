/** PARI FpV.c matrix products. Columns and rows have unused slot zero.
 * @see Deviation: PARI matrix product representation
 */
import { F2m_mul } from './F2v.js';
import { ZM_mul } from './ZV.js';
import { residue } from './_polynomial_division.js';
import {
  classical,
  fromColumns,
  toColumns,
  winograd,
  type Rows,
} from './_matrix_mul.js';

/** Native word-modulus multiplication with 64-bit Strassen–Winograd thresholds. */
export function Flm_mul(x: bigint[][], y: bigint[][], p: bigint): bigint[][] {
  if (p < 1n || p >= 1n << 64n)
    throw new RangeError('modulus must be a positive word integer');
  const e = p.toString(2).length - 1,
    cutoff = e <= 29 ? 140 : e <= 62 ? 40 : 70;
  const mul = (a: Rows, b: Rows, m: number, n: number, k: number): Rows =>
    Math.min(m + 1, n + 1, k + 1) <= cutoff
      ? classical(a, b, m, n, k, p)
      : winograd(a, b, m, n, k, mul, p);
  const a = fromColumns(x),
    b = fromColumns(y),
    k = Math.max(0, y.length - 1);
  return toColumns(mul(a, b, a.length, Math.max(0, x.length - 1), k), k);
}

/** Native FpM multiplication dispatches binary, word and arbitrary-integer backends.
 * @see Deviation: PARI signed composition and matrix boundaries
 */
export function FpM_mul(x: bigint[][], y: bigint[][], p: bigint): bigint[][] {
  if (y.length <= 1) return [[]];
  if (x.length <= 1) return [[], ...y.slice(1).map(() => [0n])];
  const pp = p < 0n ? -p : p;
  const red = (v: bigint[][]) => [
    [],
    ...v.slice(1).map((c) => [0n, ...c.slice(1).map((z) => residue(z, pp))]),
  ];
  if (pp === 2n) {
    const pack = (v: bigint[][]) =>
      v.slice(1).map((c) => c.slice(1).reduce((z, a, i) => z | ((a & 1n) << BigInt(i)), 0n));
    const rows = x[1]!.length - 1,
      result = F2m_mul(pack(x), pack(y), rows);
    return [
      [],
      ...result.map((c) => [0n, ...Array.from({ length: rows }, (_, i) => (c >> BigInt(i)) & 1n)]),
    ];
  }
  return pp > 0n && pp < 1n << 64n ? Flm_mul(red(x), red(y), pp) : red(ZM_mul(x, y));
}
