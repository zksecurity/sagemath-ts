/** PARI lll.c: QR/Cholesky precision selection for FLATTER.
 * @see Deviation: PARI adaptive Gram-Schmidt dependencies
 */
import { QR_init, type QrScalar } from './bibli1.js';
import { RgM_Cholesky } from './alglin2.js';
import { type MatrixReal } from './_real_matrix.js';
import { itor, type MpReal } from './qfb.js';
import { abscmprr } from './kernel/none/cmp.js';

type GsoMatrix = QrScalar[][];
const highExpo = 1n << 61n;
const abs = (x: bigint): bigint => (x < 0n ? -x : x);
const expo = (x: QrScalar): bigint =>
  typeof x === 'bigint' ? (x ? BigInt(abs(x).toString(2).length - 1) : -highExpo) : BigInt(x.e);
const max = (x: bigint, y: bigint): bigint => (x > y ? x : y);
function shape(M: GsoMatrix): void {
  if (!M.length || M[0]!.length < M.length || M.some((c) => c.length !== M[0]!.length))
    throw new RangeError('LLL Gram-Schmidt requires nonempty columns with rows >= columns');
}
/** Native total descending diagonal-exponent drops; diagonal types must agree. */
export function drop(R: GsoMatrix): bigint {
  shape(R);
  let s = 0n,
    m = expo(R[0]![0]!);
  for (let i = 1; i < R.length; i++) {
    const a = R[i]![i]!,
      b = R[i - 1]![i - 1]!;
    let cmp: number;
    if (typeof a === 'bigint' && typeof b === 'bigint')
      cmp = abs(a) < abs(b) ? -1 : abs(a) > abs(b) ? 1 : 0;
    else if (typeof a !== 'bigint' && typeof b !== 'bigint') cmp = abscmprr(a, b);
    else throw new RangeError('drop requires matching integer or real diagonal types');
    if (cmp >= 0) {
      s += m - expo(b);
      m = expo(a);
    }
  }
  return s + m - expo(R[R.length - 1]![R.length - 1]!);
}
/** Native weighted diagonal-exponent potential. */
export function potential(R: GsoMatrix): bigint {
  shape(R);
  let s = 0n,
    mul = BigInt(R.length - 1);
  for (let i = 0; i < R.length; i++, mul -= 2n) s += mul * expo(R[i]![i]!);
  return s;
}
/** Native difference between largest and smallest diagonal exponents. */
export function spread(R: GsoMatrix): bigint {
  shape(R);
  const values = R.map((c, i) => expo(c[i]!));
  return values.reduce(max) - values.reduce((x, y) => (x < y ? x : y));
}
/** Native exponent bound for an invertible triangular matrix's condition number. */
export function condition_bound(U: GsoMatrix, lower = false): bigint {
  shape(U);
  const n = U.length,
    y = Array<bigint>(n);
  let e = (y[n - 1] = -expo(U[n - 1]![n - 1]!));
  for (let i = n - 2; i >= 0; i--) {
    let s = 0n;
    for (let j = i + 1; j < n; j++) s = max(s, expo(lower ? U[i]![j]! : U[j]![i]!) + y[j]!);
    y[i] = s - expo(U[i]![i]!);
    e = max(e, y[i]!);
  }
  return U.flat().reduce<bigint>((s, v) => max(s, expo(v)), -highExpo) + e;
}
/** Native extra bits selected from spread and condition bound. */
export function GS_extraprec(L: GsoMatrix, lower = false): bigint {
  const C = condition_bound(L, lower),
    S = spread(L),
    n = BigInt(L.length);
  return max(2n * S + 2n * n, C - S - 2n * n);
}
const precision = (bits: bigint): number => {
  const p = ((bits + 63n) >> 6n) << 6n;
  if (p <= 0n || p > BigInt(Number.MAX_SAFE_INTEGER))
    throw new RangeError('LLL precision exceeds exact JavaScript integer range');
  return Number(p);
};
/** Native precision-selected conversion of an invertible upper integer matrix. */
export function gramschmidt_upper(M: bigint[][]): MpReal[][] {
  const p = precision(BigInt(M.length + 31) + GS_extraprec(M));
  return M.map((c) => c.map((x) => itor(x, p)));
}
/** Native adaptive QR for an integer basis with full column rank. */
export function gramschmidt_dynprec(M: bigint[][]): QrScalar[][] {
  shape(M);
  if (M[0]!.length === M.length && M.every((c, j) => c.every((x, i) => i <= j || x === 0n)))
    return gramschmidt_upper(M);
  const minprec = BigInt(M.length + 31);
  let bits = minprec;
  for (;;) {
    const p = precision(bits),
      [, , , L] = QR_init(
        M.map((c) => c.map((x) => itor(x, p))),
        p
      );
    if (
      L === null ||
      L.some((c, i) => c[i] === 0n || (typeof c[i] !== 'bigint' && !(c[i] as MpReal).s))
    ) {
      bits *= 2n;
      continue;
    }
    const required = minprec + GS_extraprec(L, true);
    if (bits >= required) return L.map((_, j) => L.map((c) => c[j]!));
    bits = max((4n * bits) / 3n, required);
  }
}
/** Native adaptive Cholesky for a positive-definite integer Gram matrix. */
export function RgM_Cholesky_dynprec(M: bigint[][]): MatrixReal[][] {
  shape(M);
  if (M.some((c) => c.length !== M.length))
    throw new RangeError('adaptive Cholesky requires square columns');
  const minprec = BigInt(M.length + 31);
  let bits = minprec;
  for (;;) {
    const p = precision(bits),
      L = RgM_Cholesky(
        M.map((c) => c.map((x) => itor(x, p))),
        p
      );
    if (L === null) {
      bits *= 2n;
      continue;
    }
    const required = minprec + GS_extraprec(L);
    if (bits >= required) return L;
    bits = max((4n * bits) / 3n, required);
  }
}
