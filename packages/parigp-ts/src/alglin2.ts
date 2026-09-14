import { type MatrixReal, matrixRealInv, matrixRealMul, matrixRealSub } from './_real_matrix.js';
import { sqrtr_abs } from './qfb.js';
import { PariError } from './errors.js';

/** Native positive Gaussian reduction for real matrices, using the upper triangle.
 * Columns contain MpReal or exact integer zero cells.
 * @see Deviation: PARI real Cholesky and triangular inverse adapters
 */
export function qfgaussred_positive(a: MatrixReal[][]): MatrixReal[][] | null {
  const n = a.length;
  if (a.some((c) => c.length !== n))
    throw new PariError('inconsistent dimensions in qfgaussred_positive');
  const b = a.map((c, j) => c.map((x, i) => (i <= j ? x : 0n)));
  for (let k = 0; k < n; k++) {
    const p = b[k]![k]!;
    if (p === 0n || p.s <= 0) return null;
    const invp = matrixRealInv(p),
      bk = b.map((c) => c[k]!);
    for (let i = k + 1; i < n; i++) b[i]![k] = matrixRealMul(bk[i]!, invp);
    for (let i = k + 1; i < n; i++) {
      const c = bk[i]!;
      for (let j = i; j < n; j++) b[j]![i] = matrixRealSub(b[j]![i]!, matrixRealMul(c, b[j]![k]!));
    }
  }
  return b;
}

/** Native Cholesky on real/zero columns; null for a nonpositive pivot.
 * Existing real cells determine square-root precision, as in native gsqrt.
 * @see Deviation: PARI real Cholesky and triangular inverse adapters
 */
export function RgM_Cholesky(M: MatrixReal[][], precision = 64): MatrixReal[][] | null {
  if (!Number.isSafeInteger(precision) || precision <= 0 || precision % 64)
    throw new RangeError('RgM_Cholesky precision must be a positive multiple of 64 bits');
  const L = qfgaussred_positive(M);
  if (L === null) return null;
  const R: MatrixReal[][] = L.map((c) => Array<MatrixReal>(c.length).fill(0n));
  for (let i = 0; i < L.length; i++) {
    const r = sqrtr_abs(L[i]![i] as Exclude<MatrixReal, 0n>);
    for (let j = 0; j < L.length; j++) R[j]![i] = i === j ? r : matrixRealMul(r, L[j]![i]!);
  }
  return R;
}

import { RgXQ_norm, RgXQ_trace } from './RgX.js';
import {
  type RationalPair,
  type RationalPolynomialData,
  normalizedPolynomial,
  rationalPair,
} from './_rational_polynomial.js';
/** Reduced QQ polynomial modulo a nonconstant rational polynomial. */
export interface RationalPolMod {
  value: RationalPolynomialData;
  modulus: RationalPolynomialData;
}
/** Native gnorm t_POLMOD branch, including scalar powering without a matrix.
 * @see Deviation: PARI rational trace and norm adapters
 */
export function gnorm(x: RationalPolMod): RationalPair {
  const a = normalizedPolynomial(x.value),
    T = normalizedPolynomial(x.modulus);
  if (a[0].length <= 1) {
    const n = BigInt(T[0].length - 1);
    return rationalPair((a[0][0] ?? 0n) ** n, a[1] ** n);
  }
  return RgXQ_norm(a, T);
}
/** Native gtrace t_POLMOD branch, including scalar degree multiplication.
 * @see Deviation: PARI rational trace and norm adapters
 */
export function gtrace(x: RationalPolMod): RationalPair {
  const a = normalizedPolynomial(x.value),
    T = normalizedPolynomial(x.modulus);
  if (a[0].length <= 1) return rationalPair((a[0][0] ?? 0n) * BigInt(T[0].length - 1), a[1]);
  return RgXQ_trace(a, T);
}
