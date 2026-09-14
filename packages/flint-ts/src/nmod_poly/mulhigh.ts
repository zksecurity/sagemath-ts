/** nmod_poly/mulhigh.c classical/KS dispatch.
 * The source leaves the low part zero for classical products and computes it
 * for KS; only coefficients at or above start have a mathematical guarantee.
 * @see Deviation: FLINT polynomial square-root kernels
 */
import { _nmod_poly_mulhigh_classical } from './mulhigh_classical.js';
import { _nmod_poly_mul_KS } from './mul.js';
import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
export function _nmod_poly_mulhigh(a: readonly bigint[], b: readonly bigint[], start: number, p: bigint): bigint[] {
  if (!Number.isSafeInteger(start) || start < 0) throw new RangeError('invalid start index');
  if (p < 2n || p >= 1n << 64n) throw new RangeError('modulus must fit an unsigned word');
  let A = k.normalized(a, p), B = k.normalized(b, p);
  if (A.length < B.length) [A, B] = [B, A];
  if (!B.length || start >= A.length + B.length - 1) return [];
  if (A.length + B.length <= 6 ||
      (2 * p.toString(2).length + A.length.toString(2).length <= 64 && A.length + B.length < 16))
    return _nmod_poly_mulhigh_classical(A, B, start, p);
  return _nmod_poly_mul_KS(A, B, p);
}
