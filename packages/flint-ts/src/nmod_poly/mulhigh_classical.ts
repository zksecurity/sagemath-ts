/** nmod_poly/mulhigh_classical.c with exact BigInt accumulation.
 * @see Deviation: FLINT polynomial square-root kernels
 */
import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
export function _nmod_poly_mulhigh_classical(a: readonly bigint[], b: readonly bigint[], start: number, p: bigint): bigint[] {
  if (!Number.isSafeInteger(start) || start < 0) throw new RangeError('invalid start index');
  if (p < 2n || p >= 1n << 64n) throw new RangeError('modulus must fit an unsigned word');
  const A = k.normalized(a, p), B = k.normalized(b, p), length = A.length + B.length - 1;
  if (!A.length || !B.length || start >= length) return [];
  const out = Array<bigint>(length).fill(0n);
  for (let i = 0; i < A.length; i++)
    for (let j = Math.max(0, start - i); j < B.length; j++) out[i + j] = out[i + j]! + A[i]! * B[j]!;
  return k.normalized(out, p);
}
