/** @see Deviation: Polynomial Roots and Truncated Series */
/** fmpz_poly/mullow.c: truncated classical/packed multiplication.
 * @see Deviation: Polynomial Integer Powers and Portable Native Products
 */
import { _fmpz_poly_resultant_kernels as k } from './gcd.js';
export function _fmpz_poly_mullow(a: readonly bigint[], b: readonly bigint[], n: number): bigint[] {
  if (!Number.isSafeInteger(n) || n < 0) throw new RangeError('invalid truncation length');
  let A = k.normalized(a.slice(0, n)),
    B = a === b ? A : k.normalized(b.slice(0, n));
  if (A.length < B.length) [A, B] = [B, A];
  if (!B.length) return [];
  n = Math.min(n, A.length + B.length - 1);
  const bits1 = k.maxBits(A),
    bits2 = k.maxBits(B);
  const tiny =
    bits1 <= 62 &&
    bits2 <= 62 &&
    (B.length < 50 || (4 * B.length >= 3 * n && n < 150 + bits1 + bits2)) &&
    bits1 + bits2 + k.bits(BigInt(B.length)) <= 127;
  if (B.length < 7 || tiny) {
    const out = Array<bigint>(n).fill(0n);
    for (let i = 0; i < A.length; i++)
      for (let j = 0; j < Math.min(B.length, n - i); j++) out[i + j] = out[i + j]! + A[i]! * B[j]!;
    return k.normalized(out);
  }
  // The portable packed/Karatsuba boundary bounds both input degrees first.
  return k.normalized(k.multiply(A, B).slice(0, n));
}
