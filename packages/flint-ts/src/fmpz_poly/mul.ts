/** fmpz_poly/mul.c: native tiny products and bounded portable KS/Karatsuba.
 * @see Deviation: Polynomial Integer Powers and Portable Native Products
 */
import { _fmpz_poly_resultant_kernels as k } from './gcd.js';
export function _fmpz_poly_mul(a: readonly bigint[], b: readonly bigint[]): bigint[] {
  let A = k.normalized(a),
    B = a === b ? A : k.normalized(b);
  if (A.length < B.length) [A, B] = [B, A];
  if (!B.length) return [];
  const bitsA = k.maxBits(A),
    bitsB = k.maxBits(B);
  const tiny =
    A === B
      ? bitsA <= 62 && A.length < 50 + 3 * bitsA && 2 * bitsA + k.bits(BigInt(A.length)) <= 127
      : bitsA <= 62 &&
        bitsB <= 62 &&
        (B.length < 40 + (bitsA + bitsB) / 2 || A.length < 70 + (bitsA + bitsB) / 2) &&
        bitsA + bitsB + k.bits(BigInt(B.length)) <= 127;
  if (B.length < 7 || tiny) {
    const out = Array<bigint>(A.length + B.length - 1).fill(0n);
    if (A === B) {
      for (let i = 0; i < A.length; i++) {
        out[2 * i] = out[2 * i]! + A[i]! * A[i]!;
        for (let j = i + 1; j < A.length; j++) out[i + j] = out[i + j]! + 2n * A[i]! * A[j]!;
      }
    } else
      for (let i = 0; i < A.length; i++)
        for (let j = 0; j < B.length; j++) out[i + j] = out[i + j]! + A[i]! * B[j]!;
    return k.normalized(out);
  }
  return k.multiply(A, B);
}
