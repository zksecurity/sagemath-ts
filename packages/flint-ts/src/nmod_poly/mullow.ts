/** @see Deviation: Polynomial Roots and Truncated Series */
/** nmod_poly/mullow.c portable dispatch (without fft_small).
 * @see Deviation: Polynomial Integer Powers and Portable Native Products
 */
import { _nmod_poly_mul, _nmod_poly_mul_KS } from './mul.js';
import { _fmpz_poly_resultant_kernels as k } from '../fmpz_poly/gcd.js';
export function _nmod_poly_mullow(
  a: readonly bigint[],
  b: readonly bigint[],
  n: number,
  p: bigint
): bigint[] {
  if (!Number.isSafeInteger(n) || n < 0) throw new RangeError('invalid truncation length');
  if (p < 1n || p >= 1n << 64n) throw new RangeError('modulus must fit a positive unsigned word');
  const normalized = (x: readonly bigint[]) => k.normalized(x.map((c) => ((c % p) + p) % p));
  let A = normalized(a.slice(0, n)),
    B = a === b ? A : normalized(b.slice(0, n));
  if (A.length < B.length) [A, B] = [B, A];
  if (!B.length) return [];
  n = Math.min(n, A.length + B.length - 1);
  const bits = p.toString(2).length;
  if (B.length <= 5 || (n !== A.length + B.length - 1 && n < 10 + Math.floor((bits * bits) / 10))) {
    const out = Array<bigint>(n).fill(0n);
    for (let i = 0; i < A.length; i++)
      for (let j = 0; j < Math.min(B.length, n - i); j++) out[i + j] = out[i + j]! + A[i]! * B[j]!;
    return normalized(out);
  }
  if ((A.length + B.length - 1) * (2 * bits + k.bits(BigInt(B.length))) > 1 << 19)
    return normalized(k.multiply(A, B).slice(0, n));
  return normalized(
    (n === A.length + B.length - 1 ? _nmod_poly_mul(A, B, p) : _nmod_poly_mul_KS(A, B, p)).slice(
      0,
      n
    )
  );
}
