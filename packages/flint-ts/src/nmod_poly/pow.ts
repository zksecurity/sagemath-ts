/** Dense-array port of nmod_poly/pow.c and pow_binexp.c.
 * @see Deviation: Native Modular Polynomial Products and Fraction Fields
 */
import { _nmod_poly_mul } from './mul.js';
export function _nmod_poly_pow(a: readonly bigint[], e: bigint, p: bigint): bigint[] {
  if (e < 0n || e >= 1n << 64n) throw new RangeError('exponent must fit an unsigned word');
  const A = _nmod_poly_mul(a, [1n], p);
  if (!A.length) return e === 0n && p !== 1n ? [1n] : [];
  if (A.length === 1) {
    let n = e,
      base = A[0]!,
      result = 1n;
    while (n) {
      if (n & 1n) result = (result * base) % p;
      n >>= 1n;
      if (n) base = (base * base) % p;
    }
    return result === 0n ? [] : [result];
  }
  if (e === 0n) return [1n];
  if (e === 1n) return A;
  // The native buffer-swap parity only changes storage; retain the exact
  // left-to-right squaring/multiplication sequence with immutable arrays.
  const bits = e.toString(2);
  let result = A;
  for (let i = 1; i < bits.length; i++) {
    result = _nmod_poly_mul(result, result, p);
    if (bits[i] === '1') result = _nmod_poly_mul(result, A, p);
  }
  return result;
}
