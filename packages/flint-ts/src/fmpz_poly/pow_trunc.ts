/** @see Deviation: Polynomial Roots and Truncated Series */
/** fmpz_poly/pow_trunc.c: unsigned-word binary powers with truncated products. */
import { _fmpz_poly_mullow } from './mullow.js';
import { _fmpz_poly_resultant_kernels as k } from './gcd.js';
export function _fmpz_poly_pow_trunc(a: readonly bigint[], e: bigint, n: number): bigint[] {
  if (e < 0n || e >= 1n << 64n) throw new RangeError('exponent must fit an unsigned word');
  if (!Number.isSafeInteger(n) || n < 0) throw new RangeError('invalid truncation length');
  if (!n) return [];
  if (!e) return [1n];
  const A = k.normalized(a.slice(0, n));
  if (!A.length) return [];
  if (A.length === 1) return [A[0]! ** e];
  let result = A;
  for (const bit of e.toString(2).slice(1)) {
    result = _fmpz_poly_mullow(result, result, n);
    if (bit === '1') result = _fmpz_poly_mullow(result, A, n);
  }
  return result;
}
