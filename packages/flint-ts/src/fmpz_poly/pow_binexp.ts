/** fmpz_poly/pow_binexp.c: left-to-right binary powering; native buffer swaps are implicit. */
import { _fmpz_poly_resultant_kernels as k } from './gcd.js';
import { _fmpz_poly_pow_small } from './pow_small.js';
export function _fmpz_poly_pow_binexp(a: readonly bigint[], e: bigint): bigint[] {
  if (e < 0n || e >= 1n << 64n) throw new RangeError('exponent must fit an unsigned word');
  const A = k.normalized([...a]);
  if (e < 3n) return _fmpz_poly_pow_small(A, e);
  if (!A.length) return [];
  if (A.length === 1) return [A[0]! ** e];
  const bits = e.toString(2);
  let result = A;
  for (let i = 1; i < bits.length; i++) {
    result = k.multiply(result, result);
    if (bits[i] === '1') result = k.multiply(result, A);
  }
  return result;
}
