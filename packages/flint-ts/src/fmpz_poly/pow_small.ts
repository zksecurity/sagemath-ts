/** fmpz_poly/pow_small.c, using the shared exact FLINT product adapter. */
import { _fmpz_poly_resultant_kernels as k } from './gcd.js';
export function _fmpz_poly_pow_small(a: readonly bigint[], e: bigint): bigint[] {
  if (e < 0n || e > 4n) throw new RangeError('small exponent must be between 0 and 4');
  const A = k.normalized([...a]);
  if (e === 0n) return [1n];
  if (e === 1n) return A;
  const square = k.multiply(A, A);
  return e === 2n ? square : e === 3n ? k.multiply(square, A) : k.multiply(square, square);
}
