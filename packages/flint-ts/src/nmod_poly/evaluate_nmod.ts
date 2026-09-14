import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
/** FLINT nmod_poly/evaluate_nmod.c, constant and zero-point shortcuts included.
 * @see Deviation: Polynomial Evaluation and Composition
 */
export function _nmod_poly_evaluate_nmod(a: readonly bigint[], x: bigint, p: bigint): bigint {
  if (p < 2n) throw new RangeError('modulus must be at least 2');
  const A = k.normalized(a, p),
    c = ((x % p) + p) % p;
  if (!A.length) return 0n;
  if (A.length === 1 || c === 0n) return A[0]!;
  let value = A[A.length - 1]!;
  for (let i = A.length - 2; i >= 0; i--) value = (value * c + A[i]!) % p;
  return value;
}
