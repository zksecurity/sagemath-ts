/** fmpq_poly/mul.c: cross-content cancellation, with the native alias square.
 * @see Deviation: Polynomial Integer Powers and Portable Native Products
 */
import { _fmpz_poly_mul } from '../fmpz_poly/mul.js';
import { _fmpz_poly_resultant_kernels as k } from '../fmpz_poly/gcd.js';
export function _fmpq_poly_mul(
  a: readonly bigint[],
  denA: bigint,
  b: readonly bigint[],
  denB: bigint
): [bigint[], bigint] {
  if (denA <= 0n || denB <= 0n) throw new RangeError('input denominators must be positive');
  const A = k.normalized(a),
    B = a === b ? A : k.normalized(b);
  if (!A.length || !B.length) return [[], 1n];
  const out = _fmpz_poly_mul(A, B),
    den = denA * denB;
  if (A === B) return [out, den];
  const g =
    (denB === 1n ? 1n : k.content([...A, denB])) * (denA === 1n ? 1n : k.content([...B, denA]));
  return [out.map((c) => c / g), den / g];
}
