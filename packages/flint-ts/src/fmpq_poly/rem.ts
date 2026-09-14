/** fmpq_poly/rem.c: integer pseudo remainder followed by denominator clearing. * @see Deviation: Number-field reduction and native array adapters
 */
import { _fmpz_poly_pseudo_rem } from '../fmpz_poly/pseudo_rem.js';
import { _fmpz_poly_resultant_kernels as k } from '../fmpz_poly/gcd.js';

export function _fmpq_poly_rem(
  a: readonly bigint[],
  denA: bigint,
  b: readonly bigint[],
  denB: bigint
): [bigint[], bigint] {
  if (denA <= 0n || denB <= 0n) throw new RangeError('input denominators must be positive');
  const [A, da] = canonical(a, denA),
    [B] = canonical(b, denB);
  if (!B.length) throw new RangeError('division by zero polynomial');
  if (A.length < B.length) return [A, da];
  if (B.length === 1) return [[], 1n];
  const [r, d] = _fmpz_poly_pseudo_rem(A, B);
  return canonical(r, da * B[B.length - 1]! ** d);
}
function canonical(a: readonly bigint[], den: bigint): [bigint[], bigint] {
  const A = k.normalized(a);
  if (!A.length) return [[], 1n];
  let common = k.content([...A, den]);
  if (den < 0n) common = -common;
  return [A.map((c) => c / common), den / common];
}
