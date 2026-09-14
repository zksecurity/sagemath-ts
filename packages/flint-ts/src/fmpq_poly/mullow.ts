/** @see Deviation: Polynomial Roots and Truncated Series */
/** fmpq_poly/mullow.c: integer numerator product followed by canonicalization. */
import { _fmpz_poly_mullow } from '../fmpz_poly/mullow.js';
import { _fmpz_poly_resultant_kernels as k } from '../fmpz_poly/gcd.js';
export function _fmpq_poly_mullow(
  a: readonly bigint[],
  denA: bigint,
  b: readonly bigint[],
  denB: bigint,
  n: number
): [bigint[], bigint] {
  if (denA <= 0n || denB <= 0n) throw new RangeError('input denominators must be positive');
  const out = _fmpz_poly_mullow(a, b, n),
    den = denA * denB,
    g = k.content([...out, den]);
  return [out.map((c) => c / g), den / g];
}
