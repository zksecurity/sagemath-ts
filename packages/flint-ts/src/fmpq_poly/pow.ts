/** fmpq_poly/pow.c: integer numerator power and common-denominator power. */
import { _fmpz_poly_resultant_kernels as k } from '../fmpz_poly/gcd.js';
import { _fmpz_poly_pow } from '../fmpz_poly/pow.js';
export function _fmpq_poly_pow(a: readonly bigint[], den: bigint, e: bigint): [bigint[], bigint] {
  if (den <= 0n) throw new RangeError('input denominator must be positive');
  const A = k.normalized([...a]);
  const numerator = _fmpz_poly_pow(A, e);
  return [numerator, numerator.length ? den ** e : 1n];
}
