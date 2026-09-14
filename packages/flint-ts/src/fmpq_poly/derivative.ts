import { _fmpz_poly_derivative } from '../fmpz_poly/derivative.js';
import { _fmpz_poly_resultant_kernels as k } from '../fmpz_poly/gcd.js';
/** Native rational derivative, with canonical numerator/denominator output. */
export function _fmpq_poly_derivative(a: readonly bigint[], den: bigint): [bigint[], bigint] {
  if (den <= 0n) throw new RangeError('input denominator must be positive');
  const coeffs = _fmpz_poly_derivative(a);
  const c = k.content([...coeffs, den]);
  return [coeffs.map((x) => x / c), den / c];
}
