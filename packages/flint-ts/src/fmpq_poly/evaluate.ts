import { _fmpz_poly_evaluate_fmpz } from '../fmpz_poly/evaluate_fmpz.js';
import { _fmpz_poly_evaluate_fmpq } from '../fmpz_poly/evaluate_fmpq.js';
/** FLINT's internal rational evaluation at an integer (unreduced output).
 * @see Deviation: Polynomial Evaluation and Composition
 */
export function _fmpq_poly_evaluate_fmpz(
  a: readonly bigint[],
  denominator: bigint,
  x: bigint
): [bigint, bigint] {
  if (denominator <= 0n) throw new RangeError('input denominator must be positive');
  return [_fmpz_poly_evaluate_fmpz(a, x), denominator];
}
/** FLINT's internal rational evaluation at a rational (unreduced output).
 * @see Deviation: Polynomial Evaluation and Composition
 */
export function _fmpq_poly_evaluate_fmpq(
  a: readonly bigint[],
  denominator: bigint,
  numerator: bigint,
  pointDenominator: bigint
): [bigint, bigint] {
  if (denominator <= 0n) throw new RangeError('input denominator must be positive');
  const [num, den] = _fmpz_poly_evaluate_fmpq(a, numerator, pointDenominator);
  return [num, den * denominator];
}
