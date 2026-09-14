/** fmpz_poly/pseudo_rem.c delegates to divide-and-conquer pseudo division. * @see Deviation: Number-field reduction and native array adapters
 */
import { _fmpz_poly_pseudo_divrem_divconquer } from './pseudo_divrem_divconquer.js';
export function _fmpz_poly_pseudo_rem(
  a: readonly bigint[],
  b: readonly bigint[]
): [bigint[], bigint] {
  const [, r, d] = _fmpz_poly_pseudo_divrem_divconquer(a, b);
  return [r, d];
}
