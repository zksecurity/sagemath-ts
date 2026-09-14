/** nmod_poly/remove.c. Return the updated polynomial and factor multiplicity.
 * @see Deviation: Native FLINT factor array kernels
 */
import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
import { _nmod_poly_divrem } from './divrem.js';
export function nmod_poly_remove(
  f: readonly bigint[],
  g: readonly bigint[],
  p: bigint
): [bigint[], number] {
  let a = k.normalized(f, p);
  const b = k.normalized(g, p);
  if (b.length <= 1) throw new RangeError('factor must have positive degree');
  let exponent = 0;
  while (a.length >= b.length) {
    const [q, r] = _nmod_poly_divrem(a, b, p);
    if (r.length) break;
    a = q;
    exponent++;
  }
  return [a, exponent];
}
