/** Single-thread nmod_poly_factor/factor_kaltofen_shoup.c.
 * @see Deviation: Native FLINT factor array kernels
 */
import { _nmod_poly_xgcd_kernels as k } from '../nmod_poly/gcd.js';
import { _nmod_poly_make_monic } from '../nmod_poly/make_monic.js';
import { nmod_poly_remove } from '../nmod_poly/remove.js';
import { nmod_poly_factor_squarefree } from './factor_squarefree.js';
import { nmod_poly_factor_distinct_deg } from './factor_distinct_deg.js';
import { nmod_poly_factor_equal_deg } from './factor_equal_deg.js';
export function nmod_poly_factor_kaltofen_shoup(
  poly: readonly bigint[],
  p: bigint
): Array<[bigint[], number]> {
  const input = k.normalized(poly, p);
  if (!input.length) throw new RangeError('polynomial division by zero');
  if (input.length === 1) return [];
  let v = _nmod_poly_make_monic(input, p);
  if (v.length === 2) return [[v, 1]];
  const out: Array<[bigint[], number]> = [];
  for (const [squarefree] of nmod_poly_factor_squarefree(v, p))
    for (const [component, degree] of nmod_poly_factor_distinct_deg(squarefree, p))
      for (const [factor] of nmod_poly_factor_equal_deg(component, degree, p)) {
        let exponent: number;
        [v, exponent] = nmod_poly_remove(v, factor, p);
        out.push([factor, exponent]);
      }
  return out;
}
