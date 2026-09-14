/** Recursive nmod_poly_factor/factor_equal_deg.c, including fresh state per split.
 * Input must be squarefree with irreducible factors all of degree d.
 * @see Deviation: Native FLINT factor array kernels
 */
import { flint_rand_init } from '../flint.js';
import { _nmod_poly_divrem } from '../nmod_poly/divrem.js';
import { _nmod_poly_xgcd_kernels as k } from '../nmod_poly/gcd.js';
import { nmod_poly_factor_equal_deg_prob } from './factor_equal_deg_prob.js';
export function nmod_poly_factor_equal_deg(
  pol: readonly bigint[],
  d: number,
  p: bigint
): Array<[bigint[], number]> {
  const f = k.normalized(pol, p);
  if (!Number.isSafeInteger(d) || d <= 0 || f.length <= 1 || (f.length - 1) % d)
    throw new RangeError('invalid equal-degree factorization input');
  if (f.length === d + 1) return [[f, 1]];
  const state = flint_rand_init();
  let factor: bigint[] | null;
  do {
    factor = nmod_poly_factor_equal_deg_prob(state, f, d, p);
  } while (!factor);
  const quotient = _nmod_poly_divrem(f, factor, p)[0];
  return [
    ...nmod_poly_factor_equal_deg(factor, d, p),
    ...nmod_poly_factor_equal_deg(quotient, d, p),
  ];
}
