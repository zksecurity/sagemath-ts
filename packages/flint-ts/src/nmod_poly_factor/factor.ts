/** nmod_poly_factor/factor.c: native deflation, squarefree and CZ/KS dispatch.
 * Returns [leading coefficient, factors in native insertion order].
 * @see Deviation: Native FLINT factor array kernels
 */
import { _nmod_poly_xgcd_kernels as k } from '../nmod_poly/gcd.js';
import { _nmod_poly_make_monic } from '../nmod_poly/make_monic.js';
import { nmod_poly_deflation } from '../nmod_poly/deflation.js';
import { nmod_poly_deflate } from '../nmod_poly/deflate.js';
import { nmod_poly_inflate } from '../nmod_poly/inflate.js';
import { nmod_poly_factor_squarefree } from './factor_squarefree.js';
import { nmod_poly_factor_cantor_zassenhaus } from './factor_cantor_zassenhaus.js';
import { nmod_poly_factor_kaltofen_shoup } from './factor_kaltofen_shoup.js';
type Factors = Array<[bigint[], number]>;
type Result = [bigint, Factors];
function factorWithoutDeflation(f: bigint[], p: bigint, ks: boolean): Result {
  if (f.length <= 1) return [f[0] ?? 0n, []];
  const unit = f[f.length - 1]!,
    monic = _nmod_poly_make_monic(f, p);
  if (f.length === 2) return [unit, [[monic, 1]]];
  const out: Factors = [];
  for (const [s, e] of nmod_poly_factor_squarefree(monic, p)) {
    const factors = ks
      ? nmod_poly_factor_kaltofen_shoup(s, p)
      : nmod_poly_factor_cantor_zassenhaus(s, p);
    for (const [g, m] of factors) out.push([g, m * e]);
  }
  return [unit, out];
}
function factorWithDeflation(input: readonly bigint[], p: bigint, ks: boolean): Result {
  const f = k.normalized(input, p);
  if (f.length <= 1) return [f[0] ?? 0n, []];
  const stride = nmod_poly_deflation(f, p);
  if (stride === 1) return factorWithoutDeflation(f, p, ks);
  const [unit, factors] = factorWithoutDeflation(nmod_poly_deflate(f, stride, p), p, ks),
    out: Factors = [];
  for (const [g, e] of factors) {
    const [, inflated] = factorWithoutDeflation(nmod_poly_inflate(g, stride, p), p, ks);
    for (const [h, m] of inflated) out.push([h, m * e]);
  }
  return [unit, out];
}
export function nmod_poly_factor(input: readonly bigint[], p: bigint): Result {
  const f = k.normalized(input, p);
  return factorWithDeflation(f, p, f.length - 1 >= 10 + Math.floor(50 / p.toString(2).length));
}
export function nmod_poly_factor_with_cantor_zassenhaus(
  input: readonly bigint[],
  p: bigint
): Result {
  return factorWithDeflation(input, p, false);
}
export function nmod_poly_factor_with_kaltofen_shoup(input: readonly bigint[], p: bigint): Result {
  return factorWithDeflation(input, p, true);
}

// Preserve the pre-existing native factor-container type alongside the function.
export type nmod_poly_factor = import('../nmod_poly.js').nmod_poly_factor;
