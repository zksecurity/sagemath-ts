/** Dense port of FLINT nmod_poly_factor/factor_squarefree.c.
 * @see Deviation: Polynomial squarefree backend dispatch
 */
import { _nmod_poly_xgcd_kernels as k } from '../nmod_poly/gcd.js';
import { _nmod_poly_gcd } from '../nmod_poly/gcd.js';
import { _nmod_poly_make_monic } from '../nmod_poly/make_monic.js';
import { _nmod_poly_divrem } from '../nmod_poly/divrem.js';
export function nmod_poly_factor_squarefree(f: readonly bigint[], p: bigint): Array<[bigint[], number]> {
  const a = k.normalized(f, p), out: Array<[bigint[], number]> = [];
  if (a.length <= 1) return out;
  const monic = (a: bigint[]) => _nmod_poly_make_monic(a, p);
  if (a.length === 2) return [[monic(a), 1]];
  const gcd = (a: bigint[], b: bigint[]) => monic(_nmod_poly_gcd(a, b, p));
  const quotient = (a: bigint[], b: bigint[]) => _nmod_poly_divrem(a, b, p)[0];
  const root = (a: bigint[]) => Array.from({length: Math.floor((a.length - 1) / Number(p)) + 1}, (_, i) => a[i * Number(p)]!);
  const derivative = k.normalized(a.slice(1).map((x, i) => x * BigInt(i + 1)), p);
  if (!derivative.length) return nmod_poly_factor_squarefree(root(a), p).map(([g, e]) => [g, e * Number(p)]);
  let g = gcd(a, derivative), g1 = quotient(a, g), i = 1;
  while (!(g1.length === 1 && g1[0] === 1n)) {
    const h = gcd(g1, g), z = quotient(g1, h);
    if (z.length > 1) out.push([monic(z), i]);
    i++; g1 = h; g = quotient(g, h);
  }
  g = monic(g);
  if (!(g.length === 1 && g[0] === 1n))
    out.push(...nmod_poly_factor_squarefree(root(g), p).map(([f, e]) => [f, e * Number(p)] as [bigint[], number]));
  return out;
}
