/** nmod_poly_factor/factor_cantor_zassenhaus.c. Input has positive degree. */
import { _nmod_poly_make_monic } from '../nmod_poly/make_monic.js';
import { _nmod_poly_gcd, _nmod_poly_xgcd_kernels as k } from '../nmod_poly/gcd.js';
import { _nmod_poly_powmod_ui_binexp } from '../nmod_poly/powmod_binexp.js';
import { nmod_poly_remove } from '../nmod_poly/remove.js';
import { nmod_poly_factor_equal_deg } from './factor_equal_deg.js';
export function nmod_poly_factor_cantor_zassenhaus(
  f: readonly bigint[],
  p: bigint
): Array<[bigint[], number]> {
  let v = _nmod_poly_make_monic(f, p),
    h = [0n, 1n],
    i = 0;
  const out: Array<[bigint[], number]> = [];
  do {
    i++;
    h = _nmod_poly_powmod_ui_binexp(h, p, v, p);
    const difference = h.length ? h.slice() : [0n];
    difference[1] = (difference[1] ?? 0n) - 1n;
    const g = _nmod_poly_make_monic(_nmod_poly_gcd(k.normalized(difference, p), v, p), p);
    if (g.length !== 1) {
      for (const [factor] of nmod_poly_factor_equal_deg(g, i, p)) {
        let exponent: number;
        [v, exponent] = nmod_poly_remove(v, factor, p);
        out.push([factor, exponent]);
      }
    }
  } while (v.length >= 2 * i + 3);
  if (v.length > 1) out.push([v, 1]);
  return out;
}
