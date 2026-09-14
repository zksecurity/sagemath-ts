/** nmod_poly_factor/factor_distinct_deg.c: baby/giant-step distinct-degree splitting.
 * Input is squarefree of positive degree. Each pair contains a factor and the
 * common degree of its irreducible factors (not its multiplicity).
 * @see Deviation: Native FLINT factor array kernels
 */
import { _nmod_poly_make_monic } from '../nmod_poly/make_monic.js';
import { _nmod_poly_gcd, _nmod_poly_xgcd_kernels as k } from '../nmod_poly/gcd.js';
import { _nmod_poly_divrem } from '../nmod_poly/divrem.js';
import { _nmod_poly_mul } from '../nmod_poly/mul.js';
import { _nmod_poly_inv_series_newton } from '../nmod_poly/inv_series_newton.js';
import { _nmod_poly_powmod_x_fmpz_preinv } from '../nmod_poly/powmod_x_preinv.js';
import {
  _nmod_poly_powmod_fmpz_binexp_preinv,
  _nmod_poly_preinv_remainder,
} from '../nmod_poly/powmod_binexp_preinv.js';
import {
  nmod_poly_precompute_matrix,
  nmod_poly_compose_mod_brent_kung_precomp_preinv,
  _nmod_poly_reduce_matrix_mod_poly,
} from '../nmod_poly/compose_mod_brent_kung_precomp_preinv.js';
import { nmod_poly_compose_mod_brent_kung_vec_preinv } from '../nmod_poly/compose_mod_brent_kung_vec_preinv.js';
import { nmod_poly_remove } from '../nmod_poly/remove.js';
export function nmod_poly_factor_distinct_deg(
  poly: readonly bigint[],
  p: bigint
): Array<[bigint[], number]> {
  let v = _nmod_poly_make_monic(poly, p);
  const n = v.length - 1;
  if (n < 1) throw new RangeError('polynomial must have positive degree');
  if (n === 1) return [[v, 1]];
  const beta = 0.5 * (1 - Math.log(2) / Math.log(n)),
    l = Math.ceil(n ** beta),
    m = Math.ceil((0.5 * n) / l);
  const inverse = (f: bigint[]) => _nmod_poly_inv_series_newton(f.slice().reverse(), f.length, p);
  let vinv = inverse(v);
  const h: bigint[][] = [[0n, 1n], _nmod_poly_powmod_x_fmpz_preinv(p, v, vinv, p)];
  if (p.toString(2).length > Math.floor(((Math.floor(Math.sqrt(n)) + 1) * 3) / 4)) {
    let i = 1;
    for (; i < BigInt(l).toString(2).length; i++) {
      const count = 2 ** (i - 1);
      h.push(
        ...nmod_poly_compose_mod_brent_kung_vec_preinv(
          h.slice(1, 1 + count),
          count,
          h[count]!,
          v,
          vinv,
          p
        )
      );
    }
    const count = 2 ** (i - 1);
    h.push(
      ...nmod_poly_compose_mod_brent_kung_vec_preinv(
        h.slice(1, 1 + count),
        l - count,
        h[count]!,
        v,
        vinv,
        p
      )
    );
  } else
    for (let i = 2; i <= l; i++)
      h.push(_nmod_poly_powmod_fmpz_binexp_preinv(h[i - 1]!, p, v, vinv, p));
  const H: bigint[][] = [h[l]!],
    I: bigint[][] = Array.from({ length: m }, () => []);
  let HH = nmod_poly_precompute_matrix(H[0]!, v, vinv, p),
    d = 1;
  const difference = (a: bigint[], b: bigint[]) =>
    k.normalized(
      Array.from({ length: Math.max(a.length, b.length) }, (_, i) => (a[i] ?? 0n) - (b[i] ?? 0n)),
      p
    );
  const gcd = (a: bigint[], b: bigint[]) => _nmod_poly_make_monic(_nmod_poly_gcd(a, b, p), p);
  for (let j = 0; j < m; j++) {
    if (j > 0) {
      let previous = H[j - 1]!;
      if (I[j - 1]!.length > 1) {
        HH = _nmod_poly_reduce_matrix_mod_poly(HH, v, p);
        previous = _nmod_poly_divrem(previous, v, p)[1];
      }
      H[j] = nmod_poly_compose_mod_brent_kung_precomp_preinv(previous, HH, v, vinv, p);
    }
    let interval = [1n];
    for (let i = l - 1; i >= 0 && 2 * d <= v.length - 1; i--, d++) {
      const tmp = difference(H[j]!, _nmod_poly_divrem(h[i]!, v, p)[1]);
      interval = _nmod_poly_preinv_remainder(_nmod_poly_mul(tmp, interval, p), v, vinv, p);
    }
    I[j] = gcd(v, interval);
    if (I[j]!.length > 1) {
      v = nmod_poly_remove(v, I[j]!, p)[0];
      vinv = inverse(v);
    }
    if (v.length - 1 < 2 * d) break;
  }
  const out: Array<[bigint[], number]> = [];
  if (v.length > 1) out.push([v, v.length - 1]);
  for (let j = 0; j < m; j++) {
    if (I[j]!.length - 1 > (j + 1) * l || j === 0) {
      let g = I[j]!;
      for (let i = l - 1; i >= 0 && g.length > 1; i--) {
        const f = gcd(g, difference(H[j]!, h[i]!));
        if (f.length > 1) {
          out.push([f, l * (j + 1) - i]);
          g = nmod_poly_remove(g, f, p)[0];
        }
      }
    } else if (I[j]!.length > 1) out.push([_nmod_poly_make_monic(I[j]!, p), I[j]!.length - 1]);
  }
  return out;
}
