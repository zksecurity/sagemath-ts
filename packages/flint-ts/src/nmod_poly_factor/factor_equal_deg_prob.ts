/** nmod_poly_factor/factor_equal_deg_prob.c. A failed split returns null.
 * @see Deviation: Native FLINT factor array kernels
 */
import type { flint_rand_t } from '../flint.js';
import { nmod_poly_randtest } from '../nmod_poly/randtest.js';
import { _nmod_poly_gcd, _nmod_poly_xgcd_kernels as k } from '../nmod_poly/gcd.js';
import { _nmod_poly_make_monic } from '../nmod_poly/make_monic.js';
import { _nmod_poly_divrem } from '../nmod_poly/divrem.js';
import { _nmod_poly_inv_series_newton } from '../nmod_poly/inv_series_newton.js';
import { _nmod_poly_powmod_fmpz_binexp_preinv as pow } from '../nmod_poly/powmod_binexp_preinv.js';
export function nmod_poly_factor_equal_deg_prob(
  state: flint_rand_t,
  pol: readonly bigint[],
  d: number,
  p: bigint
): bigint[] | null {
  const f = k.normalized(pol, p);
  if (f.length < 3 || !Number.isSafeInteger(d) || d <= 0)
    throw new RangeError('equal-degree splitting requires degree at least two and positive d');
  let a: bigint[];
  do {
    a = nmod_poly_randtest(state, f.length - 1, p);
  } while (a.length <= 1);
  const gcd = (a: bigint[], b: bigint[]) => _nmod_poly_make_monic(_nmod_poly_gcd(a, b, p), p);
  const candidate = gcd(a, f);
  if (candidate.length !== 1) return candidate;
  const inv = _nmod_poly_inv_series_newton(f.slice().reverse(), f.length, p);
  let b: bigint[];
  if (p > 2n) b = pow(a, (p ** BigInt(d) - 1n) / 2n, f, inv, p);
  else {
    b = _nmod_poly_divrem(a, f, p)[1];
    let c = b;
    for (let i = 1; i < d; i++) {
      c = pow(c, 2n, f, inv, p);
      b = k.normalized(
        Array.from({ length: Math.max(b.length, c.length) }, (_, j) => (b[j] ?? 0n) + (c[j] ?? 0n)),
        p
      );
    }
    b = _nmod_poly_divrem(b, f, p)[1];
  }
  b[0] = ((b[0] ?? 0n) - 1n + p) % p;
  const result = gcd(k.normalized(b, p), f);
  return result.length <= 1 || result.length === f.length ? null : result;
}
