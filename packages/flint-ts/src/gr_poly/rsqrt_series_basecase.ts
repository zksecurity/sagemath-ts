/** Modular instantiation of gr_poly/rsqrt_series_basecase.c.
 * @see Deviation: FLINT polynomial square-root kernels
 */
import { _gr_poly_sqrt_series_basecase } from './sqrt_series_basecase.js';
import { _gr_poly_inv_series_basecase } from './inv_series_basecase.js';
import { _nmod_poly_xgcd_kernels as k, _nmod_poly_resultant_kernels as r } from '../nmod_poly/gcd.js';
export function _gr_poly_rsqrt_series_basecase(f: readonly bigint[], n: number, p: bigint): bigint[] {
  if (!Number.isSafeInteger(n) || n < 0) throw new RangeError('invalid truncation length');
  if (n === 0) return [];
  const a = k.normalized(f.slice(0, n), p);
  if (a.length <= 1 || n === 2) {
    const root = _gr_poly_sqrt_series_basecase(a, 1, p)[0] ?? 0n;
    const inv = k.inverse(root, p);
    if (a.length <= 1) return [inv];
    return k.normalized([inv, r.mod(-inv * a[1]! * k.inverse(a[0]!, p) * k.inverse(2n, p), p)], p);
  }
  return _gr_poly_inv_series_basecase(_gr_poly_sqrt_series_basecase(a, n, p), n, p);
}
