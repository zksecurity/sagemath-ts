/** Modular instantiation of gr_poly/inv_series_basecase.c.
 * @see Deviation: FLINT polynomial square-root kernels
 */
import { _nmod_poly_xgcd_kernels as k, _nmod_poly_resultant_kernels as r } from '../nmod_poly/gcd.js';
export function _gr_poly_inv_series_basecase(f: readonly bigint[], n: number, p: bigint): bigint[] {
  if (!Number.isSafeInteger(n) || n < 0) throw new RangeError('invalid truncation length');
  if (n === 0) return [];
  const a = k.normalized(f.slice(0, n), p);
  const inv = k.inverse(a[0] ?? 0n, p), out = Array<bigint>(n).fill(0n);
  out[0] = inv;
  if (a.length === 1) return [inv];
  if (a.length === 2 || a.slice(1, -1).every(c => c === 0n)) {
    const step = a.length - 1, ratio = r.mod(-inv * a[step]!, p);
    for (let i = step; i < n; i += step) out[i] = out[i - step]! * ratio % p;
  } else {
    for (let i = 1; i < n; i++) {
      let sum = 0n;
      for (let j = 1; j <= Math.min(i, a.length - 1); j++) sum += a[j]! * out[i - j]!;
      out[i] = r.mod(-sum * inv, p);
    }
  }
  return k.normalized(out, p);
}
