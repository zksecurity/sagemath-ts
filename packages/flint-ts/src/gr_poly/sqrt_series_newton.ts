/** Modular instantiation of gr_poly/sqrt_series_newton.c, including Karp-Markstein.
 * @see Deviation: FLINT polynomial square-root kernels
 */
import { _gr_poly_sqrt_series_basecase } from './sqrt_series_basecase.js';
import { _gr_poly_rsqrt_series_basecase } from './rsqrt_series_basecase.js';
import { _nmod_poly_mullow } from '../nmod_poly/mullow.js';
import { _nmod_poly_xgcd_kernels as k, _nmod_poly_resultant_kernels as r } from '../nmod_poly/gcd.js';
export function _gr_poly_sqrt_series_newton(f: readonly bigint[], n: number, cutoff: number, p: bigint): bigint[] {
  if (!Number.isSafeInteger(n) || n < 0 || !Number.isSafeInteger(cutoff)) throw new RangeError('invalid truncation length or cutoff');
  if (n === 0) return [];
  if (n < cutoff) return _gr_poly_sqrt_series_basecase(f, n, p);
  cutoff = Math.max(cutoff, 2);
  const a = k.normalized(f.slice(0, n), p), sizes = [n], half = k.inverse(2n, p);
  while (sizes[sizes.length - 1]! >= cutoff) sizes.push(Math.ceil(sizes[sizes.length - 1]! / 2));
  let m = sizes[sizes.length - 1]!;
  let g = _gr_poly_rsqrt_series_basecase(a.slice(0, m), m, p);
  for (let i = sizes.length - 2; i >= 1; i--) {
    const length = sizes[i]!;
    const tlen = Math.min(2 * m - 1, length), ulen = Math.min(length, m + tlen - 1);
    const t = _nmod_poly_mullow(g, g, tlen, p);
    const u = _nmod_poly_mullow(g, t, ulen, p);
    const product = _nmod_poly_mullow(u, a.slice(0, length), length, p);
    while (g.length < m) g.push(0n);
    for (let j = m; j < length; j++) g[j] = r.mod(-(product[j] ?? 0n) * half, p);
    m = length;
  }
  // Convert the final inverse-root approximation to a root with one high correction.
  m = Math.ceil(n / 2);
  const v = _nmod_poly_mullow(g, a, m, p);
  const square = _nmod_poly_mullow(v, v, Math.min(2 * m - 1, n), p);
  const high = Array.from({length: n - m}, (_, i) => r.mod((a[m + i] ?? 0n) - (square[m + i] ?? 0n), p));
  const correction = _nmod_poly_mullow(g, high, n - m, p);
  const out = Array.from({length: n}, (_, i) => i < m ? v[i] ?? 0n : (correction[i - m] ?? 0n) * half % p);
  return k.normalized(out, p);
}
