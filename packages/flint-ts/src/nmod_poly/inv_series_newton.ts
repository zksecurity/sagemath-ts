/** nmod_poly/gr_poly inverse-series Newton lifting with bounded native products.
 * @see Deviation: Polynomial Modular Powers
 */
import { _nmod_poly_mullow } from './mullow.js';
import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
export function _nmod_poly_inv_series_newton(a: readonly bigint[], n: number, p: bigint): bigint[] {
  if (p < 2n) throw new RangeError('modulus must be at least 2');
  if (!Number.isSafeInteger(n) || n <= 0) throw new RangeError('invalid truncation length');
  const A = k.normalized(a.slice(0, n), p),
    inv = k.inverse(A[0] ?? 0n, p);
  if (A.length === 1) return [inv];
  let g = [inv];
  if (A.length <= 10) {
    for (let i = 1; i < n; i++) {
      let sum = 0n;
      for (let j = 1; j <= Math.min(i, A.length - 1); j++) sum += A[j]! * g[i - j]!;
      g.push((((-sum * inv) % p) + p) % p);
    }
    return k.normalized(g, p);
  }
  const sizes = [n];
  while (sizes[sizes.length - 1]! > 1) sizes.push(Math.ceil(sizes[sizes.length - 1]! / 2));
  for (let i = sizes.length - 1; i > 0; i--) {
    const m = sizes[i]!,
      length = sizes[i - 1]!,
      high = _nmod_poly_mullow(A.slice(0, length), g, length, p).slice(m);
    const correction = _nmod_poly_mullow(g, high, length - m, p);
    while (g.length < m) g.push(0n);
    for (const c of correction) g.push(c === 0n ? 0n : p - c);
  }
  return k.normalized(g, p);
}
