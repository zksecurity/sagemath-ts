/** Modular instantiation of gr_poly/sqrt_series_basecase.c.
 * @see Deviation: FLINT polynomial square-root kernels
 */
import { n_sqrtmod } from '../ulong_extras/sqrtmod.js';
import { _nmod_poly_xgcd_kernels as k, _nmod_poly_resultant_kernels as r } from '../nmod_poly/gcd.js';
export function _gr_poly_sqrt_series_basecase(f: readonly bigint[], n: number, p: bigint): bigint[] {
  if (!Number.isSafeInteger(n) || n < 0) throw new RangeError('invalid truncation length');
  if (n === 0) return [];
  const a = k.normalized(f.slice(0, n), p), constant = a[0] ?? 0n;
  const root = n_sqrtmod(constant, p);
  if (root === 0n && constant !== 0n) throw new RangeError('constant term is not a square');
  if (a.length <= 1) return root ? [root] : [];
  const half = k.inverse(2n, p), invroot = k.inverse(root, p);
  const out = Array<bigint>(n).fill(0n);
  out[0] = root;
  for (let i = 1; i < n; i++) {
    let value = a[i] ?? 0n;
    if (!(i & 1)) value -= out[i / 2]! ** 2n;
    value = r.mod(value * half, p);
    for (let j = 1; j <= Math.floor((i - 1) / 2); j++) value -= out[j]! * out[i - j]!;
    out[i] = r.mod(value * invroot, p);
  }
  return k.normalized(out, p);
}
