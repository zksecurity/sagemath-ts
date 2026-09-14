/** FLINT nmod_poly/sqrt.c: valuation reduction, series root and high-product check.
 * @see Deviation: FLINT polynomial square-root kernels
 */
import { n_sqrtmod } from '../ulong_extras/sqrtmod.js';
import { _nmod_poly_sqrt_series } from './sqrt_series.js';
import { _nmod_poly_mulhigh } from './mulhigh.js';
import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
export function _nmod_poly_sqrt(f: readonly bigint[], p: bigint): bigint[] | null {
  if (p < 2n || p >= 1n << 64n) throw new RangeError('modulus must fit an unsigned word');
  const a = k.normalized(f, p);
  if (!(a.length & 1)) return a.length === 0 ? [] : null;
  if (p === 2n) {
    for (let i = 1; i < a.length; i += 2) if (a[i] !== 0n) return null;
    return a.filter((_, i) => !(i & 1));
  }
  let offset = 0;
  while (a[offset] === 0n) {
    if (a[offset + 1] !== 0n) return null;
    offset += 2;
  }
  const b = a.slice(offset), d = b[0]!, c = d === 1n ? 1n : n_sqrtmod(d, p);
  if (c === 0n) return null;
  if (b.length === 1) return Array<bigint>(offset / 2).fill(0n).concat(c);
  const length = Math.floor(b.length / 2) + 1;
  const low = b.slice(0, length), inv = c === 1n ? 1n : k.inverse(d, p);
  const root = _nmod_poly_sqrt_series(low.map(x => x * inv % p), length, p).map(x => x * c % p);
  const square = _nmod_poly_mulhigh(root, root, length, p);
  for (let i = length; i < b.length; i++) if ((square[i] ?? 0n) !== b[i]) return null;
  return k.normalized(Array<bigint>(offset / 2).fill(0n).concat(root), p);
}
