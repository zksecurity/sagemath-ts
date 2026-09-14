/** nmod_poly/inflate.c: zero inflation evaluates at one.
 * @see Deviation: Native FLINT factor array kernels
 */
import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
import { _nmod_poly_evaluate_nmod } from './evaluate_nmod.js';
export function nmod_poly_inflate(f: readonly bigint[], inflation: number, p: bigint): bigint[] {
  if (!Number.isSafeInteger(inflation) || inflation < 0)
    throw new RangeError('inflation must be nonnegative');
  const a = k.normalized(f, p);
  if (a.length <= 1 || inflation === 1) return a;
  if (inflation === 0) return k.normalized([_nmod_poly_evaluate_nmod(a, 1n, p)], p);
  const result = Array<bigint>((a.length - 1) * inflation + 1).fill(0n);
  for (let i = 0; i < a.length; i++) result[i * inflation] = a[i]!;
  return result;
}
