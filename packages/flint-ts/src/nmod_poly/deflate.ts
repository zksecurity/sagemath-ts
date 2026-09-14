/** nmod_poly/deflate.c. The native result length is retained, even for a non-dividing stride.
 * @see Deviation: Native FLINT factor array kernels
 */
import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
export function nmod_poly_deflate(f: readonly bigint[], deflation: number, p: bigint): bigint[] {
  if (!Number.isSafeInteger(deflation) || deflation <= 0)
    throw new RangeError('deflation must be positive');
  const a = k.normalized(f, p);
  if (a.length <= 1 || deflation === 1) return a;
  return Array.from(
    { length: Math.floor((a.length - 1) / deflation) + 1 },
    (_, i) => a[i * deflation]!
  );
}
