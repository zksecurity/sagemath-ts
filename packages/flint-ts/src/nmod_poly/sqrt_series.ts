/** nmod_poly/sqrt_series.c and the portable gr/nmod.c cutoff dispatch.
 * @see Deviation: FLINT polynomial square-root kernels
 */
import { _gr_poly_sqrt_series_basecase } from '../gr_poly/sqrt_series_basecase.js';
import { _gr_poly_sqrt_series_newton } from '../gr_poly/sqrt_series_newton.js';
const CUTOFF = [32767,632,732,928,1443,1731,2364,2490,2893,3173,5316,5412,5727,6123,6613,7290,7572,8023,9114,9105,8656,10645,11290,13223,11507,15489,12328,9338,9517,9795,9596,13162,10168,8013,9949,10654,9932,12222,11287,11623,11971,12577,12207,13886,14160,12200,13207,15943,15320,14290,15933,15463,14281,15457,15302,17929,18106,17058,14844,17740,17916,18640,18093,18638];
export function _nmod_poly_sqrt_series(f: readonly bigint[], n: number, p: bigint): bigint[] {
  if (p < 2n || p >= 1n << 64n) throw new RangeError('modulus must fit an unsigned word');
  if (!Number.isSafeInteger(n) || n < 0) throw new RangeError('invalid truncation length');
  if (!f.length || n === 0) return [];
  const cutoff = CUTOFF[p.toString(2).length - 1]!;
  return Math.min(f.length, n) < cutoff
    ? _gr_poly_sqrt_series_basecase(f, n, p)
    : _gr_poly_sqrt_series_newton(f, n, cutoff, p);
}
