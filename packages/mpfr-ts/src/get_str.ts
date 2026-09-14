import type { mpfr_t, mpfr_rnd_t } from './types.js';
import { unsupported, rounding, nearest } from './round_raw.js';

/** get_str.c decimal scaled-quotient rounding; returns native digits/exponent.
 * Allocated C buffers are represented by a string and the exponent output.
 */
export function mpfr_get_str(
  base: number,
  digits: number,
  x: mpfr_t,
  mode: mpfr_rnd_t = 'RNDN'
): [string, number] {
  rounding(mode);
  if (base !== 10) unsupported('non-decimal string formatting');
  if (!Number.isSafeInteger(digits) || digits < 0) throw new RangeError('invalid MPFR digit count');
  if (digits === 0) digits = 1 + Math.ceil(x.precision * Math.LOG10E * Math.LN2);
  if (digits > 4096) unsupported('more than 4096 output digits');
  if (x.kind === 'nan') return ['@NaN@', 0];
  const sign = x.sign < 0 ? '-' : '';
  if (x.kind === 'inf') return [sign + '@Inf@', 0];
  if (x.kind === 'zero') return [sign + '0'.repeat(digits), 0];
  const shift = x.exponent - x.precision;
  let n = shift >= 0 ? x.mantissa << BigInt(shift) : x.mantissa,
    d = shift < 0 ? 1n << BigInt(-shift) : 1n;
  let e = n.toString().length - d.toString().length;
  if (e >= 0 ? n < d * 10n ** BigInt(e) : n * 10n ** BigInt(-e) < d) e--;
  e++;
  const scale = digits - e;
  if (scale >= 0) n *= 10n ** BigInt(scale);
  else d *= 10n ** BigInt(-scale);
  let q = nearest(n, d);
  if (q >= 10n ** BigInt(digits)) {
    q /= 10n;
    e++;
  }
  return [sign + q.toString().padStart(digits, '0'), e];
}
