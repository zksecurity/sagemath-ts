import type { mpfr_t, mpfr_rnd_t } from './types.js';
import { rounding } from './round_raw.js';
import { roundedDyadic, singular } from './arithmetic.js';

/** MPFR mul.c full mantissa product, normalization and nearest-even rounding.
 * @see Deviation: Native MPFR binary arithmetic
 */
export function mpfr_mul(a: mpfr_t, b: mpfr_t, c: mpfr_t, mode: mpfr_rnd_t = 'RNDN'): number {
  rounding(mode);
  if (b.kind === 'nan' || c.kind === 'nan') return singular(a, 'nan');
  const sign = (b.sign * c.sign) as 1 | -1;
  if (b.kind === 'inf' || c.kind === 'inf') {
    return b.kind === 'zero' || c.kind === 'zero' ? singular(a, 'nan') : singular(a, 'inf', sign);
  }
  if (b.kind === 'zero' || c.kind === 'zero') return singular(a, 'zero', sign);
  return roundedDyadic(
    a,
    b.mantissa * c.mantissa,
    b.exponent + c.exponent - b.precision - c.precision,
    sign
  );
}
