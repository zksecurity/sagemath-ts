import type { mpfr_t, mpfr_rnd_t } from './types.js';
import { rounding } from './round_raw.js';
/** MPFR 4.2.1 set.c/set4: copy the sign, including NaN, and round the significand.
 * @see Deviation: Native Real Fractional Parts and Comparisons
 */
export function mpfr_set(x: mpfr_t, y: mpfr_t, mode: mpfr_rnd_t = 'RNDN'): number {
  rounding(mode);
  const source = { ...y };
  x.sign = source.sign;
  x.kind = source.kind;
  x.exponent = source.exponent;
  if (source.kind !== 'finite') {
    x.mantissa = 0n;
    x.exponent = 0;
    return 0;
  }
  const shift = x.precision - source.precision;
  if (shift >= 0) {
    x.mantissa = source.mantissa << BigInt(shift);
    return 0;
  }
  const divisor = 1n << BigInt(-shift),
    remainder = source.mantissa % divisor;
  let rounded = source.mantissa / divisor;
  const up = 2n * remainder > divisor || (2n * remainder === divisor && (rounded & 1n) !== 0n);
  if (up) rounded++;
  if (rounded === 1n << BigInt(x.precision)) {
    rounded >>= 1n;
    x.exponent++;
  }
  x.mantissa = rounded;
  return remainder === 0n ? 0 : source.sign * (up ? 1 : -1);
}
