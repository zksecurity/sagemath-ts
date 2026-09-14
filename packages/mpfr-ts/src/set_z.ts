import type { mpfr_t, mpfr_rnd_t } from './types.js';
import { rounding } from './round_raw.js';

/** MPFR 4.2.1 set_z.c -> set_z_2exp.c with exponent zero: normalize the
 * integer, retain the target precision and round using the discarded bits.
 * @see Deviation: Native Real Numeric Construction (nearest-even scope).
 */
export function mpfr_set_z(x: mpfr_t, value: bigint, mode: mpfr_rnd_t = 'RNDN'): number {
  rounding(mode);
  x.sign = value < 0n ? -1 : 1;
  x.mantissa = 0n;
  x.exponent = 0;
  if (value === 0n) {
    x.kind = 'zero';
    return 0;
  }
  const magnitude = value < 0n ? -value : value;
  const length = magnitude.toString(2).length;
  x.kind = 'finite';
  x.exponent = length;
  const shift = x.precision - length;
  if (shift >= 0) {
    x.mantissa = magnitude << BigInt(shift);
    return 0;
  }
  const divisor = 1n << BigInt(-shift),
    remainder = magnitude % divisor;
  let rounded = magnitude / divisor;
  const up = 2n * remainder > divisor || (2n * remainder === divisor && (rounded & 1n) !== 0n);
  if (up) rounded++;
  if (rounded === 1n << BigInt(x.precision)) {
    rounded >>= 1n;
    x.exponent++;
  }
  x.mantissa = rounded;
  return remainder === 0n ? 0 : x.sign * (up ? 1 : -1);
}
