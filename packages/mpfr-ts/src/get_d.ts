import type { mpfr_t, mpfr_rnd_t } from './types.js';
import { rounding, nearest } from './round_raw.js';

/** get_d.c: singular values, under/overflow and normal/subnormal rounding. */
export function mpfr_get_d(x: mpfr_t, mode: mpfr_rnd_t = 'RNDN'): number {
  rounding(mode);
  if (x.kind === 'nan') return NaN;
  if (x.kind === 'inf') return x.sign * Infinity;
  if (x.kind === 'zero') return x.sign < 0 ? -0 : 0;
  if (x.exponent > 1024) return x.sign * Infinity;
  if (x.exponent < -1074) return x.sign < 0 ? -0 : 0;
  const quantum = Math.max(x.exponent - 53, -1074),
    shift = quantum - (x.exponent - x.precision);
  const q = shift >= 0 ? nearest(x.mantissa, 1n << BigInt(shift)) : x.mantissa << BigInt(-shift);
  return x.sign * Number(q) * 2 ** quantum;
}
