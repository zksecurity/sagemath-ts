import type { mpfr_t, mpfr_rnd_t } from './types.js';
import { rounding } from './round_raw.js';
import { mpfr_integer_p } from './isinteger.js';
import { mpfr_set } from './set.js';
/** MPFR 4.2.1 frac.c: isolate the fractional significand bits and normalize.
 * @see Deviation: Native Real Fractional Parts and Comparisons
 */
export function mpfr_frac(x: mpfr_t, y: mpfr_t, mode: mpfr_rnd_t = 'RNDN'): number {
  rounding(mode);
  const source = { ...y };
  if (source.kind === 'nan') {
    x.kind = 'nan';
    x.mantissa = 0n;
    x.exponent = 0;
    return 0;
  }
  if (source.kind === 'inf' || mpfr_integer_p(source)) {
    x.kind = 'zero';
    x.sign = source.sign;
    x.mantissa = 0n;
    x.exponent = 0;
    return 0;
  }
  if (source.exponent <= 0) return mpfr_set(x, source, mode);
  const bits = source.precision - source.exponent;
  const fraction = source.mantissa & ((1n << BigInt(bits)) - 1n);
  const precision = fraction.toString(2).length;
  const status = mpfr_set(
    x,
    {
      kind: 'finite',
      sign: source.sign,
      mantissa: fraction,
      precision,
      exponent: precision - bits,
    },
    mode
  );
  // frac.c's same-buffer limb path calls mpfr_round_raw, which marks exact
  // nearest-even ties with +/-2. The temporary-buffer path uses set's +/-1.
  const padding = (64 - (source.precision % 64)) % 64;
  const fractionLimb = Math.floor((precision + padding - 1) / 64);
  if (status !== 0 && Math.floor((x.precision - 1) / 64) >= fractionLimb) {
    const discarded = precision - x.precision;
    if (discarded > 0) {
      const divisor = 1n << BigInt(discarded);
      if (2n * (fraction % divisor) === divisor) return 2 * status;
    }
  }
  return status;
}
