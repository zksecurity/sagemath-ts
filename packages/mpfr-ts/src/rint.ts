import type { mpfr_t, mpfr_rnd_t } from './types.js';
import { unsupported } from './round_raw.js';
import { mpfr_integer_p } from './isinteger.js';
/** MPFR 4.2.1 rint.c: round directly on the integer/target-precision grid.
 * The inexact magnitude is 2 for nonintegral input, 1 for an unrepresentable
 * integer, 0 for exact conversion. RNDNA is the private round() tie-away mode.
 * @see Deviation: Native Real Predicates and Integer Rounding
 */
export function mpfr_rint(x: mpfr_t, y: mpfr_t, mode: mpfr_rnd_t | 'RNDNA' = 'RNDN'): number {
  if (!['RNDN', 'RNDZ', 'RNDU', 'RNDD', 'RNDA', 'RNDNA'].includes(mode))
    unsupported('mpfr_rint rounding mode ' + mode);
  const source = { ...y }; // The original supports in-place rounding.
  if (source.kind !== 'finite') {
    x.kind = source.kind;
    x.mantissa = 0n;
    x.exponent = 0;
    if (source.kind !== 'nan') x.sign = source.sign;
    return 0;
  }
  if (source.exponent <= 0) {
    // rint.c handles |u| < 1 without allocating exponent-sized discarded bits.
    const away =
      mode === 'RNDA' ||
      (mode === 'RNDU' && source.sign > 0) ||
      (mode === 'RNDD' && source.sign < 0);
    const nearestUp =
      (mode === 'RNDN' || mode === 'RNDNA') &&
      source.exponent === 0 &&
      (mode === 'RNDNA' || (source.mantissa & (source.mantissa - 1n)) !== 0n);
    const up = away || nearestUp;
    x.sign = source.sign;
    x.kind = up ? 'finite' : 'zero';
    x.exponent = up ? 1 : 0;
    x.mantissa = up ? 1n << BigInt(x.precision - 1) : 0n;
    return source.sign * (up ? 2 : -2);
  }
  const integral = mpfr_integer_p(source);
  const step = Math.max(0, source.exponent - x.precision);
  const shift = step - source.exponent + source.precision;
  let quotient: bigint,
    remainder = 0n,
    divisor = 1n;
  if (shift <= 0) quotient = source.mantissa << BigInt(-shift);
  else {
    divisor = 1n << BigInt(shift);
    quotient = source.mantissa / divisor;
    remainder = source.mantissa % divisor;
  }
  let up = false;
  if (remainder !== 0n) {
    if (mode === 'RNDA') up = true;
    else if (mode === 'RNDU') up = source.sign > 0;
    else if (mode === 'RNDD') up = source.sign < 0;
    else if (mode === 'RNDN' || mode === 'RNDNA')
      up =
        2n * remainder > divisor ||
        (2n * remainder === divisor && (mode === 'RNDNA' || (quotient & 1n) !== 0n));
  }
  // Normalize the retained bits directly, as rint.c does. Work stays bounded
  // by the source/destination precision even when the exponent is enormous.
  const rounded = quotient + (up ? 1n : 0n);
  x.sign = source.sign;
  if (rounded === 0n) {
    x.kind = 'zero';
    x.mantissa = 0n;
    x.exponent = 0;
  } else {
    const length = rounded.toString(2).length;
    const padding = x.precision - length;
    x.kind = 'finite';
    x.exponent = length + step;
    x.mantissa = padding >= 0 ? rounded << BigInt(padding) : rounded >> BigInt(-padding);
  }
  return remainder === 0n ? 0 : source.sign * (up ? 1 : -1) * (integral ? 1 : 2);
}
export function mpfr_roundeven(x: mpfr_t, y: mpfr_t): number {
  return mpfr_rint(x, y, 'RNDN');
}
export function mpfr_round(x: mpfr_t, y: mpfr_t): number {
  return mpfr_rint(x, y, 'RNDNA');
}
export function mpfr_trunc(x: mpfr_t, y: mpfr_t): number {
  return mpfr_rint(x, y, 'RNDZ');
}
export function mpfr_ceil(x: mpfr_t, y: mpfr_t): number {
  return mpfr_rint(x, y, 'RNDU');
}
export function mpfr_floor(x: mpfr_t, y: mpfr_t): number {
  return mpfr_rint(x, y, 'RNDD');
}
