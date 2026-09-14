import type { mpfr_t } from './types.js';
import { mpfr_set } from './set.js';
import { unsupported } from './round_raw.js';

/** Normalize an exact dyadic result, round within Sage's configured MPFR exponent range. */
export function roundedDyadic(
  a: mpfr_t,
  magnitude: bigint,
  exponent: number,
  sign: 1 | -1
): number {
  const bits = magnitude.toString(2).length;
  const source: mpfr_t = {
    kind: 'finite',
    sign,
    mantissa: magnitude,
    precision: bits,
    exponent: exponent + bits,
  };
  // Sage initializes MPFR's exponent range to +/- (2^62-1). Every safe
  // JavaScript integer exponent is inside that range; larger exponents need
  // an exponent representation wider than the current mpfr_t interface.
  if (!Number.isSafeInteger(source.exponent)) unsupported('exponent outside safe integer range');
  return mpfr_set(a, source);
}

/** Set a singular result; NaN keeps the destination sign, as MPFR_SET_NAN does. */
export function singular(a: mpfr_t, kind: 'nan' | 'inf' | 'zero', sign?: 1 | -1): number {
  a.kind = kind;
  a.mantissa = 0n;
  a.exponent = 0;
  if (sign !== undefined) a.sign = sign;
  return 0;
}
