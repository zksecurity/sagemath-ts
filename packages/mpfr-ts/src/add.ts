import type { mpfr_t, mpfr_rnd_t } from './types.js';
import { rounding } from './round_raw.js';
import { mpfr_set } from './set.js';
import { roundedDyadic, singular } from './arithmetic.js';

/** MPFR add.c: singular dispatch, aligned addition/subtraction and rounding.
 * @see Deviation: Native MPFR binary arithmetic
 */
export function mpfr_add(a: mpfr_t, b: mpfr_t, c: mpfr_t, mode: mpfr_rnd_t = 'RNDN'): number {
  rounding(mode);
  if (b.kind === 'nan' || c.kind === 'nan') return singular(a, 'nan');
  if (b.kind === 'inf')
    return c.kind === 'inf' && c.sign !== b.sign ? singular(a, 'nan') : singular(a, 'inf', b.sign);
  if (c.kind === 'inf') return singular(a, 'inf', c.sign);
  if (b.kind === 'zero' && c.kind === 'zero')
    return singular(a, 'zero', b.sign === -1 && c.sign === -1 ? -1 : 1);
  if (b.kind === 'zero') return mpfr_set(a, c, mode);
  if (c.kind === 'zero') return mpfr_set(a, b, mode);
  // sub1.c preserves MPFR_RNDRAW_EVEN's +/-2 on a midpoint rounded
  // toward zero when the subtrahend is too small to affect that rounding.
  // Capture this before writing a, which may alias either operand.
  let midpointStatus = 0;
  if (b.sign !== c.sign) {
    const dominant = b.exponent >= c.exponent ? b : c;
    const smaller = dominant === b ? c : b;
    const discarded = dominant.precision - a.precision;
    if (
      discarded > 0 &&
      dominant.exponent - smaller.exponent >= Math.max(a.precision, dominant.precision) + 2
    ) {
      const unit = 1n << BigInt(discarded);
      if (dominant.mantissa % unit === unit / 2n && (dominant.mantissa / unit) % 2n === 0n)
        midpointStatus = -2 * dominant.sign;
    }
  }
  const working = Math.max(a.precision, b.precision, c.precision) + 4;
  // Keep a sticky bit when exponents are far apart. Storage then depends on
  // precision, even for million-bit exponent gaps; cancellation stays exact.
  const exponent = Math.max(
    Math.min(b.exponent - b.precision, c.exponent - c.precision),
    Math.max(b.exponent, c.exponent) - 2 * working
  );
  const aligned = (x: mpfr_t): bigint => {
    const shift = x.exponent - x.precision - exponent;
    if (shift >= 0) return BigInt(x.sign) * (x.mantissa << BigInt(shift));
    if (-shift >= x.precision) return BigInt(x.sign);
    const unit = 1n << BigInt(-shift);
    return BigInt(x.sign) * (x.mantissa / unit + (x.mantissa % unit === 0n ? 0n : 1n));
  };
  const value = aligned(b) + aligned(c);
  if (value === 0n) return singular(a, 'zero', 1);
  const status = roundedDyadic(a, value < 0n ? -value : value, exponent, value < 0n ? -1 : 1);
  return midpointStatus || status;
}
