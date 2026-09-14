import type { mpfr_t } from './types.js';
import { nearest } from './round_raw.js';

/** Full-digit, exact-base-power branch of strtofr.c, with an exact quotient
 * remainder replacing the native interval/Ziv test. Input bounds keep the
 * full powers finite; huge exponents never trigger an unbounded allocation.
 */
export function setRatio(x: mpfr_t, n: bigint, d: bigint): void {
  if (n === 0n) {
    x.kind = 'zero';
    x.mantissa = 0n;
    x.exponent = 0;
    return;
  }
  let e = n.toString(2).length - d.toString(2).length;
  if (e >= 0 ? n < d << BigInt(e) : n << BigInt(-e) < d) e--;
  e++;
  const shift = x.precision - e;
  let mantissa = shift >= 0 ? nearest(n << BigInt(shift), d) : nearest(n, d << BigInt(-shift));
  if (mantissa === 1n << BigInt(x.precision)) {
    mantissa >>= 1n;
    e++;
  }
  x.kind = 'finite';
  x.mantissa = mantissa;
  x.exponent = e;
}
