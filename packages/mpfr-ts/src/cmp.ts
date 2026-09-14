import type { mpfr_t } from './types.js';
/** MPFR 4.2.1 cmp.c: compare signs, exponents, then aligned significands.
 * NaN compares as zero; MPFR's global erange flag is not exposed.
 * @see Deviation: Native Real Fractional Parts and Comparisons
 */
export function mpfr_cmp(x: mpfr_t, y: mpfr_t): number {
  if (x.kind === 'nan' || y.kind === 'nan') return 0;
  if (x.kind === 'inf') return y.kind === 'inf' && x.sign === y.sign ? 0 : x.sign;
  if (y.kind === 'inf') return -y.sign;
  if (x.kind === 'zero') return y.kind === 'zero' ? 0 : -y.sign;
  if (y.kind === 'zero') return x.sign;
  if (x.sign !== y.sign) return x.sign;
  if (x.exponent !== y.exponent) return x.sign * (x.exponent > y.exponent ? 1 : -1);
  const precision = Math.max(x.precision, y.precision);
  const left = x.mantissa << BigInt(precision - x.precision),
    right = y.mantissa << BigInt(precision - y.precision);
  return left === right ? 0 : x.sign * (left > right ? 1 : -1);
}
