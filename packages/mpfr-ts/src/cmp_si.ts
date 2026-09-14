import type { mpfr_t } from './types.js';
/** MPFR 4.2.1 cmp_si.c, exponent-zero branch. No global erange flag is exposed. @see Deviation: Native Real Predicates and Integer Rounding */
export function mpfr_cmp_si(x: mpfr_t, n: bigint): number {
  const sign = n < 0n ? -1 : 1;
  if (x.kind === 'nan') return 0;
  if (x.kind === 'inf') return x.sign;
  if (x.kind === 'zero') return n === 0n ? 0 : -sign;
  if (x.sign !== sign || n === 0n) return x.sign;
  const magnitude = n < 0n ? -n : n,
    length = magnitude.toString(2).length;
  if (x.exponent !== length) return x.sign * (x.exponent > length ? 1 : -1);
  const shift = x.precision - length;
  const left = shift >= 0 ? x.mantissa : x.mantissa << BigInt(-shift);
  const right = shift >= 0 ? magnitude << BigInt(shift) : magnitude;
  return left === right ? 0 : x.sign * (left > right ? 1 : -1);
}
