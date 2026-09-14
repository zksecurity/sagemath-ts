import type { mpfr_t } from './types.js';
/** MPFR 4.2.1 isinteger.c: inspect exponent and fractional significand bits. */
export function mpfr_integer_p(x: mpfr_t): boolean {
  if (x.kind !== 'finite') return x.kind === 'zero';
  if (x.exponent <= 0) return false;
  if (x.exponent >= x.precision) return true;
  return (x.mantissa & ((1n << BigInt(x.precision - x.exponent)) - 1n)) === 0n;
}
