import type { mpfr_t, mpfr_rnd_t } from './types.js';
import { mpfr_rint } from './rint.js';
/** MPFR 4.2.1 get_z.c: return [integer, native inexact status].
 * Singular inputs produce zero; MPFR's global erange flag is not exposed.
 * @see Deviation: Native Real Predicates and Integer Rounding
 */
export function mpfr_get_z(x: mpfr_t, mode: mpfr_rnd_t = 'RNDN'): [bigint, number] {
  if (x.kind !== 'finite') return [0n, 0];
  // Native get_z's temporary precision holds every integral bit. It may exceed
  // the public real precision limit because the result is an unbounded BigInt.
  const rounded = { ...x, precision: Math.max(1, x.exponent) };
  const status = mpfr_rint(rounded, x, mode);
  if (rounded.kind === 'zero') return [0n, status];
  const shift = rounded.exponent - rounded.precision;
  const magnitude =
    shift >= 0 ? rounded.mantissa << BigInt(shift) : rounded.mantissa >> BigInt(-shift);
  return [BigInt(rounded.sign) * magnitude, status];
}
