import type { mpfr_t, mpfr_rnd_t } from './types.js';
import { rounding } from './round_raw.js';
/** MPFR 4.2.1 set_d.c: extract IEEE significand, normalize, then set4 rounding.
 * @see Deviation: Gaussian Parameter Representation (nearest-even scope and flags).
 */
export function mpfr_set_d(x: mpfr_t, value: number, mode: mpfr_rnd_t = 'RNDN'): number {
  rounding(mode);
  if (Number.isNaN(value)) {
    x.kind = 'nan';
    x.mantissa = 0n;
    x.exponent = 0;
    return 0;
  }
  const bits = new DataView(new ArrayBuffer(8));
  bits.setFloat64(0, value, false);
  const raw = bits.getBigUint64(0, false);
  x.sign = raw >> 63n ? -1 : 1;
  x.mantissa = 0n;
  x.exponent = 0;
  if (value === 0) {
    x.kind = 'zero';
    return 0;
  }
  if (!Number.isFinite(value)) {
    x.kind = 'inf';
    return 0;
  }
  const exp = Number((raw >> 52n) & 2047n);
  const mantissa = (raw & ((1n << 52n) - 1n)) | (exp ? 1n << 52n : 0n);
  const length = mantissa.toString(2).length;
  x.exponent = (exp || 1) - 1075 + length;
  x.kind = 'finite';
  const shift = x.precision - length;
  if (shift >= 0) {
    x.mantissa = mantissa << BigInt(shift);
    return 0;
  }
  const denominator = 1n << BigInt(-shift),
    remainder = mantissa % denominator;
  let rounded = mantissa / denominator;
  const up =
    2n * remainder > denominator || (2n * remainder === denominator && (rounded & 1n) !== 0n);
  if (up) rounded++;
  if (rounded === 1n << BigInt(x.precision)) {
    rounded >>= 1n;
    x.exponent++;
  }
  x.mantissa = rounded;
  return remainder === 0n ? 0 : x.sign * (up ? 1 : -1);
}
