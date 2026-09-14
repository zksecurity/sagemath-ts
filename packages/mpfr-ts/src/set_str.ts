import type { mpfr_t, mpfr_rnd_t } from './types.js';
import { unsupported, rounding } from './round_raw.js';
import { setRatio } from './strtofr.js';

/** set_str.c wraps strtofr and returns -1 unless every character was consumed.
 * Like the original, malformed nonempty input can still change the destination.
 */
export function mpfr_set_str(
  x: mpfr_t,
  text: string,
  base = 10,
  mode: mpfr_rnd_t = 'RNDN'
): number {
  rounding(mode);
  if (base !== 10) unsupported('non-decimal string conversion');
  if (text === '') return -1;
  const prefix =
    /^[\t\n\r\f\v ]*([+-]?)(?:((?:@nan@|nan)(?:\([A-Za-z0-9_]*\))?)|(@inf@|infinity|inf)|((?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?))/i.exec(
      text
    );
  x.sign = prefix?.[1] === '-' ? -1 : 1;
  x.mantissa = 0n;
  x.exponent = 0;
  if (!prefix) {
    x.kind = 'zero';
    return -1;
  }
  const status = prefix[0].length === text.length ? 0 : -1;
  if (prefix[2]) {
    x.kind = 'nan';
    // strtofr resets the destination to positive zero before parsing NaN.
    x.sign = 1;
    return status;
  }
  if (prefix[3]) {
    x.kind = 'inf';
    return status;
  }
  const [mantissa, exponent = '0'] = prefix[4]!.split(/[eE]/);
  const decimalPlaces = mantissa!.includes('.') ? mantissa!.length - mantissa!.indexOf('.') - 1 : 0;
  const digits = mantissa!.replace('.', '').replace(/^0+/, '');
  if (!digits) {
    x.kind = 'zero';
    return status;
  }
  if (digits.length > 4096) unsupported('decimal mantissa longer than 4096 digits');
  const scale = BigInt(exponent) - BigInt(decimalPlaces);
  if (scale > 4096n || scale < -4096n) unsupported('decimal exponent outside -4096..4096');
  const n = BigInt(digits);
  setRatio(x, scale >= 0n ? n * 10n ** scale : n, scale < 0n ? 10n ** -scale : 1n);
  return status;
}
