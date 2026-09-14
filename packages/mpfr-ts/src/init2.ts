import type { mpfr_t } from './types.js';
import { unsupported } from './round_raw.js';

export function mpfr_init2(precision: number): mpfr_t {
  if (!Number.isSafeInteger(precision) || precision < 1)
    throw new RangeError('invalid MPFR precision');
  if (precision > 4096) unsupported('precision above 4096 bits');
  return { precision, sign: 1, kind: 'nan', mantissa: 0n, exponent: 0 };
}
