import type { mpfr_rnd_t } from './types.js';

export function unsupported(message: string): never {
  throw Object.assign(new Error('MPFR_NOT_IMPLEMENTED: ' + message), {
    name: 'NotImplementedError',
  });
}

export function rounding(mode: mpfr_rnd_t): void {
  if (mode !== 'RNDN') unsupported('rounding mode ' + mode);
}

export function nearest(n: bigint, d: bigint): bigint {
  const q = n / d,
    r = n % d;
  return q + (2n * r > d || (2n * r === d && (q & 1n) !== 0n) ? 1n : 0n);
}
