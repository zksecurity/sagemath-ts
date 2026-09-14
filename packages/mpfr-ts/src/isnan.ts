import type { mpfr_t } from './types.js';
/** MPFR 4.2.1 isnan.c; global exception flags are not exposed. */
export function mpfr_nan_p(x: mpfr_t): boolean {
  return x.kind === 'nan';
}
