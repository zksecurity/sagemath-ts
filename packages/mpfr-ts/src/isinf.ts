import type { mpfr_t } from './types.js';
/** MPFR 4.2.1 isinf.c; global exception flags are not exposed. */
export function mpfr_inf_p(x: mpfr_t): boolean {
  return x.kind === 'inf';
}
