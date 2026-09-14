import type { mpfr_t } from './types.js';
/** MPFR 4.2.1 isnum.c; global exception flags are not exposed. */
export function mpfr_number_p(x: mpfr_t): boolean {
  return x.kind === 'finite' || x.kind === 'zero';
}
