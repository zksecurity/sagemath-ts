import type { mpfr_t } from './types.js';
/** MPFR 4.2.1 sgn.c; global exception flags are not exposed. @see Deviation: Native Real Predicates and Integer Rounding */
export function mpfr_sgn(x: mpfr_t): number {
  return x.kind === 'zero' || x.kind === 'nan' ? 0 : x.sign;
}
