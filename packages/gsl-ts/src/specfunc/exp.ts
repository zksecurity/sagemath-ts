/** GSL specfunc/exp.c: native normal-range underflow/overflow guards.
 * @see Deviation: Polynomial Modular Powers
 */
export function gsl_sf_exp(x: number): number {
  if (x > 7.0978271289338397e2) return Number.POSITIVE_INFINITY;
  if (x < -7.0839641853226408e2) return 0;
  return Math.exp(x);
}
