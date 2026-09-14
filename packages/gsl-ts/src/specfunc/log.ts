/** GSL specfunc/log.c: value with the native error handler disabled.
 * @see Deviation: Polynomial Modular Powers
 */
export function gsl_sf_log(x: number): number {
  return x <= 0 ? Number.NaN : Math.log(x);
}
