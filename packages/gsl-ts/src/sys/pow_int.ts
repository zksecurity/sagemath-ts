/** GSL sys/pow_int.c; IEEE-754 operation order is intentional.
 * @see Deviation: Polynomial Modular Powers
 */
export function gsl_pow_int(x: number, n: number): number {
  if (!Number.isInteger(n) || n < -(2 ** 31) || n >= 2 ** 31)
    throw new RangeError('exponent must fit a signed int');
  if (n < 0) {
    x = 1 / x;
    n = -n;
  }
  let value = 1;
  do {
    if (n % 2) value *= x;
    n = Math.floor(n / 2);
    x *= x;
  } while (n);
  return value;
}
