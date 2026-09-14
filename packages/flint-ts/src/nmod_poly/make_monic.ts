/** Nonzero dense-array port of nmod_poly/make_monic.c.
 * Coefficients are reduced modulo n; the leading coefficient must be a unit.
 * @see Deviation: Dense FLINT polynomial GCD kernels
 */
export function _nmod_poly_make_monic(a: readonly bigint[], n: bigint): bigint[] {
  if (n < 2n) throw new RangeError('modulus must be at least 2');
  const mod = (x: bigint) => ((x % n) + n) % n;
  const A = a.map(mod);
  while (A.length && A[A.length - 1] === 0n) A.pop();
  if (!A.length) throw new RangeError('polynomial division by zero');
  let r = n,
    s = A[A.length - 1]!,
    x = 0n,
    y = 1n;
  while (s) {
    const q = r / s;
    [r, s] = [s, r - q * s];
    [x, y] = [y, x - q * y];
  }
  if (r !== 1n) throw new RangeError('leading coefficient is not invertible');
  return A.map((c) => mod(c * x));
}
