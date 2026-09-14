/** @see Deviation: Polynomial Roots and Truncated Series */
/** fmpz_poly/inv_series.c: binomial/basecase and high-half Newton inversion. */
import { _fmpz_poly_mullow } from './mullow.js';
import { _fmpz_poly_resultant_kernels as k } from './gcd.js';
export function _fmpz_poly_inv_series(a: readonly bigint[], n: number): bigint[] {
  if (!Number.isSafeInteger(n) || n <= 0) throw new RangeError('precision must be positive');
  const A = k.normalized(a.slice(0, n));
  if (!A.length || (A[0] !== 1n && A[0] !== -1n))
    throw new RangeError('constant coefficient must be 1 or -1');
  const sizes = [n];
  if (A.length >= 64)
    while (sizes[sizes.length - 1]! >= 64) sizes.push(Math.ceil(sizes[sizes.length - 1]! / 2));
  let precision = sizes.pop()!,
    out = Array<bigint>(precision).fill(0n);
  out[0] = A[0]!;
  const Q = A.slice(0, precision),
    step = Q.length - 1;
  if (step > 0 && Q.slice(1, -1).every((c) => c === 0n)) {
    const ratio = -Q[0]! * Q[step]!;
    for (let i = step; i < precision; i += step) out[i] = out[i - step]! * ratio;
  } else {
    for (let i = 1; i < precision; i++) {
      let sum = 0n;
      for (let j = 1; j < Math.min(i + 1, Q.length); j++) sum += Q[j]! * out[i - j]!;
      out[i] = -Q[0]! * sum;
    }
  }
  while (sizes.length) {
    const next = sizes.pop()!,
      product = _fmpz_poly_mullow(A.slice(0, next), out, next);
    const high = _fmpz_poly_mullow(out, product.slice(precision), next - precision);
    out = out.concat(Array.from({ length: next - precision }, (_, i) => -(high[i] ?? 0n)));
    precision = next;
  }
  return k.normalized(out);
}
