/** @see Deviation: Polynomial Roots and Truncated Series */
/** fmpq_poly/inv_series_newton.c: reversed-division basecase, then high-half Newton. */
import { _fmpz_poly_inv_series } from '../fmpz_poly/inv_series.js';
import { _fmpz_poly_mullow } from '../fmpz_poly/mullow.js';
import { _fmpz_poly_resultant_kernels as k } from '../fmpz_poly/gcd.js';
export function _fmpq_poly_inv_series_newton(
  a: readonly bigint[],
  den: bigint,
  n: number
): [bigint[], bigint] {
  if (!Number.isSafeInteger(n) || n <= 0) throw new RangeError('precision must be positive');
  if (den <= 0n) throw new RangeError('input denominator must be positive');
  const A = k.normalized(a.slice(0, n));
  if (!A.length || !A[0]) throw new RangeError('constant coefficient must be nonzero');
  const canonical = (v: bigint[], d: bigint): [bigint[], bigint] => {
    const g = k.content([...v, d]) * (d < 0n ? -1n : 1n);
    return [k.normalized(v.map((c) => c / g)), d / g];
  };
  if ((A[0] === 1n || A[0] === -1n) && den === 1n) return [_fmpz_poly_inv_series(A, n), 1n];
  if (A.length === 1) return canonical([den], A[0]!);
  const sizes = [n];
  while (sizes[sizes.length - 1]! > 24) sizes.push(Math.ceil(sizes[sizes.length - 1]! / 2));
  let precision = sizes.pop()!;
  // Reverse Q and divide x^(precision+len(Q)-2) by it. This is the bounded
  // pseudo-division basecase used by native rational polynomial division.
  const B = A.slice(0, precision).reverse(),
    R = Array<bigint>(precision + B.length - 1).fill(0n);
  R[R.length - 1] = 1n;
  let quotient = Array<bigint>(precision).fill(0n),
    scale = 1n;
  const lead = B[B.length - 1]!;
  for (let i = precision - 1; i >= 0; i--) {
    const c = R[i + B.length - 1]!;
    if (c % lead) {
      quotient = quotient.map((q) => q * lead);
      quotient[i] = c;
      for (let j = 0; j < R.length; j++) R[j] = R[j]! * lead;
      scale *= lead;
    } else quotient[i] = c / lead;
    for (let j = 0; j < B.length - 1; j++) R[i + j] = R[i + j]! - quotient[i]! * B[j]!;
    R[i + B.length - 1] = 0n;
  }
  let [out, d] = canonical(
    quotient.reverse().map((c) => c * den),
    scale
  );
  while (sizes.length) {
    const next = sizes.pop()!,
      product = _fmpz_poly_mullow(A.slice(0, next), out, next),
      wden = den * d;
    const high = _fmpz_poly_mullow(out, product.slice(precision), next - precision);
    const low = Array.from({ length: precision }, (_, i) => (out[i] ?? 0n) * wden);
    [out, d] = canonical(
      low.concat(Array.from({ length: next - precision }, (_, i) => -(high[i] ?? 0n))),
      d * wden
    );
    precision = next;
  }
  return [out, d];
}
