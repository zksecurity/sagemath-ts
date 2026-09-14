import { _fmpz_poly_resultant_kernels as k, _fmpz_poly_gcd } from '../fmpz_poly/gcd.js';
import { _fmpz_poly_divrem } from '../fmpz_poly/divrem.js';
import { _fmpz_poly_xgcd } from '../fmpz_poly/xgcd.js';
import { _fmpq_poly_gcd } from './gcd.js';

export type DenseRationalPolynomial = [bigint[], bigint];

/** Native rational XGCD from numerator arrays and positive common denominators.
 * @see Deviation: Polynomial Integer and Rational Extended GCD
 */
export function _fmpq_poly_xgcd(
  a: readonly bigint[],
  denA: bigint,
  b: readonly bigint[],
  denB: bigint
): [DenseRationalPolynomial, DenseRationalPolynomial, DenseRationalPolynomial] {
  if (denA <= 0n || denB <= 0n) throw new RangeError('input denominators must be positive');
  let [A, dA] = canonical(a, denA),
    [B, dB] = canonical(b, denB);
  const swapped = A.length < B.length;
  if (swapped) {
    [A, B] = [B, A];
    [dA, dB] = [dB, dA];
  }
  const finish = (
    g: DenseRationalPolynomial,
    s: DenseRationalPolynomial,
    t: DenseRationalPolynomial
  ): [DenseRationalPolynomial, DenseRationalPolynomial, DenseRationalPolynomial] =>
    swapped ? [g, t, s] : [g, s, t];
  if (!A.length)
    return [
      [[], 1n],
      [[], 1n],
      [[], 1n],
    ];
  if (!B.length) return finish(_fmpq_poly_gcd(A, []), canonical([dA], A[A.length - 1]!), [[], 1n]);
  if (B.length === 1) return finish([[1n], 1n], [[], 1n], canonical([dB], B[0]!));
  const cA = k.content(A),
    cB = k.content(B);
  const primA = A.map((c) => c / cA),
    primB = B.map((c) => c / cB);
  const G = _fmpz_poly_gcd(primA, primB);
  const C = G.length > 1 ? exact(primA, G) : primA;
  const D = G.length > 1 ? exact(primB, G) : primB;
  const [r, S, T] = _fmpz_poly_xgcd(C, D);
  const denG = G[G.length - 1]!;
  return finish(
    [G, denG],
    canonical(
      S.map((c) => c * dA),
      cA * r * denG
    ),
    canonical(
      T.map((c) => c * dB),
      cB * r * denG
    )
  );
}
function canonical(a: readonly bigint[], den: bigint): DenseRationalPolynomial {
  const A = k.normalized(a);
  if (!A.length) return [[], 1n];
  let d = k.content([...A, den]);
  if (den < 0n) d = -d;
  return [A.map((c) => c / d), den / d];
}
function exact(a: bigint[], b: bigint[]): bigint[] {
  const result = _fmpz_poly_divrem(a, b, true);
  if (!result || result[1].length) throw new Error('non-exact primitive quotient in rational XGCD');
  return result[0];
}
