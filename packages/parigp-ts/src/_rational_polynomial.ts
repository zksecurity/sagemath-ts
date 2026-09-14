/** Exact QQ coefficient/storage adapters for PARI's rational polynomial kernels. */
import { gcd } from './ff.js';
import { PariError } from './errors.js';
import { trimPolynomial } from './_polynomial_packing.js';
export type RationalPair = [bigint, bigint];
export type RationalPolynomialData = [bigint[], bigint];
export function rationalPair(n: bigint, d: bigint): RationalPair {
  // Qdivii handles numerator one through inversion; other numerators use dvmdii.
  if (!d) throw new PariError(`impossible inverse in ${n === 1n ? 'gdiv' : 'dvmdii'}: 0`);
  const c = gcd(n, d) * (d < 0n ? -1n : 1n);
  return [n / c, d / c];
}
export function rationalProduct(a: RationalPair, b: RationalPair): RationalPair {
  const g = gcd(a[0], b[1]),
    h = gcd(b[0], a[1]);
  return rationalPair((a[0] / g) * (b[0] / h), (a[1] / h) * (b[1] / g));
}
export function rationalDifference(a: RationalPair, b: RationalPair): RationalPair {
  const g = gcd(a[1], b[1]);
  return rationalPair(a[0] * (b[1] / g) - b[0] * (a[1] / g), (a[1] / g) * b[1]);
}
export function rationalQuotient(a: RationalPair, b: RationalPair): RationalPair {
  return rationalProduct(a, rationalPair(b[1], b[0]));
}
export function normalizedPolynomial([input, den]: RationalPolynomialData): RationalPolynomialData {
  if (!den) throw new PariError('impossible inverse in gdiv: 0');
  const a = trimPolynomial(input),
    g = a.reduce((c, x) => gcd(c, x), den < 0n ? -den : den) * (den < 0n ? -1n : 1n);
  return [a.map((x) => x / g), den / g];
}
/** Q_primitive_part: integral primitive coefficients and their rational content. */
export function primitivePolynomial(input: RationalPolynomialData): [bigint[], RationalPair] {
  const [a, d] = normalizedPolynomial(input),
    c = a.reduce((g, x) => gcd(g, x), 0n);
  return c ? [a.map((x) => x / c), [c, d]] : [[], [0n, 1n]];
}
export function scaledPolynomial(a: bigint[], q: RationalPair): RationalPolynomialData {
  return normalizedPolynomial([a.map((x) => x * q[0]), q[1]]);
}
export function polynomialFromFractions(a: RationalPair[]): RationalPolynomialData {
  const den = a.reduce((d, c) => (d / gcd(d, c[1])) * c[1], 1n);
  return normalizedPolynomial([a.map((c) => c[0] * (den / c[1])), den]);
}
