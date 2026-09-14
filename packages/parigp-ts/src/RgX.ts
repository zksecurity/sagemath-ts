/** PARI RgX.c:31 optimal Brent–Kung table size, retaining the first tied minimum. */
export function brent_kung_optpow(d: number, n: number, m: number): number {
  if (![d, n, m].every(Number.isSafeInteger))
    throw new RangeError('parameters must be safe integers');
  let best = 1,
    cost = BigInt(n) * BigInt(d - 1);
  for (let size = 2; size <= d; size++) {
    const next = BigInt(m) * BigInt(size - 1) + BigInt(n) * BigInt(Math.floor((d - 1) / size));
    if (next < cost) {
      best = size;
      cost = next;
    }
  }
  return best;
}

import { PariError } from './errors.js';
import { QX_mul, QX_ZX_rem } from './ZX.js';
import { QXQ_mul, QXQ_norm } from './polarit3.js';
import {
  type RationalPair,
  type RationalPolynomialData,
  normalizedPolynomial,
  rationalPair,
  rationalProduct,
  rationalDifference,
  rationalQuotient,
  polynomialFromFractions,
  primitivePolynomial,
} from './_rational_polynomial.js';
/** Native RgX_rem over QQ: monic-integer fast path or generic quotient elimination.
 * @see Deviation: PARI rational trace and norm adapters
 */
export function RgX_rem(
  input: RationalPolynomialData,
  divisor: RationalPolynomialData
): RationalPolynomialData {
  const x = normalizedPolynomial(input),
    y = normalizedPolynomial(divisor),
    dx = x[0].length - 1,
    dy = y[0].length - 1;
  if (y[1] === 1n && y[0].at(-1) === 1n) return QX_ZX_rem(x, y[0]);
  if (dy < 0) throw new PariError('impossible inverse in RgX_divrem: 0');
  if (dx < dy) return x;
  if (!dy) return [[], 1n];
  const dz = dx - dy,
    q: Array<RationalPair> = Array.from({ length: dz + 1 }, () => [0n, 1n]),
    coeff = (i: number): RationalPair => rationalPair(y[0][i]!, y[1]),
    lead = coeff(dy);
  for (let i = dx; i >= dy; i--) {
    let v = rationalPair(x[0][i]!, x[1]);
    for (let j = i - dy + 1; j <= i && j <= dz; j++)
      v = rationalDifference(v, rationalProduct(q[j]!, coeff(i - j)));
    q[i - dy] = rationalQuotient(v, lead);
  }
  const r: Array<RationalPair> = [];
  for (let i = 0; i < dy; i++) {
    let v = rationalPair(x[0][i]!, x[1]);
    for (let j = 0; j <= i && j <= dz; j++)
      v = rationalDifference(v, rationalProduct(q[j]!, coeff(i - j)));
    r.push(v);
  }
  return polynomialFromFractions(r);
}
/** Native RgXQ_mul QQ dispatch; nonintegral/nonmonic moduli use generic reduction.
 * @see Deviation: PARI rational trace and norm adapters
 */
export function RgXQ_mul(
  x: RationalPolynomialData,
  y: RationalPolynomialData,
  modulus: RationalPolynomialData
): RationalPolynomialData {
  const T = normalizedPolynomial(modulus);
  return T[1] === 1n && T[0].at(-1) === 1n ? QXQ_mul(x, y, T[0]) : RgX_rem(QX_mul(x, y), T);
}
/** Native derivative/product/remainder trace formula over QQ.
 * @see Deviation: PARI rational trace and norm adapters
 */
export function RgXQ_trace(
  x: RationalPolynomialData,
  modulus: RationalPolynomialData
): RationalPair {
  const T = normalizedPolynomial(modulus),
    n = T[0].length - 1;
  if (!n) return [0n, 1n];
  const derivative: RationalPolynomialData = [T[0].slice(1).map((c, i) => c * BigInt(i + 1)), T[1]],
    z = RgXQ_mul(x, derivative, T);
  return rationalPair((z[0][n - 1] ?? 0n) * T[1], z[1] * T[0][n]!);
}
/** Native resultant/leading-coefficient norm formula over QQ.
 * @see Deviation: PARI rational trace and norm adapters
 */
export function RgXQ_norm(
  x: RationalPolynomialData,
  modulus: RationalPolynomialData
): RationalPair {
  const [T] = primitivePolynomial(modulus);
  if (T.length === 1) return [1n, 1n];
  return QXQ_norm(x, T);
}
