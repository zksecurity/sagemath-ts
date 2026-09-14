/** PARI basemath/polarit3.c mixed coefficient polynomial operations. */
import { extensionPolynomial, type ExtensionPolynomial } from './_extension_polynomial.js';
/** Native normalization preserves integer/polynomial tags and forces integer one.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
export function FpXQX_normalize(
  x: ExtensionPolynomial,
  T: bigint[],
  p: bigint
): ExtensionPolynomial {
  return extensionPolynomial(0, 3, p, T, x);
}

import { ZX_mul, ZX_rem } from './ZX.js';
import { ZX_resultant } from './galconj.js';
import { trimPolynomial } from './_polynomial_packing.js';
import {
  type RationalPair,
  type RationalPolynomialData,
  primitivePolynomial,
  rationalPair,
  rationalProduct,
  scaledPolynomial,
} from './_rational_polynomial.js';
/** Native QXQ_mul with a monic integral modulus and rational contents.
 * @see Deviation: PARI rational trace and norm adapters
 */
export function QXQ_mul(
  x: RationalPolynomialData,
  y: RationalPolynomialData,
  T: bigint[]
): RationalPolynomialData {
  const [a, ca] = primitivePolynomial(x),
    [b, cb] = primitivePolynomial(y);
  return scaledPolynomial(ZX_rem(ZX_mul(a, b), T), rationalProduct(ca, cb));
}
/** Native QXQ_norm, with integral nonzero modulus, primitive content and resultant.
 * @see Deviation: PARI rational trace and norm adapters
 */
export function QXQ_norm(x: RationalPolynomialData, T: bigint[]): RationalPair {
  const [a, c] = primitivePolynomial(x);
  if (!a.length) return [0n, 1n];
  T = trimPolynomial(T);
  const n = BigInt(T.length - 1),
    r = ZX_resultant(T, a);
  return rationalPair(r * c[0] ** n, c[1] ** n * T.at(-1)! ** BigInt(a.length - 1));
}
