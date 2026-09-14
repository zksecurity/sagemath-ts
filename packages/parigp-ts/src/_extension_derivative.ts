/** Native outer derivatives without reduction by an inner field modulus.
 * @see Deviation: PARI extension root-count and derivative adapters
 */
import { type ExtensionPolynomial as P } from './_extension_polynomial.js';
import { trimExtension } from './_extension_field.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { FpX_red, FpX_Fp_mul } from './ffinit.js';
import { PariError } from './errors.js';
export function extensionDerivative(word: boolean, x: P, p: bigint): P {
  if (word) {
    if (p < 1n || p >= 1n << 64n) throw new RangeError('modulus must be a positive word integer');
    if (x.flat().some((c) => c < 0n || c >= p))
      throw new RangeError('word polynomial coefficients must be reduced');
  }
  x = trimExtension(x.map((c) => (typeof c === 'bigint' ? c : trimPolynomial(c))));
  if (x.length < 2) return [];
  const modulus = p < 0n ? -p : p;
  return trimExtension(
    x.slice(1).map((c, j) => {
      if (!modulus && typeof c === 'bigint' && c === 0n) return 0n;
      if (!modulus)
        throw new PariError(
          `impossible inverse in ${typeof c === 'bigint' ? 'dvmdii' : 'umodui'}: 0.`
        );
      const i = BigInt(j + 1);
      return typeof c === 'bigint' ? (FpX_red([c * i], p)[0] ?? 0n) : FpX_Fp_mul(c, i % modulus, p);
    })
  );
}
