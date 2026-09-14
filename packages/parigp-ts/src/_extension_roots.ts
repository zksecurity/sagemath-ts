/** Native distinct-root split parts via Frobenius and unscaled polynomial GCD.
 * @see Deviation: PARI extension root-count and derivative adapters
 */
import { type ExtensionPolynomial as P, validateExtensionInputs } from './_extension_polynomial.js';
import { extensionField, trimExtension } from './_extension_field.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { extensionFrobenius } from './_extension_frobenius.js';
import { extensionPolynomial } from './_extension_polynomial.js';
import { extensionGcd } from './_extension_gcd.js';
import { FpX_red } from './ffinit.js';
export function extensionSplitPart(mode: 0 | 1 | 2, f: P, T: bigint[] | bigint, p: bigint): P {
  f = trimExtension(f.map((c) => (typeof c === 'bigint' ? c : trimPolynomial(c))));
  if (Array.isArray(T)) T = trimPolynomial(T);
  const modulus = p < 0n ? -p : p;
  if (mode === 0 && modulus > 0n && modulus < 1n << 64n) {
    const word = trimExtension(f.map((c) => FpX_red(typeof c === 'bigint' ? [c] : c, modulus)));
    const r = extensionSplitPart(1, word, FpX_red(T as bigint[], modulus), modulus);
    return r.map((c) => ((c as bigint[]).length < 2 ? ((c as bigint[])[0] ?? 0n) : c));
  }
  if (mode === 1) {
    if (p < 1n || p >= 1n << 64n) throw new RangeError('modulus must be a positive word integer');
    if ([...(T as bigint[]), ...f.flat()].some((c) => c < 0n || c >= p))
      throw new RangeError('word polynomial coefficients must be reduced');
  } else if (mode === 2 && ((T as bigint) < 0n || f.some((c) => (c as bigint) < 0n)))
    throw new RangeError('polynomial bits must be nonnegative');
  if (f.length <= 2) return f;
  validateExtensionInputs(mode, p, T, f);
  f = extensionPolynomial(mode, 2, p, T, f, []);
  const X: P = mode === 1 ? [[], [1n]] : [0n, 1n];
  const z = extensionField(mode, T, p).psub(extensionFrobenius(mode, false, [], f, T, p), X);
  return extensionGcd(mode, 0, p, T, z, f) as P;
}
