/** Native extension random projection and delayed-reduction dot products.
 * @see Deviation: PARI extension projection adapters
 */
import { random_FpX } from './FpX.js';
import { random_Flx, Flx_mul } from './Flx.js';
import { ZX_mul } from './ZX.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { extensionField, trimExtension } from './_extension_field.js';
import {
  validateExtensionInputs,
  type ExtensionPolynomial as P,
  type ExtensionCoefficient as C,
} from './_extension_polynomial.js';

export function extensionProjection(
  mode: 0 | 1,
  op: 0 | 1,
  length: number,
  T: bigint[],
  p: bigint,
  a: P = [],
  b: P = [],
  innerInverse?: bigint[]
): P | C {
  T = trimPolynomial(T);
  if (op === 0) {
    if (!Number.isSafeInteger(length) || length < 0)
      throw new RangeError('length must be nonnegative');
    if (mode === 1 && (p < 0n || p >= 1n << 64n))
      throw new RangeError('modulus must be an unsigned word integer');
    if (mode === 1 && T.some((c) => c < 0n || c >= 1n << 64n))
      throw new RangeError('modulus coefficients must be unsigned words');
    const degree = T.length - 1;
    if (length && degree < 0) throw new RangeError('coefficient length must be nonnegative');
    return trimExtension(
      Array.from({ length }, () => (mode === 1 ? random_Flx(degree, p) : random_FpX(degree, p)))
    );
  }
  validateExtensionInputs(mode, p, T, a, b);
  a = trimExtension(a.map((c) => (Array.isArray(c) ? trimPolynomial(c) : c)));
  b = trimExtension(b.map((c) => (Array.isArray(c) ? trimPolynomial(c) : c)));
  const lengthDot = Math.min(a.length, b.length);
  const f = extensionField(mode, T, p, innerInverse);
  if (mode === 1) {
    let c: bigint[] = [];
    for (let i = 0; i < lengthDot; i++)
      c = f.add(c, Flx_mul(a[i] as bigint[], b[i] as bigint[], p)) as bigint[];
    return lengthDot ? f.red(c) : [];
  }
  let c: C = 0n;
  for (let i = 0; i < lengthDot; i++) {
    const x = a[i]!,
      y = b[i]!;
    const product =
      typeof x === 'bigint' && typeof y === 'bigint'
        ? x * y
        : ZX_mul(typeof x === 'bigint' ? [x] : x, typeof y === 'bigint' ? [y] : y);
    if (typeof c === 'bigint' && typeof product === 'bigint') c += product;
    else {
      const aa = typeof c === 'bigint' ? [c] : c,
        bb = typeof product === 'bigint' ? [product] : product;
      c = trimPolynomial(
        Array.from(
          { length: Math.max(aa.length, bb.length) },
          (_, i) => (aa[i] ?? 0n) + (bb[i] ?? 0n)
        )
      );
    }
  }
  return lengthDot ? f.red(c) : 0n;
}
