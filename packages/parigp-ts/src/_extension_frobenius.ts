/** Native prime/extension Frobenius phases and automorphism dispatch.
 * @see Deviation: PARI extension Frobenius adapters
 */
import { type ExtensionPolynomial as P, validateExtensionInputs } from './_extension_polynomial.js';
import { type ExtensionModulus as M, extensionGetRed } from './_extension_division.js';
import { trimExtension as trim } from './_extension_field.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { polynomialQuotient } from './_polynomial_quotient.js';
import { extensionQuotient } from './_extension_quotient.js';
import { extensionAutomorphism } from './_extension_automorphism.js';
import { FpXQ_pow, FpX_red } from './ffinit.js';
import { Flx_rem, Flx_sqr } from './Flx.js';
import { F2x_degree, F2x_Frobenius } from './F2x.js';
import { gen_powu_i } from './bb_group.js';
import { isqrt } from './ifactor.js';

export function wordFrobenius(T: bigint[], p: bigint, innerInverse?: bigint[]): bigint[] {
  T = trimPolynomial(T);
  const x = [0n, 1n];
  // Flxq_powu returns its direct square before preparing an unprovided cache.
  if (p === 2n && innerInverse === undefined) return Flx_rem(Flx_sqr(x, p), T, p);
  const c = polynomialQuotient(T, p, true, innerInverse);
  return gen_powu_i(
    x,
    p,
    (a) => c.reduce(Flx_sqr(a, p)),
    (a, b) => c.reduce(c.multiply(a, b))
  );
}

export function extensionFrobenius(
  mode: 0 | 1 | 2,
  half: boolean,
  a: P,
  S: M,
  T: bigint[] | bigint,
  p: bigint
): P {
  const canonical = (a: P) => trim(a.map((c) => (Array.isArray(c) ? trimPolynomial(c) : c)));
  a = canonical(a);
  S = Array.isArray(S) ? canonical(S) : S;
  if (Array.isArray(T)) T = trimPolynomial(T);
  const raw = canonical(Array.isArray(S) ? S : S.polynomial);
  validateExtensionInputs(mode, p, T, raw, a);
  // Native word conversion requires a coefficient modulus of positive degree.
  // A raw leading multiple of p otherwise produces an invalid extension context.
  if (mode === 0 && p < 1n << 64n && FpX_red(T as bigint[], p).length < 2)
    throw new RangeError('extension modulus must have positive degree');
  if (half && mode === 0 && p < 1n << 64n) {
    const convert = (v: P) => trim(v.map((c) => FpX_red(typeof c === 'bigint' ? [c] : c, p)));
    const s = Array.isArray(S)
      ? convert(S)
      : { polynomial: convert(S.polynomial), inverse: convert(S.inverse) };
    const r = extensionFrobenius(1, true, convert(a), s, FpX_red(T as bigint[], p), p);
    return r.map((c) => ((c as bigint[]).length < 2 ? ((c as bigint[])[0] ?? 0n) : c));
  }
  let inverse: bigint[] | undefined;
  if (half) {
    inverse = polynomialQuotient(T as bigint[], p, mode === 1).inverse;
    S = extensionGetRed(mode, S, T, p, inverse);
  }
  const X: P = mode === 1 ? [[], [1n]] : [0n, 1n];
  const xp =
    mode === 2
      ? F2x_Frobenius(T as bigint)
      : mode === 1
        ? wordFrobenius(T as bigint[], p, inverse)
        : inverse === undefined
          ? FpXQ_pow([0n, 1n], p, T as bigint[], p)
          : polynomialQuotient(T as bigint[], p, false, inverse).power([0n, 1n]);
  const Xp = extensionQuotient(
    mode,
    mode === 2 ? 1 : mode === 1 ? 7 : 5,
    p,
    p,
    T,
    X,
    [],
    S,
    inverse
  ) as P;
  const degreeT = mode === 2 ? F2x_degree(T as bigint) : (T as bigint[]).length - 1;
  if (half) {
    const ap = extensionQuotient(mode, mode === 1 ? 7 : 5, p, p >> 1n, T, a, [], S, inverse) as P;
    return extensionAutomorphism(mode, 2, [xp, Xp, ap], BigInt(degreeT), S, T, p, inverse)[2]!;
  }
  const q = mode === 2 ? 1n << BigInt(degreeT) : p ** BigInt(degreeT);
  const exponent = mode === 2 ? degreeT : q.toString(2).length - 1;
  const threshold =
    (BigInt(degreeT).toString(2).length - 1) * Number(isqrt(BigInt(raw.length - 1)));
  if (exponent >= threshold)
    return extensionAutomorphism(mode, 0, [xp, Xp], BigInt(degreeT), S, T, p)[1];
  return extensionQuotient(mode, 5, p, q, T, X, [], S) as P;
}
