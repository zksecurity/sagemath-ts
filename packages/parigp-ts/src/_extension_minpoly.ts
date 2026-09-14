/** Native Shoup transposed-power projection over a prime-field extension.
 * @see Deviation: PARI extension minimal-polynomial adapters
 */
import {
  extensionPolynomial as arithmetic,
  extensionTruncatedPolynomial,
  validateExtensionInputs,
  type ExtensionPolynomial as P,
  type ExtensionCoefficient as C,
} from './_extension_polynomial.js';
import { extensionField, trimExtension as trim } from './_extension_field.js';
import {
  extensionDivision,
  extensionGetRed,
  type ExtensionModulus as M,
} from './_extension_division.js';
import { extensionQuotient } from './_extension_quotient.js';
import { extensionComposition } from './_extension_composition.js';
import { extensionGcd } from './_extension_gcd.js';
import { extensionProjection } from './_extension_projection.js';
import { trimPolynomial } from './_polynomial_packing.js';

export function extensionMinimalPolynomial(
  mode: 0 | 1,
  x: P,
  S: M,
  T: bigint[],
  p: bigint,
  innerInverse?: bigint[]
): P {
  const canonical = (a: P) => trim(a.map((c) => (Array.isArray(c) ? trimPolynomial(c) : c)));
  T = trimPolynomial(T);
  x = canonical(x);
  const raw = canonical(Array.isArray(S) ? S : S.polynomial),
    n = raw.length - 1;
  validateExtensionInputs(mode, p, T, x, raw);
  if (n < 1) throw new RangeError('minimal polynomial modulus must have positive degree');
  if (x.length > n) throw new RangeError('power polynomial exceeds modulus degree');
  if (Array.isArray(S)) S = raw;
  const f = extensionField(mode, T, p, innerInverse);
  const lowmul = (a: P, b: P, length: number) =>
    extensionTruncatedPolynomial(mode, false, length, T, p, a, b, innerInverse);
  const dot = (a: P, b: P) => extensionProjection(mode, 1, 0, T, p, a, b, innerInverse) as C;
  const reverse = (a: P, length: number) =>
    trim(Array.from({ length }, (_, i) => a[length - 1 - i] ?? f.z));
  const shift = (a: P, k: number): P =>
    !a.length ? [] : k < 0 ? a.slice(-k) : [...Array<C>(k).fill(f.z), ...a];
  type Transposed = [P, P, P];
  // Unlike automorphism powers, native minpoly prepares S with the original T.
  // The power-table call prepares its own inner reciprocal without replacing T here.
  S = extensionGetRed(mode, S, T, p, innerInverse);
  const h = Array.isArray(S) ? undefined : S.inverse;
  const transposedInit = (tau: P): Transposed => {
    const ft = reverse(raw, n + 1),
      bt = reverse(tau, n);
    const bht = h
      ? lowmul(bt, h, n - 1)
      : reverse(extensionDivision(mode, 4, p, T, shift(tau, n - 1), S, innerInverse) as P, n - 1);
    return [bt, bht, ft];
  };
  const mul = (a: P, b: P) => arithmetic(mode, 0, p, T, a, b, innerInverse);
  const transposedMultiply = ([bt, bht, ft]: Transposed, a: P): P => {
    if (!a.length) return [];
    const t2 = shift(mul(bt, a), 1 - n);
    if (!bht.length) return t2;
    const t1 = shift(mul(ft, a), -n);
    return f.psub(t2, shift(lowmul(t1, bht, n - 1), 1));
  };
  const powers = extensionQuotient(
    mode,
    6,
    p,
    BigInt(Math.floor(Math.sqrt(2 * n))),
    T,
    x,
    [],
    S,
    innerInverse
  ) as P[];
  let g: P = [f.u],
    tau: P = [f.u];
  while (tau.length) {
    if (g.length - 1 === n) {
      tau = [f.u];
      g = [f.u];
    }
    let v = extensionProjection(mode, 0, n, T, p) as P;
    v = transposedMultiply(transposedInit(tau), v);
    const m = 2 * (n - g.length + 1),
      k = Math.floor(Math.sqrt(m)),
      tr = transposedInit(powers[k]!);
    const c: P = Array<C>(m).fill(f.z);
    for (let i = 0; i < m; i += k) {
      for (let j = 0; j < Math.min(m - i, k); j++) c[m - 1 - i - j] = dot(v, powers[j]!);
      v = transposedMultiply(tr, v);
    }
    // polxn_FpXX and polxn_FlxX both use polynomial coefficient tags, including zero.
    const monomial: P = [...Array.from({ length: m }, () => [] as bigint[]), [1n]];
    const matrix = extensionGcd(mode, 2, p, T, monomial, trim(c), innerInverse) as [[P, P], [P, P]];
    const next = matrix[1][1];
    if (next.length < 2) continue;
    g = mul(g, next);
    tau = extensionQuotient(
      mode,
      0,
      p,
      0n,
      T,
      tau,
      extensionComposition(mode, 1, p, T, next, [], powers, S, innerInverse),
      S,
      innerInverse
    ) as P;
  }
  return arithmetic(mode, 3, p, T, g, [], innerInverse);
}
