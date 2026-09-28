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


import { FpX_sub, FpX_Fp_mul, FpX_red } from './ffinit.js';
import { Flx_mul, Flx_sqr, Flx_divrem, Flx_resultant, Flv_polint } from './Flx.js';
import { FpX_resultant, FpX_eval, FpV_polint } from './FpX.js';
import { Fp_powu } from './arith1.js';

type Bivariate = bigint[][];
function trimBivariate(x: Bivariate): Bivariate {
  let length = x.length;
  while (length && !x[length - 1]!.length) length--;
  return x.slice(0, length);
}
function wordPolynomialPower(x: bigint[], n: number, p: bigint): bigint[] {
  let result = [1n];
  while (n) {
    if (n % 2) result = Flx_mul(result, x, p);
    n = Math.floor(n / 2);
    if (n) x = Flx_sqr(x, p);
  }
  return result;
}
/** polarit3.c:1779: pseudo-remainder, including skipped-degree scaling. */
function wordPseudoRemainder(u: Bivariate, v: Bivariate, p: bigint): Bivariate {
  const degree = v.length - 1,
    leading = v[degree]!;
  let r = u.map((c) => c.slice()),
    remaining = u.length - v.length + 1;
  while (r.length >= v.length) {
    const top = r.at(-1)!,
      shift = r.length - v.length;
    r = r.map((c) => Flx_mul(c, leading, p));
    for (let j = 0; j < v.length; j++)
      r[j + shift] = FpX_sub(r[j + shift]!, Flx_mul(top, v[j]!, p), p);
    r = trimBivariate(r);
    remaining--;
  }
  if (!r.length || !remaining) return r;
  const scale = wordPolynomialPower(leading, remaining, p);
  return r.map((c) => Flx_mul(c, scale, p));
}
/** polarit3.c:1822: exact x^n/y^(n-1) through Lazard's binary schedule. */
function wordLazard(x: bigint[], y: bigint[], n: number, p: bigint): bigint[] {
  if (n === 1) return x;
  let bit = 2 ** Math.floor(Math.log2(n)),
    c = x;
  n -= bit;
  while (bit > 1) {
    bit /= 2;
    c = Flx_divrem(Flx_sqr(c, p), y, p)[0];
    if (n >= bit) {
      c = Flx_divrem(Flx_mul(c, x, p), y, p)[0];
      n -= bit;
    }
  }
  return c;
}
/** polarit3.c:1840: polynomial-coefficient subresultant sequence. */
function wordSubresultant(u: Bivariate, v: Bivariate, p: bigint): bigint[] {
  let sign = 1;
  if (u.length < v.length) {
    if ((u.length - 1) % 2 && (v.length - 1) % 2) sign = -sign;
    [u, v] = [v, u];
  }
  if (!v.length) return [];
  if (v.length === 1) return wordPolynomialPower(v[0]!, u.length - 1, p);
  let g = [1n],
    h = [1n];
  for (;;) {
    const r = wordPseudoRemainder(u, v, p);
    if (!r.length) return [];
    const du = u.length - 1,
      dv = v.length - 1,
      delta = du - dv;
    u = v;
    let divisor = g;
    g = u.at(-1)!;
    if (delta === 1) {
      divisor = Flx_mul(h, divisor, p);
      h = g;
    } else if (delta > 1) {
      divisor = Flx_mul(wordPolynomialPower(h, delta, p), divisor, p);
      h = wordLazard(g, h, delta, p);
    }
    if (du % 2 && dv % 2) sign = -sign;
    v = r.map((c) => Flx_divrem(c, divisor, p)[0]);
    if (r.length === 1) {
      let z = v[0]!;
      if (dv > 1) z = wordLazard(z, h, dv, p);
      return sign < 0 ? FpX_Fp_mul(z, p - 1n, p) : z;
    }
  }
}
function resultantInterpolation(
  T: bigint[],
  Q: Bivariate,
  bound: number,
  p: bigint,
  word: boolean
): bigint[] {
  const points: bigint[] = [];
  for (let i = 0, n = 1n; i < bound; i += 2, n++) points.push(n, p - n);
  if (bound % 2 === 0) points.push(0n);
  const values = points.map((point) => {
    const evaluated = trimPolynomial(Q.map((c) => FpX_eval(c, point, p)));
    let value = (word ? Flx_resultant : FpX_resultant)(T, evaluated, p);
    const drop = Q.length - evaluated.length,
      leading = T.at(-1)!;
    if (drop && leading !== 1n) value = (value * Fp_powu(leading, BigInt(drop), p)) % p;
    return value;
  });
  return (word ? Flv_polint : FpV_polint)(points, values, p);
}
/** Native word bivariate resultant, with Q[eliminated degree][retained degree].
 * @see Deviation: PARI bivariate polynomial storage
 */
export function Flx_FlxY_resultant(T: bigint[], Q: Bivariate, p: bigint): bigint[] {
  T = trimPolynomial(T);
  Q = trimBivariate(Q.map(trimPolynomial));
  const retainedDegree = Math.max(-1, ...Q.map((c) => c.length - 1));
  const bound = (T.length - 1) * retainedDegree;
  return bound < 0 || BigInt(bound) >= p
    ? wordSubresultant(
        T.map((c) => (c === 0n ? [] : [c])),
        Q,
        p
      )
    : resultantInterpolation(T, Q, bound, p, true);
}
/** Native arbitrary-prime branch, preserving the existing port's coefficient layout.
 * @see Deviation: PARI bivariate polynomial storage
 */
export function FpX_FpXY_resultant(T: bigint[], Q: Bivariate, p: bigint): bigint[] {
  if (p < 1n << 64n)
    return Flx_FlxY_resultant(
      FpX_red(T, p),
      Q.map((c) => FpX_red(c, p)),
      p
    );
  T = trimPolynomial(T);
  Q = trimBivariate(Q.map(trimPolynomial));
  const bound = (T.length - 1) * Math.max(-1, ...Q.map((c) => c.length - 1));
  return resultantInterpolation(T, Q, bound, p, false);
}
