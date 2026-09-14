import { Flx_gcd, Flx_normalize, Flx_sqr, random_Flx } from './Flx.js';
import { FpX_gcd, FpX_normalize, FpX_sqr, random_FpX } from './FpX.js';
import { brent_kung_optpow } from './RgX.js';
import { polynomialDistinctDegrees } from './_polynomial_ddf.js';
import { polynomialMinimalPolynomial } from './_polynomial_minpoly.js';
import { type PolynomialQuotient, polynomialQuotient } from './_polynomial_quotient.js';
import { gen_pow_i, gen_powu_i } from './bb_group.js';
import { PariError } from './errors.js';
/** PARI FpX_factor.c Shoup equal-degree splitting and native random-state use. */
import { type FpX, FpX_add, FpX_sub } from './ffinit.js';
import { isPrime } from './ifactor.js';
import { random_Fl, randomi } from './random.js';

export function polynomialFactorShoup(Tp: FpX, p: bigint, word: boolean): FpX[] {
  const gcd = (a: FpX, b: FpX) => (word ? Flx_gcd(a, b, p) : FpX_gcd(a, b, p));
  const norm = (a: FpX) => (word ? Flx_normalize(a, p) : FpX_normalize(a, p));
  const randomScalar = () => (word ? random_Fl(p) : randomi(p));
  const randomPoly = (n: number) => (word ? random_Flx(n, p) : random_FpX(n, p));
  const context = (f: FpX) => polynomialQuotient(f, p, word);
  const check = (name: string, ct: number) => {
    if (!word && ct === 10 && !isPrime(p))
      throw new PariError(`not a prime number in ${name}: ${p}.`);
  };
  // All factor-splitting exponents are (p-1)/2, hence positive for odd p.
  const power = (T: PolynomialQuotient, x: FpX): FpX => {
    const e = p >> 1n;
    if (e === 1n) return x.slice();
    return gen_pow_i(
      x,
      e,
      (a) => T.reduce(word ? Flx_sqr(a, p) : FpX_sqr(a, p)),
      (a, b) => T.reduce(T.multiply(a, b))
    );
  };
  const evaluate = (T: PolynomialQuotient, x: FpX, t: FpX) =>
    T.evaluate(x, T.powers(t, Math.floor(Math.sqrt(x.length - 1))));
  const trace = (T: PolynomialQuotient, XP: FpX, g: FpX, d: number): FpX => {
    type Pair = [FpX, FpX];
    const multiply = (a: Pair, b: Pair): Pair => {
      const v = T.powers(a[0], brent_kung_optpow(Math.max(b[0].length, b[1].length) - 1, 2, 1));
      return [T.evaluate(b[0], v), FpX_add(a[1], T.evaluate(b[1], v), p)];
    };
    return gen_powu_i<Pair>([XP, g], BigInt(d), (a) => multiply(a, a), multiply)[1];
  };
  const simple = (f: FpX, XP: FpX, d: number): FpX[] => {
    const n = f.length - 1;
    if (n === d) return [f];
    const T = context(f);
    XP = T.reduce(XP);
    let part: FpX = [];
    let ct = 0;
    for (;;) {
      const t = trace(T, XP, randomPoly(n), d);
      if (!t.length) continue;
      for (let i = 0; i < 10; i++) {
        const R = power(T, FpX_add(t, [randomScalar()], p));
        part = gcd(FpX_sub(R, [1n], p), f);
        if (part.length > 1 && part.length < f.length) break;
      }
      if (part.length > 1 && part.length < f.length) break;
      check('FpX_edf_simple', ++ct);
    }
    part = norm(part);
    const other = T.divide(f, part)[0];
    return [...simple(part, XP, d), ...simple(other, XP, d)];
  };
  const recurse = (T: PolynomialQuotient, XP: FpX, hp: FpX, t: FpX, d: number): FpX[] => {
    const h = context(hp);
    t = T.reduce(t);
    let u1: FpX = [];
    let ct = 0;
    do {
      const R = power(h, [randomScalar(), 1n]);
      u1 = gcd(FpX_sub(R, [1n], p), hp);
      check('FpX_edf_rec', ++ct);
    } while (u1.length === 1 || u1.length === hp.length);
    const f1 = norm(gcd(evaluate(T, u1, t), T.T));
    const u2 = h.divide(hp, u1)[0];
    const f2 = T.divide(T.T, f1)[0];
    const branch = (f: FpX, u: FpX): FpX[] =>
      u.length === 2
        ? word && f.length - 1 !== d
          ? edf(f, XP, d)
          : [f]
        : recurse(context(f), XP, u, t, d);
    return [...branch(f1, u1), ...branch(f2, u2)];
  };
  const edf = (f: FpX, XP: FpX, d: number): FpX[] => {
    const n = f.length - 1;
    const r = n / d;
    if (word && r === 1) return [f];
    const T = context(f);
    XP = T.reduce(XP);
    let h: FpX;
    let t: FpX;
    let ct = 0;
    do {
      t = trace(T, XP, randomPoly(n), d);
      h = polynomialMinimalPolynomial(t, f, p, word, T);
      check('FpX_edf', ++ct);
    } while (word ? h.length <= 2 : h.length - 1 !== r);
    return recurse(T, XP, h, t, d);
  };
  const T = context(Tp);
  const XP = T.power([0n, 1n]);
  const parts = polynomialDistinctDegrees(T, XP);
  const out: FpX[] = [];
  const e = p.toString(2).length - 1;
  for (const [i, part] of parts.entries()) {
    const d = i + 1;
    const n = part.length - 1;
    if (!n) continue;
    const f = norm(part);
    if (n === d) {
      out.push(f);
      continue;
    }
    out.push(...(n / d <= e * Math.floor(Math.log2(e)) ? edf(f, XP, d) : simple(f, XP, d)));
  }
  return out;
}
