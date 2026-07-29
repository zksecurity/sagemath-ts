/**
 * @module sage/rings/number_field/number_field_embeddings
 * @description Archimedean embeddings of a number field: the real and complex
 * roots of its defining polynomial.
 *
 * SageMath computes `K.embeddings(RR/CC)` as `K.defining_polynomial().roots(R)`,
 * which for a rational polynomial delegates to PARI (`polroots` / `polrootsreal`)
 * or, for certified intervals, to `sage.rings.polynomial.complex_roots`.  This
 * file ports both:
 *
 * - the **real** roots by transcribing PARI's Uspensky (Vincent--Collins--Akritas)
 *   root isolation, `reference/pari/src/basemath/rootpol.c:2257` (`usp`),
 *   `:2379` (`ZX_Uspensky`) and `:2600` (`ZX_realroots_irred`), together with
 *   the Descartes counter `X2XP1` (`rootpol.c:1987`), `ZX_rescale2prim`
 *   (`rootpol.c:2241`) and `fujiwara_bound_real` (`rootpol.c:1658`);
 * - the **complex** roots by transcribing SageMath's own certified algorithm,
 *   `reference/sage/src/sage/rings/polynomial/complex_roots.py:154`
 *   (`complex_roots`), `:50` (`interval_roots`), `:93` (`intervals_disjoint`)
 *   and `reference/sage/src/sage/rings/polynomial/refine_root.pyx:27`
 *   (`refine_root`): a floating-point estimate of every root followed by an
 *   interval-arithmetic Newton certification which *proves* that each returned
 *   box contains exactly one root and that the boxes are pairwise disjoint.
 *
 * This file lives in `rings/number_field/` rather than `rings/polynomial/`
 * because the polynomial package is owned by another module; see the note on
 * `pari_nf.ts`, which hosts PARI's `nf` layer for the same reason.
 *
 * @see Deviation: Number Field Kernel Ported Locally Instead of parigp-ts
 */

import { NotImplementedError, ValueError } from '../../errors.js';
import { Rational } from '../rational.js';
import type { ZPoly } from './pari_nf.js';
import { zpDeg, zpDerivative, zpNorm } from './pari_nf.js';

/* ------------------------------------------------------------------ */
/* Dyadic rounding                                                     */
/* ------------------------------------------------------------------ */

function babs(a: bigint): bigint {
  return a < 0n ? -a : a;
}

/** Number of bits of |a| (0 for a = 0). */
export function bitLength(a: bigint): number {
  let x = babs(a);
  if (x === 0n) return 0;
  let n = 0;
  // 64 bits at a time, then bit by bit
  while (x >= 0x10000000000000000n) {
    x >>= 64n;
    n += 64;
  }
  while (x > 0n) {
    x >>= 1n;
    n++;
  }
  return n;
}

/**
 * Round `r` outward (away from zero when `up` is the sign of the desired
 * direction) to a dyadic rational with about `bits` significant bits.
 *
 * `dir = -1` rounds down (towards -infinity), `dir = +1` rounds up.  Rounding
 * *outward* is what keeps interval arithmetic sound.
 */
export function roundDyadic(r: Rational, bits: number, dir: -1 | 1): Rational {
  if (r.isZero()) return Rational.zero();
  const num = r.numerator;
  const den = r.denominator;
  // shift so that the numerator has ~bits bits: value = num/den = (num*2^s)/(den*2^s)
  const e = bitLength(num) - bitLength(den);
  const s = bits - e;
  // n = round(r * 2^s)
  let a = num;
  let b = den;
  if (s >= 0) a <<= BigInt(s);
  else b <<= BigInt(-s);
  // floor / ceil of a/b
  let q = a / b;
  const rem = a - q * b;
  if (rem !== 0n) {
    if (dir === 1) {
      if (a > 0n) q += 1n;
      // a<0: truncation already rounded towards +infinity
    } else {
      if (a < 0n) q -= 1n;
    }
  }
  if (s >= 0) return new Rational(q, 1n << BigInt(s));
  return new Rational(q << BigInt(-s), 1n);
}

/* ------------------------------------------------------------------ */
/* Real intervals over dyadic rationals                                */
/* ------------------------------------------------------------------ */

/** A closed real interval `[lo, hi]` with exact rational endpoints. */
export interface RI {
  lo: Rational;
  hi: Rational;
}

function ri(lo: Rational, hi: Rational): RI {
  return { lo, hi };
}

export function ri_exact(x: Rational): RI {
  return { lo: x, hi: x };
}

function ri_round(x: RI, prec: number): RI {
  return { lo: roundDyadic(x.lo, prec, -1), hi: roundDyadic(x.hi, prec, 1) };
}

function ri_add(a: RI, b: RI, prec: number): RI {
  return ri_round(ri(a.lo.add(b.lo), a.hi.add(b.hi)), prec);
}

function ri_sub(a: RI, b: RI, prec: number): RI {
  return ri_round(ri(a.lo.sub(b.hi), a.hi.sub(b.lo)), prec);
}

function ri_neg(a: RI): RI {
  return ri(a.hi.neg(), a.lo.neg());
}

function ri_mul(a: RI, b: RI, prec: number): RI {
  const p = [a.lo.mul(b.lo), a.lo.mul(b.hi), a.hi.mul(b.lo), a.hi.mul(b.hi)];
  let lo = p[0]!;
  let hi = p[0]!;
  for (const x of p) {
    if (x.lt(lo)) lo = x;
    if (x.gt(hi)) hi = x;
  }
  return ri_round(ri(lo, hi), prec);
}

function ri_contains_zero(a: RI): boolean {
  return a.lo.le(Rational.zero()) && a.hi.ge(Rational.zero());
}

function ri_contains(outer: RI, inner: RI): boolean {
  return outer.lo.le(inner.lo) && outer.hi.ge(inner.hi);
}

function ri_union(a: RI, b: RI): RI {
  return ri(a.lo.lt(b.lo) ? a.lo : b.lo, a.hi.gt(b.hi) ? a.hi : b.hi);
}

function ri_center(a: RI): Rational {
  return a.lo.add(a.hi).div(new Rational(2n));
}

function ri_diameter(a: RI): Rational {
  return a.hi.sub(a.lo);
}

function ri_overlaps(a: RI, b: RI): boolean {
  return a.lo.le(b.hi) && b.lo.le(a.hi);
}

/* ------------------------------------------------------------------ */
/* Complex intervals                                                   */
/* ------------------------------------------------------------------ */

/** A complex rectangle `re x im`, the analogue of Sage's `ComplexIntervalField`. */
export interface CI {
  re: RI;
  im: RI;
}

export function ci(re: RI, im: RI): CI {
  return { re, im };
}

export function ci_exact(re: Rational, im: Rational): CI {
  return { re: ri_exact(re), im: ri_exact(im) };
}

function ci_add(a: CI, b: CI, prec: number): CI {
  return ci(ri_add(a.re, b.re, prec), ri_add(a.im, b.im, prec));
}

function ci_sub(a: CI, b: CI, prec: number): CI {
  return ci(ri_sub(a.re, b.re, prec), ri_sub(a.im, b.im, prec));
}

function ci_mul(a: CI, b: CI, prec: number): CI {
  const re = ri_sub(ri_mul(a.re, b.re, prec), ri_mul(a.im, b.im, prec), prec);
  const im = ri_add(ri_mul(a.re, b.im, prec), ri_mul(a.im, b.re, prec), prec);
  return ci(re, im);
}

/** `1/a`; returns `null` when `a` may be zero. */
function ci_inv(a: CI, prec: number): CI | null {
  // |a|^2 = re^2 + im^2, computed as an interval
  const n = ri_add(ri_mul(a.re, a.re, prec), ri_mul(a.im, a.im, prec), prec);
  if (ri_contains_zero(n)) return null;
  // n > 0 here, so 1/n = [1/hi, 1/lo]
  const inv = ri_round(ri(n.hi.inv(), n.lo.inv()), prec);
  return ci(ri_mul(a.re, inv, prec), ri_mul(ri_neg(a.im), inv, prec));
}

function ci_div(a: CI, b: CI, prec: number): CI | null {
  const bi = ci_inv(b, prec);
  if (bi === null) return null;
  return ci_mul(a, bi, prec);
}

function ci_contains_zero(a: CI): boolean {
  return ri_contains_zero(a.re) && ri_contains_zero(a.im);
}

function ci_contains(outer: CI, inner: CI): boolean {
  return ri_contains(outer.re, inner.re) && ri_contains(outer.im, inner.im);
}

function ci_union(a: CI, b: CI): CI {
  return ci(ri_union(a.re, b.re), ri_union(a.im, b.im));
}

function ci_center(a: CI): CI {
  return ci_exact(ri_center(a.re), ri_center(a.im));
}

function ci_diameter(a: CI): Rational {
  const dr = ri_diameter(a.re);
  const di = ri_diameter(a.im);
  return dr.gt(di) ? dr : di;
}

function ci_overlaps(a: CI, b: CI): boolean {
  return ri_overlaps(a.re, b.re) && ri_overlaps(a.im, b.im);
}

/** Horner evaluation of an integer polynomial on a complex interval. */
function ci_polyeval(p: ZPoly, x: CI, prec: number): CI {
  const d = zpDeg(p);
  if (d < 0) return ci_exact(Rational.zero(), Rational.zero());
  let acc = ci_exact(new Rational(p[d]!), Rational.zero());
  for (let i = d - 1; i >= 0; i--) {
    acc = ci_mul(acc, x, prec);
    acc = ci_add(acc, ci_exact(new Rational(p[i]!), Rational.zero()), prec);
  }
  return acc;
}

/* ------------------------------------------------------------------ */
/* PARI rootpol.c: real root isolation                                 */
/* ------------------------------------------------------------------ */

/** `ZX_Z_translate`: `P(X + t)` by repeated synthetic division. */
export function ZX_Z_translate(P: ZPoly, t: bigint): ZPoly {
  if (t === 0n) return [...P];
  const a = [...P];
  const d = zpDeg(a);
  for (let i = 0; i < d; i++) {
    for (let j = d - 1; j >= i; j--) a[j] = a[j]! + t * a[j + 1]!;
  }
  return a;
}

/** `ZX_unscale2n`: `P(2^n X)`. */
export function ZX_unscale2n(P: ZPoly, n: number): ZPoly {
  if (n === 0) return [...P];
  const out: ZPoly = [];
  for (let i = 0; i < P.length; i++) out.push(P[i]! << BigInt(n * i));
  return out;
}

/** `ZX_z_unscale(P, -1)`: `P(-X)`. */
export function ZX_neg_unscale(P: ZPoly): ZPoly {
  return P.map((c, i) => (i % 2 === 0 ? c : -c));
}

/**
 * `ZX_rescale2prim` (rootpol.c:2241): the primitive part (w.r.t. powers of 2)
 * of `2^deg P * P(X/2)`.
 */
export function ZX_rescale2prim(P: ZPoly): ZPoly {
  const d = zpDeg(P);
  if (d < 0) return [0n];
  const val2 = (x: bigint): number => {
    if (x === 0n) return Number.MAX_SAFE_INTEGER;
    let v = 0;
    let y = babs(x);
    while ((y & 1n) === 0n) {
      y >>= 1n;
      v++;
    }
    return v;
  };
  let v = val2(P[d]!);
  for (let n = 1, i = d - 1; v > n && i >= 0; i--, n++) v = Math.min(v, val2(P[i]!) + n);
  const out: ZPoly = new Array(d + 1).fill(0n);
  const shift = (x: bigint, s: number): bigint => (s >= 0 ? x << BigInt(s) : x >> BigInt(-s));
  out[d] = v ? shift(P[d]!, -v) : P[d]!;
  for (let i = d - 1, n = 1 - v; i >= 0; i--, n++) out[i] = shift(P[i]!, n);
  return out;
}

function bsign(x: bigint): number {
  return x > 0n ? 1 : x < 0n ? -1 : 0;
}

/**
 * `X2XP1` (rootpol.c:1987): the number of sign variations of
 * `(X+1)^deg P * P(1/(X+1))`, i.e. (Descartes) an upper bound with the right
 * parity for the number of roots of `P` in `(0,1)`.  Returns 0, 1, or 2, where
 * 2 means "two or more, unknown".
 */
export function X2XP1(P: ZPoly): 0 | 1 | 2 {
  const dP = zpDeg(P);
  const a = P.slice(0, dP + 1);
  let lim = dP;
  for (let idx = 0; idx < lim; idx++) a[idx + 1] = a[idx + 1]! + a[idx]!;
  let s = -bsign(a[lim]!);
  lim--;
  let nb = 0;
  for (let i = 1; i < dP; i++) {
    const s2 = -bsign(a[0]!);
    let flag = s2 === s;
    for (let idx = 0; idx < lim; idx++) {
      a[idx + 1] = a[idx + 1]! + a[idx]!;
      if (flag) flag = s2 !== bsign(a[idx + 1]!);
    }
    if (s === bsign(a[lim]!)) {
      if (++nb >= 2) return 2;
      s = -s;
    }
    if (flag) return nb as 0 | 1;
    lim--;
  }
  if (lim >= 0 && s === bsign(a[lim]!)) nb++;
  return nb as 0 | 1 | 2;
}

/**
 * `usp` (rootpol.c:2257) with `flag = 0`: isolating intervals, with dyadic
 * endpoints, for the roots of `Q0` in `(0, 1)`.  `Q0` must have no rational
 * root (guaranteed here: it comes from an irreducible polynomial of degree > 1).
 */
export function usp(Q0in: ZPoly): RI[] {
  let Q0 = zpNorm([...Q0in]);
  const sol: RI[] = [];
  const Lc: bigint[] = [0n];
  const Lk: number[] = [0];
  let k = 0;
  let ind = 0;
  let indf = 1;
  let Q = Q0;
  let c = 0n;
  let nb_todo = 1;
  const pow2 = (e: number): Rational => new Rational(1n, 1n << BigInt(e));
  while (nb_todo) {
    const nc = Lc[ind]!;
    if (Lk[ind] === k + 1) {
      Q0 = ZX_rescale2prim(Q0);
      Q = Q0;
      c = 0n;
    }
    if (nc !== c) Q = ZX_Z_translate(Q, nc - c);
    k = Lk[ind]!;
    ind++;
    c = nc;
    nb_todo--;
    const nb = X2XP1(Q);
    if (nb === 1) {
      const h = pow2(k);
      sol.push(ri(new Rational(c).mul(h), new Rational(c + 1n).mul(h)));
    } else if (nb) {
      Lc[indf] = c * 2n;
      Lc[indf + 1] = c * 2n + 1n;
      Lk[indf] = k + 1;
      Lk[indf + 1] = k + 1;
      indf += 2;
      nb_todo += 2;
    }
  }
  return sol;
}

/**
 * An exact upper bound for `log2` of the largest positive real root of `P`
 * (`sign = 1`) or of `-P(-X)` (`sign = -1`).
 *
 * Upstream's `fujiwara_bound_real` (rootpol.c:1658) zeroes the coefficients
 * that cannot contribute and calls `fujiwara_bound` (rootpol.c:1628), which
 * works in C doubles.  We keep the same reduction but replace the double
 * logarithms with exact integer bit lengths, so the value returned is a
 * *proved* upper bound (never smaller than upstream's, at most one bit larger).
 *
 * @see Deviation: Exact bit-length bound instead of PARI's double logarithms
 */
export function fujiwara_bound_real_bits(P: ZPoly, sign: 1 | -1): number {
  const n = zpDeg(P);
  const x = P.slice(0, n + 1);
  let signeven: number;
  let signodd: number;
  if (bsign(x[n]!) > 0) {
    signeven = 1;
    signodd = sign;
  } else {
    signeven = -1;
    signodd = -sign;
  }
  for (let i = 0; i < n; i++) {
    if ((n - i) % 2) {
      if (bsign(x[i]!) === signodd) x[i] = 0n;
    } else {
      if (bsign(x[i]!) === signeven) x[i] = 0n;
    }
  }
  // Fujiwara: log2|z| <= 1 + max_i (log2|c_i| - log2|c_n|)/(n-i), with the
  // constant term getting an extra -1.  Upper bound each log2|c_i| by
  // bitLength(c_i) and lower bound log2|c_n| by bitLength(c_n) - 1.
  const lcBits = bitLength(x[n]!) - 1;
  let Lmax = -Infinity;
  if (x[0] !== 0n) Lmax = (bitLength(x[0]!) - lcBits - 1) / n;
  for (let i = 1; i < n; i++) {
    if (x[i] === 0n) continue;
    const L = (bitLength(x[i]!) - lcBits) / (n - i);
    if (L > Lmax) Lmax = L;
  }
  if (Lmax === -Infinity) return 0;
  return Math.max(0, Math.ceil(Lmax) + 1);
}

/** `RgX_deflate_max`: write `P(x) = Q(x^h)` with `h` maximal; returns `[Q, h]`. */
export function ZX_deflate_max(P: ZPoly): [ZPoly, number] {
  const d = zpDeg(P);
  if (d <= 0) return [[...P], 1];
  let h = 0;
  for (let i = 1; i <= d; i++) {
    if (P[i] !== 0n) {
      h = h === 0 ? i : gcdNum(h, i);
      if (h === 1) return [[...P], 1];
    }
  }
  if (h <= 1) return [[...P], 1];
  const Q: ZPoly = [];
  for (let i = 0; i * h <= d; i++) Q.push(P[i * h]!);
  return [Q, h];
}

function gcdNum(a: number, b: number): number {
  while (b) {
    const t = a % b;
    a = b;
    b = t;
  }
  return a;
}

/**
 * Isolating intervals for all real roots of an irreducible `ZX` of degree > 1.
 *
 * Transcription of `ZX_realroots_irred` (rootpol.c:2600) composed with
 * `ZX_Uspensky(P, ., 1|4, .)` (rootpol.c:2379) in the `flag = 0` (interval)
 * mode: `flag & 4` skips the rational-root search because `P` is irreducible.
 */
export function ZX_realroots_irred_intervals(Pin: ZPoly): RI[] {
  const [P, h] = ZX_deflate_max(zpNorm([...Pin]));
  const dP = zpDeg(P);
  if (dP < 1) return [];
  if (dP === 1) {
    const r = new Rational(-P[0]!, P[1]!);
    return liftDeflation([ri(r, r)], h);
  }
  const out: RI[] = [];
  // b = +infinity branch
  const bp = fujiwara_bound_real_bits(P, 1);
  for (const s of usp(ZX_unscale2n(P, bp))) {
    const m = new Rational(1n << BigInt(bp));
    out.push(ri(s.lo.mul(m), s.hi.mul(m)));
  }
  // a = -infinity branch, only when odd(h) (upstream restricts to [0, oo) when
  // h is even, because the deflated variable is x^h >= 0)
  if (h % 2 === 1) {
    const bm = fujiwara_bound_real_bits(P, -1);
    const neg: RI[] = [];
    for (const s of usp(ZX_unscale2n(ZX_neg_unscale(P), bm))) {
      const m = new Rational(1n << BigInt(bm));
      neg.push(ri(s.hi.mul(m).neg(), s.lo.mul(m).neg()));
    }
    neg.reverse();
    out.unshift(...neg);
  } else {
    // the deflated polynomial is evaluated at x^h > 0 only; upstream passes
    // `gen_0` as the lower bound, i.e. only positive roots are kept.
  }
  out.sort((a, b) => a.lo.cmp(b.lo));
  return liftDeflation(out, h);
}

/** Undo `x -> x^h`: replace each interval by the real `h`-th roots it yields. */
function liftDeflation(sols: RI[], h: number): RI[] {
  if (h === 1) return sols;
  const out: RI[] = [];
  for (const s of sols) {
    if (h % 2 === 0 && s.hi.le(Rational.zero())) continue;
    const r = nthRootInterval(s, h);
    out.push(r);
    if (h % 2 === 0) out.push(ri_neg(r));
  }
  out.sort((a, b) => a.lo.cmp(b.lo));
  return out;
}

/** An enclosing interval for the real `h`-th root of every point of `s`. */
function nthRootInterval(s: RI, h: number): RI {
  return ri(realNthRoot(s.lo, h, -1, 80), realNthRoot(s.hi, h, 1, 80));
}

/** A dyadic lower (`dir=-1`) or upper (`dir=1`) bound for `x^(1/h)`. */
function realNthRoot(x: Rational, h: number, dir: -1 | 1, bits: number): Rational {
  const neg = x.isNegative();
  const a = x.abs();
  // bisection on [0, max(1,a)]
  let lo = Rational.zero();
  let hi = a.gt(Rational.one()) ? a : Rational.one();
  const target = a;
  for (let i = 0; i < bits + 8; i++) {
    const m = lo.add(hi).div(new Rational(2n));
    if (m.pow(BigInt(h)).le(target)) lo = m;
    else hi = m;
    lo = roundDyadic(lo, bits, -1);
    hi = roundDyadic(hi, bits, 1);
  }
  if (neg) {
    // odd h only: (-a)^(1/h) = -(a^(1/h)); bounds swap
    return dir === -1 ? hi.neg() : lo.neg();
  }
  return dir === -1 ? lo : hi;
}

/**
 * Refine an isolating interval of a *simple* real root to `prec` bits by exact
 * bisection.  Upstream refines with `polsolve` (rootpol.c:2139), a bisection
 * followed by Newton in `t_REAL` arithmetic; we keep the bisection and drop the
 * Newton acceleration so that every step stays exact (the endpoints remain
 * dyadic and the sign of `P` at a dyadic point is an exact integer sign).
 *
 * @see Deviation: Exact bisection instead of PARI's floating-point Newton step
 */
export function refineRealRoot(P: ZPoly, iv: RI, prec: number): RI {
  const evalSign = (x: Rational): number => {
    let acc = Rational.zero();
    for (let i = zpDeg(P); i >= 0; i--) acc = acc.mul(x).add(new Rational(P[i]!));
    return acc.isZero() ? 0 : acc.isPositive() ? 1 : -1;
  };
  let lo = iv.lo;
  let hi = iv.hi;
  let slo = evalSign(lo);
  if (slo === 0) return ri(lo, lo);
  const shi = evalSign(hi);
  if (shi === 0) return ri(hi, hi);
  if (slo === shi) {
    throw new ValueError('refineRealRoot: the interval does not isolate a sign change');
  }
  // width target: 2^-prec relative to max(1,|root|)
  const scale = Math.max(1, bitLength(lo.abs().ceil()), bitLength(hi.abs().ceil()));
  const eps = new Rational(1n, 1n << BigInt(prec + scale + 2));
  let guard = 0;
  while (hi.sub(lo).gt(eps) && guard++ < prec + scale + 64) {
    const m = lo.add(hi).div(new Rational(2n));
    const sm = evalSign(m);
    if (sm === 0) return ri(m, m);
    // `slo` is invariant: we only ever move `lo` to a point of the same sign.
    if (sm === slo) lo = m;
    else hi = m;
  }
  return ri(lo, hi);
}

/**
 * All real roots of an irreducible `ZX`, as isolating intervals of width
 * `< 2^-prec`.
 */
export function ZX_realroots_irred(P: ZPoly, prec: number): RI[] {
  const iso = ZX_realroots_irred_intervals(P);
  const [Q, h] = ZX_deflate_max(zpNorm([...P]));
  void Q;
  if (h > 1) {
    // the intervals came out of an nth-root extraction; refine them against P
    return iso.map((s) => refineRealRootFromEnclosure(P, s, prec));
  }
  return iso.map((s) => refineRealRoot(P, s, prec));
}

function refineRealRootFromEnclosure(P: ZPoly, s: RI, prec: number): RI {
  // widen slightly until the endpoints have opposite signs, then bisect
  let lo = s.lo;
  let hi = s.hi;
  const evalR = (x: Rational): Rational => {
    let acc = Rational.zero();
    for (let i = zpDeg(P); i >= 0; i--) acc = acc.mul(x).add(new Rational(P[i]!));
    return acc;
  };
  let w = hi.sub(lo);
  if (w.isZero()) w = new Rational(1n, 1n << 40n);
  for (let i = 0; i < 64; i++) {
    const a = evalR(lo);
    const b = evalR(hi);
    if (a.isZero()) return ri(lo, lo);
    if (b.isZero()) return ri(hi, hi);
    if (a.isPositive() !== b.isPositive()) return refineRealRoot(P, ri(lo, hi), prec);
    lo = lo.sub(w);
    hi = hi.add(w);
    w = w.mul(new Rational(2n));
  }
  throw new ValueError('could not bracket a real root');
}

/* ------------------------------------------------------------------ */
/* Sage complex_roots.py / refine_root.pyx: certified complex roots    */
/* ------------------------------------------------------------------ */

/**
 * `refine_root` (refine_root.pyx:27): interval Newton--Raphson.  Returns an
 * interval that provably contains exactly one root of `p`, or `null`.
 */
export function refine_root(p: ZPoly, pd: ZPoly, irt0: CI, prec: number): CI | null {
  const refinement_steps = 10;
  let irt = irt0;
  let smashed_real = false;
  let smashed_imag = false;
  for (let i = 0; i < refinement_steps; i++) {
    const slope = ci_polyeval(pd, irt, prec);
    if (ci_contains_zero(slope)) return null;
    const center = ci_center(irt);
    const val = ci_polyeval(p, center, prec);
    const q = ci_div(val, slope, prec);
    if (q === null) return null;
    const nirt = ci_sub(center, q, prec);
    if (
      ci_contains(irt, nirt) &&
      (ci_diameter(nirt).mul(new Rational(8n)).ge(ci_diameter(irt)) || i >= 8)
    ) {
      return nirt;
    }
    if (i & 1) {
      irt = nirt;
    } else {
      irt = ci_union(irt, nirt);
      if (i >= 6) {
        const rD = ri_diameter(irt.re);
        const iD = ri_diameter(irt.im);
        const md = rD.gt(iD) ? rD : iD;
        const mdi = ri(md.neg(), md);
        irt = ci(ri_add(irt.re, mdi, prec), ri_add(irt.im, mdi, prec));
      }
    }
    if (!smashed_real && ri_contains_zero(irt.re)) {
      irt = ci(ri_exact(Rational.zero()), irt.im);
      smashed_real = true;
    }
    if (!smashed_imag && ri_contains_zero(irt.im)) {
      irt = ci(irt.re, ri_exact(Rational.zero()));
      smashed_imag = true;
    }
  }
  return null;
}

/** `intervals_disjoint` (complex_roots.py:93) — pairwise disjointness. */
export function intervals_disjoint(intvs: CI[]): boolean {
  for (let a = 0; a < intvs.length; a++) {
    for (let b = a + 1; b < intvs.length; b++) {
      if (ci_overlaps(intvs[a]!, intvs[b]!)) return false;
    }
  }
  return true;
}

/**
 * A floating-point estimate of every complex root, by the Durand--Kerner
 * (Weierstrass) iteration on `p/lc`.  This is the role SageMath's
 * `complex_roots` gives to `cfac.roots()` (PARI `polroots` / numpy): a *guess*
 * that the interval Newton step below has to certify.  Nothing downstream
 * trusts it.
 */
function durandKerner(p: ZPoly): Array<[number, number]> {
  const d = zpDeg(p);
  const lc = p[d]!;
  // scale the variable so the roots are O(1): x = R*y with R a Cauchy bound
  let maxAbs = 0n;
  for (let i = 0; i < d; i++) {
    const a = babs(p[i]!);
    if (a > maxAbs) maxAbs = a;
  }
  const Rbits = Math.max(0, bitLength(maxAbs) - bitLength(babs(lc)) + 1);
  const R = Math.pow(2, Rbits);
  // c[i] = p[i] * R^i / (lc * R^d), a monic polynomial in y
  const c: Array<[number, number]> = [];
  for (let i = 0; i <= d; i++) {
    const num = Number(p[i]!) / Number(lc);
    c.push([num * Math.pow(2, (i - d) * Rbits), 0]);
  }
  const cadd = (a: [number, number], b: [number, number]): [number, number] => [
    a[0] + b[0],
    a[1] + b[1],
  ];
  const csub = (a: [number, number], b: [number, number]): [number, number] => [
    a[0] - b[0],
    a[1] - b[1],
  ];
  const cmul = (a: [number, number], b: [number, number]): [number, number] => [
    a[0] * b[0] - a[1] * b[1],
    a[0] * b[1] + a[1] * b[0],
  ];
  const cdiv = (a: [number, number], b: [number, number]): [number, number] => {
    const n = b[0] * b[0] + b[1] * b[1];
    return [(a[0] * b[0] + a[1] * b[1]) / n, (a[1] * b[0] - a[0] * b[1]) / n];
  };
  const evalp = (z: [number, number]): [number, number] => {
    let acc: [number, number] = c[d]!;
    for (let i = d - 1; i >= 0; i--) acc = cadd(cmul(acc, z), c[i]!);
    return acc;
  };
  // Aberth-style spread starting points on a circle
  const z: Array<[number, number]> = [];
  for (let i = 0; i < d; i++) {
    const th = (2 * Math.PI * i) / d + 0.35;
    z.push([0.9 * Math.cos(th), 0.9 * Math.sin(th)]);
  }
  for (let iter = 0; iter < 2000; iter++) {
    let moved = 0;
    for (let i = 0; i < d; i++) {
      let den: [number, number] = [1, 0];
      for (let j = 0; j < d; j++) {
        if (i === j) continue;
        den = cmul(den, csub(z[i]!, z[j]!));
      }
      if (den[0] === 0 && den[1] === 0) continue;
      const delta = cdiv(evalp(z[i]!), den);
      if (!isFinite(delta[0]) || !isFinite(delta[1])) continue;
      z[i] = csub(z[i]!, delta);
      moved = Math.max(moved, Math.abs(delta[0]) + Math.abs(delta[1]));
    }
    if (moved < 1e-15) break;
  }
  return z.map(([a, b]) => [a * R, b * R]);
}

/**
 * `complex_roots(p)` (complex_roots.py:154) for a squarefree integer
 * polynomial: certified isolating boxes for all `deg p` complex roots.
 *
 * The returned boxes are pairwise disjoint and each provably contains exactly
 * one root (interval Newton), so returning `deg p` of them is a proof that
 * every root has been found.
 */
export function complex_roots(p: ZPoly, min_prec = 53): CI[] {
  const d = zpDeg(p);
  if (d < 1) return [];
  const pd = zpDerivative(p);
  const rts = durandKerner(p);
  let prec = 53;
  for (let attempt = 0; attempt < 12; attempt++) {
    const workPrec = Math.max(prec, min_prec) + 16;
    const irts: CI[] = [];
    let ok = true;
    for (const [a, b] of rts) {
      const start = ci_exact(dyadicFromNumber(a, prec), dyadicFromNumber(b, prec));
      const irt = refine_root(p, pd, start, workPrec);
      if (irt === null) {
        ok = false;
        break;
      }
      irts.push(irt);
    }
    if (ok && irts.length === d && intervals_disjoint(irts)) return irts;
    prec *= 2;
  }
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: complex_roots could not certify the roots of this polynomial'
  );
}

/** Exact rational value of a JS double. */
export function dyadicFromNumber(x: number, prec: number): Rational {
  if (!isFinite(x)) throw new ValueError('non-finite root estimate');
  if (x === 0) return Rational.zero();
  // decompose the double exactly
  let e = 0;
  let m = x;
  while (Math.abs(m) < 1 << 20) {
    m *= 2 ** 32;
    e -= 32;
    if (e < -2000) break;
  }
  while (!Number.isInteger(m)) {
    m *= 2;
    e -= 1;
    if (e < -3000) break;
  }
  const r = new Rational(BigInt(Math.round(m)), 1n).mul(
    e >= 0 ? new Rational(1n << BigInt(e)) : new Rational(1n, 1n << BigInt(-e))
  );
  return roundDyadic(r, prec, 1);
}

/* ------------------------------------------------------------------ */
/* Presentation                                                        */
/* ------------------------------------------------------------------ */

/** Approximate an interval by a JS double (the midpoint). */
export function ri_toNumber(x: RI): number {
  return ri_center(x).toNumber();
}

/**
 * PARI's `cmp_complex_appr` (rootpol.c:1862), the order in which `polroots`
 * returns its roots and therefore the order of SageMath's
 * `K.embeddings(CC)`: increasing real part, then increasing imaginary part.
 */
export function cmp_complex_appr(a: CI, b: CI): number {
  const ar = ri_center(a.re);
  const br = ri_center(b.re);
  const c = ar.cmp(br);
  if (c !== 0) return c;
  return ri_center(a.im).cmp(ri_center(b.im));
}

/* ------------------------------------------------------------------ */
/* Archimedean places                                                  */
/* ------------------------------------------------------------------ */

/** The minimal view of a number field that an embedding needs. */
export interface EmbeddableField {
  degree(): number;
  toString(): string;
}

/** An element that can be written out in the power basis of the generator. */
export interface EmbeddableElement {
  list(): Rational[];
}

/**
 * A ring homomorphism from a number field into `RR` or `CC`, i.e. an
 * archimedean place, determined by the image of the generator.
 *
 * The image is stored as a *certified* complex interval, so `__call__` returns
 * an enclosure rather than a rounded float: the interval provably contains the
 * true value.
 *
 * @see Reference: sage/rings/number_field/number_field.py:9375 (embeddings)
 */
export class NumberFieldEmbedding {
  private readonly _domain: EmbeddableField;
  private readonly _image: CI;
  private readonly _prec: number;
  private readonly _real: boolean;

  constructor(domain: EmbeddableField, image: CI, prec: number, isReal: boolean) {
    this._domain = domain;
    this._image = image;
    this._prec = prec;
    this._real = isReal;
  }

  domain(): EmbeddableField {
    return this._domain;
  }

  /** `'RR'` for a real place, `'CC'` for a complex one. */
  codomain(): 'RR' | 'CC' {
    return this._real ? 'RR' : 'CC';
  }

  prec(): number {
    return this._prec;
  }

  is_real(): boolean {
    return this._real;
  }

  /** The image of the field generator, as a certified complex interval. */
  im_gens(): CI[] {
    return [this._image];
  }

  /** `sigma(x)`, as a certified complex interval. */
  __call__(x: EmbeddableElement): CI {
    const coeffs = x.list();
    const p = this._prec + 16;
    let acc = ci_exact(Rational.zero(), Rational.zero());
    for (let i = coeffs.length - 1; i >= 0; i--) {
      acc = ci_mul(acc, this._image, p);
      acc = ci_add(acc, ci_exact(coeffs[i]!, Rational.zero()), p);
    }
    return acc;
  }

  /** `sigma(x)` as a pair of JS doubles. */
  evalNumber(x: EmbeddableElement): { re: number; im: number } {
    const v = this.__call__(x);
    return { re: ri_toNumber(v.re), im: ri_toNumber(v.im) };
  }

  /** `|sigma(x)|^2` as a certified real interval. */
  absSquared(x: EmbeddableElement): RI {
    const v = this.__call__(x);
    const p = this._prec + 16;
    return ri_add(ri_mul(v.re, v.re, p), ri_mul(v.im, v.im, p), p);
  }

  toString(): string {
    const re = ri_toNumber(this._image.re);
    const im = ri_toNumber(this._image.im);
    const img = this._real
      ? `${re}`
      : `${re} ${im < 0 ? '-' : '+'} ${Math.abs(im)}*I`;
    return (
      `Ring morphism:\n  From: ${this._domain.toString()}\n` +
      `  To:   ${this._real ? 'Real' : 'Complex'} Field with ${this._prec} bits of precision\n` +
      `  Defn: a |--> ${img}`
    );
  }
}

export const _internal = {
  ri,
  ri_add,
  ri_sub,
  ri_mul,
  ri_center,
  ri_diameter,
  ri_overlaps,
  ci_mul,
  ci_div,
  ci_polyeval,
  ci_center,
  ci_diameter,
  ci_overlaps,
};
