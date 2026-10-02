/** PARI odd-extension elliptic point kernels.
 * Copyright (C) The PARI group; adapted under GPL-2.0-or-later.
 * Source: FlxqE.c:40–246.
 * Inputs are reduced ascending coefficient arrays over a valid odd finite field.
 * @see Deviation: PARI odd-extension elliptic kernels
 */
import {
  oddElliptic,
  type FqEllipticPoint,
  type FqEllipticChange,
  type FlxqECoefficient,
} from './_odd_elliptic.js';
export function FlxqE_add(
  P: FqEllipticPoint,
  Q: FqEllipticPoint,
  a: FlxqECoefficient,
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(1, T, p).add(P, Q, a);
}

export function FlxqE_dbl(
  P: FqEllipticPoint,
  a: FlxqECoefficient,
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(1, T, p).dbl(P, a);
}

export function FlxqE_neg(P: FqEllipticPoint, T: bigint[], p: bigint): FqEllipticPoint {
  return oddElliptic(1, T, p).neg(P);
}

export function FlxqE_sub(
  P: FqEllipticPoint,
  Q: FqEllipticPoint,
  a: FlxqECoefficient,
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(1, T, p).sub(P, Q, a);
}

export function FlxqE_mul(
  P: FqEllipticPoint,
  n: bigint,
  a: FlxqECoefficient,
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(1, T, p).mul(P, n, a);
}

export function FlxqE_changepoint(
  P: FqEllipticPoint,
  ch: FqEllipticChange,
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(1, T, p).change(P, ch);
}

export function FlxqE_changepointinv(
  P: FqEllipticPoint,
  ch: FqEllipticChange,
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(1, T, p).changeinv(P, ch);
}

import { gen_order, type GroupOrder } from './bb_group.js';
/** FlxqE.c:313: exact order from a supplied annihilating multiple.
 * @see Deviation: PARI generic and extension-curve order adapters
 */
export function FlxqE_order(P: FqEllipticPoint, order: GroupOrder, a: FlxqECoefficient, T: bigint[], p: bigint): bigint {
  return gen_order(P, order, (Q, n) => FlxqE_mul(Q, n, a, T, p), Q => Q.isInfinity);
}

import { extensionField } from './_extension_field.js';
import { Flxq_issquare, Flxq_sqrt, Flxq_trace } from './Flx.js';
import { gen_pow_i } from './bb_group.js';
/** Supersingular characteristic-three count; FlxqE.c:1348–1375.
 * q=3^n, n=degree(T), and a4 is nonzero. Corrects the native odd-degree
 * twist-sign omission for a randomly selected nonsquare sqrt(-a4).
 * @see Deviation: PARI characteristic-three supersingular count correction
 */
export function F3xq_ellcardj(a4: bigint[], a6: bigint[], T: bigint[], q: bigint, n: number): bigint {
  const F = extensionField(1, T, 3n), q1 = q + 1n;
  const minusA = F.neg(a4) as bigint[];
  if (!Flxq_issquare(minusA, T, 3n)) return q1;
  const root = Flxq_sqrt(minusA, T, 3n)!;
  const denominator = F.mul(minusA, root);
  const c = F.mul(a6, F.inv(denominator)) as bigint[];
  const t = Flxq_trace(c, T, 3n);
  if (n % 2 === 1) {
    if (t === 0n) return q1;
    const q3 = 3n ** BigInt((n + 1) / 2);
    let trace = (t === 1n) !== (n % 4 === 1) ? q3 : -q3;
    // x = root*X introduces the quadratic twist root^3. Its character is
    // chi(root); replacing root by -root must not change the curve's count.
    if (!Flxq_issquare(root, T, 3n)) trace = -trace;
    return q1 - trace;
  }
  const q2 = 3n ** BigInt(n / 2);
  const W = gen_pow_i(a4, q >> 2n, a => F.sqr(a) as bigint[],
    (a, b) => F.mul(a, b) as bigint[]);
  const sign = (W[0] === 1n) !== (n % 4 === 2);
  if (t !== 0n) return sign ? q1 + q2 : q1 - q2;
  return sign ? q1 - 2n * q2 : q1 + 2n * q2;
}

import { constantJCard } from './_extension_elliptic_cardinality.js';
/** Native constant-j extension count for p>3; q=p^n, n=degree(T).
 * @see Deviation: PARI constant-j extension counting adapters
 */
export function Flxq_ellcardj(a4: bigint[], a6: bigint[], j: bigint, T: bigint[], q: bigint, p: bigint, n: number): bigint {
  return constantJCard(1,a4,a6,j,T,q,p,n);
}

import { random_Flx, Flxq_norm } from './Flx.js';
import { FlxqX_nbroots } from './FpXQX_factor.js';
import { isqrt } from './ifactor.js';
import { kronecker, xgcd } from './ff.js';
import { PariError } from './errors.js';

/** FlxqE.c:1038–1104: multiple of ord(f) near h, congruent to h modulo B. */
function extensionOrderMultiple(
  f: FqEllipticPoint,
  h: bigint,
  bound: bigint,
  B: bigint,
  a4: bigint[],
  T: bigint[],
  p: bigint
): bigint {
  const E = oddElliptic(1, T, p);
  let ceiling = isqrt(bound / B);
  if (ceiling * ceiling * B < bound) ceiling++;
  const s = ceiling / 2n;
  const fh = E.mul(f, h, a4);
  if (fh.isInfinity) return h;
  const F = E.mul(f, B, a4);
  let P: FqEllipticPoint = fh;
  if (s < 3n) {
    let Q: FqEllipticPoint = P;
    for (let i = 1n; ; i++) {
      P = E.add(P, F, a4);
      if (P.isInfinity) return h + i * B;
      Q = E.sub(Q, F, a4);
      if (Q.isInfinity) return h - i * B;
    }
  }
  // Native hash buckets ultimately select the first equal x-coordinate in
  // stable baby-step order. Exact coefficient keys avoid hash collisions.
  const babies = new Map<string, bigint>();
  for (let i = 1n; i <= s; i++) {
    if (P.isInfinity) throw new PariError('bug in Flxq_ellcard baby steps');
    const key = P.x.join(',');
    if (!babies.has(key)) babies.set(key, i - 1n);
    P = E.add(P, F, a4);
    if (P.isInfinity) return h + i * B;
  }
  const fg = E.sub(P, fh, a4);
  if (fg.isInfinity) return s * B;
  P = fg;
  for (let i = 1n; ; i++) {
    if (P.isInfinity) throw new PariError('bug in Flxq_ellcard giant steps');
    const j = babies.get(P.x.join(','));
    if (j !== undefined) {
      const Q = E.add(E.mul(F, j, a4), fh, a4);
      if (!Q.isInfinity) {
        const sameY = P.y.length === Q.y.length && P.y.every((x, k) => x === Q.y[k]);
        return h + (s * (sameY ? -i : i) + j) * B;
      }
    }
    P = E.add(P, fg, a4);
  }
}

/** FlxqE.c:1191–1231: native extension Shanks–Mestre counter.
 * Requires a nonsingular short model in the native Shanks dispatch domain:
 * p>3 is a word prime, q=p^degree(T), j generates the full field, and the
 * preceding tiny-field/Satoh/Kedlaya branches do not apply; expi(q)<=62.
 * Special curves outside that domain can make both native and port searches stall.
 * @see Deviation: PARI extension Shanks counting adapter
 */
export function Flxq_ellcard_Shanks(
  a4: bigint[],
  a6: bigint[],
  q: bigint,
  T: bigint[],
  p: bigint
): bigint {
  const F = extensionField(1, T, p),
    n = T.length - 1,
    q1 = q + 1n,
    q2 = 2n * q1;
  const bound = isqrt(16n * q) + 1n;
  const roots = FlxqX_nbroots([a6, a4, [], [1n]], T, p);
  let A = roots === 0 ? 1n : 0n,
    B = roots === 3 ? 4n : 2n,
    KRO = -1;
  const mod = (a: bigint, b: bigint) => ((a % b) + b) % b;
  const closest = () => A + B * ((2n * (q1 - A) + B) / (2n * B));
  for (;;) {
    let h = closest();
    KRO = -KRO;
    let point: Extract<FqEllipticPoint, { isInfinity: false }>;
    for (;;) {
      const x = random_Flx(n, p),
        x2 = F.sqr(x);
      const u = F.add(a6, F.mul(F.add(a4, x2), x)) as bigint[];
      const symbol = u.length ? kronecker(Flxq_norm(u, T, p), p) : 0;
      if (symbol === KRO) {
        point = { isInfinity: false, x: F.mul(u, x) as bigint[], y: F.sqr(u) as bigint[] };
        break;
      }
    }
    const twistA = F.mul(a4, point.y) as bigint[];
    const multiple = extensionOrderMultiple(point, h, bound, B, twistA, T, p);
    h = FlxqE_order(point, multiple, twistA, T, p);
    const [g, u] = xgcd(B, h),
      m = h / g,
      nextB = B * m;
    A = mod(A + B * mod((-A / g) * u, m), nextB);
    B = nextB;
    if (B >= bound) {
      h = closest();
      return KRO === 1 ? h : q2 - h;
    }
    A = mod(q2 - A, B);
  }
}

import { FpX_rem } from './FpX.js';
import { FpXQ_mul } from './ffinit.js';
import { ZX_mul } from './ZX.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { gen_powu_i } from './bb_group.js';
import { Fp_sqrt } from './ff.js';
import { Zp_sqrtlift } from './Zp.js';
/** FlxqE.c:547–565: Frobenius on a prime cyclotomic modulus T=Phi_l.
 * x is reduced modulo T, gcd(p,l)=1; coefficients reduce modulo q.
 * @see Deviation: PARI scalar lifts and cyclotomic counting dependencies
 */
export function ZpXQ_frob_cyc(x: bigint[], T: bigint[], q: bigint, p: bigint): bigint[] {
  const length = T.length,
    out = Array<bigint>(length).fill(0n);
  for (let i = 0; i < x.length; i++) out[Number((BigInt(i) * p) % BigInt(length))] = x[i]!;
  return FpX_rem(out, T, q);
}
/** FlxqE.c:568–581: Frobenius using precomputed powers of X^degree(T).
 * Xm=[] selects the cyclotomic shortcut; otherwise Xm has p powers.
 * T is the lifted field modulus and x is reduced modulo T.
 * @see Deviation: PARI scalar lifts and cyclotomic counting dependencies
 */
export function ZpXQ_frob(
  x: bigint[],
  Xm: bigint[][],
  T: bigint[],
  q: bigint,
  p: bigint
): bigint[] {
  if (!Xm.length) return ZpXQ_frob_cyc(x, T, q, p);
  const degree = T.length - 1;
  const blocks = Array.from({ length: Xm.length }, () => Array<bigint>(degree).fill(0n));
  for (let i = 0; i < x.length; i++) {
    const exponent = BigInt(i) * p;
    blocks[Number(exponent / BigInt(degree))]![Number(exponent % BigInt(degree))] = x[i]!;
  }
  let sum: bigint[] = [];
  for (let i = 0; i < Xm.length; i++) {
    const product = ZX_mul(blocks[i]!, Xm[i]!);
    sum = trimPolynomial(
      Array.from(
        { length: Math.max(sum.length, product.length) },
        (_, j) => (sum[j] ?? 0n) + (product[j] ?? 0n)
      )
    );
  }
  return FpX_rem(sum, T, q);
}
/** FlxqE.c:699–723: norm by semidirect powering on a prime cyclotomic field.
 * T=Phi_l irreducible modulo the small word prime p, q=p^e, x a unit.
 * Native degree one returns a polynomial copy rather than a scalar.
 * @see Deviation: PARI scalar lifts and cyclotomic counting dependencies
 */
export function ZpXQ_norm_pcyc(x: bigint[], T: bigint[], q: bigint, p: bigint): bigint | bigint[] {
  const degree = T.length - 1;
  if (degree === 1) return x.slice();
  type State = [bigint[], bigint];
  const multiply = (a: State, b: State): State => [
    FpXQ_mul(a[0], ZpXQ_frob_cyc(b[0], T, q, a[1]), T, q),
    (a[1] * b[1]) % BigInt(degree + 1),
  ];
  const z = gen_powu_i<State>([x, p], BigInt(degree), (a) => multiply(a, a), multiply);
  return z[0][0] ?? 0n;
}
/** FlxqE.c:727–731: chosen scalar square root of the cyclotomic norm.
 * T has degree>=2, x is a unit and its norm is square modulo odd p.
 * @see Deviation: PARI scalar lifts and cyclotomic counting dependencies
 */
export function ZpXQ_sqrtnorm_pcyc(
  x: bigint[],
  T: bigint[],
  q: bigint,
  p: bigint,
  e: number
): bigint {
  const z = ZpXQ_norm_pcyc(x, T, q, p) as bigint;
  const root = Fp_sqrt(z, p);
  if (root === null) throw new PariError('square cyclotomic norm required');
  return Zp_sqrtlift(z, root, p, e);
}

import { ZpXQ_log } from './Zp.js';
import { FpXQ_trace } from './FpX.js';
import { Fp_div } from './ff.js';
import { cvtop } from './gen2.js';
import { Qp_exp } from './trans1.js';

/** FlxqE.c:735: square root of the norm of a=1 mod odd p, q=p^e, e>=2.
 * T is a monic unramified quotient modulus. Preserve native relative
 * precision e-1 when exponentiating half the logarithm's trace.
 * @see Deviation: PARI p-adic precision records
 */
export function ZpXQ_sqrtnorm(a: bigint[], T: bigint[], q: bigint, p: bigint, e: number): bigint {
  const s = Fp_div(FpXQ_trace(ZpXQ_log(a, T, p, e), T, q), 2n, q);
  return Qp_exp(cvtop(s, p, e - 1)).unit % q;
}

import { gen_ZpX_Newton, gen_ZpX_Dixon } from './Zp.js';
import { FpX_red, FpX_sub, FpXQ_pow } from './ffinit.js';
import { Flxq_lroot_fast_pre } from './Flx.js';

/** FlxqE.c:746–809: lift an element satisfying Frob(x)=x^p using native
 * Newton iteration and a Dixon solve of the linearized Frobenius equation.
 * T is Teichmuller-lifted, Tp=T mod p; Xm and sqx are the forward/inverse
 * Frobenius power tables. Small word p; N>=1. At N=1 copy x unchanged.
 * @see Deviation: PARI Teichmuller element lifting adapters
 */
export function Teichmuller_lift(
  x: bigint[],
  Xm: bigint[][],
  T: bigint[],
  sqx: bigint[][],
  Tp: bigint[],
  p: bigint,
  pi: bigint,
  N: number
): bigint[] {
  return gen_ZpX_Newton<[bigint[], bigint[]]>(
    x,
    p,
    N,
    (x2, q) => {
      const TN = FpX_red(T, q),
        XN = Xm.map((x) => FpX_red(x, q));
      const y2 = ZpXQ_frob(x2, XN, TN, q, p);
      const x1 = FpXQ_pow(x2, p - 1n, TN, q);
      return [FpX_sub(y2, FpXQ_mul(x2, x1, TN, q), q), x1];
    },
    (V, v, qM, M) =>
      gen_ZpX_Dixon<[bigint[], bigint[], bigint[][]]>(
        [FpX_red(v[1], qM), FpX_red(T, qM), Xm.map((x) => FpX_red(x, qM))],
        V,
        qM,
        p,
        M,
        (F, x2, q) => {
          const y2 = ZpXQ_frob(x2, F[2], F[1], q, p);
          const term = ZX_mul(F[0], x2).map((c) => c * p);
          return FpX_rem(FpX_sub(y2, term, q), F[1], q);
        },
        (x) => Flxq_lroot_fast_pre(FpX_red(x, p), sqx, Tp, p, pi)
      )
  );
}
