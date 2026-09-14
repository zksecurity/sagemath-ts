/** Encode successive ceil-halvings of a positive precision, basemath/Zp.c. */
export function quadratic_prec_mask(n: number): bigint {
  if (!Number.isSafeInteger(n) || n < 2)
    throw new RangeError('quadratic_prec_mask requires an integer at least 2');
  let a = n,
    mask = 0n;
  for (let i = 1n; ; i++, mask <<= 1n) {
    mask |= BigInt(a % 2);
    a = Math.ceil(a / 2);
    if (a === 1) return mask | (1n << i);
  }
}

import {
  FpX_add,
  FpX_sub,
  FpX_mul,
  FpX_red,
  FpX_normalize,
  FpX_Fp_mul,
  FpXQ_mul,
} from './ffinit.js';
import {
  FpX_divrem,
  FpX_eval,
  FpX_FpXQ_eval,
  FpX_FpXQV_eval,
  FpXQ_powers,
  FpX_extgcd,
} from './FpX.js';
import { ZX_mul, ZX_deriv } from './ZX.js';
import { Fp_inv, Fp_sub, Fp_mul, Fp_sqr } from './ff.js';
import { FpX_roots, FpX_split_part } from './FpX_factor.js';
import { polynomialQuotientInverse } from './_polynomial_quotient_power.js';
import { polynomialQuotient } from './_polynomial_quotient.js';
import { residue } from './_polynomial_division.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { displayExtensionError } from './_extension_display.js';
import { PariError } from './errors.js';
import { brent_kung_optpow } from './RgX.js';
type ZX = bigint[];
const addZ = (a: ZX, b: ZX): ZX =>
  trimPolynomial(
    Array.from({ length: Math.max(a.length, b.length) }, (_, i) => (a[i] ?? 0n) + (b[i] ?? 0n))
  );
const subZ = (a: ZX, b: ZX): ZX =>
  addZ(
    a,
    b.map((c) => -c)
  );
const mulZ = (a: ZX, b: bigint): ZX => trimPolynomial(a.map((c) => c * b));
/** Native exact divisions are unchecked. Reject violations of the lifting
 * precondition instead of reproducing uninitialized/unchecked GMP results.
 * @see Deviation: PARI Hensel lifting adapters
 */
function divExact(a: ZX, d: bigint): ZX {
  return trimPolynomial(
    a.map((c) => {
      if (c % d !== 0n) throw new RangeError('Hensel lifting requires exact polynomial division');
      return c / d;
    })
  );
}
/** Zp.c:809: lift the root and reciprocal derivative along the native mask. */
export function ZpX_liftroot(f: ZX, a: bigint, p: bigint, e: number): bigint {
  let q = p;
  a = residue(a, q);
  if (e === 1) return a;
  let mask = quadratic_prec_mask(e),
    fr = FpX_red(f, q),
    W = Fp_inv(FpX_eval(ZX_deriv(fr), a, q), q);
  for (;;) {
    q *= q;
    if (mask & 1n) q /= p;
    mask >>= 1n;
    fr = FpX_red(f, q);
    a = Fp_sub(a, Fp_mul(W, FpX_eval(fr, a, q), q), q);
    if (mask === 1n) return a;
    W = Fp_sub(W << 1n, Fp_mul(Fp_sqr(W, q), FpX_eval(ZX_deriv(fr), a, q), q), q);
  }
}
interface LiftTree {
  v: ZX[];
  w: ZX[];
  link: number[];
}
/** Zp.c:383: merge the two smallest-degree factors, retaining stable ties. */
function buildTree(Q: ZX[], p: bigint): LiftTree {
  const k = Q.length,
    v: ZX[] = [[], ...Q.map((q) => trimPolynomial(q))],
    w: ZX[] = [[]],
    link = [0, ...Q.map((_, i) => -i - 1)];
  for (let j = 1, i = k + 1; j <= 2 * k - 5; j += 2, i++) {
    for (const start of [j, j + 1]) {
      let min = start;
      for (let s = start + 1; s < i; s++) if (v[s]!.length < v[min]!.length) min = s;
      [v[start], v[min]] = [v[min]!, v[start]!];
      [link[start], link[min]] = [link[min]!, link[start]!];
    }
    v[i] = FpX_mul(v[j]!, v[j + 1]!, p);
    link[i] = j;
  }
  for (let j = 1; j <= 2 * k - 3; j += 2) {
    const [d, u, vv] = FpX_extgcd(v[j]!, v[j + 1]!, p);
    if (d.length > 1)
      throw new PariError(
        'elements not coprime in BuildTree:\n    ' +
          displayExtensionError(v[j]!, 0) +
          '\n    ' +
          displayExtensionError(v[j + 1]!, 0)
      );
    if (d[0] === 1n) {
      w[j] = u;
      w[j + 1] = vv;
    } else {
      const inv = Fp_inv(d[0] ?? 0n, p);
      w[j] = FpX_Fp_mul(u, inv, p);
      w[j + 1] = FpX_Fp_mul(vv, inv, p);
    }
  }
  return { v, w, link };
}
/** Zp.c:443: lift a*b and a*u+b*v, with final inverse update optional. */
function henselLift(
  tree: LiftTree,
  j: number,
  f: ZX,
  pd: bigint,
  p0: bigint,
  noinv: boolean
): void {
  const { v: V, w: W } = tree,
    a = V[j]!,
    b = V[j + 1]!,
    u = W[j]!,
    v = W[j + 1]!;
  const correction = (g: ZX): [ZX, ZX] => {
    g = FpX_red(divExact(g, p0), pd);
    const z = FpX_mul(v, g, pd),
      [t, s] = FpX_divrem(z, a, pd);
    return [mulZ(s, p0), mulZ(FpX_red(addZ(ZX_mul(u, g), ZX_mul(t, b)), pd), p0)];
  };
  const [s, t] = correction(subZ(f, ZX_mul(a, b))),
    a2 = addZ(a, s),
    b2 = addZ(b, t);
  V[j] = a2;
  V[j + 1] = b2;
  if (noinv) return;
  const [vs, ut] = correction(subZ([1n], addZ(ZX_mul(u, a2), ZX_mul(v, b2))));
  W[j] = addZ(u, ut);
  W[j + 1] = addZ(v, vs);
}
function recTreeLift(
  tree: LiftTree,
  j: number,
  f: ZX,
  pd: bigint,
  p0: bigint,
  noinv: boolean
): void {
  if (j < 0) return;
  henselLift(tree, j, f, pd, p0, noinv);
  recTreeLift(tree, tree.link[j]!, tree.v[j]!, pd, p0, noinv);
  recTreeLift(tree, tree.link[j + 1]!, tree.v[j + 1]!, pd, p0, noinv);
}
/** Source-based factor tree and logarithmic precision schedule (Zp.c:555). */
function multiLift(f: ZX, Q: ZX[], p: bigint, e: number, inverses: boolean): LiftTree | ZX[] {
  if (Q.length < 2) throw new PariError('domain error in MultiLift: #(modular factors) < 2');
  if (e < 1) throw new PariError('domain error in MultiLift: precision < 1');
  if (e === 1) return Q.map((q) => trimPolynomial(q));
  const tree = buildTree(Q, p);
  let mask = quadratic_prec_mask(e),
    power = p;
  while (mask > 1n) {
    const pd = mask & 1n ? power / p : power;
    mask >>= 1n;
    recTreeLift(tree, tree.v.length - 2, f, pd, power, !inverses && mask === 1n);
    power *= pd;
  }
  if (inverses) return tree;
  const out: ZX[] = new Array(Q.length);
  for (let i = 1; i < tree.v.length; i++)
    if (tree.link[i]! < 0) out[-tree.link[i]! - 1] = tree.v[i]!;
  return out;
}
/** Zero-indexed factor vectors in the dependency; galconj adapts its facade. */
export function ZpX_liftfact(pol: ZX, Q: ZX[], p: bigint, e: number): ZX[] {
  pol = FpX_normalize(pol, p ** BigInt(e));
  if (Q.length === 1) return [pol];
  return multiLift(pol, Q, p, e, false) as ZX[];
}
function bezoutPropagate(tree: LiftTree, j: number, U: ZX | null, f: ZX, pe: bigint): void {
  if (j < 0) return;
  let Q = FpX_mul(tree.v[j]!, tree.w[j]!, pe),
    R: ZX;
  if (U) {
    Q = FpXQ_mul(Q, U, f, pe);
    R = FpX_sub(U, Q, pe);
  } else R = FpX_sub([1n], Q, pe);
  tree.w[j + 1] = Q;
  tree.w[j] = R;
  bezoutPropagate(tree, tree.link[j]!, R, f, pe);
  bezoutPropagate(tree, tree.link[j + 1]!, Q, f, pe);
}
export function bezout_lift_fact(pol: ZX, Q: ZX[], p: bigint, e: number): ZX[] {
  if (Q.length === 1) return [[1n]];
  const pe = p ** BigInt(e);
  pol = FpX_normalize(pol, pe);
  const tree = multiLift(pol, Q, p, e, true);
  if (Array.isArray(tree))
    throw new RangeError('bezout_lift_fact requires precision at least 2 for multiple factors');
  bezoutPropagate(tree, tree.v.length - 2, null, pol, pe);
  const out: ZX[] = new Array(Q.length);
  for (let i = 1; i < tree.v.length; i++)
    if (tree.link[i]! < 0) out[-tree.link[i]! - 1] = tree.w[i]!;
  return out;
}
function liftrootsFull(f: ZX, S: bigint[], p: bigint, e: number): bigint[] {
  const q = p ** BigInt(e),
    factors = S.map((a) => [-a, 1n]);
  return ZpX_liftfact(f, factors, p, e).map((f) => residue(-f[0]!, q));
}
export function ZpX_liftroots(f: ZX, S: bigint[], p: bigint, e: number): bigint[] {
  if (S.length === trimPolynomial(f).length - 1) return liftrootsFull(f, S, p, e);
  return S.map((a) => ZpX_liftroot(f, a, p, e));
}
export function ZpX_roots(F: ZX, p: bigint, e: number): bigint[] {
  const f = FpX_normalize(F, p),
    g = FpX_normalize(FpX_split_part(f, p), p);
  if (g.length < f.length) {
    const h = FpX_divrem(f, g, p)[0];
    F = ZpX_liftfact(F, [g, h], p, e)[0]!;
  }
  return liftrootsFull(F, FpX_roots(g, p), p, e);
}
/** Zp.c:1280: quotient Newton lift with shared evaluation powers and lifted
 * reciprocal derivative. Native n=1 returns S without coefficient reduction.
 */
export function ZpX_ZpXQ_liftroot(P: ZX, S: ZX, T: ZX, p: bigint, n: number): ZX {
  T = trimPolynomial(T);
  S = trimPolynomial(S);
  polynomialQuotient(T, p ** BigInt(n), false);
  if (n === 1) return S.slice();
  let mask = quadratic_prec_mask(n) >> 1n,
    q2 = p,
    q = p * p;
  let Tq = FpX_red(T, q),
    Tq2 = FpX_red(Tq, q2),
    Pq = FpX_red(P, q);
  let W = polynomialQuotientInverse(
    FpX_FpXQ_eval(FpX_red(ZX_deriv(P), q2), S, Tq2, q2),
    Tq2,
    q2,
    false
  );
  let Q = divExact(FpX_FpXQ_eval(Pq, S, Tq, q), q2);
  const r = brent_kung_optpow(trimPolynomial(P).length - 1, 4, 3);
  for (;;) {
    const H = FpXQ_mul(W, Q, Tq2, q2),
      Sq = FpX_sub(S, mulZ(H, q2), q);
    if (mask === 1n) return Sq;
    let qq = q * q;
    if (mask & 1n) qq /= p;
    mask >>= 1n;
    const Pqq = FpX_red(P, qq),
      Tqq = FpX_red(T, qq),
      Spow = FpXQ_powers(Sq, r, Tqq, qq);
    Q = divExact(FpX_FpXQV_eval(Pqq, Spow, Tqq, qq), q);
    const dP = FpX_FpXQV_eval(
      FpX_red(ZX_deriv(Pq), q),
      Spow.map((v) => FpX_red(v, q)),
      Tq,
      q
    );
    let Wq = divExact(FpX_sub(FpXQ_mul(W, dP, Tq, q), [1n], q), q2);
    Wq = FpX_sub(W, mulZ(FpXQ_mul(W, Wq, Tq2, q2), q2), q);
    S = Sq;
    W = Wq;
    q2 = q;
    q = qq;
    Tq2 = Tq;
    Tq = Tqq;
    Pq = Pqq;
  }
}
