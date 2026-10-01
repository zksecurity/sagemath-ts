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
  FpXQ_pow,
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
import { residue, inverseCoefficient } from './_polynomial_division.js';
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

/** Polynomial or nested polynomial vector, corresponding to native FpXT_red. */
export type ZpPolynomialTree = bigint[] | ZpPolynomialTree[];
function reduceTree<F extends ZpPolynomialTree>(F: F, q: bigint): F {
  return (
    F.length === 0 || typeof F[0] === 'bigint'
      ? FpX_red(F as ZX, q)
      : (F as ZpPolynomialTree[]).map((f) => reduceTree(f, q))
  ) as F;
}
/** Zp.c:993: balanced Dixon solve, retaining native reduction/callback order.
 * q=p^N, N>=1; invl solves the linear equation modulo p.
 * @see Deviation: PARI polynomial Newton and Dixon adapters
 */
export function gen_ZpX_Dixon<F extends ZpPolynomialTree>(
  F: F,
  V: ZX,
  q: bigint,
  p: bigint,
  N: number,
  lin: (F: F, d: ZX, q: bigint) => ZX,
  invl: (d: ZX) => ZX
): ZX {
  if (!Number.isSafeInteger(N) || N < 1)
    throw new RangeError('Dixon precision must be a positive integer');
  V = FpX_red(V, q);
  if (N === 1) return invl(V);
  const N2 = Math.ceil(N / 2),
    M = N - N2;
  F = reduceTree(F, q);
  const qM = p ** BigInt(M),
    q2 = M === N2 ? qM : qM * p;
  const VN2 = gen_ZpX_Dixon(F, V, q2, p, N2, lin, invl);
  const bil = lin(F, VN2, q);
  const V2 = divExact(subZ(V, bil), q2);
  const VM = gen_ZpX_Dixon(F, V2, qM, p, M, lin, invl);
  return FpX_red(addZ(VN2, mulZ(VM, q2)), q);
}
/** Zp.c:1098: Newton lifting with the native ceil-halving precision mask.
 * x is a solution modulo p; evaluate returns [residual, ...callback state].
 * @see Deviation: PARI polynomial Newton and Dixon adapters
 */
export function gen_ZpX_Newton<V extends [ZX, ...unknown[]]>(
  x: ZX,
  p: bigint,
  n: number,
  evaluate: (x: ZX, q: bigint) => V,
  invd: (V: ZX, v: V, q: bigint, M: number) => ZX
): ZX {
  if (n === 1) return x.slice();
  let mask = quadratic_prec_mask(n),
    N = 1,
    q = p;
  while (mask > 1n) {
    const N2 = N,
      q2 = q;
    N *= 2;
    let M: number, qM: bigint;
    if (mask & 1n) {
      N--;
      M = N2 - 1;
      qM = q2 / p;
      q = qM * q2;
    } else {
      M = N2;
      qM = q2;
      q = q2 * q2;
    }
    mask >>= 1n;
    const v = evaluate(x, q),
      V = divExact(v[0], q2);
    x = FpX_sub(x, mulZ(invd(V, v, qM, M), q2), q);
  }
  return x;
}
/** Zp.c:1161: x is an inverse of a modulo (T,p).
 * @see Deviation: PARI polynomial Newton and Dixon adapters
 */
export function ZpXQ_invlift(a: ZX, x: ZX, T: ZX, p: bigint, e: number): ZX {
  return gen_ZpX_Newton<[ZX, ZX]>(
    x,
    p,
    e,
    (x, q) => [FpX_sub(FpXQ_mul(x, FpX_red(a, q), FpX_red(T, q), q), [1n], q), x],
    (V, v, q) => FpXQ_mul(V, v[1], FpX_red(T, q), q)
  );
}
/** Zp.c:1169: word inverse initialization, then native Newton lifting.
 * @see Deviation: PARI polynomial Newton and Dixon adapters
 */
export function ZpXQ_inv(a: ZX, T: ZX, p: bigint, e: number): ZX {
  const ai = polynomialQuotientInverse(FpX_red(a, p), FpX_red(T, p), p, p > 0n && p < 1n << 64n);
  return ZpXQ_invlift(a, ai, T, p, e);
}
/** Zp.c:1183: multiply by the lifted inverse; q=p^e.
 * @see Deviation: PARI polynomial Newton and Dixon adapters
 */
export function ZpXQ_div(a: ZX, b: ZX, T: ZX, q: bigint, p: bigint, e: number): ZX {
  return FpXQ_mul(a, ZpXQ_inv(b, T, p, e), T, q);
}

import { ZX_sqr } from './ZX.js';
import { FpXQX_mul, FpXQX_red } from './FpXX.js';
import { type ExtensionPolynomial } from './_extension_polynomial.js';
import { gen_powu_i } from './bb_group.js';
function splitPolynomial(f: ZX, k: number): ZX[] {
  const parts = Array.from({ length: k }, () => [] as ZX);
  for (let i = 0; i < f.length; i++) parts[i % k]!.push(f[i]!);
  return parts.map(trimPolynomial);
}
const shiftPolynomial = (f: ZX, n: number): ZX =>
  f.length ? [...Array<bigint>(n).fill(0n), ...f] : [];
function canonicalLinear(F: ZX[], V: ZX, q: bigint): ZX {
  const parts = splitPolynomial(V, F.length);
  const dot = parts.reduce((sum, part, i) => addZ(sum, ZX_mul(part, F[i]!)), [] as ZX);
  return FpX_sub(V, dot, q);
}
/** Zp.c:1412–1455: specialized cubic Frobenius lift. */
function ternaryTeichmuller(P: ZX, n: number): ZX {
  return gen_ZpX_Newton<[ZX, ZX, ZX, ZX, ZX, ZX, ZX]>(
    P,
    3n,
    n,
    (f, q) => {
      const [h1, h2, h3] = splitPolynomial(f, 3) as [ZX, ZX, ZX];
      const h1s = ZX_sqr(h1),
        h2s = ZX_sqr(h2),
        h3s = ZX_sqr(h3);
      const h12 = ZX_mul(h1, h2),
        h13 = ZX_mul(h1, h3),
        h23 = ZX_mul(h2, h3);
      const th = ZX_mul(subZ(h2s, mulZ(h13, 3n)), h2);
      const value = addZ(
        shiftPolynomial(ZX_mul(h3, h3s), 2),
        addZ(shiftPolynomial(th, 1), ZX_mul(h1, h1s))
      );
      return [FpX_sub(f, value, q), h1s, h2s, h3s, h12, h13, h23];
    },
    (V, v, q, M) => {
      const [, h1s, h2s, h3s, h12, h13, h23] = v;
      const F = [
        subZ(h1s, shiftPolynomial(h23, 1)),
        shiftPolynomial(subZ(h2s, h13), 1),
        subZ(shiftPolynomial(h3s, 2), shiftPolynomial(h12, 1)),
      ].map((f) => mulZ(f, 3n));
      return gen_ZpX_Dixon(F, V, q, 3n, M, canonicalLinear, (d) => d);
    }
  );
}
/** Zp.c:1457–1543: canonical polynomial lift using the cyclic norm product.
 * P has reduced coefficients modulo a word prime p; n>=1.
 * @see Deviation: PARI polynomial Newton and Dixon adapters
 */
export function Flx_Teichmuller(P: ZX, p: bigint, n: number): ZX {
  if (p === 3n) return ternaryTeichmuller(P, n);
  const prime = Number(p);
  if (!Number.isSafeInteger(prime) || prime < 2)
    throw new RangeError('Teichmuller characteristic must fit a polynomial array');
  const shift = (P: ExtensionPolynomial, n: bigint): ExtensionPolynomial =>
    P.map((c, i) => {
      const s = Number((n * BigInt(i)) % p),
        r = Array<bigint>(prime).fill(0n);
      if (typeof c === 'bigint') r[s] = c;
      else for (let j = 0; j < c.length; j++) r[(j + s) % prime] = c[j]!;
      return trimPolynomial(r);
    });
  return gen_ZpX_Newton<[ZX, ZX[]]>(
    P,
    p,
    n,
    (f, q) => {
      const T = [-1n, ...Array<bigint>(prime - 1).fill(0n), 1n];
      type State = [ExtensionPolynomial, bigint];
      const multiply = (a: State, b: State): State => [
        FpXQX_mul(a[0], shift(b[0], a[1]), T, q),
        a[1] + b[1],
      ];
      const product = gen_powu_i<State>(
        [shift(f, 1n), 1n],
        p - 1n,
        (a) => multiply(a, a),
        multiply
      )[0];
      const norm = trimPolynomial(
        FpXQX_red(product, Array<bigint>(prime).fill(1n), q).map((c) =>
          typeof c === 'bigint' ? c : (c[0] ?? 0n)
        )
      );
      const value = FpX_mul(norm, f, q).filter((_, i) => i % prime === 0);
      const d = splitPolynomial(norm, prime);
      const F = [
        mulZ(d[0]!, p),
        ...Array.from({ length: prime - 1 }, (_, i) =>
          mulZ(shiftPolynomial(d[prime - 1 - i]!, 1), p)
        ),
      ];
      return [subZ(f, value), F];
    },
    (V, v, q, M) => gen_ZpX_Dixon(v[1], V, q, p, M, canonicalLinear, (d) => d)
  );
}

import { Fp_pow } from './ff.js';
import { gen_pow_i } from './bb_group.js';
/** Zp.c:133–200: lift a supplied root of X^n-b, retaining reciprocal-derivative
 * updates, the native precision mask and signed binary remainders.
 * e>=1; for e>1 the derivative must be a unit modulo p and n positive.
 * @see Deviation: PARI scalar lifts and cyclotomic counting dependencies
 */
export function Zp_sqrtnlift(b: bigint, n: bigint, a: bigint, p: bigint, e: number): bigint {
  if (e === 1) return a;
  const square = n === 2n,
    binary = p === 2n;
  let mask = quadratic_prec_mask(e),
    q = p,
    precision = 1;
  let w = Fp_inv(square ? 2n * a : Fp_mul(n, Fp_pow(a, n - 1n, p), p), p);
  // Native Fp_pow2n passes n=0 into gen_pow during the derivative update
  // for a binary linear root (n=1), causing undefined behavior. Supply the
  // multiplicative identity explicitly rather than reproducing that defect.
  const pow2 = (x: bigint, n: bigint, q: bigint): bigint =>
    n === 0n
      ? 1n
      : gen_pow_i(
          x,
          n,
          (x) => (x * x) % q,
          (x, y) => (x * y) % q
        );
  for (;;) {
    if (binary) {
      precision = 2 * precision - Number(mask & 1n);
      mask >>= 1n;
      const modulus = 1n << BigInt(precision);
      a = (a - w * (pow2(a, n, modulus) - b)) % modulus;
      if (mask === 1n) break;
      w = 2n * w - ((((w * w) % modulus) * (n * pow2(a, n - 1n, modulus))) % modulus);
      continue;
    }
    q *= q;
    if (mask & 1n) q /= p;
    mask >>= 1n;
    a = residue(a - w * (Fp_pow(a, n, q) - b), q);
    if (mask === 1n) break;
    const correction = Fp_mul(Fp_sqr(w, q), square ? a : n * Fp_pow(a, n - 1n, q), q);
    if (q < 1n << 64n && n > 0n && n < 1n << 64n)
      w = residue(square ? 2n * (w - correction) : 2n * w - correction, q);
    else w = square ? 2n * (w - correction) : 2n * w - correction;
  }
  return binary && a < 0n ? a + (1n << BigInt(precision)) : a;
}
/** Zp.c:204: square-root specialization of Zp_sqrtnlift.
 * @see Deviation: PARI scalar lifts and cyclotomic counting dependencies
 */
export function Zp_sqrtlift(b: bigint, a: bigint, p: bigint, e: number): bigint {
  return Zp_sqrtnlift(b, 2n, a, p, e);
}

/** Zp.c:55: lift an inverse, multiplying by b only at the final stage. */
function Zp_divlift(b: bigint | null, a: bigint, x: bigint, p: bigint, n: number): bigint {
  if (n === 1) return x;
  let mask = quadratic_prec_mask(n),
    q = p;
  while (mask > 1n) {
    const q2 = q;
    q *= q;
    if (mask & 1n) q /= p;
    mask >>= 1n;
    const v = residue(x * residue(a, q) - 1n, q);
    if (mask > 1n || b === null) x = residue(x - v * x, q);
    else {
      const y = residue(x * b, q);
      x = residue(y - v * residue(y, q2), q);
    }
  }
  return x;
}
/** Lift a supplied reciprocal modulo p to p^e; e=1 returns x unchanged. */
export function Zp_invlift(a: bigint, x: bigint, p: bigint, e: number): bigint {
  return Zp_divlift(null, a, x, p, e);
}
/** Zp.c:94: native word/generic inverse initialization and Newton lifting. */
export function Zp_inv(a: bigint, p: bigint, e: number): bigint {
  return Zp_invlift(a, inverseCoefficient(residue(a, p), p, p < 1n << 64n), p, e);
}
/** Zp.c:108: native quotient lift. At e=1 PARI returns 1/a, ignoring b. */
export function Zp_div(b: bigint, a: bigint, p: bigint, e: number): bigint {
  return Zp_divlift(b, a, inverseCoefficient(residue(a, p), p, p < 1n << 64n), p, e);
}

/** Zp.c:300: p-adic exponential by digit splitting and binary-split Taylor
 * sums. Requires p prime, e>=1, and p|a (4|a for p=2).
 * @see Deviation: PARI p-adic series precision adapters
 */
export function Zp_exp(a: bigint, p: bigint, e: number): bigint {
  const binary = p === 2n,
    word = p < 1n << 64n,
    pe = p ** BigInt(e);
  let N = binary ? e : e + Number(BigInt(e) / (p - 2n));
  let trunc = binary ? 4 : 2,
    truncMod = p * p,
    ans = 1n,
    denominator = 1n;
  for (;;) {
    const f = binary ? a % (1n << BigInt(trunc)) : residue(a, truncMod);
    a -= f;
    if (f !== 0n) {
      const num = Array<bigint>(N + 1).fill(1n);
      const den = Array.from({ length: N + 1 }, (_, i) => BigInt(i || 1));
      let step = 1,
        hpow = f;
      for (;;) {
        for (let i = 0; i <= N - step; i += step * 2) {
          num[i] = num[i]! * den[i + step]! + hpow * num[i + step]!;
          den[i] = den[i]! * den[i + step]!;
        }
        step *= 2;
        if (step > N) break;
        hpow *= hpow;
      }
      if (word) {
        let valuation = 0n;
        for (let n = BigInt(N) / p; n; n /= p) valuation += n;
        const d = p ** valuation;
        num[0] = num[0]! / d;
        den[0] = den[0]! / d;
      }
      ans = residue(ans * num[0]!, pe);
      denominator = residue(denominator * den[0]!, pe);
    }
    if (trunc > e) break;
    if (!binary) truncMod *= truncMod;
    trunc *= 2;
    N = Math.floor(N / 2);
  }
  return Zp_div(ans, denominator, p, e);
}

/** Zp.c:902: normalized atanh argument after raising a to p^k. */
function ZpXQ_log_to_ath(x: ZX, k: number, T: ZX, p: bigint, e: number, pe: bigint): ZX {
  const bd = addZ(x, [1n]);
  const binary = p === 2n;
  // Native shifti truncates toward zero, including negative coefficients.
  const bn = binary
    ? trimPolynomial(x.map((c) => c / (1n << BigInt(k + 1))))
    : divExact(subZ(x, [1n]), p ** BigInt(k));
  const denominator = binary ? trimPolynomial(bd.map((c) => c / 2n)) : bd;
  const bdi = ZpXQ_invlift(denominator, binary ? [1n] : [Fp_inv(2n, p)], T, p, e);
  return FpXQ_mul(bn, bdi, T, pe);
}
/** Zp.c:925: log(a) modulo (T,p^N), a=1 mod p. For p=2 require N>=2.
 * Retains native power/atanh algorithm with an integer sizing heuristic.
 * @see Deviation: PARI p-adic series precision adapters
 */
export function ZpXQ_log(a: ZX, T: ZX, p: bigint, N: number): ZX {
  const binary = p === 2n,
    word = !binary && p < 1n << 64n;
  // Replace floating log2(p) with floor(log2(p)); within a constant factor,
  // with the same cube-root balance and exact modular arithmetic throughout.
  const lp = BigInt(p.toString(2).length - 1),
    half = BigInt(Math.floor(N / 2));
  let k = 1;
  while (BigInt(k + 1) ** 3n * lp * lp <= half) k++;
  const e = binary ? N - 1 : N,
    pe = p ** BigInt(e),
    pNk = p ** BigInt(N + k);
  const l = Math.trunc((e - 2) / (2 * (k + Number(binary))));
  const ak = FpXQ_pow(a, p ** BigInt(k), FpX_red(T, pNk), pNk);
  const b = ZpXQ_log_to_ath(ak, k, T, p, e, pe);
  const pol: ZX = [];
  for (let i = 0; i <= l; i++) {
    let z = BigInt(2 * i + 1),
      w = 0;
    if (word)
      while (z % p === 0n) {
        z /= p;
        w++;
      }
    pol.push(residue(p ** BigInt(2 * i * k - w) * Fp_inv(z, pe), pe));
  }
  const s = FpX_FpXQ_eval(pol, FpXQ_mul(b, b, T, pe), T, pe);
  const result = mulZ(FpXQ_mul(b, s, T, pe), 2n);
  return binary ? result : FpX_red(result, pe);
}
