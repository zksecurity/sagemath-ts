/** PARI QX_factor.c:1020–1212: deflation and modular integer-polynomial GCD. */
import { gcd } from './ff.js';
import { nextprime } from './ifactor.js';
import { Flx_gcd } from './Flx.js';
import { inverseCoefficient, residue } from './_polynomial_division.js';
import { trimPolynomial } from './_polynomial_packing.js';

function content(a: bigint[]): bigint {
  return a.reduce((c, x) => gcd(c, x), 0n);
}
function divides(a: bigint[], b: bigint[]): boolean {
  const r = a.slice(),
    d = b.length - 1,
    lead = b[d]!;
  for (let i = r.length - 1; i >= d; i--) {
    if (r[i]! % lead !== 0n) return false;
    const q = r[i]! / lead;
    for (let j = 0; j <= d; j++) r[i - d + j] -= q * b[j]!;
  }
  return r.every((c) => c === 0n);
}
/** Content/valuation extraction, modular GCD, centered CRT and exact division.
 * @see Deviation: PARI integral-basis denominator adapters
 */
export function ZX_gcd(x: bigint[], y: bigint[]): bigint[] {
  let a = trimPolynomial(x),
    b = trimPolynomial(y);
  if (!a.length) return b;
  if (!b.length) return a;
  const ca = content(a),
    cb = content(b),
    c = gcd(ca, cb);
  a = a.map((v) => v / ca);
  b = b.map((v) => v / cb);
  let va = 0,
    vb = 0;
  while (a[va] === 0n) va++;
  while (b[vb] === 0n) vb++;
  a = a.slice(va);
  b = b.slice(vb);
  const v = Math.min(va, vb),
    prefix = new Array<bigint>(v).fill(0n);
  if (a.length === 1 || b.length === 1) return [...prefix, c];
  const g = gcd(a.at(-1)!, b.at(-1)!),
    ag = a.map((v) => v * g),
    bg = b.map((v) => v * g);
  let H: bigint[] = [],
    modulus = 1n,
    degree = Infinity,
    p = 1n << 63n;
  for (;;) {
    p = nextprime(p + 1n);
    if (residue(a.at(-1)!, p) === 0n || residue(b.at(-1)!, p) === 0n) continue;
    const h = Flx_gcd(
      a.map((v) => residue(v, p)),
      b.map((v) => residue(v, p)),
      p
    );
    if (h.length === 1) return [...prefix, c];
    if (h.length > degree) continue;
    const scale = residue(g * inverseCoefficient(h.at(-1)!, p, true), p);
    const hp = h.map((v) => residue(v * scale, p));
    if (h.length < degree) {
      degree = h.length;
      H = hp;
      modulus = p;
    } else {
      const inv = inverseCoefficient(residue(modulus, p), p, true);
      H = H.map((v, i) => v + modulus * residue((hp[i]! - v) * inv, p));
      modulus *= p;
    }
    H = H.map((v) => {
      const r = residue(v, modulus);
      return 2n * r > modulus ? r - modulus : r;
    });
    if (divides(bg, H) && divides(ag, H)) {
      const hc = content(H);
      return [...prefix, ...H.map((v) => (v / hc) * c)];
    }
  }
}
export function ZX_is_squarefree(input: bigint[]): boolean {
  let x = trimPolynomial(input);
  if (!x.length) return false;
  let m = 0n;
  for (let i = 1; i < x.length; i++)
    if (x[i]) {
      m = gcd(m, BigInt(i));
      if (m === 1n) break;
    }
  if (m > 1n) {
    if (x[0] === 0n) return false;
    x = x.filter((_, i) => i % Number(m) === 0);
  }
  const derivative = x.slice(1).map((v, i) => v * BigInt(i + 1));
  return ZX_gcd(x, derivative).length === 1;
}

import { ZM_lll_norms, LLL_INPLACE, LLL_NOFLATTER } from './lll.js';
import { PariError } from './errors.js';
import { ZX_Z_eval } from './ZX.js';
import { cmpir } from './kernel/none/level1.js';
import { vecbinomial as binomial } from './bibli2.js';
import { powruhalf } from './trans1.js';
import { ceil_safe } from './gen3.js';
import {
  itor,
  sqrr,
  sqrtr_abs,
  addrr,
  addir,
  mulir,
  divri,
  divrr,
  mulrr,
  real_0,
  type MpReal,
} from './qfb.js';
import { mppi } from './trans1.js';
import { cmprr } from './kernel/none/cmp.js';

function compare(x: bigint | MpReal, y: bigint | MpReal): number {
  if (typeof x !== 'bigint' && typeof y !== 'bigint') return cmprr(x, y);
  if (typeof x === 'bigint' && typeof y === 'bigint') return x < y ? -1 : x > y ? 1 : 0;
  return typeof x === 'bigint' ? cmpir(x, y as MpReal) : -cmpir(y as bigint, x);
}
/** Native QX_factor.c Mignotte_bound; the polynomial has positive degree and is normalized.
 * @see Deviation: PARI factorization scalar and bound adapters
 */
export function Mignotte_bound(S: bigint[]): bigint | MpReal {
  const d = S.length - 1,
    lc = S[d]!,
    bin = binomial(d - 1),
    binlc = lc === 1n || lc === -1n ? bin : bin.map((v) => v * lc);
  let s = sqrr(itor(S[0]!, 64));
  for (let i = 1; i <= d; i++) s = addrr(s, sqrr(itor(S[i]!, 64)));
  const N2 = sqrtr_abs(s);
  let C: bigint | MpReal = binlc[0]!;
  if (compare(C, N2) < 0) C = N2;
  for (let i = 1; i < d; i++) {
    const t = addir(binlc[i]!, mulir(bin[i - 1]!, N2));
    if (compare(C, t) < 0) C = t;
  }
  return C;
}
/** Native QX_factor.c Beauzamy_bound; the polynomial has positive degree and is normalized.
 * @see Deviation: PARI factorization scalar and bound adapters
 */
export function Beauzamy_bound(S: bigint[]): MpReal {
  const d = S.length - 1,
    bin = binomial(d);
  let s = real_0(64);
  for (let i = 0; i <= d; i++) {
    const c = S[i]!;
    if (c) s = addrr(s, divri(itor(c * c, 64), bin[i]!));
  }
  let C = powruhalf(itor(3n, 64), BigInt(3 + 2 * d)) as MpReal;
  C = divrr(mulrr(C, s), mulir(BigInt(4 * d), mppi(64)));
  const lc = S[d]!;
  return mulir(lc < 0n ? -lc : lc, sqrtr_abs(C));
}
/** Native QX_factor.c factor_bound; the polynomial has positive degree and is normalized.
 * @see Deviation: PARI factorization scalar and bound adapters
 */
export function factor_bound(S: bigint[]): bigint {
  const a = Mignotte_bound(S),
    b = Beauzamy_bound(S);
  return ceil_safe(compare(a, b) <= 0 ? a : b);
}
/** Positive-degree integer-polynomial root bound, following Fujiwara and exact refinement.
 * @see Deviation: PARI factorization scalar and bound adapters
 */
export function root_bound(S: bigint[]): bigint {
  const d = S.length - 1,
    lc = S[d]! < 0n ? -S[d]! : S[d]!,
    Q = S.slice(0, -1).map((c) => (c < 0n ? -c : c));
  const log = (c: bigint) => {
    const r = itor(c, 64);
    return Math.log2(Number(r.m)) + (r.e - 63);
  };
  const loglc = log(lc);
  let L = Q[0] ? (log(Q[0]!) - loglc - 1) / d : -Infinity;
  for (let i = 1; i < d; i++) if (Q[i]) L = Math.max(L, (log(Q[i]!) - loglc) / (d - i));
  let k = Math.trunc(L + 1);
  const evalAt = (x: bigint) => ZX_Z_eval(Q, x);
  const shiftEval = (k: number) => {
    let z = 0n;
    for (let i = Q.length - 1; i >= 0; i--) z = (z << BigInt(k)) + Q[i]!;
    return z;
  };
  for (; k >= 0; k--) if (shiftEval(k) >= lc << BigInt(d * k)) break;
  if (k < 0) k = 0;
  let x = 1n << BigInt(k),
    y = x << 1n;
  if (d > 2000) return y;
  for (let j = 0; ; j++) {
    const z = (x + y) >> 1n;
    if (z === x || j > 5) break;
    if (evalAt(z) < lc * z ** BigInt(d)) y = z;
    else x = z;
  }
  return y;
}

/** QX_factor.c LLL_check_progress; the native debug-timing output is omitted.
 * @see Deviation: PARI factor recombination progress adapters
 */
export function LLL_check_progress(
  Bnorm: MpReal,
  n0: number,
  m: bigint[][],
  final: boolean
): bigint[][] | null {
  const [u, norm] = ZM_lll_norms(m, final ? 0.999 : 0.75, LLL_INPLACE | LLL_NOFLATTER);
  let R = m.length;
  if (R && !norm) throw new PariError('bug in PARI/GP (Segmentation Fault), please report');
  for (; R > 0; R--) if (cmprr(norm![R - 1]!, Bnorm) < 0) break;
  if (R > u.length) throw new PariError('bug in PARI/GP (Segmentation Fault), please report');
  if (!R) throw new PariError('bug in LLL_cmbf [no factor], please report');
  if (R === 1) return null;
  return (u as bigint[][]).slice(0, R).map((c) => c.slice(0, n0));
}

import { ZX_mul } from './ZX.js';
import { centermodii } from './polarit2.js';
import { diviiround } from './gen3.js';
import { logintall } from './ispower.js';
import { cmbf_maxK } from './nffactor.js';

/** QX_factor.c exact division with an optional quotient coefficient bound.
 * Native callers require a normalized nonconstant divisor. The leading quotient
 * coefficient is deliberately not checked against B by the original routine.
 * @see Deviation: PARI bounded factor recombination adapters
 */
export function ZX_divides_i(x: bigint[], y: bigint[], B: bigint | null = null): bigint[] | null {
  const dx = x.length - 1, dy = y.length - 1, dz = dx - dy;
  if (dy < 1 || !y[dy]) throw new RangeError('ZX_divides_i requires a normalized nonconstant divisor');
  if (dz < 0) return null;
  const lead = y[dy]!, z = new Array<bigint>(dz + 1);
  if (x[dx]! % lead) return null;
  z[dz] = x[dx]! / lead;
  let i = dx - 1;
  for (; i >= dy; i--) {
    let s = x[i]!;
    for (let j = i - dy + 1; j <= i && j <= dz; j++) s -= z[j]! * y[i - j]!;
    if (s % lead) return null;
    s /= lead;
    if (B !== null && (s < 0n ? -s : s) > (B < 0n ? -B : B)) return null;
    z[i - dy] = s;
  }
  for (; i >= 0; i--) {
    let s = x[i]!;
    for (let j = 0; j <= i && j <= dz; j++) s -= z[j]! * y[i - j]!;
    if (s) return null;
  }
  return z;
}
/** Unbounded native integer-polynomial exact quotient; null when not divisible.
 * @see Deviation: PARI bounded factor recombination adapters
 */
export function ZX_divides(x: bigint[], y: bigint[]): bigint[] | null {
  return ZX_divides_i(x, y, null);
}

/** Native [flag, a, b, q^a, q^b] precision selection for the 31-bit trace window.
 * @see Deviation: PARI bounded factor recombination adapters
 */
export function cmbf_precs(q: bigint, A: bigint, B: bigint): [number, number, number, bigint, bigint] {
  if (q < 2n || q >= 1n << 31n || A <= 0n || B <= 0n)
    throw new RangeError('cmbf_precs requires 2 <= q < 2^31 and positive bounds');
  const d = Math.trunc(31 / Math.log2(Number(q)) - 1e-5);
  let b = logintall(B, q) + 1, a = b + d, flag = 1;
  let qa = q ** BigInt(a);
  if (qa <= A) {
    a = logintall(A, q) + 1;
    qa = q ** BigInt(a);
    b = a - d;
    flag = 0;
  }
  return [flag, a, b, qa, q ** BigInt(b)];
}

/** Native bounded subset recombination, including trace filters and factor order.
 * Returns [factors, associated modular factors, maximum subset size, done].
 * Inputs are copied; callers supply the native congruence and trace-bound contract.
 * @see Deviation: PARI bounded factor recombination adapters
 */
export function cmbf(
  pol: bigint[], famod: bigint[][], bound: bigint, p: bigint,
  a: number, b: number, klim: number
): [bigint[][], bigint[][][], number, boolean] {
  pol = pol.slice();
  // One-indexed temporary vectors retain the source's combination schedule.
  const fm = [[], ...famod.map(f => f.slice())];
  let count = famod.length, K = 1, maxK = cmbf_maxK(count);
  const pa = p ** BigInt(a), half = pa >> 1n, pb = p ** BigInt(b),
    window = p ** BigInt(a - b), windowHalf = window >> 1n;
  let lc: bigint | null = pol.at(-1)! < 0n ? -pol.at(-1)! : pol.at(-1)!;
  if (lc === 1n) lc = null;
  let scaled = lc === null ? pol : pol.map(c => c * lc!);
  const trace1 = [0n], trace2 = [0n], degrees = [0],
    indices = new Array<number>(count + 1), degreesBefore = new Array<number>(count + 1);
  degreesBefore[0] = 0;
  for (let i = 1; i <= count; i++) {
    const f = fm[i]!, d = f.length - 1;
    degrees[i] = d;
    let t1 = f[d - 1]!, t2 = t1 * t1;
    if (d > 1) t2 -= 2n * f[d - 2]!;
    t2 = residue(t2, pa);
    if (lc !== null) {
      t1 = residue(lc * t1, pa);
      t2 = residue(lc * lc * t2, pa);
    }
    trace1[i] = diviiround(t1, pb);
    trace2[i] = diviiround(t2, pb);
  }
  const factors: bigint[][] = [], groups: bigint[][][] = [];
  const primitive = (f: bigint[]) => { const c = content(f); return f.map(v => v / c); };
  nextK: for (;;) {
    if (K > maxK || 2 * K > count) break;
    indices[1] = 1;
    const allowedTrace = BigInt((K + 1) >> 1);
    let i = 1, degree = degrees[1]!;
    for (;;) {
      for (let j = i; j < K; j++) {
        degreesBefore[j] = degree;
        indices[j + 1] = indices[j]! + 1;
        degree += degrees[indices[j + 1]!]!;
      }
      if (degree <= klim) {
        const passesTrace = (traces: bigint[]) => {
          let t = traces[indices[1]!]!;
          for (let j = 2; j <= K; j++) {
            t += traces[indices[j]!]!;
            // Native Fl_add subtracts the modulus at most once.
            if (t >= window) t -= window;
          }
          if (t > windowHalf) t = window - t;
          return t <= allowedTrace;
        };
        if (passesTrace(trace1) && passesTrace(trace2)) {
          let trailing = lc;
          for (let j = 1; j <= K; j++) {
            let c = fm[indices[j]!]![0]!;
            if (trailing !== null) c *= trailing;
            trailing = centermodii(c, pa, half);
          }
          if (trailing && scaled[0]! % trailing === 0n) {
            let candidate: bigint[] | null = lc === null ? null : [lc];
            for (let j = 1; j <= K; j++) {
              const f = fm[indices[j]!]!;
              candidate = (candidate === null ? f : ZX_mul(candidate, f)).map(c => centermodii(c, pa, half));
            }
            const quotient = ZX_divides_i(scaled, candidate!, bound);
            if (quotient) {
              groups.push(indices.slice(1, K + 1).map(j => fm[j]!));
              const factor = primitive(candidate!);
              factors.push(factor);
              pol = lc === null ? quotient : quotient.map(c => c / factor.at(-1)!);
              for (let j = 1, dest = 1, index = 1; j <= count; j++) {
                if (index <= K && j === indices[index]) index++;
                else {
                  fm[dest] = fm[j]!;
                  trace1[dest] = trace1[j]!; trace2[dest] = trace2[j]!;
                  degrees[dest++] = degrees[j]!;
                }
              }
              count -= K;
              maxK = cmbf_maxK(count);
              if (count < 2 * K) break nextK;
              i = 1; degree = degrees[indices[1]!]!;
              bound = factor_bound(pol);
              if (lc !== null) lc = pol.at(-1)! < 0n ? -pol.at(-1)! : pol.at(-1)!;
              scaled = lc === null ? pol : pol.map(c => c * lc!);
              continue;
            }
          }
        }
      }
      i = K + 1;
      for (;;) {
        if (--i === 0) { K++; continue nextK; }
        if (++indices[i]! <= count - K + i) {
          degree = degreesBefore[i - 1]! + degrees[indices[i]!]!;
          if (degree <= klim) break;
        }
      }
    }
  }
  let done = true;
  if (pol.length > 1) {
    if (pol.at(-1)! < 0n) pol = pol.map(c => -c);
    if (count >= 2 * K) done = false;
    groups.push(fm.slice(1, count + 1));
    factors.push(pol);
  }
  return [factors, groups, maxK, done];
}

import { FpXV_prod, FpX_deriv } from './FpX.js';
import { FpX_red } from './ffinit.js';
import { Flx_factor, Flx_nbfact_by_degree, Flx_nbroots } from './FpX_factor.js';
import { Flx_normalize } from './Flx.js';
import { ZM_hnf_knapsack, ZM_hnf } from './hnf_snf.js';
import { ZM_mul } from './ZV.js';
import { FpM_image } from './alglin1.js';
import { ZpX_liftfact } from './Zp.js';
import { dbltor } from './kernel/none/mp_indep.js';
import { polsym_gen } from './polarit2.js';
import { factoru } from './ifactor.js';
import { ZX_deriv } from './ZX.js';
import { primitivePolynomial, type RationalPolynomialData } from './_rational_polynomial.js';
const abs = (x: bigint) => (x < 0n ? -x : x);
const primitive = (x: bigint[]) => {
  const c = x.reduce((g, c) => gcd(g, c), 0n);
  return x.map((v) => v / c);
};
const log2 = (x: bigint) => {
  const r = itor(x, 64);
  return Math.log2(Number(r.m)) + (r.e - 63);
};
/** Native chk_factors_get from QX_factor.c; see the documented stage preconditions.
 * @see Deviation: PARI integer polynomial factorization adapters
 */
export function chk_factors_get(
  lt: bigint | null,
  famod: bigint[][],
  c: bigint[],
  T: null,
  N: bigint
): bigint[] | bigint {
  const V = famod.filter((_, j) => c[j] !== 0n).map((f) => f.slice());
  if (lt !== null && V.length) V[0] = V[0]!.map((v) => v * lt);
  return FpXV_prod(V, N);
}
/** Native chk_factors from QX_factor.c; see the documented stage preconditions.
 * @see Deviation: PARI integer polynomial factorization adapters
 */
export function chk_factors(
  P: bigint[],
  M: bigint[][],
  bound: bigint,
  famod: bigint[][],
  pa: bigint
): bigint[][] | null {
  const piv = ZM_hnf_knapsack(M);
  if (!piv) return null;
  const r = piv.length;
  if (r < 2) throw new RangeError('chk_factors requires at least two retained factor columns');
  let pol = P,
    lt: bigint | null = abs(pol.at(-1)!);
  if (lt === 1n) lt = null;
  let scaled = lt === null ? pol : pol.map((c) => c * lt!);
  const out: bigint[][] = [],
    half = pa >> 1n;
  for (let i = 0; ; ) {
    const product = chk_factors_get(lt, famod, piv[i]!, null, pa) as bigint[];
    let y = product.map((c) => (abs(c) > half ? c - pa : c));
    const quotient = ZX_divides_i(scaled, y, bound);
    if (!quotient) return null;
    pol = quotient;
    if (lt !== null) y = primitive(y);
    out.push(y);
    if (++i >= r - 1) break;
    if (lt !== null) {
      pol = pol.map((c) => c / y.at(-1)!);
      lt = abs(pol.at(-1)!);
      scaled = pol.map((c) => c * lt!);
    } else scaled = pol;
  }
  out.push(primitive(pol));
  return out;
}
/** Native LLL_cmbf from QX_factor.c; see the documented stage preconditions.
 * @see Deviation: PARI integer polynomial factorization adapters
 */
export function LLL_cmbf(
  P: bigint[],
  famod: bigint[][],
  p: bigint,
  pa: bigint,
  bound: bigint,
  a: number,
  rec: number
): bigint[][] {
  const n0 = famod.length,
    degree = P.length - 1,
    logp = Math.log(Number(p)),
    LOGp2 = Math.LN2 / logp,
    b0 = Math.log(degree * 2) / logp,
    C = Math.ceil(Math.sqrt(n0 / 4)),
    Bnorm = dbltor(n0 * (C * C + n0 / 4) * 1.00001);
  let lc: bigint | null = abs(P.at(-1)!);
  if (lc === 1n) lc = null;
  const Br = root_bound(P) * (lc ?? 1n),
    logBr = log2(Br) * LOGp2;
  let CM = Array.from({ length: n0 }, (_, j) =>
    Array.from({ length: n0 }, (_, i) => (i === j ? BigInt(C) : 0n))
  );
  let TT: (bigint[] | null)[] = Array(n0).fill(null);
  for (let tmax = 0; ; tmax++) {
    const tnew = tmax + 1,
      r = CM.length,
      bmin = Math.ceil(b0 + tnew * logBr),
      old = CM;
    if (a <= bmin) {
      a = Math.ceil(bmin + 3 * logBr) + 1;
      let next = 1;
      while (next < a) next *= 2;
      a = next;
      pa = p ** BigInt(a);
      famod = ZpX_liftfact(P, famod, p, a);
      TT = Array(n0).fill(null);
    }
    const Tra = famod.map((f, i) => {
      const t = polsym_gen(f, TT[i]!, tnew, null, pa);
      TT[i] = t;
      return [t[tnew]! * (lc === null ? 1n : lc ** BigInt(tnew))];
    });
    let first = true,
      b = 0,
      delta = 0;
    for (;;) {
      const M = CM.map((c) => c.map((v) => v / BigInt(C))),
        half = pa >> 1n,
        T2 = ZM_mul([[], ...Tra.map((c) => [0n, ...c])], [[], ...M.map((c) => [0n, ...c])])
          .slice(1)
          .map((c) => c.slice(1).map((v) => centermodii(v, pa, half)));
      const expo = T2.flat().reduce(
        (e, v) => (v ? Math.max(e, abs(v).toString(2).length - 1) : e),
        -Infinity
      );
      if (first) {
        b = Math.max(bmin, Math.trunc((expo - Math.max(32, 0.4 * r)) * LOGp2));
        delta = a - b;
      } else {
        const b0 = Math.trunc(expo * LOGp2);
        if (b0 < b) b = b0;
        b = Math.max(b - delta, bmin);
        if (b - Math.trunc(delta / 2) < bmin) b = bmin;
      }
      const q = p ** BigInt(b),
        m = CM.map((c, j) => [...c, ...T2[j]!.map((v) => diviiround(v, q))]);
      if (first) {
        first = false;
        m.push([...Array<bigint>(n0).fill(0n), p ** BigInt(a - b)]);
      }
      const next = LLL_check_progress(Bnorm, n0, m, b === bmin);
      if (next === null) return [P.slice()];
      CM = next;
      if (b <= bmin) break;
    }
    const i = CM.length;
    if (i === r && CM.every((c, j) => c.every((v, k) => v === old[j]![k]))) {
      CM = old;
      continue;
    }
    if (FpM_image([[], ...CM.map((c) => [0n, ...c])], 27449n).length - 1 !== CM.length)
      CM = ZM_hnf(CM);
    if (i <= r && i * rec < n0) {
      const list = chk_factors(
        P,
        CM.map((c) => c.map((v) => v / BigInt(C))),
        bound,
        famod,
        pa
      );
      if (list) return list;
    }
  }
}
/** Native combine_factors from QX_factor.c; see the documented stage preconditions.
 * @see Deviation: PARI integer polynomial factorization adapters
 */
export function combine_factors(
  target: bigint[],
  famod: bigint[][],
  p: bigint,
  klim: number
): bigint[][] {
  let A = factor_bound(target);
  const B = BigInt(target.length - 1) * (abs(target.at(-1)!) * root_bound(target)) ** 2n,
    [, a, b, pa] = cmbf_precs(p, A, B);
  famod = ZpX_liftfact(target, famod, p, a);
  const [res, groups, maxK] = cmbf(target, famod, A, p, a, b, klim),
    last = res.length - 1;
  famod = groups[last]!;
  if (maxK > 0 && famod.length > 2 * maxK) {
    if (last) A = factor_bound(res[last]!);
    return [...res.slice(0, -1), ...LLL_cmbf(res[last]!, famod, p, pa, A, a, maxK)];
  }
  return res;
}
/** Native pick_prime from QX_factor.c; see the documented stage preconditions.
 * @see Deviation: PARI integer polynomial factorization adapters
 */
export function pick_prime(A: bigint[], fl: number): bigint {
  const degree = A.length - 1,
    lead = A.at(-1)!;
  let nmax = degree + 1,
    chosen = 0n,
    p = 1n;
  for (let np = 0; np < 7; ) {
    p = nextprime(p + 1n);
    if (lead % p === 0n) continue;
    const z = FpX_red(A, p);
    if (Flx_gcd(z, FpX_deriv(z, p), p).length !== 1) continue;
    let count: number;
    if (fl === 1) {
      count = Flx_nbroots(z, p);
      if (!count) return 0n;
    } else {
      const { D, nb } = Flx_nbfact_by_degree(z, p);
      count = nb;
      if (fl === 0 && count === 1) return 0n;
      if (fl !== 0 && D.findIndex((n, i) => i > 0 && n > 0) > fl) return 0n;
    }
    if (count < nmax) {
      nmax = count;
      chosen = p;
      if (degree > 100 && nmax < 5) break;
    }
    np++;
  }
  return chosen;
}
/** Native ZX_DDF_max from QX_factor.c; see the documented stage preconditions.
 * @see Deviation: PARI integer polynomial factorization adapters
 */
export function ZX_DDF_max(A: bigint[], dmax: number): bigint[][] {
  const p = pick_prime(A, dmax);
  if (!p) return [A.slice()];
  const modular = Flx_factor(Flx_normalize(FpX_red(A, p), p), p).map(([f]) => f);
  return combine_factors(A, modular, p, A.length - 2);
}
/** Native ZX_DDF from QX_factor.c; see the documented stage preconditions.
 * @see Deviation: PARI integer polynomial factorization adapters
 */
export function ZX_DDF(x: bigint[]): bigint[][] {
  let m = 0n;
  for (let i = 1; i < x.length; i++) if (x[i]) m = gcd(m, BigInt(i));
  const d = Number(m);
  let L = ZX_DDF_max(
    x.filter((_, i) => i % d === 0),
    0
  );
  if (d > 1) {
    const primes: bigint[] = [];
    for (const [p, e] of factoru(m)) for (let j = 0n; j < e; j++) primes.push(p);
    for (const p of primes.reverse())
      L = L.flatMap((f) => {
        const out = Array<bigint>((f.length - 1) * Number(p) + 1).fill(0n);
        f.forEach((v, i) => (out[i * Number(p)] = v));
        return ZX_DDF_max(out, 0);
      });
  }
  return L;
}
/** Native ZX_gcd_all from QX_factor.c; see the documented stage preconditions.
 * @see Deviation: PARI integer polynomial factorization adapters
 */
export function ZX_gcd_all(A: bigint[], B: bigint[]): [bigint[], bigint[]] {
  A = trimPolynomial(A);
  B = trimPolynomial(B);
  if (!A.length) return [B.slice(), []];
  if (!B.length) return [A.slice(), [1n]];
  const H = ZX_gcd(A, B),
    G = primitive(H),
    P = primitive(A);
  if (G.length === 1) return [H, P];
  const quotient = ZX_divides(P, G)!;
  return [H, quotient.map((c) => c * G.at(-1)!)];
}
/** Native ZX_squff from QX_factor.c; see the documented stage preconditions.
 * @see Deviation: PARI integer polynomial factorization adapters
 */
export function ZX_squff(f: bigint[]): [bigint[][], number[]] {
  if (f.at(-1)! < 0n) f = f.map((c) => -c);
  let [T, V] = ZX_gcd_all(f, ZX_deriv(f));
  const P: bigint[][] = [],
    E: number[] = [];
  for (let k = 1; ; k++) {
    const [W, next] = ZX_gcd_all(T, V);
    T = next;
    const dw = W.length - 1,
      dv = V.length - 1;
    if (!dw) {
      if (dv) {
        P.push(primitive(V));
        E.push(k);
      }
      break;
    }
    if (dw === dv) {
      for (;;) {
        const q = ZX_divides(T, V);
        if (q === null) break;
        k++;
        T = q;
      }
    } else {
      P.push(primitive(ZX_divides(V, W)!));
      E.push(k);
      V = W;
    }
  }
  return [P, E];
}
/** Native ZX_factor from QX_factor.c; see the documented stage preconditions.
 * @see Deviation: PARI integer polynomial factorization adapters
 */
export function ZX_factor(input: bigint[]): [bigint[], number][] {
  let T = trimPolynomial(input);
  if (!T.length) return [[[], 1]];
  let valuation = 0;
  while (T[valuation] === 0n) valuation++;
  T = T.slice(valuation);
  const [Q, E] = ZX_squff(T),
    out: [bigint[], number][] = [];
  Q.forEach((q, i) => ZX_DDF(q).forEach((f) => out.push([f, E[i]!])));
  if (valuation) out.push([[0n, 1n], valuation]);
  return out.sort(([a], [b]) => {
    if (a.length !== b.length) return a.length - b.length;
    for (let i = a.length - 1; i >= 0; i--) if (a[i] !== b[i]) return a[i]! < b[i]! ? -1 : 1;
    return 0;
  });
}
/** Native QX_factor from QX_factor.c; see the documented stage preconditions.
 * @see Deviation: PARI integer polynomial factorization adapters
 */
export function QX_factor(x: RationalPolynomialData): [bigint[], number][] {
  return ZX_factor(primitivePolynomial(x)[0]);
}
/** Native ZX_is_irred from QX_factor.c; see the documented stage preconditions.
 * @see Deviation: PARI integer polynomial factorization adapters
 */
export function ZX_is_irred(input: bigint[]): boolean {
  const x = trimPolynomial(input);
  if (x.length < 2) return false;
  if (x.length === 2) return true;
  if (x[0] === 0n || !ZX_is_squarefree(x)) return false;
  return ZX_DDF(x).length === 1;
}
