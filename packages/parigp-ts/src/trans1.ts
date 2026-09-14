import { atanhuu } from './trans2.js';
import { PariDomainError } from './matkermod.js';
/** Real exponential, pi and powering kernels from PARI basemath/trans1.c.
 * @see Deviation: PARI factorial real representation and transcendental dependencies
 */
import {
  itor,
  rtor,
  mulrr,
  mulir,
  divrr,
  divri,
  divir,
  divru,
  sqrr,
  sqrtr_abs as sqrtr,
  shiftr,
  addrr,
  subrr,
  addir,
  addrs,
  real_1,
  real_0_bit,
  absr,
  mulsr,
  subrs,
  real_0,
  type MpReal,
} from './qfb.js';
import { rtodbl, invr } from './kernel/none/mp_indep.js';
import { PariError } from './errors.js';
import { gen_powu_i } from './bb_group.js';
import { quadratic_prec_mask } from './Zp.js';
const nbits = (p: number) => Math.ceil(p / 64) * 64;
const truncate = (x: MpReal, p: number): MpReal => ({ ...x, p, m: x.m >> BigInt(x.p - p) });
const afffix = (x: MpReal, p: number) => rtor(x, Math.min(x.p, p));
/** Native exp(abs(x))-1 for a nonzero normalized real. Preserve its full exponent. */
export function exp1r_abs(x: MpReal): MpReal<bigint> {
  const l = x.p,
    a = x.e;
  let b = l;
  if (b + a <= 0) return shiftedReal({ ...x, s: x.s ? 1 : 0 }, 0n);
  const B = Math.floor(b / 3) + 64 + Math.floor(4096 / b);
  let d = a / 2,
    m = Math.trunc(d + Math.sqrt(d * d + B));
  if (m < -a * 0.1) m = 0;
  const log2 = Math.log2(Number(x.m >> BigInt(x.p - 64))) + (x.e - 63);
  d = m - log2 - 1 / Math.LN2;
  while (d <= 0) {
    d++;
    m++;
  }
  const L = l + nbits(m);
  b += m;
  let n = Math.trunc(b / d);
  if (n === 1) n = Math.trunc(b / (d + Math.log2(n + 1)));
  while (n * (d + Math.log2(n + 1)) < b) n++;
  const X = { ...shiftr(rtor(x, L), -m), s: 1 as const };
  let p2: MpReal;
  if (n === 1) p2 = X;
  else {
    let s = 0,
      l1 = nbits(Math.trunc(d + n + 16));
    p2 = real_1(L);
    for (let i = n; i >= 2; i--) {
      const p3 = divru(truncate(X, l1), i);
      const v = s - p3.e,
        q = Math.floor(v / 64);
      s = v - q * 64;
      l1 = Math.min(L, l1 + q * 64);
      const p1 = addrr(real_1(l1), i === n ? p3 : mulrr(p3, p2));
      p2 = rtor(p1, l1);
    }
    p2 = mulrr(X, p2);
  }
  let scale = 0n;
  for (let i = 1; i <= m; i++) {
    if (p2.p > L) p2 = truncate(p2, L);
    if (scale) {
      // At this exponent, native addsr(2,p2) returns a copy of p2. Keep its
      // binary scale separately before doubling can exceed Number precision.
      const squared = mulrr(p2, { ...p2 });
      scale = 2n * scale + BigInt(squared.e);
      p2 = { ...squared, e: 0 };
      if (scale >= 1n << 61n) throw new PariError('overflow in expo()');
    } else {
      p2 = p2.e < -L ? shiftr(p2, 1) : mulrr(p2, addrs(p2, 2));
      if (p2.e >= 2 ** 40 && p2.e > p2.p + 64) {
        scale = BigInt(p2.e);
        p2 = { ...p2, e: 0 };
      }
    }
  }
  return shiftedReal(afffix(p2, l), scale);
}
/** Native modlog2 keeps the C double quotient estimate, including above 2^53. */
function modlog2(x: MpReal): { residual: MpReal; shift: bigint } {
  const d = rtodbl(x),
    qd = (Math.abs(d) + Math.LN2 / 2) / Math.LN2;
  if (qd >= 2 ** 63) throw new PariError('overflow in expo()');
  const shift = BigInt(Math.trunc(qd)) * (d < 0 ? -1n : 1n);
  const residual = shift ? subrr(rtor(x, x.p + 64), mulir(shift, mplog2(x.p + 64))) : x;
  return { residual, shift };
}
function shiftedReal(x: MpReal, shift: bigint): MpReal<bigint> {
  const e = BigInt(x.e) + shift;
  if (e >= 1n << 61n || e < -(1n << 61n)) throw new PariError('overflow in expo()');
  return { ...x, e };
}
function mpexp_basecase(x: MpReal): MpReal<bigint> {
  const { residual: y, shift } = modlog2(x);
  if (!y.s) return shiftedReal(real_1(x.p), shift);
  const small = exp1r_abs(y);
  // Range reduction bounds this exponent well inside the safe Number range.
  let z = addrs({ ...small, e: Number(small.e) }, 1);
  if (y.s < 0) z = invr(z);
  if (shift && z.p > x.p) z = rtor(z, x.p);
  return shiftedReal(z, shift);
}
/** Native series/Newton branches; retain the full signed PARI exponent. */
export function mpexp(x: MpReal): MpReal<bigint> {
  if (!x.s) return shiftedReal(x.e >= 0 ? real_0_bit(x.e) : real_1(nbits(-x.e)), 0n);
  if (x.p <= 4224) return mpexp_basecase(x);
  const { residual, shift } = modlog2(x),
    l = x.p;
  if (!residual.s) return shiftedReal(real_1(l), shift);
  // The native cache is prepared before the logarithms in the Newton steps.
  mppi(l);
  let mask = quadratic_prec_mask(l + 64),
    p = 1;
  for (let i = 0; i < 12; i++) {
    p *= 2;
    if (mask & 1n) p--;
    mask >>= 1n;
  }
  const initial = mpexp_basecase(rtor(residual, nbits(p)));
  let a = rtor({ ...initial, e: Number(initial.e) }, l + 64);
  let input = addrs(residual, 1);
  if (input.p < l + 64) input = rtor(input, l + 64);
  let t: MpReal;
  for (;;) {
    p *= 2;
    if (mask & 1n) p--;
    mask >>= 1n;
    const prec = nbits(p);
    t = mulrr(truncate(a, prec), subrr(truncate(input, prec), logr_abs(truncate(a, prec))));
    if (mask === 1n) break;
    // affrr writes the current logical precision; later setprec re-exposes the
    // existing allocation. Preserve the unchanged trailing buffer words.
    const rounded = rtor(t, prec);
    const tailBits = a.p - prec,
      tailMask = (1n << BigInt(tailBits)) - 1n;
    a = { ...rounded, p: a.p, m: (rounded.m << BigInt(tailBits)) | (a.m & tailMask) };
  }
  return shiftedReal(rtor(t!, l), shift);
}
let cachedPi: MpReal | undefined;

/** Native constpi cache followed by precision conversion. */
export function mppi(p: number): MpReal {
  if (cachedPi && cachedPi.p >= p) return rtor(cachedPi, p);
  const n = Math.trunc(1 + p / 47.11041314);
  const split = (a: number, b: number): [bigint, bigint, bigint] => {
    if (b - a === 1) {
      const k = BigInt(a),
        P = a ? (6n * k - 5n) * (2n * k - 1n) * (1n - 6n * k) : 1n,
        Q = a ? k * k * k * 10939058860032000n : 1n,
        T = (13591409n + 545140134n * k) * P;
      return [P, Q, T];
    }
    const c = Math.floor((a + b) / 2),
      [lp, lq, lt] = split(a, c),
      [rp, rq, rt] = split(c, b);
    return [lp * rp, lq * rq, rq * lt + lp * rt];
  };
  const [P, Q, T] = split(0, n),
    C = 640320n;
  cachedPi = rtor(mulrr(divri(itor(Q * (C / 12n), p + 64), T), sqrtr(itor(C, p + 64))), p);
  return cachedPi;
}

/** PARI powru, with binary and sliding-window operation schedules. */
export function powru(x: MpReal, n: bigint): MpReal {
  if (n === 0n) return x.s ? real_1(x.p) : x.e >= 0 ? real_0_bit(x.e) : real_1(nbits(-x.e));
  return gen_powu_i(x, n, sqrr, mulrr);
}

/** Arithmetic/geometric mean of 1 and a positive real, with native gap cutoff. */
export function agm1r_abs(x: MpReal): MpReal {
  let a = shiftr(addrr(real_1(x.p), x), -1),
    b = sqrtr(x);
  for (;;) {
    const gap = subrr(b, a);
    if (!gap.s || gap.e - b.e < 5 - x.p) break;
    const previous = a;
    a = shiftr(addrr(previous, b), -1);
    b = sqrtr(mulrr(previous, b));
  }
  return afffix(a, x.p);
}

/** Native high-precision logarithm using Pi/(2 AGM), then binary rescaling. */
export function logagmr_abs(q: MpReal): MpReal {
  const p = q.p + 64,
    limit = Math.floor(p / 2);
  const Q = { ...shiftr(rtor(q, p), limit - q.e), s: 1 as const };
  const inverse = shiftr(invr(Q), 2);
  const y = divrr(shiftr(mppi(p), -1), agm1r_abs(inverse));
  return afffix(addrr(y, mulir(BigInt(q.e - limit), mplog2(p))), q.p);
}


export interface Abpq {
  a: bigint[];
  b: bigint[];
  p: bigint[];
  q: bigint[];
}
export interface AbpqResult { P: bigint; Q: bigint; B: bigint; T: bigint }

/** PARI abpq_init, with zero-based coefficient arrays including the final slot. */
export function abpq_init(n: number): Abpq {
  return { a: new Array(n + 1), b: new Array(n + 1), p: new Array(n + 1), q: new Array(n + 1) };
}
function T2(A: Abpq, n: number): bigint {
  return A.p[n]! * (A.a[n]! * A.b[n + 1]! * A.q[n + 1]! + A.b[n]! * A.a[n + 1]! * A.p[n + 1]!);
}
/** Native Haible/Papanikolaou binary split on [n1,n2), including small bases. */
export function abpq_sum(n1: number, n2: number, A: Abpq): AbpqResult {
  if (n2 - n1 === 1) return { P: A.p[n1]!, Q: A.q[n1]!, B: A.b[n1]!, T: A.a[n1]! * A.p[n1]! };
  if (n2 - n1 === 2) return {
    P: A.p[n1]! * A.p[n1 + 1]!, Q: A.q[n1]! * A.q[n1 + 1]!,
    B: A.b[n1]! * A.b[n1 + 1]!, T: T2(A, n1),
  };
  if (n2 - n1 === 3) {
    const q = A.q[n1 + 1]! * A.q[n1 + 2]!, b = A.b[n1 + 1]! * A.b[n1 + 2]!;
    return {
      P: A.p[n1]! * A.p[n1 + 1]! * A.p[n1 + 2]!, Q: A.q[n1]! * q,
      B: A.b[n1]! * b, T: A.p[n1]! * (b * q * A.a[n1]! + A.b[n1]! * T2(A, n1 + 1)),
    };
  }
  const n = Math.floor((n1 + n2) / 2), L = abpq_sum(n1, n, A), R = abpq_sum(n, n2, A);
  return { P: L.P * R.P, Q: L.Q * R.Q, B: L.B * R.B, T: R.B * R.Q * L.T + L.B * L.P * R.T };
}

function log2_split(p: number): MpReal {
  const u = atanhuu(1n, 26n, p), v = shiftr(atanhuu(1n, 4801n, p), 1), w = shiftr(atanhuu(1n, 8749n, p), 3);
  return addrr(mulir(18n, u), addrr({ ...v, s: -1 }, w));
}
let cachedLog2: MpReal | undefined;
/** Native logarithm cache, populated with one guard word. */
export function constlog2(p: number): MpReal {
  if (!cachedLog2 || cachedLog2.p < p) cachedLog2 = rtor(log2_split(p + 64), p);
  return cachedLog2;
}
export function mplog2(p: number): MpReal {
  return rtor(constlog2(p), p);
}
/** Native double estimate reads and converts the complete leading unsigned word. */
function dbllog2r(x: MpReal): number {
  return Math.log2(Number(x.m >> BigInt(x.p - 64))) + (x.e - 63);
}
function logr_aux(y: MpReal): MpReal {
  const L = y.p, d = -2 * dbllog2r(y);
  let k = Math.trunc(2 * (L / d));
  if (k % 2 === 0) k++;
  if (k < 3) return y;
  const y2 = sqrr(y), incs = Math.trunc(d);
  let s = 0, l1 = nbits(Math.trunc(d));
  let S = rtor(divru(real_1(l1), k), l1), T: MpReal;
  for (k -= 2; ; k -= 2) {
    T = mulrr(S, truncate(y2, l1));
    if (k === 1) break;
    const words = Math.floor((s + incs) / 64);
    s = (s + incs) % 64;
    l1 = Math.min(L, l1 + words * 64);
    S = rtor(addrr(divru(real_1(l1), k), T), l1);
  }
  return mulrr(y, addrs(T, 1));
}
/** Native logarithm series/AGM dispatch and progressive working precision.
 * @see Deviation: Shared native logarithms and Buchmann transcendental results
 */
export function logr_abs(X: MpReal): MpReal {
  // The original low-level helper assumes nonzero input; retain the typed guard.
  if (!X.s) throw new PariDomainError('logr_abs', 'argument', '=', '0');
  const p = X.p, u = X.m >> BigInt(p - 64);
  let EX = X.e, D: bigint;
  if (u > 12297829382473034410n) { EX++; D = (1n << BigInt(p)) - 1n - X.m; }
  else D = X.m - (1n << BigInt(p - 1));
  if (!D) return EX ? mulsr(EX, mplog2(p)) : real_0(p);
  const a = p - D.toString(2).length, lost = 64 * Math.floor(a / 64);
  let L = p + 64;
  const b = L - lost;
  if (b > 24 * a * Math.log2(L / 64 + 2) && p > 384) return logagmr_abs(X);
  const d = -a / 2;
  let m = Math.trunc(d + Math.sqrt(d * d + Math.trunc(b / 6)));
  if (m > b - a) m = b - a;
  if (m < 0.2 * a) m = 0;
  else L += nbits(m);
  let x = shiftr(rtor(absr(X), L), -EX);
  for (let k = 1; k <= m; k++) x = sqrtr(x);
  let y = logr_aux(divrr(subrs(x, 1), addrs(x, 1)));
  y = shiftr(y, m + 1);
  if (EX) y = addrr(y, mulsr(EX, mplog2(p + 64)));
  return afffix(y, EX ? p : p - lost);
}

import { sqrtr as nativeSqrtr, type MpComplex } from './qfb.js';
/** Native unsigned half-integral power; odd powers of negative reals become complex.
 * @see Deviation: PARI factorization scalar and bound adapters
 */
export function powruhalf(x: MpReal & { s: 0 | 1 }, n: bigint): MpReal;
export function powruhalf(x: MpReal, n: bigint): MpReal | MpComplex;
export function powruhalf(x: MpReal, n: bigint): MpReal | MpComplex {
  return n & 1n ? nativeSqrtr(powru(x, n)) : powru(x, n >> 1n);
}
