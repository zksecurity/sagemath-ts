/**
 * @module parigp-ts/qfrep
 * @description Representation numbers (theta series) of a positive definite
 * integral quadratic form: PARI's `qfrep0` / GP's `qfrep`.
 *
 * Direct port of the PARI sources vendored under `reference/pari/`:
 *
 * - `src/basemath/bibli1.c:1649-1650`   `qfrep0`
 * - `src/basemath/bibli1.c:1299-1462`   `minim0_dolll` (the Fincke-Pohst
 *                                       enumeration; we port the
 *                                       `min_VECSMALL` / `min_VECSMALL2`
 *                                       branches, which are the only ones
 *                                       `qfrep0` uses)
 * - `src/basemath/bibli1.c:1126-1134`   `minim_lll`  (LLL-reduce the Gram
 *                                       matrix first)
 * - `src/basemath/bibli1.c:1136-1155`   `forqfvec_init_dolll`
 * - `src/basemath/bibli1.c:1121-1125`   `err_minim`
 * - `src/basemath/alglin2.c:1566-1607`  `qfgaussred_positive` (Cholesky /
 *                                       LDL^T of the form; [GTM138, Algo 2.7.6])
 * - `src/basemath/lll.c:2683-2692`      `lllgramint` = `ZM_lll(x, 0.99,
 *                                       LLL_IM | LLL_GRAM)`
 *
 * GP documentation (`src/functions/linear_algebra/qfrep`):
 *
 * > `q` being a square and symmetric matrix with integer entries representing a
 * > positive definite quadratic form, count the vectors representing successive
 * > integers.
 * >  - If `flag = 0`, count all vectors.  Outputs the vector whose `i`-th entry,
 * >    `1 <= i <= B`, is half the number of vectors `v` such that `q(v) = i`.
 * >  - If `flag = 1`, count vectors of even norm.  Outputs the vector whose
 * >    `i`-th entry, `1 <= i <= B`, is half the number of vectors such that
 * >    `q(v) = 2i`.
 *
 * ## Exact arithmetic (deliberate deviation from PARI)
 *
 * PARI runs the Fincke-Pohst tree search in C `double`s (`minim0_dolll` keeps
 * `v`, `y`, `z`, `q` as `double*`) and compensates with a fudge factor
 * `BOUND = borne * (1 + 1e-10)`; the norm of each vector found is then recovered
 * by rounding, `(ulong)(p + 0.5)`.  That is *inexact*: for badly conditioned
 * forms a vector sitting exactly on the boundary can be gained or lost.
 *
 * This port performs the identical enumeration, in the identical order, with the
 * identical pruning conditions, but in **exact integer arithmetic**.  The
 * Cholesky data is rescaled so that every quantity in the inner loop is a
 * `bigint`; see {@link cholesky_scaled}.  Consequences:
 *
 *  - counts are exactly right, never off by the boundary cases `eps` is meant to
 *    absorb;
 *  - the `pari_err_PREC` that `minim0_dolll` raises when the *Cholesky* loses
 *    precision cannot happen here.  The `pari_err_PREC` checks that depend only
 *    on the size of `B` (`is_bigint(BORNE)` and `(long)BOUND != sBORNE`) ARE
 *    reproduced, so the accepted range of `B` is exactly PARI's.
 */

import { PariDomainError, PariTypeError, type ZM } from './matkermod.js';

export { PariDomainError, PariTypeError, type ZM } from './matkermod.js';

/* ------------------------------------------------------------------ */
/* Errors                                                              */
/* ------------------------------------------------------------------ */

/**
 * PARI `pari_err_PREC`.
 *
 * NOTE: the other PARI error kinds used here (`PariTypeError`,
 * `PariDomainError`, ...) currently live in `matkermod.ts` because `parigp-ts`
 * has no shared errors module; `pari_err_PREC` had no user yet.  When
 * `src/errors.ts` appears, this class should move there and be re-exported from
 * here for compatibility.
 */
export class PariPrecError extends Error {
  constructor(fun: string) {
    super(`precision too low in ${fun}`);
    this.name = 'PariPrecError';
  }
}

/** `bibli1.c:1121-1125` `err_minim`. */
function err_minim(): never {
  throw new PariDomainError('minim0', 'form', 'is not', 'positive definite');
}

/* ------------------------------------------------------------------ */
/* Small integer helpers                                               */
/* ------------------------------------------------------------------ */

function abs(a: bigint): bigint {
  return a < 0n ? -a : a;
}

function gcd(a: bigint, b: bigint): bigint {
  a = abs(a);
  b = abs(b);
  while (b) {
    const t = a % b;
    a = b;
    b = t;
  }
  return a;
}

function lcm(a: bigint, b: bigint): bigint {
  if (a === 0n || b === 0n) return 0n;
  return abs((a / gcd(a, b)) * b);
}

/** floor(a / b) for b > 0, correct for negative `a`. */
function floorDiv(a: bigint, b: bigint): bigint {
  const q = a / b;
  return a % b < 0n ? q - 1n : q;
}

/**
 * floor(sqrt(n)) for n >= 0 (PARI's `sqrtint`).
 *
 * NOTE: this duplicates `sqrti` in `qfb.ts`; it is kept local so that this
 * module does not have to depend on `qfb.ts` (which drags in `ifactor.ts` and
 * `ff.ts`).  The two should be merged once `parigp-ts` grows a shared integer
 * module.  Exported for its unit test only - do not re-export it from the
 * package index, `sqrti` is already there.
 */
export function isqrt(n: bigint): bigint {
  if (n < 0n) throw new RangeError('isqrt: negative argument');
  if (n < 2n) return n;
  // Newton, seeded from the bit length.
  let x = 1n << ((BigInt(n.toString(2).length) >> 1n) + 1n);
  for (;;) {
    const y = (x + n / x) >> 1n;
    if (y >= x) break;
    x = y;
  }
  return x;
}

/* ------------------------------------------------------------------ */
/* Rationals (only used outside the inner loop)                        */
/* ------------------------------------------------------------------ */

/** Normalised rational: `d > 0`, `gcd(|n|, d) = 1`. */
export interface Frac {
  n: bigint;
  d: bigint;
}

function frac(n: bigint, d: bigint): Frac {
  if (d === 0n) throw new RangeError('frac: zero denominator');
  if (d < 0n) {
    n = -n;
    d = -d;
  }
  const g = gcd(n, d);
  if (g > 1n) {
    n /= g;
    d /= g;
  }
  return { n, d };
}
const FZERO: Frac = { n: 0n, d: 1n };
function fromZ(n: bigint): Frac {
  return { n, d: 1n };
}
function fadd(a: Frac, b: Frac): Frac {
  return frac(a.n * b.d + b.n * a.d, a.d * b.d);
}
function fsub(a: Frac, b: Frac): Frac {
  return frac(a.n * b.d - b.n * a.d, a.d * b.d);
}
function fmul(a: Frac, b: Frac): Frac {
  return frac(a.n * b.n, a.d * b.d);
}
function fdiv(a: Frac, b: Frac): Frac {
  if (b.n === 0n) throw new RangeError('fdiv: division by zero');
  return frac(a.n * b.d, a.d * b.n);
}
function fsign(a: Frac): number {
  return a.n === 0n ? 0 : a.n < 0n ? -1 : 1;
}
/** sign(a - b) */
function fcmp(a: Frac, b: Frac): number {
  const l = a.n * b.d;
  const r = b.n * a.d;
  return l === r ? 0 : l < r ? -1 : 1;
}
/** Nearest integer, ties away from zero (PARI's `ground` on rationals). */
function fround(a: Frac): bigint {
  const twice = 2n * a.n;
  const d2 = 2n * a.d;
  // round(n/d) = floor((2n + d) / (2d)) rounds halves up; PARI rounds
  // half away from zero, which only differs on exact halves.
  if (a.n >= 0n) return floorDiv(twice + a.d, d2);
  return -floorDiv(-twice + a.d, d2);
}

/* ------------------------------------------------------------------ */
/* Matrix helpers.  Layout is PARI's: `A[j][i]` is the (i, j) entry.    */
/* ------------------------------------------------------------------ */

function zm_copy(A: ZM): ZM {
  return A.map((c) => c.slice());
}

function zm_identity(n: number): ZM {
  const u: ZM = [];
  for (let j = 0; j < n; j++) {
    const col = new Array<bigint>(n).fill(0n);
    col[j] = 1n;
    u.push(col);
  }
  return u;
}

/** `u^T * G * u` for square integer matrices (PARI's `qf_ZM_apply`). */
export function qf_ZM_apply(G: ZM, u: ZM): ZM {
  const n = G.length;
  // t = G * u  (column j of t = G * column j of u)
  const t: ZM = [];
  for (let j = 0; j < n; j++) {
    const uj = u[j]!;
    const col = new Array<bigint>(n).fill(0n);
    for (let k = 0; k < n; k++) {
      const c = uj[k]!;
      if (c === 0n) continue;
      const Gk = G[k]!;
      for (let i = 0; i < n; i++) col[i]! += c * Gk[i]!;
    }
    t.push(col);
  }
  // out = u^T * t
  const out: ZM = [];
  for (let j = 0; j < n; j++) {
    const tj = t[j]!;
    const col = new Array<bigint>(n).fill(0n);
    for (let i = 0; i < n; i++) {
      const ui = u[i]!;
      let s = 0n;
      for (let k = 0; k < n; k++) s += ui[k]! * tj[k]!;
      col[i] = s;
    }
    out.push(col);
  }
  return out;
}

/** Determinant of a square integer matrix (Bareiss, exact). */
export function ZM_det(A: ZM): bigint {
  const n = A.length;
  if (n === 0) return 1n;
  // work on a row-major copy m[i][j]
  const m: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const row = new Array<bigint>(n);
    for (let j = 0; j < n; j++) row[j] = A[j]![i]!;
    m.push(row);
  }
  let sign = 1n;
  let prev = 1n;
  for (let k = 0; k < n - 1; k++) {
    if (m[k]![k]! === 0n) {
      let p = -1;
      for (let i = k + 1; i < n; i++) {
        if (m[i]![k]! !== 0n) {
          p = i;
          break;
        }
      }
      if (p < 0) return 0n;
      const t = m[k]!;
      m[k] = m[p]!;
      m[p] = t;
      sign = -sign;
    }
    for (let i = k + 1; i < n; i++) {
      for (let j = k + 1; j < n; j++) {
        m[i]![j] = (m[i]![j]! * m[k]![k]! - m[i]![k]! * m[k]![j]!) / prev;
      }
      m[i]![k] = 0n;
    }
    prev = m[k]![k]!;
  }
  return sign * m[n - 1]![n - 1]!;
}

/* ------------------------------------------------------------------ */
/* qfgaussred_positive  (alglin2.c:1566-1607)                          */
/* ------------------------------------------------------------------ */

/**
 * Cholesky (LDL^T) decomposition of the positive definite matrix `a`,
 * `alglin2.c:1566-1607`.
 *
 * Returns the upper triangular matrix `b` in row-major, 0-indexed form:
 * `b[k][k]` is the `k`-th pivot `v_k` and `b[k][j]` (`j > k`) is the
 * coefficient `q_{kj}`, so that
 *
 *     q(x) = sum_k v_k * (x_k + sum_{j>k} q_{kj} x_j)^2 .
 *
 * Returns `null` if `a` is not positive definite (PARI returns `NULL`).
 *
 * PARI first attempts this in floating point (`forqfvec_init_dolll`) and only
 * falls back to the exact computation when that fails; we always take the
 * exact path.  Only the upper triangle of `a` is read, exactly as PARI does.
 */
export function qfgaussred_positive(a: ZM): Frac[][] | null {
  const n = a.length;
  const b: Frac[][] = [];
  for (let i = 0; i < n; i++) b.push(new Array<Frac>(n).fill(FZERO));
  for (let j = 0; j < n; j++) {
    for (let i = 0; i <= j; i++) b[i]![j] = fromZ(a[j]![i]!);
  }
  for (let k = 0; k < n; k++) {
    const p = b[k]![k]!;
    if (fsign(p) <= 0) return null; /* not positive definite */
    const invp = fdiv(fromZ(1n), p);
    const bk = b[k]!.slice(); /* PARI's `bk = row(b,k)`: a *copy* */
    for (let i = k + 1; i < n; i++) b[k]![i] = fmul(bk[i]!, invp);
    for (let i = k + 1; i < n; i++) {
      const c = bk[i]!;
      if (c.n === 0n) continue;
      for (let j = i; j < n; j++) {
        b[i]![j] = fsub(b[i]![j]!, fmul(c, b[k]![j]!));
      }
    }
  }
  return b;
}

/* ------------------------------------------------------------------ */
/* lllgramint  (lll.c:2690-2692)                                       */
/* ------------------------------------------------------------------ */

/**
 * LLL reduction of the *Gram matrix* `G` of a positive definite integral
 * quadratic form: PARI's `lllgramint(G) = ZM_lll(G, 0.99, LLL_IM | LLL_GRAM)`,
 * `lll.c:2690-2692`.
 *
 * Returns the unimodular `u` with `u^T G u` LLL-reduced (columns of `u` are the
 * coordinates of the new basis in terms of the old one), or `null` when the
 * form is not positive definite.
 *
 * DEVIATION: PARI's `ZM_lll` is a heavily engineered floating-point/`flatter`
 * hybrid (`lll.c`, ~2600 lines).  We use the textbook exact-rational LLL
 * (Cohen, *A Course in Computational Algebraic Number Theory*, Algorithm 2.6.3,
 * which is what `ZM_lll` computes, with the same default `delta = 0.99`,
 * `lll.c:474`).  Only the *speed* of {@link qfrep0} depends on `u`: the
 * representation numbers are invariant under any unimodular change of basis, so
 * this cannot change the result.  `qfrep0` additionally verifies that the `u`
 * returned here really is unimodular and falls back to the identity otherwise.
 */
export function lllgramint(G: ZM, delta: Frac = frac(99n, 100n)): ZM | null {
  const n = G.length;
  if (n === 0) return [];
  if (n === 1) return G[0]![0]! > 0n ? zm_identity(1) : null;

  const A = zm_copy(G); /* A = u^T G u, kept up to date */
  const u = zm_identity(n);
  const mu: Frac[][] = [];
  for (let i = 0; i < n; i++) mu.push(new Array<Frac>(n).fill(FZERO));
  const B = new Array<Frac>(n).fill(FZERO);

  const dot = (i: number, j: number): bigint => A[j]![i]!;

  /* b_k <- b_k - q b_l  (columns of u, and rows/cols k of A) */
  const transform = (k: number, l: number, q: bigint) => {
    if (q === 0n) return;
    const uk = u[k]!;
    const ul = u[l]!;
    for (let i = 0; i < n; i++) uk[i]! -= q * ul[i]!;
    /* A_{k,j} -= q A_{l,j} for all j, then A_{j,k} -= q A_{j,l} */
    for (let j = 0; j < n; j++) A[j]![k]! -= q * A[j]![l]!;
    for (let i = 0; i < n; i++) A[k]![i]! -= q * A[l]![i]!;
  };

  const swap = (k: number) => {
    const t = u[k]!;
    u[k] = u[k - 1]!;
    u[k - 1] = t;
    const tc = A[k]!;
    A[k] = A[k - 1]!;
    A[k - 1] = tc;
    for (let j = 0; j < n; j++) {
      const c = A[j]!;
      const s = c[k]!;
      c[k] = c[k - 1]!;
      c[k - 1] = s;
    }
  };

  B[0] = fromZ(dot(0, 0));
  if (fsign(B[0]!) <= 0) return null;
  let k = 1;
  let kmax = 0;
  const half = frac(1n, 2n);
  const mhalf = frac(-1n, 2n);

  const RED = (kk: number, l: number) => {
    if (fcmp(mu[kk]![l]!, half) <= 0 && fcmp(mu[kk]![l]!, mhalf) >= 0) return;
    const q = fround(mu[kk]![l]!);
    transform(kk, l, q);
    mu[kk]![l] = fsub(mu[kk]![l]!, fromZ(q));
    for (let i = 0; i < l; i++) {
      mu[kk]![i] = fsub(mu[kk]![i]!, fmul(fromZ(q), mu[l]![i]!));
    }
  };

  let guard = 0;
  const maxSteps = 1000000;
  while (k < n) {
    /* Safety valve: give up on reduction (u stays unimodular, so the caller
     * still gets a correct - merely less reduced - basis). */
    if (++guard > maxSteps) break;
    if (k > kmax) {
      kmax = k;
      for (let j = 0; j <= k; j++) {
        let s = fromZ(dot(k, j));
        for (let i = 0; i < j; i++) {
          s = fsub(s, fmul(fmul(mu[j]![i]!, mu[k]![i]!), B[i]!));
        }
        if (j < k) mu[k]![j] = fdiv(s, B[j]!);
        else {
          if (fsign(s) <= 0) return null; /* not positive definite */
          B[k] = s;
        }
      }
    }
    RED(k, k - 1);
    /* B[k] < (delta - mu[k][k-1]^2) B[k-1] ? */
    const m = mu[k]![k - 1]!;
    const lhs = B[k]!;
    const rhs = fmul(fsub(delta, fmul(m, m)), B[k - 1]!);
    if (fcmp(lhs, rhs) < 0) {
      /* SWAP(k) */
      const MU = mu[k]![k - 1]!;
      const BB = fadd(B[k]!, fmul(fmul(MU, MU), B[k - 1]!));
      if (fsign(BB) <= 0) return null;
      mu[k]![k - 1] = fdiv(fmul(MU, B[k - 1]!), BB);
      B[k] = fdiv(fmul(B[k - 1]!, B[k]!), BB);
      B[k - 1] = BB;
      swap(k);
      for (let j = 0; j <= k - 2; j++) {
        const t = mu[k - 1]![j]!;
        mu[k - 1]![j] = mu[k]![j]!;
        mu[k]![j] = t;
      }
      for (let i = k + 1; i <= kmax; i++) {
        const t = mu[i]![k]!;
        mu[i]![k] = fsub(mu[i]![k - 1]!, fmul(MU, t));
        mu[i]![k - 1] = fadd(t, fmul(mu[k]![k - 1]!, mu[i]![k]!));
      }
      k = Math.max(1, k - 1);
    } else {
      for (let l = k - 2; l >= 0; l--) RED(k, l);
      k++;
    }
  }
  return u;
}

/* ------------------------------------------------------------------ */
/* Scaled (all-integer) Cholesky data for the enumeration              */
/* ------------------------------------------------------------------ */

/**
 * Integer rescaling of the Cholesky data of `a`.
 *
 * With `v_k`, `q_{kj}` as in {@link qfgaussred_positive}, put
 * `d_0 = 1`, `d_k = v_1 ... v_k` (the leading principal minors of `a`, hence
 * integers) and
 *
 *     P[k][j] = d_k q_{kj}   (integers, by Cramer's rule)
 *     W_k     = d_k (x_k + sum_{j>k} q_{kj} x_j) = d_k x_k + sum_{j>k} P[k][j] x_j
 *
 * so that `v_k (x_k + z_k)^2 = W_k^2 / (d_{k-1} d_k)` and, with
 * `C = lcm_k(d_{k-1} d_k)` and `c_k = C / (d_{k-1} d_k)`,
 *
 *     C * q(x) = sum_k c_k W_k^2 .
 *
 * Every quantity in the enumeration is then an integer: PARI's `y[l]`, `z[l]`
 * and `BOUND` become `Y[l] = C*y[l]`, `Z[l] = d_l*z[l]` and `C*borne`.
 */
interface ScaledCholesky {
  n: number;
  /** `d[k]` for `k = 0..n`, `d[0] = 1`. */
  d: bigint[];
  /** `c[k]` for `k = 1..n`. */
  c: bigint[];
  /** `P[k][j] = d_k q_{kj}` for `1 <= k < j <= n`. */
  P: bigint[][];
  /** the common denominator `C`. */
  C: bigint;
}

function cholesky_scaled(a: ZM): ScaledCholesky | null {
  const n = a.length;
  const b = qfgaussred_positive(a);
  if (!b) return null;
  const d = new Array<bigint>(n + 1).fill(1n);
  let dk: Frac = fromZ(1n);
  for (let k = 1; k <= n; k++) {
    dk = fmul(dk, b[k - 1]![k - 1]!);
    if (dk.d !== 1n) {
      /* cannot happen for an integral form: d_k is a principal minor */
      throw new Error('qfrep: internal error, non-integral principal minor');
    }
    d[k] = dk.n;
  }
  const P: bigint[][] = [];
  for (let k = 0; k <= n; k++) P.push(new Array<bigint>(n + 1).fill(0n));
  for (let k = 1; k <= n; k++) {
    for (let j = k + 1; j <= n; j++) {
      const t = fmul(fromZ(d[k]!), b[k - 1]![j - 1]!);
      if (t.d !== 1n) {
        throw new Error('qfrep: internal error, non-integral Cholesky coefficient');
      }
      P[k]![j] = t.n;
    }
  }
  let C = 1n;
  for (let k = 1; k <= n; k++) C = lcm(C, d[k - 1]! * d[k]!);
  const c = new Array<bigint>(n + 1).fill(0n);
  for (let k = 1; k <= n; k++) c[k] = C / (d[k - 1]! * d[k]!);
  return { n, d, c, P, C };
}

/* ------------------------------------------------------------------ */
/* minim0 (min_VECSMALL / min_VECSMALL2)   bibli1.c:1299-1462          */
/* ------------------------------------------------------------------ */

/** Check the input is a square matrix of integers (`RgM_is_ZM`). */
function check_ZM_square(a: ZM, fun: string): number {
  if (!Array.isArray(a)) throw new PariTypeError(fun, 'not a t_MAT');
  const n = a.length;
  for (const col of a) {
    if (!Array.isArray(col) || col.length !== n) {
      throw new PariTypeError(fun, 'not a square t_MAT');
    }
    for (const e of col) {
      if (typeof e !== 'bigint') throw new PariTypeError(fun, 'not a t_INT');
    }
  }
  return n;
}

/**
 * `minim0(a, BORNE, 0, min_VECSMALL | min_VECSMALL2)`, `bibli1.c:1299-1462`.
 *
 * @param a     Gram matrix, PARI column-major (`a[j][i]` = entry `(i,j)`).
 * @param sBORNE bound `B`, already floored, `> 0`.
 * @param even  `true` for `min_VECSMALL2` (`flag & 1`).
 * @returns `L`, 0-indexed: `L[i-1]` = half the number of `v` with `q(v) = i`
 *          (resp. `q(v) = 2i`), `1 <= i <= B`.
 */
function minim0_vecsmall(a: ZM, sBORNE: bigint, even: boolean): bigint[] {
  const n = check_ZM_square(a, 'qfminim');

  /* bibli1.c:1324-1330: allocate L of length B, then double the search bound
   * for min_VECSMALL2. */
  let bound = sBORNE;
  if (even) bound <<= 1n;

  /* bibli1.c:1370-1371: PARI computes BOUND = sBORNE * (1 + eps) as a double
   * and bails out when the rounded value no longer represents sBORNE.  We keep
   * the check verbatim so the accepted range of B is exactly PARI's, even
   * though our own arithmetic is exact.  (PARI runs it after allocating L;
   * for any B large enough to fail here, allocating L is hopeless anyway, and
   * running the check first yields the more informative error.) */
  {
    const eps = 1e-10;
    const BOUND = Number(bound) * (1 + eps);
    if (!Number.isFinite(BOUND) || BigInt(Math.trunc(BOUND)) !== bound) {
      throw new PariPrecError('qfminim');
    }
  }

  const L = new Array<bigint>(Number(sBORNE)).fill(0n);
  if (n === 0) return L;

  /* bibli1.c:1345 forqfvec_init_dolll: LLL-reduce first (dolll = 1). */
  let A = a;
  {
    const u = lllgramint(a);
    if (u === null) err_minim(); /* minim_lll / qfgaussred_positive failure */
    /* Safety net: the enumeration is only correct if u is unimodular.  It is
     * by construction, but a wrong u would silently produce wrong counts, so
     * we check and fall back to the unreduced form rather than risk that. */
    const det = ZM_det(u);
    if (det === 1n || det === -1n) A = qf_ZM_apply(a, u);
  }

  const S = cholesky_scaled(A);
  if (!S) err_minim();

  const { d, c, P, C } = S;
  const CB = C * bound;

  /* PARI's arrays are 1-indexed; so are ours. */
  const x = new Array<bigint>(n + 2).fill(0n);
  const Y = new Array<bigint>(n + 2).fill(0n); /* Y[k] = C * y[k] */
  const Z = new Array<bigint>(n + 2).fill(0n); /* Z[k] = d_k * z[k] */

  let k = n;
  Y[n] = 0n;
  Z[n] = 0n;
  /* x[n] = floor(sqrt(BOUND / v[n])) */
  x[n] = isqrt(CB / c[n]!) / d[n]!;

  for (;;) {
    do {
      if (k > 1) {
        const l = k - 1;
        /* z[l] = sum_{j>=k} q[l][j] x[j] */
        let zl = 0n;
        const Pl = P[l]!;
        for (let j = k; j <= n; j++) {
          const xj = x[j]!;
          if (xj !== 0n) zl += Pl[j]! * xj;
        }
        Z[l] = zl;
        const W = d[k]! * x[k]! + Z[k]!;
        Y[l] = Y[k]! + c[k]! * W * W;
        /* x[l] = floor(sqrt((BOUND - y[l]) / v[l]) - z[l]) */
        const T = CB - Y[l]!;
        const Wmax = isqrt(T / c[l]!);
        x[l] = floorDiv(Wmax - zl, d[l]!);
        k = l;
      }
      for (;;) {
        const W = d[k]! * x[k]! + Z[k]!;
        if (Y[k]! + c[k]! * W * W <= CB) break;
        k++;
        if (k > n) {
          /* unreachable: the enumeration always stops at the zero vector */
          throw new Error('qfrep: internal error, enumeration overflow');
        }
        x[k]! -= 1n;
      }
    } while (k > 1);

    /* bibli1.c:1398: if (!x[1] && y[1] <= eps) break; */
    if (x[1]! === 0n && Y[1]! === 0n) break;

    const W = d[1]! * x[1]! + Z[1]!;
    const p = Y[1]! + c[1]! * W * W; /* = C * norm(x) */
    if (p % C !== 0n) throw new Error('qfrep: internal error, non-integral norm');
    const norm = p / C;
    /* bibli1.c:1439-1446 */
    if (even) {
      if ((norm & 1n) === 0n) {
        const i = Number(norm >> 1n) - 1;
        L[i]! += 1n;
      }
    } else {
      const i = Number(norm) - 1;
      L[i]! += 1n;
    }

    x[1]! -= 1n; /* the `for(;;x[1]--)` step */
  }
  return L;
}

/* ------------------------------------------------------------------ */
/* Public entry points                                                 */
/* ------------------------------------------------------------------ */

/** Anything we accept where PARI wants a `t_INT` bound. */
export type BoundLike = bigint | number;

/**
 * `qfrep0(a, borne, flag)`, `bibli1.c:1649-1650`:
 *
 *     GEN qfrep0(GEN a, GEN borne, long flag)
 *     { return minim0(a, borne, gen_0, (flag & 1)? min_VECSMALL2: min_VECSMALL); }
 *
 * `a` must be a square symmetric integral matrix representing a positive
 * definite quadratic form (PARI reads only its upper triangle in the Cholesky
 * step, exactly as we do).
 *
 * @param a     Gram matrix in PARI column-major layout: `a[j][i]` is the
 *              `(i, j)` entry.  Use `zm_from_rows` from `matkermod.ts` to build
 *              one from rows (identical for symmetric input).
 * @param borne the bound `B`; non-integral values are floored (`gfloor`).
 * @param flag  bit 0: count only vectors of even norm and halve the norms.
 * @returns a 0-indexed array of length `max(B, 0)`; entry `i-1` is **half** the
 *          number of `v != 0` with `q(v) = i` (`flag & 1 == 0`) resp.
 *          `q(v) = 2i` (`flag & 1 == 1`).
 *
 * @example
 * // GP: q = [2,1;1,3]; qfrep(q, 5) -> Vecsmall([0, 1, 2, 0, 0])
 * qfrep0([[2n, 1n], [1n, 3n]], 5) // => [0n, 1n, 2n, 0n, 0n]
 */
export function qfrep0(a: ZM, borne: BoundLike, flag = 0): bigint[] {
  let sBORNE: bigint;
  if (typeof borne === 'number') {
    if (!Number.isFinite(borne)) throw new PariTypeError('minim0', 'not a t_INT');
    sBORNE = BigInt(Math.floor(borne));
  } else if (typeof borne === 'bigint') {
    sBORNE = borne;
  } else {
    throw new PariTypeError('minim0', 'not a t_INT');
  }
  /* bibli1.c:1315: is_bigint(BORNE) => pari_err_PREC("qfminim") */
  if (sBORNE >= 1n << 63n || sBORNE < -(1n << 63n)) throw new PariPrecError('qfminim');
  if (sBORNE < 0n) sBORNE = 0n;
  /* bibli1.c:1326: sBORNE <= 0 => empty Vecsmall (no type checks are run) */
  if (sBORNE <= 0n) return [];
  return minim0_vecsmall(a, sBORNE, (flag & 1) !== 0);
}

/**
 * GP's `qfrep(q, B, {flag = 0})`.
 *
 * Bit 0 of `flag` is `qfrep0`'s flag.  cypari2 (`cypari2/gen.pyx:4222-4251`)
 * additionally uses bit 1 to choose between a `t_VEC` and a `t_VECSMALL`
 * result; both are the same JavaScript array here, so bit 1 is accepted and
 * ignored.  Any other bit is a `pari_err_FLAG` in GP, but PARI itself only ever
 * masks with 1, which is what we do.
 */
export function qfrep(q: ZM, B: BoundLike, flag = 0): bigint[] {
  return qfrep0(q, B, flag & 1);
}
