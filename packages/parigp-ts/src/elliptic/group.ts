/**
 * @module parigp-ts/elliptic/group
 * @description Group structure functions for elliptic curves over finite fields
 *
 * Port of PARI/GP functions from:
 * - elliptic.c:6332-6388 (ellcard)
 * - FpE.c:1470-1475 (Fp_ellgroup)
 * - FpE.c:406-423 (FpE_order)
 * - bb_group.c:986-1048 (gen_ellgroup)
 *
 * This module implements:
 * - ellcard(E): Compute the cardinality (number of points) on E(Fq)
 * - ellgroup(E): Determine the group structure [d1] or [d1, d2]
 * - ellgenerators(E): Find generators for the group
 * - ellorder(E, P): Compute the order of a point P
 * - trace_of_frobenius(E): Compute a_p = p + 1 - #E(Fp)
 */

import {
  Fp_add,
  Fp_center,
  Fp_div,
  Fp_double,
  Fp_inv,
  Fp_mul,
  Fp_mulu,
  Fp_neg,
  Fp_pow,
  Fp_red,
  Fp_sqr,
  Fp_sqrt,
  Fp_sub,
  gcd,
  kronecker,
  xgcd,
} from '../ff.js';
import { Z_factor } from '../ifactor.js';
import { ellweilpairing } from './advanced.js';

/**
 * Elliptic curve representation for PARI-style functions.
 *
 * We use short Weierstrass form: y^2 = x^3 + a4*x + a6 over Fp
 */
export interface EllipticCurveFp {
  /** Coefficient a4 in y^2 = x^3 + a4*x + a6 */
  a4: bigint;
  /** Coefficient a6 in y^2 = x^3 + a4*x + a6 */
  a6: bigint;
  /** Prime field characteristic p */
  p: bigint;
  /** Cached cardinality (if computed) */
  _card?: bigint;
  /** Cached group structure (if computed) */
  _group?: bigint[];
  /**
   * Cached group exponent `m` produced by `gen_ellgroup` (PARI's `*pm`
   * out-parameter); `gen_ellgens` needs it.
   */
  _m?: bigint;
  /** Cached generators (if computed) */
  _generators?: EllipticPointFp[];
}

/**
 * Point on an elliptic curve over Fp.
 *
 * PARI represents points as [x, y] for finite points
 * and [0] (or a special flag) for the point at infinity.
 */
export interface EllipticPointFp {
  /** x-coordinate (null for point at infinity) */
  x: bigint | null;
  /** y-coordinate (null for point at infinity) */
  y: bigint | null;
  /** Whether this is the point at infinity */
  isInfinity: boolean;
}

/**
 * Create the point at infinity.
 *
 * Reference: PARI pariinl.h:2753 - ellinf()
 */
export function ellinf(): EllipticPointFp {
  return { x: null, y: null, isInfinity: true };
}

/**
 * Create a finite point.
 *
 * @param x - x-coordinate
 * @param y - y-coordinate
 */
export function ellpoint(x: bigint, y: bigint): EllipticPointFp {
  return { x, y, isInfinity: false };
}

/**
 * Check if a point is the point at infinity.
 *
 * Reference: PARI pariinl.h:2754 - ell_is_inf()
 */
export function ell_is_inf(P: EllipticPointFp): boolean {
  return P.isInfinity;
}

/**
 * Check equality of two points.
 */
export function ellequal(P: EllipticPointFp, Q: EllipticPointFp): boolean {
  if (P.isInfinity && Q.isInfinity) return true;
  if (P.isInfinity || Q.isInfinity) return false;
  return P.x === Q.x && P.y === Q.y;
}

// ============================================================================
// Modular Arithmetic Helpers
//
// All Fp_* primitives are the ones from `../ff.ts` (the package's single
// port of pariinl.h/arith1.c).  This module used to carry a second, subtly
// different copy -- in particular a `Fp_sqrt` without PARI's smallest-root
// normalisation (arith1.c:1277).
// ============================================================================

/**
 * Modular reduction ensuring positive result (PARI `Fp_red`/`modii`).
 */
function mod(a: bigint, p: bigint): bigint {
  return Fp_red(a, p);
}

/**
 * Modular tripling.
 */
function Fp_triple(a: bigint, p: bigint): bigint {
  return mod(a * 3n, p);
}

/**
 * Integer square root (floor). PARI `sqrti`.
 */
function isqrt(n: bigint): bigint {
  if (n < 0n) throw new Error('Square root of negative number');
  if (n < 2n) return n;

  let x = n;
  let y = (x + 1n) / 2n;

  while (y < x) {
    x = y;
    y = (x + n / x) / 2n;
  }

  return x;
}

/**
 * Floor division (BigInt `/` truncates towards zero).
 */
function floorDiv(a: bigint, b: bigint): bigint {
  const q = a / b;
  return a % b !== 0n && a < 0n !== b < 0n ? q - 1n : q;
}

/**
 * LCM of two integers.
 */
function lcm(a: bigint, b: bigint): bigint {
  if (a === 0n || b === 0n) return 0n;
  const absA = a < 0n ? -a : a;
  const absB = b < 0n ? -b : b;
  return (absA / gcd(absA, absB)) * absB;
}

/**
 * p-adic valuation of n (PARI `Z_pval`).
 */
function valuation(n: bigint, q: bigint): bigint {
  if (n === 0n) return -1n; // Convention
  let v = 0n;
  let m = n < 0n ? -n : n;
  while (m % q === 0n) {
    v++;
    m /= q;
  }
  return v;
}

/**
 * Solve x = a1 (mod b1), x = a2 (mod b2).
 *
 * Reference: PARI arith2.c - Z_chinese_all (returns the new modulus too)
 */
function Z_chinese_all(a1: bigint, a2: bigint, b1: bigint, b2: bigint): [bigint, bigint] {
  const d = gcd(b1, b2);
  const L = (b1 / d) * b2;
  if (mod(a2 - a1, d) !== 0n) {
    throw new Error('Z_chinese_all: inconsistent congruences');
  }
  if (d === b2) return [mod(a1, L), L];
  const m = b2 / d;
  const [, u] = xgcd(mod(b1 / d, m), m);
  const t = mod(((a2 - a1) / d) * u, m);
  return [mod(a1 + b1 * t, L), L];
}

// ============================================================================
// Point Arithmetic
// ============================================================================

/**
 * Negate a point.
 *
 * For y^2 = x^3 + a4*x + a6: -P = (x, -y)
 *
 * Reference: PARI FpE.c - FpE_neg
 */
export function FpE_neg(P: EllipticPointFp, p: bigint): EllipticPointFp {
  if (P.isInfinity) return ellinf();
  return ellpoint(P.x!, Fp_neg(P.y!, p));
}

/**
 * Add two points on an elliptic curve over Fp.
 *
 * Uses the standard chord-and-tangent formulas for short Weierstrass form.
 *
 * Reference: PARI FpE.c - FpE_add
 */
export function FpE_add(
  P: EllipticPointFp,
  Q: EllipticPointFp,
  a4: bigint,
  p: bigint
): EllipticPointFp {
  // P + O = P
  if (P.isInfinity) return Q;
  // O + Q = Q
  if (Q.isInfinity) return P;

  const x1 = P.x!;
  const y1 = P.y!;
  const x2 = Q.x!;
  const y2 = Q.y!;

  let slope: bigint;

  if (x1 === x2) {
    // Same x-coordinate
    if (Fp_add(y1, y2, p) === 0n) {
      // P + (-P) = O (y1 = -y2)
      return ellinf();
    }
    // P = Q, use doubling formula
    // slope = (3*x1^2 + a4) / (2*y1)
    const num = Fp_add(Fp_triple(Fp_sqr(x1, p), p), mod(a4, p), p);
    const den = Fp_double(y1, p);
    slope = Fp_div(num, den, p);
  } else {
    // General case
    // slope = (y2 - y1) / (x2 - x1)
    slope = Fp_div(Fp_sub(y2, y1, p), Fp_sub(x2, x1, p), p);
  }

  // x3 = slope^2 - x1 - x2
  const x3 = Fp_sub(Fp_sub(Fp_sqr(slope, p), x1, p), x2, p);

  // y3 = slope * (x1 - x3) - y1
  const y3 = Fp_sub(Fp_mul(slope, Fp_sub(x1, x3, p), p), y1, p);

  return ellpoint(x3, y3);
}

/**
 * Double a point on an elliptic curve over Fp.
 *
 * Reference: PARI FpE.c - FpE_dbl
 */
export function FpE_dbl(P: EllipticPointFp, a4: bigint, p: bigint): EllipticPointFp {
  if (P.isInfinity) return ellinf();

  const x1 = P.x!;
  const y1 = P.y!;

  // Check for point of order 2
  if (y1 === 0n) return ellinf();

  // slope = (3*x1^2 + a4) / (2*y1)
  const num = Fp_add(Fp_triple(Fp_sqr(x1, p), p), mod(a4, p), p);
  const den = Fp_double(y1, p);
  const slope = Fp_div(num, den, p);

  // x3 = slope^2 - 2*x1
  const x3 = Fp_sub(Fp_sqr(slope, p), Fp_double(x1, p), p);

  // y3 = slope * (x1 - x3) - y1
  const y3 = Fp_sub(Fp_mul(slope, Fp_sub(x1, x3, p), p), y1, p);

  return ellpoint(x3, y3);
}

// ============================================================================
// Jacobian Coordinate Operations (internal to FpE_mul)
// ============================================================================

interface JacobianPointFp {
  X: bigint;
  Y: bigint;
  Z: bigint;
}

function ellinf_FpJ(): JacobianPointFp {
  return { X: 1n, Y: 1n, Z: 0n };
}

function FpE_to_FpJ(P: EllipticPointFp): JacobianPointFp {
  if (P.isInfinity) {
    return ellinf_FpJ();
  }
  return { X: P.x!, Y: P.y!, Z: 1n };
}

function FpJ_to_FpE(P: JacobianPointFp, p: bigint): EllipticPointFp {
  if (P.Z === 0n) {
    return ellinf();
  }
  const zInv = Fp_inv(P.Z, p);
  const z2 = Fp_sqr(zInv, p);
  const z3 = Fp_mul(zInv, z2, p);
  return ellpoint(Fp_mul(P.X, z2, p), Fp_mul(P.Y, z3, p));
}

function FpJ_dbl(P: JacobianPointFp, a4: bigint, p: bigint): JacobianPointFp {
  if (P.Z === 0n) {
    return ellinf_FpJ();
  }

  const X1 = P.X;
  const Y1 = P.Y;
  const Z1 = P.Z;

  const XX = Fp_sqr(X1, p);
  const YY = Fp_sqr(Y1, p);
  const YYYY = Fp_sqr(YY, p);
  const ZZ = Fp_sqr(Z1, p);

  const S = Fp_double(Fp_sub(Fp_sqr(Fp_add(X1, YY, p), p), Fp_add(XX, YYYY, p), p), p);
  const M = Fp_add(Fp_triple(XX, p), Fp_mul(a4, Fp_sqr(ZZ, p), p), p);
  const T = Fp_sub(Fp_sqr(M, p), Fp_double(S, p), p);

  const X3 = T;
  const Y3 = Fp_sub(Fp_mul(M, Fp_sub(S, T, p), p), Fp_mulu(YYYY, 8, p), p);
  const Z3 = Fp_sub(Fp_sqr(Fp_add(Y1, Z1, p), p), Fp_add(YY, ZZ, p), p);

  return { X: X3, Y: Y3, Z: Z3 };
}

function FpJ_add(P: JacobianPointFp, Q: JacobianPointFp, a4: bigint, p: bigint): JacobianPointFp {
  if (Q.Z === 0n) {
    return { X: P.X, Y: P.Y, Z: P.Z };
  }
  if (P.Z === 0n) {
    return { X: Q.X, Y: Q.Y, Z: Q.Z };
  }

  const X1 = P.X;
  const Y1 = P.Y;
  const Z1 = P.Z;
  const X2 = Q.X;
  const Y2 = Q.Y;
  const Z2 = Q.Z;

  const Z1Z1 = Fp_sqr(Z1, p);
  const Z2Z2 = Fp_sqr(Z2, p);
  const U1 = Fp_mul(X1, Z2Z2, p);
  const U2 = Fp_mul(X2, Z1Z1, p);
  const S1 = Fp_mul(Y1, Fp_mul(Z2, Z2Z2, p), p);
  const S2 = Fp_mul(Y2, Fp_mul(Z1, Z1Z1, p), p);
  const H = Fp_sub(U2, U1, p);
  const r = Fp_double(Fp_sub(S2, S1, p), p);

  if (H === 0n) {
    if (r === 0n) {
      return FpJ_dbl(P, a4, p);
    }
    return ellinf_FpJ();
  }

  const I = Fp_sqr(Fp_double(H, p), p);
  const J = Fp_mul(H, I, p);
  const V = Fp_mul(U1, I, p);
  const W = Fp_sub(Fp_sqr(r, p), Fp_add(J, Fp_double(V, p), p), p);

  const X3 = W;
  const Y3 = Fp_sub(Fp_mul(r, Fp_sub(V, W, p), p), Fp_double(Fp_mul(S1, J, p), p), p);
  const Z3 = Fp_mul(Fp_sub(Fp_sqr(Fp_add(Z1, Z2, p), p), Fp_add(Z1Z1, Z2Z2, p), p), H, p);

  return { X: X3, Y: Y3, Z: Z3 };
}

function gen_pow_FpJ(P: JacobianPointFp, n: bigint, a4: bigint, p: bigint): JacobianPointFp {
  if (n === 0n || P.Z === 0n) {
    return ellinf_FpJ();
  }

  let R = ellinf_FpJ();
  let Q: JacobianPointFp = { X: P.X, Y: P.Y, Z: P.Z };

  let remaining = n;
  while (remaining > 0n) {
    if ((remaining & 1n) === 1n) {
      R = FpJ_add(R, Q, a4, p);
    }
    Q = FpJ_dbl(Q, a4, p);
    remaining >>= 1n;
  }

  return R;
}

/**
 * Scalar multiplication [n]P using double-and-add.
 *
 * Reference: PARI FpE.c - FpE_mul
 */
export function FpE_mul(P: EllipticPointFp, n: bigint, a4: bigint, p: bigint): EllipticPointFp {
  if (n === 0n || P.isInfinity) return ellinf();

  if (n < 0n) {
    n = -n;
    P = FpE_neg(P, p);
  }

  if (n === 1n) {
    return P.isInfinity ? ellinf() : ellpoint(P.x!, P.y!);
  }

  if (n === 2n) {
    return FpE_dbl(P, a4, p);
  }

  const Q = gen_pow_FpJ(FpE_to_FpJ(P), n, a4, p);
  return FpJ_to_FpE(Q, p);
}

// ============================================================================
// Factorization
// ============================================================================

/**
 * Factor a positive integer.
 *
 * Delegates to the package's `Z_factor` (ifactor.ts), exactly as PARI's
 * elliptic-curve code calls `Z_factor`.  The previous local trial-division
 * copy had no primality short-circuit and made `ellorder` on a large prime
 * order take seconds.
 *
 * Reference: PARI ifactor1.c - Z_factor
 */
function factor(n: bigint): [bigint, bigint][] {
  if (n <= 1n) return [];
  return Z_factor(n).filter(([q]) => q > 0n) as [bigint, bigint][];
}

// ============================================================================
// FpX helpers (only what Fp_ellcard_Shanks needs: roots of a cubic)
// ============================================================================

/** Drop trailing zero coefficients (little-endian representation). */
function FpX_trim(a: bigint[]): bigint[] {
  let i = a.length;
  while (i > 0 && a[i - 1] === 0n) i--;
  return a.slice(0, i);
}

/** Polynomial product over Fp. */
function FpX_mul(a: bigint[], b: bigint[], p: bigint): bigint[] {
  if (a.length === 0 || b.length === 0) return [];
  const r = new Array<bigint>(a.length + b.length - 1).fill(0n);
  for (let i = 0; i < a.length; i++) {
    if (a[i] === 0n) continue;
    for (let j = 0; j < b.length; j++) {
      r[i + j] = mod(r[i + j]! + a[i]! * b[j]!, p);
    }
  }
  return FpX_trim(r);
}

/** Remainder of a modulo b over Fp (b nonzero). */
function FpX_rem(a: bigint[], b: bigint[], p: bigint): bigint[] {
  a = FpX_trim(a);
  b = FpX_trim(b);
  if (b.length === 0) throw new Error('FpX_rem: division by zero polynomial');
  if (a.length < b.length) return a;
  const inv = Fp_inv(b[b.length - 1]!, p);
  const r = [...a];
  for (let i = r.length - 1; i >= b.length - 1; i--) {
    const c = Fp_mul(r[i]!, inv, p);
    if (c === 0n) continue;
    const off = i - b.length + 1;
    for (let j = 0; j < b.length; j++) {
      r[off + j] = Fp_sub(r[off + j]!, Fp_mul(c, b[j]!, p), p);
    }
  }
  return FpX_trim(r);
}

/** GCD of two polynomials over Fp. */
function FpX_gcd(a: bigint[], b: bigint[], p: bigint): bigint[] {
  a = FpX_trim(a);
  b = FpX_trim(b);
  while (b.length !== 0) {
    const r = FpX_rem(a, b, p);
    a = b;
    b = r;
  }
  return a;
}

/** x^n mod T over Fp. */
function FpXQ_pow(x: bigint[], n: bigint, T: bigint[], p: bigint): bigint[] {
  let result: bigint[] = [1n];
  let base = FpX_rem(x, T, p);
  let e = n;
  while (e > 0n) {
    if ((e & 1n) === 1n) result = FpX_rem(FpX_mul(result, base, p), T, p);
    base = FpX_rem(FpX_mul(base, base, p), T, p);
    e >>= 1n;
  }
  return result;
}

/**
 * Number of distinct roots in Fp of the 2-division polynomial
 * x^3 + c4*x + c6 (i.e. the number of rational points of order 2, plus the
 * convention used by `Fp_ellcard_Shanks`).
 *
 * Computed as deg gcd(x^p - x, T), which for a squarefree cubic is 0, 1 or 3.
 *
 * Reference: PARI FpX.c - FpX_nbroots, used at FpE.c:938-943
 */
function FpX_nbroots_cubic(c4: bigint, c6: bigint, p: bigint): number {
  const T = FpX_trim([mod(c6, p), mod(c4, p), 0n, 1n]);
  if (T.length <= 1) return 0;
  if (p === 2n || p === 3n) {
    // tiny fields: direct enumeration is exact and cheaper
    let n = 0;
    for (let x = 0n; x < p; x++) {
      if (mod(x * x * x + c4 * x + c6, p) === 0n) n++;
    }
    return n;
  }
  // g = x^p - x (mod T)
  const g = FpXQ_pow([0n, 1n], p, T, p);
  while (g.length < 2) g.push(0n);
  g[1] = Fp_sub(g[1]!, 1n, p);
  const XpMinusX = FpX_trim(g);

  // T | x^p - x: every root of T is rational
  if (XpMinusX.length === 0) return T.length - 1;
  return FpX_gcd(T, XpMinusX, p).length - 1;
}

// ============================================================================
// Random Point Generation
// ============================================================================

/**
 * Generate a random field element in [0, p-1].
 */
function Fp_random(p: bigint): bigint {
  const bits = p.toString(2).length;
  const bytes = Math.ceil(bits / 8);

  for (;;) {
    let result = 0n;
    for (let i = 0; i < bytes; i++) {
      result = (result << 8n) | BigInt(Math.floor(Math.random() * 256));
    }
    result >>= BigInt(bytes * 8 - bits);
    if (result < p) return result;
  }
}

/**
 * Find a random point on the curve.
 *
 * Port of PARI FpE.c:369-385 - random_FpE:
 *
 * ```c
 *   do {
 *     x   = randomi(p);
 *     x2  = Fp_sqr(x, p);
 *     rhs = Fp_add(Fp_mul(x, Fp_add(x2, a4, p), p), a6, p);
 *   } while ((!signe(rhs) && !signe(Fp_add(Fp_mulu(x2,3,p),a4,p)))
 *           || kronecker(rhs, p) < 0);
 *   y = Fp_sqrt(rhs, p);
 * ```
 *
 * Like PARI we return the canonical square root `Fp_sqrt(rhs, p)`; the sign
 * is not randomized (`<P>` = `<-P>`, so no caller is affected).
 */
export function FpE_random(E: EllipticCurveFp): EllipticPointFp {
  const { a4, a6, p } = E;

  for (;;) {
    const x = Fp_random(p);
    const x2 = Fp_sqr(x, p);
    const rhs = Fp_add(Fp_mul(x, Fp_add(x2, a4, p), p), a6, p);

    // singular point of a singular curve
    if (rhs === 0n && Fp_add(Fp_mulu(x2, 3, p), a4, p) === 0n) continue;
    if (kronecker(rhs, p) < 0) continue;

    const y = Fp_sqrt(rhs, p);
    if (y === null) throw new Error(`FpE_random: ${p} is not prime`);
    return ellpoint(x, y);
  }
}

// ============================================================================
// Cardinality Computation
// ============================================================================

/**
 * Compute the trace of an elliptic curve over Fp using naive enumeration.
 *
 * The trace a_p satisfies: #E(Fp) = p + 1 - a_p
 *
 * Reference: PARI FpE.c:811-831 - Fl_elltrace_naive
 *
 * @param a4 - Coefficient a4 in y^2 = x^3 + a4*x + a6
 * @param a6 - Coefficient a6
 * @param p - Prime field characteristic
 * @returns The trace a_p
 */
export function Fp_elltrace_naive(a4: bigint, a6: bigint, p: bigint): bigint {
  let trace = 0n;

  for (let x = 0n; x < p; x++) {
    // y^2 = x^3 + a4*x + a6
    const ySquared = mod(Fp_pow(x, 3n, p) + a4 * x + a6, p);

    // Add the Legendre symbol to trace
    trace -= BigInt(kronecker(ySquared, p));
  }

  return trace;
}

/**
 * Compute elliptic curve cardinality using exhaustive enumeration.
 *
 * Reference: PARI FpE.c:811-831
 */
function ellcard_exhaustive(E: EllipticCurveFp): bigint {
  const { a4, a6, p } = E;
  return p + 1n - Fp_elltrace_naive(a4, a6, p);
}

/**
 * Find the lift of a (mod b) which is closest to c.
 *
 * `x = round((c-a)/b) = floor((2(c-a) + b) / 2b)`; return `a + b*x`.
 *
 * Reference: PARI FpE.c:862-867 - closest_lift
 */
function closest_lift(a: bigint, b: bigint, c: bigint): bigint {
  return a + b * floorDiv(2n * (c - a) + b, 2n * b);
}

/**
 * Table size for the baby-step/giant-step search.
 *
 * PARI: `t = ceil(sqrt(pordmin / B)); return t >> 1;` -- computed exactly
 * here as the smallest t with t^2 * B >= pordmin.
 *
 * Reference: PARI FpE.c:869-878 - get_table_size
 */
function get_table_size(pordmin: bigint, B: bigint): bigint {
  let t = isqrt(pordmin / B);
  while (t * t * B < pordmin) t++;
  while (t > 0n && (t - 1n) * (t - 1n) * B >= pordmin) t--;
  /* `if (is_bigint(t)) pari_err_OVERFLOW(...)` */
  if (t >= 1n << 63n) {
    throw new Error("ellap [large prime: install the 'seadata' package]");
  }
  return t >> 1n;
}

/**
 * Find x such that kronecker(u = x^3 + c4*x + c6, p) is KRO.
 * Return the point [x*u, u^2] on E_u (KRO = 1) / E^twist (KRO = -1).
 *
 * Reference: PARI FpE.c:880-894 - Fp_ellpoint
 */
function Fp_ellpoint(
  KRO: number,
  startX: bigint,
  c4: bigint,
  c6: bigint,
  p: bigint
): [EllipticPointFp, bigint] {
  let x = startX;

  for (;;) {
    x++; /* u = x^3 + c4 x + c6 */
    const u = mod(c6 + x * (c4 + x * x), p);
    if (kronecker(u, p) === KRO) {
      return [ellpoint(mod(x * u, p), Fp_sqr(u, p)), x];
    }
  }
}

/**
 * Solve x^2 + y^2 = p for a prime p = 1 (mod 4).
 *
 * This is the `cornacchia2(4, p)` instance used by `ap_j1728`
 * (x^2 + 4y^2 = 4p with x = 2u reduces to u^2 + y^2 = p).
 *
 * Reference: PARI Qfb.c:2028-2077 - cornacchia2
 */
function sum_of_two_squares(p: bigint): [bigint, bigint] {
  const i = Fp_sqrt(Fp_neg(1n, p), p);
  if (i === null) throw new Error(`sum_of_two_squares: ${p} is not 1 mod 4`);

  let a = p;
  let b = i;
  while (b * b > p) {
    const t = a % b;
    a = b;
    b = t;
  }
  const s2 = p - b * b;
  const s = isqrt(s2);
  if (s * s !== s2) throw new Error(`sum_of_two_squares failed for ${p}`);
  return [b, s];
}

/**
 * Trace of Frobenius for the CM curve y^2 = x^3 + a4*x (j = 1728).
 *
 * ```c
 * static GEN ap_j1728(GEN a4, GEN p)
 * {
 *   if (mod4(p) != 1) return gen_0;
 *   (void)cornacchia2(utoipos(4), p, &a, &b);   // a^2 + 4b^2 = 4p
 *   if (Mod4(a)==0) a = b;
 *   if (Mod2(a)==1) a = shifti(a,1);
 *   if (Mod8(a)==6) a = negi(a);
 *   e = shifti(p,-2);                            // (p-1)/4
 *   return centermod(mulii(a, Fp_pow(a4, e, p)), p);
 * }
 * ```
 *
 * The three normalisation steps amount to: `a = 2c` where c is the odd
 * member of the pair (u, v) with u^2 + v^2 = p, signed so that c = 1 (mod 4).
 *
 * Reference: PARI FpE.c:1293-1305 - ap_j1728
 */
function ap_j1728(a4: bigint, p: bigint): bigint {
  if (mod(p, 4n) !== 1n) return 0n;

  const [u, v] = sum_of_two_squares(p);
  let c = (u & 1n) === 1n ? u : v; /* Mod4(a)==0 -> a = b; Mod2(a)==1 -> a = 2a */
  if (mod(c, 4n) !== 1n) c = -c; /* Mod8(a)==6 -> a = -a */
  const a = 2n * c;

  const e = (p - 1n) / 4n;
  return Fp_center(mod(a * Fp_pow(mod(a4, p), e, p), p), p);
}

/**
 * Exact order of the point z on E, knowing that [o]z = O.
 *
 * Reference: PARI FpE.c:405-423 - FpE_order (-> bb_group.c gen_order)
 */
function FpE_order(z: EllipticPointFp, o: bigint, a4: bigint, p: bigint): bigint {
  if (o <= 0n) throw new Error(`FpE_order: invalid bound ${o}`);
  let order = o;
  for (const [q] of factor(o)) {
    while (order % q === 0n) {
      const t = order / q;
      if (ell_is_inf(FpE_mul(z, t, a4, p))) order = t;
      else break;
    }
  }
  return order;
}

/**
 * Multiplicative order of a in Fp*, knowing that a^N = 1.
 *
 * Reference: PARI arith1.c - Fp_order
 */
function Fp_order(a: bigint, N: bigint, p: bigint): bigint {
  let order = N;
  for (const [q] of factor(N)) {
    while (order % q === 0n) {
      const t = order / q;
      if (Fp_pow(a, t, p) === 1n) order = t;
      else break;
    }
  }
  return order;
}

/**
 * Baby-step/giant-step search for a multiple of ord(f) of the form h + k*B.
 *
 * Mirrors the `s < 3` naive branch and the BSGS branch of
 * `Fp_ellcard_Shanks` (FpE.c:975-1120).  PARI hashes the low word of the
 * coordinates and re-verifies; since we compare exact bigints no
 * verification step is needed.
 *
 * @returns a nonzero integer h' with [h']f = O
 */
function ellcard_bsgs_search(
  f: EllipticPointFp,
  fh: EllipticPointFp,
  h0: bigint,
  B: bigint,
  s: bigint,
  a4: bigint,
  p: bigint
): bigint {
  const F = FpE_mul(f, B, a4, p);

  /* h + k*B annihilates f; 0 carries no information, and a negative value is
   * equivalent to its absolute value. */
  const accept = (v: bigint): bigint | null => (v === 0n ? null : v < 0n ? -v : v);

  if (s < 3n) {
    /* we're nearly done: naive search (FpE.c:983-994) */
    let P = fh;
    let q1 = fh;
    const mF = FpE_neg(F, p);
    for (let i = 1n; ; i++) {
      P = FpE_add(P, F, a4, p);
      if (ell_is_inf(P)) {
        const v = accept(h0 + i * B);
        if (v !== null) return v;
      }
      q1 = FpE_add(q1, mF, a4, p);
      if (ell_is_inf(q1)) {
        const v = accept(h0 - i * B);
        if (v !== null) return v;
      }
    }
  }

  /* baby steps: table of h.f + j.F for j = 0 .. s-1 */
  const tbl = new Map<bigint, [bigint, bigint][]>();
  let P = fh;
  for (let j = 0n; j < s; j++) {
    const key = P.x!;
    const bucket = tbl.get(key);
    if (bucket) bucket.push([P.y!, j]);
    else tbl.set(key, [[P.y!, j]]);
    P = FpE_add(P, F, a4, p);
    if (ell_is_inf(P)) {
      const v = accept(h0 + (j + 1n) * B);
      if (v !== null) return v;
    }
  }

  /* fg = s.F = (h.f + s.F) - h.f */
  const sF = FpE_add(P, FpE_neg(fh, p), a4, p);
  if (ell_is_inf(sF)) return s * B;

  /* giant steps: ftest = (s*i).F */
  let ftest = sF;
  for (let i = 1n; ; i++) {
    if (ell_is_inf(ftest)) return s * i * B;
    const bucket = tbl.get(ftest.x!);
    if (bucket) {
      const negY = Fp_neg(ftest.y!, p);
      for (const [y, j2] of bucket) {
        /* [h + j2*B] f == +/- [s*i*B] f */
        let v: bigint | null = null;
        if (y === ftest.y!) v = accept(h0 + (j2 - s * i) * B);
        else if (y === negY) v = accept(h0 + (j2 + s * i) * B);
        if (v !== null) return v;
      }
    }
    ftest = FpE_add(ftest, sF, a4, p);
  }
}

/**
 * Compute #E(Fp) with Shanks/Mestre baby-step giant-step.
 *
 * Faithful port of PARI FpE.c:920-1135 - Fp_ellcard_Shanks:
 * - `pordmin = ceil(4 sqrt(p))`, so that #E is determined once it is known
 *   modulo B >= pordmin;
 * - (A, B) is seeded from the number of rational 2-torsion points;
 * - each round works on E_u for a point u of Legendre symbol +/-1, i.e. it
 *   alternates between the curve and its quadratic twist, using
 *   `#E(Fp) + #E'(Fp) = 2p + 2` to transport the congruence;
 * - the loop only exits once B >= pordmin (the answer is then unique in the
 *   Hasse interval);
 * - c6 = 0 (j = 1728) is handled by `ap_j1728` since the seed point
 *   [0, c6^2] degenerates there.
 *
 * Assumes p > 457.  Exported (PARI keeps it `static`) so that the test suite
 * can exercise the BSGS branch for primes below `ellcard`'s 2048 cutoff.
 */
export function Fp_ellcard_Shanks(c4: bigint, c6: bigint, p: bigint): bigint {
  if (c6 === 0n) {
    return p + 1n - ap_j1728(c4, p);
  }

  /* once #E(Fp) is known mod B >= pordmin, it is completely determined */
  const pordmin = isqrt(16n * p) + 1n; /* ceil( 4 sqrt(p) ) */
  const p1p = p + 1n;
  const p2p = 2n * p1p;

  let x = 0n;
  let KRO = 0;
  let A: bigint;
  let B: bigint;

  /* how many 2-torsion points ? */
  switch (FpX_nbroots_cubic(c4, c6, p)) {
    case 3:
      A = 0n;
      B = 4n;
      break;
    case 1:
      A = 0n;
      B = 2n;
      break;
    default:
      A = 1n;
      B = 2n;
      break; /* 0 */
  }

  for (;;) {
    let h = closest_lift(A, B, p1p);

    let f: EllipticPointFp;
    if (KRO === 0) {
      /* first time, initialize */
      KRO = kronecker(c6, p);
      f = ellpoint(0n, Fp_sqr(c6, p));
    } else {
      KRO = -KRO;
      [f, x] = Fp_ellpoint(KRO, x, c4, c6, p);
    }

    /* [ux, u^2] is on E_u: y^2 = x^3 + c4 u^2 x + c6 u^3
     * E_u isomorphic to E (resp. E') iff KRO = 1 (resp. -1)
     * #E(F_p) = p+1 - a_p, #E'(F_p) = p+1 + a_p */
    const a4 = Fp_mul(c4, f.y!, p); /* c4 for E_u */
    const fh = FpE_mul(f, h, a4, p);

    if (!ell_is_inf(fh)) {
      const s = get_table_size(pordmin, B);
      h = ellcard_bsgs_search(f, fh, h, B, s, a4, p);
    }

    /* found a point of exponent h on E_u */
    const ord = FpE_order(f, h, a4, p);
    /* ord | #E_u(Fp) = A (mod B) */
    [A, B] = Z_chinese_all(A, 0n, B, ord);
    if (B >= pordmin) break;

    /* not done: update A mod B for the _next_ curve, isomorphic to
     * the quadratic twist of this one; #E(Fp)+#E'(Fp) = 2p+2 */
    A = mod(p2p - A, B);
  }

  const h = closest_lift(A, B, p1p);
  return KRO === 1 ? h : p2p - h;
}

/**
 * Compute the cardinality (number of points) of E(Fp).
 *
 * PARI dispatches on the size of p (FpE.c:1424-1437 - Fp_ellcard):
 * - expi(p) < 11 (p < 2048): naive trace enumeration;
 * - CM curves: closed formulas (we implement the j = 1728 case, which
 *   `Fp_ellcard_Shanks` needs anyway);
 * - expi(p) >= 56: SEA (not ported -- we fall through to Shanks);
 * - otherwise: Shanks/Mestre baby-step giant-step.
 *
 * @param E - The elliptic curve
 * @returns The number of points on E(Fp)
 *
 * @see Deviation: PARI Elliptic Curve Advanced Algorithms Missing (parigp-ts)
 */
export function ellcard(E: EllipticCurveFp): bigint {
  if (E._card !== undefined) {
    return E._card;
  }

  const { a4, a6, p } = E;

  /* expi(p) < 11 <=> p < 2^11 */
  const card = p < 2048n ? ellcard_exhaustive(E) : Fp_ellcard_Shanks(mod(a4, p), mod(a6, p), p);

  E._card = card;
  return card;
}

// ============================================================================
// Point Order
// ============================================================================

/**
 * Compute the order of a point P on E(Fp).
 *
 * The order is the smallest positive integer n such that [n]P = O.
 * Uses factorization of the curve order.
 *
 * Reference: PARI FpE.c:405-423 - FpE_order
 *
 * @param E - The elliptic curve
 * @param P - A point on E
 * @param curveOrder - Optional known curve order (will be computed if not provided)
 * @returns The order of P
 */
export function ellorder(E: EllipticCurveFp, P: EllipticPointFp, curveOrder?: bigint): bigint {
  if (P.isInfinity) {
    return 1n;
  }

  const N = curveOrder ?? ellcard(E);
  return FpE_order(P, N, E.a4, E.p);
}

// ============================================================================
// Group Structure
// ============================================================================

/**
 * Order of the Weil pairing of P and Q, both of order dividing m.
 *
 * ```c
 * static GEN _FpE_pairorder(void *E, GEN P, GEN Q, GEN m, GEN F)
 * { return Fp_order(FpE_weilpairing(P,Q,m,e->a4,e->p), F, e->p); }
 * ```
 *
 * Reference: PARI FpE.c:1462-1467 - _FpE_pairorder
 */
function pairorder(
  E: EllipticCurveFp,
  P: EllipticPointFp,
  Q: EllipticPointFp,
  m: bigint,
  F: bigint
): bigint {
  const w = ellweilpairing(E, P, Q, m);
  return Fp_order(w, F, E.p);
}

/**
 * `c = prod_{q^2 | (N, d^2)} q^{v_q(N)}` together with its factorization;
 * `c` is a multiple of d2.  Returns null when the group must be cyclic.
 *
 * ```c
 * static GEN d2_multiple(GEN N, GEN d)
 * {
 *   GEN Q = gel(Z_factor(gcdii(N,d)), 1);
 *   for (i = 1, j = 1; i < l; i++) {
 *     long v = Z_pval(N, gel(Q,i));
 *     if (v <= 1) continue;
 *     gel(P,j) = gel(Q,i); gel(E,j) = utoipos(v); j++;
 *   }
 *   if (j == 1) return NULL;
 *   return mkvec2(factorback2(P,E), mkmat2(P,E));
 * }
 * ```
 *
 * Note that the exponent kept is `v_q(N)`, **not** `v_q(gcd(N,d))`, and that
 * primes with `v_q(N) <= 1` are dropped entirely.
 *
 * Reference: PARI bb_group.c:967-986 - d2_multiple
 */
function d2_multiple(N: bigint, d: bigint): { N0: bigint; fa: [bigint, bigint][] } | null {
  const Q = factor(gcd(N, d));
  const fa: [bigint, bigint][] = [];
  let N0 = 1n;
  for (const [q] of Q) {
    const v = valuation(N, q);
    if (v <= 1n) continue;
    fa.push([q, v]);
    N0 *= q ** v;
  }
  if (fa.length === 0) return null;
  return { N0, fa };
}

/**
 * Elementary divisors [d1, d2] (d2 | d1) of a group of order N whose
 * "second" invariant divides d, together with the exponent `m` handed to
 * `gen_ellgens`.
 *
 * Reference: PARI bb_group.c:988-1046 - gen_ellgroup
 *
 * @see Deviation: `m` is `g1` (i.e. `dm = N1`), not the last iteration's
 * `lcm(ord P, ord Q)`.  The vendored PARI 2.18-dev writes `*pm = g1` and then
 * immediately overwrites it with `*pm = m`; when the primes of N0 are not all
 * settled in a single iteration, that `m` need not be a multiple of d2, and
 * `gen_ellgens`'s `do ... while (d != d2)` can then never terminate (the
 * pairing of two m-torsion points has order dividing m).  Concrete case:
 * E/F_43: y^2 = x^3+7x+8, group [12, 2^2*3] -- roughly 0.5 % of runs produce
 * m = 4 with d2 = 3.  `g1` always satisfies d2 | g1 | d1 (the l-parts are
 * l^a and l^b with a >= b), and PARI 2.15.4 (Sage 10.3) never hangs on that
 * curve over 4000 runs, so `g1` is what the shipping behaviour requires.
 */
function gen_ellgroup(E: EllipticCurveFp, N: bigint, d: bigint): { D: bigint[]; m: bigint } {
  const { a4, p } = E;

  if (N === 1n) return { D: [], m: 1n };

  const F = d2_multiple(N, d);
  if (F === null) return { D: [N], m: 1n };

  const N0 = F.N0; /* a multiple of d2 */
  const N1 = N / N0; /* N1 | d1 */
  const L0 = F.fa.map(([q]) => q); /* primes dividing N0 */
  const E0 = F.fa.map(([, e]) => e); /* ... and their exponents */
  const n0 = L0.length;

  let g1 = 1n;
  let g2 = 1n;
  let n = 0;

  for (let guard = 0; guard < 10000; guard++) {
    /* g1 | (d1/N1), g2 | d2 */
    const P = FpE_mul(FpE_random(E), N1, a4, p);
    const s = ell_is_inf(P) ? 1n : FpE_order(P, N0, a4, p); /* s | N0 */
    if (s === N0) return { D: [N], m: 1n };

    const Q = FpE_mul(FpE_random(E), N1, a4, p);
    const t = ell_is_inf(Q) ? 1n : FpE_order(Q, N0, a4, p); /* t | N0 */
    if (t === N0) return { D: [N], m: 1n };

    const m = lcm(s, t); /* m | N0 */
    const mo = m * pairorder(E, P, Q, m, N0);

    /* For each prime l dividing N0, check whether P and Q
     * generate all rational points of order a power of l */
    for (let j = 0; j < n0; j++) {
      const e = E0[j]!;
      if (e === 0n) continue;
      const l = L0[j]!;
      if (valuation(mo, l) === e) {
        const vm = valuation(m, l);
        g1 *= l ** vm;
        g2 *= l ** (e - vm);
        if (++n === n0) {
          /* done with all primes l */
          if (g2 === 1n) return { D: [N], m: 1n };
          return { D: [g1 * N1, g2], m: g1 };
        }
        E0[j] = 0n; /* done with this prime l */
      }
    }
  }

  throw new Error('gen_ellgroup: failed to determine the group structure');
}

/**
 * Determine the group structure of E(Fp).
 *
 * Returns [] for the trivial group, [d1] for cyclic groups Z/d1Z,
 * or [d1, d2] for non-cyclic groups Z/d1Z x Z/d2Z with d2 | d1.
 *
 * Reference: PARI FpE.c:1469-1476 - Fp_ellgroup
 *   `gen_ellgroup(N, p-1, &m, E, &FpE_group, _FpE_pairorder)`
 *
 * @param E - The elliptic curve
 * @returns Group structure as [], [d1] or [d1, d2]
 */
export function ellgroup(E: EllipticCurveFp): bigint[] {
  if (E._group !== undefined) {
    return E._group;
  }

  const N = ellcard(E);
  const { D, m } = gen_ellgroup(E, N, E.p - 1n);

  E._group = D;
  E._m = m;
  return D;
}

// ============================================================================
// Generators
// ============================================================================

/**
 * Find generators for the group E(Fp).
 *
 * Reference: PARI FpE.c:1478-1497 - Fp_ellgens
 *   - cyclic: `gen_gener(d1, ...)`
 *   - otherwise: `gen_ellgens(d1, d2, m, ...)` (bb_group.c:1048-1072)
 *
 * ```c
 * do { P = grp->rand(E); s = gen_order(P, F, E, grp); } while (!equalii(s,d1));
 * do { Q = grp->rand(E);
 *      d = pairorder(E, grp->pow(E,P,dm), grp->pow(E,Q,dm), m, F);
 * } while (!equalii(d, d2));
 * ```
 *
 * @param E - The elliptic curve
 * @returns Array of generator points
 */
export function ellgenerators(E: EllipticCurveFp): EllipticPointFp[] {
  if (E._generators !== undefined) {
    return E._generators;
  }

  const { a4, p } = E;
  const N = ellcard(E);
  const D = ellgroup(E);

  if (D.length === 0) {
    E._generators = [];
    return [];
  }

  const d1 = D[0]!;

  /* gen_gener: a random point of order exactly d1 generates */
  let P: EllipticPointFp | null = null;
  for (let attempts = 0; attempts < 10000; attempts++) {
    const R = FpE_random(E);
    if (FpE_order(R, N, a4, p) === d1) {
      P = R;
      break;
    }
  }
  if (P === null) throw new Error('ellgenerators: failed to find a generator');

  if (D.length === 1) {
    E._generators = [P];
    return [P];
  }

  /* gen_ellgens: pick Q so that the Weil pairing of [dm]P and [dm]Q has
   * order exactly d2.  This is PARI's independence test -- enumerating the
   * multiples of a point would be O(d2) = O(sqrt p). */
  const d2 = D[1]!;
  const m = E._m ?? d1;
  const dm = d1 / m;

  for (let attempts = 0; attempts < 10000; attempts++) {
    const Q = FpE_random(E);
    const d = pairorder(E, FpE_mul(P, dm, a4, p), FpE_mul(Q, dm, a4, p), m, d1);
    if (d === d2) {
      E._generators = [P, Q];
      return [P, Q];
    }
  }

  throw new Error('ellgenerators: failed to find a second generator');
}

// ============================================================================
// Trace of Frobenius
// ============================================================================

/**
 * Compute the trace of Frobenius for E/Fp.
 *
 * The trace a_p satisfies: #E(Fp) = p + 1 - a_p
 *
 * By the Hasse bound: |a_p| <= 2*sqrt(p)
 *
 * Reference: PARI elliptic.c:5241-5250
 *
 * @param E - The elliptic curve
 * @returns The trace of Frobenius
 */
export function trace_of_frobenius(E: EllipticCurveFp): bigint {
  const { p } = E;
  const card = ellcard(E);
  return p + 1n - card;
}

// ============================================================================
// Convenience Functions
// ============================================================================

/**
 * Create an elliptic curve in short Weierstrass form.
 *
 * @param a4 - Coefficient a4 in y^2 = x^3 + a4*x + a6
 * @param a6 - Coefficient a6
 * @param p - Prime field characteristic
 */
export function ellinit_Fp(a4: bigint, a6: bigint, p: bigint): EllipticCurveFp {
  // Reduce coefficients mod p
  a4 = mod(a4, p);
  a6 = mod(a6, p);

  // Check discriminant: -4a4^3 - 27a6^2 != 0
  const disc = mod(-4n * Fp_pow(a4, 3n, p) - 27n * Fp_pow(a6, 2n, p), p);
  if (disc === 0n) {
    throw new Error('Singular curve: discriminant is zero');
  }

  return { a4, a6, p };
}

/**
 * Check if a point is on the curve.
 *
 * @param E - The elliptic curve
 * @param P - The point to check
 */
export function ellisoncurve(E: EllipticCurveFp, P: EllipticPointFp): boolean {
  if (P.isInfinity) return true;

  const { a4, a6, p } = E;
  const x = P.x!;
  const y = P.y!;

  // Check: y^2 = x^3 + a4*x + a6
  const lhs = Fp_sqr(y, p);
  const rhs = mod(Fp_pow(x, 3n, p) + a4 * x + a6, p);

  return lhs === rhs;
}

/**
 * Lift an x-coordinate to a point on the curve.
 *
 * @param E - The elliptic curve
 * @param x - The x-coordinate
 * @returns A point (x, y) on E, or null if no such point exists
 */
export function elllift_x(E: EllipticCurveFp, x: bigint): EllipticPointFp | null {
  const { a4, a6, p } = E;
  x = mod(x, p);

  // y^2 = x^3 + a4*x + a6
  const ySquared = mod(Fp_pow(x, 3n, p) + a4 * x + a6, p);
  const y = Fp_sqrt(ySquared, p);

  if (y === null) return null;
  return ellpoint(x, y);
}
