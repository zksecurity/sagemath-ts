/**
 * Elliptic Curve Point Operations - Jacobian Coordinates
 *
 * TypeScript port of PARI/GP elliptic curve point arithmetic using
 * Jacobian coordinates for efficient scalar multiplication.
 *
 * Reference: reference/pari/src/basemath/FpE.c
 *
 * This module implements Jacobian coordinate operations for elliptic curves
 * in short Weierstrass form: y^2 = x^3 + a4*x + a6 over Fp.
 *
 * PARI/GP is free software under the GNU GPL v2+.
 */

import { type EllipticCurve, EllCurveType, EllipticCurveError, ell_to_a4a6_bc, ellinit } from './init.js';
import { FpE_mul } from './group.js';
import {
  type FFEllipticCurve, type FFEllipticInputPoint, type FFEllipticPoint, FF_ellmul,
} from '../_elliptic_finite_field.js';

import { Fp_inv, Fp_add, Fp_double, Fp_mul, Fp_mulu, Fp_neg, Fp_sqr, Fp_sub } from '../ff.js';

import {
  type EllipticPoint,
  FpE_isoncurve,
  FpE_to_FpJ,
  FpJ_is_inf,
  FpJ_to_FpE,
  type EllipticPointFinite,
  type JacobianPoint,
  type ShortWeierstrassCurve,
  ell_is_inf,
  ellinf,
  ellinf_FpJ,
  mkpoint,
} from './points.js';

// Re-export types and utilities from points.ts for convenience
export type { EllipticPoint, EllipticPointFinite, JacobianPoint, ShortWeierstrassCurve };

export { ellinf, ell_is_inf, ellinf_FpJ, FpJ_is_inf, FpE_to_FpJ, FpJ_to_FpE, mkpoint };

const scalarModels = new WeakMap<
  ShortWeierstrassCurve,
  {
    a4: bigint;
    a6: bigint;
    p: bigint;
    model: EllipticCurve | FFEllipticCurve;
  }
>();

/** Adapt the short record (implicit a1=a2=a3=0) to native ellinit invariants. */
function scalarModel(E: ShortWeierstrassCurve): EllipticCurve | FFEllipticCurve {
  const cached = scalarModels.get(E);
  if (cached && cached.a4 === E.a4 && cached.a6 === E.a6 && cached.p === E.p) return cached.model;
  const model = ellinit([E.a4, E.a6], E.p);
  if (model === null) throw new EllipticCurveError('Curve is singular (discriminant is zero mod p)');
  scalarModels.set(E, { a4: E.a4, a6: E.a6, p: E.p, model });
  return model;
}

// =============================================================================
// Jacobian Coordinate Operations
// =============================================================================
// Source: FpE.c:26-105

/**
 * Negate a Jacobian point
 * Source: FpE.c:108-111, FpJ_neg
 */
export function FpJ_neg(P: JacobianPoint, p: bigint): JacobianPoint {
  return { X: P.X, Y: Fp_neg(P.Y, p), Z: P.Z };
}

/**
 * Jacobian point doubling using dbl-2007-bl formula
 * Cost: 1M + 8S + 1*a + 10add + 1*8 + 2*2 + 1*3
 * Source: FpE.c:39-61, FpJ_dbl
 * Reference: http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-dbl-2007-bl
 */
export function FpJ_dbl(P: JacobianPoint, a4: bigint, p: bigint): JacobianPoint {
  // Point at infinity doubles to itself
  if (P.Z === 0n) {
    return ellinf_FpJ();
  }

  const X1 = P.X;
  const Y1 = P.Y;
  const Z1 = P.Z;

  // XX = X1^2
  const XX = Fp_sqr(X1, p);
  // YY = Y1^2
  const YY = Fp_sqr(Y1, p);
  // YYYY = YY^2
  const YYYY = Fp_sqr(YY, p);
  // ZZ = Z1^2
  const ZZ = Fp_sqr(Z1, p);

  // S = 2*((X1+YY)^2 - XX - YYYY)
  const S = Fp_double(Fp_sub(Fp_sqr(Fp_add(X1, YY, p), p), Fp_add(XX, YYYY, p), p), p);

  // M = 3*XX + a4*ZZ^2
  const M = Fp_add(Fp_mulu(XX, 3, p), Fp_mul(a4, Fp_sqr(ZZ, p), p), p);

  // T = M^2 - 2*S
  const T = Fp_sub(Fp_sqr(M, p), Fp_double(S, p), p);

  // X3 = T
  const X3 = T;

  // Y3 = M*(S-T) - 8*YYYY
  const Y3 = Fp_sub(Fp_mul(M, Fp_sub(S, T, p), p), Fp_mulu(YYYY, 8, p), p);

  // Z3 = (Y1+Z1)^2 - YY - ZZ
  const Z3 = Fp_sub(Fp_sqr(Fp_add(Y1, Z1, p), p), Fp_add(YY, ZZ, p), p);

  return { X: X3, Y: Y3, Z: Z3 };
}

/**
 * Jacobian point addition using add-2007-bl formula
 * Cost: 11M + 5S + 9add + 4*2
 * Source: FpE.c:65-105, FpJ_add
 * Reference: http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-add-2007-bl
 */
export function FpJ_add(P: JacobianPoint, Q: JacobianPoint, a4: bigint, p: bigint): JacobianPoint {
  // Handle point at infinity
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

  // Z1Z1 = Z1^2
  const Z1Z1 = Fp_sqr(Z1, p);
  // Z2Z2 = Z2^2
  const Z2Z2 = Fp_sqr(Z2, p);
  // U1 = X1*Z2Z2
  const U1 = Fp_mul(X1, Z2Z2, p);
  // U2 = X2*Z1Z1
  const U2 = Fp_mul(X2, Z1Z1, p);
  // S1 = Y1*Z2*Z2Z2
  const S1 = Fp_mul(Y1, Fp_mul(Z2, Z2Z2, p), p);
  // S2 = Y2*Z1*Z1Z1
  const S2 = Fp_mul(Y2, Fp_mul(Z1, Z1Z1, p), p);
  // H = U2 - U1
  const H = Fp_sub(U2, U1, p);
  // r = 2*(S2 - S1)
  const r = Fp_double(Fp_sub(S2, S1, p), p);

  // If points are equal we must double
  if (H === 0n) {
    if (r === 0n) {
      // Points are equal, so double
      return FpJ_dbl(P, a4, p);
    } else {
      // Points are inverses, return infinity
      return ellinf_FpJ();
    }
  }

  // I = (2*H)^2
  const I = Fp_sqr(Fp_double(H, p), p);
  // J = H*I
  const J = Fp_mul(H, I, p);
  // V = U1*I
  const V = Fp_mul(U1, I, p);

  // X3 = r^2 - J - 2*V
  const X3 = Fp_sub(Fp_sqr(r, p), Fp_add(J, Fp_double(V, p), p), p);

  // Y3 = r*(V - X3) - 2*S1*J
  const Y3 = Fp_sub(Fp_mul(r, Fp_sub(V, X3, p), p), Fp_double(Fp_mul(S1, J, p), p), p);

  // Z3 = ((Z1+Z2)^2 - Z1Z1 - Z2Z2)*H
  const Z3 = Fp_mul(Fp_sub(Fp_sqr(Fp_add(Z1, Z2, p), p), Fp_add(Z1Z1, Z2Z2, p), p), H, p);

  return { X: X3, Y: Y3, Z: Z3 };
}

// =============================================================================
// Affine Point Operations using points.ts types
// =============================================================================

/**
 * Check if a point lies on the curve y^2 = x^3 + a4*x + a6 (mod p)
 * Source: elliptic.c:2003-2037, ellisoncurve
 *
 * This is a wrapper around FpE_isoncurve from points.ts.
 */
export function ellisoncurve(E: ShortWeierstrassCurve, P: EllipticPoint): boolean {
  return FpE_isoncurve(E, P);
}

/**
 * Point negation: -P = (x, -y) for short Weierstrass form
 * Source: elliptic.c:2113-2126, ellneg
 * Source: FpE.c:316-320, FpE_neg
 */
export function ellneg(E: ShortWeierstrassCurve, P: EllipticPoint): EllipticPoint {
  if (ell_is_inf(P)) {
    return ellinf();
  }

  const { x, y } = P;
  const { p } = E;

  return mkpoint(x, Fp_neg(y, p));
}

/**
 * Affine point doubling for short Weierstrass curves
 * Source: FpE.c:257-268, FpE_dbl_slope
 */
function FpE_dbl(P: EllipticPoint, a4: bigint, p: bigint): EllipticPoint {
  if (ell_is_inf(P)) {
    return ellinf();
  }

  const { x, y } = P;

  // If y = 0, the tangent is vertical and 2P = O
  if (y === 0n) {
    return ellinf();
  }

  // slope = (3*x^2 + a4) / (2*y)
  const num = Fp_add(Fp_mulu(Fp_sqr(x, p), 3, p), a4, p);
  const den = Fp_mulu(y, 2, p);
  const slope = Fp_mul(num, modInv(den, p), p);

  // x' = slope^2 - 2*x
  const xr = Fp_sub(Fp_sqr(slope, p), Fp_mulu(x, 2, p), p);

  // y' = slope*(x - x') - y
  const yr = Fp_sub(Fp_mul(slope, Fp_sub(x, xr, p), p), y, p);

  return mkpoint(xr, yr);
}

/**
 * Modular inverse helper
 */
function modInv(a: bigint, p: bigint): bigint {
  let [oldR, r] = [a % p, p];
  let [oldS, s] = [1n, 0n];

  while (r !== 0n) {
    const quotient = oldR / r;
    [oldR, r] = [r, oldR - quotient * r];
    [oldS, s] = [s, oldS - quotient * s];
  }

  if (oldR !== 1n && oldR !== -1n) {
    throw new Error(`No inverse: gcd(${a}, ${p}) = ${oldR}`);
  }

  let result = oldS % p;
  if (result < 0n) result += p;
  return result;
}

/**
 * Affine point addition for short Weierstrass curves
 * Source: FpE.c:279-306, FpE_add_slope
 */
function FpE_add(P: EllipticPoint, Q: EllipticPoint, a4: bigint, p: bigint): EllipticPoint {
  if (ell_is_inf(P)) {
    if (ell_is_inf(Q)) {
      return ellinf();
    }
    return mkpoint(Q.x, Q.y);
  }
  if (ell_is_inf(Q)) {
    return mkpoint(P.x, P.y);
  }

  const Px = P.x;
  const Py = P.y;
  const Qx = Q.x;
  const Qy = Q.y;

  // Check if x-coordinates are equal
  if (Px === Qx) {
    // If y-coordinates are equal, use doubling formula
    if (Py === Qy) {
      return FpE_dbl(P, a4, p);
    }
    // Otherwise points are inverses, return infinity
    return ellinf();
  }

  // slope = (Qy - Py) / (Qx - Px)
  const slope = Fp_mul(Fp_sub(Py, Qy, p), modInv(Fp_sub(Px, Qx, p), p), p);

  // x' = slope^2 - Px - Qx
  const xr = Fp_sub(Fp_sub(Fp_sqr(slope, p), Px, p), Qx, p);

  // y' = slope*(Px - x') - Py
  const yr = Fp_sub(Fp_mul(slope, Fp_sub(Px, xr, p), p), Py, p);

  return mkpoint(xr, yr);
}

/**
 * Point addition: P + Q on curve E
 * Source: elliptic.c:2061-2094, elladd
 * Source: FpE.c:301-306, FpE_add
 */
export function elladd(
  E: ShortWeierstrassCurve,
  P: EllipticPoint,
  Q: EllipticPoint
): EllipticPoint {
  return FpE_add(P, Q, E.a4, E.p);
}

/**
 * Point subtraction: P - Q = P + (-Q)
 * Source: elliptic.c:2129-2136, ellsub
 * Source: FpE.c:323-328, FpE_sub
 */
export function ellsub(
  E: ShortWeierstrassCurve,
  P: EllipticPoint,
  Q: EllipticPoint
): EllipticPoint {
  return elladd(E, P, ellneg(E, Q));
}

// =============================================================================
// Scalar Multiplication
// =============================================================================

/**
 * Scalar multiplication: [n]P
 * Source: elliptic.c:2306-2316, ellmul_Z
 * Source: FpE.c:345-365, _FpE_mul and FpE_mul
 *
 * Uses Jacobian coordinates internally for efficiency. General ellinit records
 * over p > 3, including the legacy short-record adapter, follow ellffmul's
 * coordinate-change/FpE_mul/inverse-change route. Short records in characteristic
 * three use the native FF model and FlxqE kernel.
 * @see Deviation: General-model PARI scalar multiplication
 */
export function ellmul(
  E: FFEllipticCurve, P: FFEllipticInputPoint, n: bigint
): FFEllipticPoint;
export function ellmul(
  E: ShortWeierstrassCurve | EllipticCurve, P: EllipticPoint, n: bigint
): EllipticPoint;
export function ellmul(
  E: ShortWeierstrassCurve | EllipticCurve | FFEllipticCurve,
  P: EllipticPoint | FFEllipticInputPoint, n: bigint
): EllipticPoint | FFEllipticPoint {
  if ('field' in E) return P.isInfinity ? ellinf() : FF_ellmul(E, P, n);
  return ellmulPrime(E, P as EllipticPoint, n);
}

function ellmulPrime(
  E: ShortWeierstrassCurve | EllipticCurve,
  P: EllipticPoint,
  n: bigint
): EllipticPoint {
  if (!('type' in E)) {
    const model = scalarModel(E);
    if ('field' in model) {
      const Q = FF_ellmul(model, P, n);
      if (Q.isInfinity) return ellinf();
      const x = Q.x.value, y = Q.y.value;
      return mkpoint(typeof x === 'bigint' ? x : (x[0] ?? 0n), typeof y === 'bigint' ? y : (y[0] ?? 0n));
    }
    E = model;
  }
  if (ell_is_inf(P)) return ellinf();
  if (E.type !== EllCurveType.t_ELL_Fp || E.p === undefined || E.p <= 3n)
    throw new Error('PARI_NOT_IMPLEMENTED: ellmul over this base field');
  // elliptic.c:2289 ellffmul: transform into PARI's short model and back.
  const [a4, , ch] = ell_to_a4a6_bc(E, E.p);
  const T = FpE_changepointinv(P, ch, E.p);
  const Q = FpE_mul(T.isInfinity ? { isInfinity: true, x: null, y: null } : T, n, a4, E.p);
  return FpE_changepoint(Q.isInfinity ? ellinf() : mkpoint(Q.x!, Q.y!), ch, E.p);
}

/** PARI FpE.c:190, change from transformed to original coordinates. */
export function FpE_changepoint(
  P: EllipticPoint,
  ch: readonly [bigint, bigint, bigint, bigint],
  p: bigint
): EllipticPoint {
  if (ell_is_inf(P)) return P;
  const [u, r, s, t] = ch;
  const v = Fp_inv(u, p),
    v2 = Fp_sqr(v, p),
    v3 = Fp_mul(v, v2, p);
  const c = Fp_sub(P.x, r, p);
  return mkpoint(Fp_mul(v2, c, p), Fp_mul(v3, Fp_sub(P.y, Fp_add(Fp_mul(s, c, p), t, p), p), p));
}

/** PARI FpE.c:211, change from original to transformed coordinates. */
export function FpE_changepointinv(
  P: EllipticPoint,
  ch: readonly [bigint, bigint, bigint, bigint],
  p: bigint
): EllipticPoint {
  if (ell_is_inf(P)) return P;
  const [u, r, s, t] = ch;
  const u2 = Fp_sqr(u, p),
    u3 = Fp_mul(u, u2, p),
    c = Fp_mul(u2, P.x, p);
  return mkpoint(Fp_add(c, r, p), Fp_add(Fp_mul(u3, P.y, p), Fp_add(Fp_mul(s, c, p), t, p), p));
}
