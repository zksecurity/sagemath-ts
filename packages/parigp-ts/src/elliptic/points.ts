/**
 * @module parigp-ts/elliptic/points
 * @description Elliptic curve point operations
 *
 * Port of PARI/GP functions from:
 * - reference/pari/src/basemath/elliptic.c (ellordinate)
 * - reference/pari/src/basemath/FpE.c (random_FpE, coordinate conversions)
 * - reference/pari/src/headers/pariinl.h (ellinf, ell_is_inf)
 */

import { Fp_halve, Fp_sqrt } from '../ff.js';

/**
 * Representation of a point on an elliptic curve.
 *
 * PARI/GP represents:
 * - Point at infinity as a single-element vector [0]
 * - Finite points as [x, y]
 *
 * We use a discriminated union type for type safety.
 */
export type EllipticPoint = { isInfinity: true } | { isInfinity: false; x: bigint; y: bigint };

/**
 * The finite branch of {@link EllipticPoint}, i.e. an affine point with coordinates.
 *
 * Constructors that can never produce the point at infinity return this instead of the
 * full union, so callers reach `.x`/`.y` without re-narrowing. It remains assignable to
 * `EllipticPoint` everywhere.
 */
export type EllipticPointFinite = Extract<EllipticPoint, { isInfinity: false }>;

/**
 * Jacobian coordinates representation.
 *
 * From FpE.c: Jacobian coordinates (X : Y : Z) represent affine (X/Z^2, Y/Z^3).
 * Point at infinity is represented as Z = 0, conventionally (1 : 1 : 0).
 */
export interface JacobianPoint {
  X: bigint;
  Y: bigint;
  Z: bigint;
}

/**
 * Short Weierstrass curve parameters: y^2 = x^3 + a4*x + a6 over Fp.
 */
export interface ShortWeierstrassCurve {
  a4: bigint;
  a6: bigint;
  p: bigint;
}

/**
 * Return the point at infinity.
 *
 * Port of PARI/GP's ellinf() from pariinl.h:2754:
 *   INLINE GEN ellinf(void) { return mkvec(gen_0); }
 *
 * @returns The point at infinity
 */
export function ellinf(): EllipticPoint {
  return { isInfinity: true };
}

/**
 * Check if a point is the point at infinity.
 *
 * Port of PARI/GP's ell_is_inf() from pariinl.h:2753:
 *   INLINE int ell_is_inf(GEN z) { return lg(z) == 2; }
 *
 * @param P - Point to check
 * @returns true if P is the point at infinity
 */
export function ell_is_inf(P: EllipticPoint): P is { isInfinity: true } {
  return P.isInfinity;
}

/**
 * Create a finite point.
 *
 * @param x - x-coordinate
 * @param y - y-coordinate
 * @returns The point (x, y)
 */
export function mkpoint(x: bigint, y: bigint): EllipticPointFinite {
  return { isInfinity: false, x, y };
}

/**
 * Compute a mod n, ensuring the result is in [0, n).
 */
function mod(a: bigint, n: bigint): bigint {
  const result = a % n;
  return result < 0n ? result + n : result;
}

/**
 * Compute the modular inverse of a modulo m.
 *
 * @throws Error if gcd(a, m) != 1
 */
function modInverse(a: bigint, m: bigint): bigint {
  a = mod(a, m);
  if (a === 0n) {
    throw new Error('modular inverse of zero does not exist');
  }

  let oldR = a;
  let r = m;
  let oldS = 1n;
  let s = 0n;

  while (r !== 0n) {
    const quotient = oldR / r;
    const tempR = r;
    r = oldR - quotient * r;
    oldR = tempR;
    const tempS = s;
    s = oldS - quotient * s;
    oldS = tempS;
  }

  if (oldR !== 1n) {
    throw new Error(`modular inverse does not exist: gcd(${a}, ${m}) = ${oldR}`);
  }

  return mod(oldS, m);
}

/**
 * Fast modular exponentiation: a^n mod m.
 */
function powerMod(a: bigint, n: bigint, m: bigint): bigint {
  if (n < 0n) {
    a = modInverse(a, m);
    n = -n;
  }

  a = mod(a, m);
  if (n === 0n) return 1n;

  let result = 1n;
  let base = a;

  while (n > 0n) {
    if ((n & 1n) === 1n) {
      result = mod(result * base, m);
    }
    base = mod(base * base, m);
    n >>= 1n;
  }

  return result;
}

/**
 * Compute the Kronecker/Legendre symbol (a/p).
 *
 * @param a - Integer
 * @param p - Odd prime
 * @returns -1, 0, or 1
 */
function kronecker(a: bigint, p: bigint): bigint {
  a = mod(a, p);
  if (a === 0n) return 0n;
  if (p === 2n) return 1n;

  // Euler's criterion: (a/p) = a^((p-1)/2) mod p
  const result = powerMod(a, (p - 1n) / 2n, p);
  if (result === p - 1n) return -1n;
  return result;
}

/**
 * Generate a random bigint in range [0, n).
 */
function randomi(n: bigint): bigint {
  if (n <= 0n) {
    throw new Error('randomi requires positive argument');
  }

  if (n <= BigInt(Number.MAX_SAFE_INTEGER)) {
    return BigInt(Math.floor(Math.random() * Number(n)));
  }

  // For larger values, generate random bytes
  const bits = n.toString(2).length;
  const bytes = Math.ceil(bits / 8);
  let result: bigint;

  do {
    result = 0n;
    for (let i = 0; i < bytes; i++) {
      result = (result << 8n) | BigInt(Math.floor(Math.random() * 256));
    }
    result = result % n;
  } while (result >= n);

  return result;
}

/**
 * Find y-coordinates for a given x on an elliptic curve.
 *
 * Port of PARI/GP's ellordinate() from elliptic.c:2140-2211.
 *
 * For a short Weierstrass curve y^2 = x^3 + a4*x + a6 over Fp,
 * this computes the possible y values for a given x.
 *
 * Algorithm:
 * 1. Compute rhs = x^3 + a4*x + a6
 * 2. Check if rhs is a quadratic residue (using Kronecker symbol)
 * 3. If not, return empty array (no points with this x)
 * 4. If rhs = 0, return [0] (single point)
 * 5. Otherwise, compute sqrt(rhs) and return both square roots
 *
 * @param E - Curve parameters (a4, a6, p)
 * @param x - x-coordinate to find y for
 * @returns Array of 0, 1, or 2 y-values
 *
 * @example
 * ```typescript
 * // Curve: y^2 = x^3 + x + 1 over F_23
 * const E = { a4: 1n, a6: 1n, p: 23n };
 * ellordinate(E, 0n)  // [1n, 22n] since 1^2 = 1 = 0 + 0 + 1
 * ellordinate(E, 2n)  // [] if 2^3 + 2 + 1 = 11 is not a QR mod 23
 * ```
 */
export function ellordinate(E: ShortWeierstrassCurve, x: bigint): bigint[] {
  const { a4, a6, p } = E;

  // Normalize x to [0, p)
  x = mod(x, p);

  // PARI solves y*(y+b) = a with a = ec_f_evalx(E,x), b = ec_h_evalx(E,x).
  // In short Weierstrass form (a1 = a2 = a3 = 0):
  //   a = x^3 + a4*x + a6,  b = 0,  D = b^2 + 4a = 4a
  const x2 = mod(x * x, p);
  const a = mod(mod(x * mod(x2 + a4, p), p) + a6, p);
  const D = mod(4n * a, p);

  // D == 0: the single root y = -b/2 = 0 (elliptic.c:233-239)
  if (D === 0n) {
    return [0n];
  }

  // t_ELL_Fp branch (elliptic.c:245-251): no point when D is a non-residue
  const kr = kronecker(D, p);
  if (kr < 0n) {
    return [];
  }

  // d = Fp_sqrt(D); y1 = (d - b)/2 = d/2; y2 = y1 - d  (elliptic.c:279-282)
  const d = Fp_sqrt(D, p);
  if (d === null) {
    return [];
  }
  const y1 = Fp_halve(d, p);
  const y2 = mod(y1 - d, p);

  return [y1, y2];
}

/**
 * Generate a random nonsingular point on an elliptic curve.
 *
 * Port of PARI/GP's random_FpE() from FpE.c:370-385:
 *
 * ```c
 * GEN random_FpE(GEN a4, GEN a6, GEN p)
 * {
 *   pari_sp ltop = avma;
 *   GEN x, x2, y, rhs;
 *   do
 *   {
 *     set_avma(ltop);
 *     x   = randomi(p);
 *     x2  = Fp_sqr(x, p);
 *     rhs = Fp_add(Fp_mul(x, Fp_add(x2, a4, p), p), a6, p);
 *   } while ((!signe(rhs) && !signe(Fp_add(Fp_mulu(x2,3,p),a4,p)))
 *           || kronecker(rhs, p) < 0);
 *   y = Fp_sqrt(rhs, p);
 *   if (!y) pari_err_PRIME("random_FpE", p);
 *   return gc_GEN(ltop, mkvec2(x, y));
 * }
 * ```
 *
 * We return `Fp_sqrt(rhs, p)` exactly like PARI (the canonical, smallest
 * root); the sign is *not* randomized. `<P>` and `<-P>` are the same
 * subgroup, so every consumer (order, group structure, pairings) is
 * unaffected, and matching PARI keeps the two `random_FpE` in this package
 * consistent.
 *
 * The condition `(!signe(rhs) && !signe(Fp_add(Fp_mulu(x2,3,p),a4,p)))`
 * checks if the point would be singular:
 * - rhs = 0 means y = 0
 * - 3*x^2 + a4 = 0 means the curve has a singular point at (x, 0)
 *
 * @param E - Curve parameters (a4, a6, p)
 * @returns A random point on the curve (never the point at infinity)
 *
 * @example
 * ```typescript
 * const E = { a4: 1n, a6: 1n, p: 23n };
 * const P = random_FpE(E);
 * // P is a random point on y^2 = x^3 + x + 1 over F_23
 * ```
 */
export function random_FpE(E: ShortWeierstrassCurve): EllipticPoint {
  const { a4, a6, p } = E;

  // PARI's do/while: retry while the candidate is the singular point of a
  // singular curve, or while rhs is a non-residue.  A `continue` inside a
  // do/while jumps to the *condition*, which would accept rhs = 0, so the
  // loop is written with an explicit `while (true)` here.
  let x: bigint;
  let rhs: bigint;
  for (;;) {
    x = randomi(p);
    const x2 = mod(x * x, p);
    // rhs = x*(x^2 + a4) + a6
    rhs = mod(mod(x * mod(x2 + a4, p), p) + a6, p);

    // (!signe(rhs) && !signe(Fp_add(Fp_mulu(x2,3,p),a4,p))): singular point
    if (rhs === 0n && mod(3n * x2 + a4, p) === 0n) continue;

    if (kronecker(rhs, p) < 0n) continue;
    break;
  }

  // y = Fp_sqrt(rhs, p); pari_err_PRIME if it fails
  const y = Fp_sqrt(rhs, p);
  if (y === null) {
    throw new Error(`random_FpE: ${p} is not prime`);
  }

  return mkpoint(x, y);
}

/**
 * Convert an affine point to Jacobian coordinates.
 *
 * Port of PARI/GP's FpE_to_FpJ() from FpE.c:114-117:
 *
 * ```c
 * GEN FpE_to_FpJ(GEN P)
 * {
 *   return ell_is_inf(P) ? ellinf_FpJ()
 *        : mkvec3(icopy(gel(P,1)),icopy(gel(P,2)), gen_1);
 * }
 * ```
 *
 * Jacobian coordinates (X : Y : Z) represent affine (X/Z^2, Y/Z^3).
 * A finite point (x, y) maps to (x, y, 1).
 * The point at infinity maps to (1, 1, 0).
 *
 * @param P - Affine point
 * @returns Jacobian coordinates
 */
export function FpE_to_FpJ(P: EllipticPoint): JacobianPoint {
  if (ell_is_inf(P)) {
    // Point at infinity in Jacobian: (1 : 1 : 0)
    return { X: 1n, Y: 1n, Z: 0n };
  }

  return { X: P.x, Y: P.y, Z: 1n };
}

/**
 * Convert Jacobian coordinates to an affine point.
 *
 * Port of PARI/GP's FpJ_to_FpE() from FpE.c:121-130:
 *
 * ```c
 * GEN FpJ_to_FpE(GEN P, GEN p)
 * {
 *   if (signe(gel(P,3)) == 0) return ellinf();
 *   else
 *   {
 *     GEN Z = Fp_inv(gel(P,3), p);
 *     GEN Z2 = Fp_sqr(Z, p), Z3 = Fp_mul(Z, Z2, p);
 *     retmkvec2(Fp_mul(gel(P,1), Z2, p), Fp_mul(gel(P,2), Z3, p));
 *   }
 * }
 * ```
 *
 * Jacobian coordinates (X : Y : Z) represent affine (X/Z^2, Y/Z^3).
 * If Z = 0, return the point at infinity.
 *
 * @param P - Jacobian point
 * @param p - Prime modulus
 * @returns Affine point
 */
export function FpJ_to_FpE(P: JacobianPoint, p: bigint): EllipticPoint {
  if (P.Z === 0n) {
    return ellinf();
  }

  const Z = modInverse(P.Z, p);
  const Z2 = mod(Z * Z, p);
  const Z3 = mod(Z * Z2, p);

  return mkpoint(mod(P.X * Z2, p), mod(P.Y * Z3, p));
}

/**
 * Return the Jacobian point at infinity.
 *
 * Port of PARI/GP's ellinf_FpJ() from FpE.c:34-35:
 *
 * ```c
 * static GEN ellinf_FpJ(void)
 * { return mkvec3(gen_1, gen_1, gen_0); }
 * ```
 *
 * @returns Jacobian point at infinity (1 : 1 : 0)
 */
export function ellinf_FpJ(): JacobianPoint {
  return { X: 1n, Y: 1n, Z: 0n };
}

/**
 * Check if a Jacobian point is the point at infinity.
 *
 * @param P - Jacobian point
 * @returns true if P is at infinity (Z = 0)
 */
export function FpJ_is_inf(P: JacobianPoint): boolean {
  return P.Z === 0n;
}

/**
 * Verify that a point is on the curve.
 *
 * For curve y^2 = x^3 + a4*x + a6, checks that y^2 - x^3 - a4*x - a6 = 0.
 *
 * @param E - Curve parameters
 * @param P - Point to check
 * @returns true if P is on the curve
 */
export function FpE_isoncurve(E: ShortWeierstrassCurve, P: EllipticPoint): boolean {
  if (ell_is_inf(P)) {
    return true;
  }

  const { a4, a6, p } = E;
  const { x, y } = P;

  // Check y^2 = x^3 + a4*x + a6
  const lhs = mod(y * y, p);
  const x2 = mod(x * x, p);
  const x3 = mod(x2 * x, p);
  const rhs = mod(x3 + mod(a4 * x, p) + a6, p);

  return lhs === rhs;
}
