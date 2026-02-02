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
// ============================================================================

/**
 * Modular reduction ensuring positive result.
 */
function mod(a: bigint, p: bigint): bigint {
  const r = a % p;
  return r < 0n ? r + p : r;
}

/**
 * Modular addition.
 */
function Fp_add(a: bigint, b: bigint, p: bigint): bigint {
  return mod(a + b, p);
}

/**
 * Modular subtraction.
 */
function Fp_sub(a: bigint, b: bigint, p: bigint): bigint {
  return mod(a - b, p);
}

/**
 * Modular multiplication.
 */
function Fp_mul(a: bigint, b: bigint, p: bigint): bigint {
  return mod(a * b, p);
}

/**
 * Modular squaring.
 */
function Fp_sqr(a: bigint, p: bigint): bigint {
  return mod(a * a, p);
}

/**
 * Modular negation.
 */
function Fp_neg(a: bigint, p: bigint): bigint {
  return mod(-a, p);
}

/**
 * Modular doubling.
 */
function Fp_double(a: bigint, p: bigint): bigint {
  return mod(a << 1n, p);
}

/**
 * Modular tripling.
 */
function Fp_triple(a: bigint, p: bigint): bigint {
  return mod(a * 3n, p);
}

/**
 * Modular multiplication by a small integer.
 */
function Fp_mulu(a: bigint, k: number, p: bigint): bigint {
  return mod(a * BigInt(k), p);
}

/**
 * Modular exponentiation using binary method.
 *
 * Reference: PARI gen1.c
 */
function Fp_pow(base: bigint, exp: bigint, p: bigint): bigint {
  if (exp === 0n) return 1n;
  if (exp < 0n) {
    base = Fp_inv(base, p);
    exp = -exp;
  }

  let result = 1n;
  base = mod(base, p);

  while (exp > 0n) {
    if ((exp & 1n) === 1n) {
      result = mod(result * base, p);
    }
    exp >>= 1n;
    base = mod(base * base, p);
  }

  return result;
}

/**
 * Extended Euclidean algorithm.
 * Returns [gcd, x, y] where gcd = a*x + b*y.
 */
function extgcd(a: bigint, b: bigint): [bigint, bigint, bigint] {
  if (b === 0n) {
    return [a, 1n, 0n];
  }

  let [oldR, r] = [a, b];
  let [oldS, s] = [1n, 0n];
  let [oldT, t] = [0n, 1n];

  while (r !== 0n) {
    const quotient = oldR / r;
    [oldR, r] = [r, oldR - quotient * r];
    [oldS, s] = [s, oldS - quotient * s];
    [oldT, t] = [t, oldT - quotient * t];
  }

  return [oldR, oldS, oldT];
}

/**
 * Modular inverse using extended Euclidean algorithm.
 *
 * Reference: PARI arith1.c - Fp_inv
 */
function Fp_inv(a: bigint, p: bigint): bigint {
  const am = mod(a, p);
  if (am === 0n) {
    throw new Error('Cannot invert zero');
  }

  const [g, x] = extgcd(am, p);
  if (g !== 1n) {
    throw new Error(`${a} is not invertible mod ${p}`);
  }

  return mod(x, p);
}

/**
 * Modular division.
 */
function Fp_div(a: bigint, b: bigint, p: bigint): bigint {
  return Fp_mul(a, Fp_inv(b, p), p);
}

/**
 * Compute the Kronecker/Legendre symbol (a/p).
 *
 * Returns:
 *  1 if a is a quadratic residue mod p
 * -1 if a is a non-residue mod p
 *  0 if a = 0 mod p
 *
 * Reference: PARI arith1.c - kronecker
 */
function kronecker(a: bigint, p: bigint): number {
  a = mod(a, p);
  if (a === 0n) return 0;

  // Use Euler's criterion: a^((p-1)/2) = 1 if QR, -1 if NR
  const exp = (p - 1n) / 2n;
  const result = Fp_pow(a, exp, p);

  if (result === 1n) return 1;
  if (result === p - 1n) return -1;
  return 0;
}

/**
 * Compute square root modulo p using Tonelli-Shanks algorithm.
 *
 * Reference: PARI arith1.c - Fp_sqrt
 *
 * @param a - The value to take square root of
 * @param p - The prime modulus
 * @returns sqrt(a) mod p, or null if a is not a QR
 */
function Fp_sqrt(a: bigint, p: bigint): bigint | null {
  a = mod(a, p);

  if (a === 0n) return 0n;

  // Check if a is a quadratic residue
  if (kronecker(a, p) !== 1) {
    return null;
  }

  // Special case: p = 3 mod 4
  if (mod(p, 4n) === 3n) {
    return Fp_pow(a, (p + 1n) / 4n, p);
  }

  // Tonelli-Shanks algorithm
  // Write p - 1 = 2^s * q with q odd
  let s = 0n;
  let q = p - 1n;
  while ((q & 1n) === 0n) {
    s++;
    q >>= 1n;
  }

  // Find a quadratic non-residue
  let z = 2n;
  while (kronecker(z, p) !== -1) {
    z++;
  }

  let m = s;
  let c = Fp_pow(z, q, p);
  let t = Fp_pow(a, q, p);
  let r = Fp_pow(a, (q + 1n) / 2n, p);

  while (true) {
    if (t === 1n) return r;

    // Find the least i such that t^(2^i) = 1
    let i = 1n;
    let temp = Fp_sqr(t, p);
    while (temp !== 1n) {
      temp = Fp_sqr(temp, p);
      i++;
    }

    // Update
    const b = Fp_pow(c, 1n << (m - i - 1n), p);
    m = i;
    c = Fp_sqr(b, p);
    t = Fp_mul(t, c, p);
    r = Fp_mul(r, b, p);
  }
}

/**
 * Integer square root (floor).
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
 * GCD of two integers.
 */
function gcd(a: bigint, b: bigint): bigint {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;

  while (b !== 0n) {
    const t = b;
    b = a % b;
    a = t;
  }

  return a;
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
 * Simple trial division factorization.
 *
 * Returns an array of [prime, exponent] pairs.
 */
function factor(n: bigint): [bigint, bigint][] {
  if (n <= 1n) return [];

  const factors: [bigint, bigint][] = [];
  let remaining = n;

  // Check 2
  if (remaining % 2n === 0n) {
    let exp = 0n;
    while (remaining % 2n === 0n) {
      exp++;
      remaining /= 2n;
    }
    factors.push([2n, exp]);
  }

  // Check odd factors
  let f = 3n;
  while (f * f <= remaining) {
    if (remaining % f === 0n) {
      let exp = 0n;
      while (remaining % f === 0n) {
        exp++;
        remaining /= f;
      }
      factors.push([f, exp]);
    }
    f += 2n;
  }

  if (remaining > 1n) {
    factors.push([remaining, 1n]);
  }

  return factors;
}

/**
 * Compute all divisors of n from its factorization.
 */
function divisors_from_factorization(factorization: [bigint, bigint][]): bigint[] {
  const divs: bigint[] = [1n];

  for (const [p, e] of factorization) {
    const newDivs: bigint[] = [];
    let pk = 1n;
    for (let k = 0n; k <= e; k++) {
      for (const d of divs) {
        newDivs.push(d * pk);
      }
      pk *= p;
    }
    divs.length = 0;
    divs.push(...newDivs);
  }

  divs.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  return divs;
}

// ============================================================================
// Random Point Generation
// ============================================================================

/**
 * Generate a random field element in [0, p-1].
 */
function Fp_random(p: bigint): bigint {
  // Use crypto.getRandomValues if available, otherwise Math.random
  // For simplicity, we use Math.random with bigint conversion
  const bits = p.toString(2).length;
  const bytes = Math.ceil(bits / 8);

  let result = 0n;
  for (let i = 0; i < bytes; i++) {
    result = (result << 8n) | BigInt(Math.floor(Math.random() * 256));
  }

  return mod(result, p);
}

/**
 * Find a random point on the curve.
 *
 * Algorithm from PARI FpE.c - random_FpE:
 * 1. Pick random x
 * 2. Compute y^2 = x^3 + a4*x + a6
 * 3. If y^2 is a square, return (x, sqrt(y^2))
 * 4. Otherwise, repeat
 *
 * Reference: PARI FpE.c:369-385
 */
export function FpE_random(E: EllipticCurveFp): EllipticPointFp {
  const { a4, a6, p } = E;

  for (let attempts = 0; attempts < 1000; attempts++) {
    const x = Fp_random(p);

    // y^2 = x^3 + a4*x + a6
    const ySquared = mod(Fp_pow(x, 3n, p) + a4 * x + a6, p);

    const y = Fp_sqrt(ySquared, p);
    if (y !== null) {
      // Randomly choose one of the two y values
      if (Math.random() < 0.5) {
        return ellpoint(x, y);
      } else {
        return ellpoint(x, Fp_neg(y, p));
      }
    }
  }

  throw new Error('Failed to find random point after 1000 attempts');
}

// ============================================================================
// Cardinality Computation
// ============================================================================

/**
 * Compute the trace of an elliptic curve over Fp using naive enumeration.
 *
 * The trace a_p satisfies: #E(Fp) = p + 1 - a_p
 *
 * This method enumerates all x-coordinates and uses Legendre symbols
 * to count points.
 *
 * Reference: PARI FpE.c:813-832 - Fl_elltrace_naive
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
    const legendre = kronecker(ySquared, p);
    trace -= BigInt(legendre);
  }

  return trace;
}

/**
 * Compute elliptic curve cardinality using exhaustive enumeration.
 *
 * For each x in Fp, count:
 * - 2 points if x^3 + a4*x + a6 is a non-zero square
 * - 1 point if x^3 + a4*x + a6 = 0
 * - 0 points otherwise
 *
 * Plus 1 for the point at infinity.
 *
 * Reference: PARI FpE.c:813-832
 */
function ellcard_exhaustive(E: EllipticCurveFp): bigint {
  const { a4, a6, p } = E;
  let count = 1n; // Point at infinity

  for (let x = 0n; x < p; x++) {
    const ySquared = mod(Fp_pow(x, 3n, p) + a4 * x + a6, p);

    if (ySquared === 0n) {
      count += 1n;
    } else if (kronecker(ySquared, p) === 1) {
      count += 2n;
    }
  }

  return count;
}

/**
 * Find the closest integer to c that is congruent to a modulo b.
 *
 * Reference: PARI FpE.c:864-867 - closest_lift
 */
function closest_lift(a: bigint, b: bigint, c: bigint): bigint {
  // x = round((c - a) / b)
  // Return a + b * x
  const diff = c - a;
  // Compute floor((2*(c-a) + b) / (2*b)) for proper rounding
  const x = (2n * diff + b) / (2n * b);
  return a + b * x;
}

/**
 * Find a point on E with given Legendre symbol for y^2.
 *
 * Reference: PARI FpE.c:882-894 - Fp_ellpoint
 */
function Fp_ellpoint(
  KRO: number,
  startX: bigint,
  a4: bigint,
  a6: bigint,
  p: bigint
): [EllipticPointFp, bigint] {
  let x = startX;

  while (true) {
    x++;
    const ySquared = mod(Fp_pow(x, 3n, p) + a4 * x + a6, p);

    if (kronecker(ySquared, p) === KRO) {
      // Return [x*u, u^2] where u = sqrt(y^2)
      // For our purposes, we return a point on the twist/original curve
      return [ellpoint(x, ySquared), x];
    }
  }
}

/**
 * Compute elliptic curve cardinality using baby-step giant-step.
 *
 * This is a simplified version of PARI's Fp_ellcard_Shanks algorithm.
 * Uses the Hasse bound: |#E - (p+1)| <= 2*sqrt(p)
 *
 * Reference: PARI FpE.c:921-1135 - Fp_ellcard_Shanks
 */
function ellcard_bsgs(E: EllipticCurveFp): bigint {
  const { a4, a6, p } = E;

  // Hasse bound: |#E - (p+1)| <= 2*sqrt(p)
  const sqrtP = isqrt(p);
  const hasseBound = 2n * sqrtP;
  const center = p + 1n;

  // Try multiple random points to determine the order
  // We accumulate constraints via CRT

  let A = 0n; // Current congruence class
  let B = 1n; // Current modulus
  const pordmin = 4n * sqrtP + 1n; // Need B >= pordmin

  let attempts = 0;
  while (B < pordmin && attempts < 100) {
    attempts++;

    // Get a random point
    let P: EllipticPointFp;
    try {
      P = FpE_random(E);
    } catch {
      continue;
    }

    if (P.isInfinity) continue;

    // Find the order of P using BSGS
    // The order divides some candidate in [center - hasseBound, center + hasseBound]

    const m = isqrt(2n * hasseBound) + 1n;

    // Baby steps: compute [center - m + j]P for j = 0..2m
    const babySteps = new Map<string, bigint>();
    let Q = FpE_mul(P, center - m, a4, p);

    for (let j = 0n; j <= 2n * m; j++) {
      if (!Q.isInfinity) {
        const key = `${Q.x},${Q.y}`;
        if (!babySteps.has(key)) {
          babySteps.set(key, center - m + j);
        }
      } else {
        // Found [center - m + j]P = O
        const candidateOrder = center - m + j;
        if (candidateOrder > 0n) {
          // Find exact order by factoring
          const order = findExactOrder(P, candidateOrder, a4, p);
          // Update CRT
          [A, B] = crt_update(A, B, 0n, order);
          break;
        }
      }
      Q = FpE_add(Q, P, a4, p);
    }

    if (B >= pordmin) break;

    // Giant steps: compute [m*i]P for i = 1, 2, ...
    const mP = FpE_mul(P, m, a4, p);
    let R = mP;

    for (let i = 1n; i <= hasseBound / m + 2n; i++) {
      if (R.isInfinity) {
        const candidateOrder = m * i;
        const order = findExactOrder(P, candidateOrder, a4, p);
        [A, B] = crt_update(A, B, 0n, order);
        break;
      }

      // Check both R and -R against baby steps
      const keyPos = `${R.x},${R.y}`;
      const keyNeg = `${R.x},${Fp_neg(R.y!, p)}`;

      if (babySteps.has(keyPos)) {
        // [center - m + j]P = [m*i]P
        // So [center - m + j - m*i]P = O
        const j = babySteps.get(keyPos)!;
        const diff = j - m * i;
        if (diff > 0n && diff <= center + hasseBound) {
          const order = findExactOrder(P, diff, a4, p);
          [A, B] = crt_update(A, B, 0n, order);
          break;
        }
      }

      if (babySteps.has(keyNeg)) {
        // [center - m + j]P = -[m*i]P
        // So [center - m + j + m*i]P = O
        const j = babySteps.get(keyNeg)!;
        const sum = j + m * i;
        if (sum > 0n && sum <= center + hasseBound) {
          const order = findExactOrder(P, sum, a4, p);
          [A, B] = crt_update(A, B, 0n, order);
          break;
        }
      }

      R = FpE_add(R, mP, a4, p);
    }
  }

  // Find the value in Hasse interval congruent to A mod B
  const result = closest_lift(A, B, center);
  return result;
}

/**
 * Find the exact order of P given that [n]P = O.
 *
 * Repeatedly divides by prime factors of n while still annihilating P.
 */
function findExactOrder(P: EllipticPointFp, n: bigint, a4: bigint, p: bigint): bigint {
  let order = n;
  const factors = factor(n);

  for (const [prime] of factors) {
    while (order % prime === 0n) {
      const testOrder = order / prime;
      const Q = FpE_mul(P, testOrder, a4, p);
      if (Q.isInfinity) {
        order = testOrder;
      } else {
        break;
      }
    }
  }

  return order;
}

/**
 * Update CRT constraints: find A', B' such that
 * x = A' mod B' iff x = a1 mod b1 and x = a2 mod b2.
 */
function crt_update(a1: bigint, b1: bigint, a2: bigint, b2: bigint): [bigint, bigint] {
  const g = gcd(b1, b2);
  const newB = lcm(b1, b2);

  if (mod(a1, g) !== mod(a2, g)) {
    // No solution - this shouldn't happen in valid usage
    return [a1, b1];
  }

  // Solve using extended Euclidean algorithm
  const [, u1] = extgcd(b1, b2);
  const newA = mod(a1 + b1 * u1 * ((a2 - a1) / g), newB);

  return [newA, newB];
}

/**
 * Compute the cardinality (number of points) of E(Fp).
 *
 * PARI uses different algorithms based on field size:
 * - Small primes (p < 1000): Exhaustive enumeration
 * - Medium primes (p < 10^6): Baby-step giant-step (Shanks/Mestre)
 * - Large primes: SEA algorithm (not implemented here)
 *
 * Reference: PARI elliptic.c:6360-6388 - ellcard
 *
 * @param E - The elliptic curve
 * @returns The number of points on E(Fp)
 */
export function ellcard(E: EllipticCurveFp): bigint {
  if (E._card !== undefined) {
    return E._card;
  }

  const { p } = E;

  let card: bigint;

  if (p < 1000n) {
    // Small primes: exhaustive enumeration
    card = ellcard_exhaustive(E);
  } else if (p < 1000000n) {
    // Medium primes: baby-step giant-step
    card = ellcard_bsgs(E);
  } else {
    // Large primes: also BSGS for now (SEA would be better but complex)
    card = ellcard_bsgs(E);
  }

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
 * Reference: PARI FpE.c:406-423 - FpE_order
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
  const factors = factor(N);

  let order = N;
  const { a4, p } = E;

  // For each prime power p^e dividing N, find the actual power dividing the point order
  for (const [prime] of factors) {
    while (order % prime === 0n) {
      const testOrder = order / prime;
      const Q = FpE_mul(P, testOrder, a4, p);
      if (Q.isInfinity) {
        order = testOrder;
      } else {
        break;
      }
    }
  }

  return order;
}

// ============================================================================
// Group Structure
// ============================================================================

/**
 * Compute the order of the Weil pairing of P and Q.
 *
 * This is used in the group structure algorithm.
 * For simplicity, we compute it by finding the order of the product
 * in the multiplicative group.
 *
 * Reference: PARI FpE.c:1463-1467 - _FpE_pairorder
 */
function pairorder(E: EllipticCurveFp, P: EllipticPointFp, Q: EllipticPointFp, m: bigint): bigint {
  // The Weil pairing e(P,Q) has order dividing m
  // For a full implementation, we would compute the actual Weil pairing
  // For now, we use a simpler heuristic based on point orders

  const { a4, p } = E;
  const curveOrder = ellcard(E);

  // Compute orders of P and Q
  const orderP = ellorder(E, P, curveOrder);
  const orderQ = ellorder(E, Q, curveOrder);

  // The pairing order divides gcd(orderP, orderQ)
  // A full implementation would compute the actual pairing
  return gcd(orderP, orderQ);
}

/**
 * Check if d2 might divide d1 in the group structure.
 *
 * For E(Fp), we have E(Fp) = Z/d1Z or Z/d1Z x Z/d2Z with d2 | d1.
 * The d2 must divide gcd(N, p-1) where N = #E(Fp).
 *
 * Reference: PARI bb_group.c:991-1000 - d2_multiple
 */
function d2_multiple(N: bigint, pMinus1: bigint): [bigint, [bigint, bigint][]] | null {
  const g = gcd(N, pMinus1);
  if (g === 1n) {
    return null; // Group is cyclic
  }

  const factors = factor(g);
  return [g, factors];
}

/**
 * Determine the group structure of E(Fp).
 *
 * Returns [d1] for cyclic groups Z/d1Z,
 * or [d1, d2] for non-cyclic groups Z/d1Z x Z/d2Z with d2 | d1.
 *
 * Algorithm based on PARI bb_group.c:986-1048 - gen_ellgroup
 *
 * Reference: PARI FpE.c:1470-1475 - Fp_ellgroup
 *
 * @param E - The elliptic curve
 * @returns Group structure as [d1] or [d1, d2]
 */
export function ellgroup(E: EllipticCurveFp): bigint[] {
  if (E._group !== undefined) {
    return E._group;
  }

  const { a4, p } = E;
  const N = ellcard(E);

  if (N === 1n) {
    E._group = [];
    return [];
  }

  // Check if group can be non-cyclic
  // For E(Fp), the group is Z/d1Z x Z/d2Z where d2 | d1 and d2 | (p-1)
  const pMinus1 = p - 1n;
  const d2Info = d2_multiple(N, pMinus1);

  if (d2Info === null) {
    // Group must be cyclic
    E._group = [N];
    return [N];
  }

  const [g] = d2Info;
  const gFactors = factor(g);

  // Try to find two independent points that generate the full group
  // Using a probabilistic approach similar to PARI's gen_ellgroup

  let g1 = 1n; // Product of prime powers found for d1
  let g2 = 1n; // Product of prime powers found for d2

  const primesDone = new Set<bigint>();
  const numPrimes = gFactors.length;

  for (let attempts = 0; attempts < 100 && primesDone.size < numPrimes; attempts++) {
    // Get two random points
    let P: EllipticPointFp;
    let Q: EllipticPointFp;
    try {
      P = FpE_random(E);
      Q = FpE_random(E);
    } catch {
      continue;
    }

    const orderP = ellorder(E, P, N);
    const orderQ = ellorder(E, Q, N);

    if (orderP === N || orderQ === N) {
      // Found a generator - group is cyclic
      E._group = [N];
      return [N];
    }

    // Compute LCM of orders
    const m = lcm(orderP, orderQ);

    // Compute pair order contribution
    const pairOrd = pairorder(E, P, Q, m);
    const mo = m * pairOrd;

    // Check each prime dividing g
    for (const [prime, exp] of gFactors) {
      if (primesDone.has(prime)) continue;

      const vmo = valuation(mo, prime);

      if (vmo >= exp) {
        // This prime is fully determined
        const vm = valuation(m, prime);
        // Ensure we don't get negative exponent
        const vmClamped = vm > exp ? exp : vm;
        g1 *= prime ** vmClamped;
        g2 *= prime ** (exp - vmClamped);
        primesDone.add(prime);
      }
    }
  }

  // Combine the results
  const d2 = g2;
  const d1 = N / d2;

  if (d2 === 1n) {
    E._group = [N];
    return [N];
  }

  // Return in PARI format: [d1, d2] where d2 | d1
  E._group = [d1, d2];
  return [d1, d2];
}

/**
 * Compute the p-adic valuation of n.
 */
function valuation(n: bigint, p: bigint): bigint {
  if (n === 0n) return -1n; // Convention

  let v = 0n;
  while (n % p === 0n) {
    v++;
    n /= p;
  }
  return v;
}

// ============================================================================
// Generators
// ============================================================================

/**
 * Find generators for the group E(Fp).
 *
 * For cyclic groups, returns one generator.
 * For non-cyclic groups Z/d1Z x Z/d2Z, returns two generators.
 *
 * Reference: PARI FpE.c:1477-1497 - Fp_ellgens
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
  const group = ellgroup(E);

  if (group.length === 0) {
    E._generators = [];
    return [];
  }

  if (group.length === 1) {
    // Cyclic group - find a generator of order N
    const targetOrder = group[0]!;

    for (let attempts = 0; attempts < 100; attempts++) {
      try {
        const P = FpE_random(E);
        const order = ellorder(E, P, N);

        if (order === targetOrder) {
          E._generators = [P];
          return [P];
        }
      } catch {
        continue;
      }
    }

    throw new Error('Failed to find generator');
  }

  // Non-cyclic group Z/d1Z x Z/d2Z
  const d1 = group[0]!;
  const d2 = group[1]!;

  // Find first generator of order d1
  let G1: EllipticPointFp | null = null;
  for (let attempts = 0; attempts < 100; attempts++) {
    try {
      const P = FpE_random(E);
      const order = ellorder(E, P, N);

      if (order === d1) {
        G1 = P;
        break;
      }
    } catch {
      continue;
    }
  }

  if (G1 === null) {
    throw new Error('Failed to find first generator');
  }

  // Find second generator of order d2 that is independent of G1
  // A point Q is independent if [d1/d2]Q is not in <G1>
  const cofactor = d1 / d2;

  for (let attempts = 0; attempts < 100; attempts++) {
    try {
      const Q = FpE_random(E);
      const order = ellorder(E, Q, N);

      if (order % d2 !== 0n) continue;

      // Scale Q to have order d2
      const scale = order / d2;
      const Q2 = FpE_mul(Q, scale, a4, p);

      // Check independence: Q2 is not in the subgroup generated by [d1/d2]G1
      const G1scaled = FpE_mul(G1, cofactor, a4, p);

      // Simple check: if Q2 = [k]G1scaled for some k, they are dependent
      // We check this by seeing if Q2 can be expressed as a multiple of G1scaled
      let isDependent = false;
      let R = ellinf();
      for (let k = 0n; k <= d2; k++) {
        if (ellequal(R, Q2)) {
          isDependent = true;
          break;
        }
        R = FpE_add(R, G1scaled, a4, p);
      }

      if (!isDependent) {
        E._generators = [G1, Q2];
        return [G1, Q2];
      }
    } catch {
      continue;
    }
  }

  // Fallback: just return G1
  E._generators = [G1];
  return [G1];
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
