/**
 * @module parigp-ts/elliptic/advanced
 * @description Advanced elliptic curve functions (stubs for future implementation)
 *
 * This module provides stubs for advanced PARI/GP elliptic curve operations:
 * - Discrete logarithm (elllog)
 * - Pairings (elltatepairing, ellweilpairing)
 * - SEA algorithm (ellcard_sea)
 * - Isogenies (ellisogeny, ellisogenyapply)
 * - Frobenius endomorphism (ellfrobenius)
 *
 * Port of PARI/GP functions from:
 * - reference/pari/src/basemath/elliptic.c:6244-6299 (elllog, pairings)
 * - reference/pari/src/basemath/FpE.c:426-443 (FpE_log)
 * - reference/pari/src/basemath/FpE.c:589-621 (FpE_weilpairing, FpE_tatepairing)
 * - reference/pari/src/basemath/ellsea.c:2112-2113 (Fp_ellcard_SEA)
 * - reference/pari/src/basemath/ellisog.c:631-669 (ellisogeny)
 * - reference/pari/src/basemath/ellisog.c:187-213 (ellisogenyapply)
 */

import {
  type EllipticCurveFp,
  type EllipticPointFp,
  ellorder,
  ellcard,
  ellinf,
  ellequal,
  FpE_add,
  FpE_mul,
  FpE_neg,
  ell_is_inf,
} from './group.js';

// ============================================================================
// Types
// ============================================================================

/**
 * Isogeny representation.
 *
 * PARI/GP represents an isogeny as [f, g, h] where:
 * - f(x) = phi_x(x) / h(x)^2 is the x-coordinate map
 * - g(x, y) = phi_y(x, y) / h(x)^3 is the y-coordinate map
 * - h(x) is a polynomial whose roots are the x-coordinates of the kernel
 *
 * Reference: ellisog.c:215-218
 */
export interface Isogeny {
  /** Numerator polynomial for x-coordinate: x' = f(x) / h(x)^2 */
  f: bigint[];
  /** Numerator polynomial for y-coordinate (in x and y): y' = g(x,y) / h(x)^3 */
  g: bigint[][];
  /** Kernel polynomial: h(x) */
  h: bigint[];
  /** Source curve */
  domain: EllipticCurveFp;
  /** Target curve (codomain) */
  codomain: EllipticCurveFp;
}

/**
 * Simplified isogeny for point application.
 *
 * For applying an isogeny to a point, we use the evaluated form.
 */
export interface IsogenyMap {
  /** Polynomial f(x) */
  f: bigint[];
  /** Polynomial g(x, y) represented as coefficients */
  g: bigint[][];
  /** Polynomial h(x) */
  h: bigint[];
}

// ============================================================================
// Discrete Logarithm
// ============================================================================

// ============================================================================
// Helper Functions for Discrete Logarithm
// ============================================================================

/**
 * Modular reduction ensuring positive result.
 */
function mod(a: bigint, p: bigint): bigint {
  const r = a % p;
  return r < 0n ? r + p : r;
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
 * Trial division factorization.
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
 * Baby-step giant-step algorithm for discrete logarithm in a group of prime order.
 *
 * Computes n such that P = [n]G, where G has order q (prime).
 *
 * Reference: PARI/GP gen_Shanks_log in bb_group.c:494-542
 *
 * @param P - Target point
 * @param G - Generator point
 * @param q - Order of G (must be prime)
 * @param E - Elliptic curve
 * @returns n such that P = [n]G, or throws if not found
 */
function bsgs_log(
  P: EllipticPointFp,
  G: EllipticPointFp,
  q: bigint,
  E: EllipticCurveFp
): bigint {
  const { a4, p } = E;

  // Handle trivial cases
  if (ell_is_inf(P)) {
    return 0n;
  }
  if (ellequal(P, G)) {
    return 1n;
  }

  // Compute m = ceil(sqrt(q))
  const m = isqrt(q) + 1n;

  // Baby steps: compute G^(-j) * P for j = 0, 1, ..., m-1
  // Store in a hash table: hash(point) -> j
  const babySteps = new Map<string, bigint>();
  const Ginv = FpE_neg(G, p);

  let current = P;
  for (let j = 0n; j < m; j++) {
    if (ell_is_inf(current)) {
      // Found: P = [j]G
      return j;
    }
    const key = `${current.x},${current.y}`;
    if (!babySteps.has(key)) {
      babySteps.set(key, j);
    }
    current = FpE_add(current, Ginv, a4, p);
  }

  // Giant steps: compute [m*i]G for i = 1, 2, ..., m
  const mG = FpE_mul(G, m, a4, p);
  let giant = mG;

  for (let i = 1n; i <= m + 1n; i++) {
    if (ell_is_inf(giant)) {
      // Check if 0 is in baby steps (P = [j]G where j = m*i)
      const infKey = 'null,null';
      if (babySteps.has(infKey)) {
        return babySteps.get(infKey)!;
      }
    } else {
      const key = `${giant.x},${giant.y}`;
      if (babySteps.has(key)) {
        // Found: [m*i]G = G^(-j) * P
        // So P = [m*i + j]G
        const j = babySteps.get(key)!;
        const result = mod(m * i + j, q);
        return result;
      }
    }
    giant = FpE_add(giant, mG, a4, p);
  }

  throw new Error('elllog: P is not in the subgroup generated by G');
}

/**
 * Discrete logarithm in a prime power order subgroup using the
 * lifting method (Pohlig-Hellman for prime powers).
 *
 * For a group of order q^e, computes n = n_0 + n_1*q + n_2*q^2 + ... + n_{e-1}*q^{e-1}
 * where each n_i is found using BSGS on the q-torsion subgroup.
 *
 * Reference: PARI/GP gen_PH_log in bb_group.c:605-660
 *
 * @param P - Target point
 * @param G - Generator point
 * @param q - Prime base
 * @param e - Exponent (order of G is q^e)
 * @param E - Elliptic curve
 * @returns n such that P = [n]G
 */
function pohlig_hellman_prime_power(
  P: EllipticPointFp,
  G: EllipticPointFp,
  q: bigint,
  e: bigint,
  E: EllipticCurveFp
): bigint {
  const { a4, p } = E;

  // Handle trivial cases
  if (ell_is_inf(P)) {
    return 0n;
  }
  if (e === 1n) {
    return bsgs_log(P, G, q, E);
  }

  // Precompute powers of q: q^0, q^1, ..., q^e
  const qPowers: bigint[] = [1n];
  for (let i = 1n; i <= e; i++) {
    qPowers.push(qPowers[Number(i - 1n)]! * q);
  }

  // Compute Gq = [q^(e-1)]G, which has order q
  const Gq = FpE_mul(G, qPowers[Number(e - 1n)]!, a4, p);

  // Compute the discrete log digit by digit
  let n = 0n;
  let currentP = P;
  const Ginv = FpE_neg(G, p);

  for (let j = 0n; j < e; j++) {
    // Compute [q^(e-1-j)]currentP, which has the form [n_j]Gq
    const Pj = FpE_mul(currentP, qPowers[Number(e - 1n - j)]!, a4, p);

    // Find n_j using BSGS: Pj = [n_j]Gq
    const nj = bsgs_log(Pj, Gq, q, E);

    // Update n = n + n_j * q^j
    n = n + nj * qPowers[Number(j)]!;

    // Update currentP = currentP - [n_j * q^j]G
    if (nj !== 0n) {
      const correction = FpE_mul(Ginv, nj * qPowers[Number(j)]!, a4, p);
      currentP = FpE_add(currentP, correction, a4, p);
    }
  }

  return n;
}

/**
 * Chinese Remainder Theorem for combining residues.
 *
 * @param residues - Array of residues
 * @param moduli - Array of pairwise coprime moduli
 * @returns x such that x = residues[i] mod moduli[i] for all i
 */
function crt(residues: bigint[], moduli: bigint[]): bigint {
  if (residues.length === 0) {
    return 0n;
  }
  if (residues.length === 1) {
    return mod(residues[0]!, moduli[0]!);
  }

  let result = residues[0]!;
  let modProduct = moduli[0]!;

  for (let i = 1; i < residues.length; i++) {
    const mi = moduli[i]!;
    const ri = residues[i]!;

    // Solve: result + modProduct * x = ri (mod mi)
    // x = (ri - result) * modProduct^(-1) (mod mi)
    const diff = mod(ri - result, mi);
    const inv = Fp_inv(modProduct, mi);
    const x = mod(diff * inv, mi);

    result = result + modProduct * x;
    modProduct = modProduct * mi;
  }

  return mod(result, modProduct);
}

/**
 * Compute the discrete logarithm of P with respect to base G.
 *
 * Given points P and G on an elliptic curve E, computes n such that P = [n]G
 * (scalar multiplication), if such n exists.
 *
 * PARI/GP uses the Pohlig-Hellman algorithm for this computation,
 * which exploits the factorization of the group order to reduce
 * the discrete log problem to smaller subgroups.
 *
 * @param E - Elliptic curve over Fp
 * @param P - Point whose discrete log we want to compute
 * @param G - Base point (generator)
 * @param order - Optional: Order of G (improves efficiency if known)
 * @returns n such that P = [n]G
 * @throws Error if P is not in the subgroup generated by G
 *
 * @see Reference: pari/src/basemath/elliptic.c:6244-6259
 * @see Reference: pari/src/basemath/FpE.c:426-443 (FpE_log)
 *
 * @example
 * ```typescript
 * const E = { a4: 0n, a6: 7n, p: 23n };
 * const G = { x: 1n, y: 7n, isInfinity: false };
 * const P = FpE_mul(E, G, 5n); // P = [5]G
 * const n = elllog(E, P, G); // returns 5n
 * ```
 */
export function elllog(
  E: EllipticCurveFp,
  P: EllipticPointFp,
  G: EllipticPointFp,
  order?: bigint
): bigint {
  const { a4, p } = E;

  // Handle trivial cases
  if (ell_is_inf(P)) {
    return 0n;
  }
  if (ellequal(P, G)) {
    return 1n;
  }
  if (ell_is_inf(G)) {
    // If G is infinity, P can only be infinity (which we handled above)
    throw new Error('elllog: P is not in the subgroup generated by G (G is identity)');
  }

  // Get the order of G
  const ord = order ?? ellorder(E, G, ellcard(E));

  if (ord === 1n) {
    // G is the identity, which we handled above
    throw new Error('elllog: P is not in the subgroup generated by G');
  }

  // Factor the order
  const factors = factor(ord);

  if (factors.length === 0) {
    throw new Error('elllog: invalid order');
  }

  // If order is prime, use BSGS directly
  if (factors.length === 1 && factors[0]![1] === 1n) {
    return bsgs_log(P, G, ord, E);
  }

  // Pohlig-Hellman: compute discrete log modulo each prime power
  const residues: bigint[] = [];
  const moduli: bigint[] = [];

  for (const [q, e] of factors) {
    const qe = q ** e;
    const cofactor = ord / qe;

    // Project P and G to the q^e torsion subgroup
    // [cofactor]P and [cofactor]G have order dividing q^e
    const Pq = FpE_mul(P, cofactor, a4, p);
    const Gq = FpE_mul(G, cofactor, a4, p);

    // Check if Gq is identity (shouldn't happen if ord is correct)
    if (ell_is_inf(Gq)) {
      // If Pq is not identity, no solution exists
      if (!ell_is_inf(Pq)) {
        throw new Error('elllog: P is not in the subgroup generated by G');
      }
      // Both are identity, so any value works; contribute 0
      residues.push(0n);
      moduli.push(qe);
      continue;
    }

    // Compute discrete log in the q^e subgroup
    const nq = pohlig_hellman_prime_power(Pq, Gq, q, e, E);
    residues.push(nq);
    moduli.push(qe);
  }

  // Combine using CRT
  const result = crt(residues, moduli);

  // Verify the result (optional but helps catch bugs)
  const check = FpE_mul(G, result, a4, p);
  if (!ellequal(check, P)) {
    throw new Error('elllog: internal error - result verification failed');
  }

  return result;
}

// ============================================================================
// Pairings - Helper Functions
// ============================================================================

/**
 * Create a finite point.
 */
function ellpoint(x: bigint, y: bigint): EllipticPointFp {
  return { x, y, isInfinity: false };
}

/**
 * Point doubling with slope output.
 *
 * Reference: FpE.c:257-276 (FpE_dbl_slope)
 */
function FpE_dbl_slope(
  P: EllipticPointFp,
  a4: bigint,
  p: bigint
): { point: EllipticPointFp; slope: bigint | null } {
  if (ell_is_inf(P)) {
    return { point: ellinf(), slope: null };
  }

  const x = P.x!;
  const y = P.y!;

  // Point of order 2: y = 0
  if (y === 0n) {
    return { point: ellinf(), slope: null };
  }

  // slope = (3*x^2 + a4) / (2*y)
  const x2 = mod(x * x, p);
  const num = mod(3n * x2 + a4, p);
  const den = mod(2n * y, p);
  const slope = mod(num * Fp_inv(den, p), p);

  // x' = slope^2 - 2*x
  const slopeSquared = mod(slope * slope, p);
  const x3 = mod(slopeSquared - 2n * x, p);
  const x3Normalized = x3 < 0n ? x3 + p : x3;

  // y' = slope*(x - x') - y
  const y3 = mod(slope * mod(x - x3Normalized, p) - y, p);
  const y3Normalized = y3 < 0n ? y3 + p : y3;

  return { point: ellpoint(x3Normalized, y3Normalized), slope };
}

/**
 * Point addition with slope output.
 *
 * Reference: FpE.c:279-306 (FpE_add_slope)
 */
function FpE_add_slope(
  P: EllipticPointFp,
  R: EllipticPointFp,
  a4: bigint,
  p: bigint
): { point: EllipticPointFp; slope: bigint | null } {
  if (ell_is_inf(P)) {
    return { point: R.isInfinity ? ellinf() : ellpoint(R.x!, R.y!), slope: null };
  }
  if (ell_is_inf(R)) {
    return { point: ellpoint(P.x!, P.y!), slope: null };
  }

  const x1 = P.x!;
  const y1 = P.y!;
  const x2 = R.x!;
  const y2 = R.y!;

  if (x1 === x2) {
    if (y1 === y2) {
      // Same point, use doubling
      return FpE_dbl_slope(P, a4, p);
    }
    // Inverse points
    return { point: ellinf(), slope: null };
  }

  // slope = (y2 - y1) / (x2 - x1)
  const dy = mod(y2 - y1, p);
  const dx = mod(x2 - x1, p);
  const slope = mod(dy * Fp_inv(dx, p), p);

  // x' = slope^2 - x1 - x2
  const slopeSquared = mod(slope * slope, p);
  const x3 = mod(slopeSquared - x1 - x2, p);
  const x3Normalized = x3 < 0n ? x3 + p : x3;

  // y' = slope*(x1 - x') - y1
  const y3 = mod(slope * mod(x1 - x3Normalized, p) - y1, p);
  const y3Normalized = y3 < 0n ? y3 + p : y3;

  return { point: ellpoint(x3Normalized, y3Normalized), slope };
}

/**
 * Evaluate the vertical line through P at point Q.
 *
 * The vertical line through P is x - x_P = 0.
 * Evaluating at Q gives x_Q - x_P.
 *
 * Special cases from FpE.c:453-462:
 * - If P is infinity, return 1
 * - If Q.x != P.x, return Q.x - P.x
 * - If P.y != 0, return 1
 * - Otherwise return 1/(3*P.x^2 + a4) (the tangent slope denominator)
 *
 * Reference: FpE.c:453-462 (FpE_vert)
 */
function FpE_vert(
  P: EllipticPointFp,
  Q: EllipticPointFp,
  a4: bigint,
  p: bigint
): bigint {
  if (ell_is_inf(P)) {
    return 1n;
  }

  const Px = P.x!;
  const Py = P.y!;
  const Qx = Q.x!;

  if (Qx !== Px) {
    const result = mod(Qx - Px, p);
    return result < 0n ? result + p : result;
  }

  if (Py !== 0n) {
    return 1n;
  }

  // P.y = 0: return 1/(3*x^2 + a4)
  const denom = mod(3n * mod(Px * Px, p) + a4, p);
  return Fp_inv(denom, p);
}

/**
 * Evaluate the line through R with given slope at point Q.
 *
 * The line is: y - y_R = slope * (x - x_R)
 * Evaluating at Q: y_Q - y_R - slope * (x_Q - x_R)
 *                = y_Q - (slope * (x_Q - x_R) + y_R)
 *
 * Reference: FpE.c:464-483 (FpE_Miller_line)
 */
function FpE_Miller_line(
  R: EllipticPointFp,
  Q: EllipticPointFp,
  slope: bigint,
  a4: bigint,
  p: bigint
): bigint {
  const Qx = Q.x!;
  const Qy = Q.y!;
  const Rx = R.x!;
  const Ry = R.y!;

  // tmp1 = Qx - Rx
  const tmp1 = mod(Qx - Rx, p);
  // tmp2 = slope * tmp1 + Ry
  const tmp2 = mod(mod(slope * tmp1, p) + Ry, p);

  // If Qy != tmp2, return Qy - tmp2
  if (Qy !== tmp2) {
    const result = mod(Qy - tmp2, p);
    return result < 0n ? result + p : result;
  }

  // Q is on the line, need special handling
  if (Qy === 0n) {
    return 1n;
  }

  // Compute tangent slope at Q and check if it equals the line slope
  const y2i = Fp_inv(mod(2n * Qy, p), p);
  const Qx2 = mod(Qx * Qx, p);
  const s1 = mod(mod(3n * Qx2 + a4, p) * y2i, p);

  if (s1 !== slope) {
    const result = mod(s1 - slope, p);
    return result < 0n ? result + p : result;
  }

  // Double tangency
  const s2 = mod(mod(3n * Qx - mod(s1 * s1, p), p) * y2i, p);
  return s2 !== 0n ? s2 : y2i;
}

/**
 * Evaluate the tangent line at R at point Q, and update R to 2*R.
 *
 * Reference: FpE.c:490-506 (FpE_tangent_update)
 */
function FpE_tangent_update(
  R: EllipticPointFp,
  Q: EllipticPointFp,
  a4: bigint,
  p: bigint
): { line: bigint; R: EllipticPointFp } {
  if (ell_is_inf(R)) {
    return { line: 1n, R: ellinf() };
  }

  if (R.y === 0n) {
    return { line: FpE_vert(R, Q, a4, p), R: ellinf() };
  }

  const { point: newR, slope } = FpE_dbl_slope(R, a4, p);
  const line = FpE_Miller_line(R, Q, slope!, a4, p);
  return { line, R: newR };
}

/**
 * Evaluate the chord through R and P at point Q, and update R to R + P.
 *
 * Reference: FpE.c:513-538 (FpE_chord_update)
 */
function FpE_chord_update(
  R: EllipticPointFp,
  P: EllipticPointFp,
  Q: EllipticPointFp,
  a4: bigint,
  p: bigint
): { line: bigint; R: EllipticPointFp } {
  if (ell_is_inf(R)) {
    return { line: FpE_vert(P, Q, a4, p), R: P.isInfinity ? ellinf() : ellpoint(P.x!, P.y!) };
  }

  if (ell_is_inf(P)) {
    return { line: FpE_vert(R, Q, a4, p), R: ellpoint(R.x!, R.y!) };
  }

  if (P.x === R.x) {
    if (P.y === R.y) {
      // Same point, use tangent
      return FpE_tangent_update(R, Q, a4, p);
    }
    // Inverse points
    return { line: FpE_vert(R, Q, a4, p), R: ellinf() };
  }

  const { point: newR, slope } = FpE_add_slope(P, R, a4, p);
  const line = FpE_Miller_line(R, Q, slope!, a4, p);
  return { line, R: newR };
}

// ============================================================================
// Pairings
// ============================================================================

/**
 * Compute the Tate pairing of points P and Q.
 *
 * The Tate pairing is a bilinear map e: E[m] x E(Fq)/mE(Fq) -> Fq* / (Fq*)^m
 * For m-torsion points P in E[m] and any point Q.
 *
 * Uses Miller's algorithm for efficient computation.
 *
 * Note: This implementation returns the "raw" Miller function value f_P(Q).
 * For the full Tate pairing with final exponentiation, the result should be
 * raised to the power (q^k - 1) / m where k is the embedding degree.
 *
 * @param E - Elliptic curve over Fp
 * @param P - First point (m-torsion point)
 * @param Q - Second point
 * @param m - Torsion order (P must satisfy [m]P = O)
 * @returns Element of Fp* representing the Tate pairing value
 *
 * @see Reference: pari/src/basemath/elliptic.c:6282-6299
 * @see Reference: pari/src/basemath/FpE.c:609-621 (FpE_tatepairing)
 *
 * @example
 * ```typescript
 * const E = { a4: 0n, a6: 1n, p: 101n };
 * const P = { x: 2n, y: 3n, isInfinity: false };
 * const Q = { x: 5n, y: 7n, isInfinity: false };
 * const m = 17n; // assuming [17]P = O
 * const result = elltatepairing(E, P, Q, m);
 * ```
 */
export function elltatepairing(
  E: EllipticCurveFp,
  P: EllipticPointFp,
  Q: EllipticPointFp,
  m: bigint
): bigint {
  // Handle trivial cases
  if (ell_is_inf(P) || ell_is_inf(Q)) {
    return 1n;
  }

  return _FpE_Miller(P, Q, m, E);
}

/**
 * Compute the Weil pairing of points P and Q.
 *
 * The Weil pairing is a bilinear, alternating, non-degenerate map
 * e: E[m] x E[m] -> mu_m (m-th roots of unity)
 *
 * For two m-torsion points P and Q, the Weil pairing satisfies:
 * - e(P, P) = 1 (alternating)
 * - e(P, Q) * e(Q, P) = 1 (antisymmetric)
 * - e(P, Q)^m = 1 (values are m-th roots of unity)
 *
 * Computed as: e(P, Q) = (-1)^m * f_P(Q) / f_Q(P) where f_P, f_Q are Miller functions.
 *
 * @param E - Elliptic curve over Fp (or extension field)
 * @param P - First m-torsion point
 * @param Q - Second m-torsion point
 * @param m - Torsion order ([m]P = [m]Q = O)
 * @returns Element of Fp* (m-th root of unity)
 *
 * @see Reference: pari/src/basemath/elliptic.c:6263-6279
 * @see Reference: pari/src/basemath/FpE.c:589-606 (FpE_weilpairing)
 *
 * @example
 * ```typescript
 * const E = { a4: 0n, a6: 1n, p: 101n };
 * const P = { x: 2n, y: 3n, isInfinity: false };
 * const Q = { x: 5n, y: 7n, isInfinity: false };
 * const m = 17n; // order of torsion subgroup
 * const zeta = ellweilpairing(E, P, Q, m); // m-th root of unity
 * ```
 */
export function ellweilpairing(
  E: EllipticCurveFp,
  P: EllipticPointFp,
  Q: EllipticPointFp,
  m: bigint
): bigint {
  const { p } = E;

  // Handle trivial cases: e(O, Q) = e(P, O) = e(P, P) = 1
  if (ell_is_inf(P) || ell_is_inf(Q) || ellequal(P, Q)) {
    return 1n;
  }

  // Compute f_P(Q) and f_Q(P)
  const N = _FpE_Miller(P, Q, m, E);
  const D = _FpE_Miller(Q, P, m, E);

  // w = f_P(Q) / f_Q(P)
  let w = mod(N * Fp_inv(D, p), p);

  // If m is odd, multiply by -1
  // This comes from the definition: e(P,Q) = (-1)^m * f_P(Q) / f_Q(P)
  if ((m & 1n) === 1n) {
    w = mod(-w, p);
    if (w < 0n) w += p;
  }

  return w;
}

// ============================================================================
// SEA Algorithm (Schoof-Elkies-Atkin)
// ============================================================================

/**
 * Compute the cardinality of E(Fp) using the SEA algorithm.
 *
 * The Schoof-Elkies-Atkin (SEA) algorithm is an efficient method for
 * computing the number of points on an elliptic curve over a prime field.
 * It improves on Schoof's algorithm by using modular polynomials to
 * distinguish "Elkies primes" (where the modular polynomial factors)
 * from "Atkin primes".
 *
 * Time complexity: O((log p)^4) with fast arithmetic
 *
 * This is typically used for large primes where baby-step giant-step
 * becomes impractical (roughly p > 10^9).
 *
 * @param E - Elliptic curve over Fp in short Weierstrass form
 * @param smallfact - Optional: bound for small factor base (affects performance)
 * @returns Number of points #E(Fp) = p + 1 - t where t is the trace of Frobenius
 *
 * @see Reference: pari/src/basemath/ellsea.c:2112-2113 (Fp_ellcard_SEA)
 * @see Reference: pari/src/basemath/elliptic.c:6332-6353 (ellsea)
 *
 * @example
 * ```typescript
 * // Large prime example
 * const p = 2n ** 256n - 2n ** 32n - 977n; // secp256k1 prime
 * const E = { a4: 0n, a6: 7n, p };
 * const card = ellcard_sea(E);
 * ```
 */
export function ellcard_sea(
  E: EllipticCurveFp,
  smallfact?: number
): bigint {
  // Parameters intentionally unused in stub
  void E;
  void smallfact;
  throw new Error('PARI_NOT_IMPLEMENTED: ellcard_sea - Schoof-Elkies-Atkin point counting');
}

// ============================================================================
// Isogenies
// ============================================================================

/**
 * Compute an isogeny with a given kernel.
 *
 * Given an elliptic curve E and a subgroup G (specified either as a
 * generator point or as a polynomial whose roots are x-coordinates of
 * kernel points), computes the quotient curve E/G and the isogeny
 * phi: E -> E/G.
 *
 * Uses Velu's formulas for computing isogenies from kernel points.
 *
 * @param E - Source elliptic curve
 * @param G - Kernel specification: either a point generating the kernel,
 *            or a polynomial whose roots are x-coordinates of kernel points
 * @param onlyImage - If true, only return the codomain curve (not the isogeny maps)
 * @returns If onlyImage is true: the codomain curve E/G
 *          Otherwise: [codomain, isogeny] where isogeny can be applied via ellisogenyapply
 *
 * @see Reference: pari/src/basemath/ellisog.c:631-669 (ellisogeny)
 *
 * @example
 * ```typescript
 * // Compute 3-isogeny
 * const E = { a4: 1n, a6: 1n, p: 101n };
 * const P = { x: 5n, y: 10n, isInfinity: false }; // point of order 3
 * const [E2, phi] = ellisogeny(E, P);
 * // E2 is the codomain curve
 * // phi can be used with ellisogenyapply
 * ```
 */
export function ellisogeny(
  E: EllipticCurveFp,
  G: EllipticPointFp | bigint[],
  onlyImage?: boolean
): EllipticCurveFp | [EllipticCurveFp, IsogenyMap] {
  // Parameters intentionally unused in stub
  void E;
  void G;
  void onlyImage;
  throw new Error('PARI_NOT_IMPLEMENTED: ellisogeny - Velu\'s isogeny formulas');
}

/**
 * Apply an isogeny to a point.
 *
 * Given an isogeny phi (as returned by ellisogeny) and a point P
 * on the domain curve, computes phi(P) on the codomain curve.
 *
 * The isogeny is represented as rational functions:
 * - x' = f(x) / h(x)^2
 * - y' = g(x, y) / h(x)^3
 *
 * @param phi - Isogeny map from ellisogeny
 * @param P - Point on the domain curve
 * @returns Image point phi(P) on the codomain curve
 *
 * @see Reference: pari/src/basemath/ellisog.c:187-213 (ellisogenyapply)
 *
 * @example
 * ```typescript
 * const E = { a4: 1n, a6: 1n, p: 101n };
 * const kernel = { x: 5n, y: 10n, isInfinity: false };
 * const [E2, phi] = ellisogeny(E, kernel);
 * const P = { x: 20n, y: 30n, isInfinity: false };
 * const Q = ellisogenyapply(phi, P); // Q is on E2
 * ```
 */
export function ellisogenyapply(
  phi: IsogenyMap,
  P: EllipticPointFp
): EllipticPointFp {
  // Parameters intentionally unused in stub
  void phi;
  void P;
  throw new Error('PARI_NOT_IMPLEMENTED: ellisogenyapply - Apply isogeny to point');
}

/**
 * Compose two isogenies.
 *
 * Given isogenies phi: E1 -> E2 and psi: E2 -> E3,
 * computes the composed isogeny psi o phi: E1 -> E3.
 *
 * @param psi - Second isogeny (applied second)
 * @param phi - First isogeny (applied first)
 * @returns Composed isogeny psi o phi
 *
 * @see Reference: pari/src/basemath/ellisog.c:192 (ellcompisog branch)
 */
export function ellisogenycompose(
  psi: IsogenyMap,
  phi: IsogenyMap
): IsogenyMap {
  // Parameters intentionally unused in stub
  void psi;
  void phi;
  throw new Error('PARI_NOT_IMPLEMENTED: ellisogenycompose - Compose isogenies');
}

// ============================================================================
// Frobenius Endomorphism
// ============================================================================

/**
 * Compute the n-th power of the Frobenius endomorphism on a point.
 *
 * For a curve E over Fq where q = p^k, the Frobenius endomorphism
 * pi: E -> E is defined by pi(x, y) = (x^p, y^p).
 *
 * This function computes pi^n(P) = (x^(p^n), y^(p^n)).
 *
 * For curves over Fp (prime field), Frobenius is the identity,
 * so this is mainly useful for extension fields Fq.
 *
 * @param E - Elliptic curve over Fq
 * @param P - Point on E
 * @param n - Power of Frobenius (default: 1)
 * @returns The point pi^n(P) = (x^(q^n), y^(q^n))
 *
 * @example
 * ```typescript
 * // For curves over Fp, Frobenius is identity
 * const E = { a4: 0n, a6: 7n, p: 101n };
 * const P = { x: 5n, y: 10n, isInfinity: false };
 * const piP = ellfrobenius(E, P, 1); // equals P for prime fields
 * ```
 */
export function ellfrobenius(
  E: EllipticCurveFp,
  P: EllipticPointFp,
  n?: number
): EllipticPointFp {
  // Parameters intentionally unused in stub
  void E;
  void P;
  void n;
  throw new Error('PARI_NOT_IMPLEMENTED: ellfrobenius - Frobenius endomorphism');
}

// ============================================================================
// Miller's Algorithm (Internal helper for pairings)
// ============================================================================

/**
 * Miller's algorithm for computing the Miller function.
 *
 * This is an internal function used by both Tate and Weil pairings.
 * Computes f_P(Q) where f_P is the Miller function associated with
 * point P and evaluated at Q.
 *
 * The Miller function f_P satisfies:
 * div(f_P) = m(P) - m(O) - ([m]P) + (O)
 *
 * Algorithm: Double-and-add method following PARI/GP's gen_pow_i pattern.
 * For each bit of m:
 * - Double: Update accumulator by squaring and multiplying by tangent line / vertical
 * - Add: If bit is 1, multiply by chord line / vertical
 *
 * @param P - Base point for Miller function
 * @param Q - Evaluation point
 * @param m - Order (P is m-torsion)
 * @param E - Elliptic curve
 * @returns f_P(Q) as an element of Fp*
 *
 * @see Reference: pari/src/basemath/FpE.c:565-586 (FpE_Miller)
 * @internal
 */
export function _FpE_Miller(
  P: EllipticPointFp,
  Q: EllipticPointFp,
  m: bigint,
  E: EllipticCurveFp
): bigint {
  const { a4, p } = E;

  // Handle trivial cases
  if (m === 0n) {
    return 1n;
  }

  if (ell_is_inf(P) || ell_is_inf(Q)) {
    return 1n;
  }

  // Get the binary representation of m
  // We use a double-and-add algorithm similar to PARI's gen_pow_i
  // State: [N, D, point] where the result is N/D and point is the running sum
  let N = 1n;
  let D = 1n;
  let R = P.isInfinity ? ellinf() : ellpoint(P.x!, P.y!);

  // Process bits from second-highest to lowest
  const bits: boolean[] = [];
  let temp = m;
  while (temp > 0n) {
    bits.push((temp & 1n) === 1n);
    temp >>= 1n;
  }
  bits.reverse();

  // Skip the highest bit (always 1 for m > 0)
  for (let i = 1; i < bits.length; i++) {
    // Double step
    N = mod(N * N, p);
    D = mod(D * D, p);

    const tangent = FpE_tangent_update(R, Q, a4, p);
    N = mod(N * tangent.line, p);
    R = tangent.R;

    const v = FpE_vert(R, Q, a4, p);
    D = mod(D * v, p);

    // Add step if bit is 1
    if (bits[i]) {
      const chord = FpE_chord_update(R, P, Q, a4, p);
      N = mod(N * chord.line, p);
      R = chord.R;

      const v2 = FpE_vert(R, Q, a4, p);
      D = mod(D * v2, p);
    }
  }

  return mod(N * Fp_inv(D, p), p);
}

// ============================================================================
// Embedding Degree
// ============================================================================

/**
 * Compute the embedding degree of a curve.
 *
 * The embedding degree k is the smallest positive integer such that
 * the m-th roots of unity are contained in Fq^k, i.e., m | q^k - 1.
 *
 * This is important for pairing-based cryptography as it determines
 * the extension field where the pairing takes values.
 *
 * @param E - Elliptic curve
 * @param m - Subgroup order (typically the largest prime factor of #E(Fq))
 * @returns The embedding degree k
 *
 * @example
 * ```typescript
 * const E = { a4: 0n, a6: 1n, p: 101n };
 * const m = 17n; // prime dividing #E(Fp)
 * const k = ellembeddingdegree(E, m); // embedding degree
 * ```
 */
export function ellembeddingdegree(
  E: EllipticCurveFp,
  m: bigint
): number {
  const { p } = E;

  // Find smallest k such that m | p^k - 1
  // This is equivalent to p^k = 1 (mod m)

  // Special case: if m = 1, any k works
  if (m === 1n) {
    return 1;
  }

  // Compute p mod m
  let pk = mod(p, m);

  // Find the multiplicative order of p in (Z/mZ)*
  for (let k = 1; k <= 1000; k++) {
    if (pk === 1n) {
      return k;
    }
    pk = mod(pk * p, m);
  }

  throw new Error(`Embedding degree exceeds 1000 for m = ${m}`);
}

// ============================================================================
// Division Polynomials
// ============================================================================

/**
 * Polynomial arithmetic over Fp.
 * Polynomials are represented as arrays where poly[i] is the coefficient of x^i.
 * Empty array represents the zero polynomial.
 */

/**
 * Normalize a polynomial by removing trailing zeros.
 */
function polyNormalize(poly: bigint[], p: bigint): bigint[] {
  // Reduce all coefficients mod p
  const result = poly.map((c) => mod(c, p));

  // Remove trailing zeros
  let len = result.length;
  while (len > 0 && result[len - 1] === 0n) {
    len--;
  }
  return result.slice(0, len);
}

/**
 * Return the zero polynomial.
 */
function polyZero(): bigint[] {
  return [];
}

/**
 * Return the constant polynomial 1.
 */
function polyOne(): bigint[] {
  return [1n];
}

/**
 * Return the polynomial x.
 */
function polyX(): bigint[] {
  return [0n, 1n];
}

/**
 * Subtract two polynomials mod p.
 */
function polySub(a: bigint[], b: bigint[], p: bigint): bigint[] {
  const maxLen = Math.max(a.length, b.length);
  const result: bigint[] = [];

  for (let i = 0; i < maxLen; i++) {
    const ai = i < a.length ? a[i]! : 0n;
    const bi = i < b.length ? b[i]! : 0n;
    result.push(mod(ai - bi, p));
  }

  return polyNormalize(result, p);
}

/**
 * Multiply two polynomials mod p.
 */
function polyMul(a: bigint[], b: bigint[], p: bigint): bigint[] {
  if (a.length === 0 || b.length === 0) {
    return [];
  }

  const resultLen = a.length + b.length - 1;
  const result: bigint[] = new Array(resultLen).fill(0n);

  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < b.length; j++) {
      result[i + j] = mod(result[i + j]! + a[i]! * b[j]!, p);
    }
  }

  return polyNormalize(result, p);
}

/**
 * Multiply polynomial by a scalar mod p.
 */
function polyScalarMul(poly: bigint[], c: bigint, p: bigint): bigint[] {
  const cm = mod(c, p);
  if (cm === 0n) return [];
  return polyNormalize(poly.map((coeff) => mod(coeff * cm, p)), p);
}

/**
 * Square a polynomial mod p.
 */
function polySqr(a: bigint[], p: bigint): bigint[] {
  return polyMul(a, a, p);
}

/**
 * Compute a^3 for a polynomial mod p.
 */
function polyCube(a: bigint[], p: bigint): bigint[] {
  return polyMul(polyMul(a, a, p), a, p);
}

/**
 * Shift polynomial by n (multiply by x^n).
 * polyShift([a, b, c], 2) = [0, 0, a, b, c] = (a + bx + cx^2) * x^2
 */
function polyShift(poly: bigint[], n: number): bigint[] {
  if (poly.length === 0 || n === 0) return poly;
  const result: bigint[] = new Array(n).fill(0n);
  return result.concat(poly);
}

/**
 * Compute the b-model polynomial: 4x^3 + b2*x^2 + 2*b4*x + b6
 * For short Weierstrass y^2 = x^3 + a4*x + a6:
 *   b2 = 0, b4 = 2*a4, b6 = 4*a6
 * So: 4x^3 + 0*x^2 + 4*a4*x + 4*a6 = 4*(x^3 + a4*x + a6)
 *
 * This represents (2y)^2 = 4y^2 on the curve.
 */
function ecBModel(E: EllipticCurveFp): bigint[] {
  const { a4, a6, p } = E;
  // 4x^3 + 4*a4*x + 4*a6
  return polyNormalize([mod(4n * a6, p), mod(4n * a4, p), 0n, 4n], p);
}

/**
 * Compute phi_2: x^4 - b4*x^2 - 2*b6*x - b8
 * For short Weierstrass y^2 = x^3 + a4*x + a6:
 *   b4 = 2*a4, b6 = 4*a6, b8 = -a4^2
 * So: x^4 - 2*a4*x^2 - 8*a6*x + a4^2
 */
function ecPhi2(E: EllipticCurveFp): bigint[] {
  const { a4, a6, p } = E;
  const a4sq = mod(a4 * a4, p);
  // x^4 - 2*a4*x^2 - 8*a6*x + a4^2
  return polyNormalize([a4sq, mod(-8n * a6, p), mod(-2n * a4, p), 0n, 1n], p);
}

/**
 * Compute base division polynomials for n <= 4.
 *
 * For short Weierstrass form y^2 = x^3 + a4*x + a6:
 * - psi_0 = 0
 * - psi_1 = 1
 * - psi_2 = 2y (we store 1, multiply by d2 = 4y^2 later for even indices)
 * - psi_3 = 3x^4 + 6*a4*x^2 + 12*a6*x - a4^2
 * - psi_4 = 4y*(x^6 + 5*a4*x^4 + 20*a6*x^3 - 5*a4^2*x^2 - 4*a4*a6*x - 8*a6^2 - a4^3)
 *
 * We return "normalized" versions:
 * - For odd n: psi_n directly (polynomial in x)
 * - For even n: psi_n / (2y) (polynomial in x, will be multiplied by d2 = 4y^2 later)
 *
 * PARI's approach: store f_n where:
 * - f_n = psi_n for odd n
 * - f_n = psi_n / 2y for even n (i.e., divide out one factor of 2y)
 */
function elldivpol4(E: EllipticCurveFp, n: number): bigint[] {
  const { a4, a6, p } = E;

  if (n === 0) return polyZero();
  if (n === 1) return polyOne();
  if (n === 2) return polyOne(); // psi_2 / 2y = 1

  if (n === 3) {
    // psi_3 = 3x^4 + 6*a4*x^2 + 12*a6*x - a4^2
    const a4sq = mod(a4 * a4, p);
    return polyNormalize(
      [mod(-a4sq, p), mod(12n * a6, p), mod(6n * a4, p), 0n, 3n],
      p
    );
  }

  if (n === 4) {
    // psi_4 / 2y = 2*(x^6 + 5*a4*x^4 + 20*a6*x^3 - 5*a4^2*x^2 - 4*a4*a6*x - 8*a6^2 - a4^3)
    const a4sq = mod(a4 * a4, p);
    const a4cube = mod(a4sq * a4, p);
    const a6sq = mod(a6 * a6, p);
    return polyNormalize(
      [
        mod(-2n * (8n * a6sq + a4cube), p), // constant term
        mod(-2n * 4n * a4 * a6, p), // x^1
        mod(-2n * 5n * a4sq, p), // x^2
        mod(2n * 20n * a6, p), // x^3
        mod(2n * 5n * a4, p), // x^4
        0n, // x^5
        2n, // x^6
      ],
      p
    );
  }

  throw new Error(`elldivpol4: n must be <= 4, got ${n}`);
}

/**
 * Memoized recursive computation of division polynomials.
 *
 * Uses the recurrence relations:
 * - psi_{2m+1} = psi_{m+2}*psi_m^3 - psi_{m-1}*psi_{m+1}^3
 * - psi_{2m} = (psi_m / 2y) * (psi_{m+2}*psi_{m-1}^2 - psi_{m-2}*psi_{m+1}^2)
 *
 * PARI normalizes by powers of d2 = (2y)^2 = 4y^2:
 * - For odd n, f_n = psi_n
 * - For even n, f_n = psi_n / 2y
 *
 * The recurrence in terms of f_n (with T = d2^2 = 16y^4):
 * - f_{2m+1} when m is odd: f_{m+2}*f_m^3 - T*f_{m-1}*f_{m+1}^3
 * - f_{2m+1} when m is even: T*f_{m+2}*f_m^3 - f_{m-1}*f_{m+1}^3
 * - f_{2m} = f_m * (f_{m+2}*f_{m-1}^2 - f_{m-2}*f_{m+1}^2)
 */
function elldivpol0(
  E: EllipticCurveFp,
  cache: Map<number, bigint[]>,
  T: bigint[], // d2^2 = (4y^2)^2 = 16y^4 = 16(x^3+a4*x+a6)^2
  n: number
): bigint[] {
  const { p } = E;

  // Check cache
  const cached = cache.get(n);
  if (cached !== undefined) return cached;

  // Base cases
  if (n <= 4) {
    const result = elldivpol4(E, n);
    cache.set(n, result);
    return result;
  }

  const m = Math.floor(n / 2);
  let result: bigint[];

  if (n % 2 === 1) {
    // Odd case: f_{2m+1}
    const fm2 = elldivpol0(E, cache, T, m + 2);
    const fm = elldivpol0(E, cache, T, m);
    const fm1minus = elldivpol0(E, cache, T, m - 1);
    const fm1plus = elldivpol0(E, cache, T, m + 1);

    const t1 = polyMul(fm2, polyCube(fm, p), p);
    const t2 = polyMul(fm1minus, polyCube(fm1plus, p), p);

    if (m % 2 === 1) {
      // f_{4l+3} = f_{2l+3}*f_{2l+1}^3 - T*f_{2l}*f_{2l+2}^3, where m = 2l+1
      result = polySub(t1, polyMul(T, t2, p), p);
    } else {
      // f_{4l+1} = T*f_{2l+2}*f_{2l}^3 - f_{2l-1}*f_{2l+1}^3, where m = 2l
      result = polySub(polyMul(T, t1, p), t2, p);
    }
  } else {
    // Even case: f_{2m} = f_m * (f_{m+2}*f_{m-1}^2 - f_{m-2}*f_{m+1}^2)
    const fm = elldivpol0(E, cache, T, m);
    const fm2plus = elldivpol0(E, cache, T, m + 2);
    const fm1minus = elldivpol0(E, cache, T, m - 1);
    const fm2minus = elldivpol0(E, cache, T, m - 2);
    const fm1plus = elldivpol0(E, cache, T, m + 1);

    const t1 = polyMul(fm2plus, polySqr(fm1minus, p), p);
    const t2 = polyMul(fm2minus, polySqr(fm1plus, p), p);
    result = polyMul(fm, polySub(t1, t2, p), p);
  }

  cache.set(n, result);
  return result;
}

/**
 * Compute the n-th division polynomial of an elliptic curve.
 *
 * The n-th division polynomial psi_n is defined recursively and
 * has the property that its roots are the x-coordinates of the
 * non-trivial n-torsion points.
 *
 * For a short Weierstrass curve y^2 = x^3 + a4*x + a6:
 * - For odd n: returns psi_n(x), a polynomial in x only
 * - For even n: returns psi_n(x) / (2y), a polynomial in x only
 *
 * The degree is:
 * - For odd n: (n^2 - 1) / 2
 * - For even n: (n^2 - 4) / 2 (after dividing by 2y)
 *
 * @param E - Elliptic curve in short Weierstrass form over Fp
 * @param n - Division polynomial index (must be non-negative)
 * @returns Coefficients of the division polynomial as an array where result[i] is the coefficient of x^i
 *
 * @see Reference: pari/src/basemath/elliptic.c:6757-6782 (elldivpol)
 *
 * @example
 * ```typescript
 * const E = { a4: 1n, a6: 1n, p: 23n };
 * const psi3 = elldivpol(E, 3);  // 3x^4 + 6x^2 + 12x - 1
 * ```
 */
export function elldivpol(E: EllipticCurveFp, n: number): bigint[] {
  const { p } = E;
  const n0 = n;
  n = Math.abs(n);

  // Base cases that don't need d2
  if (n === 1 || n === 3) {
    const result = elldivpol4(E, n);
    return n0 < 0 ? polyScalarMul(result, -1n, p) : result;
  }

  // d2 = 4y^2 = 4*(x^3 + a4*x + a6) for short Weierstrass
  const d2 = ecBModel(E);

  if (n <= 4) {
    let f = elldivpol4(E, n);
    // For even n, multiply by d2 to get the full psi_n
    // But we want to return psi_n / 2y for even n, so multiply by d2^(1/2) factor
    // Actually, for elldivpol PARI returns: for even n, f * d2
    if (n % 2 === 0) {
      f = polyMul(f, d2, p);
    }
    return n0 < 0 ? polyScalarMul(f, -1n, p) : f;
  }

  // For n > 4, use the recursive formula
  const T = polySqr(d2, p); // T = d2^2 = 16y^4
  const cache = new Map<number, bigint[]>();
  let f = elldivpol0(E, cache, T, n);

  // For even n, multiply by d2
  if (n % 2 === 0) {
    f = polyMul(f, d2, p);
  }

  return n0 < 0 ? polyScalarMul(f, -1n, p) : f;
}

/**
 * Compute [phi_n, psi_n^2] for the x-coordinate of [n]P.
 *
 * The x-coordinate of [n]P can be expressed as:
 * x([n]P) = phi_n(x(P)) / psi_n(x(P))^2
 *
 * This function returns both polynomials as coefficient arrays.
 *
 * For short Weierstrass y^2 = x^3 + a4*x + a6:
 * - phi_n = x * psi_n^2 - psi_{n-1} * psi_{n+1}
 * - psi_n^2 is the square of the n-th division polynomial
 *
 * @param E - Elliptic curve in short Weierstrass form over Fp
 * @param n - Multiplication factor (non-negative)
 * @returns [phi_n, psi_n^2] as coefficient arrays
 *
 * @see Reference: pari/src/basemath/elliptic.c:6785-6828 (ellxn)
 *
 * @example
 * ```typescript
 * const E = { a4: 1n, a6: 1n, p: 23n };
 * const [phi2, psi2sq] = ellxn(E, 2);
 * // For a point P = (x, y), x([2]P) = phi2(x) / psi2sq(x)
 * ```
 */
export function ellxn(E: EllipticCurveFp, n: number): [bigint[], bigint[]] {
  const { p } = E;
  n = Math.abs(n);

  // d2 = 4y^2 = 4*(x^3 + a4*x + a6)
  const d2 = ecBModel(E);

  if (n === 0) {
    return [polyZero(), polyZero()];
  }

  if (n === 1) {
    // [1]P = P, so phi_1 = x, psi_1^2 = 1
    return [polyX(), polyOne()];
  }

  if (n === 2) {
    // psi_2 = 2y, so psi_2^2 = 4y^2 = d2
    // phi_2 = x^4 - 2*a4*x^2 - 8*a6*x + a4^2
    return [ecPhi2(E), d2];
  }

  // For n > 2, use the recursive division polynomials
  const T = polySqr(d2, p);
  const cache = new Map<number, bigint[]>();

  // f_n, f_{n-1}, f_{n+1}
  const fn = elldivpol0(E, cache, T, n);
  const fnMinus1 = elldivpol0(E, cache, T, n - 1);
  const fnPlus1 = elldivpol0(E, cache, T, n + 1);

  // Compute psi_n^2
  // For odd n: psi_n = f_n, so psi_n^2 = f_n^2
  // For even n: psi_n = f_n * 2y, so psi_n^2 = f_n^2 * 4y^2 = f_n^2 * d2
  let psiNSq = polySqr(fn, p);
  if (n % 2 === 0) {
    psiNSq = polyMul(psiNSq, d2, p);
  }

  // Compute u = psi_{n-1} * psi_{n+1}
  // psi_{n-1} * psi_{n+1}:
  // - If n is odd: n-1 even, n+1 even, so psi_{n-1}*psi_{n+1} = f_{n-1}*f_{n+1} * d2^2 = f_{n-1}*f_{n+1}*T
  // - If n is even: n-1 odd, n+1 odd, so psi_{n-1}*psi_{n+1} = f_{n-1}*f_{n+1}
  let u = polyMul(fnMinus1, fnPlus1, p);
  if (n % 2 === 1) {
    // n odd means n-1 and n+1 are even
    u = polyMul(u, d2, p);
  }

  // phi_n = x * psi_n^2 - psi_{n-1} * psi_{n+1}
  // x * psi_n^2 is polyShift(psiNSq, 1)
  const phiN = polySub(polyShift(psiNSq, 1), u, p);

  return [phiN, psiNSq];
}
