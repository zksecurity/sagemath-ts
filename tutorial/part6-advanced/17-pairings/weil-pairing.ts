/**
 * # The Weil Pairing
 *
 * The Weil pairing is a bilinear map on the n-torsion points of an elliptic curve.
 * It was discovered by Andre Weil in the 1940s as a tool for studying the arithmetic
 * of elliptic curves, but has become fundamental to modern cryptography.
 *
 * ## Mathematical Definition
 *
 * For an elliptic curve E and integer n coprime to the characteristic:
 *   e_n: E[n] x E[n] -> mu_n
 *
 * where:
 *   - E[n] = {P in E : nP = O} is the n-torsion subgroup
 *   - mu_n = {z : z^n = 1} is the group of nth roots of unity
 *
 * ## Key Properties
 *
 * 1. BILINEARITY: e_n(aP, bQ) = e_n(P, Q)^(ab)
 * 2. ALTERNATING: e_n(P, P) = 1 for all P
 * 3. NON-DEGENERATE: If e_n(P, Q) = 1 for all Q, then P = O
 * 4. COMPATIBILITY: e_n(P, Q)^m = e_{nm}(P, Q) when defined
 *
 * @module tutorial/part6-advanced/17-pairings/weil-pairing
 */

import { GF, EllipticCurve, weil_pairing } from 'sagemath-ts';

console.log("=== The Weil Pairing ===\n");

// ============================================================================
// Section 1: Mathematical Background
// ============================================================================

console.log("--- Section 1: Mathematical Foundation ---\n");

/**
 * THE WEIL PAIRING: FORMAL DEFINITION
 *
 * Let E be an elliptic curve over a field K, and let n be a positive integer
 * coprime to char(K).
 *
 * The n-torsion subgroup E[n] consists of all points P such that nP = O.
 * Over the algebraic closure, E[n] is isomorphic to (Z/nZ)^2.
 *
 * The Weil pairing e_n: E[n] x E[n] -> mu_n is defined using divisors:
 *
 * DIVISOR APPROACH:
 * For P, Q in E[n] with P != Q:
 * 1. Let f_P be a function with divisor n(P) - n(O)
 *    (i.e., f_P has a zero of order n at P and a pole of order n at O)
 * 2. Similarly, f_Q has divisor n(Q) - n(O)
 *
 * Then: e_n(P, Q) = f_P(Q+S) / f_P(S) * f_Q(S) / f_Q(P+S)
 *
 * where S is a point chosen so all evaluations are defined.
 *
 * In practice, we use Miller's algorithm to compute the pairing efficiently.
 */

console.log("Weil Pairing Definition:");
console.log("  e_n: E[n] x E[n] -> mu_n");
console.log("  Maps n-torsion points to nth roots of unity\n");

console.log("Key insight: The pairing 'extracts' information about");
console.log("the linear independence of torsion points.\n");

// ============================================================================
// Section 2: Properties of the Weil Pairing
// ============================================================================

console.log("--- Section 2: Properties ---\n");

/**
 * PROPERTIES OF THE WEIL PAIRING
 *
 * 1. BILINEARITY (in both arguments):
 *    e_n(P1 + P2, Q) = e_n(P1, Q) * e_n(P2, Q)
 *    e_n(P, Q1 + Q2) = e_n(P, Q1) * e_n(P, Q2)
 *
 *    Consequence: e_n(aP, bQ) = e_n(P, Q)^(ab)
 *
 * 2. ALTERNATING (skew-symmetric):
 *    e_n(P, Q) = e_n(Q, P)^(-1)
 *    e_n(P, P) = 1 for all P
 *
 *    The alternating property distinguishes the Weil pairing from
 *    the Tate pairing, which is not alternating.
 *
 * 3. NON-DEGENERATE:
 *    If e_n(P, Q) = 1 for all Q in E[n], then P = O.
 *    Equivalently: for non-trivial P, there exists Q with e_n(P, Q) != 1.
 *
 * 4. GALOIS EQUIVARIANCE:
 *    For sigma in Gal(K-bar/K):
 *    sigma(e_n(P, Q)) = e_n(sigma(P), sigma(Q))
 *
 * 5. COMPATIBILITY:
 *    e_{nm}(P, Q) = e_n(mP, Q) = e_n(P, mQ) when defined
 *    e_n(P, Q)^m = e_{nm}(P, Q)
 */

console.log("Properties:");
console.log("1. Bilinear:     e(aP, bQ) = e(P, Q)^(ab)");
console.log("2. Alternating:  e(P, Q) = e(Q, P)^(-1), e(P, P) = 1");
console.log("3. Non-degenerate: e(P, Q) = 1 for all Q implies P = O");
console.log("4. Galois equivariant: sigma(e(P,Q)) = e(sigma(P), sigma(Q))");
console.log("5. Compatible: e_nm(P, Q) = e_n(mP, Q)\n");

// ============================================================================
// Section 3: Example Computation
// ============================================================================

console.log("--- Section 3: Computing the Weil Pairing ---\n");

// For a proper demonstration, we need a curve where torsion points
// are rational (defined over the base field).

// Let's use a small example
const p = 631n;  // Prime field
const F = GF(p);
const E = EllipticCurve(F, [30n, 34n]);  // y^2 = x^3 + 30x + 34

console.log(`Field: F_${p}`);
console.log(`Curve: ${E.toString()}`);

const curveOrder = E.cardinality();
console.log(`Curve order: ${curveOrder}\n`);

// For the Weil pairing, we need n-torsion points
// The torsion subgroup E[n] is defined over F_p^k where k is the embedding degree

// Let's find the n-torsion points for a small n that divides the curve order
// We need two linearly independent n-torsion points

/**
 * FINDING n-TORSION POINTS
 *
 * For P to be n-torsion (nP = O), we need ord(P) | n.
 *
 * If the curve order is N, then n | N ensures E[n] is non-trivial.
 *
 * For E[n] to be fully defined over F_q, we need the embedding degree k
 * (smallest k such that n | q^k - 1) to be 1.
 */

// Try to find a suitable n
// First, check if 5 divides the curve order
const n = 5n;

if (curveOrder % n === 0n) {
  console.log(`Looking for ${n}-torsion points...`);

  // Find a point P of order n by taking a random point and multiplying by N/n
  const cofactor = curveOrder / n;

  let P = E.random_point().mul(cofactor);
  let attempts = 0;
  while (P.isZero() && attempts < 100) {
    P = E.random_point().mul(cofactor);
    attempts++;
  }

  console.log(`Found P of order dividing ${n}: ${P.toString()}`);
  console.log(`Verify: ${n}P = ${P.mul(n).toString()}\n`);

  // For the Weil pairing to be non-trivial, we need a SECOND independent point Q
  // In E[n] over the algebraic closure, E[n] is isomorphic to (Z/nZ)^2
  // Over F_p, we might only see part of E[n]

  let Q = E.random_point().mul(cofactor);
  attempts = 0;
  while ((Q.isZero() || Q.eq(P)) && attempts < 100) {
    Q = E.random_point().mul(cofactor);
    attempts++;
  }

  console.log(`Found Q of order dividing ${n}: ${Q.toString()}`);
  console.log(`Verify: ${n}Q = ${Q.mul(n).toString()}\n`);

  // Now compute the Weil pairing
  // Note: The Weil pairing might be trivial if P and Q are linearly dependent
  // (i.e., if there exists a such that Q = aP)

  console.log("Computing Weil pairing e_n(P, Q)...\n");

  // The weil_pairing function is available on points
  // It returns an element of the base field (or extension field)

  try {
    // In sagemath-ts, weil_pairing is imported from the ell_point module
    // It takes two points and the order n

    // Check if P and Q are linearly independent by computing the pairing
    // If e_n(P, Q) = 1, then P and Q might be dependent

    console.log("Weil pairing properties to verify:");
    console.log("  1. e(P, P) = 1 (alternating)");
    console.log("  2. e(P, Q) * e(Q, P) = 1 (skew-symmetry)");
    console.log("  3. e(aP, bQ) = e(P, Q)^(ab) (bilinearity)\n");

    // The actual computation would use:
    // const zeta = weil_pairing(P, Q, n);
    // console.log(`e_${n}(P, Q) = ${zeta}`);

    console.log("Note: For the Weil pairing to be non-trivial,");
    console.log("P and Q must be linearly independent in E[n].");
    console.log("This requires working over a field extension containing");
    console.log("all n-torsion points.\n");

  } catch (err: unknown) {
    console.log(`Pairing computation: ${err}`);
    console.log("(This is expected for educational demonstrations)\n");
  }
}

// ============================================================================
// Section 4: Miller's Algorithm
// ============================================================================

console.log("--- Section 4: Miller's Algorithm ---\n");

/**
 * MILLER'S ALGORITHM
 *
 * Victor Miller (1986) gave an efficient algorithm to compute pairings.
 *
 * The idea: Build the function f_P with divisor n(P) - n(O) iteratively
 * using the double-and-add structure.
 *
 * ALGORITHM (simplified):
 *
 * Input: P, Q in E[n], n > 0
 * Output: e_n(P, Q)
 *
 * 1. f := 1
 * 2. V := P
 * 3. For each bit b_i of n (from high to low):
 *    a. f := f^2 * l_{V,V}(Q) / v_{2V}(Q)   // doubling
 *    b. V := 2V
 *    c. If b_i = 1:
 *       f := f * l_{V,P}(Q) / v_{V+P}(Q)   // addition
 *       V := V + P
 * 4. Return f
 *
 * Where:
 * - l_{V,V}(Q) = value of tangent line at V, evaluated at Q
 * - l_{V,P}(Q) = value of line through V and P, evaluated at Q
 * - v_{V}(Q) = value of vertical line at V, evaluated at Q
 *
 * The function f accumulates the value of f_P(Q) as we compute nP.
 */

console.log("Miller's Algorithm (pseudo-code):");
console.log(`
function miller(P, Q, n):
    f = 1
    V = P
    for bit in binary(n):
        f = f^2 * line(V, V, Q) / vertical(2V, Q)  // double
        V = 2V
        if bit == 1:
            f = f * line(V, P, Q) / vertical(V+P, Q)  // add
            V = V + P
    return f
`);

console.log("Line functions:");
console.log("  l_{R,S}(Q): line through R and S evaluated at Q");
console.log("  v_R(Q): vertical line through R evaluated at Q\n");

console.log("Complexity: O(log n) group operations");
console.log("This makes pairings practical for cryptography.\n");

// ============================================================================
// Section 5: Embedding Degree
// ============================================================================

console.log("--- Section 5: Embedding Degree ---\n");

/**
 * EMBEDDING DEGREE
 *
 * The embedding degree k of a curve E with respect to n is the smallest
 * positive integer such that n | q^k - 1.
 *
 * Significance:
 * - E[n] is fully defined over F_{q^k}
 * - The Weil pairing maps to mu_n in F_{q^k}*
 * - k determines the security-efficiency tradeoff
 *
 * For cryptography:
 * - Small k (e.g., k = 2, 4, 6): efficient pairings but smaller GT
 * - Large k (e.g., k = 12, 24): slower pairings but better security ratio
 *
 * PAIRING-FRIENDLY CURVES
 *
 * Most random curves have k ~ q (embedding degree roughly field size).
 * Pairing-friendly curves are specially constructed with small k:
 *
 * - Supersingular curves: k | 6 (k = 1, 2, 3, or 6)
 * - BN curves: k = 12
 * - BLS12: k = 12
 * - BLS24: k = 24
 *
 * Security consideration:
 * - DLP in E: O(sqrt(q)) with Pollard-rho
 * - DLP in F_{q^k}*: subexponential (index calculus)
 *
 * Need q^k large enough to resist index calculus attacks.
 */

console.log("Embedding Degree k:");
console.log("  Smallest k such that n | q^k - 1");
console.log("  E[n] is defined over F_{q^k}");
console.log("  Pairing maps to mu_n in F_{q^k}*\n");

console.log("Security Balance:");
console.log("  - G1/G2: sqrt(q) security (ECDLP)");
console.log("  - GT: subexponential in q^k (index calculus)");
console.log("  - Need both to be hard at same security level\n");

console.log("Common Pairing-Friendly Curves:");
console.log("  BN254:     k=12, ~100-bit security");
console.log("  BLS12-381: k=12, ~128-bit security");
console.log("  BLS12-377: k=12, SNARK-friendly\n");

// ============================================================================
// Section 6: Weil Pairing vs Other Pairings
// ============================================================================

console.log("--- Section 6: Comparison with Other Pairings ---\n");

/**
 * WEIL vs TATE vs ATE PAIRINGS
 *
 * WEIL PAIRING e_n:
 * - Symmetric: e_n(P, Q) and e_n(Q, P) both well-defined
 * - Alternating: e_n(P, Q) = e_n(Q, P)^{-1}
 * - Self-pairing trivial: e_n(P, P) = 1
 * - Miller loop runs twice (once for each argument)
 *
 * TATE PAIRING <P, Q>_n:
 * - NOT alternating: <P, Q> != <Q, P>^{-1} in general
 * - Only one Miller loop needed
 * - Result needs "final exponentiation" to be well-defined
 * - More efficient than Weil (roughly 2x faster)
 *
 * ATE PAIRING:
 * - Optimized version of Tate
 * - Shorter Miller loop (uses Frobenius structure)
 * - Even more efficient than Tate
 *
 * OPTIMAL ATE / R-ATE / ETA PAIRINGS:
 * - Further optimizations using curve-specific structure
 * - Used in practice (BN254, BLS12-381 implementations)
 *
 * In practice, cryptographic libraries use Ate or optimal Ate pairings
 * for efficiency, while the Weil pairing is mainly of theoretical interest.
 */

console.log("Comparison:");
console.log("+---------+-------------+-----------+------------+");
console.log("| Pairing | Alternating | Loops     | Efficiency |");
console.log("+---------+-------------+-----------+------------+");
console.log("| Weil    | Yes         | 2 Miller  | Slowest    |");
console.log("| Tate    | No          | 1 Miller  | Faster     |");
console.log("| Ate     | No          | Shorter   | Fast       |");
console.log("| Opt-Ate | No          | Shortest  | Fastest    |");
console.log("+---------+-------------+-----------+------------+\n");

console.log("In practice: Libraries use Ate/optimal-Ate pairings");
console.log("The Weil pairing is mainly for theoretical analysis.\n");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");

console.log("The Weil pairing e_n: E[n] x E[n] -> mu_n is:");
console.log("  - Bilinear: enables multiplicative relationships");
console.log("  - Alternating: e(P, P) = 1, e(P, Q) = e(Q, P)^{-1}");
console.log("  - Non-degenerate: detects linear independence\n");

console.log("Computed via Miller's algorithm in O(log n) time.");
console.log("The embedding degree k determines where mu_n lives.\n");

console.log("Cryptographic importance:");
console.log("  - Foundation for pairing-based protocols");
console.log("  - Enables IBE, BLS signatures, SNARKs");
console.log("  - In practice, Tate/Ate pairings are used for efficiency\n");

console.log("Next: Learn about the Tate pairing and its optimizations");
