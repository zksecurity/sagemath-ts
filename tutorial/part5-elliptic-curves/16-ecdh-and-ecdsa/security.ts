/**
 * # Elliptic Curve Security Considerations
 *
 * The security of elliptic curve cryptography rests on the hardness of
 * the Elliptic Curve Discrete Logarithm Problem (ECDLP). This tutorial
 * explores what makes ECDLP hard and how curves can be weak.
 *
 * ## The ECDLP
 *
 * Given: E, P (base point), Q (target point where Q = nP for some n)
 * Find: n
 *
 * Best generic attack: O(sqrt(n)) using Pollard's rho algorithm.
 * For 256-bit n, this requires 2^128 operations - infeasible.
 *
 * ## Why Not Index Calculus?
 *
 * Unlike classical DLP in Z_p*, there is no known sub-exponential
 * index calculus attack on general elliptic curves. This is why ECC
 * can use smaller keys than RSA for equivalent security.
 *
 * @module tutorial/elliptic-curves/security
 */

import {
  EllipticCurve,
  is_supersingular,
  is_ordinary,
  discrete_log,
  embedding_degree,
} from 'sagemath-ts';
import { FiniteFieldPrime, factor } from 'sagemath-ts';

// ============================================================================
// The ECDLP Problem
// ============================================================================

console.log("=== The Elliptic Curve Discrete Logarithm Problem ===\n");

/**
 * ECDLP: Given P and Q = nP, find n.
 *
 * This is the foundation of ECC security.
 * - ECDH security: Can't find ab from aP and bP
 * - ECDSA security: Can't find d from Q = dG
 */

const F101 = new FiniteFieldPrime(101n);
const E = EllipticCurve(F101, [2n, 3n]);
const G = E.gens()[0]!;
const n = G.order();

console.log(`Curve: ${E.toString()}`);
console.log(`Generator G: ${G.toString()}`);
console.log(`Order n: ${n}`);
console.log("");

// Create a DLP instance
const secretK = 37n;
const Q = G.mul(secretK);

console.log("ECDLP Instance:");
console.log(`  Given: P = ${G.toString()}`);
console.log(`  Given: Q = ${Q.toString()}`);
console.log(`  Find: k such that Q = kP`);
console.log(`  (Secret answer: k = ${secretK})`);
console.log("");

// ============================================================================
// Generic Attacks
// ============================================================================

console.log("=== Generic Attacks on ECDLP ===\n");

/**
 * Generic attacks work on any group:
 *
 * 1. BRUTE FORCE: O(n)
 *    Try k = 1, 2, 3, ... until kP = Q
 *
 * 2. BABY-STEP GIANT-STEP: O(sqrt(n)) time and space
 *    - Baby steps: Store {iP : i = 0, 1, ..., m} where m = ceil(sqrt(n))
 *    - Giant steps: Check if Q - jmP matches any baby step
 *
 * 3. POLLARD'S RHO: O(sqrt(n)) time, O(1) space
 *    - Random walk in the group
 *    - Detect cycle using Floyd's algorithm
 *    - Most practical generic attack
 *
 * For 256-bit curves: sqrt(2^256) = 2^128 operations
 * At 10^12 ops/second, this takes 10^25 years
 */

console.log("1. BRUTE FORCE");
console.log("   Complexity: O(n)");
console.log("   For n = 2^256: 2^256 operations (infeasible)");
console.log("");

console.log("2. BABY-STEP GIANT-STEP (BSGS)");
console.log("   Complexity: O(sqrt(n)) time and space");
console.log("   For n = 2^256: 2^128 operations and storage");
console.log("");

console.log("3. POLLARD'S RHO");
console.log("   Complexity: O(sqrt(n)) time, O(1) space");
console.log("   Most practical for large groups");
console.log("   For n = 2^256: 2^128 operations");
console.log("");

// Demonstrate BSGS on small example
console.log("Demo: Solving ECDLP with baby-step giant-step:");

try {
  const startTime = Date.now();
  const recovered = discrete_log(G, Q, n);
  const endTime = Date.now();

  console.log(`  Recovered k = ${recovered}`);
  console.log(`  Correct? ${recovered === secretK}`);
  console.log(`  Time: ${endTime - startTime}ms`);
} catch (e) {
  console.log(`  Error: ${(e as Error).message}`);
}
console.log("");

// ============================================================================
// Pohlig-Hellman Attack
// ============================================================================

console.log("=== Pohlig-Hellman Attack ===\n");

/**
 * If the group order n has small prime factors, ECDLP is easier.
 *
 * Pohlig-Hellman:
 * 1. Factor n = p1^e1 * p2^e2 * ...
 * 2. Solve DLP in each prime-power subgroup
 * 3. Combine with CRT
 *
 * Complexity: O(sum(e_i * (sqrt(p_i) + log(n))))
 *
 * This is why we need LARGE PRIME ORDER subgroups!
 */

console.log("Pohlig-Hellman reduces DLP to subgroup DLPs.");
console.log("");

// Demonstrate with composite order
const F97 = new FiniteFieldPrime(97n);
const E_composite = EllipticCurve(F97, [1n, 1n]);
const order_composite = E_composite.cardinality();
const factorization = factor(order_composite);

console.log(`Curve with composite order over F_97:`);
console.log(`  Order: ${order_composite}`);
console.log(`  Factorization: ${factorization.map(([p, e]) => e > 1 ? `${p}^${e}` : `${p}`).join(' * ')}`);
console.log("");

// Security impact
const maxPrimeFactor = factorization.reduce(
  (max, [p]) => (p > max ? p : max),
  1n
);

console.log("Security analysis:");
console.log(`  Largest prime factor: ${maxPrimeFactor}`);
console.log(`  Effective security: ~${Math.floor(Math.log2(Number(maxPrimeFactor)) / 2)} bits`);
console.log("");

console.log("LESSON: Use curves where #E (or the subgroup order) is prime!");
console.log("");

// ============================================================================
// MOV Attack (Pairing-Based)
// ============================================================================

console.log("=== MOV Attack (for low embedding degree) ===\n");

/**
 * The MOV attack transfers ECDLP to DLP in a finite field extension.
 *
 * Using the Weil or Tate pairing:
 *   e(P, Q) is an element of F_{q^k}*
 *
 * If k (embedding degree) is small, we can:
 * 1. Compute e(G, Q) and e(G, G)^n in F_{q^k}
 * 2. Solve DLP in F_{q^k}* using index calculus
 *
 * Index calculus in F_{q^k}* has sub-exponential complexity!
 *
 * For security:
 * - k must be large enough that q^k is too big for index calculus
 * - Standard curves have astronomical k (safe)
 * - Pairing-friendly curves have moderate k by design
 */

console.log("MOV attack uses pairings to transfer ECDLP to F_{q^k}*");
console.log("");

// Check embedding degree for our curve
console.log("Checking embedding degree:");

const P = E.gens()[0]!;
const ord = P.order();

try {
  const k = embedding_degree(E, ord);
  console.log(`  Subgroup order: ${ord}`);
  console.log(`  Embedding degree k: ${k}`);

  if (k < 20n) {
    console.log(`  WARNING: Low embedding degree - potentially vulnerable to MOV!`);
  } else {
    console.log(`  Safe: Embedding degree is large`);
  }
} catch (e) {
  console.log(`  Could not compute embedding degree: ${(e as Error).message}`);
}
console.log("");

console.log("For cryptographic curves:");
console.log("  - secp256k1: k is astronomically large (~10^77)");
console.log("  - P-256: k is astronomically large");
console.log("  - BLS12-381: k = 12 (intentionally small for pairings)");
console.log("");

// ============================================================================
// Supersingular Curves
// ============================================================================

console.log("=== Supersingular vs Ordinary Curves ===\n");

/**
 * Supersingular curves have special properties:
 * - Trace t = 0 (mod p)
 * - Embedding degree k divides 6 (SMALL!)
 * - Endomorphism ring is non-commutative (quaternion algebra)
 *
 * For cryptographic ECDLP, avoid supersingular curves!
 * (But they're useful for isogeny-based crypto like SIDH/SIKE)
 */

// Check if our curve is supersingular
console.log(`Curve over F_101:`);
console.log(`  Is supersingular: ${is_supersingular(E)}`);
console.log(`  Is ordinary: ${is_ordinary(E)}`);
console.log("");

// Find and analyze a supersingular curve
console.log("Finding a supersingular curve over F_101...");

for (let a = 0n; a < 101n; a++) {
  for (let b = 0n; b < 101n; b++) {
    try {
      const E_test = EllipticCurve(F101, [a, b]);
      if (is_supersingular(E_test)) {
        const order = E_test.cardinality();
        console.log(`  Found: y^2 = x^3 + ${a}x + ${b}`);
        console.log(`  Order: ${order} (= p + 1 = 102 for supersingular)`);
        console.log(`  Trace: ${101n + 1n - order}`);
        console.log(`  DO NOT use for ECDLP-based crypto!`);
        console.log(`  (Embedding degree divides 6, vulnerable to MOV)`);
        break;
      }
    } catch {
      // Singular
    }
  }
}
console.log("");

// ============================================================================
// Anomalous Curves
// ============================================================================

console.log("=== Anomalous Curves (CRITICALLY WEAK) ===\n");

/**
 * An ANOMALOUS curve has #E(F_p) = p.
 *
 * Smart's attack (1999) solves ECDLP on anomalous curves in polynomial time!
 *
 * The attack:
 * 1. Lift curve and points to Q_p (p-adic numbers)
 * 2. Use the p-adic logarithm
 * 3. Reduces to linear algebra
 *
 * NEVER use anomalous curves for cryptography!
 */

console.log("Anomalous curves have #E = p (trace t = 1)");
console.log("Smart's attack solves ECDLP in POLYNOMIAL TIME!");
console.log("");

// Search for an anomalous curve
console.log("Finding an anomalous curve (DO NOT USE FOR CRYPTO):");

let anomalousFound = false;
for (const p of [7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n]) {
  if (anomalousFound) break;

  const Fp = new FiniteFieldPrime(p);

  for (let a = 0n; a < p && !anomalousFound; a++) {
    for (let b = 0n; b < p && !anomalousFound; b++) {
      try {
        const E_anom = EllipticCurve(Fp, [a, b]);
        if (E_anom.cardinality() === p) {
          console.log(`  Found over F_${p}: y^2 = x^3 + ${a}x + ${b}`);
          console.log(`  Order: ${p} = p (ANOMALOUS!)`);
          console.log(`  Trace: 1`);
          console.log(`  Smart's attack breaks this in O(log^3 p) time!`);
          anomalousFound = true;
        }
      } catch {
        // Singular
      }
    }
  }
}
console.log("");

// ============================================================================
// Invalid Curve Attacks
// ============================================================================

console.log("=== Invalid Curve Attacks ===\n");

/**
 * If implementations don't validate points, attackers can:
 *
 * 1. Send points NOT on the curve
 * 2. These points might be on a different curve with weak properties
 * 3. Recover private key through multiple queries
 *
 * Mitigation:
 * - Always verify received points are on the correct curve
 * - Check x^3 + ax + b is a quadratic residue
 * - Check point is in the correct subgroup
 */

console.log("Invalid curve attack: attacker sends points on wrong curve");
console.log("");

console.log("Example:");
console.log(`  Valid curve: y^2 = x^3 + 2x + 3 over F_101`);

// Create a point on a DIFFERENT curve
const E_invalid = EllipticCurve(F101, [2n, 4n]); // Different b
const P_invalid = E_invalid.lift_x(16n); // x=16 is valid on this curve

console.log(`  Invalid point (on y^2 = x^3 + 2x + 4): ${P_invalid.toString()}`);

// Check if this point is on the original curve
const x_inv = P_invalid.x!.value;
const y_inv = P_invalid.y!.value;
const lhs = (y_inv * y_inv) % 101n;
const rhs = (x_inv ** 3n + 2n * x_inv + 3n) % 101n;

console.log(`  On original curve? ${lhs === rhs} (LHS=${lhs}, RHS=${rhs})`);
console.log("");

console.log("Mitigation: ALWAYS validate received points!");
console.log("  1. Check (x, y) satisfies y^2 = x^3 + ax + b");
console.log("  2. Check point is not the identity");
console.log("  3. For cofactor > 1, check point order");
console.log("");

// ============================================================================
// Side-Channel Attacks
// ============================================================================

console.log("=== Side-Channel Attacks ===\n");

/**
 * Even with secure curves, implementations can leak information:
 *
 * 1. TIMING ATTACKS
 *    - Branching on secret data reveals information
 *    - Mitigation: Constant-time implementations
 *
 * 2. POWER ANALYSIS
 *    - Power consumption reveals operations
 *    - Mitigation: Balanced operations, masking
 *
 * 3. ELECTROMAGNETIC EMANATIONS
 *    - Similar to power analysis
 *    - Mitigation: Shielding, masking
 *
 * 4. CACHE ATTACKS
 *    - Memory access patterns reveal secrets
 *    - Mitigation: Constant-time table lookups
 */

console.log("Side-channel attacks exploit implementation details:");
console.log("");
console.log("1. TIMING ATTACKS");
console.log("   - Double-and-add leaks bit pattern of scalar");
console.log("   - Mitigation: Montgomery ladder, constant-time code");
console.log("");
console.log("2. POWER/EM ANALYSIS");
console.log("   - Different operations have different power signatures");
console.log("   - Mitigation: Unified addition formulas, masking");
console.log("");
console.log("3. FAULT ATTACKS");
console.log("   - Induce errors during computation");
console.log("   - Mitigation: Verification, redundancy");
console.log("");

// ============================================================================
// Security Level Summary
// ============================================================================

console.log("=== Security Level Summary ===\n");

console.log("Curve/Key Size | Security Level | Equivalent Symmetric");
console.log("---------------|----------------|---------------------");
console.log("160-bit curve  | 80 bits        | Not recommended");
console.log("224-bit curve  | 112 bits       | 3DES equivalent");
console.log("256-bit curve  | 128 bits       | AES-128 equivalent");
console.log("384-bit curve  | 192 bits       | AES-192 equivalent");
console.log("521-bit curve  | 256 bits       | AES-256 equivalent");
console.log("");

console.log("Minimum recommendations (2024):");
console.log("  - General use: 256-bit curves (128-bit security)");
console.log("  - Long-term security: 384-bit curves");
console.log("  - Post-quantum concerns: Consider hybrid approaches");
console.log("");

// ============================================================================
// Checklist for Secure ECC
// ============================================================================

console.log("=== ECC Security Checklist ===\n");

console.log("CURVE SELECTION:");
console.log("  [ ] Use standard, well-vetted curve");
console.log("  [ ] Prime order or small cofactor");
console.log("  [ ] Not anomalous (#E != p)");
console.log("  [ ] Not supersingular (for ECDLP)");
console.log("  [ ] Large embedding degree");
console.log("");

console.log("IMPLEMENTATION:");
console.log("  [ ] Constant-time scalar multiplication");
console.log("  [ ] Validate all received points");
console.log("  [ ] Clear cofactor when needed");
console.log("  [ ] Use cryptographic RNG for nonces");
console.log("  [ ] Never reuse nonces (ECDSA!)");
console.log("");

console.log("PROTOCOL:");
console.log("  [ ] Authenticate public keys");
console.log("  [ ] Use proven constructions (ECDH, ECDSA, etc.)");
console.log("  [ ] Apply key derivation functions");
console.log("  [ ] Consider forward secrecy (ephemeral keys)");
console.log("");

// ============================================================================
// Exercises
// ============================================================================

console.log("=== Exercises ===\n");

console.log("1. Implement the baby-step giant-step algorithm:");
console.log("   - Compare performance to brute force");
console.log("   - Verify it finds correct discrete log");
console.log("");

console.log("2. Research the MOV attack:");
console.log("   - How does the Weil pairing transfer ECDLP to DLP?");
console.log("   - Why is BLS12-381 safe despite k=12?");
console.log("");

console.log("3. Investigate Smart's attack on anomalous curves:");
console.log("   - What are p-adic numbers?");
console.log("   - Why does lifting to Q_p help?");
console.log("");

console.log("4. Find the embedding degree for:");
console.log("   - E: y^2 = x^3 + x + 1 over F_23 (with respect to a generator)");
console.log("   - Determine if MOV attack is practical");
console.log("");
