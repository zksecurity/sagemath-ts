/**
 * # Examples of Elliptic Curves
 *
 * This tutorial explores various elliptic curves and their properties,
 * building intuition for different curve types and their applications.
 *
 * ## Curve Families
 *
 * Different curves have different properties:
 * - **Supersingular curves**: Special structure, used in isogeny crypto
 * - **Ordinary curves**: Most common, used in ECDSA/ECDH
 * - **Curves with j = 0 or 1728**: Extra automorphisms
 * - **Anomalous curves**: #E = p (insecure for ECDLP!)
 *
 * @module tutorial/elliptic-curves/examples
 */

import {
  EllipticCurve,
  is_supersingular,
  is_ordinary,
} from 'sagemath-ts';
import { FiniteFieldPrime } from 'sagemath-ts';

// ============================================================================
// Example 1: Small Field Exploration
// ============================================================================

console.log("=== Example 1: Curves Over Small Fields ===\n");

// Over F_7, let's explore all possible curves
const F7 = new FiniteFieldPrime(7n);

console.log("All non-singular curves y^2 = x^3 + ax + b over F_7:");
console.log("(a, b)  -> order,  j-invariant");

for (let a = 0n; a < 7n; a++) {
  for (let b = 0n; b < 7n; b++) {
    // Check if curve is non-singular: 4a^3 + 27b^2 != 0
    const disc = (4n * a ** 3n + 27n * b ** 2n) % 7n;
    if (disc === 0n) continue;

    try {
      const E = EllipticCurve(F7, [a, b]);
      const order = E.cardinality();
      const j = E.j_invariant();
      console.log(`  (${a}, ${b}) -> order=${order}, j=${j}`);
    } catch {
      // Singular curve, skip
    }
  }
}
console.log("");

// ============================================================================
// Example 2: The Hasse Bound
// ============================================================================

console.log("=== Example 2: The Hasse Bound ===\n");

/**
 * Hasse's theorem: For an elliptic curve E over F_q,
 *
 *   |#E(F_q) - (q + 1)| <= 2*sqrt(q)
 *
 * This means:
 *   q + 1 - 2*sqrt(q) <= #E <= q + 1 + 2*sqrt(q)
 *
 * For F_97: 97 + 1 - 2*sqrt(97) = 78.3...
 *           97 + 1 + 2*sqrt(97) = 117.7...
 * So #E is between 79 and 117
 */

const F97 = new FiniteFieldPrime(97n);
const sqrtQ = Math.floor(Math.sqrt(97));
const lower = 97n + 1n - 2n * BigInt(sqrtQ) - 1n; // Round down
const upper = 97n + 1n + 2n * BigInt(sqrtQ) + 1n; // Round up

console.log(`Hasse bound for F_97: ${lower} <= #E <= ${upper}`);
console.log("");

// Try several curves and verify they satisfy Hasse bound
console.log("Curve orders over F_97:");
const orders = new Set<bigint>();

for (let a = 1n; a <= 10n; a++) {
  for (let b = 1n; b <= 10n; b++) {
    try {
      const E = EllipticCurve(F97, [a, b]);
      const order = E.cardinality();
      if (!orders.has(order)) {
        orders.add(order);
        const inBound = order >= lower && order <= upper;
        console.log(`  E(a=${a}, b=${b}): order=${order}, in Hasse bound: ${inBound}`);
      }
    } catch {
      // Singular
    }
  }
}
console.log("");

// ============================================================================
// Example 3: Supersingular vs Ordinary Curves
// ============================================================================

console.log("=== Example 3: Supersingular vs Ordinary Curves ===\n");

/**
 * A curve E over F_p is SUPERSINGULAR if p | trace(Frobenius)
 * Equivalently, if #E = p + 1 (trace = 0)
 *
 * Properties of supersingular curves:
 * - Embedding degree divides 6 (small!)
 * - j-invariant in F_p^2
 * - Used in isogeny-based cryptography (SIDH, SIKE)
 *
 * Most cryptographic curves are ORDINARY (not supersingular)
 */

// Find a supersingular curve over F_97
// Supersingular means #E = p + 1 = 98
console.log("Searching for supersingular curves over F_97...");
console.log("(Supersingular means #E = p + 1 = 98)");

let supersingularFound = false;
for (let a = 0n; a < 97n && !supersingularFound; a++) {
  for (let b = 0n; b < 97n && !supersingularFound; b++) {
    try {
      const E = EllipticCurve(F97, [a, b]);
      const order = E.cardinality();
      if (order === 98n) {
        console.log(`  Found: y^2 = x^3 + ${a}x + ${b}`);
        console.log(`  Order: ${order}`);
        console.log(`  Trace: ${97n + 1n - order}`);
        console.log(`  Supersingular: ${is_supersingular(E)}`);
        supersingularFound = true;
      }
    } catch {
      // Singular
    }
  }
}

if (!supersingularFound) {
  console.log("  No supersingular curve found over F_97");
}
console.log("");

// Example of an ordinary curve
const E_ord = EllipticCurve(F97, [2n, 3n]);
console.log("Ordinary curve example:");
console.log(`  E: y^2 = x^3 + 2x + 3`);
console.log(`  Order: ${E_ord.cardinality()}`);
console.log(`  Trace: ${97n + 1n - E_ord.cardinality()}`);
console.log(`  Ordinary: ${is_ordinary(E_ord)}`);
console.log("");

// ============================================================================
// Example 4: Special j-Invariants
// ============================================================================

console.log("=== Example 4: Special j-Invariants ===\n");

const F101 = new FiniteFieldPrime(101n);

/**
 * Curves with j = 0:
 * - Extra automorphisms (order 6 instead of 2)
 * - Form: y^2 = x^3 + b (a = 0)
 */
console.log("Curves with j = 0 (form: y^2 = x^3 + b):");
const E_j0 = EllipticCurve(F101, [0n, 1n]);
console.log(`  E: y^2 = x^3 + 1 over F_101`);
console.log(`  j-invariant: ${E_j0.j_invariant()}`);
console.log(`  Order: ${E_j0.cardinality()}`);
console.log("");

/**
 * Curves with j = 1728:
 * - Extra automorphisms (order 4 instead of 2)
 * - Form: y^2 = x^3 + ax (b = 0)
 */
console.log("Curves with j = 1728 (form: y^2 = x^3 + ax):");
const E_j1728 = EllipticCurve(F101, [1n, 0n]);
console.log(`  E: y^2 = x^3 + x over F_101`);
console.log(`  j-invariant: ${E_j1728.j_invariant()}`);
console.log(`  Order: ${E_j1728.cardinality()}`);
console.log("");

// ============================================================================
// Example 5: Anomalous Curves (INSECURE!)
// ============================================================================

console.log("=== Example 5: Anomalous Curves (INSECURE!) ===\n");

/**
 * An ANOMALOUS curve over F_p has #E = p.
 *
 * WARNING: These curves are INSECURE for ECDLP!
 * The Smart-ASS attack solves ECDLP in polynomial time.
 *
 * Never use anomalous curves for cryptography!
 */

console.log("An anomalous curve has #E = p (INSECURE for crypto!)");
console.log("Searching for anomalous curves over small primes...\n");

for (const p of [11n, 13n, 17n, 19n, 23n, 29n, 31n]) {
  const Fp = new FiniteFieldPrime(p);

  for (let a = 0n; a < p; a++) {
    for (let b = 0n; b < p; b++) {
      try {
        const E = EllipticCurve(Fp, [a, b]);
        const order = E.cardinality();
        if (order === p) {
          console.log(`  F_${p}: y^2 = x^3 + ${a}x + ${b} has order ${order} = p`);
          console.log(`    THIS IS ANOMALOUS - DO NOT USE FOR CRYPTO!`);
          break;
        }
      } catch {
        // Singular
      }
    }
  }
}
console.log("");

// ============================================================================
// Example 6: Trace Distribution
// ============================================================================

console.log("=== Example 6: Trace Distribution ===\n");

/**
 * The trace t = p + 1 - #E determines the curve order.
 * By Hasse: |t| <= 2*sqrt(p)
 *
 * The distribution of traces is related to the Sato-Tate conjecture
 * (now a theorem, proved by Taylor and collaborators).
 */

const F53 = new FiniteFieldPrime(53n);
const traceCount = new Map<bigint, number>();

console.log("Trace distribution for curves over F_53:");

for (let a = 1n; a < 53n; a++) {
  for (let b = 1n; b < 53n; b++) {
    try {
      const E = EllipticCurve(F53, [a, b]);
      const trace = E.trace_of_frobenius();
      traceCount.set(trace, (traceCount.get(trace) || 0) + 1);
    } catch {
      // Singular
    }
  }
}

// Sort and display
const sortedTraces = Array.from(traceCount.entries()).sort((a, b) =>
  a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0
);

for (const [trace, count] of sortedTraces) {
  const bar = '#'.repeat(Math.min(count / 10, 50));
  console.log(`  t=${String(trace).padStart(3)}: ${String(count).padStart(4)} ${bar}`);
}
console.log("");

// ============================================================================
// Example 7: Isomorphic Curves
// ============================================================================

console.log("=== Example 7: Isomorphic Curves ===\n");

/**
 * Two curves over F_q with the same j-invariant are isomorphic
 * over the algebraic closure. They may or may not be isomorphic
 * over F_q itself (quadratic twist).
 */

const F103 = new FiniteFieldPrime(103n);

const E1 = EllipticCurve(F103, [1n, 1n]);
const E2 = EllipticCurve(F103, [4n, 8n]); // (a*u^4, b*u^6) for some u

console.log("Comparing two curves:");
console.log(`  E1: y^2 = x^3 + x + 1,  j = ${E1.j_invariant()}, #E = ${E1.cardinality()}`);
console.log(`  E2: y^2 = x^3 + 4x + 8, j = ${E2.j_invariant()}, #E = ${E2.cardinality()}`);

const sameJ = E1.j_invariant().eq(E2.j_invariant());
console.log(`  Same j-invariant? ${sameJ}`);
if (sameJ) {
  console.log("  These curves are isomorphic over the algebraic closure.");
  if (E1.cardinality() === E2.cardinality()) {
    console.log("  Same order -> isomorphic over F_103");
  } else {
    console.log("  Different order -> quadratic twists of each other");
  }
}
console.log("");

// ============================================================================
// Summary Table
// ============================================================================

console.log("=== Summary of Curve Properties ===\n");

console.log("Property                | Cryptographic Impact");
console.log("------------------------|-------------------------------------------");
console.log("Large prime order       | REQUIRED - prevents small subgroup attacks");
console.log("Ordinary (not supersg.) | RECOMMENDED - larger embedding degree");
console.log("j != 0, 1728            | RECOMMENDED - smaller automorphism group");
console.log("#E != p (not anomalous) | REQUIRED - prevents Smart-ASS attack");
console.log("Small cofactor          | RECOMMENDED - efficient, secure");
console.log("");

// ============================================================================
// Exercises
// ============================================================================

console.log("=== Exercises ===\n");

console.log("1. Find all curves over F_11 with prime order");
console.log("");

console.log("2. Compute the j-invariants that appear over F_13");
console.log("   How many distinct j-invariants are there?");
console.log("");

console.log("3. Find a supersingular curve over F_p for p = 2, 3, 5, 7, 11");
console.log("   What pattern do you notice about their orders?");
console.log("");

console.log("4. For curves E: y^2 = x^3 + 1 and E': y^2 = x^3 + 2 over F_43,");
console.log("   verify they have the same j-invariant but different orders");
console.log("   (they are quadratic twists)");
console.log("");
