/**
 * # Elliptic Curves Over Extension Fields F_q
 *
 * While most cryptographic applications use prime fields F_p, elliptic curves
 * can also be defined over extension fields F_q where q = p^n for some n > 1.
 *
 * ## Extension Fields
 *
 * An extension field F_q = F_{p^n} is constructed by:
 * 1. Taking a prime field F_p
 * 2. Adjoining a root of an irreducible polynomial of degree n
 *
 * Example: F_4 = F_2[x]/(x^2 + x + 1)
 * Elements: {0, 1, alpha, alpha+1} where alpha^2 + alpha + 1 = 0
 *
 * ## Characteristic 2 Fields
 *
 * Binary extension fields F_{2^n} were popular for:
 * - Hardware efficiency (XOR and shift operations)
 * - NIST B-curves (B-163, B-233, B-283, B-409, B-571)
 *
 * However, recent index calculus advances have made them less attractive.
 *
 * ## Tower Fields
 *
 * For pairing-based cryptography, we often work with tower extensions:
 * - F_p -> F_{p^2} -> F_{p^6} -> F_{p^12}
 * - Used in BLS signatures, BN curves, etc.
 *
 * @module tutorial/elliptic-curves/ec-over-fq
 */

import { EllipticCurve, FiniteFieldPrime } from 'sagemath-ts';

// ============================================================================
// Conceptual Introduction
// ============================================================================

console.log("=== Elliptic Curves Over Extension Fields ===\n");

/**
 * While sagemath-ts currently focuses on prime fields, understanding
 * extension fields is important for:
 * - Pairing-based cryptography
 * - Characteristic 2 curves (historical)
 * - Torsion point analysis
 * - Complex multiplication theory
 */

console.log("Extension fields F_q = F_{p^n} extend prime fields F_p");
console.log("");

console.log("Key examples:");
console.log("  - F_4 = F_2[x]/(x^2 + x + 1): Binary extension");
console.log("  - F_9 = F_3[x]/(x^2 + 1): Ternary extension");
console.log("  - F_{p^12}: Tower field for BN254 pairings");
console.log("");

// ============================================================================
// Simulating Extension Field Concepts with Prime Fields
// ============================================================================

console.log("=== Understanding Through Prime Field Analogies ===\n");

// We'll use a prime field but explain concepts that generalize
const F101 = new FiniteFieldPrime(101n);
const E = EllipticCurve(F101, [1n, 1n]);

console.log(`Working with E: y^2 = x^3 + x + 1 over F_101`);
console.log(`This is a prime field, but concepts generalize to F_q`);
console.log("");

// ============================================================================
// Frobenius Endomorphism
// ============================================================================

console.log("=== The Frobenius Endomorphism ===\n");

/**
 * For a curve E over F_q, the Frobenius map is:
 *   phi: (x, y) -> (x^q, y^q)
 *
 * Properties:
 * 1. phi is an endomorphism of E
 * 2. phi fixes exactly E(F_q) (the F_q-rational points)
 * 3. phi satisfies: phi^2 - t*phi + q = 0 (in the endomorphism ring)
 *    where t = trace of Frobenius
 *
 * For prime fields F_p:
 *   (x, y) -> (x^p, y^p) = (x, y) for all points in E(F_p)
 *   (since x, y are in F_p and x^p = x by Fermat's little theorem)
 */

const trace = E.trace_of_frobenius();
console.log(`Trace of Frobenius: t = ${trace}`);
console.log(`Frobenius satisfies: phi^2 - ${trace}*phi + 101 = 0`);
console.log("");

// Verify: for points in E(F_p), Frobenius is identity
const P = E.random_point();
if (!P.isZero()) {
  const x = P.x!;
  const y = P.y!;

  // x^p mod p = x for x in F_p
  const xFrobenius = x.pow(101n);
  const yFrobenius = y.pow(101n);

  console.log(`Point P = (${x}, ${y})`);
  console.log(`Frobenius(P) = (${x}^101, ${y}^101) = (${xFrobenius}, ${yFrobenius})`);
  console.log(`Frobenius fixes F_p points: ${x.eq(xFrobenius) && y.eq(yFrobenius)}`);
}
console.log("");

// ============================================================================
// Points Over Field Extensions
// ============================================================================

console.log("=== Points Over Extension Fields ===\n");

/**
 * Given E/F_p, we can consider points with coordinates in extension fields:
 *
 * E(F_p) subset E(F_{p^2}) subset E(F_{p^3}) subset ...
 *
 * Example: E might have 100 points over F_p but 10000 points over F_{p^2}.
 *
 * The TORSION POINTS E[n] (points P with nP = O) often require
 * extension fields to exist fully.
 *
 * For E[n] to be fully defined over F_q:
 * - We need q^k - 1 divisible by n for some k (the embedding degree)
 */

const cardP = E.cardinality();
console.log(`#E(F_101) = ${cardP}`);

// The cardinality over F_{p^2} can be computed from the trace
// #E(F_{p^2}) = p^2 + 1 - (t^2 - 2p) where t is trace over F_p
const cardP2 = 101n ** 2n + 1n - (trace ** 2n - 2n * 101n);
console.log(`#E(F_{101^2}) = ${cardP2} (computed from trace)`);
console.log("");

// ============================================================================
// Embedding Degree
// ============================================================================

console.log("=== Embedding Degree ===\n");

/**
 * The EMBEDDING DEGREE k of E with respect to n is the smallest k such that
 * n divides q^k - 1.
 *
 * This is the extension degree needed for the n-torsion to embed into
 * the multiplicative group F_{q^k}*.
 *
 * For pairing-based crypto:
 * - k should be moderate (6-48) for efficiency
 * - k must provide adequate security in F_{q^k}*
 *
 * BN curves: k = 12 (popular for SNARKs)
 * BLS curves: k = 12 or 24
 */

// Find embedding degree for a subgroup of E
const curveOrder = E.cardinality();
const pointOrder = P.order();

console.log(`Looking for embedding degree with respect to n = ${pointOrder}`);

// Find smallest k such that pointOrder | (101^k - 1)
let k = 1n;
let qPowK = 101n;

while (k <= 20n) {
  if ((qPowK - 1n) % pointOrder === 0n) {
    console.log(`  Embedding degree k = ${k}`);
    console.log(`  101^${k} - 1 = ${qPowK - 1n} = ${pointOrder} * ${(qPowK - 1n) / pointOrder}`);
    break;
  }
  k++;
  qPowK *= 101n;
}

if (k > 20n) {
  console.log(`  Embedding degree > 20`);
}
console.log("");

// ============================================================================
// Characteristic 2 Curves (Historical Context)
// ============================================================================

console.log("=== Characteristic 2 Curves (Historical) ===\n");

/**
 * Binary fields F_{2^n} were popular because:
 * - Addition is XOR (very fast in hardware)
 * - Squaring is efficient (Frobenius map)
 * - Field multiplication can use shift-and-XOR
 *
 * Curve equation in characteristic 2 (non-supersingular):
 *   y^2 + xy = x^3 + ax^2 + b
 *
 * NIST B-curves: B-163, B-233, B-283, B-409, B-571
 *
 * HOWEVER: Recent advances in index calculus for binary fields
 * (Joux et al., 2014) have made them less attractive for new deployments.
 * Modern recommendations favor prime fields.
 */

console.log("Binary curves y^2 + xy = x^3 + ax^2 + b were popular for:");
console.log("  - Hardware implementation (XOR operations)");
console.log("  - Constrained devices (smart cards, IoT)");
console.log("");
console.log("Caution: Index calculus attacks have improved significantly.");
console.log("Modern recommendation: Use prime field curves instead.");
console.log("");

// ============================================================================
// Tower Fields for Pairings
// ============================================================================

console.log("=== Tower Fields for Pairing Crypto ===\n");

/**
 * Pairing-based cryptography uses bilinear maps:
 *   e: G1 x G2 -> GT
 *
 * where G1, G2 are subgroups of E and GT is a subgroup of F_{q^k}*.
 *
 * For efficiency, F_{q^k} is constructed as a tower:
 *   F_p -> F_{p^2} -> F_{p^6} -> F_{p^12}
 *
 * Example for BN254 (k=12):
 *   F_{p^2} = F_p[u]/(u^2 + 1)
 *   F_{p^6} = F_{p^2}[v]/(v^3 - (9+u))
 *   F_{p^12} = F_{p^6}[w]/(w^2 - v)
 *
 * This tower structure enables efficient arithmetic and the
 * "optimal ate pairing" computation.
 */

console.log("Tower field construction for BN254:");
console.log("  Level 1: F_p (prime field)");
console.log("  Level 2: F_{p^2} = F_p[u]/(u^2 + 1)");
console.log("  Level 3: F_{p^6} = F_{p^2}[v]/(v^3 - (9+u))");
console.log("  Level 4: F_{p^12} = F_{p^6}[w]/(w^2 - v)");
console.log("");
console.log("This enables efficient pairing computation for:");
console.log("  - BLS signatures");
console.log("  - Groth16 SNARKs");
console.log("  - Identity-based encryption");
console.log("");

// ============================================================================
// Weil and Tate Pairings
// ============================================================================

console.log("=== Pairings on Elliptic Curves ===\n");

/**
 * The Weil pairing e_n: E[n] x E[n] -> mu_n is a bilinear map where
 * mu_n is the group of n-th roots of unity in F_{q^k}*.
 *
 * Properties:
 * 1. Bilinear: e_n(aP, bQ) = e_n(P, Q)^(ab)
 * 2. Non-degenerate: e_n(P, Q) != 1 for independent P, Q
 * 3. Alternating: e_n(P, P) = 1
 *
 * The Tate pairing and its variants (ate, optimal ate) are used in practice.
 *
 * Applications:
 * - Identity-based encryption
 * - Short signatures (BLS)
 * - Zero-knowledge proofs (Groth16, PLONK)
 */

console.log("The Weil pairing e_n maps pairs of n-torsion points to");
console.log("n-th roots of unity in the multiplicative group.");
console.log("");
console.log("Bilinearity: e_n(aP, bQ) = e_n(P, Q)^(ab)");
console.log("");
console.log("This enables:");
console.log("  - Decisional Diffie-Hellman in G1 is EASY (via pairing)");
console.log("  - But Computational Diffie-Hellman in G1 is still HARD");
console.log("  - Opening new cryptographic primitives (IBE, ABE, etc.)");
console.log("");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");

console.log("Extension fields F_q = F_{p^n} generalize prime fields:");
console.log("");
console.log("1. Construction: F_p[x]/(f(x)) where f is irreducible");
console.log("2. Elements: polynomials of degree < n with coefficients in F_p");
console.log("3. Arithmetic: polynomial arithmetic modulo f(x)");
console.log("");
console.log("Cryptographic relevance:");
console.log("  - Pairing-based crypto requires extension fields");
console.log("  - Torsion points may live in extensions");
console.log("  - Binary fields (char 2) are now discouraged");
console.log("");

// ============================================================================
// Exercises
// ============================================================================

console.log("=== Exercises ===\n");

console.log("1. Compute #E(F_{p^2}) for E: y^2 = x^3 + x + 1 over F_7");
console.log("   using the formula #E(F_{p^2}) = p^2 + 1 - (t^2 - 2p)");
console.log("");

console.log("2. Find the embedding degree of E: y^2 = x^3 + 2x + 3 over F_37");
console.log("   with respect to a point of order 19");
console.log("");

console.log("3. (Advanced) Explain why the embedding degree affects");
console.log("   the security of pairing-based cryptography");
console.log("");
