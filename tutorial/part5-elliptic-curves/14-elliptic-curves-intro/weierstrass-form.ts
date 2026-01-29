/**
 * # Elliptic Curves: The Weierstrass Form
 *
 * Elliptic curves are algebraic curves defined by cubic equations that have been
 * fundamental to number theory since the 17th century. In modern cryptography,
 * they provide the most efficient public-key cryptosystems known.
 *
 * ## What is an Elliptic Curve?
 *
 * An elliptic curve is NOT an ellipse! The name comes from elliptic integrals,
 * which arise when computing the arc length of an ellipse.
 *
 * An elliptic curve over a field K is the set of solutions (x, y) to a cubic
 * equation, together with a special "point at infinity" denoted O.
 *
 * ## The Weierstrass Equation
 *
 * The most general form is the **long Weierstrass form**:
 *
 *   y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
 *
 * The coefficients [a1, a2, a3, a4, a6] are called the **a-invariants**.
 * The strange numbering comes from a weighted degree system.
 *
 * ## Short Weierstrass Form
 *
 * For fields with characteristic not equal to 2 or 3, we can simplify to:
 *
 *   y^2 = x^3 + a*x + b
 *
 * This is the **short Weierstrass form**, used in most cryptographic applications.
 * Here a = a4 and b = a6, with a1 = a2 = a3 = 0.
 *
 * ## Non-Singularity Condition
 *
 * For the curve to be an elliptic curve (and not singular), the discriminant
 * must be non-zero:
 *
 *   Delta = -16(4a^3 + 27b^2) != 0
 *
 * A singular curve has a "cusp" or "node" and does NOT form a group.
 *
 * ## The j-Invariant
 *
 * The j-invariant classifies elliptic curves up to isomorphism:
 *
 *   j = 1728 * 4a^3 / (4a^3 + 27b^2)
 *
 * Two curves over an algebraically closed field are isomorphic iff they have
 * the same j-invariant.
 *
 * @module tutorial/elliptic-curves/weierstrass-form
 */

import { EllipticCurve, GF, FiniteFieldPrime } from 'sagemath-ts';

// ============================================================================
// Creating Elliptic Curves
// ============================================================================

console.log("=== Creating Elliptic Curves ===\n");

// Create a prime field F_23
const F23 = new FiniteFieldPrime(23n);

// Short Weierstrass form: y^2 = x^3 + x + 1 over F_23
// Using coefficients [a, b] = [1, 1]
const E1 = EllipticCurve(F23, [1n, 1n]);

console.log("Curve E1:", E1.toString());
console.log("  Coefficients [a, b] = [1, 1]");
console.log("  Equation: y^2 = x^3 + x + 1");
console.log("");

// Another curve: y^2 = x^3 + 2x + 3
const E2 = EllipticCurve(F23, [2n, 3n]);

console.log("Curve E2:", E2.toString());
console.log("  Coefficients [a, b] = [2, 3]");
console.log("  Equation: y^2 = x^3 + 2x + 3");
console.log("");

// ============================================================================
// Curve Invariants
// ============================================================================

console.log("=== Curve Invariants ===\n");

// The discriminant tells us if the curve is non-singular
const disc1 = E1.discriminant();
console.log(`E1 discriminant: ${disc1}`);
console.log(`  Non-singular: ${!disc1.isZero()}`);

const disc2 = E2.discriminant();
console.log(`E2 discriminant: ${disc2}`);
console.log(`  Non-singular: ${!disc2.isZero()}`);
console.log("");

// The j-invariant classifies curves up to isomorphism
const j1 = E1.j_invariant();
const j2 = E2.j_invariant();

console.log(`E1 j-invariant: ${j1}`);
console.log(`E2 j-invariant: ${j2}`);
console.log(`Same isomorphism class? ${j1.eq(j2)}`);
console.log("");

// ============================================================================
// Special j-Invariants
// ============================================================================

console.log("=== Special j-Invariants ===\n");

/**
 * Special values of j have extra symmetry:
 * - j = 0: Curves with automorphism group of order 6
 * - j = 1728: Curves with automorphism group of order 4
 * - Other j: Automorphism group of order 2 (just +P and -P)
 */

// Find a curve with j = 0
// For j = 0, we need 4a^3 = 0, so a = 0
// Equation: y^2 = x^3 + b
const F101 = new FiniteFieldPrime(101n);
const E_j0 = EllipticCurve(F101, [0n, 1n]); // y^2 = x^3 + 1

console.log("Curve with j = 0:");
console.log(`  E: ${E_j0.toString()}`);
console.log(`  j-invariant: ${E_j0.j_invariant()}`);
console.log("");

// Find a curve with j = 1728
// For j = 1728, we need 27b^2 = 0, so b = 0
// Equation: y^2 = x^3 + ax
const E_j1728 = EllipticCurve(F101, [1n, 0n]); // y^2 = x^3 + x

console.log("Curve with j = 1728:");
console.log(`  E: ${E_j1728.toString()}`);
console.log(`  j-invariant: ${E_j1728.j_invariant()}`);
console.log("");

// ============================================================================
// Singular Curves (What NOT to Use)
// ============================================================================

console.log("=== Singular Curves (Invalid for Cryptography) ===\n");

/**
 * A singular curve has discriminant = 0 and does NOT form a group.
 * The library will throw an error if you try to create one.
 *
 * For y^2 = x^3 + ax + b, singularity occurs when 4a^3 + 27b^2 = 0.
 *
 * Example: y^2 = x^3 (both a=0 and b=0) has a cusp at (0,0).
 * Example: y^2 = x^3 + x^2 has a node at (0,0).
 */

console.log("Attempting to create singular curve y^2 = x^3...");
try {
  const singular = EllipticCurve(F23, [0n, 0n]);
  console.log("ERROR: Should have thrown an exception!");
} catch (e) {
  console.log(`  Correctly rejected: ${(e as Error).message}`);
}
console.log("");

// ============================================================================
// Why Weierstrass Form?
// ============================================================================

console.log("=== Why Use Weierstrass Form? ===\n");

/**
 * The Weierstrass form is universal:
 *
 * 1. Every elliptic curve over a field can be written in Weierstrass form
 *    (possibly the long form for char 2 or 3)
 *
 * 2. Standard form makes comparison easy (via j-invariant)
 *
 * 3. Formulas for point addition are well-known and optimized
 *
 * 4. Cryptographic standards (NIST, SEC, etc.) use this form
 *
 * Other forms exist and can be more efficient:
 * - Montgomery form: By^2 = x^3 + Ax^2 + x (used in Curve25519)
 * - Edwards form: x^2 + y^2 = 1 + dx^2y^2 (used in Ed25519)
 * - But these can all be converted to Weierstrass form
 */

console.log("Key takeaways:");
console.log("  1. Elliptic curves are cubic equations with non-zero discriminant");
console.log("  2. Short Weierstrass form y^2 = x^3 + ax + b is standard");
console.log("  3. The j-invariant classifies curves up to isomorphism");
console.log("  4. Singular curves (Delta = 0) cannot be used for cryptography");
console.log("");

// ============================================================================
// Exercises
// ============================================================================

console.log("=== Exercises ===\n");

console.log("1. Create the curve y^2 = x^3 - 3x + 1 over F_97 and compute its j-invariant");
console.log("");

console.log("2. Find a curve over F_31 with j-invariant equal to 5");
console.log("   Hint: Solve j = 1728 * 4a^3 / (4a^3 + 27b^2) for a and b");
console.log("");

console.log("3. Verify that two curves with the same j-invariant are isomorphic:");
console.log("   E1: y^2 = x^3 + x + 1 over F_101");
console.log("   E2: y^2 = x^3 + 4x + 4 over F_101 (this is a twist)");
console.log("");
