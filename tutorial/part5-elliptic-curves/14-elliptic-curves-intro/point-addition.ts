/**
 * # Point Addition on Elliptic Curves
 *
 * The magic of elliptic curves lies in the fact that points on the curve
 * can be "added" together to produce another point on the curve.
 *
 * ## Geometric Intuition
 *
 * Consider two points P and Q on an elliptic curve y^2 = x^3 + ax + b.
 *
 * **To add P + Q:**
 * 1. Draw a line through P and Q
 * 2. This line intersects the curve at exactly one more point, R
 * 3. Reflect R across the x-axis to get P + Q
 *
 * Why does this work? Bezout's theorem says a degree-2 curve (the line)
 * and a degree-3 curve (the cubic) intersect at exactly 2*3 = 6 points
 * (counting multiplicity and points at infinity). Since the line has
 * degree 1 in y, it contributes 3 intersection points: P, Q, and R.
 *
 * **To double P (compute P + P = 2P):**
 * 1. Draw the tangent line to the curve at P
 * 2. This line intersects the curve at one more point, R
 * 3. Reflect R across the x-axis to get 2P
 *
 * ## The Point at Infinity
 *
 * When P = -Q (same x-coordinate, opposite y), the line through P and Q
 * is vertical and doesn't intersect the curve again. We say it intersects
 * at the "point at infinity" O.
 *
 * The point at infinity O serves as the identity element: P + O = P for all P.
 *
 * @module tutorial/elliptic-curves/point-addition
 */

import { EllipticCurve, EllipticCurvePoint, FiniteFieldPrime } from 'sagemath-ts';

// ============================================================================
// Basic Point Operations
// ============================================================================

console.log("=== Basic Point Operations ===\n");

// Create a curve over F_23
const F23 = new FiniteFieldPrime(23n);
const E = EllipticCurve(F23, [1n, 1n]); // y^2 = x^3 + x + 1

console.log(`Working with curve: ${E.toString()}`);
console.log("");

// Create some points on the curve
// First, find valid points by checking which x-coordinates work
console.log("Finding points on the curve...");

// Point at infinity (identity element)
const O = E.zero();
console.log(`Point at infinity O: ${O.toString()}`);

// Find a point by lifting an x-coordinate
try {
  const P = E.lift_x(0n);
  console.log(`Point P with x=0: ${P.toString()}`);

  const Q = E.lift_x(3n);
  console.log(`Point Q with x=3: ${Q.toString()}`);
  console.log("");

  // ============================================================================
  // Point Addition
  // ============================================================================

  console.log("=== Point Addition ===\n");

  // Add two different points
  const R = P.add(Q);
  console.log(`P + Q = ${R.toString()}`);

  // Addition is commutative
  const R2 = Q.add(P);
  console.log(`Q + P = ${R2.toString()}`);
  console.log(`P + Q == Q + P? ${R.eq(R2)}`);
  console.log("");

  // ============================================================================
  // Point Doubling
  // ============================================================================

  console.log("=== Point Doubling ===\n");

  // Double a point (P + P = 2P)
  const twoP = P.add(P);
  console.log(`2P = P + P = ${twoP.toString()}`);

  // Using the double() method (more efficient internally)
  const twoPfast = P.double();
  console.log(`2P (using double) = ${twoPfast.toString()}`);
  console.log(`Same result? ${twoP.eq(twoPfast)}`);
  console.log("");

  // ============================================================================
  // The Identity Element
  // ============================================================================

  console.log("=== The Identity Element O ===\n");

  // Adding the identity does nothing
  const PplusO = P.add(O);
  console.log(`P + O = ${PplusO.toString()}`);
  console.log(`P + O == P? ${PplusO.eq(P)}`);

  // Identity works on both sides
  const OplusP = O.add(P);
  console.log(`O + P = ${OplusP.toString()}`);
  console.log(`O + P == P? ${OplusP.eq(P)}`);
  console.log("");

  // ============================================================================
  // Point Negation and Inverse
  // ============================================================================

  console.log("=== Point Negation ===\n");

  // Negate a point: -P = (x, -y)
  const negP = P.neg();
  console.log(`P = ${P.toString()}`);
  console.log(`-P = ${negP.toString()}`);

  // P + (-P) = O
  const shouldBeO = P.add(negP);
  console.log(`P + (-P) = ${shouldBeO.toString()}`);
  console.log(`P + (-P) == O? ${shouldBeO.isZero()}`);
  console.log("");

  // ============================================================================
  // Addition Formulas (The Math)
  // ============================================================================

  console.log("=== Addition Formulas ===\n");

  /**
   * For points P = (x1, y1) and Q = (x2, y2) on y^2 = x^3 + ax + b:
   *
   * If P != Q (general addition):
   *   lambda = (y2 - y1) / (x2 - x1)      (slope of line PQ)
   *   x3 = lambda^2 - x1 - x2
   *   y3 = lambda * (x1 - x3) - y1
   *
   * If P = Q (doubling):
   *   lambda = (3*x1^2 + a) / (2*y1)      (slope of tangent)
   *   x3 = lambda^2 - 2*x1
   *   y3 = lambda * (x1 - x3) - y1
   *
   * Let's verify these formulas manually:
   */

  console.log("Manual verification of addition formula:");
  const x1 = P.x;
  const y1 = P.y;
  const x2 = Q.x;
  const y2 = Q.y;

  if (x1 && y1 && x2 && y2) {
    console.log(`  P = (${x1}, ${y1})`);
    console.log(`  Q = (${x2}, ${y2})`);

    // Compute lambda = (y2 - y1) / (x2 - x1)
    const dy = y2.sub(y1);
    const dx = x2.sub(x1);
    const lambda = dy.div(dx);
    console.log(`  lambda = (${y2} - ${y1}) / (${x2} - ${x1}) = ${lambda}`);

    // Compute x3 = lambda^2 - x1 - x2
    const x3 = lambda.mul(lambda).sub(x1).sub(x2);
    console.log(`  x3 = ${lambda}^2 - ${x1} - ${x2} = ${x3}`);

    // Compute y3 = lambda * (x1 - x3) - y1
    const y3 = lambda.mul(x1.sub(x3)).sub(y1);
    console.log(`  y3 = ${lambda} * (${x1} - ${x3}) - ${y1} = ${y3}`);

    console.log(`  Computed P + Q = (${x3}, ${y3})`);
    console.log(`  Library result = ${R.toString()}`);
  }
  console.log("");

  // ============================================================================
  // Subtraction
  // ============================================================================

  console.log("=== Point Subtraction ===\n");

  // Subtraction: P - Q = P + (-Q)
  const PminusQ = P.sub(Q);
  console.log(`P - Q = ${PminusQ.toString()}`);

  // Verify: (P - Q) + Q should equal P
  const verify = PminusQ.add(Q);
  console.log(`(P - Q) + Q = ${verify.toString()}`);
  console.log(`(P - Q) + Q == P? ${verify.eq(P)}`);
  console.log("");

  // ============================================================================
  // Associativity
  // ============================================================================

  console.log("=== Associativity ===\n");

  // Point addition is associative: (P + Q) + R = P + (Q + R)
  const S = E.lift_x(5n);
  console.log(`S = ${S.toString()}`);

  const leftAssoc = P.add(Q).add(S);
  const rightAssoc = P.add(Q.add(S));

  console.log(`(P + Q) + S = ${leftAssoc.toString()}`);
  console.log(`P + (Q + S) = ${rightAssoc.toString()}`);
  console.log(`Associative? ${leftAssoc.eq(rightAssoc)}`);
  console.log("");

} catch (e) {
  console.log(`Error: ${(e as Error).message}`);
}

// ============================================================================
// Why This Matters for Cryptography
// ============================================================================

console.log("=== Why This Matters ===\n");

console.log("Point addition gives elliptic curves a GROUP structure:");
console.log("  - Closure: P + Q is always on the curve");
console.log("  - Associativity: (P + Q) + R = P + (Q + R)");
console.log("  - Identity: O (point at infinity)");
console.log("  - Inverses: Every P has -P such that P + (-P) = O");
console.log("");
console.log("This group structure is the foundation of elliptic curve cryptography!");
console.log("The difficulty of the 'discrete log problem' in this group");
console.log("(finding n given P and nP) is what makes ECC secure.");
console.log("");

// ============================================================================
// Exercises
// ============================================================================

console.log("=== Exercises ===\n");

console.log("1. On E: y^2 = x^3 + 2x + 3 over F_97:");
console.log("   - Find three points P, Q, R on the curve");
console.log("   - Verify P + (Q + R) = (P + Q) + R");
console.log("");

console.log("2. Compute 2P, 3P, 4P, ... until you get back to O or a repeated point");
console.log("   This gives you the order of P in the group");
console.log("");

console.log("3. Implement point addition manually using the formulas above");
console.log("   and verify your results match the library");
console.log("");
