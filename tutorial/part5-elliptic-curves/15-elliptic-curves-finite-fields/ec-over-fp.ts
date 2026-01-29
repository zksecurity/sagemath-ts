/**
 * # Elliptic Curves Over Prime Fields F_p
 *
 * Elliptic curves over prime fields F_p are the workhorses of modern cryptography.
 * Every standard curve (secp256k1, P-256, P-384, etc.) is defined over a prime field.
 *
 * ## Why Prime Fields?
 *
 * 1. **Simple arithmetic**: Operations are just modular arithmetic
 * 2. **Well-understood**: Decades of cryptanalysis and optimization
 * 3. **Hardware support**: Many CPUs have optimized big integer ops
 * 4. **Standards compliance**: NIST, SEC, Brainpool all use prime fields
 *
 * ## The Curve Equation
 *
 * Over F_p (p > 3), elliptic curves are given in short Weierstrass form:
 *
 *   y^2 = x^3 + ax + b  (mod p)
 *
 * where a, b are elements of F_p and 4a^3 + 27b^2 != 0 (non-singular).
 *
 * @module tutorial/elliptic-curves/ec-over-fp
 */

import { EllipticCurve, FiniteFieldPrime } from 'sagemath-ts';

// ============================================================================
// Creating Curves Over Prime Fields
// ============================================================================

console.log("=== Creating Curves Over Prime Fields ===\n");

// A small prime field for demonstration
const F101 = new FiniteFieldPrime(101n);

// Create a curve: y^2 = x^3 + 2x + 3 over F_101
const E = EllipticCurve(F101, [2n, 3n]);

console.log(`Field: F_${F101.characteristic}`);
console.log(`Curve: ${E.toString()}`);
console.log("");

// ============================================================================
// Field Arithmetic Review
// ============================================================================

console.log("=== Field Arithmetic in F_p ===\n");

/**
 * In F_p, all arithmetic is done modulo p:
 * - Addition: (a + b) mod p
 * - Subtraction: (a - b) mod p
 * - Multiplication: (a * b) mod p
 * - Division: a * b^(-1) mod p (using Fermat's little theorem: b^(-1) = b^(p-2))
 */

const a = F101.__call__(47n);
const b = F101.__call__(73n);

console.log(`In F_101:`);
console.log(`  a = ${a},  b = ${b}`);
console.log(`  a + b = ${a.add(b)} (mod 101)`);
console.log(`  a - b = ${a.sub(b)} (mod 101)`);
console.log(`  a * b = ${a.mul(b)} (mod 101)`);
console.log(`  a / b = ${a.div(b)} (mod 101)`);
console.log(`  a^(-1) = ${a.inv()} (verify: ${a.mul(a.inv())})`);
console.log("");

// ============================================================================
// Points on the Curve
// ============================================================================

console.log("=== Points on the Curve ===\n");

// The point at infinity (identity)
const O = E.zero();
console.log(`Point at infinity: ${O.toString()}`);

// Find points by lifting x-coordinates
console.log("\nFinding points on E: y^2 = x^3 + 2x + 3:");

for (let x = 0n; x <= 10n; x++) {
  // Check if there's a point with this x-coordinate
  if (E.is_x_coord(x)) {
    const points = E.lift_x(x, true);
    for (const P of points) {
      console.log(`  x=${x}: ${P.toString()}`);
    }
  } else {
    console.log(`  x=${x}: no point (x^3 + 2x + 3 is not a QR)`);
  }
}
console.log("");

// ============================================================================
// Verifying Points Are On the Curve
// ============================================================================

console.log("=== Verifying Points ===\n");

const P = E.lift_x(3n);
console.log(`Point P = ${P.toString()}`);

// Manually verify: y^2 should equal x^3 + 2x + 3
if (!P.isZero()) {
  const x = P.x!;
  const y = P.y!;

  const lhs = y.mul(y);  // y^2
  const rhs = x.mul(x).mul(x).add(F101.__call__(2n).mul(x)).add(F101.__call__(3n));  // x^3 + 2x + 3

  console.log(`  LHS: y^2 = ${y}^2 = ${lhs}`);
  console.log(`  RHS: x^3 + 2x + 3 = ${x}^3 + 2*${x} + 3 = ${rhs}`);
  console.log(`  LHS == RHS? ${lhs.eq(rhs)}`);
}
console.log("");

// ============================================================================
// Random Points
// ============================================================================

console.log("=== Generating Random Points ===\n");

/**
 * Random point generation is crucial for:
 * - Key generation (private key -> random scalar -> public key)
 * - Randomized protocols (e.g., ElGamal encryption)
 * - Testing and verification
 */

console.log("Five random points on E:");
for (let i = 0; i < 5; i++) {
  const randP = E.random_point();
  console.log(`  Random point ${i + 1}: ${randP.toString()}`);
}
console.log("");

// ============================================================================
// Point Operations
// ============================================================================

console.log("=== Point Operations Over F_p ===\n");

// Find Q by using a different valid point (x=5 works: 5^3 + 2*5 + 3 = 138 ≡ 37 mod 101, sqrt exists)
const Q = E.lift_x(5n);
console.log(`P = ${P.toString()}`);
console.log(`Q = ${Q.toString()}`);
console.log("");

// Addition
const sum = P.add(Q);
console.log(`P + Q = ${sum.toString()}`);

// Doubling
const doubled = P.double();
console.log(`2P = ${doubled.toString()}`);

// Scalar multiplication
const scalar = 17n;
const scalarMul = P.mul(scalar);
console.log(`${scalar}P = ${scalarMul.toString()}`);
console.log("");

// ============================================================================
// Large Prime Fields (Cryptographic Size)
// ============================================================================

console.log("=== Cryptographic-Size Prime Fields ===\n");

/**
 * Real cryptographic curves use primes of 256+ bits.
 * Common choices:
 * - secp256k1: p = 2^256 - 2^32 - 977
 * - P-256: p = 2^256 - 2^224 + 2^192 + 2^96 - 1
 * - P-384: p = 2^384 - 2^128 - 2^96 + 2^32 - 1
 *
 * These primes are chosen for:
 * 1. Efficient reduction (special form)
 * 2. Security (large enough to resist attacks)
 * 3. Compatibility (standard bit sizes)
 */

// Example with a medium-sized prime (for demonstration)
const p256_like = 2n ** 32n - 5n;  // A 32-bit prime for demonstration
const F_large = new FiniteFieldPrime(p256_like);
const E_large = EllipticCurve(F_large, [2n, 3n]);

console.log(`Larger field: F_${p256_like}`);
console.log(`Bit size: ${p256_like.toString(2).length} bits`);
console.log(`Curve: ${E_large.toString()}`);

const P_large = E_large.random_point();
console.log(`Random point: ${P_large.toString()}`);
console.log(`Point order: ${P_large.order()}`);
console.log("");

// ============================================================================
// Quadratic Residues and Point Existence
// ============================================================================

console.log("=== Quadratic Residues ===\n");

/**
 * A point (x, y) exists on y^2 = x^3 + ax + b iff x^3 + ax + b is a
 * quadratic residue (QR) in F_p.
 *
 * In F_p with p odd:
 * - Exactly (p-1)/2 non-zero elements are QRs
 * - Testing: a is QR iff a^((p-1)/2) = 1 (Euler's criterion)
 * - Square root: Use Tonelli-Shanks algorithm
 */

console.log("Testing which x-values give points on E over F_101:");

let numPoints = 1;  // Count point at infinity
let numWithY = 0;
let numWithTwoY = 0;

for (let x = 0n; x < 101n; x++) {
  const xElem = F101.__call__(x);
  // Compute x^3 + 2x + 3
  const rhs = xElem.mul(xElem).mul(xElem)
    .add(F101.__call__(2n).mul(xElem))
    .add(F101.__call__(3n));

  if (rhs.isZero()) {
    numPoints += 1;  // y = 0, one point
    numWithY++;
  } else if (rhs.is_square()) {
    numPoints += 2;  // y = +/- sqrt(rhs), two points
    numWithTwoY++;
  }
  // else: no point with this x
}

console.log(`  Total points on E(F_101): ${numPoints}`);
console.log(`  x-values with y=0: ${numWithY}`);
console.log(`  x-values with two y values: ${numWithTwoY}`);
console.log(`  x-values with no point: ${101 - numWithY - numWithTwoY}`);

// Verify with cardinality
const cardE = E.cardinality();
console.log(`  Cardinality from library: ${cardE}`);
console.log(`  Match? ${BigInt(numPoints) === cardE}`);
console.log("");

// ============================================================================
// Curve Order and Hasse Bound
// ============================================================================

console.log("=== Curve Order and Hasse Bound ===\n");

/**
 * For E/F_p, the Hasse bound says:
 *
 *   |#E(F_p) - (p + 1)| <= 2*sqrt(p)
 *
 * Write #E = p + 1 - t, where t is the "trace of Frobenius".
 * Then |t| <= 2*sqrt(p).
 */

const p = 101n;
const sqrtP = BigInt(Math.floor(Math.sqrt(Number(p))));
const hasseUpper = p + 1n + 2n * sqrtP;
const hasseLower = p + 1n - 2n * sqrtP;

console.log(`For F_${p}:`);
console.log(`  Hasse bound: ${hasseLower} <= #E <= ${hasseUpper}`);
console.log(`  Actual #E = ${cardE}`);
console.log(`  Trace t = ${p + 1n - cardE}`);
console.log(`  Within bound? ${cardE >= hasseLower && cardE <= hasseUpper}`);
console.log("");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Key Points for Elliptic Curves Over F_p ===\n");

console.log("1. Curve equation: y^2 = x^3 + ax + b (mod p)");
console.log("2. Non-singularity: 4a^3 + 27b^2 != 0 (mod p)");
console.log("3. Points: solutions (x,y) plus point at infinity O");
console.log("4. Arithmetic: all computations mod p");
console.log("5. Order: #E(F_p) is between p+1-2sqrt(p) and p+1+2sqrt(p)");
console.log("");

// ============================================================================
// Exercises
// ============================================================================

console.log("=== Exercises ===\n");

console.log("1. Enumerate all points on y^2 = x^3 + x + 1 over F_13");
console.log("   Verify the count matches the cardinality()");
console.log("");

console.log("2. Find all x in F_23 such that x^3 + 5x + 7 is a quadratic residue");
console.log("");

console.log("3. For p = 127, find a curve E with #E = 128 (supersingular)");
console.log("   and another with #E = 127 (anomalous - DO NOT USE)");
console.log("");
