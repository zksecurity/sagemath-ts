/**
 * # Standard Elliptic Curves for Cryptography
 *
 * Choosing the right elliptic curve is critical for security and performance.
 * This tutorial covers standard curves used in practice and their properties.
 *
 * ## Why Standard Curves?
 *
 * Using well-vetted, standardized curves provides:
 * - Security assurance (extensive cryptanalysis)
 * - Interoperability (everyone uses the same parameters)
 * - Optimized implementations (hardware and software)
 * - Regulatory compliance (NIST, FIPS requirements)
 *
 * ## Curve Families
 *
 * 1. **NIST curves** (P-256, P-384, P-521)
 * 2. **Koblitz curves** (secp256k1 used by Bitcoin)
 * 3. **Modern curves** (Curve25519, Ed25519)
 * 4. **Pairing-friendly curves** (BN254, BLS12-381)
 *
 * @module tutorial/elliptic-curves/curve-selection
 */

import { EllipticCurve, FiniteFieldPrime } from 'sagemath-ts';

// ============================================================================
// secp256k1: Bitcoin's Curve
// ============================================================================

console.log("=== secp256k1: The Bitcoin Curve ===\n");

/**
 * secp256k1 Parameters:
 *
 * Prime p = 2^256 - 2^32 - 977
 *         = 0xFFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFE FFFFFC2F
 *
 * Curve: y^2 = x^3 + 7 (a = 0, b = 7)
 *
 * Order n = 0xFFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFE BAAEDCE6 AF48A03B BFD25E8C D0364141
 *
 * Generator G:
 *   Gx = 0x79BE667E F9DCBBAC 55A06295 CE870B07 029BFCDB 2DCE28D9 59F2815B 16F81798
 *   Gy = 0x483ADA77 26A3C465 5DA4FBFC 0E1108A8 FD17B448 A6855419 9C47D08F FB10D4B8
 *
 * Cofactor h = 1 (prime order - ideal!)
 */

console.log("secp256k1 (used by Bitcoin, Ethereum):");
console.log("");
console.log("Parameters:");
console.log("  p = 2^256 - 2^32 - 977");
console.log("  a = 0");
console.log("  b = 7");
console.log("  Equation: y^2 = x^3 + 7");
console.log("");
console.log("Properties:");
console.log("  - Order is prime (cofactor h = 1)");
console.log("  - Koblitz curve (a = 0, efficient endomorphism)");
console.log("  - NOT a NIST standard (independently chosen)");
console.log("  - Efficient implementation due to special form");
console.log("");

// Demonstrate with smaller prime that's actually prime
const p_small = 4294966337n; // Next prime after 2^32 - 977
const F_secp = new FiniteFieldPrime(p_small);
const E_secp = EllipticCurve(F_secp, [0n, 7n]);

console.log("Small demonstration (secp256k1-like with smaller prime):");
console.log(`  p = ${p_small} (a prime near 2^32)`);
console.log(`  Curve order: ${E_secp.cardinality()}`);

const G_secp = E_secp.gens()[0]!;
console.log(`  Generator: ${G_secp.toString()}`);
console.log(`  Generator order: ${G_secp.order()}`);
console.log("");

// ============================================================================
// NIST P-256 (secp256r1)
// ============================================================================

console.log("=== NIST P-256 (secp256r1): The TLS Standard ===\n");

/**
 * P-256 Parameters:
 *
 * Prime p = 2^256 - 2^224 + 2^192 + 2^96 - 1
 *
 * Curve: y^2 = x^3 - 3x + b (a = -3 for efficiency)
 *
 * The parameter b is chosen using SHA-1 of a seed (verifiably random).
 *
 * Properties:
 * - Random curve (not Koblitz)
 * - a = -3 optimization for point doubling
 * - Prime order (cofactor h = 1)
 * - Widely deployed in TLS, code signing, etc.
 */

console.log("NIST P-256 (secp256r1):");
console.log("");
console.log("Parameters:");
console.log("  p = 2^256 - 2^224 + 2^192 + 2^96 - 1");
console.log("  a = -3 (optimizes doubling formula)");
console.log("  b = [verifiably random from SHA-1 of seed]");
console.log("");
console.log("Properties:");
console.log("  - NIST standard (FIPS 186-4)");
console.log("  - Most widely used curve in TLS");
console.log("  - Prime order (cofactor h = 1)");
console.log("  - ~128-bit security level");
console.log("");
console.log("Usage:");
console.log("  - TLS 1.2/1.3 (most HTTPS connections)");
console.log("  - ECDSA certificates");
console.log("  - Apple's Secure Enclave");
console.log("  - AWS KMS, Google Cloud KMS");
console.log("");

// ============================================================================
// Curve25519 and Ed25519
// ============================================================================

console.log("=== Curve25519 & Ed25519: Modern Curves ===\n");

/**
 * Curve25519 (Montgomery form):
 *   y^2 = x^3 + 486662x^2 + x over F_p where p = 2^255 - 19
 *
 * Ed25519 (twisted Edwards form):
 *   -x^2 + y^2 = 1 - (121665/121666)x^2y^2 over the same field
 *
 * These are birationally equivalent curves designed by Daniel J. Bernstein.
 *
 * Advantages:
 * - Complete addition formulas (no special cases)
 * - Fast, constant-time implementations
 * - Designed to resist implementation errors
 * - Smaller code size
 *
 * Curve25519: Used for ECDH (X25519)
 * Ed25519: Used for signatures (EdDSA)
 */

console.log("Curve25519 / Ed25519:");
console.log("");
console.log("Parameters:");
console.log("  p = 2^255 - 19 (efficient Mersenne-like prime)");
console.log("  Security: ~128-bit");
console.log("");
console.log("Design philosophy:");
console.log("  - Safe curves: resistant to common implementation mistakes");
console.log("  - Constant-time: prevent timing side-channel attacks");
console.log("  - Compact: small keys and signatures");
console.log("");
console.log("Forms:");
console.log("  - Curve25519: Montgomery form for fast ECDH (X25519)");
console.log("  - Ed25519: Edwards form for fast signatures (EdDSA)");
console.log("  - Both represent the same underlying curve!");
console.log("");
console.log("Usage:");
console.log("  - Signal Protocol (encrypted messaging)");
console.log("  - SSH (modern key exchange)");
console.log("  - TLS 1.3 (X25519 is mandatory to implement)");
console.log("  - WireGuard VPN");
console.log("  - Tor");
console.log("");

// ============================================================================
// BLS12-381: Pairing-Friendly Curve
// ============================================================================

console.log("=== BLS12-381: Pairing-Friendly Curve ===\n");

/**
 * BLS12-381 is designed for efficient pairings:
 *   e: G1 x G2 -> GT
 *
 * Parameters:
 * - Embedding degree k = 12
 * - ~128-bit security for both ECDLP and pairing attacks
 * - Supports BLS signatures and SNARKs
 *
 * G1: Points on E(F_p)
 * G2: Points on E'(F_{p^2}) (twisted curve)
 * GT: Subgroup of F_{p^12}*
 */

console.log("BLS12-381:");
console.log("");
console.log("Parameters:");
console.log("  - Embedding degree k = 12");
console.log("  - Prime p: 381 bits");
console.log("  - Subgroup order r: 255 bits");
console.log("  - Security: ~128-bit");
console.log("");
console.log("Groups:");
console.log("  - G1: Subgroup of E(F_p), order r");
console.log("  - G2: Subgroup of E'(F_{p^2}), order r");
console.log("  - GT: Subgroup of F_{p^12}*, order r");
console.log("");
console.log("Applications:");
console.log("  - BLS signatures (short, aggregatable)");
console.log("  - zk-SNARKs (Groth16, PLONK)");
console.log("  - Ethereum 2.0 (validator signatures)");
console.log("  - Zcash (shielded transactions)");
console.log("");

// ============================================================================
// Curve Comparison Table
// ============================================================================

console.log("=== Curve Comparison Table ===\n");

console.log("Curve       | Bits | Cofactor | a=-3? | Special Form | Use Case");
console.log("------------|------|----------|-------|--------------|------------------");
console.log("secp256k1   | 256  | 1        | No    | Koblitz      | Bitcoin, Ethereum");
console.log("P-256       | 256  | 1        | Yes   | Random       | TLS, certificates");
console.log("P-384       | 384  | 1        | Yes   | Random       | High-security TLS");
console.log("Curve25519  | 255  | 8        | N/A   | Montgomery   | Modern ECDH");
console.log("Ed25519     | 255  | 8        | N/A   | Edwards      | Modern signatures");
console.log("BLS12-381   | 381  | varies   | N/A   | Pairing      | zk-SNARKs, BLS");
console.log("");

// ============================================================================
// Curve Selection Criteria
// ============================================================================

console.log("=== How to Choose a Curve ===\n");

console.log("1. SECURITY LEVEL");
console.log("   - 128-bit: P-256, secp256k1, Curve25519");
console.log("   - 192-bit: P-384");
console.log("   - 256-bit: P-521");
console.log("");

console.log("2. INTEROPERABILITY");
console.log("   - Web/TLS: P-256 (widest support)");
console.log("   - Blockchain: secp256k1 (Bitcoin/Ethereum ecosystem)");
console.log("   - Modern protocols: Curve25519/Ed25519");
console.log("");

console.log("3. PERFORMANCE");
console.log("   - Fastest: Curve25519/Ed25519 (designed for speed)");
console.log("   - Fast: secp256k1 (Koblitz endomorphism)");
console.log("   - Standard: P-256 (a=-3 optimization)");
console.log("");

console.log("4. SPECIAL REQUIREMENTS");
console.log("   - Pairings needed: BLS12-381");
console.log("   - Constant-time required: Curve25519/Ed25519");
console.log("   - FIPS compliance: NIST curves (P-256/384/521)");
console.log("");

console.log("5. IMPLEMENTATION SAFETY");
console.log("   - Complete formulas: Edwards curves (no edge cases)");
console.log("   - Montgomery ladder: Curve25519 (inherently constant-time)");
console.log("   - Weierstrass: Need careful implementation");
console.log("");

// ============================================================================
// Curve Security Considerations
// ============================================================================

console.log("=== Security Considerations ===\n");

console.log("AVOID these curve properties:");
console.log("");
console.log("1. ANOMALOUS (#E = p)");
console.log("   - Smart's attack solves ECDLP in polynomial time");
console.log("   - All standard curves avoid this");
console.log("");

console.log("2. SUPERSINGULAR");
console.log("   - Low embedding degree (k <= 6)");
console.log("   - MOV attack reduces ECDLP to DLP in extension field");
console.log("   - Only use for isogeny-based crypto");
console.log("");

console.log("3. LOW EMBEDDING DEGREE");
console.log("   - MOV/Frey-Ruck attack if k is small");
console.log("   - Standard curves have astronomically large k");
console.log("");

console.log("4. TWIST ATTACKS");
console.log("   - Invalid curve attacks on implementations");
console.log("   - Solution: Validate received points are on curve");
console.log("");

console.log("5. SMALL SUBGROUP ATTACKS");
console.log("   - If cofactor h > 1, adversary can exploit subgroups");
console.log("   - Solution: Multiply points by h (cofactor clearing)");
console.log("");

// ============================================================================
// Implementing with Standard Parameters
// ============================================================================

console.log("=== Using Standard Parameters ===\n");

/**
 * When implementing EC crypto, use standardized parameters:
 *
 * 1. For general TLS/HTTPS: P-256
 * 2. For Bitcoin/Ethereum: secp256k1
 * 3. For modern protocols: X25519 (ECDH), Ed25519 (signatures)
 * 4. For zk-proofs: BLS12-381
 *
 * NEVER generate your own curve parameters!
 * - Curve generation is subtle and easy to get wrong
 * - Standard curves have been extensively analyzed
 * - Custom curves raise red flags in security audits
 */

console.log("Golden rules:");
console.log("  1. USE STANDARD CURVES - never generate your own");
console.log("  2. USE STANDARD IMPLEMENTATIONS - don't roll your own crypto");
console.log("  3. VALIDATE INPUTS - check points are on the curve");
console.log("  4. HANDLE COFACTORS - multiply by h for subgroup safety");
console.log("");

// Example: Demonstrate parameter validation
console.log("Example: Validating a point is on secp256k1-like curve:");

const x_test = 12345n;
const E_test = E_secp;

if (E_test.is_x_coord(x_test)) {
  const P_test = E_test.lift_x(x_test);
  console.log(`  x = ${x_test}`);
  console.log(`  Valid point: ${P_test.toString()}`);

  // Verify it's on the curve
  const px = P_test.x!.value;
  const py = P_test.y!.value;
  const lhs = (py * py) % p_small;
  const rhs = (px * px * px + 7n) % p_small;
  console.log(`  y^2 mod p = ${lhs}`);
  console.log(`  x^3 + 7 mod p = ${rhs}`);
  console.log(`  On curve: ${lhs === rhs}`);
} else {
  console.log(`  x = ${x_test} does not correspond to a curve point`);
}
console.log("");

// ============================================================================
// Exercises
// ============================================================================

console.log("=== Exercises ===\n");

console.log("1. Compare P-256 and secp256k1:");
console.log("   - Why does P-256 use a = -3?");
console.log("   - Why does secp256k1 use a = 0?");
console.log("   - Which is faster and why?");
console.log("");

console.log("2. Research Curve25519:");
console.log("   - What is the Montgomery ladder?");
console.log("   - Why is cofactor 8 acceptable here?");
console.log("   - What makes it 'safe' against implementation errors?");
console.log("");

console.log("3. Investigate a curve not covered here:");
console.log("   - Brainpool curves (used in Europe)");
console.log("   - SM2 curve (Chinese standard)");
console.log("   - NUMS curves (Microsoft)");
console.log("");
