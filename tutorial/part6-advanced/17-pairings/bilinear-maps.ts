/**
 * # Bilinear Maps: The Mathematical Foundation of Pairings
 *
 * Bilinear maps are one of the most powerful tools in modern cryptography.
 * They enable constructions that were previously thought impossible, including:
 * - Identity-Based Encryption (IBE)
 * - Short signatures (BLS)
 * - Tripartite key exchange
 * - Attribute-Based Encryption
 * - Zero-Knowledge Proofs (SNARKs)
 *
 * ## What is a Bilinear Map?
 *
 * A bilinear map (or pairing) is a function e: G1 x G2 -> GT where:
 * - G1 and G2 are additive groups (typically points on elliptic curves)
 * - GT is a multiplicative group (typically in a finite field extension)
 *
 * The key property is BILINEARITY:
 *   e(aP, bQ) = e(P, Q)^(ab)  for all scalars a, b
 *
 * This seemingly simple property is what makes pairings so powerful.
 *
 * ## Why Bilinearity Matters
 *
 * Without pairings, we can compute:
 *   - aP (scalar multiplication in G1)
 *   - bQ (scalar multiplication in G2)
 *
 * But there's NO efficient way to compute:
 *   - Given P, aP, Q, bQ -> compute abP (this would break Diffie-Hellman!)
 *
 * WITH pairings, we CAN verify relationships:
 *   - e(aP, bQ) = e(P, Q)^(ab) = e(abP, Q) = e(P, abQ)
 *
 * This enables checking multiplicative relationships between discrete logs!
 *
 * @module tutorial/part6-advanced/17-pairings/bilinear-maps
 */

import { GF, EllipticCurve, weil_pairing } from 'sagemath-ts';

console.log("=== Bilinear Maps: Foundation of Pairing-Based Cryptography ===\n");

// ============================================================================
// Section 1: Understanding Bilinearity
// ============================================================================

console.log("--- Section 1: Understanding Bilinearity ---\n");

/**
 * MATHEMATICAL DEFINITION
 *
 * A bilinear map e: G1 x G2 -> GT satisfies:
 *
 * 1. BILINEARITY (the key property):
 *    - e(P1 + P2, Q) = e(P1, Q) * e(P2, Q)  (linear in first argument)
 *    - e(P, Q1 + Q2) = e(P, Q1) * e(P, Q2)  (linear in second argument)
 *
 * 2. NON-DEGENERACY:
 *    - If P generates G1 and Q generates G2, then e(P, Q) generates GT
 *    - Equivalently: e(P, Q) != 1 for generators P, Q
 *
 * 3. COMPUTABILITY:
 *    - e(P, Q) can be computed efficiently (polynomial time)
 *
 * From bilinearity, we can derive:
 *    e(aP, bQ) = e(P, Q)^(ab)
 *
 * Proof:
 *    e(aP, bQ) = e(P + P + ... + P, Q + Q + ... + Q)  (a copies, b copies)
 *             = e(P, Q) * e(P, Q) * ... * e(P, Q)     (ab times, by bilinearity)
 *             = e(P, Q)^(ab)
 */

console.log("Bilinearity: e(aP, bQ) = e(P, Q)^(ab)");
console.log("This connects addition in G1/G2 with multiplication in GT.\n");

// Let's demonstrate with a concrete example using a small curve
// Note: We use small parameters for educational purposes.
// In practice, you would use pairing-friendly curves like BN254 or BLS12-381.

const p = 59n;  // Small prime for demonstration
const F = GF(p);

// Create a curve E over F_p
// For pairings, we need specific curve constructions (pairing-friendly curves)
// Here we use a simple example to illustrate the concepts
const E = EllipticCurve(F, [1n, 0n]);  // y^2 = x^3 + x

console.log(`Working over F_${p}`);
console.log(`Curve: ${E.toString()}\n`);

// ============================================================================
// Section 2: Types of Pairings
// ============================================================================

console.log("--- Section 2: Types of Pairings ---\n");

/**
 * SYMMETRIC vs ASYMMETRIC PAIRINGS
 *
 * Type 1 (Symmetric): G1 = G2 = G
 *   - e: G x G -> GT
 *   - Typically on supersingular curves
 *   - Simpler but less efficient
 *
 * Type 2: G1 != G2, but there exists an efficient map phi: G2 -> G1
 *   - e: G1 x G2 -> GT
 *   - Used in some protocols
 *
 * Type 3 (Asymmetric): G1 != G2, no efficient map between them
 *   - e: G1 x G2 -> GT
 *   - Most common in modern protocols
 *   - BLS12-381, BN254 are Type 3
 *
 * The embedding degree k determines:
 *   - GT is a subgroup of F_{q^k}*
 *   - Security: discrete log in GT should be hard
 *   - Efficiency: smaller k means faster pairings but lower security
 */

console.log("Pairing Types:");
console.log("  Type 1: G1 = G2 (symmetric, e.g., supersingular curves)");
console.log("  Type 2: G1 != G2 with map G2 -> G1");
console.log("  Type 3: G1 != G2 with no efficient map (most common)\n");

// ============================================================================
// Section 3: Properties and Their Cryptographic Uses
// ============================================================================

console.log("--- Section 3: Properties and Their Cryptographic Uses ---\n");

/**
 * BILINEARITY ENABLES NEW CRYPTOGRAPHIC PRIMITIVES
 *
 * 1. THREE-PARTY KEY EXCHANGE (Joux's Protocol)
 *    Alice: a, broadcasts aP
 *    Bob:   b, broadcasts bP
 *    Carol: c, broadcasts cP
 *
 *    Shared key = e(aP, bP)^c = e(aP, cP)^b = e(bP, cP)^a = e(P, P)^(abc)
 *
 *    Without pairings, three-party key exchange requires two rounds.
 *    With pairings, it's done in ONE round!
 *
 * 2. BLS SIGNATURES
 *    Sign: sigma = sk * H(m)  where H: message -> G1
 *    Verify: e(sigma, g) == e(H(m), pk)  where pk = sk * g
 *
 *    Why it works:
 *    e(sk * H(m), g) = e(H(m), sk * g) = e(H(m), pk)
 *
 * 3. IDENTITY-BASED ENCRYPTION
 *    Master secret: s
 *    User's private key: s * H(identity)
 *    Anyone can encrypt to "user@example.com" without pre-shared keys!
 */

console.log("Cryptographic Applications of Bilinearity:\n");

console.log("1. Three-Party Key Exchange (Joux):");
console.log("   K = e(P, P)^(abc) computed by each party in ONE round\n");

console.log("2. BLS Signatures:");
console.log("   Sign: sigma = sk * H(m)");
console.log("   Verify: e(sigma, g) ?= e(H(m), pk)");
console.log("   Uses bilinearity: e(sk*H(m), g) = e(H(m), sk*g)\n");

console.log("3. Identity-Based Encryption:");
console.log("   Public key IS the identity (email, name, etc.)");
console.log("   No certificates needed!\n");

// ============================================================================
// Section 4: The DDH Gap Problem
// ============================================================================

console.log("--- Section 4: The DDH Gap Problem ---\n");

/**
 * DECISIONAL DIFFIE-HELLMAN (DDH) AND PAIRINGS
 *
 * DDH Problem: Given (g, g^a, g^b, Z), determine if Z = g^(ab)
 *
 * In groups WITHOUT pairings:
 *   DDH is believed to be hard (basis of ElGamal security)
 *
 * In groups WITH pairings:
 *   DDH is EASY to solve!
 *   Check: e(g^a, g^b) ?= e(g, Z)
 *   If Z = g^(ab): e(g, g)^(ab) = e(g, g^(ab)) = e(g, Z) CHECK
 *   If Z != g^(ab): the check fails
 *
 * This creates a "GAP GROUP":
 *   - CDH (Computational DH) is still hard: can't COMPUTE g^(ab)
 *   - DDH (Decisional DH) is easy: can CHECK if Z = g^(ab)
 *
 * Gap groups enable:
 *   - Signature schemes where verification is efficient
 *   - Encryption schemes with special properties
 *   - Zero-knowledge proofs with efficient verification
 */

console.log("DDH in Pairing Groups:");
console.log("  Given (P, aP, bP, Z), check if Z = abP:");
console.log("  e(aP, bP) ?= e(P, Z)");
console.log("  This makes DDH EASY but CDH remains HARD.\n");

console.log("Gap Group Properties:");
console.log("  - Can VERIFY multiplicative relationships");
console.log("  - Cannot COMPUTE the actual value");
console.log("  - Perfect for signatures: sign is hard, verify is easy\n");

// ============================================================================
// Section 5: Pairing Security Assumptions
// ============================================================================

console.log("--- Section 5: Security Assumptions ---\n");

/**
 * STANDARD ASSUMPTIONS IN PAIRING-BASED CRYPTOGRAPHY
 *
 * 1. Bilinear Diffie-Hellman (BDH):
 *    Given (P, aP, bP, cP), compute e(P, P)^(abc)
 *    This is hard if CDH is hard in G1.
 *
 * 2. Decisional BDH (DBDH):
 *    Given (P, aP, bP, cP, Z), decide if Z = e(P, P)^(abc)
 *    Stronger assumption (implies BDH).
 *
 * 3. q-Strong DH (q-SDH):
 *    Given (g, g^s, g^(s^2), ..., g^(s^q)), compute (c, g^(1/(s+c)))
 *    Used in some signature schemes and SNARKs.
 *
 * 4. Knowledge of Exponent (KoE):
 *    If you can produce (g^a, h^a), you "know" a
 *    Used in SNARKs (Groth16, etc.)
 *
 * SECURITY LEVELS
 *
 * For 128-bit security:
 *   - G1, G2: ~256-bit groups (elliptic curve points)
 *   - GT: ~3072-bit groups (finite field extension)
 *
 * Common choices:
 *   - BN254: 254-bit, embedding degree 12, ~100-bit security
 *   - BLS12-381: 381-bit, embedding degree 12, ~128-bit security
 *   - BLS12-377: 377-bit, embedding degree 12, used in SNARKs
 */

console.log("Security Assumptions:");
console.log("  BDH: Given P, aP, bP, cP, compute e(P,P)^(abc)");
console.log("  DBDH: Decide if Z = e(P,P)^(abc)");
console.log("  q-SDH: Given powers of s, compute (c, g^(1/(s+c)))\n");

console.log("Practical Curves for ~128-bit Security:");
console.log("  BLS12-381: 381-bit, k=12 (Ethereum 2.0, Zcash Sapling)");
console.log("  BN254: 254-bit, k=12 (Ethereum contracts, older)");
console.log("  BLS12-377: 377-bit, k=12 (Aleo, SNARK-friendly)\n");

// ============================================================================
// Section 6: Demonstration with Code
// ============================================================================

console.log("--- Section 6: Bilinearity Demonstration ---\n");

// Find the curve order to locate torsion points
const curveOrder = E.cardinality();
console.log(`Curve order: ${curveOrder}`);

// For the Weil pairing, we need n-torsion points
// Let's find a small torsion subgroup
// We'll use a divisor n of the curve order

// Factor the order to find a suitable n
// For demonstration, we'll look for small order points
const n = 5n;  // We'll look for 5-torsion points if they exist

// Note: The Weil pairing is defined on E[n] x E[n] -> mu_n
// where E[n] is the n-torsion subgroup and mu_n are nth roots of unity

// For a complete demonstration, we would:
// 1. Find two independent n-torsion points P and Q
// 2. Compute e(aP, bQ)
// 3. Verify that e(aP, bQ) = e(P, Q)^(ab)

console.log("\nBilinearity property:");
console.log("  e(aP, bQ) = e(P, Q)^(ab)");
console.log("  where a and b are any integers");
console.log("\nThis enables:");
console.log("  - Checking if discrete logs are equal");
console.log("  - Verifying polynomial relationships");
console.log("  - Building SNARKs and other zero-knowledge systems\n");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");
console.log("Bilinear maps provide a bridge between:");
console.log("  - Addition in elliptic curve groups (G1, G2)");
console.log("  - Multiplication in finite field extensions (GT)\n");

console.log("Key property: e(aP, bQ) = e(P, Q)^(ab)");
console.log("This enables checking multiplicative relationships");
console.log("without revealing the underlying discrete logarithms.\n");

console.log("Applications:");
console.log("  - BLS signatures (short, aggregatable)");
console.log("  - Identity-based encryption");
console.log("  - Tripartite key exchange");
console.log("  - SNARKs (Groth16, PLONK)");
console.log("  - Attribute-based encryption\n");

console.log("Next: Learn about specific pairing constructions (Weil, Tate)");
