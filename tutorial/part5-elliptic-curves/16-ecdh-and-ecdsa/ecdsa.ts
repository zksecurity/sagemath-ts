/**
 * # Elliptic Curve Digital Signature Algorithm (ECDSA)
 *
 * ECDSA is the elliptic curve variant of DSA, providing digital signatures
 * with much smaller key sizes than RSA for equivalent security.
 *
 * ## What is a Digital Signature?
 *
 * A digital signature provides:
 * - **Authentication**: Verifies the signer's identity
 * - **Integrity**: Detects if the message was modified
 * - **Non-repudiation**: Signer cannot deny having signed
 *
 * ## ECDSA Overview
 *
 * - Key generation: Private key d, Public key Q = dG
 * - Signing: Uses private key and message hash to create (r, s)
 * - Verification: Uses public key and message hash to verify (r, s)
 *
 * ## Security
 *
 * Security is based on:
 * - Hardness of ECDLP (recovering d from Q = dG)
 * - Quality of random number k used in signing
 *
 * @module tutorial/elliptic-curves/ecdsa
 */

import { EllipticCurve, EllipticCurvePoint, FiniteFieldPrime, inverse_mod, gcd } from 'sagemath-ts';

// ============================================================================
// ECDSA Setup
// ============================================================================

console.log("=== ECDSA: Elliptic Curve Digital Signature Algorithm ===\n");

// Setup: Choose a curve and base point
// Using a curve with PRIME order (crucial for ECDSA!)
// E: y^2 = x^3 + 3x + 2 over F_127 has prime order 139
const F127 = new FiniteFieldPrime(127n);
const E = EllipticCurve(F127, [3n, 2n]);

// Find a generator
const G = E.gens()[0]!;
const n = G.order(); // Order of base point (should be prime)

console.log("Domain Parameters:");
console.log(`  Curve: ${E.toString()}`);
console.log(`  Base point G: ${G.toString()}`);
console.log(`  Order n: ${n}`);
console.log(`  (In practice, n would be a 256-bit prime)`);
console.log("");

// ============================================================================
// Key Generation
// ============================================================================

console.log("=== Key Generation ===\n");

/**
 * Key Generation:
 * 1. Pick random d in [1, n-1]
 * 2. Compute Q = dG
 * 3. Private key: d
 * 4. Public key: Q
 */

function generateKeyPair(G: EllipticCurvePoint, n: bigint): { d: bigint; Q: EllipticCurvePoint } {
  // Random private key in [1, n-1]
  const d = BigInt(Math.floor(Math.random() * Number(n - 1n))) + 1n;
  const Q = G.mul(d);
  return { d, Q };
}

const { d: privateKey, Q: publicKey } = generateKeyPair(G, n);

console.log("Key Pair Generated:");
console.log(`  Private key (d): ${privateKey} (KEEP SECRET!)`);
console.log(`  Public key (Q):  ${publicKey.toString()}`);
console.log("");

// ============================================================================
// Signing Algorithm
// ============================================================================

console.log("=== ECDSA Signing ===\n");

/**
 * ECDSA Signing Algorithm:
 *
 * Input: Private key d, message hash e (as integer)
 * Output: Signature (r, s)
 *
 * 1. Pick random k in [1, n-1] (CRITICAL: must be unique and secret!)
 * 2. Compute R = kG
 * 3. Let r = x_R mod n. If r = 0, go to step 1.
 * 4. Compute s = k^(-1) * (e + d*r) mod n. If s = 0, go to step 1.
 * 5. Return (r, s)
 */

function sign(
  message: bigint,
  privateKey: bigint,
  G: EllipticCurvePoint,
  n: bigint
): { r: bigint; s: bigint; k: bigint } {
  let r = 0n;
  let s = 0n;
  let k = 0n;
  let R: EllipticCurvePoint;

  // Keep trying until we get valid r and s
  while (r === 0n || s === 0n) {
    // Step 1: Pick random k that is coprime to n (required for inverse)
    do {
      k = BigInt(Math.floor(Math.random() * Number(n - 1n))) + 1n;
    } while (gcd(k, n) !== 1n);

    // Step 2: Compute R = kG
    R = G.mul(k);

    // Step 3: r = x_R mod n
    if (!R.isZero() && R.x !== null) {
      r = R.x.value % n;
      if (r < 0n) r += n;
    }

    if (r === 0n) continue;

    // Step 4: s = k^(-1) * (e + d*r) mod n
    const kInv = inverse_mod(k, n);
    s = (kInv * (message + privateKey * r)) % n;
    if (s < 0n) s += n;
  }

  return { r, s, k };
}

// Example: Sign a message
// In practice, message would be hash(actual_message) truncated to fit in n
const messageHash = 42n; // Simulated hash of "Hello, World!"

console.log(`Message hash (e): ${messageHash}`);
console.log(`(In practice: e = SHA-256(message) truncated to bit-length of n)`);
console.log("");

const { r, s, k } = sign(messageHash, privateKey, G, n);

console.log("Signing Process:");
console.log(`  1. Generate random k: ${k} (MUST be unique and secret!)`);
console.log(`  2. Compute R = kG: ${G.mul(k).toString()}`);
console.log(`  3. r = x_R mod n: ${r}`);
console.log(`  4. s = k^(-1) * (e + d*r) mod n: ${s}`);
console.log("");
console.log(`Signature: (r, s) = (${r}, ${s})`);
console.log("");

// ============================================================================
// Verification Algorithm
// ============================================================================

console.log("=== ECDSA Verification ===\n");

/**
 * ECDSA Verification Algorithm:
 *
 * Input: Public key Q, message hash e, signature (r, s)
 * Output: VALID or INVALID
 *
 * 1. Verify r, s are in [1, n-1]
 * 2. Compute w = s^(-1) mod n
 * 3. Compute u1 = e*w mod n and u2 = r*w mod n
 * 4. Compute R' = u1*G + u2*Q
 * 5. If R' = O, reject
 * 6. Let v = x_{R'} mod n
 * 7. Accept if v = r, reject otherwise
 */

function verify(
  message: bigint,
  signature: { r: bigint; s: bigint },
  publicKey: EllipticCurvePoint,
  G: EllipticCurvePoint,
  n: bigint
): boolean {
  const { r, s } = signature;

  // Step 1: Check range
  if (r < 1n || r >= n || s < 1n || s >= n) {
    console.log("  Verification failed: r or s out of range");
    return false;
  }

  // Step 2: w = s^(-1) mod n
  const w = inverse_mod(s, n);

  // Step 3: u1 = e*w mod n, u2 = r*w mod n
  const u1 = (message * w) % n;
  const u2 = (r * w) % n;

  // Step 4: R' = u1*G + u2*Q
  const R_prime = G.mul(u1).add(publicKey.mul(u2));

  // Step 5: Check R' != O
  if (R_prime.isZero()) {
    console.log("  Verification failed: R' is point at infinity");
    return false;
  }

  // Step 6: v = x_{R'} mod n
  let v = R_prime.x!.value % n;
  if (v < 0n) v += n;

  // Step 7: Accept if v = r
  console.log(`  w = s^(-1) mod n = ${w}`);
  console.log(`  u1 = e*w mod n = ${u1}`);
  console.log(`  u2 = r*w mod n = ${u2}`);
  console.log(`  R' = u1*G + u2*Q = ${R_prime.toString()}`);
  console.log(`  v = x_{R'} mod n = ${v}`);
  console.log(`  r = ${r}`);
  console.log(`  v == r? ${v === r}`);

  return v === r;
}

console.log("Verification Process:");
const isValid = verify(messageHash, { r, s }, publicKey, G, n);
console.log("");
console.log(`Signature is ${isValid ? "VALID" : "INVALID"}`);
console.log("");

// ============================================================================
// Why Verification Works
// ============================================================================

console.log("=== Why Verification Works ===\n");

/**
 * Mathematical justification:
 *
 * We have: s = k^(-1) * (e + d*r) mod n
 * So: k = s^(-1) * (e + d*r) mod n
 *     k = e*s^(-1) + d*r*s^(-1) mod n
 *     k = e*w + d*r*w mod n
 *     k = u1 + d*u2 mod n
 *
 * Therefore:
 *   R' = u1*G + u2*Q
 *      = u1*G + u2*(d*G)
 *      = (u1 + d*u2)*G
 *      = k*G
 *      = R
 *
 * So x_{R'} = x_R = r (mod n)
 */

console.log("Why R' = R:");
console.log(`  Given: s = k^(-1) * (e + d*r) mod n`);
console.log(`  We get: k = e*s^(-1) + d*r*s^(-1) = u1 + d*u2 mod n`);
console.log(`  So: R' = u1*G + u2*Q = u1*G + u2*d*G = (u1 + d*u2)*G = k*G = R`);
console.log(`  Therefore: x_{R'} = x_R, and v = r`);
console.log("");

// ============================================================================
// Security Critical: The Nonce k
// ============================================================================

console.log("=== CRITICAL: The Nonce k ===\n");

/**
 * THE MOST IMPORTANT SECURITY RULE IN ECDSA:
 *
 * The nonce k MUST be:
 * 1. Randomly generated
 * 2. Unique for every signature
 * 3. Kept completely secret
 *
 * If k is reused or leaked, THE PRIVATE KEY IS COMPROMISED!
 *
 * Attack if k is known:
 *   s = k^(-1) * (e + d*r) mod n
 *   s*k = e + d*r mod n
 *   d = (s*k - e) * r^(-1) mod n
 *
 * Attack if k is reused on two messages e1, e2:
 *   s1 = k^(-1) * (e1 + d*r) mod n
 *   s2 = k^(-1) * (e2 + d*r) mod n
 *   s1 - s2 = k^(-1) * (e1 - e2) mod n
 *   k = (e1 - e2) * (s1 - s2)^(-1) mod n
 *   Then recover d as above.
 *
 * Real-world examples:
 * - PlayStation 3 hack (2010): Sony used same k for all signatures
 * - Bitcoin thefts: Poor random number generators leaked k
 */

console.log("THE k VALUE IS CRITICAL FOR SECURITY!");
console.log("");
console.log("If k is known, private key can be recovered:");
console.log(`  d = (s*k - e) * r^(-1) mod n`);

// Demonstrate the attack
const d_recovered = (((s * k - messageHash) % n) * inverse_mod(r, n)) % n;
const d_normalized = d_recovered < 0n ? d_recovered + n : d_recovered;

console.log(`  With k = ${k}:`);
console.log(`    d_recovered = (${s} * ${k} - ${messageHash}) * ${r}^(-1) mod ${n}`);
console.log(`    d_recovered = ${d_normalized}`);
console.log(`    Actual d    = ${privateKey}`);
console.log(`    Attack worked? ${d_normalized === privateKey}`);
console.log("");

console.log("If k is reused on two messages:");
console.log("  k = (e1 - e2) * (s1 - s2)^(-1) mod n");
console.log("  Then private key can be recovered as above.");
console.log("");

console.log("MITIGATION: Use RFC 6979 deterministic nonce generation:");
console.log("  k = HMAC-DRBG(private_key, message_hash)");
console.log("  This eliminates random number generator failures.");
console.log("");

// ============================================================================
// Complete ECDSA Example
// ============================================================================

console.log("=== Complete ECDSA Protocol ===\n");

function ecdsaDemo() {
  // Key generation
  console.log("1. KEY GENERATION (done once):");
  const { d, Q } = generateKeyPair(G, n);
  console.log(`   Private key d = ${d}`);
  console.log(`   Public key Q = ${Q.toString()}`);
  console.log("");

  // Message
  const message = "Transfer 100 BTC to Alice";
  const messageHashDemo = 57n; // In practice: SHA-256(message)
  console.log(`2. MESSAGE: "${message}"`);
  console.log(`   Hash e = ${messageHashDemo} (simulated)`);
  console.log("");

  // Sign
  console.log("3. SIGNING (uses private key d):");
  const sig = sign(messageHashDemo, d, G, n);
  console.log(`   Signature: (r, s) = (${sig.r}, ${sig.s})`);
  console.log("");

  // Transmit (message, signature, public key)
  console.log("4. TRANSMISSION:");
  console.log(`   Send: message, (r=${sig.r}, s=${sig.s}), Q=${Q.toString()}`);
  console.log("");

  // Verify
  console.log("5. VERIFICATION (uses public key Q):");
  const valid = verify(messageHashDemo, sig, Q, G, n);
  console.log(`   Result: ${valid ? "ACCEPT" : "REJECT"}`);
  console.log("");

  // Tampered message
  console.log("6. TAMPER DETECTION:");
  const tamperedHash = 58n; // Different message
  console.log(`   Verifying with different message hash: ${tamperedHash}`);
  const tamperedValid = verify(tamperedHash, sig, Q, G, n);
  console.log(`   Result: ${tamperedValid ? "ACCEPT" : "REJECT (as expected)"}`);

  return valid;
}

ecdsaDemo();
console.log("");

// ============================================================================
// ECDSA in Practice
// ============================================================================

console.log("=== ECDSA in Practice ===\n");

console.log("Real-world usage:");
console.log("  - Bitcoin and Ethereum transactions (secp256k1)");
console.log("  - TLS certificates (P-256, P-384)");
console.log("  - Code signing (iOS, Android apps)");
console.log("  - Secure boot verification");
console.log("");

console.log("Standard curves for ECDSA:");
console.log("  - secp256k1: Used by Bitcoin, 256-bit");
console.log("  - P-256 (secp256r1): NIST, widely used in TLS");
console.log("  - P-384 (secp384r1): NIST, higher security");
console.log("  - Ed25519: Modern curve, used in EdDSA (different algorithm)");
console.log("");

console.log("Best practices:");
console.log("  1. Use RFC 6979 for deterministic k generation");
console.log("  2. Use a well-audited cryptographic library");
console.log("  3. Never implement ECDSA from scratch in production");
console.log("  4. Hash messages before signing (don't sign raw data)");
console.log("  5. Verify public key is on the curve before verification");
console.log("");

// ============================================================================
// Exercises
// ============================================================================

console.log("=== Exercises ===\n");

console.log("1. Implement ECDSA verification that rejects invalid signatures:");
console.log("   - Wrong message hash");
console.log("   - Wrong public key");
console.log("   - Corrupted r or s");
console.log("");

console.log("2. Demonstrate the k-reuse attack:");
console.log("   - Sign two different messages with the same k");
console.log("   - Recover the private key from the two signatures");
console.log("");

console.log("3. Compare ECDSA signature size with RSA-2048:");
console.log("   - ECDSA with P-256: signature is 64 bytes");
console.log("   - RSA-2048: signature is 256 bytes");
console.log("   - Both provide ~128-bit security");
console.log("");
