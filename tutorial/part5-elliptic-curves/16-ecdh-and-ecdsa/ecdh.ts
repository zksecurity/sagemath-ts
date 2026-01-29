/**
 * # Elliptic Curve Diffie-Hellman (ECDH)
 *
 * ECDH is a key agreement protocol that allows two parties to establish
 * a shared secret over an insecure channel. It is the elliptic curve
 * analog of the classic Diffie-Hellman protocol.
 *
 * ## The Protocol
 *
 * Setup: A curve E over F_p, a base point G of prime order n.
 *
 * 1. Alice picks random a in [1, n-1], computes A = aG (public key)
 * 2. Bob picks random b in [1, n-1], computes B = bG (public key)
 * 3. Alice and Bob exchange A and B publicly
 * 4. Alice computes S = a * B = a * (bG) = abG
 * 5. Bob computes S = b * A = b * (aG) = abG
 *
 * Both arrive at the same shared secret S = abG!
 *
 * ## Security
 *
 * An eavesdropper sees G, A = aG, B = bG but cannot compute S = abG.
 * This is the Elliptic Curve Computational Diffie-Hellman (ECCDH) problem,
 * which is believed to be as hard as ECDLP.
 *
 * @module tutorial/elliptic-curves/ecdh
 */

import { EllipticCurve, EllipticCurvePoint, FiniteFieldPrime } from 'sagemath-ts';

// ============================================================================
// ECDH Key Exchange
// ============================================================================

console.log("=== Elliptic Curve Diffie-Hellman (ECDH) ===\n");

// Setup: Choose a curve and base point
// Using a small prime for demonstration (real curves use 256+ bit primes)
const F101 = new FiniteFieldPrime(101n);
const E = EllipticCurve(F101, [2n, 3n]);

// Find a generator of large prime order
const curveOrder = E.cardinality();
console.log(`Curve: ${E.toString()}`);
console.log(`Curve order: ${curveOrder}`);

// Get a generator
const G = E.gens()[0]!;
const n = G.order();

console.log(`Base point G: ${G.toString()}`);
console.log(`Order of G: ${n}`);
console.log("");

// ============================================================================
// Key Generation
// ============================================================================

console.log("=== Key Generation ===\n");

/**
 * Each party generates a key pair:
 * - Private key: random integer a in [1, n-1]
 * - Public key: A = aG
 */

// Helper to generate random private key
function randomPrivateKey(order: bigint): bigint {
  // In practice, use cryptographically secure random numbers
  return BigInt(Math.floor(Math.random() * Number(order - 1n))) + 1n;
}

// Alice generates her key pair
const alicePrivate = randomPrivateKey(n);
const alicePublic = G.mul(alicePrivate);

console.log("Alice's key pair:");
console.log(`  Private key (a): ${alicePrivate} (SECRET - never share!)`);
console.log(`  Public key (A): ${alicePublic.toString()}`);
console.log("");

// Bob generates his key pair
const bobPrivate = randomPrivateKey(n);
const bobPublic = G.mul(bobPrivate);

console.log("Bob's key pair:");
console.log(`  Private key (b): ${bobPrivate} (SECRET - never share!)`);
console.log(`  Public key (B): ${bobPublic.toString()}`);
console.log("");

// ============================================================================
// Key Exchange
// ============================================================================

console.log("=== Key Exchange ===\n");

/**
 * Public information (visible to eavesdroppers):
 * - Curve E, base point G, order n
 * - Alice's public key A = aG
 * - Bob's public key B = bG
 *
 * Secret information:
 * - Alice's private key a
 * - Bob's private key b
 * - Shared secret S = abG
 */

console.log("Public information (visible to Eve):");
console.log(`  Curve: y^2 = x^3 + 2x + 3 over F_101`);
console.log(`  Base point G: ${G.toString()}`);
console.log(`  Alice's public key A: ${alicePublic.toString()}`);
console.log(`  Bob's public key B: ${bobPublic.toString()}`);
console.log("");

// ============================================================================
// Shared Secret Computation
// ============================================================================

console.log("=== Shared Secret Computation ===\n");

// Alice computes: S = a * B
const aliceSharedSecret = bobPublic.mul(alicePrivate);
console.log("Alice computes:");
console.log(`  S = a * B = ${alicePrivate} * ${bobPublic.toString()}`);
console.log(`  S = ${aliceSharedSecret.toString()}`);
console.log("");

// Bob computes: S = b * A
const bobSharedSecret = alicePublic.mul(bobPrivate);
console.log("Bob computes:");
console.log(`  S = b * A = ${bobPrivate} * ${alicePublic.toString()}`);
console.log(`  S = ${bobSharedSecret.toString()}`);
console.log("");

// Verify they match
console.log("Verification:");
console.log(`  Alice's S: ${aliceSharedSecret.toString()}`);
console.log(`  Bob's S:   ${bobSharedSecret.toString()}`);
console.log(`  Match? ${aliceSharedSecret.eq(bobSharedSecret)}`);
console.log("");

// Explain why they match
console.log("Why they match:");
console.log(`  Alice: a * B = a * (bG) = (ab)G`);
console.log(`  Bob:   b * A = b * (aG) = (ba)G = (ab)G`);
console.log(`  Scalar multiplication is commutative!`);
console.log("");

// ============================================================================
// Key Derivation
// ============================================================================

console.log("=== Key Derivation ===\n");

/**
 * The shared secret S is a point on the curve.
 * To use it as a symmetric key, we need to derive bytes from it.
 *
 * Common approaches:
 * 1. Use the x-coordinate: key = Hash(x-coordinate of S)
 * 2. Use both coordinates: key = Hash(x || y)
 * 3. Use a KDF: key = HKDF(S, info, length)
 *
 * Standards like ECIES, X25519, ECDH in TLS specify exact methods.
 */

if (!aliceSharedSecret.isZero()) {
  const sharedX = aliceSharedSecret.x!;
  const sharedY = aliceSharedSecret.y!;

  console.log("Deriving symmetric key from shared secret:");
  console.log(`  S = (${sharedX}, ${sharedY})`);
  console.log(`  x-coordinate: ${sharedX}`);
  console.log("");
  console.log("  In practice, apply a hash function:");
  console.log("    key = SHA-256(x-coordinate) or");
  console.log("    key = HKDF(x || y, salt, info, length)");
}
console.log("");

// ============================================================================
// Security Analysis
// ============================================================================

console.log("=== Security Analysis ===\n");

/**
 * The ECDH Problem:
 * Given G, A = aG, B = bG, compute S = abG.
 *
 * This is the Computational Diffie-Hellman (CDH) problem on elliptic curves.
 *
 * Reductions:
 * - CDH <= ECDLP: If you can solve ECDLP, recover a from (G, aG),
 *   then compute S = a * B.
 *
 * It's unknown whether CDH = ECDLP, but they're believed equivalent.
 *
 * Best attacks:
 * - Generic: Pollard rho, O(sqrt(n)) time
 * - No index calculus attack for properly chosen curves!
 */

console.log("What an eavesdropper (Eve) sees:");
console.log(`  G = ${G.toString()}`);
console.log(`  A = ${alicePublic.toString()}`);
console.log(`  B = ${bobPublic.toString()}`);
console.log("");

console.log("To compute S = abG, Eve would need to:");
console.log("  1. Solve ECDLP: Find a given G and A = aG, OR");
console.log("  2. Solve CDH directly: Compute abG from G, aG, bG");
console.log("");

console.log("Security assumptions:");
console.log("  - ECDLP is hard (no efficient algorithm known)");
console.log("  - CDH is hard (believed equivalent to ECDLP)");
console.log("  - For 256-bit curves: ~128-bit security");
console.log("");

// ============================================================================
// Man-in-the-Middle Attack
// ============================================================================

console.log("=== Man-in-the-Middle Attack ===\n");

/**
 * ECDH alone does NOT provide authentication!
 *
 * A man-in-the-middle (Mallory) can:
 * 1. Intercept Alice's A, send own M_A to Bob
 * 2. Intercept Bob's B, send own M_B to Alice
 * 3. Establish separate keys with Alice and Bob
 * 4. Relay and read all messages
 *
 * Solution: Authenticated Key Exchange
 * - ECDH + Digital Signatures (e.g., TLS)
 * - ECDH + MACs with pre-shared key
 * - ECDH + Identity-based cryptography
 */

console.log("WARNING: Basic ECDH is vulnerable to man-in-the-middle attacks!");
console.log("");
console.log("Attack scenario:");
console.log("  1. Mallory intercepts Alice's A, generates own keys (m, M=mG)");
console.log("  2. Mallory sends M to Bob (pretending to be Alice)");
console.log("  3. Mallory intercepts Bob's B, sends own M' to Alice");
console.log("  4. Now: Alice-Mallory share key S1, Bob-Mallory share key S2");
console.log("  5. Mallory can decrypt, read, re-encrypt all messages");
console.log("");
console.log("Solutions:");
console.log("  - Use authenticated ECDH (sign public keys)");
console.log("  - Use ECDHE + certificates (like TLS)");
console.log("  - Use authenticated encryption after key exchange");
console.log("");

// ============================================================================
// Practical Considerations
// ============================================================================

console.log("=== Practical Considerations ===\n");

console.log("1. PRIVATE KEY GENERATION");
console.log("   - Use cryptographically secure random number generator");
console.log("   - Ensure uniform distribution in [1, n-1]");
console.log("   - Never reuse private keys in ephemeral ECDHE");
console.log("");

console.log("2. PUBLIC KEY VALIDATION");
console.log("   - Check received point is on the curve");
console.log("   - Check point is not the identity O");
console.log("   - For curves with cofactor h > 1, multiply by h");
console.log("   - Prevents small subgroup attacks");
console.log("");

console.log("3. EPHEMERAL vs STATIC");
console.log("   - ECDH: Static keys, used multiple times");
console.log("   - ECDHE: Ephemeral keys, fresh for each exchange");
console.log("   - ECDHE provides forward secrecy");
console.log("");

// ============================================================================
// Complete ECDH Example
// ============================================================================

console.log("=== Complete ECDH Protocol ===\n");

function ecdhProtocol() {
  console.log("Protocol execution:");
  console.log("");

  // Setup (public parameters)
  console.log("SETUP (public parameters):");
  console.log(`  Curve: E over F_${F101.characteristic}`);
  console.log(`  Base point: G = ${G.toString()}`);
  console.log(`  Order: n = ${n}`);
  console.log("");

  // Alice
  const a = randomPrivateKey(n);
  const A = G.mul(a);
  console.log("ALICE generates keys:");
  console.log(`  a = ${a} (private, kept secret)`);
  console.log(`  A = aG = ${A.toString()} (public, sent to Bob)`);
  console.log("");

  // Bob
  const b = randomPrivateKey(n);
  const B = G.mul(b);
  console.log("BOB generates keys:");
  console.log(`  b = ${b} (private, kept secret)`);
  console.log(`  B = bG = ${B.toString()} (public, sent to Alice)`);
  console.log("");

  // Exchange
  console.log("KEY EXCHANGE:");
  console.log(`  Alice receives B = ${B.toString()}`);
  console.log(`  Bob receives A = ${A.toString()}`);
  console.log("");

  // Compute shared secrets
  const S_alice = B.mul(a);
  const S_bob = A.mul(b);

  console.log("SHARED SECRET COMPUTATION:");
  console.log(`  Alice: S = a*B = ${S_alice.toString()}`);
  console.log(`  Bob:   S = b*A = ${S_bob.toString()}`);
  console.log("");

  console.log("RESULT:");
  console.log(`  Shared secret established: ${S_alice.eq(S_bob)}`);
  if (!S_alice.isZero()) {
    console.log(`  S = ${S_alice.toString()}`);
    console.log(`  Derive symmetric key from x-coordinate: ${S_alice.x}`);
  }

  return { a, A, b, B, S: S_alice };
}

ecdhProtocol();
console.log("");

// ============================================================================
// Exercises
// ============================================================================

console.log("=== Exercises ===\n");

console.log("1. Implement ECDH over the curve y^2 = x^3 + x + 1 over F_257");
console.log("   Generate key pairs for Alice and Bob, compute shared secret");
console.log("");

console.log("2. Demonstrate a man-in-the-middle attack:");
console.log("   - Mallory generates two key pairs");
console.log("   - Show that Alice-Mallory and Bob-Mallory have different secrets");
console.log("");

console.log("3. Implement public key validation:");
console.log("   - Check point is on curve");
console.log("   - Check point is not identity");
console.log("   - Check point has correct order (for cofactor > 1)");
console.log("");
