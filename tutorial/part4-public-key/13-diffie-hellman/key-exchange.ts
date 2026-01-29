/**
 * # Diffie-Hellman Key Exchange
 *
 * The Diffie-Hellman (DH) key exchange protocol allows two parties to establish
 * a shared secret over an insecure channel without prior shared secrets.
 *
 * ## Historical Significance
 *
 * Published by Whitfield Diffie and Martin Hellman in 1976, this was the first
 * practical method for establishing a shared key over an insecure channel.
 * It revolutionized cryptography and laid the foundation for public-key crypto.
 *
 * ## Protocol Overview
 *
 * 1. Alice and Bob agree on public parameters: prime p and generator g
 * 2. Alice picks secret a, sends A = g^a mod p
 * 3. Bob picks secret b, sends B = g^b mod p
 * 4. Alice computes K = B^a mod p = g^(ab) mod p
 * 5. Bob computes K = A^b mod p = g^(ab) mod p
 * 6. Both have the same shared secret K!
 *
 * ## Security Basis
 *
 * Security relies on the Discrete Logarithm Problem (DLP):
 * Given g, p, and g^x mod p, it's hard to find x.
 */

import { gcd, is_prime, power_mod, factor, euler_phi, inverse_mod, Mod, GF } from 'sagemath-ts';

console.log('='.repeat(60));
console.log('Diffie-Hellman Key Exchange Tutorial');
console.log('='.repeat(60));

// ============================================================================
// Step 1: Public Parameters
// ============================================================================

console.log('\n=== Step 1: Public Parameters ===\n');

/**
 * The public parameters are:
 * - A large prime p (the modulus)
 * - A generator g of the multiplicative group (Z/pZ)*
 *
 * These are shared publicly and reused across many key exchanges.
 */

// For educational purposes, use a small prime
// In production: p should be at least 2048 bits
const p = 23n;

console.log(`Prime modulus: p = ${p}`);
console.log(`This is MUCH smaller than real-world parameters!`);
console.log(`Production systems use 2048+ bit primes.`);

// Verify p is prime
console.log(`\nVerification: is_prime(${p}) = ${is_prime(p)}`);

// The multiplicative group (Z/pZ)* has order p-1
const groupOrder = p - 1n;
console.log(`\nMultiplicative group order: phi(p) = p - 1 = ${groupOrder}`);
console.log(`Factorization of p-1: ${p - 1n} = ${factor(p - 1n).map(([f, e]) => e === 1n ? `${f}` : `${f}^${e}`).join(' * ')}`);

// ============================================================================
// Step 2: Finding a Generator
// ============================================================================

console.log('\n=== Step 2: Finding a Generator ===\n');

/**
 * A generator g of (Z/pZ)* satisfies:
 * - {g^1, g^2, ..., g^(p-1)} = {1, 2, ..., p-1} mod p
 * - The order of g is exactly p-1
 *
 * To find a generator:
 * 1. For each candidate g in {2, 3, ...}
 * 2. Check that g^((p-1)/q) != 1 for all prime factors q of p-1
 * 3. If all checks pass, g is a generator
 */

function findGenerator(p: bigint): bigint {
  const order = p - 1n;
  const factors = factor(order).filter(([f]) => f > 0n);

  for (let g = 2n; g < p; g++) {
    let isGenerator = true;

    for (const [q] of factors) {
      const exp = order / q;
      if (power_mod(g, exp, p) === 1n) {
        isGenerator = false;
        break;
      }
    }

    if (isGenerator) {
      return g;
    }
  }

  throw new Error('No generator found');
}

const g = findGenerator(p);
console.log(`Generator found: g = ${g}`);

// Verify it's a generator by computing its order
function computeOrder(g: bigint, p: bigint): bigint {
  const maxOrder = p - 1n;
  let order = 1n;
  let current = g;

  while (current !== 1n && order <= maxOrder) {
    current = (current * g) % p;
    order++;
  }

  return order;
}

const gOrder = computeOrder(g, p);
console.log(`Order of g: ord(${g}) = ${gOrder}`);
console.log(`Is generator? ${gOrder === p - 1n}`);

// Show powers of g
console.log(`\nPowers of g = ${g}:`);
let current = 1n;
const powers: string[] = [];
for (let i = 0n; i <= p - 1n; i++) {
  powers.push(`${g}^${i} = ${current}`);
  current = (current * g) % p;
}
console.log(powers.slice(0, 12).join(', ') + (powers.length > 12 ? ', ...' : ''));

// ============================================================================
// Step 3: The Key Exchange Protocol
// ============================================================================

console.log('\n=== Step 3: The Key Exchange Protocol ===\n');

/**
 * DIFFIE-HELLMAN KEY EXCHANGE
 *
 * Public parameters: p (prime), g (generator)
 *
 * Alice:                          Bob:
 * ------                          ----
 * Choose secret a                 Choose secret b
 * Compute A = g^a mod p           Compute B = g^b mod p
 *
 *              Alice sends A to Bob
 *              Bob sends B to Alice
 *
 * Compute K = B^a mod p           Compute K = A^b mod p
 *         K = g^(ab) mod p               K = g^(ab) mod p
 *
 * Both now share secret K!
 */

console.log('Protocol Execution:');
console.log('-------------------');

// Alice's side
const a = 6n; // Alice's private key
console.log(`\nAlice chooses secret: a = ${a}`);

const A = Mod(g, p).pow(a).value;
console.log(`Alice computes: A = g^a mod p = ${g}^${a} mod ${p} = ${A}`);
console.log(`Alice sends A = ${A} to Bob (publicly)`);

// Bob's side
const b = 15n; // Bob's private key
console.log(`\nBob chooses secret: b = ${b}`);

const B = Mod(g, p).pow(b).value;
console.log(`Bob computes: B = g^b mod p = ${g}^${b} mod ${p} = ${B}`);
console.log(`Bob sends B = ${B} to Alice (publicly)`);

// Shared secret computation
console.log('\n--- Computing Shared Secret ---');

const K_alice = Mod(B, p).pow(a).value;
console.log(`\nAlice computes: K = B^a mod p = ${B}^${a} mod ${p} = ${K_alice}`);

const K_bob = Mod(A, p).pow(b).value;
console.log(`Bob computes: K = A^b mod p = ${A}^${b} mod ${p} = ${K_bob}`);

console.log(`\nShared secret: K = ${K_alice}`);
console.log(`Both computed same value: ${K_alice === K_bob}`);

// Verify: both equal g^(ab)
const K_verify = Mod(g, p).pow(a * b).value;
console.log(`\nVerification: g^(ab) = ${g}^${a * b} mod ${p} = ${K_verify}`);

// ============================================================================
// Why It Works: The Mathematics
// ============================================================================

console.log('\n=== Why It Works ===\n');

console.log('Mathematical explanation:');
console.log('');
console.log('  Alice receives B = g^b mod p');
console.log('  Alice computes: B^a = (g^b)^a = g^(ba) mod p');
console.log('');
console.log('  Bob receives A = g^a mod p');
console.log('  Bob computes: A^b = (g^a)^b = g^(ab) mod p');
console.log('');
console.log('  Since ab = ba, both compute g^(ab) mod p');
console.log('');
console.log(`  In our example: g^(ab) = ${g}^(${a}*${b}) = ${g}^${a * b} mod ${p} = ${K_alice}`);

// ============================================================================
// What an Eavesdropper Sees
// ============================================================================

console.log('\n=== What an Eavesdropper Sees ===\n');

/**
 * An eavesdropper (Eve) observes:
 * - Public parameters: p, g
 * - Alice's public value: A = g^a mod p
 * - Bob's public value: B = g^b mod p
 *
 * To compute K = g^(ab), Eve would need to find a from A (or b from B).
 * This is the Discrete Logarithm Problem (DLP).
 */

console.log('Eve observes:');
console.log(`  p = ${p} (public parameter)`);
console.log(`  g = ${g} (public parameter)`);
console.log(`  A = ${A} (from Alice)`);
console.log(`  B = ${B} (from Bob)`);

console.log('\nTo compute K, Eve needs to solve:');
console.log(`  Find a such that ${g}^a = ${A} (mod ${p})`);
console.log('  OR');
console.log(`  Find b such that ${g}^b = ${B} (mod ${p})`);

console.log('\nThis is the Discrete Logarithm Problem (DLP).');
console.log('For small p, it\'s easy. For 2048-bit p, it\'s infeasible.');

// For our small example, we can solve it by brute force
console.log('\nBrute force solution (only works for small p):');
for (let x = 0n; x < p; x++) {
  if (Mod(g, p).pow(x).value === A) {
    console.log(`  Found: ${g}^${x} = ${A} (mod ${p}), so a = ${x}`);
    break;
  }
}

// ============================================================================
// Computational Diffie-Hellman Problem
// ============================================================================

console.log('\n=== The CDH Problem ===\n');

/**
 * The Computational Diffie-Hellman (CDH) Problem:
 *
 * Given: g, p, A = g^a, B = g^b
 * Compute: g^(ab) mod p
 *
 * Note: CDH doesn't require finding a or b!
 * It only asks for the shared secret.
 *
 * CDH is believed to be as hard as DLP, but this is not proven.
 */

console.log('Computational Diffie-Hellman (CDH) Problem:');
console.log('');
console.log('  Given: g, p, A = g^a mod p, B = g^b mod p');
console.log('  Compute: g^(ab) mod p (without knowing a or b)');
console.log('');
console.log('  CDH security is slightly weaker than DLP:');
console.log('  - DLP: Given A, find a');
console.log('  - CDH: Given A and B, find g^(ab)');
console.log('');
console.log('  If you can solve DLP, you can solve CDH.');
console.log('  Reverse implication is unknown!');

// ============================================================================
// Multiple Key Exchanges
// ============================================================================

console.log('\n=== Multiple Key Exchanges ===\n');

/**
 * DH parameters can be reused. Each party keeps one private key
 * and can exchange keys with multiple parties.
 */

console.log('Alice can establish keys with multiple parties:');
console.log('');

// Alice keeps the same a
console.log(`Alice's private key: a = ${a}`);
console.log(`Alice's public key: A = ${A}`);

// Bob (original)
console.log(`\nWith Bob (b = ${b}):`);
console.log(`  Shared secret: ${K_alice}`);

// Carol
const c_carol = 9n;
const C_carol = Mod(g, p).pow(c_carol).value;
const K_carol = Mod(C_carol, p).pow(a).value;
console.log(`\nWith Carol (c = ${c_carol}):`);
console.log(`  Carol's public key: C = ${C_carol}`);
console.log(`  Shared secret: ${K_carol}`);

// David
const d_david = 18n;
const D_david = Mod(g, p).pow(d_david).value;
const K_david = Mod(D_david, p).pow(a).value;
console.log(`\nWith David (d = ${d_david}):`);
console.log(`  David's public key: D = ${D_david}`);
console.log(`  Shared secret: ${K_david}`);

// ============================================================================
// Ephemeral vs Static Keys
// ============================================================================

console.log('\n=== Ephemeral vs Static Keys ===\n');

/**
 * Static DH: Same key pair used for many exchanges
 * - Pro: Can pre-compute public key
 * - Con: If private key is compromised, all past sessions are compromised
 *
 * Ephemeral DH (DHE): New key pair for each exchange
 * - Pro: Forward secrecy (past sessions safe even if key leaked)
 * - Con: Must generate new keys each time
 *
 * Best practice: Use ephemeral keys (DHE, ECDHE)
 */

console.log('Key Types:');
console.log('');
console.log('  Static DH:');
console.log('    - Same private key for many exchanges');
console.log('    - Risk: compromise reveals all past sessions');
console.log('');
console.log('  Ephemeral DH (DHE):');
console.log('    - New private key for each session');
console.log('    - Forward secrecy: past sessions stay secure');
console.log('    - Recommended for TLS and modern protocols');
console.log('');

console.log('Forward Secrecy Example:');

// Simulate ephemeral key exchange
console.log('\n  Session 1:');
const a1 = 3n,
  b1 = 7n;
const A1 = Mod(g, p).pow(a1).value;
const B1 = Mod(g, p).pow(b1).value;
const K1 = Mod(B1, p).pow(a1).value;
console.log(`    Ephemeral keys: a1=${a1}, b1=${b1}`);
console.log(`    Shared secret: K1=${K1}`);
console.log(`    (Keys destroyed after session)`);

console.log('\n  Session 2:');
const a2 = 11n,
  b2 = 13n;
const A2 = Mod(g, p).pow(a2).value;
const B2 = Mod(g, p).pow(b2).value;
const K2 = Mod(B2, p).pow(a2).value;
console.log(`    Ephemeral keys: a2=${a2}, b2=${b2}`);
console.log(`    Shared secret: K2=${K2}`);
console.log(`    (Completely independent of Session 1)`);

// ============================================================================
// Deriving Encryption Keys
// ============================================================================

console.log('\n=== Deriving Encryption Keys ===\n');

/**
 * The shared secret K should NOT be used directly as an encryption key.
 * Instead, derive keys using a Key Derivation Function (KDF).
 *
 * Common KDFs:
 * - HKDF (HMAC-based KDF)
 * - SHA-256(K || context)
 * - PBKDF2 (for password-based scenarios)
 */

console.log('Key Derivation (conceptual):');
console.log('');
console.log('  Shared secret: K = g^(ab) mod p');
console.log('');
console.log('  Encryption key: enc_key = KDF(K, "encryption")');
console.log('  MAC key:        mac_key = KDF(K, "authentication")');
console.log('  IV/Nonce:       iv = KDF(K, "nonce")');
console.log('');
console.log('  Using raw K as encryption key is dangerous!');
console.log('  - K has algebraic structure');
console.log('  - K may not be uniformly distributed');
console.log('  - KDF removes these weaknesses');

// Simulate simple key derivation
function simpleKDF(K: bigint, context: string): bigint {
  // THIS IS NOT A SECURE KDF - just for demonstration
  let hash = K;
  for (let i = 0; i < context.length; i++) {
    hash = (hash * 31n + BigInt(context.charCodeAt(i))) % (2n ** 128n);
  }
  return hash;
}

console.log('\nSimulated key derivation:');
const encKey = simpleKDF(K_alice, 'encryption');
const macKey = simpleKDF(K_alice, 'authentication');
console.log(`  Shared secret K = ${K_alice}`);
console.log(`  Encryption key = ${encKey}`);
console.log(`  MAC key = ${macKey}`);

// ============================================================================
// Summary
// ============================================================================

console.log('\n=== Summary ===\n');

console.log('DIFFIE-HELLMAN KEY EXCHANGE:');
console.log('');
console.log('  Protocol:');
console.log('    1. Agree on public p (prime) and g (generator)');
console.log('    2. Alice: secret a, sends A = g^a mod p');
console.log('    3. Bob: secret b, sends B = g^b mod p');
console.log('    4. Both compute: K = g^(ab) mod p');
console.log('');
console.log('  Security:');
console.log('    - Based on Discrete Logarithm Problem');
console.log('    - Parameters must be large (2048+ bits)');
console.log('    - Use ephemeral keys for forward secrecy');
console.log('');
console.log('  Best Practices:');
console.log('    - Use standard parameters (RFC 3526, etc.)');
console.log('    - Derive keys with KDF, don\'t use K directly');
console.log('    - Consider ECDH for better performance');

console.log('\n' + '='.repeat(60));
console.log('Diffie-Hellman Key Exchange Tutorial Complete');
console.log('='.repeat(60));
