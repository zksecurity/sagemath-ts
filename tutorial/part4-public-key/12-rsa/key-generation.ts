/**
 * # RSA Key Generation
 *
 * RSA (Rivest-Shamir-Adleman) is one of the first practical public-key cryptosystems
 * and is widely used for secure data transmission. This tutorial walks through the
 * mathematical foundations and implementation of RSA key generation.
 *
 * ## Mathematical Foundation
 *
 * RSA security relies on the computational difficulty of factoring the product of
 * two large prime numbers (the "factoring problem").
 *
 * ## Key Generation Algorithm
 *
 * 1. Choose two distinct large prime numbers p and q
 * 2. Compute n = p * q (the modulus)
 * 3. Compute phi(n) = (p-1)(q-1) (Euler's totient)
 * 4. Choose e such that 1 < e < phi(n) and gcd(e, phi(n)) = 1
 * 5. Compute d = e^(-1) mod phi(n) (the private exponent)
 *
 * Public key: (n, e)
 * Private key: (n, d) or equivalently (p, q, d)
 *
 * ## Security Considerations
 *
 * - p and q should be of similar bit length but not too close
 * - |p - q| should not be too small (prevents Fermat factorization)
 * - p - 1 and q - 1 should have large prime factors (prevents Pollard's p-1)
 * - e = 65537 (0x10001) is commonly used as it's prime and has low Hamming weight
 */

import {
  gcd,
  is_prime,
  euler_phi,
  inverse_mod,
  next_prime,
  factor,
  formatFactorization,
} from 'sagemath-ts';
import { Mod } from 'sagemath-ts';

console.log('='.repeat(60));
console.log('RSA Key Generation Tutorial');
console.log('='.repeat(60));

// ============================================================================
// Step 1: Choose Two Large Primes
// ============================================================================

console.log('\n=== Step 1: Choose Two Large Primes ===\n');

/**
 * In practice, we use cryptographically secure random number generators
 * to find primes of 1024+ bits. For educational purposes, we use smaller primes.
 *
 * Key insight: The primes should be "random" and of similar size, but not too close.
 * If |p - q| is small, Fermat's factorization can quickly find p and q.
 */

// For demonstration, we use small primes
// In production: use 512-bit or larger primes for 1024-bit RSA (minimum recommended today: 2048-bit)
const p = 61n;
const q = 53n;

console.log(`Selected primes:`);
console.log(`  p = ${p}`);
console.log(`  q = ${q}`);

// Verify they are prime
console.log(`\nVerification:`);
console.log(`  is_prime(${p}) = ${is_prime(p)}`);
console.log(`  is_prime(${q}) = ${is_prime(q)}`);

// Security check: ensure p and q are not too close
const difference = p > q ? p - q : q - p;
console.log(`  |p - q| = ${difference}`);

// ============================================================================
// Step 2: Compute the Modulus n
// ============================================================================

console.log('\n=== Step 2: Compute Modulus n = p * q ===\n');

const n = p * q;
console.log(`n = p * q = ${p} * ${q} = ${n}`);

/**
 * The modulus n is public. Its security relies on the difficulty of
 * factoring it back into p and q.
 *
 * Bit length of n determines the security level:
 * - 1024 bits: ~80 bits of security (deprecated)
 * - 2048 bits: ~112 bits of security (current minimum)
 * - 3072 bits: ~128 bits of security
 * - 4096 bits: ~140 bits of security
 */

const bitLength = n.toString(2).length;
console.log(`Bit length of n: ${bitLength} bits`);

// ============================================================================
// Step 3: Compute Euler's Totient phi(n)
// ============================================================================

console.log("\n=== Step 3: Compute Euler's Totient phi(n) ===\n");

/**
 * For n = p * q where p, q are distinct primes:
 *   phi(n) = (p - 1)(q - 1)
 *
 * This is because:
 * - phi(p) = p - 1 for prime p
 * - phi(q) = q - 1 for prime q
 * - phi(p*q) = phi(p) * phi(q) for coprime p, q
 *
 * phi(n) counts integers in [1, n-1] that are coprime to n.
 * These are exactly the elements that have multiplicative inverses mod n.
 */

const phi_n = (p - 1n) * (q - 1n);
console.log(`phi(n) = (p-1)(q-1) = (${p}-1)(${q}-1) = ${p - 1n} * ${q - 1n} = ${phi_n}`);

// Verify using the library function
const phi_n_direct = euler_phi(n);
console.log(`euler_phi(${n}) = ${phi_n_direct} (verification)`);

/**
 * IMPORTANT: phi(n) must be kept secret!
 * If an attacker knows phi(n), they can compute:
 *   d = e^(-1) mod phi(n)
 * which breaks the cryptosystem.
 *
 * In fact, knowing phi(n) is equivalent to knowing the factorization:
 * Given n and phi(n), we can solve for p and q:
 *   p + q = n - phi(n) + 1
 *   p * q = n
 * This is a quadratic equation with solutions p and q.
 */

// Demonstrate recovering p and q from n and phi(n)
console.log(`\nRecovering primes from n and phi(n):`);
const sum_pq = n - phi_n + 1n;
console.log(`  p + q = n - phi(n) + 1 = ${sum_pq}`);
console.log(`  p * q = n = ${n}`);
// p and q are roots of: x^2 - (p+q)x + pq = 0

// ============================================================================
// Step 4: Choose Public Exponent e
// ============================================================================

console.log('\n=== Step 4: Choose Public Exponent e ===\n');

/**
 * The public exponent e must satisfy:
 * 1. 1 < e < phi(n)
 * 2. gcd(e, phi(n)) = 1 (so e has an inverse mod phi(n))
 *
 * Common choices:
 * - e = 3: Fast encryption but vulnerable to some attacks
 * - e = 17: Another small prime
 * - e = 65537 (2^16 + 1): Standard choice, good balance
 *
 * The value 65537 is popular because:
 * - It's prime, so likely coprime to phi(n)
 * - It has only two 1-bits (Hamming weight 2), making exponentiation fast
 * - It's large enough to resist small exponent attacks
 */

// Common choice: e = 65537
let e = 65537n;

// Check if e is coprime to phi(n)
if (gcd(e, phi_n) !== 1n) {
  // If not coprime, find next suitable e
  console.log(`gcd(${e}, ${phi_n}) != 1, finding alternative...`);
  e = 3n;
  while (gcd(e, phi_n) !== 1n) {
    e += 2n;
  }
}

console.log(`Public exponent e = ${e}`);
console.log(`gcd(e, phi(n)) = gcd(${e}, ${phi_n}) = ${gcd(e, phi_n)}`);

// For our small example, let's use e = 17 which is more instructive
e = 17n;
console.log(`\nFor this tutorial, using e = ${e}`);
console.log(`gcd(${e}, ${phi_n}) = ${gcd(e, phi_n)}`);

// ============================================================================
// Step 5: Compute Private Exponent d
// ============================================================================

console.log('\n=== Step 5: Compute Private Exponent d ===\n');

/**
 * The private exponent d is the modular inverse of e modulo phi(n):
 *   d * e ≡ 1 (mod phi(n))
 *
 * This means: d * e = 1 + k * phi(n) for some integer k
 *
 * d is computed using the Extended Euclidean Algorithm.
 */

const d = inverse_mod(e, phi_n);
console.log(`d = e^(-1) mod phi(n)`);
console.log(`d = ${e}^(-1) mod ${phi_n}`);
console.log(`d = ${d}`);

// Verify: d * e mod phi(n) should be 1
const verification = (d * e) % phi_n;
console.log(`\nVerification: d * e mod phi(n) = ${d} * ${e} mod ${phi_n} = ${verification}`);

// ============================================================================
// Summary: The RSA Key Pair
// ============================================================================

console.log('\n=== RSA Key Pair Summary ===\n');

console.log('PUBLIC KEY (can be shared):');
console.log(`  Modulus n = ${n}`);
console.log(`  Public exponent e = ${e}`);

console.log('\nPRIVATE KEY (must be kept secret):');
console.log(`  Modulus n = ${n}`);
console.log(`  Private exponent d = ${d}`);

console.log('\nSECRET VALUES (destroy after key generation):');
console.log(`  Prime p = ${p}`);
console.log(`  Prime q = ${q}`);
console.log(`  phi(n) = ${phi_n}`);

// ============================================================================
// Understanding Why RSA Works
// ============================================================================

console.log('\n=== Why RSA Works: The Mathematical Foundation ===\n');

/**
 * RSA relies on Euler's theorem (a generalization of Fermat's little theorem):
 *
 * For any a coprime to n:
 *   a^phi(n) ≡ 1 (mod n)
 *
 * Since e * d ≡ 1 (mod phi(n)), we have:
 *   e * d = 1 + k * phi(n) for some integer k
 *
 * Therefore, for any message m coprime to n:
 *   (m^e)^d = m^(e*d) = m^(1 + k*phi(n)) = m * (m^phi(n))^k ≡ m * 1^k ≡ m (mod n)
 *
 * Similarly:
 *   (m^d)^e ≡ m (mod n)
 *
 * This shows that encryption and decryption are inverse operations.
 */

console.log("Euler's theorem: For gcd(a, n) = 1:");
console.log(`  a^phi(n) ≡ 1 (mod n)`);
console.log(`  a^${phi_n} ≡ 1 (mod ${n})`);

// Demonstrate with an example
const testMessage = 42n;
const testMod = Mod(testMessage, n);
const eulerResult = testMod.pow(phi_n);
console.log(`\nExample: ${testMessage}^${phi_n} mod ${n} = ${eulerResult.value}`);

// Show the encryption-decryption cycle
console.log('\nEncryption-Decryption cycle:');
const m = 42n;
const mMod = Mod(m, n);
const encrypted = mMod.pow(e);
const decrypted = encrypted.pow(d);
console.log(`  Original message m = ${m}`);
console.log(`  Encrypted: m^e mod n = ${m}^${e} mod ${n} = ${encrypted.value}`);
console.log(`  Decrypted: (m^e)^d mod n = ${encrypted.value}^${d} mod ${n} = ${decrypted.value}`);

// ============================================================================
// Larger Example with More Realistic Primes
// ============================================================================

console.log('\n=== Larger Example ===\n');

// Use larger (but still educational) primes
const p2 = 104729n; // 100,000th prime
const q2 = 104723n; // Another prime

console.log(`Larger primes:`);
console.log(`  p = ${p2}`);
console.log(`  q = ${q2}`);

const n2 = p2 * q2;
const phi_n2 = (p2 - 1n) * (q2 - 1n);
const e2 = 65537n;
const d2 = inverse_mod(e2, phi_n2);

console.log(`\nKey parameters:`);
console.log(`  n = ${n2}`);
console.log(`  phi(n) = ${phi_n2}`);
console.log(`  e = ${e2}`);
console.log(`  d = ${d2}`);

// Test encryption/decryption
const message2 = 123456789n;
const c2 = Mod(message2, n2).pow(e2).value;
const m2_recovered = Mod(c2, n2).pow(d2).value;

console.log(`\nEncryption test:`);
console.log(`  Message: ${message2}`);
console.log(`  Ciphertext: ${c2}`);
console.log(`  Decrypted: ${m2_recovered}`);
console.log(`  Success: ${message2 === m2_recovered}`);

// ============================================================================
// Key Generation Best Practices
// ============================================================================

console.log('\n=== Key Generation Best Practices ===\n');

console.log('1. PRIME SELECTION:');
console.log('   - Use cryptographically secure random number generator');
console.log('   - Primes should be "provable" primes or "probable" primes with high confidence');
console.log('   - |p - q| should be large (at least 2^(n/2 - 100) for n-bit modulus)');
console.log('   - p - 1 and q - 1 should have large prime factors');

console.log('\n2. KEY SIZES (as of 2024):');
console.log('   - 2048 bits: Minimum for general use');
console.log('   - 3072 bits: Recommended for medium-term security');
console.log('   - 4096 bits: Recommended for long-term security');

console.log('\n3. PUBLIC EXPONENT:');
console.log('   - e = 65537 is the standard choice');
console.log('   - Avoid e = 3 without proper padding');

console.log('\n4. PRIVATE KEY STORAGE:');
console.log('   - Store in secure hardware when possible (HSM, TPM)');
console.log('   - Use key encryption for storage');
console.log('   - Consider CRT representation for faster decryption');

console.log('\n' + '='.repeat(60));
console.log('RSA Key Generation Tutorial Complete');
console.log('='.repeat(60));
