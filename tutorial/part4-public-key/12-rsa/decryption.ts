/**
 * # RSA Decryption
 *
 * This tutorial covers RSA decryption: recovering plaintext from ciphertext
 * using the private key.
 *
 * ## The Decryption Function
 *
 * Given a private key (n, d) and ciphertext c:
 *   m = c^d mod n
 *
 * ## Why Decryption Works
 *
 * Since e * d ≡ 1 (mod phi(n)), we have e * d = 1 + k * phi(n).
 * By Euler's theorem, for gcd(m, n) = 1:
 *   m^(phi(n)) ≡ 1 (mod n)
 *
 * Therefore:
 *   c^d = (m^e)^d = m^(e*d) = m^(1 + k*phi(n)) = m * (m^phi(n))^k ≡ m * 1 ≡ m (mod n)
 *
 * ## Performance Considerations
 *
 * Decryption (c^d) is much slower than encryption (m^e) because d is large.
 * The CRT optimization (covered in crt-optimization.ts) speeds up decryption.
 */

import { gcd, inverse_mod, power_mod, crt, Mod } from 'sagemath-ts';

console.log('='.repeat(60));
console.log('RSA Decryption Tutorial');
console.log('='.repeat(60));

// ============================================================================
// Setup: RSA Key Pair
// ============================================================================

console.log('\n=== RSA Key Setup ===\n');

const p = 61n;
const q = 53n;
const n = p * q;
const phi_n = (p - 1n) * (q - 1n);
const e = 17n;
const d = inverse_mod(e, phi_n);

console.log('Key Parameters:');
console.log(`  p = ${p}, q = ${q}`);
console.log(`  n = ${n}`);
console.log(`  phi(n) = ${phi_n}`);
console.log(`  e = ${e}`);
console.log(`  d = ${d}`);

// Verify e * d ≡ 1 (mod phi(n))
console.log(`\nVerification: e * d mod phi(n) = ${(e * d) % phi_n}`);

// ============================================================================
// Basic RSA Decryption
// ============================================================================

console.log('\n=== Basic RSA Decryption ===\n');

/**
 * RSA Decryption:
 *   m = c^d mod n
 *
 * This reverses encryption (c = m^e mod n) by the mathematical
 * relationship between e and d.
 */

const originalMessage = 123n;
console.log(`Original message: m = ${originalMessage}`);

// Encrypt
const ciphertext = Mod(originalMessage, n).pow(e).value;
console.log(`Ciphertext: c = m^e mod n = ${ciphertext}`);

// Decrypt
const decrypted = Mod(ciphertext, n).pow(d).value;
console.log(`Decrypted: m' = c^d mod n = ${decrypted}`);
console.log(`Recovery successful: ${decrypted === originalMessage}`);

// ============================================================================
// Step-by-Step Decryption
// ============================================================================

console.log('\n=== Step-by-Step Decryption ===\n');

/**
 * Let's trace through the decryption calculation to understand
 * the modular exponentiation process.
 */

const c = 855n; // A ciphertext
console.log(`Ciphertext c = ${c}`);
console.log(`Private exponent d = ${d}`);
console.log(`Modulus n = ${n}`);

console.log(`\nComputing c^d mod n using square-and-multiply:`);

// Show binary representation of d
const dBinary = d.toString(2);
console.log(`d in binary: ${dBinary}`);
console.log(`d has ${dBinary.length} bits`);

// The square-and-multiply algorithm processes d from most significant bit
let result = 1n;
let base = c;
let dCopy = d;
const steps: string[] = [];

for (let i = 0; dCopy > 0n; i++) {
  const bit = dCopy & 1n;

  if (bit === 1n) {
    const oldResult = result;
    result = (result * base) % n;
    steps.push(`Bit ${i} = 1: ${oldResult} * ${base} mod ${n} = ${result}`);
  }

  dCopy = dCopy >> 1n;

  if (dCopy > 0n) {
    const oldBase = base;
    base = (base * base) % n;
    // Note: We don't log every squaring to keep output manageable
  }
}

console.log(`\nFirst few multiplication steps:`);
steps.slice(0, 5).forEach((step) => console.log(`  ${step}`));
console.log(`  ... (${steps.length - 5} more steps)`);

console.log(`\nFinal result: ${result}`);

// Verify with library
const libraryResult = Mod(c, n).pow(d).value;
console.log(`Library result: ${libraryResult}`);
console.log(`Match: ${result === libraryResult}`);

// ============================================================================
// Why Decryption is Slow
// ============================================================================

console.log('\n=== Why Decryption is Slow ===\n');

/**
 * The private exponent d is typically very large (close to n in size).
 * This is because:
 *   d = e^(-1) mod phi(n)
 *
 * With e = 65537 (small), d ends up being approximately the same
 * size as phi(n), which is close to n.
 *
 * Number of modular multiplications:
 * - Encryption: ~17 for e = 65537 (only 2 one-bits)
 * - Decryption: ~bit_length(d) = ~bit_length(n) operations
 */

console.log('Encryption vs Decryption complexity:');
console.log(`  e = ${e} has ${e.toString(2).length} bits`);
console.log(`  d = ${d} has ${d.toString(2).length} bits`);
console.log('');
console.log(`  Encryption uses ~${e.toString(2).split('').filter((b) => b === '1').length + e.toString(2).length - 1} modular operations`);
console.log(`  Decryption uses ~${d.toString(2).length * 2} modular operations`);
console.log('');
console.log('  Decryption is roughly 100x slower than encryption!');

// With realistic key sizes
console.log('\nWith 2048-bit RSA:');
console.log('  e = 65537: ~17 modular operations');
console.log('  d (2048 bits): ~3000 modular operations');
console.log('  Plus, each operation is on 2048-bit numbers');

// ============================================================================
// Decryption with Different Message Sizes
// ============================================================================

console.log('\n=== Decrypting Various Messages ===\n');

const testMessages = [1n, 42n, 100n, 1000n, 3000n, n - 1n];

console.log('Testing decryption for various message values:');
console.log('-'.repeat(50));

for (const m of testMessages) {
  if (m >= n) {
    console.log(`m = ${m}: Too large (>= n = ${n})`);
    continue;
  }

  const c_test = Mod(m, n).pow(e).value;
  const m_recovered = Mod(c_test, n).pow(d).value;

  console.log(`m = ${m}:`);
  console.log(`  Encrypted: c = ${c_test}`);
  console.log(`  Decrypted: m' = ${m_recovered}`);
  console.log(`  Success: ${m === m_recovered}`);
}

// ============================================================================
// Handling Special Cases
// ============================================================================

console.log('\n=== Special Cases in Decryption ===\n');

/**
 * Edge cases to consider:
 * 1. m = 0: c = 0, decrypts to 0
 * 2. m = 1: c = 1, decrypts to 1
 * 3. gcd(m, n) > 1: Still works, but reveals information!
 */

console.log('Case 1: m = 0');
let m_special = 0n;
let c_special = Mod(m_special, n).pow(e).value;
let m_recovered = Mod(c_special, n).pow(d).value;
console.log(`  c = ${c_special}, decrypted = ${m_recovered}`);

console.log('\nCase 2: m = 1');
m_special = 1n;
c_special = Mod(m_special, n).pow(e).value;
m_recovered = Mod(c_special, n).pow(d).value;
console.log(`  c = ${c_special}, decrypted = ${m_recovered}`);

console.log('\nCase 3: gcd(m, n) > 1 (DANGEROUS!)');
// m that shares a factor with n
m_special = p; // m = 61 = p
c_special = Mod(m_special, n).pow(e).value;
m_recovered = Mod(c_special, n).pow(d).value;
console.log(`  m = ${m_special} = p`);
console.log(`  gcd(m, n) = ${gcd(m_special, n)}`);
console.log(`  c = ${c_special}`);
console.log(`  decrypted = ${m_recovered}`);
console.log(`  Decryption still works!`);

/**
 * WARNING: If an attacker can get you to decrypt a message where
 * gcd(m, n) > 1, they can factor n!
 *
 * If they know m shares a factor with n:
 *   gcd(m, n) = p or q
 *   This completely breaks the key!
 */

console.log('\n  SECURITY NOTE:');
console.log('  If attacker knows m where gcd(m, n) > 1:');
console.log(`    gcd(${m_special}, ${n}) = ${gcd(m_special, n)} = p!`);
console.log('  This would completely break the RSA key!');

// ============================================================================
// Decryption Without Knowing the Factorization
// ============================================================================

console.log('\n=== Why Attackers Cannot Decrypt ===\n');

/**
 * To decrypt, you need d.
 * To compute d = e^(-1) mod phi(n), you need phi(n).
 * To compute phi(n) = (p-1)(q-1), you need to factor n.
 *
 * The security of RSA rests on the difficulty of factoring n.
 */

console.log('Attack scenarios and why they fail:');
console.log('');

console.log('1. Compute d directly?');
console.log(`   d = e^(-1) mod phi(n)`);
console.log(`   But phi(n) = ${phi_n} requires knowing p and q`);
console.log('');

console.log('2. Factor n to find p and q?');
console.log(`   n = ${n} = ${p} * ${q}`);
console.log('   For small n, this is easy by trial division');
console.log('   For 2048-bit n, factoring takes ~2^112 operations');
console.log('');

console.log('3. Compute phi(n) without factoring?');
console.log('   Knowing phi(n) is equivalent to knowing the factorization');
console.log('   From n and phi(n):');
console.log(`     p + q = n - phi(n) + 1 = ${n} - ${phi_n} + 1 = ${n - phi_n + 1n}`);
console.log(`     p * q = n = ${n}`);
console.log('   Solving: p and q are roots of x^2 - (p+q)x + pq = 0');
console.log('');

console.log('4. Other attacks require:');
console.log('   - Side-channel information (timing, power analysis)');
console.log('   - Chosen ciphertext attacks (padding oracle)');
console.log('   - Implementation flaws');

// ============================================================================
// Partial Decryption Information Leakage
// ============================================================================

console.log('\n=== Information Leakage Considerations ===\n');

/**
 * Even partial information about decryption can be dangerous:
 *
 * 1. Knowing whether m is even or odd (1 bit of m)
 * 2. Knowing m's most significant bit
 * 3. Knowing if decryption "succeeded" (padding oracle)
 *
 * These can lead to attacks that recover the full message.
 */

console.log('Potential information leaks:');
console.log('');

console.log('1. Parity Oracle Attack:');
console.log('   If attacker learns whether m is even/odd after decryption,');
console.log('   they can recover m bit by bit using the malleability of RSA.');
console.log('');

console.log('2. Padding Oracle Attack (Bleichenbacher):');
console.log('   If server reveals whether PKCS#1 v1.5 padding is valid,');
console.log('   attacker can decrypt arbitrary ciphertexts.');
console.log('   This attack requires ~1 million queries.');
console.log('');

console.log('3. Timing Attacks:');
console.log('   If decryption time depends on the plaintext,');
console.log('   attacker can deduce information about d or m.');
console.log('   Solution: constant-time implementations.');

// ============================================================================
// Correct Decryption Implementation
// ============================================================================

console.log('\n=== Secure Decryption Practices ===\n');

/**
 * A production RSA decryption implementation should:
 * 1. Use constant-time modular exponentiation
 * 2. Validate padding before returning success/failure
 * 3. Use blinding to prevent timing attacks
 * 4. Return only valid/invalid, never partial information
 */

console.log('Production decryption checklist:');
console.log('');
console.log('  [x] Constant-time modular exponentiation');
console.log('  [x] RSA blinding (multiply by random r^e before decryption)');
console.log('  [x] Verify OAEP padding structure');
console.log('  [x] Return only success/failure');
console.log('  [x] Use CRT optimization for speed');
console.log('  [x] Clear sensitive values from memory');

// Demonstrate RSA blinding concept
console.log('\nRSA Blinding concept:');
console.log('  1. Choose random r, compute r_e = r^e mod n');
console.log('  2. Blind ciphertext: c_blind = c * r_e mod n');
console.log('  3. Decrypt: m_blind = c_blind^d mod n');
console.log('  4. Unblind: m = m_blind * r^(-1) mod n');
console.log('');
console.log('  This prevents timing attacks because the value being');
console.log('  exponentiated (c * r^e) is randomized.');

// ============================================================================
// Summary
// ============================================================================

console.log('\n=== Summary ===\n');

console.log('RSA DECRYPTION KEY POINTS:');
console.log('');
console.log('  1. Operation: m = c^d mod n');
console.log('');
console.log('  2. Why it works:');
console.log('     e * d ≡ 1 (mod phi(n))');
console.log('     (m^e)^d = m^(ed) ≡ m (mod n)');
console.log('');
console.log('  3. Performance:');
console.log('     - Much slower than encryption (d is large)');
console.log('     - CRT optimization provides ~4x speedup');
console.log('');
console.log('  4. Security:');
console.log('     - Never leak partial decryption information');
console.log('     - Use constant-time implementations');
console.log('     - Use blinding to prevent timing attacks');
console.log('     - Validate padding to prevent oracle attacks');

console.log('\n' + '='.repeat(60));
console.log('RSA Decryption Tutorial Complete');
console.log('='.repeat(60));
