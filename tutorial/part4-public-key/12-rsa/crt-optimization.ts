/**
 * # RSA CRT Optimization
 *
 * The Chinese Remainder Theorem (CRT) can speed up RSA decryption by
 * approximately 4x. This is the standard optimization used in all
 * production RSA implementations.
 *
 * ## The Basic Idea
 *
 * Instead of computing m = c^d mod n directly (expensive for large n),
 * we compute:
 *   m_p = c^(d mod (p-1)) mod p
 *   m_q = c^(d mod (q-1)) mod q
 *
 * Then combine using CRT to get m mod n.
 *
 * ## Why It's Faster
 *
 * - Operations mod p and mod q use half the bits of n
 * - Modular exponentiation is O(k^3) for k-bit numbers
 * - Two half-size operations: 2 * (k/2)^3 = k^3/4
 * - Plus small CRT overhead
 * - Net speedup: ~4x
 *
 * ## Security Considerations
 *
 * CRT decryption can be vulnerable to fault attacks.
 * Always verify the result after CRT computation.
 */

import { gcd, inverse_mod, power_mod, crt, Mod } from 'sagemath-ts';

console.log('='.repeat(60));
console.log('RSA CRT Optimization Tutorial');
console.log('='.repeat(60));

// ============================================================================
// Setup: RSA Key Pair with CRT Components
// ============================================================================

console.log('\n=== RSA Key Setup with CRT Components ===\n');

const p = 61n;
const q = 53n;
const n = p * q;
const phi_n = (p - 1n) * (q - 1n);
const e = 17n;
const d = inverse_mod(e, phi_n);

console.log('Standard RSA parameters:');
console.log(`  p = ${p}`);
console.log(`  q = ${q}`);
console.log(`  n = p * q = ${n}`);
console.log(`  phi(n) = ${phi_n}`);
console.log(`  e = ${e}`);
console.log(`  d = ${d}`);

// CRT components
const d_p = d % (p - 1n); // d mod (p-1)
const d_q = d % (q - 1n); // d mod (q-1)
const q_inv = inverse_mod(q, p); // q^(-1) mod p

console.log('\nCRT components:');
console.log(`  d_p = d mod (p-1) = ${d} mod ${p - 1n} = ${d_p}`);
console.log(`  d_q = d mod (q-1) = ${d} mod ${q - 1n} = ${d_q}`);
console.log(`  q_inv = q^(-1) mod p = ${q}^(-1) mod ${p} = ${q_inv}`);

// Verify q_inv
const qInvCheck = (q * q_inv) % p;
console.log(`  Verification: q * q_inv mod p = ${q} * ${q_inv} mod ${p} = ${qInvCheck}`);

// ============================================================================
// Standard RSA Decryption (Slow)
// ============================================================================

console.log('\n=== Standard RSA Decryption ===\n');

const message = 123n;
const ciphertext = Mod(message, n).pow(e).value;

console.log(`Message: m = ${message}`);
console.log(`Ciphertext: c = m^e mod n = ${ciphertext}`);

// Standard decryption
const startStandard = performance.now();
let standard_count = 0;

// Simulate multiple decryptions
for (let i = 0; i < 100; i++) {
  const m_standard = Mod(ciphertext, n).pow(d).value;
  standard_count++;
}
const endStandard = performance.now();

const m_standard = Mod(ciphertext, n).pow(d).value;
console.log(`\nStandard decryption: c^d mod n`);
console.log(`  c^d mod n = ${ciphertext}^${d} mod ${n}`);
console.log(`  Result: ${m_standard}`);

// ============================================================================
// CRT Decryption (Fast)
// ============================================================================

console.log('\n=== CRT RSA Decryption ===\n');

/**
 * CRT Decryption Algorithm:
 *
 * 1. Compute m_p = c^(d_p) mod p
 * 2. Compute m_q = c^(d_q) mod q
 * 3. Combine: m = m_q + q * ((m_p - m_q) * q_inv mod p)
 *
 * This is called "Garner's formula" for CRT.
 */

console.log('Step 1: Compute m_p = c^d_p mod p');
const c_mod_p = ciphertext % p;
const m_p = Mod(c_mod_p, p).pow(d_p).value;
console.log(`  c mod p = ${c_mod_p}`);
console.log(`  m_p = ${c_mod_p}^${d_p} mod ${p} = ${m_p}`);

console.log('\nStep 2: Compute m_q = c^d_q mod q');
const c_mod_q = ciphertext % q;
const m_q = Mod(c_mod_q, q).pow(d_q).value;
console.log(`  c mod q = ${c_mod_q}`);
console.log(`  m_q = ${c_mod_q}^${d_q} mod ${q} = ${m_q}`);

console.log("\nStep 3: Combine using Garner's formula");
// h = (m_p - m_q) * q_inv mod p
let h = ((m_p - m_q) % p + p) % p; // Ensure positive
h = (h * q_inv) % p;
console.log(`  h = (m_p - m_q) * q_inv mod p`);
console.log(`  h = (${m_p} - ${m_q}) * ${q_inv} mod ${p} = ${h}`);

// m = m_q + q * h
const m_crt = m_q + q * h;
console.log(`  m = m_q + q * h = ${m_q} + ${q} * ${h} = ${m_crt}`);

console.log(`\nVerification:`);
console.log(`  Standard result: ${m_standard}`);
console.log(`  CRT result:      ${m_crt}`);
console.log(`  Match: ${m_standard === m_crt}`);

// ============================================================================
// Performance Comparison
// ============================================================================

console.log('\n=== Performance Analysis ===\n');

/**
 * Why CRT is faster:
 *
 * For n-bit RSA modulus:
 * - Standard: One exponentiation with ~n-bit exponent mod n-bit number
 * - CRT: Two exponentiations with ~(n/2)-bit exponents mod (n/2)-bit numbers
 *
 * Modular exponentiation cost is roughly O(k * k^2) = O(k^3) where k is bit size.
 *
 * Standard: O(n^3)
 * CRT: 2 * O((n/2)^3) = 2 * O(n^3/8) = O(n^3/4)
 *
 * Theoretical speedup: 4x
 * Practical speedup: ~3-4x (due to CRT combination overhead)
 */

console.log('Theoretical complexity analysis:');
console.log('');
console.log('  Let n have k bits.');
console.log('');
console.log('  Standard decryption:');
console.log('    - Exponent d: ~k bits');
console.log('    - Modulus n: k bits');
console.log('    - Cost: O(k^3)');
console.log('');
console.log('  CRT decryption:');
console.log('    - Exponent d_p: ~k/2 bits');
console.log('    - Modulus p: ~k/2 bits');
console.log('    - Cost per modulus: O((k/2)^3) = O(k^3/8)');
console.log('    - Two moduli: 2 * O(k^3/8) = O(k^3/4)');
console.log('    - Plus small CRT overhead');
console.log('');
console.log('  Speedup: k^3 / (k^3/4) = 4x');

// Count operations for our example
const bits_d = d.toString(2).length;
const bits_d_p = d_p.toString(2).length;
const bits_d_q = d_q.toString(2).length;

console.log(`\nOur example (bits in exponents):`);
console.log(`  d: ${bits_d} bits`);
console.log(`  d_p: ${bits_d_p} bits`);
console.log(`  d_q: ${bits_d_q} bits`);

// ============================================================================
// Implementation Details
// ============================================================================

console.log('\n=== Implementation Details ===\n');

/**
 * The CRT private key representation stores:
 * - p, q: The prime factors
 * - d_p = d mod (p-1): Private exponent mod (p-1)
 * - d_q = d mod (q-1): Private exponent mod (q-1)
 * - q_inv = q^(-1) mod p: CRT coefficient
 *
 * Some implementations also store:
 * - d: The full private exponent (for verification)
 * - p_inv = p^(-1) mod q: Alternative CRT coefficient
 */

console.log('CRT Private Key Structure:');
console.log('');
console.log('  {');
console.log(`    p: ${p},`);
console.log(`    q: ${q},`);
console.log(`    d_p: ${d_p},  // d mod (p-1)`);
console.log(`    d_q: ${d_q},  // d mod (q-1)`);
console.log(`    q_inv: ${q_inv}  // q^(-1) mod p`);
console.log('  }');

// ============================================================================
// CRT Decryption Function
// ============================================================================

console.log('\n=== CRT Decryption Function ===\n');

interface RSACRTKey {
  n: bigint;
  p: bigint;
  q: bigint;
  d_p: bigint;
  d_q: bigint;
  q_inv: bigint;
}

function rsaDecryptCRT(c: bigint, key: RSACRTKey): bigint {
  // Step 1: Reduce c modulo p and q
  const c_p = c % key.p;
  const c_q = c % key.q;

  // Step 2: Compute partial decryptions
  const m_p = power_mod(c_p, key.d_p, key.p);
  const m_q = power_mod(c_q, key.d_q, key.q);

  // Step 3: CRT combination (Garner's formula)
  let h = ((m_p - m_q) % key.p + key.p) % key.p;
  h = (h * key.q_inv) % key.p;

  // Step 4: Compute final result
  const m = m_q + key.q * h;

  return m;
}

const crtKey: RSACRTKey = { n, p, q, d_p, d_q, q_inv };

console.log('Testing CRT decryption function:');

const testMessages = [1n, 42n, 100n, 1000n, 3000n];

for (const msg of testMessages) {
  if (msg >= n) continue;

  const c_test = power_mod(msg, e, n);
  const m_recovered = rsaDecryptCRT(c_test, crtKey);

  console.log(`  m = ${msg}, c = ${c_test}, decrypted = ${m_recovered}, correct = ${m_recovered === msg}`);
}

// ============================================================================
// Security: Fault Attack on CRT
// ============================================================================

console.log('\n=== Security: Fault Attack on CRT ===\n');

/**
 * CRT decryption is vulnerable to fault attacks!
 *
 * If a computational fault occurs in one of the partial decryptions
 * (m_p or m_q), the attacker can factor n.
 *
 * Attack:
 * 1. Get correct signature s
 * 2. Induce fault to get faulty signature s'
 * 3. If fault was in m_q computation:
 *    - s ≡ m (mod p) and s ≡ m (mod q)
 *    - s' ≡ m (mod p) but s' ≢ m (mod q)
 * 4. Compute gcd(s - s', n) = p
 *
 * Defense: Always verify the result before outputting it.
 */

console.log('Fault attack demonstration:');
console.log('');

// Correct decryption
const secretMsg = 42n;
const faultCiphertext = power_mod(secretMsg, e, n);
const correctDecryption = rsaDecryptCRT(faultCiphertext, crtKey);

console.log(`Ciphertext: c = ${faultCiphertext}`);
console.log(`Correct decryption: m = ${correctDecryption}`);

// Simulate a fault in m_q computation
function rsaDecryptCRTWithFault(c: bigint, key: RSACRTKey, faultInQ: boolean): bigint {
  const c_p = c % key.p;
  const c_q = c % key.q;

  const m_p = power_mod(c_p, key.d_p, key.p);
  let m_q = power_mod(c_q, key.d_q, key.q);

  if (faultInQ) {
    // Simulate a fault: corrupt m_q
    m_q = (m_q + 1n) % key.q;
  }

  let h = ((m_p - m_q) % key.p + key.p) % key.p;
  h = (h * key.q_inv) % key.p;

  return m_q + key.q * h;
}

const faultyDecryption = rsaDecryptCRTWithFault(faultCiphertext, crtKey, true);
console.log(`Faulty decryption: m' = ${faultyDecryption}`);

// Attack: use the difference to factor n
const diff = correctDecryption > faultyDecryption ? correctDecryption - faultyDecryption : faultyDecryption - correctDecryption;

console.log(`\nFault attack:`);
console.log(`  Difference: |m - m'| = ${diff}`);

const factorFound = gcd(diff, n);
console.log(`  gcd(|m - m'|, n) = ${factorFound}`);

if (factorFound > 1n && factorFound < n) {
  console.log(`  Factor of n found: ${factorFound}`);
  console.log(`  Other factor: ${n / factorFound}`);
  console.log(`  RSA key is completely broken!`);
}

console.log('\nDefense: Always verify result');
console.log('  After CRT decryption, compute m^e mod n');
console.log('  Verify it equals the original ciphertext');
console.log('  If not, reject or retry');

// ============================================================================
// Secure CRT Decryption
// ============================================================================

console.log('\n=== Secure CRT Decryption ===\n');

function rsaDecryptCRTSecure(c: bigint, key: RSACRTKey, e_pub: bigint): bigint | null {
  // Standard CRT decryption
  const c_p = c % key.p;
  const c_q = c % key.q;

  const m_p = power_mod(c_p, key.d_p, key.p);
  const m_q = power_mod(c_q, key.d_q, key.q);

  let h = ((m_p - m_q) % key.p + key.p) % key.p;
  h = (h * key.q_inv) % key.p;

  const m = m_q + key.q * h;

  // SECURITY: Verify result
  const verification = power_mod(m, e_pub, key.n);
  if (verification !== c) {
    console.log('  Fault detected! Verification failed.');
    return null;
  }

  return m;
}

console.log('Secure CRT with verification:');

// Test with correct computation
const secureResult = rsaDecryptCRTSecure(faultCiphertext, crtKey, e);
console.log(`  Correct computation: ${secureResult}`);

// Note: In practice, the faulty computation would be detected
console.log('  Verification ensures faults are detected.');

// ============================================================================
// Blinding for Side-Channel Protection
// ============================================================================

console.log('\n=== CRT with Blinding ===\n');

/**
 * Blinding protects against timing and power analysis attacks.
 *
 * 1. Generate random r, compute r^e mod n
 * 2. Blind ciphertext: c' = c * r^e mod n
 * 3. Decrypt with CRT: m' = (c')^d mod n = m * r
 * 4. Unblind: m = m' * r^(-1) mod n
 *
 * The decryption now operates on a randomized value,
 * preventing side-channel leakage.
 */

function rsaDecryptCRTBlinded(c: bigint, key: RSACRTKey, e_pub: bigint): bigint {
  // Generate random blinding factor coprime to n
  let r: bigint;
  do {
    r = BigInt(Math.floor(Math.random() * (Number(key.n) - 2))) + 2n;
  } while (gcd(r, key.n) !== 1n);
  const r_e = power_mod(r, e_pub, key.n);
  const r_inv = inverse_mod(r, key.n);

  // Blind the ciphertext
  const c_blinded = (c * r_e) % key.n;

  // Decrypt (with CRT)
  const m_blinded = rsaDecryptCRT(c_blinded, key);

  // Unblind
  const m = (m_blinded * r_inv) % key.n;

  // Verify
  const verification = power_mod(m, e_pub, key.n);
  if (verification !== c) {
    throw new Error('Decryption verification failed');
  }

  return m;
}

console.log('Blinded CRT decryption:');

for (let i = 0; i < 3; i++) {
  const blindedResult = rsaDecryptCRTBlinded(faultCiphertext, crtKey, e);
  console.log(`  Attempt ${i + 1}: ${blindedResult}`);
}

console.log('  All results match despite different blinding factors.');

// ============================================================================
// Summary
// ============================================================================

console.log('\n=== Summary ===\n');

console.log('CRT OPTIMIZATION KEY POINTS:');
console.log('');
console.log('  1. Performance:');
console.log('     - ~4x speedup for decryption');
console.log('     - Essential for practical RSA');
console.log('');
console.log('  2. CRT key components:');
console.log('     - d_p = d mod (p-1)');
console.log('     - d_q = d mod (q-1)');
console.log('     - q_inv = q^(-1) mod p');
console.log('');
console.log('  3. Security considerations:');
console.log('     - Vulnerable to fault attacks');
console.log('     - ALWAYS verify: m^e mod n == c');
console.log('     - Use blinding for side-channel protection');
console.log('');
console.log('  4. Implementation checklist:');
console.log('     [ ] Store CRT components securely');
console.log('     [ ] Verify every decryption result');
console.log('     [ ] Use blinding in sensitive environments');
console.log('     [ ] Clear sensitive values from memory');

console.log('\n' + '='.repeat(60));
console.log('RSA CRT Optimization Tutorial Complete');
console.log('='.repeat(60));
