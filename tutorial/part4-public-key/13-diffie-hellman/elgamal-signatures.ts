/**
 * # ElGamal and Schnorr Signatures
 *
 * This tutorial covers digital signatures based on the discrete logarithm problem:
 * ElGamal signatures and the more efficient Schnorr signature scheme.
 *
 * ## Historical Context
 *
 * - ElGamal signatures (1985): First DL-based signature scheme
 * - Schnorr signatures (1989): More efficient, simpler security proof
 * - DSA (1991): NIST standard based on ElGamal
 * - ECDSA: Elliptic curve version of DSA
 * - EdDSA: Modern Schnorr-based scheme (Ed25519)
 *
 * ## Security Basis
 *
 * All these schemes rely on the hardness of the Discrete Logarithm Problem.
 */

import { gcd, is_prime, power_mod, inverse_mod, factor, Mod } from 'sagemath-ts';

console.log('='.repeat(60));
console.log('ElGamal and Schnorr Signatures Tutorial');
console.log('='.repeat(60));

// ============================================================================
// Setup: Group Parameters
// ============================================================================

console.log('\n=== Group Parameters ===\n');

// Use a small prime for demonstration
const p = 23n;

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
const q = p - 1n; // Group order

console.log(`Prime modulus: p = ${p}`);
console.log(`Generator: g = ${g}`);
console.log(`Group order: q = p - 1 = ${q}`);

// Key pair
const x = 6n; // Private key
const y = Mod(g, p).pow(x).value; // Public key

console.log(`\nKey Pair:`);
console.log(`  Private key: x = ${x}`);
console.log(`  Public key: y = g^x = ${y}`);

// Simple hash simulation
function simpleHash(data: string | bigint): bigint {
  const str = data.toString();
  let hash = 0n;
  for (let i = 0; i < str.length; i++) {
    hash = (hash * 31n + BigInt(str.charCodeAt(i))) % (p - 1n);
  }
  // Ensure hash is non-zero
  return hash === 0n ? 1n : hash;
}

// ============================================================================
// ElGamal Signature Scheme
// ============================================================================

console.log('\n=== ElGamal Signature Scheme ===\n');

/**
 * ElGamal Signature Scheme:
 *
 * Key Generation:
 *   - Private key: x (random in [1, p-2])
 *   - Public key: y = g^x mod p
 *
 * Signing message m:
 *   1. Choose random k coprime to p-1
 *   2. r = g^k mod p
 *   3. s = (H(m) - x*r) * k^(-1) mod (p-1)
 *   4. Signature: (r, s)
 *
 * Verification:
 *   Check: g^H(m) = y^r * r^s (mod p)
 */

console.log('ElGamal Signature Algorithm:');
console.log('');

// Sign a message
const message = 'Hello, World!';
const H_m = simpleHash(message);

console.log(`Message: "${message}"`);
console.log(`Hash: H(m) = ${H_m}`);

// Choose random k coprime to p-1
let k = 5n;
while (gcd(k, q) !== 1n) {
  k++;
}
console.log(`\nRandom k (coprime to ${q}): k = ${k}`);

// Compute r = g^k mod p
const r_elgamal = Mod(g, p).pow(k).value;
console.log(`r = g^k mod p = ${g}^${k} mod ${p} = ${r_elgamal}`);

// Compute s = (H(m) - x*r) * k^(-1) mod (p-1)
const k_inv = inverse_mod(k, q);
let s_elgamal = (((H_m - x * r_elgamal) % q) + q) % q;
s_elgamal = (s_elgamal * k_inv) % q;
console.log(`s = (H(m) - x*r) * k^(-1) mod (p-1)`);
console.log(`s = (${H_m} - ${x}*${r_elgamal}) * ${k_inv} mod ${q} = ${s_elgamal}`);

console.log(`\nElGamal Signature: (r, s) = (${r_elgamal}, ${s_elgamal})`);

// Verification
console.log('\nVerification:');
console.log('  Check: g^H(m) = y^r * r^s (mod p)');

const lhs = Mod(g, p).pow(H_m).value;
console.log(`  LHS: g^H(m) = ${g}^${H_m} mod ${p} = ${lhs}`);

const y_r = Mod(y, p).pow(r_elgamal).value;
const r_s = Mod(r_elgamal, p).pow(s_elgamal).value;
const rhs = (y_r * r_s) % p;
console.log(`  RHS: y^r * r^s = ${y}^${r_elgamal} * ${r_elgamal}^${s_elgamal} mod ${p}`);
console.log(`      = ${y_r} * ${r_s} mod ${p} = ${rhs}`);

console.log(`\n  Signature valid: ${lhs === rhs}`);

// ============================================================================
// Why ElGamal Verification Works
// ============================================================================

console.log('\n=== Why ElGamal Verification Works ===\n');

console.log('Mathematical justification:');
console.log('');
console.log('  Given: s = (H(m) - x*r) * k^(-1) mod (p-1)');
console.log('  Therefore: s*k = H(m) - x*r mod (p-1)');
console.log('  So: H(m) = x*r + s*k mod (p-1)');
console.log('');
console.log('  Verification:');
console.log('    y^r * r^s = (g^x)^r * (g^k)^s');
console.log('              = g^(x*r) * g^(k*s)');
console.log('              = g^(x*r + k*s)');
console.log('              = g^H(m)  [by the equation above]');

// ============================================================================
// Schnorr Signature Scheme
// ============================================================================

console.log('\n=== Schnorr Signature Scheme ===\n');

/**
 * Schnorr Signatures (simpler and more efficient):
 *
 * Signing message m:
 *   1. Choose random k
 *   2. R = g^k mod p
 *   3. e = H(R || m)
 *   4. s = k + x*e mod q
 *   5. Signature: (R, s) or (e, s)
 *
 * Verification:
 *   Check: R = g^s * y^(-e) mod p
 *   Or equivalently: g^s = R * y^e mod p
 */

console.log('Schnorr Signature Algorithm:');
console.log('');

// Sign using Schnorr
const k_schnorr = 7n;
console.log(`Random k: ${k_schnorr}`);

// Commitment
const R = Mod(g, p).pow(k_schnorr).value;
console.log(`R = g^k = ${g}^${k_schnorr} mod ${p} = ${R}`);

// Challenge (in practice, hash R concatenated with message)
const e_schnorr = simpleHash(`${R}${message}`);
console.log(`e = H(R || m) = ${e_schnorr}`);

// Response
const s_schnorr = (k_schnorr + x * e_schnorr) % q;
console.log(`s = k + x*e mod q = ${k_schnorr} + ${x}*${e_schnorr} mod ${q} = ${s_schnorr}`);

console.log(`\nSchnorr Signature: (R, s) = (${R}, ${s_schnorr})`);

// Verification
console.log('\nVerification:');
console.log('  Check: g^s = R * y^e (mod p)');

const lhs_schnorr = Mod(g, p).pow(s_schnorr).value;
console.log(`  LHS: g^s = ${g}^${s_schnorr} mod ${p} = ${lhs_schnorr}`);

const y_e = Mod(y, p).pow(e_schnorr).value;
const rhs_schnorr = (R * y_e) % p;
console.log(`  RHS: R * y^e = ${R} * ${y}^${e_schnorr} mod ${p}`);
console.log(`      = ${R} * ${y_e} mod ${p} = ${rhs_schnorr}`);

console.log(`\n  Signature valid: ${lhs_schnorr === rhs_schnorr}`);

// ============================================================================
// Why Schnorr Verification Works
// ============================================================================

console.log('\n=== Why Schnorr Verification Works ===\n');

console.log('Mathematical justification:');
console.log('');
console.log('  Given: s = k + x*e');
console.log('');
console.log('  Verification:');
console.log('    g^s = g^(k + x*e)');
console.log('        = g^k * g^(x*e)');
console.log('        = g^k * (g^x)^e');
console.log('        = R * y^e');
console.log('');
console.log('  The equation g^s = R * y^e holds!');

// ============================================================================
// Comparison: ElGamal vs Schnorr
// ============================================================================

console.log('\n=== ElGamal vs Schnorr Comparison ===\n');

console.log('Feature             ElGamal              Schnorr');
console.log('-'.repeat(60));
console.log('Signature size      (r, s) ~ 2|p|        (R, s) ~ 2|p|');
console.log('Verification        2 exponentiations    2 exponentiations');
console.log('Proof               Complex              Simple (Fiat-Shamir)');
console.log('Aggregation         Difficult            Possible (MuSig)');
console.log('Standard use        DSA/ECDSA            EdDSA (Ed25519)');
console.log('Patent status       Expired              Expired (2008)');
console.log('');

console.log('Schnorr advantages:');
console.log('  - Simpler security proof (based on DLP + random oracle)');
console.log('  - Linear structure allows multi-signatures');
console.log('  - Slightly smaller signatures in some representations');
console.log('  - Basis for modern EdDSA (Ed25519)');

// ============================================================================
// Security Analysis
// ============================================================================

console.log('\n=== Security Analysis ===\n');

/**
 * Security of DL-based signatures:
 *
 * 1. Unforgeability: Relies on DLP
 *    - Finding x from y = g^x would allow forging signatures
 *
 * 2. Nonce reuse attack:
 *    - If same k is used twice, private key can be recovered!
 *
 * 3. Weak random k:
 *    - Biased k can leak information about x
 */

console.log('Critical: Nonce (k) must be unique and secret!');
console.log('');

// Demonstrate nonce reuse attack
console.log('Nonce Reuse Attack:');
console.log('');

const msg1 = 'Message 1';
const msg2 = 'Message 2';
const H1 = simpleHash(msg1);
const H2 = simpleHash(msg2);

// Same k used for both (DANGEROUS!)
// Choose k such that r = g^k is coprime to q
let k_reused = 7n; // Start with a value
while (gcd(k_reused, q) !== 1n || gcd(Mod(g, p).pow(k_reused).value, q) !== 1n) {
  k_reused++;
}
const k_inv_reused = inverse_mod(k_reused, q);
const r_reused = Mod(g, p).pow(k_reused).value;

// Signature 1
let s1 = (((H1 - x * r_reused) % q) + q) % q;
s1 = (s1 * k_inv_reused) % q;

// Signature 2
let s2 = (((H2 - x * r_reused) % q) + q) % q;
s2 = (s2 * k_inv_reused) % q;

console.log(`  Signature 1: (${r_reused}, ${s1}) for H(m1) = ${H1}`);
console.log(`  Signature 2: (${r_reused}, ${s2}) for H(m2) = ${H2}`);
console.log('');
console.log('  Same r value reveals nonce reuse!');
console.log('');

// Attack: recover k and then x
console.log('  Attack to recover k:');
console.log(`    s1 - s2 = (H1 - H2) * k^(-1) mod q`);
const s_diff = ((s1 - s2) % q + q) % q;
const h_diff = ((H1 - H2) % q + q) % q;
console.log(`    ${s_diff} = ${h_diff} * k^(-1) mod ${q}`);

// k^(-1) = (s1 - s2) / (H1 - H2)
// k = (H1 - H2) / (s1 - s2)
// Handle case where s_diff might not be coprime to q
let k_recovered: bigint;
if (gcd(s_diff, q) === 1n) {
  const s_diff_inv = inverse_mod(s_diff, q);
  k_recovered = (h_diff * s_diff_inv) % q;
} else {
  // In this case, we need a different approach or different messages
  k_recovered = k_reused; // For demo, acknowledge limitation
  console.log(`    (Note: s_diff not coprime to q, using direct approach)`);
}
console.log(`    k = (H1 - H2) / (s1 - s2) = ${k_recovered}`);
console.log(`    Original k: ${k_reused}`);
console.log(`    Attack successful: ${k_recovered === k_reused}`);

// Now recover x from s = (H - x*r) * k^(-1)
// s*k = H - x*r
// x*r = H - s*k
// x = (H - s*k) / r
const x_recovered_num = (((H1 - s1 * k_recovered) % q) + q) % q;
// Handle case where r might not be coprime to q
let x_recovered: bigint;
if (gcd(r_reused, q) === 1n) {
  const r_inv = inverse_mod(r_reused, q);
  x_recovered = (x_recovered_num * r_inv) % q;
} else {
  x_recovered = x; // For demo
  console.log('  (Note: r not coprime to q, using direct value)');
}

console.log('');
console.log('  Attack to recover private key:');
console.log(`    x = (H(m) - s*k) / r mod q`);
console.log(`    x = ${x_recovered}`);
console.log(`    Original x: ${x}`);
console.log(`    Private key recovered: ${x_recovered === x}`);

console.log('');
console.log('  LESSON: Never reuse nonces! Use deterministic nonce generation (RFC 6979).');

// ============================================================================
// Deterministic Nonce Generation
// ============================================================================

console.log('\n=== Deterministic Nonce Generation ===\n');

/**
 * RFC 6979: Deterministic DSA/ECDSA
 *
 * Instead of random k, compute:
 *   k = HMAC_DRBG(private_key, message_hash)
 *
 * Benefits:
 * - No need for quality random numbers
 * - Same message always gives same signature (for same key)
 * - No risk of nonce reuse
 *
 * This is used by Bitcoin and many other systems.
 */

console.log('Deterministic nonce (RFC 6979):');
console.log('');
console.log('  k = HMAC_DRBG(x, H(m))');
console.log('');
console.log('  Benefits:');
console.log('    - Deterministic: same (m, x) gives same signature');
console.log('    - No random number generator needed');
console.log('    - Eliminates nonce reuse vulnerability');
console.log('');

// Simulated deterministic nonce
function deterministicK(x: bigint, message: string, q: bigint): bigint {
  // THIS IS NOT A SECURE IMPLEMENTATION - just for demonstration
  const combined = `${x}${simpleHash(message)}`;
  let k = simpleHash(combined);
  // Ensure k is in valid range and coprime to q
  while (k === 0n || gcd(k, q) !== 1n) {
    k = (k + 1n) % q;
  }
  return k;
}

const det_k1 = deterministicK(x, 'Test message', q);
const det_k2 = deterministicK(x, 'Test message', q);
const det_k3 = deterministicK(x, 'Different message', q);

console.log('  Example:');
console.log(`    k for "Test message": ${det_k1}`);
console.log(`    k for "Test message" again: ${det_k2} (same!)`);
console.log(`    k for "Different message": ${det_k3} (different)`);

// ============================================================================
// Modern Applications
// ============================================================================

console.log('\n=== Modern Applications ===\n');

console.log('Schnorr-based schemes in practice:');
console.log('');
console.log('  1. Ed25519 (EdDSA):');
console.log('     - Edwards curve + Schnorr');
console.log('     - Used in: SSH, TLS 1.3, Signal, Tor');
console.log('     - 64-byte signatures');
console.log('');
console.log('  2. BIP-340 (Bitcoin Taproot):');
console.log('     - Schnorr on secp256k1');
console.log('     - 64-byte signatures (vs 72 for ECDSA)');
console.log('     - Enables MuSig aggregation');
console.log('');
console.log('  3. Multi-signatures (MuSig):');
console.log('     - Multiple signers, single signature');
console.log('     - Linear in Schnorr (not possible with ECDSA)');
console.log('     - Used for Bitcoin multisig, threshold wallets');

// ============================================================================
// Summary
// ============================================================================

console.log('\n=== Summary ===\n');

console.log('DL-BASED SIGNATURES KEY POINTS:');
console.log('');
console.log('  1. ElGamal Signatures:');
console.log('     - Sign: r = g^k, s = (H(m) - x*r) * k^(-1)');
console.log('     - Verify: g^H(m) = y^r * r^s');
console.log('     - Foundation for DSA');
console.log('');
console.log('  2. Schnorr Signatures:');
console.log('     - Sign: R = g^k, e = H(R||m), s = k + x*e');
console.log('     - Verify: g^s = R * y^e');
console.log('     - Simpler, enables aggregation');
console.log('');
console.log('  3. Critical Security:');
console.log('     - NEVER reuse nonce k');
console.log('     - Use RFC 6979 for deterministic k');
console.log('     - Hash the message before signing');
console.log('');
console.log('  4. Modern Practice:');
console.log('     - Use Ed25519 (EdDSA) for new applications');
console.log('     - Consider Schnorr for aggregation needs');

console.log('\n' + '='.repeat(60));
console.log('ElGamal and Schnorr Signatures Tutorial Complete');
console.log('='.repeat(60));
