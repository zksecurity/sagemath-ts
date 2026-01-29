/**
 * # CRT Optimization in RSA
 *
 * The Chinese Remainder Theorem provides a significant speedup
 * for RSA decryption and signing operations.
 *
 * ## The Problem
 *
 * RSA decryption computes: m = c^d mod n
 * where n = p*q and d can be very large (typically 2048+ bits).
 *
 * ## The CRT Speedup
 *
 * Instead of one large exponentiation mod n:
 * 1. Compute m_p = c^(d mod (p-1)) mod p
 * 2. Compute m_q = c^(d mod (q-1)) mod q
 * 3. Combine using CRT: m = CRT(m_p, m_q, p, q)
 *
 * This is approximately 4x faster because:
 * - Two exponentiations with half-size moduli
 * - Exponentiation is O(k^3) for k-bit numbers
 * - Two operations at k/2 bits: 2 * (k/2)^3 = k^3/4
 *
 * @module tutorial/part2/crt-rsa
 */

import {
  crt,
  power_mod,
  inverse_mod,
  gcd,
  euler_phi,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("CRT OPTIMIZATION IN RSA");
console.log("=".repeat(70));

// ============================================================================
// Setup: RSA Key Generation
// ============================================================================
console.log("\n=== RSA Key Setup ===\n");

// Using small primes for demonstration (real RSA uses ~1024-bit primes)
const p = 61n;
const q = 53n;
const n = p * q;
const phiN = (p - 1n) * (q - 1n);
const e = 17n;
const d = inverse_mod(e, phiN);

console.log("RSA Parameters:");
console.log(`  p = ${p}`);
console.log(`  q = ${q}`);
console.log(`  n = p * q = ${n}`);
console.log(`  phi(n) = ${phiN}`);
console.log(`  e = ${e} (public exponent)`);
console.log(`  d = ${d} (private exponent)`);

// ============================================================================
// Standard RSA Decryption
// ============================================================================
console.log("\n=== Standard RSA Decryption ===\n");

const message = 42n;
const ciphertext = power_mod(message, e, n);

console.log(`Original message: m = ${message}`);
console.log(`Ciphertext: c = m^e mod n = ${message}^${e} mod ${n} = ${ciphertext}`);

// Standard decryption: m = c^d mod n
const decryptedStandard = power_mod(ciphertext, d, n);
console.log(`\nStandard decryption: m = c^d mod n`);
console.log(`  m = ${ciphertext}^${d} mod ${n} = ${decryptedStandard}`);

// ============================================================================
// CRT-Optimized RSA Decryption
// ============================================================================
console.log("\n=== CRT-Optimized RSA Decryption ===\n");

// Precompute CRT parameters (done once during key generation)
const dp = d % (p - 1n);  // d mod (p-1)
const dq = d % (q - 1n);  // d mod (q-1)
const qInv = inverse_mod(q, p);  // q^(-1) mod p

console.log("Precomputed CRT parameters:");
console.log(`  d_p = d mod (p-1) = ${d} mod ${p - 1n} = ${dp}`);
console.log(`  d_q = d mod (q-1) = ${d} mod ${q - 1n} = ${dq}`);
console.log(`  q^(-1) mod p = ${qInv}`);

// CRT decryption
console.log("\nCRT Decryption Steps:");

// Step 1: Compute m_p = c^d_p mod p
const mp = power_mod(ciphertext, dp, p);
console.log(`  1. m_p = c^d_p mod p = ${ciphertext}^${dp} mod ${p} = ${mp}`);

// Step 2: Compute m_q = c^d_q mod q
const mq = power_mod(ciphertext, dq, q);
console.log(`  2. m_q = c^d_q mod q = ${ciphertext}^${dq} mod ${q} = ${mq}`);

// Step 3: Combine using CRT (Garner's formula)
// m = m_q + q * ((m_p - m_q) * q^(-1) mod p)
const h = (((mp - mq) * qInv) % p + p) % p;
const decryptedCRT = mq + q * h;

console.log(`  3. h = (m_p - m_q) * q^(-1) mod p`);
console.log(`       = (${mp} - ${mq}) * ${qInv} mod ${p}`);
console.log(`       = ${h}`);
console.log(`  4. m = m_q + q * h = ${mq} + ${q} * ${h} = ${decryptedCRT}`);

console.log(`\nResults:`);
console.log(`  Standard decryption: ${decryptedStandard}`);
console.log(`  CRT decryption:      ${decryptedCRT}`);
console.log(`  Match: ${decryptedStandard === decryptedCRT}`);

// ============================================================================
// Why CRT is Faster
// ============================================================================
console.log("\n=== Why CRT is Faster ===\n");

console.log("Operation Comparison:");
console.log("-".repeat(60));

console.log("\nStandard: c^d mod n");
console.log(`  - Exponent size: ${d.toString(2).length} bits (d)`);
console.log(`  - Modulus size: ${n.toString(2).length} bits (n)`);
console.log(`  - Work: O(|d| * |n|^2) = O(k * k^2) = O(k^3)`);

console.log("\nCRT: c^d_p mod p and c^d_q mod q");
console.log(`  - Exponent sizes: ${dp.toString(2).length} bits (d_p), ${dq.toString(2).length} bits (d_q)`);
console.log(`  - Modulus sizes: ${p.toString(2).length} bits (p), ${q.toString(2).length} bits (q)`);
console.log(`  - Work: 2 * O((k/2) * (k/2)^2) = 2 * O(k^3/8) = O(k^3/4)`);

console.log("\nSpeedup: Approximately 4x");
console.log("(In practice, 3-4x speedup is typical)");

// ============================================================================
// Mathematical Justification
// ============================================================================
console.log("\n=== Mathematical Justification ===\n");

console.log("Why does c^d_p mod p equal c^d mod n (mod p)?");
console.log("-".repeat(50));

console.log("\nBy Fermat's Little Theorem:");
console.log(`  c^(p-1) = 1 (mod p)  [since p is prime and gcd(c, p) = 1]`);

console.log("\nSince d = d_p + k*(p-1) for some integer k:");
console.log(`  c^d = c^(d_p + k*(p-1))`);
console.log(`      = c^d_p * (c^(p-1))^k`);
console.log(`      = c^d_p * 1^k`);
console.log(`      = c^d_p (mod p)`);

console.log("\nSimilarly:");
console.log(`  c^d = c^d_q (mod q)`);

console.log("\nBy CRT, the unique solution mod n is:");
console.log(`  m = CRT(c^d_p mod p, c^d_q mod q)`);

// Verify this reasoning
console.log("\n Verification:");
console.log(`  c^d mod p = ${power_mod(ciphertext, d, p)}`);
console.log(`  c^d_p mod p = ${power_mod(ciphertext, dp, p)}`);
console.log(`  Equal: ${power_mod(ciphertext, d, p) === power_mod(ciphertext, dp, p)}`);

// ============================================================================
// Garner's Algorithm
// ============================================================================
console.log("\n=== Garner's Algorithm for CRT Reconstruction ===\n");

console.log("The standard CRT formula requires computing the full product M.");
console.log("Garner's algorithm computes the result more efficiently.\n");

console.log("For two moduli p and q:");
console.log("  m = m_q + q * ((m_p - m_q) * q^(-1) mod p)\n");

console.log("This avoids computing n = p*q in intermediate steps,");
console.log("which is important when p and q are large.\n");

function garnerCRT(mp: bigint, mq: bigint, p: bigint, q: bigint, qInvModP: bigint): bigint {
  const h = (((mp - mq) * qInvModP) % p + p) % p;
  return mq + q * h;
}

// Verify Garner's formula matches standard CRT
const garnerResult = garnerCRT(mp, mq, p, q, qInv);
const standardCRTResult = crt(mp, mq, p, q);

console.log(`Garner's algorithm: ${garnerResult}`);
console.log(`Standard CRT:       ${standardCRTResult}`);
console.log(`Match: ${garnerResult === standardCRTResult}`);

// ============================================================================
// CRT for RSA Signing
// ============================================================================
console.log("\n=== CRT for RSA Signing ===\n");

console.log("RSA signing: s = m^d mod n (same operation as decryption)");
console.log("CRT provides the same 4x speedup for signing.\n");

const messageToSign = 123n;
console.log(`Message hash: h = ${messageToSign}`);

// Standard signing
const signatureStandard = power_mod(messageToSign, d, n);

// CRT signing
const sp = power_mod(messageToSign, dp, p);
const sq = power_mod(messageToSign, dq, q);
const signatureCRT = garnerCRT(sp, sq, p, q, qInv);

console.log(`Standard signature: ${signatureStandard}`);
console.log(`CRT signature:      ${signatureCRT}`);
console.log(`Match: ${signatureStandard === signatureCRT}`);

// Verify signature
const verified = power_mod(signatureCRT, e, n);
console.log(`\nVerification: s^e mod n = ${verified}`);
console.log(`Original message: ${messageToSign}`);
console.log(`Valid signature: ${verified === messageToSign}`);

// ============================================================================
// Fault Attack Consideration
// ============================================================================
console.log("\n=== Security Consideration: Fault Attacks ===\n");

console.log("CRT-RSA is vulnerable to fault attacks!");
console.log("If a computation error occurs in one branch:\n");

// Simulate a fault in the m_q computation
const mpCorrect = power_mod(ciphertext, dp, p);
const mqFaulty = (power_mod(ciphertext, dq, q) + 1n) % q; // Simulated fault

console.log(`Correct m_p: ${mpCorrect}`);
console.log(`Faulty m_q:  ${mqFaulty} (should be ${power_mod(ciphertext, dq, q)})`);

const faultySignature = garnerCRT(mpCorrect, mqFaulty, p, q, qInv);
console.log(`Faulty decryption: ${faultySignature}`);

// The attacker can use this to factor n!
// gcd(m - m_faulty, n) often reveals p or q
const diff = decryptedCRT - faultySignature;
const factor = gcd(diff < 0n ? -diff : diff, n);

console.log(`\nAttack: gcd(m_correct - m_faulty, n)`);
console.log(`      = gcd(${decryptedCRT} - ${faultySignature}, ${n})`);
console.log(`      = gcd(${diff}, ${n})`);
console.log(`      = ${factor}`);

if (factor !== 1n && factor !== n) {
  console.log(`\nFactor found: ${factor} (this is ${factor === p ? "p" : "q"}!)`);
  console.log("This allows complete key recovery!");
}

console.log("\nCountermeasures:");
console.log("  1. Verify signature after generation: check s^e = m");
console.log("  2. Blinding: compute (r^e * m)^d then divide by r");
console.log("  3. Hardware protections against fault injection");

// ============================================================================
// Multi-Prime RSA
// ============================================================================
console.log("\n=== Multi-Prime RSA ===\n");

console.log("Using more than 2 primes provides additional speedup.");
console.log("n = p1 * p2 * ... * pk\n");

// Three-prime RSA example
const p1 = 29n;
const p2 = 31n;
const p3 = 37n;
const n3 = p1 * p2 * p3;
const phi3 = (p1 - 1n) * (p2 - 1n) * (p3 - 1n);
const e3 = 11n;
const d3 = inverse_mod(e3, phi3);

console.log(`Three-prime RSA: n = ${p1} * ${p2} * ${p3} = ${n3}`);
console.log(`phi(n) = ${phi3}, d = ${d3}\n`);

const m3 = 100n;
const c3 = power_mod(m3, e3, n3);

// CRT decryption
const d1 = d3 % (p1 - 1n);
const d2 = d3 % (p2 - 1n);
const d3mod = d3 % (p3 - 1n);

const m1 = power_mod(c3, d1, p1);
const m2 = power_mod(c3, d2, p2);
const m3r = power_mod(c3, d3mod, p3);

console.log("CRT decryption:");
console.log(`  m_1 = c^d_1 mod p1 = ${m1}`);
console.log(`  m_2 = c^d_2 mod p2 = ${m2}`);
console.log(`  m_3 = c^d_3 mod p3 = ${m3r}`);

// Reconstruct
const decrypted3 = ((m1: bigint, m2: bigint, m3: bigint, p1: bigint, p2: bigint, p3: bigint) => {
  // First combine m1 and m2
  const partial = crt(m1, m2, p1, p2);
  // Then combine with m3
  return crt(partial, m3, p1 * p2, p3);
})(m1, m2, m3r, p1, p2, p3);

console.log(`  Reconstructed: ${decrypted3}`);
console.log(`  Original: ${m3}`);
console.log(`  Correct: ${decrypted3 === m3}`);

console.log(`
Speedup for k primes vs 2 primes:
  2 primes: 2 * (n/2)^3 / n^3 = 1/4 (4x speedup)
  3 primes: 3 * (n/3)^3 / n^3 = 1/9 (9x speedup)
  k primes: k * (n/k)^3 / n^3 = 1/k^2 (k^2 speedup)

But more primes = more fault attack surface!
`);

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
CRT Optimization in RSA:

Standard decryption: m = c^d mod n

CRT decryption:
  1. Precompute: d_p = d mod (p-1), d_q = d mod (q-1), q^(-1) mod p
  2. Compute: m_p = c^d_p mod p, m_q = c^d_q mod q
  3. Combine: m = m_q + q * ((m_p - m_q) * q^(-1) mod p)

Speedup: ~4x for 2-prime RSA

Why it works:
  - Fermat's theorem: c^(p-1) = 1 mod p
  - So c^d = c^(d mod (p-1)) mod p
  - CRT reconstructs the unique answer mod n

Security Considerations:
  - Fault attacks can factor n
  - Countermeasures: verification, blinding

Multi-prime RSA:
  - k primes gives k^2 speedup
  - Trade-off with fault attack surface
`);
