/**
 * # Factorization: Overview and Cryptographic Implications
 *
 * Integer factorization is one of the fundamental hard problems in
 * cryptography. The security of RSA and many other cryptosystems
 * depends on the difficulty of factoring large integers.
 *
 * ## The Factoring Hierarchy
 *
 * From simplest to most sophisticated:
 * 1. Trial division: O(sqrt(n))
 * 2. Pollard's rho: O(n^(1/4))
 * 3. Pollard's p-1: O(B log B) when p-1 is B-smooth
 * 4. Elliptic curve method (ECM): O(exp(sqrt(2 log p log log p)))
 * 5. Quadratic sieve: O(exp(sqrt(log n log log n)))
 * 6. Number field sieve (NFS): O(exp(c * (log n)^(1/3) * (log log n)^(2/3)))
 *
 * ## Records and Estimates
 *
 * - RSA-250 (829 bits, 250 digits): Factored in 2020, ~2700 CPU years
 * - RSA-2048: Estimated >10^12 CPU years with current algorithms
 *
 * This is why RSA-2048 is considered secure (for now).
 */

import {
  factor,
  formatFactorization,
  is_prime,
  isqrt,
  gcd,
  power_mod,
  prime_range,
} from 'sagemath-ts';

console.log("=== Factorization Overview ===\n");

// Example 1: Algorithm complexity comparison
console.log("Example 1: Algorithm Complexity Comparison");
console.log("-".repeat(60));

console.log(`
  Comparing factoring algorithms for n = p * q where p ~ q ~ sqrt(n):

  Algorithm          | Complexity                    | 20 digits | 50 digits
  -------------------|-------------------------------|-----------|----------
  Trial division     | O(n^(1/2))                    | 10^10     | 10^25
  Pollard's rho      | O(n^(1/4))                    | 10^5      | 10^12.5
  Pollard's p-1      | O(B log B) [if p-1 smooth]    | varies    | varies
  ECM                | O(exp(sqrt(2 log p log log p)))| 10^4      | 10^8
  Quadratic sieve    | O(exp(sqrt(log n log log n))) | 10^3      | 10^6
  Number field sieve | O(exp((log n)^(1/3)...))      | 10^2      | 10^5

  Note: These are rough estimates. Actual times depend on many factors.
`);

// Example 2: The RSA Challenge numbers
console.log("Example 2: RSA Challenge History");
console.log("-".repeat(60));

const rsaChallenges = [
  { name: "RSA-100", bits: 330, digits: 100, year: 1991, factored: 1991 },
  { name: "RSA-129", bits: 426, digits: 129, year: 1991, factored: 1994 },
  { name: "RSA-155", bits: 512, digits: 155, year: 1991, factored: 1999 },
  { name: "RSA-200", bits: 663, digits: 200, year: 1991, factored: 2005 },
  { name: "RSA-250", bits: 829, digits: 250, year: 1991, factored: 2020 },
  { name: "RSA-270", bits: 896, digits: 270, year: 1991, factored: null },
  { name: "RSA-2048", bits: 2048, digits: 617, year: 1991, factored: null },
];

console.log("  Number   | Bits | Digits | Status");
console.log("  " + "-".repeat(50));
for (const r of rsaChallenges) {
  const status = r.factored ? `Factored ${r.factored}` : "Unfactored";
  console.log(`  ${r.name.padEnd(10)} | ${r.bits.toString().padStart(4)} | ${r.digits.toString().padStart(6)} | ${status}`);
}

console.log(`
  Key insight: Each 50 additional digits takes about 10x more effort.
  RSA-2048 (~617 digits) is estimated to require >10^12 CPU years.
`);

// Example 3: Special number forms
console.log("\nExample 3: Special Number Forms");
console.log("-".repeat(60));

console.log(`
  Some numbers have special structure that makes factoring easier
  or harder:

  Fermat numbers: F_n = 2^(2^n) + 1
  - F_0 through F_4 are prime: 3, 5, 17, 257, 65537
  - F_5 = 2^32 + 1 is composite (Euler found a factor)
  - No other Fermat primes are known

  Mersenne numbers: M_p = 2^p - 1
  - Many are prime: M_2, M_3, M_5, M_7, M_13, ...
  - Special factoring methods exist (Lucas-Lehmer)

  Cunningham project: Factors a^n +/- 1 for small a
  - Extensive tables maintained since 1925
`);

// Check small Fermat numbers
console.log("  Fermat numbers:");
for (let n = 0; n <= 5; n++) {
  const F_n = (1n << (1n << BigInt(n))) + 1n;
  const isPr = is_prime(F_n);
  const digits = F_n.toString().length;
  console.log(`    F_${n} = 2^(2^${n}) + 1 = ${F_n} (${digits} digits) - ${isPr ? "prime" : "composite"}`);
}

// Example 4: Why factoring is hard
console.log("\nExample 4: Why Factoring is Hard");
console.log("-".repeat(60));

console.log(`
  Asymmetry of multiplication vs factoring:

  Multiplication: O(n^2) or O(n log n) with FFT
  - Multiplying two 1000-digit numbers: microseconds

  Factoring: Best known is sub-exponential but super-polynomial
  - Factoring a 2000-digit number: effectively impossible

  This asymmetry is the foundation of RSA security.

  The Factoring Assumption:
  "For random n-bit primes p, q, given n = p*q, finding p or q
  requires time exponential in n."

  This is an assumption, not a proven fact!
  - No one has proven factoring is hard
  - No one has found a polynomial-time algorithm
  - Quantum computers could change everything (Shor's algorithm)
`);

// Example 5: How sagemath-ts factors numbers
console.log("\nExample 5: Factoring with sagemath-ts");
console.log("-".repeat(60));

const testNumbers = [
  360n,                    // Small, smooth
  1234567890123456789n,    // Medium
  2n ** 64n + 1n,          // Fermat-like
  1000000007n * 1000000009n, // Product of two primes
];

for (const n of testNumbers) {
  const start = performance.now();
  const f = factor(n);
  const time = performance.now() - start;
  console.log(`  factor(${n.toString().slice(0, 30)}${n.toString().length > 30 ? "..." : ""})`);
  console.log(`    = ${formatFactorization(f)}`);
  console.log(`    Time: ${time.toFixed(2)}ms\n`);
}

// Example 6: Factoring as an oracle
console.log("Example 6: Factoring as a Cryptographic Oracle");
console.log("-".repeat(60));

console.log(`
  If you can factor n = p * q, you can:

  1. Break RSA:
     - Compute phi(n) = (p-1)(q-1)
     - Find d = e^(-1) mod phi(n)
     - Decrypt any message

  2. Forge signatures:
     - With the private key d, sign any message
     - Impersonate the key owner

  3. Break related schemes:
     - Rabin cryptosystem
     - Blum-Blum-Shub PRNG
     - Some zero-knowledge proofs

  The security reduction:
  "Breaking RSA is at least as hard as factoring n"
  (Actually, we only know: "RSA <= factoring", not equality)
`);

// Example 7: Computing phi(n) from factors
console.log("\nExample 7: From Factors to Private Key");
console.log("-".repeat(60));

// Small RSA example
const p = 61n;
const q = 53n;
const n = p * q;
const e = 17n; // Public exponent

console.log(`  RSA key generation:`);
console.log(`    p = ${p}, q = ${q}`);
console.log(`    n = p * q = ${n}`);
console.log(`    e = ${e} (public exponent)`);

// Compute private key
const phi_n = (p - 1n) * (q - 1n);
console.log(`    phi(n) = (p-1)(q-1) = ${phi_n}`);

// Extended Euclidean algorithm for modular inverse
function modInverse(a: bigint, m: bigint): bigint {
  let [old_r, r] = [a, m];
  let [old_s, s] = [1n, 0n];

  while (r !== 0n) {
    const q = old_r / r;
    [old_r, r] = [r, old_r - q * r];
    [old_s, s] = [s, old_s - q * s];
  }

  return ((old_s % m) + m) % m;
}

const d = modInverse(e, phi_n);
console.log(`    d = e^(-1) mod phi(n) = ${d}`);

// Verify: e * d = 1 mod phi(n)
console.log(`    Verify: e * d mod phi(n) = ${(e * d) % phi_n}`);

// Encrypt and decrypt a message
const message = 42n;
const ciphertext = power_mod(message, e, n);
const decrypted = power_mod(ciphertext, d, n);

console.log(`\n  Encryption/Decryption:`);
console.log(`    Message:    ${message}`);
console.log(`    Encrypted:  ${ciphertext}`);
console.log(`    Decrypted:  ${decrypted}`);

// Example 8: Safe prime generation
console.log("\nExample 8: Generating Factoring-Resistant Primes");
console.log("-".repeat(60));

console.log(`
  Historical guidelines for RSA primes:

  1. Strong primes: p-1 has a large prime factor
     - Resists Pollard's p-1 attack

  2. p+1 should have a large prime factor
     - Resists Williams' p+1 attack

  3. (p-1)/2 should be prime (safe prime)
     - p = 2q + 1 where q is prime (Sophie Germain)
     - Useful for some protocols (DH)

  Modern view:
  - For RSA-2048+, random primes are secure enough
  - The probability of a random 1024-bit prime being
    weak is negligibly small
  - Strong prime generation is unnecessarily slow
`);

// Example 9: Quantum threat
console.log("Example 9: The Quantum Threat");
console.log("-".repeat(60));

console.log(`
  Shor's Algorithm (1994):
  - Factors n in O((log n)^3) using a quantum computer
  - Would break RSA, DSA, DH, ECDSA, etc.

  Current status:
  - Largest number factored by Shor: 21 = 3 * 7 (2012)
  - Need millions of qubits for RSA-2048
  - Currently: ~100 noisy qubits available

  Timeline estimates vary wildly:
  - Optimistic: 10-15 years
  - Pessimistic: 50+ years
  - Some say never at scale

  Post-quantum cryptography:
  - Lattice-based (NTRU, Kyber, Dilithium)
  - Code-based (McEliece)
  - Hash-based signatures (SPHINCS+)
  - Isogeny-based (SIKE - broken in 2022!)

  NIST PQC standardization: Kyber and Dilithium selected (2024)
`);

// Example 10: Summary of factoring methods
console.log("\nExample 10: Factoring Method Selection");
console.log("-".repeat(60));

console.log(`
  When to use each method:

  Trial Division:
  - Small numbers (< 10^12)
  - Finding small factors of large numbers
  - Preprocessing step

  Pollard's Rho:
  - Medium-sized factors (10-25 digits)
  - Unknown factor size
  - Simple implementation

  Pollard's p-1:
  - When p-1 might be smooth
  - Very fast when it works
  - Quick to try and abandon

  ECM (Elliptic Curve Method):
  - Factors up to 50+ digits
  - Time depends on factor size, not n
  - Highly parallelizable

  Quadratic Sieve:
  - Numbers up to ~110 digits
  - Simpler than NFS
  - Good for medium-sized numbers

  Number Field Sieve:
  - Largest factorizations (100+ digits)
  - Most complex algorithm
  - State of the art
`);

// Final summary
console.log("\n=== Summary ===");
console.log(`
  Integer Factorization:

  Core problem: Given n = p * q, find p and q.

  Complexity status:
  - Not known to be in P (no polynomial algorithm)
  - Not known to be NP-complete
  - In NP and coNP (witness exists for both outcomes)
  - Solvable in quantum polynomial time (Shor)

  Cryptographic implications:
  - RSA security relies on factoring hardness
  - Key sizes must account for algorithm advances
  - Quantum computers pose an existential threat

  Algorithm landscape:
  - General-purpose: NFS (best for large n)
  - Special-purpose: ECM (best for finding small factors)
  - Historical: QS, p-1, rho (still useful in practice)

  Current recommendations:
  - RSA: Use 2048+ bit moduli (3072 for longer term)
  - Consider post-quantum alternatives for new systems
  - Monitor advances in quantum computing
`);
