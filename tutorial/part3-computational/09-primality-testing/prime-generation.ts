/**
 * # Generating Random Primes
 *
 * Generating large random primes is crucial for cryptography:
 * - RSA needs two large primes p and q
 * - DSA/ElGamal need a large prime modulus
 * - Diffie-Hellman needs a safe prime or Sophie Germain prime
 *
 * ## The Prime Number Theorem
 *
 * The number of primes less than n is approximately n / ln(n).
 *
 * This means the probability that a random n-bit number is prime is
 * about 1 / (n * ln(2)) ~ 1.44 / n.
 *
 * For a 1024-bit prime: ~1 in 710 odd numbers is prime.
 *
 * ## Algorithm
 *
 * 1. Generate a random n-bit odd number
 * 2. Test for primality (Miller-Rabin)
 * 3. If not prime, try the next odd number (or generate new random)
 * 4. Repeat until a prime is found
 *
 * Expected number of trials: O(n) where n is the bit length.
 */

import { is_prime, next_prime, previous_prime, prime_range, is_prime_power, factor } from 'sagemath-ts';

console.log("=== Generating Random Primes ===\n");

// Generate a random bigint with specified number of bits
function randomBigInt(bits: number): bigint {
  // Ensure the high bit is set (for exact bit length)
  let result = 1n << BigInt(bits - 1);

  // Fill remaining bits randomly
  for (let i = 0; i < bits - 1; i++) {
    if (Math.random() < 0.5) {
      result |= 1n << BigInt(i);
    }
  }

  return result;
}

// Generate a random odd number with specified bits
function randomOddNumber(bits: number): bigint {
  const n = randomBigInt(bits);
  return n | 1n; // Ensure odd
}

// Example 1: Prime density
console.log("Example 1: Prime Density (Prime Number Theorem)");
console.log("-".repeat(50));

const ranges: [bigint, bigint, string][] = [
  [1n, 100n, "1-100"],
  [1n, 1000n, "1-1000"],
  [1n, 10000n, "1-10000"],
  [10000n, 11000n, "10000-11000"],
  [100000n, 101000n, "100000-101000"],
];

console.log("  Range          | Primes | Density | 1/ln(mid)");
console.log("  " + "-".repeat(50));

for (const [start, end, label] of ranges) {
  const primes = prime_range(start, end);
  const count = primes.length;
  const rangeSize = Number(end - start);
  const density = (count / rangeSize * 100).toFixed(2);
  const mid = Number((start + end) / 2n);
  const theoretical = (100 / Math.log(mid)).toFixed(2);
  console.log(`  ${label.padEnd(14)} | ${count.toString().padStart(6)} | ${density.padStart(6)}% | ${theoretical}%`);
}

// Example 2: Simple prime generation
console.log("\nExample 2: Simple Prime Generation");
console.log("-".repeat(50));

function generatePrime(bits: number): { prime: bigint; attempts: number } {
  let attempts = 0;
  let candidate: bigint;

  do {
    candidate = randomOddNumber(bits);
    attempts++;
  } while (!is_prime(candidate));

  return { prime: candidate, attempts };
}

console.log("  Generating random primes of various sizes:");
for (const bits of [16, 20, 24, 28, 32]) {
  const results: number[] = [];
  for (let i = 0; i < 5; i++) {
    const { prime, attempts } = generatePrime(bits);
    results.push(attempts);
  }
  const avgAttempts = (results.reduce((a, b) => a + b, 0) / results.length).toFixed(1);
  console.log(`    ${bits}-bit prime: ~${avgAttempts} attempts (expected: ~${(bits * 0.693).toFixed(1)})`);
}

// Example 3: Incremental search vs random generation
console.log("\nExample 3: Incremental Search");
console.log("-".repeat(50));

function generatePrimeIncremental(bits: number): { prime: bigint; attempts: number } {
  let candidate = randomOddNumber(bits);
  let attempts = 0;

  while (!is_prime(candidate)) {
    candidate += 2n; // Next odd number
    attempts++;
  }

  return { prime: candidate, attempts };
}

console.log("  Comparing random vs incremental search:");
for (const bits of [24, 28, 32]) {
  // Random approach
  let randomTotal = 0;
  for (let i = 0; i < 10; i++) {
    randomTotal += generatePrime(bits).attempts;
  }

  // Incremental approach
  let incrTotal = 0;
  for (let i = 0; i < 10; i++) {
    incrTotal += generatePrimeIncremental(bits).attempts;
  }

  console.log(`    ${bits}-bit: Random avg ${(randomTotal / 10).toFixed(1)}, Incremental avg ${(incrTotal / 10).toFixed(1)}`);
}

// Example 4: Using next_prime
console.log("\nExample 4: Using next_prime");
console.log("-".repeat(50));

const startPoints = [100n, 1000n, 10000n, 100000n];

console.log("  Finding next prime after various starting points:");
for (const start of startPoints) {
  const prime = next_prime(start);
  const gap = prime - start;
  console.log(`    next_prime(${start}) = ${prime} (gap: ${gap})`);
}

// Example 5: Safe primes and Sophie Germain primes
console.log("\nExample 5: Safe Primes and Sophie Germain Primes");
console.log("-".repeat(50));

console.log(`
  Sophie Germain prime: p such that 2p + 1 is also prime
  Safe prime: q = 2p + 1 where p is also prime

  These are important for:
  - Diffie-Hellman (safe primes ensure large prime subgroup)
  - RSA (some implementations)
  - Pohlig-Hellman resistance
`);

function isSophieGermain(p: bigint): boolean {
  return is_prime(p) && is_prime(2n * p + 1n);
}

function isSafePrime(q: bigint): boolean {
  if (!is_prime(q)) return false;
  const p = (q - 1n) / 2n;
  return is_prime(p);
}

console.log("  Sophie Germain primes under 100:");
const sgPrimes: bigint[] = [];
for (let p = 2n; p < 100n; p = next_prime(p)) {
  if (isSophieGermain(p)) {
    sgPrimes.push(p);
  }
}
console.log(`    ${sgPrimes.join(', ')}`);

console.log("\n  Corresponding safe primes:");
for (const p of sgPrimes) {
  const q = 2n * p + 1n;
  console.log(`    p = ${p} -> q = ${q}`);
}

// Example 6: Strong primes for RSA
console.log("\nExample 6: Strong Primes for RSA");
console.log("-".repeat(50));

console.log(`
  A "strong prime" p has:
  1. p - 1 has a large prime factor r
  2. p + 1 has a large prime factor s
  3. r - 1 has a large prime factor t

  This helps resist:
  - Pollard p-1 attack (property 1)
  - Pollard p+1 attack (property 2)
  - Cycling attacks (property 3)

  Modern practice: Use primes large enough that these
  attacks are infeasible anyway (2048+ bits).
`);

// Example 7: Prime generation for specific crypto systems
console.log("Example 7: Primes for Different Crypto Systems");
console.log("-".repeat(50));

console.log(`
  RSA (2048-bit modulus):
  - Need two ~1024-bit primes p and q
  - p != q (obviously)
  - |p - q| should be large (prevents Fermat factoring)
  - Modern: random primes are fine

  DSA/ElGamal:
  - Need prime p and prime q | (p-1)
  - p is typically 2048-3072 bits
  - q is typically 256 bits
  - p = kq + 1 for some k

  Diffie-Hellman:
  - Option 1: Safe prime p = 2q + 1
  - Option 2: Schnorr groups with small prime-order subgroup
`);

// Example 8: Generate a prime with specific properties
console.log("Example 8: DSA-style Prime Generation");
console.log("-".repeat(50));

function generateDSAPrimes(qBits: number, pBits: number): { p: bigint; q: bigint } | null {
  // Generate q first (the smaller prime)
  const { prime: q } = generatePrime(qBits);

  // Now find p = kq + 1 that is prime
  const minK = (1n << BigInt(pBits - 1)) / q;
  const maxK = (1n << BigInt(pBits)) / q;

  for (let k = minK; k <= maxK; k++) {
    const p = k * q + 1n;

    // Check bit length
    const pBitsActual = p.toString(2).length;
    if (pBitsActual < pBits) continue;
    if (pBitsActual > pBits) break;

    if (is_prime(p)) {
      return { p, q };
    }
  }

  return null;
}

console.log("  Generating DSA-style primes (q | p-1):");
for (let i = 0; i < 3; i++) {
  const result = generateDSAPrimes(16, 32);
  if (result) {
    const { p, q } = result;
    const k = (p - 1n) / q;
    console.log(`    q = ${q} (${q.toString(2).length} bits)`);
    console.log(`    p = ${p} = ${k} * q + 1 (${p.toString(2).length} bits)`);
    console.log(`    Verify: (p-1) mod q = ${(p - 1n) % q}`);
    console.log();
  }
}

// Example 9: Primorial primes
console.log("Example 9: Special Prime Forms");
console.log("-".repeat(50));

console.log(`
  Some special forms of primes:

  Mersenne primes: 2^p - 1 (p prime)
  - Examples: 3, 7, 31, 127, 8191...
  - Efficient for modular reduction

  Proth primes: k * 2^n + 1 (k < 2^n)
  - Fast primality testing available
  - Used in some crypto implementations

  Generalized Fermat primes: a^(2^n) + 1
  - F_n = 2^(2^n) + 1 are Fermat numbers
  - Only F_0 through F_4 are known to be prime
`);

// Check some Mersenne numbers
console.log("  Mersenne numbers M_p = 2^p - 1:");
for (const p of [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n]) {
  const mp = (1n << p) - 1n;
  const isPr = is_prime(mp);
  console.log(`    M_${p} = ${mp} ${isPr ? '(prime)' : '(composite)'}`);
}

console.log("\n=== Summary ===");
console.log(`
  Prime Generation in Practice:

  1. Prime Number Theorem: ~1 in n*ln(2) n-bit odd numbers is prime
     - Expected O(n) attempts for n-bit prime

  2. Methods:
     - Random + Miller-Rabin: Simple and secure
     - Incremental search: Slightly biased but often faster
     - Sieving: Pre-filter obvious composites

  3. Special primes for cryptography:
     - Safe primes: For DH, ensure large subgroup order
     - Sophie Germain primes: p where 2p+1 is also prime
     - DSA primes: q | (p-1) structure
     - Strong primes: Resist specific factoring attacks

  4. Security considerations:
     - Use cryptographic RNG (not Math.random in production!)
     - Sufficient bit length (2048+ for RSA)
     - Check primality with high confidence (20+ MR rounds)
`);
