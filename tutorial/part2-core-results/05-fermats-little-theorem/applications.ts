/**
 * # Cryptographic Applications of Fermat's Little Theorem
 *
 * Fermat's Little Theorem has numerous applications in cryptography:
 *
 * 1. Fast modular inverse computation
 * 2. Primality testing (Fermat test)
 * 3. RSA encryption/decryption
 * 4. Pseudorandom number generation
 * 5. Hash functions and checksums
 *
 * This tutorial explores these applications with practical examples.
 *
 * @module tutorial/part2/fermat-applications
 */

import {
  power_mod,
  inverse_mod,
  gcd,
  is_prime,
  factor,
  euler_phi,
  Zmod,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("CRYPTOGRAPHIC APPLICATIONS OF FERMAT'S THEOREM");
console.log("=".repeat(70));

// ============================================================================
// Application 1: Fast Modular Inverse
// ============================================================================
console.log("\n=== Application 1: Fast Modular Inverse ===\n");

console.log("For prime p and a coprime to p:");
console.log("  a^(-1) = a^(p-2) mod p\n");

// Compare methods
const p = 1000000007n; // Large prime
const a = 123456789n;

console.log(`Computing ${a}^(-1) mod ${p}:\n`);

// Method 1: Extended Euclidean Algorithm
const invEuclid = inverse_mod(a, p);
console.log(`Method 1 (Extended Euclidean): ${invEuclid}`);

// Method 2: Fermat's Little Theorem
const invFermat = power_mod(a, p - 2n, p);
console.log(`Method 2 (Fermat's theorem):   ${invFermat}`);

// Verify
const verify = (a * invFermat) % p;
console.log(`\nVerification: ${a} * ${invFermat} mod ${p} = ${verify}`);
console.log(`Both methods agree: ${invEuclid === invFermat}`);

console.log(`
Note: Extended Euclidean is generally faster, but Fermat's method:
  - Is simpler to implement
  - Can be parallelized more easily
  - Works well when modular exponentiation is already optimized
`);

// ============================================================================
// Application 2: Fermat Primality Test
// ============================================================================
console.log("\n=== Application 2: Fermat Primality Test ===\n");

/**
 * Fermat primality test.
 * If n is prime and gcd(a,n)=1, then a^(n-1) = 1 (mod n).
 * Contrapositive: If a^(n-1) != 1 (mod n), then n is composite.
 */
function fermatPrimalityTest(n: bigint, k: number = 10): {
  isProbablePrime: boolean;
  witnesses: bigint[];
  failedWitness?: bigint;
} {
  if (n < 2n) return { isProbablePrime: false, witnesses: [] };
  if (n === 2n || n === 3n) return { isProbablePrime: true, witnesses: [] };
  if (n % 2n === 0n) return { isProbablePrime: false, witnesses: [] };

  const witnesses: bigint[] = [];

  for (let i = 0; i < k; i++) {
    // Random witness in [2, n-2]
    const a = 2n + BigInt(Math.floor(Math.random() * Number(n - 4n > 1000n ? 1000n : n - 4n)));
    witnesses.push(a);

    if (gcd(a, n) !== 1n) {
      // Found a factor!
      return { isProbablePrime: false, witnesses, failedWitness: a };
    }

    const result = power_mod(a, n - 1n, n);
    if (result !== 1n) {
      return { isProbablePrime: false, witnesses, failedWitness: a };
    }
  }

  return { isProbablePrime: true, witnesses };
}

console.log("Testing various numbers with Fermat primality test:\n");

const testCases: Array<{ n: bigint; description: string }> = [
  { n: 7n, description: "small prime" },
  { n: 15n, description: "composite (3*5)" },
  { n: 97n, description: "prime" },
  { n: 561n, description: "Carmichael number" },
  { n: 1009n, description: "prime" },
  { n: 1729n, description: "Carmichael (Hardy-Ramanujan)" },
  { n: 104729n, description: "10000th prime" },
];

for (const { n, description } of testCases) {
  const test = fermatPrimalityTest(n, 5);
  const actual = is_prime(n);
  const status = test.isProbablePrime === actual ? "OK" :
                 (test.isProbablePrime && !actual ? "FOOLED" : "OK");
  console.log(`  n = ${n.toString().padStart(8)} (${description.padEnd(25)}) ` +
              `Fermat: ${test.isProbablePrime.toString().padEnd(5)} ` +
              `Actual: ${actual.toString().padEnd(5)} [${status}]`);
}

console.log(`
Important: Carmichael numbers fool the Fermat test!
They satisfy a^(n-1) = 1 (mod n) for all a coprime to n, but are composite.
Examples: 561 = 3*11*17, 1729 = 7*13*19 (taxicab number)

For production use, prefer Miller-Rabin test (available in is_prime).
`);

// ============================================================================
// Application 3: Diffie-Hellman Key Exchange
// ============================================================================
console.log("\n=== Application 3: Diffie-Hellman Key Exchange ===\n");

console.log("Diffie-Hellman allows two parties to establish a shared secret");
console.log("over an insecure channel.\n");

// Small numbers for demonstration (use much larger in practice!)
const dhP = 23n; // Prime modulus
const dhG = 5n;  // Generator

console.log(`Public parameters: p = ${dhP}, g = ${dhG}`);

// Alice's private and public key
const alicePrivate = 6n;
const alicePublic = power_mod(dhG, alicePrivate, dhP);

// Bob's private and public key
const bobPrivate = 15n;
const bobPublic = power_mod(dhG, bobPrivate, dhP);

console.log(`\nAlice: private key a = ${alicePrivate}, public key A = g^a = ${alicePublic}`);
console.log(`Bob:   private key b = ${bobPrivate}, public key B = g^b = ${bobPublic}`);

// Shared secret computation
const aliceShared = power_mod(bobPublic, alicePrivate, dhP);   // B^a = (g^b)^a
const bobShared = power_mod(alicePublic, bobPrivate, dhP);     // A^b = (g^a)^b

console.log(`\nShared secret computation:`);
console.log(`  Alice computes: B^a = ${bobPublic}^${alicePrivate} mod ${dhP} = ${aliceShared}`);
console.log(`  Bob computes:   A^b = ${alicePublic}^${bobPrivate} mod ${dhP} = ${bobShared}`);
console.log(`  Shared secret:  ${aliceShared} (both compute g^(ab) = g^${alicePrivate * bobPrivate} mod p)`);

console.log(`
Security relies on the Discrete Logarithm Problem (DLP):
Given g, p, and g^x mod p, it's hard to find x.
Fermat's theorem ensures g^(p-1) = 1, so exponents cycle with period dividing p-1.
`);

// ============================================================================
// Application 4: ElGamal Encryption
// ============================================================================
console.log("\n=== Application 4: ElGamal Encryption ===\n");

console.log("ElGamal encryption builds on Diffie-Hellman:\n");

// Parameters (small for demonstration)
const elP = 43n;
const elG = 3n;

// Key generation
const elPrivate = 7n; // Private key x
const elPublic = power_mod(elG, elPrivate, elP); // Public key h = g^x

console.log(`Key generation:`);
console.log(`  p = ${elP}, g = ${elG}`);
console.log(`  Private key x = ${elPrivate}`);
console.log(`  Public key h = g^x = ${elPublic}`);

// Encryption
const message = 15n;
const ephemeral = 11n; // Random ephemeral key k
const c1 = power_mod(elG, ephemeral, elP);              // c1 = g^k
const sharedSecret = power_mod(elPublic, ephemeral, elP); // s = h^k = g^(xk)
const c2 = (message * sharedSecret) % elP;              // c2 = m * s

console.log(`\nEncryption of message m = ${message}:`);
console.log(`  Ephemeral key k = ${ephemeral}`);
console.log(`  c1 = g^k = ${c1}`);
console.log(`  s = h^k = ${sharedSecret} (shared secret)`);
console.log(`  c2 = m * s = ${c2}`);
console.log(`  Ciphertext: (${c1}, ${c2})`);

// Decryption
const decryptShared = power_mod(c1, elPrivate, elP);    // s = c1^x = g^(kx)
const decryptSharedInv = power_mod(decryptShared, elP - 2n, elP); // s^(-1) using Fermat
const decrypted = (c2 * decryptSharedInv) % elP;        // m = c2 * s^(-1)

console.log(`\nDecryption:`);
console.log(`  s = c1^x = ${decryptShared}`);
console.log(`  s^(-1) = ${decryptSharedInv} (using Fermat: s^(p-2))`);
console.log(`  m = c2 * s^(-1) = ${decrypted}`);
console.log(`  Correct: ${decrypted === message}`);

// ============================================================================
// Application 5: Pseudorandom Number Generation
// ============================================================================
console.log("\n=== Application 5: Blum Blum Shub PRNG ===\n");

console.log("The BBS generator uses modular squaring:");
console.log("  x_{n+1} = x_n^2 mod M, where M = p*q (p,q = 3 mod 4)\n");

// Small primes for demonstration
const bbsP = 11n; // 11 = 3 mod 4
const bbsQ = 19n; // 19 = 3 mod 4
const bbsM = bbsP * bbsQ;
let bbsState = 3n; // Seed (must be coprime to M)

console.log(`Parameters: p = ${bbsP}, q = ${bbsQ}, M = ${bbsM}`);
console.log(`Seed: ${bbsState}`);
console.log(`\nGenerating sequence (showing x_n mod 2 as random bits):`);

const bits: number[] = [];
for (let i = 0; i < 20; i++) {
  bbsState = power_mod(bbsState, 2n, bbsM);
  const bit = Number(bbsState % 2n);
  bits.push(bit);
  if (i < 10) {
    console.log(`  x_${(i + 1).toString().padStart(2)} = ${bbsState.toString().padStart(3)}, bit = ${bit}`);
  }
}

console.log(`  ...`);
console.log(`\nBit sequence: ${bits.join("")}`);
console.log(`
Security: Predicting BBS output is as hard as factoring M.
Fermat's theorem ensures the sequence eventually cycles.
Period divides lambda(M) = lcm(p-1, q-1).
`);

// ============================================================================
// Application 6: Polynomial Evaluation in Finite Fields
// ============================================================================
console.log("\n=== Application 6: Polynomial Evaluation in Finite Fields ===\n");

console.log("Evaluating polynomials over Z/pZ uses Fermat's theorem:");
console.log("Since a^p = a (mod p), we have a^n = a^(n mod (p-1)) for a != 0.\n");

const polyP = 7n;
const Zp = Zmod(polyP);

// Polynomial: f(x) = x^100 + 3x^50 + 2x + 1
function evaluatePolynomial(x: bigint): bigint {
  // Reduce exponents using Fermat
  // x^100 mod 7: 100 mod 6 = 4, so x^100 = x^4
  // x^50 mod 7: 50 mod 6 = 2, so x^50 = x^2
  const x100 = power_mod(x, 100n, polyP); // Library handles optimization
  const x50 = power_mod(x, 50n, polyP);
  return ((x100 + 3n * x50 + 2n * x + 1n) % polyP + polyP) % polyP;
}

console.log(`Polynomial f(x) = x^100 + 3x^50 + 2x + 1 over Z/${polyP}Z:`);
console.log(`(Note: x^100 mod p uses x^(100 mod (p-1)) = x^${100n % (polyP - 1n)} when x != 0)\n`);

console.log("x | f(x)");
console.log("--+-----");
for (let x = 0n; x < polyP; x++) {
  const fx = evaluatePolynomial(x);
  console.log(`${x} | ${fx}`);
}

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
Cryptographic Applications of Fermat's Little Theorem:

1. Modular Inverse: a^(-1) = a^(p-2) mod p
   - Simple alternative to extended Euclidean algorithm

2. Primality Testing: If a^(n-1) != 1 mod n, then n is composite
   - Fast but fooled by Carmichael numbers
   - Use Miller-Rabin for production

3. Diffie-Hellman: Shared secret via g^(ab) mod p
   - Security relies on discrete log problem
   - Fermat ensures cyclic structure

4. ElGamal Encryption: Public-key encryption
   - Uses Fermat for decryption (computing inverse)

5. BBS PRNG: Cryptographically secure randomness
   - Period analysis uses Fermat/Euler

6. Polynomial Evaluation: Exponent reduction
   - x^n = x^(n mod (p-1)) for x != 0

Key Insight: Fermat's theorem creates predictable cyclic structure
in the multiplicative group (Z/pZ)*, which is both useful
(for computation) and potentially dangerous (for security analysis).
`);
