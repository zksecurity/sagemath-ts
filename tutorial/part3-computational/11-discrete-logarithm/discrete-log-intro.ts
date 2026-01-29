/**
 * # The Discrete Logarithm Problem (DLP)
 *
 * The discrete logarithm problem is fundamental to many cryptographic
 * systems including Diffie-Hellman key exchange, ElGamal encryption,
 * DSA signatures, and elliptic curve cryptography.
 *
 * ## Definition
 *
 * Given a cyclic group G with generator g and an element h in G,
 * find x such that g^x = h.
 *
 * We write x = log_g(h) or x = dlog_g(h).
 *
 * ## Examples of Groups
 *
 * 1. Multiplicative group (Z/pZ)*: g^x mod p = h
 * 2. Elliptic curve group: [x]P = Q (scalar multiplication)
 * 3. Any cyclic group of prime order
 *
 * ## Why It's Hard
 *
 * Exponentiation is easy: g^x mod p takes O(log x) multiplications
 * Discrete log is hard: Best known takes sub-exponential time
 *
 * This asymmetry is the foundation of DL-based cryptography.
 */

import { power_mod, gcd, is_prime, euler_phi, factor, formatFactorization, Mod, Zmod, GF } from 'sagemath-ts';

console.log("=== The Discrete Logarithm Problem ===\n");

// Example 1: Definition and small examples
console.log("Example 1: What is the Discrete Logarithm?");
console.log("-".repeat(50));

const p = 23n;
const g = 5n; // Generator of (Z/23Z)*

console.log(`  Working in (Z/${p}Z)* with generator g = ${g}`);
console.log(`  The group has order ${p - 1n} (since p is prime)`);
console.log(`\n  Powers of g:`);

const powers = new Map<bigint, bigint>();
let power = 1n;
for (let x = 0n; x < p - 1n; x++) {
  powers.set(power, x);
  if (x < 12n) {
    console.log(`    ${g}^${x} = ${power} mod ${p}`);
  } else if (x === 12n) {
    console.log(`    ...`);
  }
  power = (power * g) % p;
}

console.log(`\n  Discrete logarithm examples:`);
for (const h of [8n, 10n, 16n, 1n]) {
  const x = powers.get(h);
  console.log(`    log_${g}(${h}) = ${x} (since ${g}^${x} = ${h} mod ${p})`);
}

// Example 2: Why DLP is hard
console.log("\nExample 2: The Asymmetry");
console.log("-".repeat(50));

console.log(`
  Forward direction (exponentiation):
  - Compute g^x mod p using square-and-multiply
  - Time: O(log x) multiplications
  - Example: Computing 5^1000000 mod 23 is instant

  Reverse direction (discrete log):
  - Given h, find x such that g^x = h mod p
  - Naive: Try x = 0, 1, 2, ... until match
  - Time: O(p) in worst case

  For cryptographic sizes:
  - p is a 256-bit to 2048-bit prime
  - Naive search would take 2^256 steps
  - Better algorithms exist, but still hard
`);

// Demonstrate fast exponentiation
console.log("  Fast exponentiation demonstration:");
const bigExp = 1000000n;
const result = power_mod(g, bigExp, p);
console.log(`    ${g}^${bigExp} mod ${p} = ${result} (computed instantly)`);

// Example 3: Groups and generators
console.log("\nExample 3: Groups and Generators");
console.log("-".repeat(50));

function isGenerator(g: bigint, p: bigint): boolean {
  // g is a generator of (Z/pZ)* iff g^((p-1)/q) != 1 for all prime q | (p-1)
  const order = p - 1n;
  const primeFactors = factor(order).filter(([pr]) => pr > 0n).map(([pr]) => pr);

  for (const q of primeFactors) {
    const exp = order / q;
    if (power_mod(g, exp, p) === 1n) {
      return false;
    }
  }
  return true;
}

console.log(`  For prime p, (Z/pZ)* is cyclic of order p-1.`);
console.log(`  g is a generator (primitive root) if g has order p-1.`);

console.log(`\n  Checking generators mod ${p}:`);
for (let g = 2n; g <= 10n; g++) {
  const isGen = isGenerator(g, p);
  const ord = multiplicativeOrder(g, p);
  console.log(`    g = ${g}: order = ${ord}, generator: ${isGen}`);
}

// Helper function for multiplicative order
function multiplicativeOrder(a: bigint, n: bigint): bigint {
  if (gcd(a, n) !== 1n) return 0n;

  let order = 1n;
  let current = a % n;

  while (current !== 1n) {
    current = (current * a) % n;
    order++;
    if (order > n) return -1n; // Shouldn't happen for valid inputs
  }

  return order;
}

// Example 4: DLP in different groups
console.log("\nExample 4: DLP in Different Groups");
console.log("-".repeat(50));

console.log(`
  The DLP can be defined in any group:

  1. (Z/pZ)* - Multiplicative group mod prime p
     - g^x = h (mod p)
     - Used in: Original DH, DSA, ElGamal

  2. (Z/nZ)* - Multiplicative group mod composite n
     - Harder to analyze (related to factoring)
     - Less commonly used

  3. Elliptic curve groups E(F_p)
     - [x]P = Q (scalar multiplication)
     - Used in: ECDH, ECDSA, EdDSA
     - Same security with smaller keys

  4. Subgroups of prime order
     - Work in a subgroup of order q where q | p-1
     - Used in: Schnorr groups, DSA
`);

// Example 5: Schnorr groups (prime order subgroups)
console.log("Example 5: Schnorr Groups (Prime Order Subgroups)");
console.log("-".repeat(50));

console.log(`
  For efficiency and security, we often work in a prime-order
  subgroup of (Z/pZ)*.

  Construction:
  1. Choose a prime q (e.g., 256 bits)
  2. Find p = kq + 1 that is also prime
  3. Find g such that g has order q (not order p-1)

  The subgroup {1, g, g^2, ..., g^(q-1)} has prime order q.
`);

// Small example
const q_small = 11n;
const p_small = 2n * q_small + 1n; // p = 23

console.log(`  Example: q = ${q_small}, p = 2q + 1 = ${p_small}`);
console.log(`  (Z/${p_small}Z)* has order ${p_small - 1n} = 2 * ${q_small}`);

// Find a generator of order q
let g_q = 2n;
while (power_mod(g_q, q_small, p_small) !== 1n) {
  g_q++;
}
// Actually, we need g^((p-1)/q) = g^2 to be the generator of subgroup
const g_subgroup = power_mod(g_q, 2n, p_small);

console.log(`\n  Generator of subgroup of order ${q_small}: g = ${g_subgroup}`);
console.log(`  Verification: g^${q_small} = ${power_mod(g_subgroup, q_small, p_small)} (should be 1)`);

console.log(`\n  Subgroup elements: {1, g, g^2, ..., g^(q-1)}`);
let elem = 1n;
const subgroup: bigint[] = [];
for (let i = 0n; i < q_small; i++) {
  subgroup.push(elem);
  elem = (elem * g_subgroup) % p_small;
}
console.log(`    ${subgroup.join(", ")}`);

// Example 6: DLP difficulty depends on group order
console.log("\nExample 6: Group Order and DLP Difficulty");
console.log("-".repeat(50));

console.log(`
  The hardness of DLP depends on the GROUP ORDER, not the modulus.

  For (Z/pZ)*:
  - Group order is p - 1
  - If p - 1 = 2^k * (small primes), DLP can be solved efficiently
    via Pohlig-Hellman algorithm

  For security:
  - Want group order to have a LARGE PRIME FACTOR
  - Safe primes p = 2q + 1 give order 2q (q is large prime)
  - Schnorr groups work in subgroup of prime order q

  Example of weak group:
`);

// Weak example: p - 1 is smooth
const p_weak = 127n; // p - 1 = 126 = 2 * 3^2 * 7
console.log(`  p = ${p_weak}, p - 1 = ${p_weak - 1n} = ${formatFactorization(factor(p_weak - 1n))}`);
console.log(`  Largest prime factor: 7 (very weak!)`);
console.log(`  DLP can be solved by Pohlig-Hellman in O(sqrt(7)) ~ 3 steps.`);

// Strong example
const p_strong = 23n; // p - 1 = 22 = 2 * 11
console.log(`\n  p = ${p_strong}, p - 1 = ${p_strong - 1n} = ${formatFactorization(factor(p_strong - 1n))}`);
console.log(`  Largest prime factor: 11 (better)`);
console.log(`  DLP requires O(sqrt(11)) ~ 4 steps.`);

// Example 7: DLP vs factoring
console.log("\nExample 7: DLP vs Factoring");
console.log("-".repeat(50));

console.log(`
  Both are believed to be hard, but they're different problems:

  Factoring:
  - Input: n = p * q
  - Output: p and q
  - Best algorithm: Number Field Sieve

  Discrete Log in (Z/pZ)*:
  - Input: g, h, p where g^x = h mod p
  - Output: x
  - Best algorithm: Number Field Sieve (different variant)

  Relationship:
  - If you can factor, you CAN'T necessarily solve DLP
  - If you can solve DLP mod p, you CAN'T necessarily factor
  - They're not known to be equivalent!

  For elliptic curves:
  - No sub-exponential algorithm known
  - Makes ECDLP harder than DLP in (Z/pZ)*
  - 256-bit EC ~ 3072-bit RSA/DH security
`);

// Example 8: Computing discrete logs for small examples
console.log("\nExample 8: Solving Small DLP Instances");
console.log("-".repeat(50));

function discreteLogNaive(g: bigint, h: bigint, p: bigint): bigint | null {
  const order = p - 1n;
  let current = 1n;

  for (let x = 0n; x < order; x++) {
    if (current === h) {
      return x;
    }
    current = (current * g) % p;
  }

  return null;
}

console.log(`  Solving discrete logs mod ${p} with generator ${g}:`);

const targets = [1n, 2n, 5n, 8n, 10n, 16n, 21n];
for (const h of targets) {
  const x = discreteLogNaive(g, h, p);
  if (x !== null) {
    const verify = power_mod(g, x, p);
    console.log(`    log_${g}(${h}) = ${x} (verify: ${g}^${x} = ${verify})`);
  } else {
    console.log(`    log_${g}(${h}) = not found (${h} not in subgroup)`);
  }
}

// Example 9: Applications
console.log("\nExample 9: Cryptographic Applications");
console.log("-".repeat(50));

console.log(`
  DLP hardness enables many cryptographic constructions:

  1. Diffie-Hellman Key Exchange:
     - Alice: a, A = g^a
     - Bob: b, B = g^b
     - Shared: K = A^b = B^a = g^(ab)
     - Attacker sees g, A, B but can't find K without DLP

  2. ElGamal Encryption:
     - Public key: (p, g, h = g^x)
     - Encrypt m: (c1, c2) = (g^k, m * h^k)
     - Decrypt: m = c2 / c1^x
     - Breaking requires DLP

  3. DSA/Schnorr Signatures:
     - Private key: x
     - Public key: y = g^x
     - Signature proves knowledge of x

  4. Zero-Knowledge Proofs:
     - Many ZK constructions rely on DLP
     - Schnorr identification protocol
`);

// Demonstrate DH key exchange
console.log("  Simple Diffie-Hellman demonstration:");
const alice_secret = 6n;
const bob_secret = 15n;

const alice_public = power_mod(g, alice_secret, p);
const bob_public = power_mod(g, bob_secret, p);

const alice_shared = power_mod(bob_public, alice_secret, p);
const bob_shared = power_mod(alice_public, bob_secret, p);

console.log(`    Alice: a = ${alice_secret}, A = g^a = ${alice_public}`);
console.log(`    Bob:   b = ${bob_secret}, B = g^b = ${bob_public}`);
console.log(`    Alice computes: B^a = ${alice_shared}`);
console.log(`    Bob computes:   A^b = ${bob_shared}`);
console.log(`    Shared secret: ${alice_shared} (matches!)`);

console.log("\n=== Summary ===");
console.log(`
  The Discrete Logarithm Problem:

  Definition: Given g, h in a group G, find x such that g^x = h.

  Key properties:
  - Exponentiation is efficient (O(log x) operations)
  - Discrete log is hard (best: sub-exponential)
  - This asymmetry enables cryptography

  Important groups:
  - (Z/pZ)*: Traditional DH, DSA, ElGamal
  - Prime-order subgroups: Schnorr groups
  - Elliptic curves: ECDH, ECDSA, EdDSA

  Security depends on:
  - Group order having large prime factors
  - Sufficient group size (256+ bits for EC, 2048+ for (Z/pZ)*)

  Next: Algorithms for solving DLP (brute force, baby-giant, Pohlig-Hellman)
`);
