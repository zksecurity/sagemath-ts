/**
 * # Diffie-Hellman Security
 *
 * This tutorial explores the security foundations of Diffie-Hellman:
 * the Discrete Logarithm Problem (DLP), Computational Diffie-Hellman (CDH),
 * and Decisional Diffie-Hellman (DDH) assumptions.
 *
 * ## Security Assumptions Hierarchy
 *
 * DLP <= CDH <= DDH
 *
 * If DLP is easy, then CDH is easy.
 * If CDH is easy, then DDH is easy.
 * (Reverse implications are not known to hold in general.)
 *
 * ## Practical Implications
 *
 * Different cryptographic constructions rely on different assumptions.
 * Understanding these helps choose appropriate parameters and protocols.
 */

import { gcd, is_prime, power_mod, factor, euler_phi, Mod } from 'sagemath-ts';

console.log('='.repeat(60));
console.log('Diffie-Hellman Security Tutorial');
console.log('='.repeat(60));

// ============================================================================
// Setup: Group Parameters
// ============================================================================

console.log('\n=== Group Parameters ===\n');

// Use a prime where p-1 has nice structure for demonstration
const p = 23n;
const g = 5n; // Generator of (Z/pZ)*

console.log(`Prime: p = ${p}`);
console.log(`Generator: g = ${g}`);
console.log(`Group order: p - 1 = ${p - 1n}`);
console.log(`Factorization of p-1: ${factor(p - 1n).map(([f, e]) => e === 1n ? `${f}` : `${f}^${e}`).join(' * ')}`);

// Verify g is a generator
function orderOfElement(g: bigint, p: bigint): bigint {
  let order = 1n;
  let current = g;
  while (current !== 1n) {
    current = (current * g) % p;
    order++;
  }
  return order;
}

const gOrder = orderOfElement(g, p);
console.log(`\nOrder of g: ${gOrder}`);
console.log(`g is a generator: ${gOrder === p - 1n}`);

// ============================================================================
// The Discrete Logarithm Problem (DLP)
// ============================================================================

console.log('\n=== Discrete Logarithm Problem (DLP) ===\n');

/**
 * DLP Definition:
 *
 * Given: Prime p, generator g, and h = g^x mod p
 * Find: x (the discrete logarithm of h base g)
 *
 * This is believed to be computationally hard for large p.
 */

console.log('DLP Definition:');
console.log('  Given: p, g, and h = g^x mod p');
console.log('  Find: x');
console.log('');

// Example DLP instance
const x_secret = 7n;
const h = Mod(g, p).pow(x_secret).value;

console.log('DLP Instance:');
console.log(`  p = ${p}`);
console.log(`  g = ${g}`);
console.log(`  h = g^x = ${h}`);
console.log(`  Challenge: Find x such that ${g}^x = ${h} (mod ${p})`);
console.log('');

// Solve by brute force (only works for small p)
console.log('Brute force solution:');
for (let x = 0n; x < p; x++) {
  if (Mod(g, p).pow(x).value === h) {
    console.log(`  Found: x = ${x}`);
    console.log(`  Verification: ${g}^${x} mod ${p} = ${Mod(g, p).pow(x).value}`);
    break;
  }
}

// ============================================================================
// DLP Algorithms
// ============================================================================

console.log('\n=== DLP Algorithms ===\n');

/**
 * Methods for solving DLP (with increasing sophistication):
 *
 * 1. Brute Force: O(p) - try all exponents
 * 2. Baby-Step Giant-Step: O(sqrt(p)) time and space
 * 3. Pollard's Rho: O(sqrt(p)) time, O(1) space
 * 4. Index Calculus: Sub-exponential for prime fields
 * 5. Number Field Sieve: Best known for large primes
 */

console.log('Algorithm Complexities:');
console.log('');
console.log('  Algorithm               Time           Space');
console.log('  ---------               ----           -----');
console.log('  Brute Force             O(p)           O(1)');
console.log('  Baby-Step Giant-Step    O(sqrt(p))     O(sqrt(p))');
console.log("  Pollard's Rho           O(sqrt(p))     O(1)");
console.log('  Index Calculus          L[1/2, c]      Polynomial');
console.log('  Number Field Sieve      L[1/3, c]      Polynomial');
console.log('');
console.log('  where L[a, c] = exp(c * (log p)^a * (log log p)^(1-a))');

// Demonstrate Baby-Step Giant-Step
console.log('\nBaby-Step Giant-Step Algorithm:');

function babyStepGiantStep(g: bigint, h: bigint, p: bigint): bigint {
  // Compute m = ceil(sqrt(p-1))
  let m = 1n;
  while (m * m < p - 1n) m++;

  console.log(`  m = ceil(sqrt(p-1)) = ${m}`);

  // Baby steps: compute g^j for j = 0, 1, ..., m-1
  const babySteps = new Map<bigint, bigint>();
  let gj = 1n;
  for (let j = 0n; j < m; j++) {
    babySteps.set(gj, j);
    gj = (gj * g) % p;
  }
  console.log(`  Baby steps computed: ${babySteps.size} values`);

  // Giant steps: compute h * (g^(-m))^i for i = 0, 1, ...
  // g^(-m) = (g^(-1))^m
  const gInv = Mod(g, p).inv().value;
  const gMinusM = Mod(gInv, p).pow(m).value;

  let gamma = h;
  for (let i = 0n; i < m; i++) {
    if (babySteps.has(gamma)) {
      const j = babySteps.get(gamma)!;
      const x = i * m + j;
      console.log(`  Match found at i=${i}, j=${j}`);
      return x;
    }
    gamma = (gamma * gMinusM) % p;
  }

  throw new Error('No solution found');
}

const x_found = babyStepGiantStep(g, h, p);
console.log(`  Solution: x = ${x_found}`);
console.log(`  Verification: ${g}^${x_found} mod ${p} = ${Mod(g, p).pow(x_found).value}`);

// ============================================================================
// Computational Diffie-Hellman (CDH)
// ============================================================================

console.log('\n=== Computational Diffie-Hellman (CDH) ===\n');

/**
 * CDH Problem:
 *
 * Given: Prime p, generator g, A = g^a mod p, B = g^b mod p
 * Compute: g^(ab) mod p
 *
 * Note: CDH doesn't require finding a or b!
 *
 * CDH is "weaker" than DLP:
 * - If you can solve DLP, you can solve CDH
 * - Reverse is not known (in general)
 */

console.log('CDH Problem:');
console.log('  Given: p, g, A = g^a, B = g^b');
console.log('  Compute: g^(ab) (without finding a or b)');
console.log('');

// CDH instance
const a = 3n;
const b = 11n;
const A = Mod(g, p).pow(a).value;
const B = Mod(g, p).pow(b).value;
const K = Mod(g, p).pow(a * b).value; // The CDH challenge answer

console.log('CDH Instance:');
console.log(`  p = ${p}, g = ${g}`);
console.log(`  A = g^a = ${A}`);
console.log(`  B = g^b = ${B}`);
console.log(`  Challenge: Compute g^(ab) = ${K}`);
console.log('');

// Show that solving DLP solves CDH
console.log('Solving CDH via DLP:');
console.log('  1. Solve DLP for A: Find a such that g^a = A');

let a_found = 0n;
for (let x = 0n; x < p; x++) {
  if (Mod(g, p).pow(x).value === A) {
    a_found = x;
    break;
  }
}
console.log(`     Found a = ${a_found}`);

console.log(`  2. Compute B^a = g^(ba) = g^(ab)`);
const cdh_solution = Mod(B, p).pow(a_found).value;
console.log(`     B^a = ${B}^${a_found} mod ${p} = ${cdh_solution}`);
console.log(`     Correct: ${cdh_solution === K}`);

// ============================================================================
// Decisional Diffie-Hellman (DDH)
// ============================================================================

console.log('\n=== Decisional Diffie-Hellman (DDH) ===\n');

/**
 * DDH Problem:
 *
 * Given: Prime p, generator g, A = g^a, B = g^b, and C
 * Decide: Is C = g^(ab) or is C random?
 *
 * DDH is "weaker" than CDH:
 * - If you can solve CDH, you can solve DDH (compute g^(ab) and compare)
 * - Reverse is not known
 *
 * DDH is the weakest (easiest) of the three problems.
 */

console.log('DDH Problem:');
console.log('  Given: p, g, A = g^a, B = g^b, and C');
console.log('  Decide: Is C = g^(ab) (real) or C random (fake)?');
console.log('');

// DDH challenger presents tuples
console.log('DDH Challenge Examples:');
console.log('');

// Real tuple (C = g^(ab))
const C_real = K;
console.log(`Real tuple: (g=${g}, A=${A}, B=${B}, C=${C_real})`);
console.log(`  This is a DH tuple: C = g^(ab)`);

// Fake tuple (C is random)
const c_random = 17n;
const C_fake = Mod(g, p).pow(c_random).value;
console.log(`\nFake tuple: (g=${g}, A=${A}, B=${B}, C=${C_fake})`);
console.log(`  This is NOT a DH tuple: C = g^${c_random} is random`);

console.log('\nSolving DDH via CDH:');
console.log('  1. Compute g^(ab) using CDH solver');
console.log('  2. Compare with given C');
console.log('  3. If equal: real tuple. If not: fake tuple.');
console.log(`\n  For real tuple: computed g^(ab) = ${K}, C = ${C_real}, match: ${K === C_real}`);
console.log(`  For fake tuple: computed g^(ab) = ${K}, C = ${C_fake}, match: ${K === C_fake}`);

// ============================================================================
// Security Assumption Hierarchy
// ============================================================================

console.log('\n=== Security Assumption Hierarchy ===\n');

/**
 * The hierarchy from hardest to easiest:
 *
 * DLP >= CDH >= DDH
 *
 * - DLP hard => CDH hard => DDH hard
 * - DDH easy does NOT imply CDH easy
 * - CDH easy does NOT imply DLP easy
 *
 * In some groups, DDH is easy but CDH is hard!
 * Example: Bilinear groups (used in pairing-based crypto)
 */

console.log('Hierarchy:');
console.log('');
console.log('  DLP (hardest)');
console.log('   |');
console.log('   v "If DLP is easy, then CDH is easy"');
console.log('  CDH');
console.log('   |');
console.log('   v "If CDH is easy, then DDH is easy"');
console.log('  DDH (easiest)');
console.log('');
console.log('Implications:');
console.log('  - Protocols based on DDH have strongest security assumption');
console.log('  - Protocols based on DLP have weakest security assumption');
console.log('  - "Stronger assumption" = more ways it could be false');

// ============================================================================
// Which Protocols Use Which Assumptions
// ============================================================================

console.log('\n=== Protocols and Their Assumptions ===\n');

console.log('Protocols relying on DLP:');
console.log('  - Schnorr signatures');
console.log('  - DSA/ECDSA (uses DLP for signing)');
console.log('');

console.log('Protocols relying on CDH:');
console.log('  - Basic Diffie-Hellman key exchange');
console.log('  - Hashed ElGamal encryption');
console.log('');

console.log('Protocols relying on DDH:');
console.log('  - ElGamal encryption (IND-CPA security)');
console.log('  - ECDH key exchange (in prime-order groups)');
console.log('  - Pedersen commitments (hiding property)');
console.log('');

console.log('Note: Using DDH assumption is "risky" in some groups.');
console.log('In groups with pairings, DDH is easy but CDH is still hard.');

// ============================================================================
// Attacks on Weak Parameters
// ============================================================================

console.log('\n=== Attacks on Weak Parameters ===\n');

/**
 * Security depends heavily on parameter choice:
 *
 * 1. Small p: DLP easily solved by baby-step giant-step
 * 2. Smooth p-1: Pohlig-Hellman attack reduces DLP to small subgroups
 * 3. Non-prime-order subgroup: Small subgroup attacks
 * 4. Weak generator: May generate small subgroup
 */

console.log('Attack 1: Small p (DLP is easy)');
console.log(`  Our p = ${p} allows brute force in ${p - 1n} steps`);
console.log('  Solution: Use p with 2048+ bits');
console.log('');

// Pohlig-Hellman for smooth p-1
console.log('Attack 2: Pohlig-Hellman (smooth p-1)');
console.log(`  p - 1 = ${p - 1n} = ${factor(p - 1n).map(([f, e]) => `${f}^${e}`).join(' * ')}`);
console.log('  If p-1 has only small factors, DLP is easy.');
console.log('  Solution: Use p where p-1 = 2q for large prime q');
console.log('');

// Demonstrate Pohlig-Hellman concept
console.log('Pohlig-Hellman concept:');
console.log(`  Let p - 1 = ${p - 1n} = 2 * 11`);
console.log('  We can solve DLP mod 2 and mod 11 separately,');
console.log('  then combine with CRT.');
console.log('  Each smaller DLP is easier than the full one.');

// ============================================================================
// Safe Primes
// ============================================================================

console.log('\n=== Safe Primes ===\n');

/**
 * A safe prime p satisfies: p = 2q + 1 where q is also prime.
 * The prime q is called a Sophie Germain prime.
 *
 * Benefits:
 * - p - 1 = 2q has only one large factor (q)
 * - Pohlig-Hellman gives no advantage
 * - The subgroup of order q has strong security
 */

function isSafePrime(p: bigint): boolean {
  if (!is_prime(p)) return false;
  const q = (p - 1n) / 2n;
  return is_prime(q);
}

console.log('Safe Prime Definition:');
console.log('  p is safe prime if p = 2q + 1 where q is prime');
console.log('');

// Find some safe primes
console.log('Small safe primes:');
const safePrimes: bigint[] = [];
for (let n = 5n; n < 200n && safePrimes.length < 10; n += 2n) {
  if (isSafePrime(n)) {
    safePrimes.push(n);
    const q = (n - 1n) / 2n;
    console.log(`  p = ${n} = 2 * ${q} + 1`);
  }
}

console.log('\nStandard safe primes (RFC 3526):');
console.log('  1536-bit, 2048-bit, 3072-bit, 4096-bit groups defined');
console.log('  These should be used in practice');

// ============================================================================
// Subgroup Security
// ============================================================================

console.log('\n=== Subgroup Security ===\n');

/**
 * When using safe prime p = 2q + 1:
 * - (Z/pZ)* has order p-1 = 2q
 * - Subgroups have order 1, 2, q, or 2q
 * - We work in the order-q subgroup for security
 *
 * Working in order-q subgroup:
 * - Use g = r^2 mod p for random r (squares form order-q subgroup)
 * - Or verify g^q = 1 mod p
 */

// Use safe prime p = 23 (q = 11)
const safePrime = 23n;
const q_sophie = (safePrime - 1n) / 2n;

console.log(`Safe prime: p = ${safePrime}`);
console.log(`Sophie Germain prime: q = ${q_sophie}`);
console.log(`Group order: p - 1 = ${safePrime - 1n} = 2 * ${q_sophie}`);
console.log('');

// Find element of order q
console.log('Finding element of order q:');
for (let r = 2n; r < safePrime; r++) {
  const g_candidate = Mod(r, safePrime).pow(2n).value; // Use square
  const order = orderOfElement(g_candidate, safePrime);
  if (order === q_sophie) {
    console.log(`  r = ${r}, g = r^2 = ${g_candidate}, order = ${order}`);
    console.log(`  Verification: g^q = ${Mod(g_candidate, safePrime).pow(q_sophie).value}`);
    break;
  }
}

// ============================================================================
// Small Subgroup Attack
// ============================================================================

console.log('\n=== Small Subgroup Attack ===\n');

/**
 * If an attacker can force the DH computation into a small subgroup,
 * they can recover the private key modulo the subgroup order.
 *
 * Attack:
 * 1. Send A' = h where h has small order k
 * 2. Victim computes K = A'^b = h^b = h^(b mod k)
 * 3. Attacker can brute-force b mod k in k steps
 * 4. Repeat with different small-order elements to recover full b
 *
 * Defense: Validate received values are in the correct subgroup
 */

console.log('Small Subgroup Attack:');
console.log('');

// In (Z/23Z)*, element of order 2 is -1 = 22
const smallOrderElem = safePrime - 1n; // order 2

console.log(`  Attacker sends A' = ${smallOrderElem} (has order 2)`);
console.log(`  Victim has secret b, computes K = A'^b`);
console.log('');

for (let b_victim = 0n; b_victim < 4n; b_victim++) {
  const K_computed = Mod(smallOrderElem, safePrime).pow(b_victim).value;
  console.log(`    If b = ${b_victim}: K = ${smallOrderElem}^${b_victim} = ${K_computed}`);
}

console.log('');
console.log('  Attacker learns: K = 1 means b is even, K = 22 means b is odd');
console.log('  This leaks 1 bit of information about b!');
console.log('');
console.log('Defense:');
console.log('  - Verify received values: Check A^q = 1');
console.log('  - Use prime-order groups');
console.log('  - Cofactor multiplication');

// ============================================================================
// Summary
// ============================================================================

console.log('\n=== Summary ===\n');

console.log('DH SECURITY KEY POINTS:');
console.log('');
console.log('  1. Problem Hierarchy (hardest to easiest):');
console.log('     DLP >= CDH >= DDH');
console.log('');
console.log('  2. Parameter Requirements:');
console.log('     - Use safe primes: p = 2q + 1');
console.log('     - Minimum 2048 bits');
console.log('     - Work in prime-order subgroup');
console.log('');
console.log('  3. Implementation Requirements:');
console.log('     - Validate public values');
console.log('     - Check subgroup membership');
console.log('     - Use standard parameters (RFC 3526)');
console.log('');
console.log('  4. Algorithm Awareness:');
console.log('     - Baby-step giant-step: O(sqrt(p))');
console.log('     - Index calculus: Sub-exponential');
console.log('     - Pohlig-Hellman: Exploits smooth group order');

console.log('\n' + '='.repeat(60));
console.log('Diffie-Hellman Security Tutorial Complete');
console.log('='.repeat(60));
