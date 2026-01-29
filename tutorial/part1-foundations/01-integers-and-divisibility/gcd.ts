/**
 * # Greatest Common Divisor (GCD) and the Euclidean Algorithm
 *
 * The Greatest Common Divisor is arguably the most important function in
 * computational number theory. It's used everywhere in cryptography, from
 * computing modular inverses to determining if numbers are coprime.
 *
 * ## Mathematical Background
 *
 * The GCD of two integers a and b, denoted gcd(a, b), is the largest
 * positive integer that divides both a and b.
 *
 * Key properties:
 *   - gcd(a, b) = gcd(b, a)                    (commutative)
 *   - gcd(a, 0) = |a|                          (identity)
 *   - gcd(a, b) = gcd(a - b, b) if a > b       (subtraction property)
 *   - gcd(a, b) = gcd(a mod b, b)              (Euclidean property)
 *
 * ## The Euclidean Algorithm
 *
 * The Euclidean Algorithm (c. 300 BCE) is one of the oldest known algorithms.
 * It computes gcd(a, b) by repeatedly applying: gcd(a, b) = gcd(b, a mod b)
 * until one operand becomes zero.
 *
 * Time complexity: O(log(min(a, b)))
 *
 * ## Why This Matters for Cryptography
 *
 *   - RSA: Finding e coprime to phi(n)
 *   - Extended GCD: Computing modular inverses
 *   - Pollard's rho: Factoring based on GCD computations
 *   - Lattice cryptography: Basis reduction algorithms
 */

import { gcd, GCD, factor, is_prime } from 'sagemath-ts';

// =============================================================================
// Example 1: Basic GCD Computation
// =============================================================================
console.log("=== Example 1: Basic GCD Computation ===\n");

/**
 * The gcd() function computes the greatest common divisor.
 * Results are always non-negative.
 */
const pairs: [bigint, bigint][] = [
  [48n, 18n],
  [100n, 35n],
  [17n, 23n],      // Two primes - GCD is 1
  [0n, 15n],       // GCD with 0
  [-12n, 8n],      // Negative numbers
];

console.log("Computing GCD of various pairs:");
for (const [a, b] of pairs) {
  console.log(`  gcd(${a}, ${b}) = ${gcd(a, b)}`);
}

console.log();

// =============================================================================
// Example 2: Understanding the Euclidean Algorithm
// =============================================================================
console.log("=== Example 2: Tracing the Euclidean Algorithm ===\n");

/**
 * Let's trace through the Euclidean Algorithm step by step.
 * This helps understand why the algorithm works and is so efficient.
 */
function euclideanTrace(a: bigint, b: bigint): void {
  console.log(`Computing gcd(${a}, ${b}) step by step:`);

  // Ensure a >= b for cleaner output
  if (a < b) [a, b] = [b, a];

  let step = 1;
  while (b !== 0n) {
    const quotient = a / b;
    const remainder = a % b;
    console.log(`  Step ${step}: ${a} = ${quotient} * ${b} + ${remainder}`);
    a = b;
    b = remainder;
    step++;
  }
  console.log(`  Result: gcd = ${a}\n`);
}

euclideanTrace(252n, 198n);
euclideanTrace(1071n, 462n);

// =============================================================================
// Example 3: GCD and Prime Factorization
// =============================================================================
console.log("=== Example 3: GCD via Prime Factorization ===\n");

/**
 * The GCD can also be computed by finding common prime factors.
 * For each prime p: the power of p in gcd(a,b) is min(power in a, power in b)
 *
 * This method is educational but inefficient for large numbers.
 * The Euclidean algorithm is much faster in practice.
 */
const num1 = 360n;  // 2^3 * 3^2 * 5
const num2 = 150n;  // 2 * 3 * 5^2

console.log(`Factorization of ${num1}:`);
const factors1 = factor(num1);
for (const [p, e] of factors1) {
  console.log(`  ${p}^${e}`);
}

console.log(`\nFactorization of ${num2}:`);
const factors2 = factor(num2);
for (const [p, e] of factors2) {
  console.log(`  ${p}^${e}`);
}

console.log(`\ngcd(${num1}, ${num2}) = ${gcd(num1, num2)}`);
// GCD = 2^1 * 3^1 * 5^1 = 30

console.log();

// =============================================================================
// Example 4: GCD of Multiple Numbers
// =============================================================================
console.log("=== Example 4: GCD of Multiple Numbers ===\n");

/**
 * The gcd() function can take an array to compute the GCD of multiple numbers.
 * gcd(a, b, c) = gcd(gcd(a, b), c)
 */
const numbers = [120n, 90n, 150n, 180n];
console.log(`Numbers: ${numbers.join(', ')}`);
console.log(`GCD of all: ${gcd(numbers)}`);

// Verify step by step
let runningGcd = numbers[0]!;
console.log("\nStep-by-step verification:");
for (let i = 1; i < numbers.length; i++) {
  const prev = runningGcd;
  runningGcd = gcd(runningGcd, numbers[i]!);
  console.log(`  gcd(${prev}, ${numbers[i]}) = ${runningGcd}`);
}

console.log();

// =============================================================================
// Example 5: Coprimality Testing
// =============================================================================
console.log("=== Example 5: Coprimality (Relative Primeness) ===\n");

/**
 * Two numbers are coprime if their GCD is 1.
 * This is crucial in cryptography:
 *   - RSA: e must be coprime to phi(n)
 *   - Modular inverse exists iff gcd(a, n) = 1
 */
function areCoprime(a: bigint, b: bigint): boolean {
  return gcd(a, b) === 1n;
}

const testPairs: [bigint, bigint, string][] = [
  [15n, 28n, "consecutive composites"],
  [8n, 15n, "no common factors"],
  [12n, 18n, "share factor 6"],
  [7n, 11n, "two primes"],
  [100n, 21n, "100 and 21"],
];

console.log("Testing coprimality:");
for (const [a, b, desc] of testPairs) {
  const coprime = areCoprime(a, b);
  console.log(`  gcd(${a}, ${b}) = ${gcd(a, b)} -> ${coprime ? "coprime" : "NOT coprime"} (${desc})`);
}

console.log();

// =============================================================================
// Example 6: RSA-Style Key Selection
// =============================================================================
console.log("=== Example 6: RSA Key Selection Simulation ===\n");

/**
 * In RSA, we need to find e coprime to phi(n) = (p-1)(q-1).
 * Let's simulate finding such an e.
 */
const p = 61n;
const q = 53n;
const n = p * q;
const phiN = (p - 1n) * (q - 1n);

console.log(`Simulating RSA key selection:`);
console.log(`  p = ${p}, q = ${q}`);
console.log(`  n = p * q = ${n}`);
console.log(`  phi(n) = (p-1)(q-1) = ${phiN}`);

// Find a suitable public exponent e
// Traditionally, e = 65537 (0x10001) is used, but let's find options
console.log("\nFinding valid public exponents e (coprime to phi(n)):");
let count = 0;
for (let e = 3n; e < 100n && count < 10; e += 2n) {
  if (gcd(e, phiN) === 1n) {
    console.log(`  e = ${e} is valid (gcd(${e}, ${phiN}) = 1)`);
    count++;
  }
}

// Check the standard exponent
const standardE = 65537n;
console.log(`\nStandard exponent: e = ${standardE}`);
console.log(`  gcd(${standardE}, ${phiN}) = ${gcd(standardE, phiN)}`);

console.log();

// =============================================================================
// Example 7: Efficiency of the Euclidean Algorithm
// =============================================================================
console.log("=== Example 7: Efficiency Analysis ===\n");

/**
 * The Euclidean Algorithm is remarkably efficient.
 * Worst case: consecutive Fibonacci numbers (requires the most steps).
 */
// Generate Fibonacci numbers
const fibs: bigint[] = [1n, 1n];
for (let i = 2; i < 15; i++) {
  fibs.push(fibs[i-1]! + fibs[i-2]!);
}

console.log("GCD of consecutive Fibonacci numbers (worst case):");
for (let i = 5; i < 12; i++) {
  const a = fibs[i]!;
  const b = fibs[i-1]!;

  // Count steps
  let steps = 0;
  let x = a, y = b;
  while (y !== 0n) {
    [x, y] = [y, x % y];
    steps++;
  }

  console.log(`  gcd(F_${i}=${a}, F_${i-1}=${b}) = ${x} in ${steps} steps`);
}

// The number of steps is approximately log_phi(min(a,b)) where phi = golden ratio

console.log();

// =============================================================================
// Example 8: Binary GCD (Stein's Algorithm)
// =============================================================================
console.log("=== Example 8: Binary GCD Concept ===\n");

/**
 * The Binary GCD (Stein's Algorithm) uses only subtraction, halving,
 * and comparisons - no division. This can be faster on some hardware.
 *
 * Key observations:
 *   - gcd(2a, 2b) = 2 * gcd(a, b)
 *   - gcd(2a, b) = gcd(a, b) if b is odd
 *   - gcd(a, b) = gcd(|a-b|, min(a,b)) if both odd
 */
function binaryGcd(a: bigint, b: bigint): bigint {
  // Handle signs
  if (a < 0n) a = -a;
  if (b < 0n) b = -b;

  if (a === 0n) return b;
  if (b === 0n) return a;

  // Find common factors of 2
  let shift = 0n;
  while (((a | b) & 1n) === 0n) {
    a >>= 1n;
    b >>= 1n;
    shift++;
  }

  // Remove remaining factors of 2 from a
  while ((a & 1n) === 0n) {
    a >>= 1n;
  }

  while (b !== 0n) {
    // Remove factors of 2 from b
    while ((b & 1n) === 0n) {
      b >>= 1n;
    }

    // Ensure a <= b
    if (a > b) [a, b] = [b, a];

    b -= a;
  }

  return a << shift;
}

console.log("Comparing standard GCD with binary GCD:");
const testCases = [[48n, 18n], [1071n, 462n], [12345n, 67890n]];
for (const [a, b] of testCases) {
  const standard = gcd(a, b);
  const binary = binaryGcd(a, b);
  console.log(`  gcd(${a}, ${b}): standard=${standard}, binary=${binary}, match=${standard === binary}`);
}

console.log();

// =============================================================================
// Exercise: Implement GCD from Scratch
// =============================================================================
console.log("=== Exercise: Implement GCD ===\n");

/**
 * EXERCISE: Implement the Euclidean algorithm yourself.
 * Then verify it matches the library's gcd() function.
 */

// TODO: Implement this function
function myGcd(a: bigint, b: bigint): bigint {
  // Make inputs non-negative
  if (a < 0n) a = -a;
  if (b < 0n) b = -b;

  // Euclidean algorithm
  while (b !== 0n) {
    const temp = b;
    b = a % b;
    a = temp;
  }

  return a;
}

console.log("Testing your GCD implementation:");
const testValues: [bigint, bigint][] = [[252n, 105n], [1234567890n, 987654321n], [0n, 42n]];
for (const [a, b] of testValues) {
  const yours = myGcd(a, b);
  const library = gcd(a, b);
  console.log(`  gcd(${a}, ${b}): yours=${yours}, library=${library}, match=${yours === library}`);
}

console.log();

// =============================================================================
// Cryptographic Application: Detecting Weak RSA Keys
// =============================================================================
console.log("=== Cryptographic Application: Detecting Weak RSA Keys ===\n");

/**
 * If two RSA moduli share a prime factor, we can factor both instantly
 * using GCD! This has been used to break real-world RSA keys with poor
 * random number generation.
 */
const p1 = 104729n;  // A prime
const q1 = 104743n;  // Another prime
const q2 = 104759n;  // Yet another prime

const n1 = p1 * q1;
const n2 = p1 * q2;  // Oops! Shares p1 with n1

console.log(`RSA Modulus 1: n1 = ${n1}`);
console.log(`RSA Modulus 2: n2 = ${n2}`);

const commonFactor = gcd(n1, n2);
console.log(`\ngcd(n1, n2) = ${commonFactor}`);

if (commonFactor > 1n && commonFactor < n1) {
  console.log("\nSECURITY BREACH: Common factor found!");
  console.log(`  n1 = ${commonFactor} * ${n1 / commonFactor}`);
  console.log(`  n2 = ${commonFactor} * ${n2 / commonFactor}`);
  console.log("Both RSA keys are now completely broken!");
}

console.log("\nThis attack has broken millions of real RSA keys with weak RNGs!");
