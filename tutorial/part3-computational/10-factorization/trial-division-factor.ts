/**
 * # Factorization by Trial Division
 *
 * Trial division is the simplest factorization algorithm. It finds the
 * prime factors of n by systematically dividing by potential prime factors.
 *
 * ## Algorithm Description
 *
 * 1. While n is even, divide by 2
 * 2. For each odd number d from 3 to sqrt(n):
 *    - While d divides n, divide n by d
 * 3. If n > 1 at the end, n is a prime factor
 *
 * ## Complexity Analysis
 *
 * Time complexity: O(sqrt(n)) = O(2^(k/2)) for k-bit number
 *
 * This is fine for small numbers but becomes impractical for:
 * - 40-digit numbers: millions of years
 * - RSA-sized numbers: longer than the age of the universe
 *
 * ## When It's Useful
 *
 * - Factoring small numbers quickly
 * - Finding small factors of large numbers
 * - As a preprocessing step before other algorithms
 */

import {
  factor,
  formatFactorization,
  is_prime,
  isqrt,
  gcd,
  prime_range,
} from 'sagemath-ts';
import { trial_division } from 'sagemath-ts/arith';

console.log("=== Factorization by Trial Division ===\n");

// Custom trial division implementation
function trialDivisionFactor(n: bigint): Array<[bigint, bigint]> {
  const factors: Array<[bigint, bigint]> = [];

  if (n === 0n) throw new Error("Cannot factor 0");

  // Handle negative numbers
  if (n < 0n) {
    factors.push([-1n, 1n]);
    n = -n;
  }

  if (n === 1n) return factors;

  // Factor out 2s
  let exp2 = 0n;
  while (n % 2n === 0n) {
    exp2++;
    n /= 2n;
  }
  if (exp2 > 0n) factors.push([2n, exp2]);

  // Try odd divisors up to sqrt(n)
  let d = 3n;
  while (d * d <= n) {
    let exp = 0n;
    while (n % d === 0n) {
      exp++;
      n /= d;
    }
    if (exp > 0n) factors.push([d, exp]);
    d += 2n;
  }

  // If n > 1, it's a prime factor
  if (n > 1n) {
    factors.push([n, 1n]);
  }

  return factors;
}

// Example 1: Basic factorization
console.log("Example 1: Basic Factorization");
console.log("-".repeat(50));

const testNumbers = [12n, 60n, 100n, 360n, 1000n, 1001n, 1024n];

for (const n of testNumbers) {
  const factors = trialDivisionFactor(n);
  const factorStr = factors
    .map(([p, e]) => (e === 1n ? `${p}` : `${p}^${e}`))
    .join(" * ");
  console.log(`  ${n} = ${factorStr || "1"}`);
}

// Example 2: Using sagemath-ts factor function
console.log("\nExample 2: Using sagemath-ts factor()");
console.log("-".repeat(50));

const moreNumbers = [
  123456789n,
  987654321n,
  1000000007n, // A prime
  2n ** 20n - 1n, // Mersenne-like
];

for (const n of moreNumbers) {
  const f = factor(n);
  console.log(`  factor(${n})`);
  console.log(`    = ${formatFactorization(f)}`);
}

// Example 3: Smooth numbers
console.log("\nExample 3: Smooth Numbers");
console.log("-".repeat(50));

console.log(`
  A number is B-smooth if all its prime factors are <= B.

  Smooth numbers are important in:
  - Quadratic sieve and number field sieve
  - Pollard's p-1 algorithm
  - Index calculus for discrete log
`);

function smoothness(n: bigint): bigint {
  const factors = factor(n < 0n ? -n : n);
  let maxPrime = 1n;
  for (const [p] of factors) {
    if (p > maxPrime && p > 0n) maxPrime = p;
  }
  return maxPrime;
}

const smoothExamples = [
  { n: 60n, expected: "2^2 * 3 * 5" },
  { n: 100n, expected: "2^2 * 5^2" },
  { n: 2310n, expected: "2 * 3 * 5 * 7 * 11" },
  { n: 720720n, expected: "2^4 * 3^2 * 5 * 7 * 11 * 13" },
];

for (const { n } of smoothExamples) {
  const s = smoothness(n);
  console.log(`  ${n} is ${s}-smooth (factors: ${formatFactorization(factor(n))})`);
}

// Example 4: Finding smallest factor
console.log("\nExample 4: Finding Smallest Factor");
console.log("-".repeat(50));

console.log(`  trial_division(n) returns the smallest prime factor of n.`);

const largish = [91n, 143n, 221n, 323n, 1001n, 10001n, 100001n];

for (const n of largish) {
  const smallestFactor = trial_division(n);
  const isPr = is_prime(n);
  console.log(`  trial_division(${n}) = ${smallestFactor}${isPr ? " (n is prime)" : ""}`);
}

// Example 5: Wheel factorization optimization
console.log("\nExample 5: Wheel Factorization");
console.log("-".repeat(50));

console.log(`
  Basic optimization: Only check numbers coprime to small primes.

  Wheel for 2: Check only odds (skip every 2nd number)
  Wheel for 2,3: Check 6k +/- 1 (skip 4 out of 6)
  Wheel for 2,3,5: Check numbers coprime to 30 (skip 22 out of 30)

  This doesn't change the O(sqrt(n)) complexity but reduces
  the constant factor significantly.
`);

function wheelFactor(n: bigint): Array<[bigint, bigint]> {
  const factors: Array<[bigint, bigint]> = [];

  if (n <= 1n) return factors;

  // Check 2, 3, 5
  for (const p of [2n, 3n, 5n]) {
    let exp = 0n;
    while (n % p === 0n) {
      exp++;
      n /= p;
    }
    if (exp > 0n) factors.push([p, exp]);
  }

  // Wheel increments for numbers coprime to 30
  const wheel = [4n, 2n, 4n, 2n, 4n, 6n, 2n, 6n];
  let d = 7n;
  let wheelIdx = 0;

  while (d * d <= n) {
    let exp = 0n;
    while (n % d === 0n) {
      exp++;
      n /= d;
    }
    if (exp > 0n) factors.push([d, exp]);

    d += wheel[wheelIdx];
    wheelIdx = (wheelIdx + 1) % wheel.length;
  }

  if (n > 1n) factors.push([n, 1n]);

  return factors;
}

// Verify wheel factorization
console.log("  Verifying wheel factorization:");
for (const n of [12345678901234567n, 98765432109876543n]) {
  const factorsNormal = factor(n);
  const factorsWheel = wheelFactor(n);

  const str1 = factorsNormal.map(([p, e]) => `${p}^${e}`).join(" * ");
  const str2 = factorsWheel.map(([p, e]) => `${p}^${e}`).join(" * ");

  console.log(`    ${n}:`);
  console.log(`      factor():      ${str1}`);
  console.log(`      wheelFactor(): ${str2}`);
}

// Example 6: Performance comparison
console.log("\nExample 6: Performance of Trial Division");
console.log("-".repeat(50));

function measureFactor(n: bigint): number {
  const start = performance.now();
  factor(n);
  return performance.now() - start;
}

console.log("  Factoring products of two primes of similar size:");

// Create semiprimes (products of two primes)
const semiprimes: Array<{ n: bigint; label: string }> = [
  { n: 143n, label: "11 * 13 (3 digits)" },
  { n: 10403n, label: "101 * 103 (5 digits)" },
  { n: 1020407n, label: "1009 * 1013 (7 digits)" },
  { n: 102030509n, label: "10007 * 10009 (9 digits)" },
  { n: 10004600209n, label: "100003 * 100003 (11 digits)" },
];

for (const { n, label } of semiprimes) {
  const time = measureFactor(n);
  const f = factor(n);
  console.log(`    ${label}: ${time.toFixed(2)}ms`);
  console.log(`      = ${formatFactorization(f)}`);
}

// Example 7: Limits of trial division
console.log("\nExample 7: Limits of Trial Division");
console.log("-".repeat(50));

console.log(`
  Time to factor n by trial division ~ O(sqrt(n)) divisions

  For a semiprime p * q where p ~ q ~ sqrt(n):
    - We need to check up to p ~ sqrt(n) divisors
    - At 10^9 divisions/second:

  Digits | Time
  -------|------------------
  10     | milliseconds
  15     | seconds
  20     | hours
  25     | years
  30     | millennia
  40     | age of universe

  RSA-2048 has 617-digit factors. Trial division would
  take 10^290 years. We need better algorithms!
`);

// Example 8: GCD-based factorization
console.log("Example 8: GCD as a Factoring Tool");
console.log("-".repeat(50));

console.log(`
  If we somehow find a number m with 1 < gcd(m, n) < n,
  we've factored n!

  This is the basis of:
  - Pollard's rho algorithm
  - Pollard's p-1 algorithm
  - Fermat's factorization
  - Quadratic sieve
  - Number field sieve
`);

// Demo: If we know some relation
const n_demo = 15n;
const m_demo = 21n; // Happens to share factor 3 with 15
console.log(`  n = ${n_demo}`);
console.log(`  Suppose we find m = ${m_demo} somehow`);
console.log(`  gcd(${m_demo}, ${n_demo}) = ${gcd(m_demo, n_demo)}`);
console.log(`  We found a factor!`);

// Example 9: Preprocessing with small primes
console.log("\nExample 9: Preprocessing with Small Primes");
console.log("-".repeat(50));

function factorWithPreprocessing(n: bigint, bound: bigint = 10000n): {
  factors: Array<[bigint, bigint]>;
  cofactor: bigint;
} {
  const factors: Array<[bigint, bigint]> = [];
  const primes = prime_range(bound);

  for (const p of primes) {
    if (p * p > n) break;
    let exp = 0n;
    while (n % p === 0n) {
      exp++;
      n /= p;
    }
    if (exp > 0n) factors.push([p, exp]);
  }

  return { factors, cofactor: n };
}

const largeWithSmallFactors = 2n ** 10n * 3n ** 5n * 5n ** 3n * 1000000007n;
console.log(`  n = 2^10 * 3^5 * 5^3 * 1000000007`);
console.log(`    = ${largeWithSmallFactors}`);

const { factors: smallFactors, cofactor } = factorWithPreprocessing(largeWithSmallFactors);
console.log(`  Small factors found: ${smallFactors.map(([p, e]) => `${p}^${e}`).join(" * ")}`);
console.log(`  Remaining cofactor: ${cofactor}`);
console.log(`  Cofactor is prime: ${is_prime(cofactor)}`);

console.log("\n=== Summary ===");
console.log(`
  Trial Division:
  - Simplest factorization algorithm
  - Time: O(sqrt(n)) = O(2^(k/2)) for k-bit number
  - Practical for ~20-digit numbers

  Optimizations:
  - Wheel factorization (skip non-prime candidates)
  - Precompute small primes (sieve)
  - Preprocessing: remove small factors first

  Key insight: Trial division works but is fundamentally
  limited by the sqrt(n) barrier. For cryptographic-size
  numbers, we need algorithms that exploit special structure:
  - Pollard rho (random walk)
  - Pollard p-1 (smooth p-1)
  - Quadratic sieve (smooth relations)
  - Number field sieve (algebraic number theory)
`);
