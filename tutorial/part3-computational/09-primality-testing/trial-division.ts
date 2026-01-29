/**
 * # Trial Division: The Simplest Primality Test
 *
 * Trial division is the most straightforward method to test whether a number
 * is prime. The idea is simple: if n is composite, it must have a factor
 * less than or equal to sqrt(n).
 *
 * ## Algorithm Description
 *
 * To test if n is prime:
 * 1. Check if n <= 1 (not prime by definition)
 * 2. Check if n == 2 (prime)
 * 3. Check if n is even (not prime if n > 2)
 * 4. For each odd d from 3 to sqrt(n), check if d divides n
 *
 * If no divisor is found, n is prime.
 *
 * ## Complexity Analysis
 *
 * Time complexity: O(sqrt(n))
 * - We only need to check divisors up to sqrt(n)
 * - If n has k bits, this is O(2^(k/2)) operations
 *
 * This becomes impractical for large numbers (e.g., 512-bit primes used
 * in cryptography), which is why we need probabilistic tests.
 *
 * ## Why sqrt(n)?
 *
 * If n = a * b where a <= b, then a <= sqrt(n).
 * Proof: If a > sqrt(n) and b > sqrt(n), then a*b > n, contradiction.
 */

import { isqrt, is_prime } from 'sagemath-ts';
import { trial_division } from 'sagemath-ts/arith';

console.log("=== Trial Division Primality Test ===\n");

// Implementation of trial division from scratch for educational purposes
function isPrimeTrialDivision(n: bigint): boolean {
  // Handle edge cases
  if (n <= 1n) return false;
  if (n === 2n) return true;
  if (n % 2n === 0n) return false;
  if (n === 3n) return true;
  if (n % 3n === 0n) return false;

  // Check divisors of form 6k +/- 1 up to sqrt(n)
  // All primes > 3 are of this form
  const sqrtN = isqrt(n);
  let i = 5n;

  while (i <= sqrtN) {
    if (n % i === 0n) return false;        // Check 6k - 1
    if (n % (i + 2n) === 0n) return false;  // Check 6k + 1
    i += 6n;
  }

  return true;
}

// Example 1: Testing small primes
console.log("Example 1: Testing small numbers");
console.log("-".repeat(40));

const smallTests = [1n, 2n, 3n, 4n, 5n, 6n, 7n, 11n, 13n, 15n, 17n, 19n, 23n];
for (const n of smallTests) {
  const result = isPrimeTrialDivision(n);
  console.log(`  ${n} is ${result ? 'prime' : 'composite'}`);
}

// Example 2: Using sagemath-ts's is_prime function
console.log("\nExample 2: Using sagemath-ts is_prime");
console.log("-".repeat(40));

const testNumbers = [
  97n,              // Prime
  100n,             // 2^2 * 5^2
  101n,             // Prime
  561n,             // 3 * 11 * 17 (Carmichael number)
  1009n,            // Prime
  1024n,            // 2^10
];

for (const n of testNumbers) {
  console.log(`  is_prime(${n}) = ${is_prime(n)}`);
}

// Example 3: Finding the smallest factor
console.log("\nExample 3: Trial division to find smallest factor");
console.log("-".repeat(40));

const composites = [15n, 21n, 35n, 91n, 143n, 561n, 1001n];
for (const n of composites) {
  const smallestFactor = trial_division(n);
  console.log(`  trial_division(${n}) = ${smallestFactor}`);
}

// Example 4: Performance demonstration
console.log("\nExample 4: Performance of trial division");
console.log("-".repeat(40));

// For small primes, trial division is fast
const measureTime = (fn: () => boolean): [boolean, number] => {
  const start = performance.now();
  const result = fn();
  const duration = performance.now() - start;
  return [result, duration];
};

const sizesAndPrimes: [string, bigint][] = [
  ["6-digit", 999983n],           // Largest 6-digit prime
  ["9-digit", 999999937n],        // A 9-digit prime
  ["12-digit", 999999999989n],    // A 12-digit prime
];

for (const [label, p] of sizesAndPrimes) {
  const [result, time] = measureTime(() => isPrimeTrialDivision(p));
  console.log(`  ${label} prime: ${result} (${time.toFixed(3)}ms)`);
}

// Example 5: Why we need better algorithms
console.log("\nExample 5: Why trial division isn't enough for cryptography");
console.log("-".repeat(40));

console.log(`
  For RSA, we need primes with hundreds of digits.
  A 256-bit number has about 77 decimal digits.

  Trial division would need to check up to 2^128 divisors.
  At 10^9 checks/second, this would take:
    2^128 / 10^9 seconds = 10^29 seconds
    = 3 * 10^21 years

  This is why we need probabilistic primality tests!
`);

// Example 6: The 6k +/- 1 optimization explained
console.log("Example 6: The 6k +/- 1 optimization");
console.log("-".repeat(40));

console.log(`
  All integers can be written as 6k + r where r in {0, 1, 2, 3, 4, 5}.

  - 6k + 0 = 6k is divisible by 2 and 3
  - 6k + 1 is coprime to 6 (potential prime)
  - 6k + 2 = 2(3k + 1) is even
  - 6k + 3 = 3(2k + 1) is divisible by 3
  - 6k + 4 = 2(3k + 2) is even
  - 6k + 5 = 6(k+1) - 1 is coprime to 6 (potential prime)

  So we only need to check numbers of the form 6k + 1 and 6k + 5.
  This reduces work by a factor of 3!
`);

// Verify the optimization
console.log("  Verifying primes > 3 are 6k +/- 1:");
const primesUnder50 = [5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n, 41n, 43n, 47n];
for (const p of primesUnder50) {
  const mod6 = p % 6n;
  console.log(`    ${p} mod 6 = ${mod6} (${mod6 === 1n ? '6k+1' : '6k-1'})`);
}

console.log("\n=== Summary ===");
console.log(`
  Trial division is:
  - Simple to understand and implement
  - Deterministic (always correct)
  - O(sqrt(n)) time complexity
  - Practical for numbers up to about 10^12
  - Too slow for cryptographic-sized primes

  Key insight: If n is composite, it has a factor <= sqrt(n).

  Next: We'll look at Fermat's test, which trades certainty for speed.
`);
