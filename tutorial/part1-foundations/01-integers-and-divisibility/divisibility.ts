/**
 * # Divisibility: The Foundation of Number Theory
 *
 * Divisibility is one of the most fundamental concepts in number theory and
 * forms the basis for understanding modular arithmetic, which is essential
 * for modern cryptography.
 *
 * ## Mathematical Background
 *
 * We say that an integer `a` divides an integer `b` (written as `a | b`) if
 * there exists an integer `k` such that `b = a * k`.
 *
 * For example:
 *   - 3 | 12 because 12 = 3 * 4
 *   - 5 | 0 because 0 = 5 * 0
 *   - 1 | n for any integer n
 *
 * Key properties of divisibility:
 *   1. Reflexivity: a | a (every number divides itself)
 *   2. Transitivity: if a | b and b | c, then a | c
 *   3. Linear combinations: if a | b and a | c, then a | (bx + cy) for any x, y
 *
 * ## Why This Matters for Cryptography
 *
 * Divisibility relationships underpin:
 *   - RSA key generation (finding large primes)
 *   - Modular arithmetic operations
 *   - Group structure in elliptic curve cryptography
 *   - Hash functions and message authentication
 */

import { ZZ, Integer, gcd, factor, divisors, is_prime } from 'sagemath-ts';

// =============================================================================
// Example 1: Basic Divisibility Tests
// =============================================================================
console.log("=== Example 1: Basic Divisibility Tests ===\n");

/**
 * The Integer class provides a `divides` method that checks if the integer
 * divides another number. Note: a.divides(b) returns true if a | b.
 */
const a = new Integer(6);
const b = new Integer(18);
const c = new Integer(20);

console.log(`Does ${a} divide ${b}? ${a.divides(b.value)}`);  // true: 18 = 6 * 3
console.log(`Does ${a} divide ${c}? ${a.divides(c.value)}`);  // false: 20/6 is not an integer

// Any number divides 0
const zero = new Integer(0);
console.log(`Does ${a} divide 0? ${a.divides(zero.value)}`);  // true

// 1 divides everything
const one = new Integer(1);
console.log(`Does 1 divide ${c}? ${one.divides(c.value)}`);   // true

console.log();

// =============================================================================
// Example 2: Finding All Divisors
// =============================================================================
console.log("=== Example 2: Finding All Divisors ===\n");

/**
 * The divisors of a number are all positive integers that divide it evenly.
 * Understanding divisors is crucial for:
 *   - Factoring large numbers (the basis of RSA security)
 *   - Computing Euler's totient function
 *   - Analyzing group structures
 */
const n = 60n;
const divs = divisors(n);

console.log(`The divisors of ${n} are: ${divs.join(', ')}`);
console.log(`Number of divisors: ${divs.length}`);

// The divisors come from prime factorization
const factorization = factor(n);
console.log(`\nPrime factorization of ${n}:`);
for (const [prime, exp] of factorization) {
  console.log(`  ${prime}^${exp}`);
}

// Formula: if n = p1^a1 * p2^a2 * ... * pk^ak
// then number of divisors = (a1+1)(a2+1)...(ak+1)
// For 60 = 2^2 * 3^1 * 5^1, divisors = (2+1)(1+1)(1+1) = 3*2*2 = 12

console.log();

// =============================================================================
// Example 3: Perfect Numbers and Divisor Sums
// =============================================================================
console.log("=== Example 3: Perfect Numbers and Divisor Sums ===\n");

/**
 * A perfect number equals the sum of its proper divisors (divisors excluding itself).
 * The first few perfect numbers are: 6, 28, 496, 8128
 *
 * Perfect numbers have fascinating connections to Mersenne primes:
 * If 2^p - 1 is prime (a Mersenne prime), then 2^(p-1) * (2^p - 1) is perfect.
 */
function sumOfProperDivisors(num: bigint): bigint {
  const divs = divisors(num);
  // Sum all divisors except the number itself
  return divs.slice(0, -1).reduce((sum, d) => sum + d, 0n);
}

function isPerfect(num: bigint): boolean {
  return sumOfProperDivisors(num) === num;
}

const candidates = [6n, 28n, 12n, 496n, 100n];
for (const num of candidates) {
  const properSum = sumOfProperDivisors(num);
  console.log(`${num}: sum of proper divisors = ${properSum}, perfect? ${isPerfect(num)}`);
}

console.log();

// =============================================================================
// Example 4: Common Divisors and GCD Preview
// =============================================================================
console.log("=== Example 4: Common Divisors ===\n");

/**
 * Common divisors of two numbers are integers that divide both.
 * The greatest common divisor (GCD) is the largest such divisor.
 *
 * This concept is fundamental to:
 *   - Reducing fractions to lowest terms
 *   - Computing modular inverses
 *   - Key generation in cryptographic protocols
 */
const x = 48n;
const y = 36n;

const divsX = divisors(x);
const divsY = divisors(y);
const commonDivs = divsX.filter(d => divsY.includes(d));

console.log(`Divisors of ${x}: ${divsX.join(', ')}`);
console.log(`Divisors of ${y}: ${divsY.join(', ')}`);
console.log(`Common divisors: ${commonDivs.join(', ')}`);
console.log(`Greatest common divisor: ${gcd(x, y)}`);

console.log();

// =============================================================================
// Example 5: Prime Factorization
// =============================================================================
console.log("=== Example 5: Prime Factorization ===\n");

/**
 * The Fundamental Theorem of Arithmetic states that every integer > 1
 * can be uniquely expressed as a product of primes.
 *
 * This is the mathematical foundation for:
 *   - RSA encryption (difficulty of factoring large semiprimes)
 *   - Pollard's rho algorithm
 *   - Index calculus methods
 */
const numbers = [84n, 100n, 997n, 2310n];

for (const num of numbers) {
  const factors = factor(num);
  const primeCheck = is_prime(num);

  if (primeCheck) {
    console.log(`${num} is prime`);
  } else {
    const formatted = factors
      .map(([p, e]) => e === 1n ? `${p}` : `${p}^${e}`)
      .join(' * ');
    console.log(`${num} = ${formatted}`);
  }
}

console.log();

// =============================================================================
// Example 6: Divisibility Rules
// =============================================================================
console.log("=== Example 6: Divisibility Rules ===\n");

/**
 * There are efficient divisibility rules for small primes that don't require
 * full division. These are useful for quick checks.
 */
const testNumber = new Integer(12345678901234567890n);

console.log(`Testing divisibility of ${testNumber}:`);

// Divisible by 2 iff last digit is even (ends in 0, 2, 4, 6, 8)
const divBy2 = testNumber.value % 2n === 0n;
console.log(`  By 2? ${divBy2} (last digit: ${testNumber.value % 10n})`);

// Divisible by 3 iff digit sum is divisible by 3
const digitSum = testNumber.digit_sum();
console.log(`  By 3? ${digitSum % 3n === 0n} (digit sum: ${digitSum})`);

// Divisible by 5 iff last digit is 0 or 5
const lastDigit = testNumber.value % 10n;
console.log(`  By 5? ${lastDigit === 0n || lastDigit === 5n} (last digit: ${lastDigit})`);

// Divisible by 9 iff digit sum is divisible by 9
console.log(`  By 9? ${digitSum % 9n === 0n} (digit sum: ${digitSum})`);

console.log();

// =============================================================================
// Exercise: Implement your own divisibility check
// =============================================================================
console.log("=== Exercise: DIY Divisibility ===\n");

/**
 * EXERCISE: Without using the divides() method or the % operator,
 * implement a function that checks if a divides b using only
 * addition and subtraction.
 *
 * Hint: Keep subtracting a from b until you reach 0 or go negative.
 */

// TODO: Implement this function
function dividesBySubtraction(a: bigint, b: bigint): boolean {
  // Handle edge cases
  if (a === 0n) return b === 0n;
  if (b === 0n) return true;  // Everything divides 0

  // Work with absolute values
  const absA = a < 0n ? -a : a;
  let absB = b < 0n ? -b : b;

  // Keep subtracting until we hit 0 or go negative
  while (absB > 0n) {
    absB -= absA;
  }

  return absB === 0n;
}

// Test your implementation
console.log("Testing dividesBySubtraction:");
console.log(`  12 divides 48? ${dividesBySubtraction(12n, 48n)}`);   // true
console.log(`  12 divides 50? ${dividesBySubtraction(12n, 50n)}`);   // false
console.log(`  7 divides 49? ${dividesBySubtraction(7n, 49n)}`);     // true

console.log();

// =============================================================================
// Cryptographic Application: Finding Coprime Numbers
// =============================================================================
console.log("=== Cryptographic Application: Finding Coprimes ===\n");

/**
 * Two numbers are coprime (relatively prime) if their GCD is 1.
 * This is essential in cryptography:
 *   - RSA requires e and phi(n) to be coprime
 *   - Modular inverses only exist for coprime numbers
 */
const modulus = 15n;
console.log(`Finding numbers coprime to ${modulus}:`);

const coprimes: bigint[] = [];
for (let i = 1n; i < modulus; i++) {
  if (gcd(i, modulus) === 1n) {
    coprimes.push(i);
  }
}

console.log(`Numbers coprime to ${modulus}: ${coprimes.join(', ')}`);
console.log(`Count: ${coprimes.length}`);

// This count is exactly Euler's totient function phi(15)!
// phi(15) = phi(3) * phi(5) = 2 * 4 = 8
