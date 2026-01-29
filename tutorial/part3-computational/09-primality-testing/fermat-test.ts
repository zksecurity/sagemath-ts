/**
 * # Fermat Primality Test and Pseudoprimes
 *
 * Fermat's Little Theorem states: If p is prime and gcd(a, p) = 1, then
 *   a^(p-1) = 1 (mod p)
 *
 * The contrapositive gives us a primality test: If a^(n-1) != 1 (mod n)
 * for some a coprime to n, then n is definitely composite.
 *
 * ## The Problem: Fermat Liars and Pseudoprimes
 *
 * Unfortunately, the converse isn't true. Some composite numbers n satisfy
 * a^(n-1) = 1 (mod n) for certain bases a. These are called:
 * - "Fermat pseudoprimes" to base a
 * - a is called a "Fermat liar" for n
 *
 * Even worse, Carmichael numbers are composite but satisfy the Fermat
 * condition for ALL bases coprime to n!
 *
 * ## Complexity Analysis
 *
 * Time: O(log n) multiplications, each with O((log n)^2) bit operations
 * Overall: O((log n)^3) bit operations
 *
 * This is MUCH faster than trial division's O(sqrt(n))!
 */

import { power_mod, gcd, is_prime, factor, euler_phi } from 'sagemath-ts';

console.log("=== Fermat Primality Test ===\n");

// Fermat primality test implementation
function fermatTest(n: bigint, a: bigint): boolean {
  // Returns true if n is probably prime (or a pseudoprime to base a)
  // Returns false if n is definitely composite

  if (n <= 1n) return false;
  if (n === 2n) return true;
  if (n % 2n === 0n) return false;

  // Ensure a is in valid range
  a = ((a % n) + n) % n;
  if (a === 0n || a === 1n || a === n - 1n) {
    // Trivial cases - need a different witness
    return true;
  }

  // Check if gcd(a, n) > 1 - this would prove n composite
  const g = gcd(a, n);
  if (g > 1n) {
    return false; // Found a factor!
  }

  // Fermat's test: a^(n-1) should equal 1 mod n
  const result = power_mod(a, n - 1n, n);
  return result === 1n;
}

// Example 1: Fermat's Little Theorem in action
console.log("Example 1: Fermat's Little Theorem for primes");
console.log("-".repeat(40));

const p = 17n;
console.log(`  For prime p = ${p}, testing a^(p-1) = 1 (mod p):`);
for (let a = 2n; a <= 5n; a++) {
  const result = power_mod(a, p - 1n, p);
  console.log(`    ${a}^${p-1n} mod ${p} = ${result}`);
}

// Example 2: Detecting composites
console.log("\nExample 2: Fermat test detecting composites");
console.log("-".repeat(40));

const composite = 15n; // 3 * 5
console.log(`  For composite n = ${composite}:`);
for (let a = 2n; a <= 8n; a++) {
  if (gcd(a, composite) > 1n) {
    console.log(`    a=${a}: gcd(${a}, ${composite}) = ${gcd(a, composite)} (factor found!)`);
  } else {
    const result = power_mod(a, composite - 1n, composite);
    const passes = result === 1n;
    console.log(`    a=${a}: ${a}^${composite-1n} mod ${composite} = ${result} ${passes ? '(LIAR!)' : '(witness)'}`);
  }
}

// Example 3: Fermat pseudoprimes
console.log("\nExample 3: Fermat Pseudoprimes");
console.log("-".repeat(40));

// 341 = 11 * 31 is a pseudoprime to base 2
const n341 = 341n;
console.log(`  n = ${n341} = 11 * 31 (composite)`);
console.log(`  2^${n341-1n} mod ${n341} = ${power_mod(2n, n341 - 1n, n341)}`);
console.log(`  This means 341 is a pseudoprime to base 2!`);
console.log(`  But: 3^${n341-1n} mod ${n341} = ${power_mod(3n, n341 - 1n, n341)} (not 1, so composite)`);

// Example 4: Carmichael numbers - the worst case
console.log("\nExample 4: Carmichael Numbers");
console.log("-".repeat(40));

// 561 = 3 * 11 * 17 is the smallest Carmichael number
const carmichael = 561n;
console.log(`  ${carmichael} = 3 * 11 * 17 is a Carmichael number`);
console.log(`  Testing Fermat's condition for various bases:`);

let liars = 0;
let witnesses = 0;
for (let a = 2n; a <= 20n; a++) {
  const g = gcd(a, carmichael);
  if (g > 1n) {
    console.log(`    a=${a}: gcd = ${g} (trivial witness)`);
    witnesses++;
  } else {
    const result = power_mod(a, carmichael - 1n, carmichael);
    if (result === 1n) {
      liars++;
    } else {
      witnesses++;
    }
  }
}
console.log(`  For bases 2-20: ${liars} liars, ${witnesses} witnesses`);
console.log(`  All bases coprime to 561 are liars!`);

// Example 5: Why Carmichael numbers fool the Fermat test
console.log("\nExample 5: Understanding Carmichael Numbers");
console.log("-".repeat(40));

console.log(`
  Korselt's Criterion: n is a Carmichael number if and only if:
  1. n is squarefree
  2. For each prime p dividing n: (p - 1) | (n - 1)

  For 561 = 3 * 11 * 17:
  - 561 is squarefree (no repeated prime factors)
  - 561 - 1 = 560 = 16 * 35 = 2^4 * 5 * 7
  - (3 - 1) = 2 divides 560
  - (11 - 1) = 10 divides 560
  - (17 - 1) = 16 divides 560

  This is why a^560 = 1 (mod 561) for all a coprime to 561!
`);

// Verify the criterion
const n = 561n;
const factors = factor(n);
console.log(`  Factorization of ${n}: ${factors.map(([p, e]) => `${p}^${e}`).join(' * ')}`);
console.log(`  n - 1 = ${n - 1n}`);
for (const [p] of factors) {
  if (p > 0n) {
    const divides = (n - 1n) % (p - 1n) === 0n;
    console.log(`  (${p} - 1) = ${p - 1n} divides ${n - 1n}? ${divides}`);
  }
}

// Example 6: Improving the Fermat test
console.log("\nExample 6: Probabilistic Fermat Test");
console.log("-".repeat(40));

function probabilisticFermat(n: bigint, iterations: number = 10): boolean {
  if (n <= 1n) return false;
  if (n === 2n || n === 3n) return true;
  if (n % 2n === 0n) return false;

  for (let i = 0; i < iterations; i++) {
    // Random base between 2 and n-2
    const range = n - 4n;
    const a = 2n + BigInt(Math.floor(Math.random() * Number(range < BigInt(Number.MAX_SAFE_INTEGER) ? range : BigInt(Number.MAX_SAFE_INTEGER))));

    if (!fermatTest(n, a)) {
      return false; // Definitely composite
    }
  }
  return true; // Probably prime (unless it's a Carmichael number!)
}

const testCases = [
  { n: 97n, expected: true },
  { n: 100n, expected: false },
  { n: 561n, expected: false }, // Carmichael - might be fooled!
  { n: 1009n, expected: true },
];

console.log("  Testing with 10 random bases:");
for (const { n, expected } of testCases) {
  const result = probabilisticFermat(n, 10);
  const actual = is_prime(n);
  const correct = result === actual;
  console.log(`    n=${n}: Fermat says ${result ? 'prime' : 'composite'}, actually ${actual ? 'prime' : 'composite'} ${correct ? '' : '(WRONG!)'}`);
}

// Example 7: Small Carmichael numbers
console.log("\nExample 7: First Few Carmichael Numbers");
console.log("-".repeat(40));

// Known small Carmichael numbers
const carmichaelNumbers = [561n, 1105n, 1729n, 2465n, 2821n, 6601n, 8911n];

console.log("  Carmichael numbers fool the Fermat test:");
for (const cn of carmichaelNumbers.slice(0, 5)) {
  const factors = factor(cn);
  const factorStr = factors.filter(([p]) => p > 0n).map(([p]) => p.toString()).join(' * ');
  console.log(`    ${cn} = ${factorStr}`);
}

// Example 8: 1729 - the Hardy-Ramanujan number
console.log("\nExample 8: 1729 - A Famous Carmichael Number");
console.log("-".repeat(40));

console.log(`
  1729 = 7 * 13 * 19 is a Carmichael number.

  It's also famous as the Hardy-Ramanujan number:
  1729 = 1^3 + 12^3 = 9^3 + 10^3

  It's the smallest number expressible as the sum of two cubes
  in two different ways!

  Ramanujan mentioned this when Hardy visited him in the hospital,
  arriving in taxi number 1729.
`);

console.log("\n=== Summary ===");
console.log(`
  Fermat's Primality Test:
  - Based on Fermat's Little Theorem
  - Very fast: O((log n)^3) vs O(sqrt(n)) for trial division
  - Can be fooled by Fermat pseudoprimes
  - Carmichael numbers fool it for ALL coprime bases

  Key concepts:
  - Fermat liar: a base a where composite n passes the test
  - Fermat witness: a base a that proves n is composite
  - Carmichael number: composite where all coprime bases are liars

  The lesson: We need a stronger test. Enter Miller-Rabin!
`);
