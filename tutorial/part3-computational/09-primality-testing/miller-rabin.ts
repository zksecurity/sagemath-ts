/**
 * # Miller-Rabin Primality Test
 *
 * The Miller-Rabin test is a refinement of the Fermat test that is immune
 * to Carmichael numbers. It's based on a key observation about primes.
 *
 * ## Mathematical Foundation
 *
 * For any odd prime p, we can write p - 1 = 2^s * d where d is odd.
 *
 * By Fermat's Little Theorem, a^(p-1) = 1 (mod p) for gcd(a,p) = 1.
 *
 * The key insight: In a field, the only square roots of 1 are 1 and -1.
 *
 * So for the sequence: a^d, a^(2d), a^(4d), ..., a^(2^s * d) = a^(p-1)
 * - The final term must be 1
 * - Each term is the square of the previous
 * - Before we hit 1, we must see -1 (= p-1)
 *
 * If we ever see a square root of 1 that isn't +/- 1, n is composite!
 *
 * ## Complexity Analysis
 *
 * Time: O(k * (log n)^3) where k is the number of witnesses tested
 * For a single witness: same as Fermat test
 *
 * Error probability: At most (1/4)^k for k random witnesses
 * With 20 witnesses, error probability < 10^(-12)
 */

import { power_mod, gcd, is_prime, factor } from 'sagemath-ts';

console.log("=== Miller-Rabin Primality Test ===\n");

// Miller-Rabin test for a single witness
function millerRabinWitness(n: bigint, a: bigint): boolean {
  // Returns true if n passes the Miller-Rabin test for base a
  // Returns false if n is definitely composite (a is a witness)

  if (n < 2n) return false;
  if (n === 2n) return true;
  if (n % 2n === 0n) return false;

  // Write n - 1 = 2^s * d where d is odd
  let d = n - 1n;
  let s = 0n;
  while (d % 2n === 0n) {
    d /= 2n;
    s++;
  }

  // Reduce a mod n
  a = a % n;
  if (a === 0n) return true;

  // Check if gcd(a, n) > 1
  if (gcd(a, n) > 1n) return false;

  // Compute a^d mod n
  let x = power_mod(a, d, n);

  // If a^d = 1 or -1, the test passes
  if (x === 1n || x === n - 1n) return true;

  // Square repeatedly, looking for -1
  for (let i = 1n; i < s; i++) {
    x = (x * x) % n;

    // If we hit 1 without seeing -1 first, n is composite
    if (x === 1n) return false;

    // If we hit -1, the test passes
    if (x === n - 1n) return true;
  }

  // We reached a^(n-1) without ever seeing -1, so n is composite
  return false;
}

// Example 1: The s and d decomposition
console.log("Example 1: Writing n-1 = 2^s * d");
console.log("-".repeat(40));

function decompose(n: bigint): [bigint, bigint] {
  let d = n;
  let s = 0n;
  while (d % 2n === 0n) {
    d /= 2n;
    s++;
  }
  return [s, d];
}

const primes = [7n, 13n, 17n, 31n, 97n, 127n];
for (const p of primes) {
  const [s, d] = decompose(p - 1n);
  console.log(`  ${p} - 1 = ${p - 1n} = 2^${s} * ${d}`);
}

// Example 2: The Miller-Rabin sequence
console.log("\nExample 2: The Miller-Rabin sequence for a prime");
console.log("-".repeat(40));

const p = 97n;
const [s_p, d_p] = decompose(p - 1n);
const a = 2n;

console.log(`  Prime p = ${p}, a = ${a}`);
console.log(`  p - 1 = ${p - 1n} = 2^${s_p} * ${d_p}`);
console.log(`  Sequence of powers:`);

let x = power_mod(a, d_p, p);
console.log(`    a^d = ${a}^${d_p} = ${x} mod ${p}`);

for (let i = 1n; i <= s_p; i++) {
  x = (x * x) % p;
  console.log(`    a^(2^${i}*d) = ${x} mod ${p}${x === p - 1n ? ' = -1' : ''}${x === 1n ? ' = 1' : ''}`);
}

// Example 3: Detecting composites that fool Fermat
console.log("\nExample 3: Miller-Rabin vs Carmichael numbers");
console.log("-".repeat(40));

const carmichael = 561n;
const [s_c, d_c] = decompose(carmichael - 1n);

console.log(`  Carmichael number n = ${carmichael} = 3 * 11 * 17`);
console.log(`  n - 1 = ${carmichael - 1n} = 2^${s_c} * ${d_c}`);

for (const base of [2n, 3n, 5n, 7n]) {
  console.log(`\n  Testing base a = ${base}:`);

  let val = power_mod(base, d_c, carmichael);
  let foundWitness = false;
  let sequence = [`a^d = ${val}`];

  if (val === 1n || val === carmichael - 1n) {
    console.log(`    ${sequence[0]} (passes)`);
    continue;
  }

  for (let i = 1n; i < s_c; i++) {
    const prev = val;
    val = (val * val) % carmichael;
    sequence.push(`${val}`);

    if (val === 1n && prev !== carmichael - 1n) {
      console.log(`    Sequence: ${sequence.join(' -> ')}`);
      console.log(`    WITNESS! Found sqrt(1) = ${prev} != +/-1`);
      foundWitness = true;
      break;
    }
    if (val === carmichael - 1n) {
      console.log(`    Sequence: ${sequence.join(' -> ')}`);
      console.log(`    Found -1, test passes for this base`);
      break;
    }
  }

  if (!foundWitness && val !== carmichael - 1n) {
    val = (val * val) % carmichael;
    sequence.push(`${val}`);
    console.log(`    Sequence: ${sequence.join(' -> ')}`);
    if (val !== 1n) {
      console.log(`    WITNESS! Final value ${val} != 1`);
    } else {
      console.log(`    Never saw -1, so this is a witness`);
    }
  }
}

// Example 4: Why Miller-Rabin works
console.log("\nExample 4: Why Miller-Rabin Works");
console.log("-".repeat(40));

console.log(`
  Key insight: In Z/nZ where n is an odd prime:
  - The equation x^2 = 1 has exactly 2 solutions: x = 1 and x = -1
  - This is because Z/nZ is a field

  For composite n:
  - Z/nZ is NOT a field
  - x^2 = 1 may have MORE than 2 solutions

  Example: n = 15 = 3 * 5
  By CRT, Z/15Z is isomorphic to Z/3Z x Z/5Z.
  Solutions to x^2 = 1:
`);

const n15 = 15n;
console.log(`  Square roots of 1 mod ${n15}:`);
for (let x = 1n; x < n15; x++) {
  if ((x * x) % n15 === 1n) {
    console.log(`    ${x}^2 = ${(x * x) % n15} mod ${n15}   (x = ${x})`);
  }
}
console.log(`  Four solutions! In a field, there would only be two.`);

// Example 5: Witness counts
console.log("\nExample 5: Witnesses vs Liars for Miller-Rabin");
console.log("-".repeat(40));

function countWitnesses(n: bigint): [number, number] {
  let witnesses = 0;
  let liars = 0;
  for (let a = 2n; a < n - 1n; a++) {
    if (millerRabinWitness(n, a)) {
      liars++;
    } else {
      witnesses++;
    }
  }
  return [witnesses, liars];
}

const testComposites = [
  { n: 9n, label: "9 = 3^2" },
  { n: 15n, label: "15 = 3*5" },
  { n: 21n, label: "21 = 3*7" },
  { n: 25n, label: "25 = 5^2" },
  { n: 561n, label: "561 = 3*11*17 (Carmichael)" },
];

console.log("  For composite n, counting witnesses in [2, n-2]:");
for (const { n, label } of testComposites) {
  if (n <= 50n) {
    const [witnesses, liars] = countWitnesses(n);
    const total = Number(n) - 3;
    const witnessRatio = ((witnesses / total) * 100).toFixed(1);
    console.log(`    ${label}: ${witnesses}/${total} witnesses (${witnessRatio}%)`);
  }
}

// Example 6: Deterministic Miller-Rabin with known witnesses
console.log("\nExample 6: Deterministic Miller-Rabin");
console.log("-".repeat(40));

console.log(`
  For numbers below certain bounds, specific witness sets
  guarantee a correct answer:

  n < 2,047: test with {2}
  n < 1,373,653: test with {2, 3}
  n < 9,080,191: test with {31, 73}
  n < 3,215,031,751: test with {2, 3, 5, 7}
  n < 341,550,071,728,321: test with {2, 3, 5, 7, 11, 13, 17}

  For n < 3,317,044,064,679,887,385,961,981:
  test with {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}

  This is what sagemath-ts uses internally!
`);

// Example 7: Probabilistic analysis
console.log("Example 7: Error Probability");
console.log("-".repeat(40));

function millerRabinProbabilistic(n: bigint, k: number): boolean {
  if (n < 2n) return false;
  if (n === 2n || n === 3n) return true;
  if (n % 2n === 0n) return false;

  // Test with k random witnesses
  for (let i = 0; i < k; i++) {
    const range = n - 4n;
    const randomVal = BigInt(Math.floor(Math.random() * Number(range < BigInt(Number.MAX_SAFE_INTEGER) ? range : BigInt(Number.MAX_SAFE_INTEGER))));
    const a = 2n + randomVal;

    if (!millerRabinWitness(n, a)) {
      return false; // Definitely composite
    }
  }
  return true; // Probably prime
}

console.log(`
  For a random base a and composite odd n:
  - At least 3/4 of bases are witnesses (prove n is composite)
  - At most 1/4 of bases are liars (say n might be prime)

  With k independent random bases:
  - P(all liars) <= (1/4)^k
  - k=10: P(error) < 10^-6
  - k=20: P(error) < 10^-12
  - k=40: P(error) < 10^-24

  This is good enough for cryptographic applications!
`);

// Example 8: Comparison with sagemath-ts is_prime
console.log("Example 8: Verification against is_prime");
console.log("-".repeat(40));

const testNumbers = [
  4n, 9n, 15n, 17n, 25n, 31n, 91n, 97n,
  341n, // Pseudoprime to base 2
  561n, // Carmichael
  1009n,
  1729n, // Carmichael (taxicab number)
];

console.log("  Testing Miller-Rabin with witnesses {2, 3, 5, 7}:");
for (const n of testNumbers) {
  const witnesses = [2n, 3n, 5n, 7n];
  const mrResult = witnesses.every(a => millerRabinWitness(n, a));
  const actual = is_prime(n);
  const match = mrResult === actual;
  console.log(`    n=${n}: MR says ${mrResult ? 'prime' : 'composite'}, is_prime says ${actual ? 'prime' : 'composite'} ${match ? '' : '(MISMATCH!)'}`);
}

console.log("\n=== Summary ===");
console.log(`
  Miller-Rabin is the workhorse of primality testing:

  1. Based on Fermat + the structure of (Z/nZ)*
  2. Immune to Carmichael numbers
  3. Error probability: at most (1/4)^k with k witnesses
  4. Can be made deterministic with known witness sets
  5. O(k * (log n)^3) time complexity

  Key insight: The existence of "nontrivial square roots of 1"
  (square roots other than +/- 1) proves a number is composite.

  In practice:
  - Used in all major crypto libraries
  - 20-40 random witnesses give negligible error probability
  - Deterministic version used for numbers below certain bounds
`);
