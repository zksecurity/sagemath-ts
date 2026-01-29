/**
 * # Pollard's p-1 Factorization Algorithm
 *
 * Pollard's p-1 method exploits a weakness in primes p where p-1 is
 * "smooth" (has only small prime factors). It's based on Fermat's
 * Little Theorem.
 *
 * ## Mathematical Foundation
 *
 * Fermat's Little Theorem: If p is prime and gcd(a, p) = 1, then
 *   a^(p-1) = 1 (mod p)
 *
 * If p-1 | M (M is a multiple of p-1), then a^M = 1 (mod p).
 *
 * Key insight: If p | n but p-1 | M and q-1 does not divide M,
 * then gcd(a^M - 1, n) might be p (not n or 1).
 *
 * ## Algorithm
 *
 * 1. Choose a smoothness bound B
 * 2. Compute M = lcm(1, 2, ..., B) or M = product of prime powers <= B
 * 3. Compute a^M mod n for some base a
 * 4. Compute gcd(a^M - 1, n)
 * 5. If 1 < gcd < n, we found a factor!
 *
 * ## Complexity
 *
 * Time: O(B * log(B) * log(n)^2) to check B-smoothness
 *
 * Works when p-1 is B-smooth. Fails if p-1 has a large prime factor.
 */

import { gcd, power_mod, is_prime, factor, formatFactorization, prime_range, lcm } from 'sagemath-ts';

console.log("=== Pollard's p-1 Factorization Algorithm ===\n");

// Stage 1 of Pollard's p-1
function pollardPMinus1Stage1(n: bigint, B: bigint, base: bigint = 2n): bigint {
  // Compute a^M mod n where M = lcm(1, 2, ..., B)
  // More efficient: compute a^(p^k) for each prime power p^k <= B

  let a = base;
  const primes = prime_range(B + 1n);

  for (const p of primes) {
    // Find largest k such that p^k <= B
    let pk = p;
    while (pk <= B) {
      a = power_mod(a, p, n);
      pk *= p;
    }
  }

  return gcd(a - 1n, n);
}

// Example 1: Why p-1 smoothness matters
console.log("Example 1: p-1 Smoothness");
console.log("-".repeat(50));

console.log(`
  A number is B-smooth if all its prime factors are <= B.

  For Pollard's p-1 to work on n = p * q, we need p-1 or q-1
  to be smooth for some reasonable bound B.
`);

const primeExamples = [17n, 101n, 1009n, 10007n, 100003n];

console.log("  Prime p | Factors of p-1");
console.log("  " + "-".repeat(40));

for (const p of primeExamples) {
  const factors = factor(p - 1n);
  const factorStr = formatFactorization(factors);
  const smoothness = factors.reduce((max, [prime]) => prime > max ? prime : max, 1n);
  console.log(`  ${p.toString().padEnd(8)} | ${factorStr.padEnd(25)} (${smoothness}-smooth)`);
}

// Example 2: Basic p-1 example
console.log("\nExample 2: Basic Pollard p-1 Example");
console.log("-".repeat(50));

// Choose n where one factor has smooth p-1
const p1 = 101n;   // p1 - 1 = 100 = 2^2 * 5^2 (25-smooth)
const q1 = 10007n; // q1 - 1 = 10006 = 2 * 5003 (5003-smooth)
const n1 = p1 * q1;

console.log(`  n = ${n1} = ${p1} * ${q1}`);
console.log(`  ${p1} - 1 = ${p1 - 1n} = ${formatFactorization(factor(p1 - 1n))} (25-smooth)`);
console.log(`  ${q1} - 1 = ${q1 - 1n} = ${formatFactorization(factor(q1 - 1n))} (5003-smooth)`);

for (const B of [10n, 25n, 50n, 100n]) {
  const result = pollardPMinus1Stage1(n1, B);
  console.log(`  B = ${B}: gcd = ${result}${result !== 1n && result !== n1 ? ` (found ${p1}!)` : ""}`);
}

// Example 3: Step-by-step p-1
console.log("\nExample 3: Step-by-Step Pollard p-1");
console.log("-".repeat(50));

function pollardPMinus1Verbose(n: bigint, B: bigint): bigint {
  console.log(`  Factoring n = ${n} with bound B = ${B}`);
  console.log(`  Computing a^M mod n where M = lcm(1, ..., ${B})`);

  let a = 2n;
  const primes = prime_range(B + 1n);

  for (const p of primes) {
    let pk = p;
    while (pk <= B) {
      const oldA = a;
      a = power_mod(a, p, n);
      console.log(`    a = ${oldA}^${p} mod ${n} = ${a}`);
      pk *= p;
    }
  }

  const g = gcd(a - 1n, n);
  console.log(`  gcd(${a} - 1, ${n}) = gcd(${a - 1n}, ${n}) = ${g}`);

  return g;
}

const n_demo = 493n; // 17 * 29; 17-1 = 16 = 2^4 (2-smooth!)
console.log(`\n  n = ${n_demo} = 17 * 29`);
console.log(`  17 - 1 = 16 = 2^4 (very smooth!)`);
console.log(`  29 - 1 = 28 = 2^2 * 7`);
pollardPMinus1Verbose(n_demo, 5n);

// Example 4: Stage 2 (for when p-1 has one large prime factor)
console.log("\nExample 4: Stage 2 of Pollard p-1");
console.log("-".repeat(50));

console.log(`
  If p-1 = m * q where m is B1-smooth and q is a prime in (B1, B2],
  Stage 2 can find the factor.

  After Stage 1, we have a = 2^M mod n where M = lcm(1, ..., B1).
  If p-1 = m * q with m | M, then a = 2^M = 2^(km) for some k.

  Now 2^(m*q) = 1 (mod p), so 2^m is a q-th root of unity mod p.

  Stage 2: For each prime q in (B1, B2], check gcd(a^q - 1, n).
`);

function pollardPMinus1Stage2(
  n: bigint,
  B1: bigint,
  B2: bigint,
  base: bigint = 2n
): bigint {
  // Stage 1
  let a = base;
  const primes1 = prime_range(B1 + 1n);

  for (const p of primes1) {
    let pk = p;
    while (pk <= B1) {
      a = power_mod(a, p, n);
      pk *= p;
    }
  }

  let g = gcd(a - 1n, n);
  if (g > 1n && g < n) return g;

  // Stage 2: Check primes between B1 and B2
  const primes2 = prime_range(B1 + 1n, B2 + 1n);

  // Use baby-step giant-step style optimization
  // For simplicity, we'll do the naive approach here
  for (const q of primes2) {
    const aq = power_mod(a, q, n);
    g = gcd(aq - 1n, n);
    if (g > 1n && g < n) return g;
  }

  return g;
}

// Example where Stage 2 helps
const p2 = 1307n;   // p2 - 1 = 1306 = 2 * 653, where 653 is prime
const q2 = 1009n;   // q2 - 1 = 1008 = 2^4 * 3^2 * 7 (7-smooth)
const n2 = p2 * q2;

console.log(`\n  n = ${n2} = ${p2} * ${q2}`);
console.log(`  ${p2} - 1 = ${formatFactorization(factor(p2 - 1n))} (has large prime 653)`);
console.log(`  ${q2} - 1 = ${formatFactorization(factor(q2 - 1n))} (7-smooth)`);

console.log("\n  Stage 1 only (B1 = 100):");
let result = pollardPMinus1Stage1(n2, 100n);
console.log(`    Result: ${result}${result !== 1n && result !== n2 ? ` (success!)` : " (failed)"}`);

console.log("\n  With Stage 2 (B1 = 10, B2 = 700):");
result = pollardPMinus1Stage2(n2, 10n, 700n);
console.log(`    Result: ${result}${result !== 1n && result !== n2 ? ` (success! Found ${q2})` : " (failed)"}`);

// Example 5: When p-1 fails
console.log("\nExample 5: When Pollard p-1 Fails");
console.log("-".repeat(50));

console.log(`
  Pollard p-1 fails when BOTH p-1 and q-1 have large prime factors.

  This is why strong primes (p-1 has a large prime factor) were
  historically recommended for RSA.

  Modern RSA (2048+ bits) uses random primes because:
  1. Primes are large enough that p-1 methods are impractical anyway
  2. Finding strong primes is expensive
  3. The probability of a random prime being weak is negligible
`);

// Two primes with large factors in p-1
const p3 = 1000000007n; // p3 - 1 = 2 * 500000003
const q3 = 1000000009n; // q3 - 1 = 2^3 * 125000001 (has large factor)
const n3 = p3 * q3;

console.log(`  n = ${n3}`);
console.log(`  Both factors have p-1 with large prime divisors.`);

const start = performance.now();
result = pollardPMinus1Stage1(n3, 10000n);
const time = performance.now() - start;
console.log(`  Stage 1 with B = 10000: ${result} (${time.toFixed(2)}ms)`);

// Example 6: Comparison with Pollard's rho
console.log("\nExample 6: Pollard p-1 vs Pollard rho");
console.log("-".repeat(50));

console.log(`
  Pollard's p-1:
  - Works when p-1 is smooth
  - Time depends on smoothness of p-1, not size of p
  - Can be very fast or completely fail

  Pollard's rho:
  - Works for any composite (given enough time)
  - Time depends on size of smallest factor
  - More predictable: O(sqrt(p)) for factor p
`);

// Compare on numbers with smooth p-1
const smoothP = 16411n;  // p - 1 = 16410 = 2 * 3 * 5 * 547 (547-smooth)
const normalQ = 10007n;  // q - 1 = 10006 = 2 * 5003 (not very smooth)
const testN = smoothP * normalQ;

console.log(`\n  n = ${testN} = ${smoothP} * ${normalQ}`);
console.log(`  ${smoothP} - 1 = ${formatFactorization(factor(smoothP - 1n))}`);

// p-1 method
const startPM1 = performance.now();
const resultPM1 = pollardPMinus1Stage1(testN, 1000n);
const timePM1 = performance.now() - startPM1;

// Pollard rho (simple version for comparison)
function pollardRhoBasic(n: bigint): bigint {
  if (n % 2n === 0n) return 2n;
  let x = 2n, y = 2n, d = 1n;
  const f = (x: bigint) => (x * x + 1n) % n;
  while (d === 1n) {
    x = f(x);
    y = f(f(y));
    d = gcd(x > y ? x - y : y - x, n);
  }
  return d;
}

const startRho = performance.now();
const resultRho = pollardRhoBasic(testN);
const timeRho = performance.now() - startRho;

console.log(`\n  Pollard p-1 (B=1000): ${resultPM1} in ${timePM1.toFixed(2)}ms`);
console.log(`  Pollard rho:          ${resultRho} in ${timeRho.toFixed(2)}ms`);

// Example 7: Choosing the bound B
console.log("\nExample 7: Choosing the Bound B");
console.log("-".repeat(50));

console.log(`
  Trade-off in choosing B:
  - Larger B: More likely to find factors (if p-1 is smooth)
  - Larger B: More computation time

  Time complexity: O(B * log B * log^2 n)

  Practical bounds:
  - B1 = 10^6 to 10^8 for Stage 1
  - B2 = 10^9 to 10^12 for Stage 2

  If p-1 is B-smooth, we succeed with high probability.
`);

// Example 8: RSA and smooth p-1
console.log("Example 8: Security Implications");
console.log("-".repeat(50));

console.log(`
  Historical concern: RSA modulus n = p * q might be factorable
  if p-1 or q-1 is smooth.

  Defense 1: "Strong primes"
  - Choose p such that p-1 has a large prime factor r
  - Also: r-1 should have a large prime factor
  - Makes p-1 methods impractical

  Defense 2: Large enough primes (modern approach)
  - For 1024-bit primes, the probability of being
    10^8-smooth is negligible (~10^-30)
  - Random primes are safe in practice

  Related attack: Williams' p+1 method
  - Works when p+1 is smooth instead of p-1
  - Similar algorithm using Lucas sequences

  Bottom line: Use sufficiently large random primes (2048+ bits for RSA)
`);

// Example 9: Computing lcm(1, 2, ..., B) efficiently
console.log("Example 9: The LCM Computation");
console.log("-".repeat(50));

console.log(`
  We need M = lcm(1, 2, 3, ..., B).

  Naive: lcm of all integers up to B
  Better: Product of prime powers p^k <= B for each prime p

  For prime p, the largest power <= B is p^floor(log_p(B)).

  lcm(1, ..., 10) = lcm(2^3, 3^2, 5, 7) = 8 * 9 * 5 * 7 = 2520
`);

function computeLCM(B: bigint): bigint {
  let result = 1n;
  const primes = prime_range(B + 1n);

  for (const p of primes) {
    let pk = p;
    while (pk <= B) {
      pk *= p;
    }
    pk /= p; // Largest power <= B
    result = lcm(result, pk);
  }

  return result;
}

for (const B of [5n, 10n, 20n, 50n, 100n]) {
  const L = computeLCM(B);
  console.log(`  lcm(1, ..., ${B}) = ${L}`);
}

console.log("\n=== Summary ===");
console.log(`
  Pollard's p-1 Method:

  Key idea: If p-1 is B-smooth, then p | gcd(a^M - 1, n)
  where M = lcm(1, ..., B).

  Algorithm:
  1. Stage 1: Compute a^M mod n
  2. Check gcd(a^M - 1, n)
  3. Stage 2: Handle one large prime factor in p-1

  Complexity: O(B * log B * log^2 n)

  When it works:
  - p-1 is B-smooth for reasonable B
  - Very fast compared to general methods

  When it fails:
  - Both p-1 and q-1 have large prime factors
  - Need Pollard rho or more advanced methods

  Security implication:
  - Avoid primes with smooth p-1 (historical concern)
  - Modern: sufficiently large random primes are safe
`);
