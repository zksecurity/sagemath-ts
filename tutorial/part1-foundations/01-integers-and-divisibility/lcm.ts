/**
 * # Least Common Multiple (LCM)
 *
 * The Least Common Multiple complements the GCD. While GCD finds the largest
 * common divisor, LCM finds the smallest common multiple.
 *
 * ## Mathematical Background
 *
 * The LCM of two integers a and b, denoted lcm(a, b), is the smallest
 * positive integer that is divisible by both a and b.
 *
 * Key relationship with GCD:
 *   lcm(a, b) * gcd(a, b) = |a * b|
 *
 * Therefore: lcm(a, b) = |a * b| / gcd(a, b)
 *
 * Properties:
 *   - lcm(a, b) = lcm(b, a)                    (commutative)
 *   - lcm(a, 1) = |a|                          (identity)
 *   - lcm(a, 0) = 0                            (convention)
 *   - lcm(a, a) = |a|
 *   - lcm(a, b) >= max(|a|, |b|)
 *
 * ## Why This Matters for Cryptography
 *
 *   - Carmichael's lambda function uses LCM
 *   - Key scheduling in block ciphers
 *   - Order calculations in groups
 *   - Chinese Remainder Theorem applications
 */

import { gcd, lcm, LCM, factor, euler_phi } from 'sagemath-ts';

// =============================================================================
// Example 1: Basic LCM Computation
// =============================================================================
console.log("=== Example 1: Basic LCM Computation ===\n");

/**
 * The lcm() function computes the least common multiple.
 * Results are always non-negative.
 */
const pairs: [bigint, bigint][] = [
  [4n, 6n],        // lcm = 12
  [12n, 18n],      // lcm = 36
  [7n, 13n],       // Coprime: lcm = 91
  [0n, 5n],        // lcm with 0 is 0
  [-6n, 8n],       // Works with negatives
];

console.log("Computing LCM of various pairs:");
for (const [a, b] of pairs) {
  console.log(`  lcm(${a}, ${b}) = ${lcm(a, b)}`);
}

console.log();

// =============================================================================
// Example 2: LCM-GCD Relationship
// =============================================================================
console.log("=== Example 2: The Fundamental LCM-GCD Relationship ===\n");

/**
 * The key identity: lcm(a, b) * gcd(a, b) = |a * b|
 *
 * This is why we compute: lcm(a, b) = |a * b| / gcd(a, b)
 */
console.log("Verifying lcm(a,b) * gcd(a,b) = |a*b|:\n");

const testPairs: [bigint, bigint][] = [
  [12n, 18n],
  [35n, 49n],
  [100n, 75n],
  [13n, 17n],
];

for (const [a, b] of testPairs) {
  const g = gcd(a, b);
  const l = lcm(a, b);
  const product = a * b;
  const absProduct = product < 0n ? -product : product;

  console.log(`  a = ${a}, b = ${b}`);
  console.log(`    gcd(a, b) = ${g}`);
  console.log(`    lcm(a, b) = ${l}`);
  console.log(`    gcd * lcm = ${g * l}`);
  console.log(`    |a * b| = ${absProduct}`);
  console.log(`    Equal? ${g * l === absProduct}\n`);
}

// =============================================================================
// Example 3: LCM via Prime Factorization
// =============================================================================
console.log("=== Example 3: LCM via Prime Factorization ===\n");

/**
 * For each prime p: the power of p in lcm(a,b) is max(power in a, power in b)
 * Compare this to GCD where we take the min.
 */
const num1 = 360n;  // 2^3 * 3^2 * 5^1
const num2 = 450n;  // 2^1 * 3^2 * 5^2

console.log(`Computing lcm(${num1}, ${num2}) via factorization:\n`);

const factors1 = factor(num1);
const factors2 = factor(num2);

console.log(`${num1} = ${factors1.map(([p, e]) => `${p}^${e}`).join(' * ')}`);
console.log(`${num2} = ${factors2.map(([p, e]) => `${p}^${e}`).join(' * ')}`);

// For LCM, we take max power for each prime:
// 2: max(3, 1) = 3
// 3: max(2, 2) = 2
// 5: max(1, 2) = 2
// lcm = 2^3 * 3^2 * 5^2 = 8 * 9 * 25 = 1800

console.log(`\nLCM takes max powers: 2^3 * 3^2 * 5^2 = ${8n * 9n * 25n}`);
console.log(`Library result: ${lcm(num1, num2)}`);

console.log();

// =============================================================================
// Example 4: LCM of Multiple Numbers
// =============================================================================
console.log("=== Example 4: LCM of Multiple Numbers ===\n");

/**
 * The lcm() function can take an array to compute the LCM of multiple numbers.
 * lcm(a, b, c) = lcm(lcm(a, b), c)
 */
const numbers = [4n, 6n, 8n, 10n, 12n];
console.log(`Numbers: ${numbers.join(', ')}`);
console.log(`LCM of all: ${lcm(numbers)}`);

// Verify step by step
let runningLcm = numbers[0]!;
console.log("\nStep-by-step verification:");
for (let i = 1; i < numbers.length; i++) {
  const prev = runningLcm;
  runningLcm = lcm(runningLcm, numbers[i]!);
  console.log(`  lcm(${prev}, ${numbers[i]}) = ${runningLcm}`);
}

console.log();

// =============================================================================
// Example 5: LCM and Coprime Numbers
// =============================================================================
console.log("=== Example 5: LCM and Coprime Numbers ===\n");

/**
 * When gcd(a, b) = 1 (coprime), then lcm(a, b) = |a * b|
 * This is because there are no shared prime factors.
 */
console.log("For coprime numbers, lcm = product:\n");

const coprimePairs: [bigint, bigint][] = [
  [3n, 5n],
  [7n, 11n],
  [8n, 9n],   // 2^3 and 3^2 have no common factors
  [14n, 15n], // 2*7 and 3*5 have no common factors
];

for (const [a, b] of coprimePairs) {
  const g = gcd(a, b);
  const l = lcm(a, b);
  console.log(`  gcd(${a}, ${b}) = ${g}, lcm(${a}, ${b}) = ${l}, product = ${a * b}`);
}

console.log();

// =============================================================================
// Example 6: Carmichael's Lambda Function
// =============================================================================
console.log("=== Example 6: Carmichael's Lambda Function ===\n");

/**
 * Carmichael's lambda function lambda(n) is defined using LCM.
 * It's the smallest positive integer m such that a^m = 1 (mod n)
 * for all a coprime to n.
 *
 * For n = p1^e1 * p2^e2 * ... * pk^ek:
 *   lambda(n) = lcm(lambda(p1^e1), lambda(p2^e2), ..., lambda(pk^ek))
 *
 * Where:
 *   lambda(2) = 1
 *   lambda(4) = 2
 *   lambda(2^e) = 2^(e-2) for e >= 3
 *   lambda(p^e) = p^(e-1) * (p-1) = phi(p^e) for odd prime p
 */
function carmichaelLambdaPrimePower(p: bigint, e: bigint): bigint {
  if (p === 2n) {
    if (e === 1n) return 1n;
    if (e === 2n) return 2n;
    return 2n ** (e - 2n);
  }
  // For odd prime: lambda(p^e) = phi(p^e) = p^(e-1) * (p - 1)
  return (p ** (e - 1n)) * (p - 1n);
}

function carmichaelLambda(n: bigint): bigint {
  if (n <= 0n) throw new Error("n must be positive");
  if (n === 1n) return 1n;

  const factors = factor(n);
  const lambdaParts: bigint[] = [];

  for (const [p, e] of factors) {
    if (p === -1n) continue;
    lambdaParts.push(carmichaelLambdaPrimePower(p, e));
  }

  return lcm(lambdaParts);
}

console.log("Carmichael lambda function values:\n");

const testValues = [1n, 2n, 3n, 4n, 5n, 6n, 7n, 8n, 9n, 10n, 15n, 16n, 100n];
for (const n of testValues) {
  const phi = euler_phi(n);
  const lambda = carmichaelLambda(n);
  console.log(`  n = ${n.toString().padStart(3)}: phi(n) = ${phi.toString().padStart(3)}, lambda(n) = ${lambda.toString().padStart(3)}, ratio = ${phi/lambda}`);
}

console.log("\nNote: lambda(n) always divides phi(n), but can be much smaller!");

console.log();

// =============================================================================
// Example 7: LCM in Fraction Operations
// =============================================================================
console.log("=== Example 7: LCM for Finding Common Denominators ===\n");

/**
 * When adding fractions, we need a common denominator.
 * The LCM of denominators gives the smallest such denominator.
 */
interface Fraction {
  num: bigint;
  den: bigint;
}

function addFractions(a: Fraction, b: Fraction): Fraction {
  const commonDen = lcm(a.den, b.den);
  const newNumA = a.num * (commonDen / a.den);
  const newNumB = b.num * (commonDen / b.den);
  const sumNum = newNumA + newNumB;

  // Reduce the result
  const g = gcd(sumNum < 0n ? -sumNum : sumNum, commonDen);
  return { num: sumNum / g, den: commonDen / g };
}

const frac1: Fraction = { num: 5n, den: 12n };
const frac2: Fraction = { num: 7n, den: 18n };

console.log(`Adding fractions: ${frac1.num}/${frac1.den} + ${frac2.num}/${frac2.den}`);
console.log(`  lcm(${frac1.den}, ${frac2.den}) = ${lcm(frac1.den, frac2.den)}`);

const sum = addFractions(frac1, frac2);
console.log(`  Result: ${sum.num}/${sum.den}`);

// Verify
const verify = Number(frac1.num)/Number(frac1.den) + Number(frac2.num)/Number(frac2.den);
console.log(`  Verification: ${verify.toFixed(6)} = ${Number(sum.num)/Number(sum.den)}`);

console.log();

// =============================================================================
// Example 8: LCM for Periodic Events
// =============================================================================
console.log("=== Example 8: LCM for Synchronization ===\n");

/**
 * LCM is used to find when periodic events coincide.
 * This is important in scheduling, signal processing, and cryptography.
 */
console.log("Problem: Three clocks chime every 4, 6, and 10 hours respectively.");
console.log("When will they all chime together again?\n");

const periods = [4n, 6n, 10n];
const coincidence = lcm(periods);

console.log(`Periods: ${periods.join(', ')} hours`);
console.log(`They coincide every ${coincidence} hours (= ${Number(coincidence)/24} days)`);

console.log("\nBreaking down why:");
for (const period of periods) {
  console.log(`  ${coincidence} / ${period} = ${coincidence / period} (whole number)`);
}

console.log();

// =============================================================================
// Example 9: Order of Elements in Groups
// =============================================================================
console.log("=== Example 9: LCM and Group Orders ===\n");

/**
 * In a group, if elements a and b have orders m and n respectively,
 * and they commute, then the order of ab divides lcm(m, n).
 *
 * For the group (Z/nZ)*, the exponent (smallest e where a^e = 1 for all a)
 * is exactly the Carmichael function lambda(n).
 */
console.log("In the group (Z/15Z)*, finding element orders:\n");

function orderInZnStar(a: bigint, n: bigint): bigint {
  if (gcd(a, n) !== 1n) throw new Error("a must be coprime to n");

  let power = a % n;
  let order = 1n;

  while (power !== 1n) {
    power = (power * a) % n;
    order++;
  }

  return order;
}

const modN = 15n;
const unitsOf15: bigint[] = [];
for (let i = 1n; i < modN; i++) {
  if (gcd(i, modN) === 1n) {
    unitsOf15.push(i);
  }
}

console.log(`Units of Z/${modN}Z: ${unitsOf15.join(', ')}\n`);

const orders: bigint[] = [];
for (const a of unitsOf15) {
  const ord = orderInZnStar(a, modN);
  orders.push(ord);
  console.log(`  ord(${a}) = ${ord}`);
}

console.log(`\nLCM of all orders: ${lcm(orders)}`);
console.log(`Carmichael lambda(${modN}): ${carmichaelLambda(modN)}`);
console.log("These are equal - the group exponent!");

console.log();

// =============================================================================
// Exercise: Implement LCM from scratch
// =============================================================================
console.log("=== Exercise: Implement LCM ===\n");

/**
 * EXERCISE: Implement LCM using the GCD-LCM relationship.
 * Remember: lcm(a, b) = |a * b| / gcd(a, b)
 */

// TODO: Implement this function
function myLcm(a: bigint, b: bigint): bigint {
  // Handle special cases
  if (a === 0n || b === 0n) return 0n;

  // Work with absolute values
  const absA = a < 0n ? -a : a;
  const absB = b < 0n ? -b : b;

  // Use the identity: lcm * gcd = product
  return (absA / gcd(absA, absB)) * absB;
}

console.log("Testing your LCM implementation:");
const testLcm: [bigint, bigint][] = [[12n, 8n], [35n, 49n], [0n, 5n], [-6n, 9n]];
for (const [a, b] of testLcm) {
  const yours = myLcm(a, b);
  const library = lcm(a, b);
  console.log(`  lcm(${a}, ${b}): yours=${yours}, library=${library}, match=${yours === library}`);
}

console.log();

// =============================================================================
// Cryptographic Application: Key Scheduling
// =============================================================================
console.log("=== Cryptographic Application: Key Scheduling ===\n");

/**
 * Some key schedules use LCM-related properties.
 * Here's a simple example: finding when two LFSR sequences align.
 */
console.log("LFSR Period Analysis:\n");

// Two LFSRs with different periods
const lfsr1Period = 31n;   // 2^5 - 1 (maximal length 5-bit LFSR)
const lfsr2Period = 63n;   // 2^6 - 1 (maximal length 6-bit LFSR)

console.log(`LFSR 1 period: ${lfsr1Period}`);
console.log(`LFSR 2 period: ${lfsr2Period}`);

const combinedPeriod = lcm(lfsr1Period, lfsr2Period);
console.log(`Combined period (via XOR): ${combinedPeriod}`);

// For coprime periods, combined period = product
console.log(`\nAre they coprime? gcd = ${gcd(lfsr1Period, lfsr2Period)}`);
console.log(`Product of periods: ${lfsr1Period * lfsr2Period}`);

console.log("\nIn stream ciphers, combining LFSRs with coprime periods");
console.log("gives a combined sequence with maximum period = product of periods.");
