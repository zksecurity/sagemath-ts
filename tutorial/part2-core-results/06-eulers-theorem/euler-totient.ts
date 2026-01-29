/**
 * # Euler's Totient Function
 *
 * Euler's totient function phi(n) counts the number of integers from 1 to n
 * that are coprime to n (i.e., gcd(k, n) = 1).
 *
 * ## Definition
 *
 *   phi(n) = |{k : 1 <= k <= n and gcd(k, n) = 1}|
 *
 * Equivalently, phi(n) is the order of the multiplicative group (Z/nZ)*.
 *
 * ## Key Formula
 *
 * For n = p1^e1 * p2^e2 * ... * pk^ek (prime factorization):
 *
 *   phi(n) = n * (1 - 1/p1) * (1 - 1/p2) * ... * (1 - 1/pk)
 *
 * Or equivalently:
 *
 *   phi(n) = product over i of p_i^(e_i - 1) * (p_i - 1)
 *
 * @module tutorial/part2/euler-totient
 */

import {
  euler_phi,
  gcd,
  factor,
  is_prime,
  Zmod,
  formatFactorization,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("EULER'S TOTIENT FUNCTION phi(n)");
console.log("=".repeat(70));

// ============================================================================
// Example 1: Definition and basic examples
// ============================================================================
console.log("\n=== Example 1: Definition of phi(n) ===\n");

console.log("phi(n) = number of integers k with 1 <= k <= n and gcd(k,n) = 1\n");

function computePhiDirectly(n: bigint): { phi: bigint; coprimes: bigint[] } {
  const coprimes: bigint[] = [];
  for (let k = 1n; k <= n; k++) {
    if (gcd(k, n) === 1n) {
      coprimes.push(k);
    }
  }
  return { phi: BigInt(coprimes.length), coprimes };
}

const examples = [1n, 2n, 6n, 8n, 10n, 12n];

for (const n of examples) {
  const { phi, coprimes } = computePhiDirectly(n);
  const coprimesStr = coprimes.length <= 12
    ? coprimes.join(", ")
    : coprimes.slice(0, 10).join(", ") + "...";
  console.log(`phi(${n.toString().padStart(2)}) = ${phi.toString().padStart(2)}  ` +
              `Coprimes: {${coprimesStr}}`);
}

// ============================================================================
// Example 2: Special cases
// ============================================================================
console.log("\n=== Example 2: Special Cases ===\n");

console.log("1. For prime p:");
console.log("   phi(p) = p - 1");
console.log("   (All numbers from 1 to p-1 are coprime to p)\n");

const primes = [2n, 3n, 5n, 7n, 11n, 13n];
for (const p of primes) {
  console.log(`   phi(${p.toString().padStart(2)}) = ${euler_phi(p).toString().padStart(2)} = ${p} - 1`);
}

console.log("\n2. For prime power p^k:");
console.log("   phi(p^k) = p^(k-1) * (p - 1) = p^k - p^(k-1)\n");

const primePowers = [
  { p: 2n, k: 3n },
  { p: 3n, k: 2n },
  { p: 5n, k: 2n },
  { p: 2n, k: 4n },
];

for (const { p, k } of primePowers) {
  const pk = p ** k;
  const phi = euler_phi(pk);
  const formula = `${p}^${k - 1n} * (${p} - 1) = ${p ** (k - 1n)} * ${p - 1n}`;
  console.log(`   phi(${p}^${k}) = phi(${pk.toString().padStart(2)}) = ${phi.toString().padStart(2)} = ${formula}`);
}

console.log("\n3. For n = 1:");
console.log(`   phi(1) = ${euler_phi(1n)} (by convention, gcd(1,1) = 1)`);

// ============================================================================
// Example 3: The multiplicative property
// ============================================================================
console.log("\n=== Example 3: Multiplicative Property ===\n");

console.log("If gcd(m, n) = 1, then phi(m*n) = phi(m) * phi(n)");
console.log("(phi is a multiplicative function)\n");

const pairs = [
  { m: 3n, n: 4n },
  { m: 5n, n: 6n },
  { m: 7n, n: 9n },
  { m: 8n, n: 15n },
];

for (const { m, n } of pairs) {
  const mn = m * n;
  const phiM = euler_phi(m);
  const phiN = euler_phi(n);
  const phiMN = euler_phi(mn);
  const product = phiM * phiN;
  const check = phiMN === product ? "OK" : "FAIL";
  console.log(`  phi(${m} * ${n}) = phi(${mn.toString().padStart(2)}) = ${phiMN.toString().padStart(2)}  ` +
              `phi(${m}) * phi(${n}) = ${phiM} * ${phiN} = ${product.toString().padStart(2)} [${check}]`);
}

console.log("\nCounterexample when gcd(m,n) > 1:");
const badM = 6n;
const badN = 4n;
console.log(`  gcd(${badM}, ${badN}) = ${gcd(badM, badN)} (not coprime)`);
console.log(`  phi(${badM} * ${badN}) = phi(${badM * badN}) = ${euler_phi(badM * badN)}`);
console.log(`  phi(${badM}) * phi(${badN}) = ${euler_phi(badM)} * ${euler_phi(badN)} = ${euler_phi(badM) * euler_phi(badN)}`);
console.log(`  These are NOT equal!`);

// ============================================================================
// Example 4: Computing phi from prime factorization
// ============================================================================
console.log("\n=== Example 4: Computing phi from Prime Factorization ===\n");

console.log("Formula: phi(n) = product of (p^(e-1) * (p - 1)) for each prime power p^e in n\n");

function computePhiFromFactors(n: bigint): { phi: bigint; explanation: string } {
  const factors = factor(n);
  let phi = 1n;
  const terms: string[] = [];

  for (const [p, e] of factors) {
    if (p === -1n) continue;
    const contribution = p ** (e - 1n) * (p - 1n);
    phi *= contribution;
    if (e === 1n) {
      terms.push(`(${p} - 1)`);
    } else {
      terms.push(`${p}^${e - 1n} * (${p} - 1)`);
    }
  }

  return { phi, explanation: terms.join(" * ") };
}

const computeExamples = [12n, 30n, 36n, 100n, 360n];

for (const n of computeExamples) {
  const factors = factor(n);
  const factorStr = formatFactorization(factors);
  const { phi, explanation } = computePhiFromFactors(n);
  const verify = euler_phi(n);

  console.log(`n = ${n.toString().padStart(3)} = ${factorStr.padEnd(15)}`);
  console.log(`    phi(${n}) = ${explanation} = ${phi}`);
  console.log(`    Verification: ${verify} [${phi === verify ? "OK" : "FAIL"}]\n`);
}

// ============================================================================
// Example 5: Relationship with the multiplicative group
// ============================================================================
console.log("\n=== Example 5: phi(n) as Order of (Z/nZ)* ===\n");

console.log("phi(n) equals the number of units (invertible elements) in Z/nZ");
console.log("This is the order of the multiplicative group (Z/nZ)*\n");

const n5 = 12n;
const Zn = Zmod(n5);
const units = Zn.units();

console.log(`Z/${n5}Z:`);
console.log(`  Elements: {0, 1, 2, ..., ${n5 - 1n}}`);
console.log(`  Units (invertible): {${units.map(u => u.value).join(", ")}}`);
console.log(`  |Units| = ${units.length}`);
console.log(`  phi(${n5}) = ${euler_phi(n5)}`);
console.log(`  Match: ${BigInt(units.length) === euler_phi(n5)}`);

console.log("\nAn element a is a unit iff gcd(a, n) = 1:");
console.log("-".repeat(50));
for (let a = 0n; a < n5; a++) {
  const g = gcd(a, n5);
  const isUnit = g === 1n;
  console.log(`  a = ${a.toString().padStart(2)}: gcd(${a}, ${n5}) = ${g.toString().padStart(2)}, unit: ${isUnit}`);
}

// ============================================================================
// Example 6: Gauss's divisor sum formula
// ============================================================================
console.log("\n=== Example 6: Gauss's Divisor Sum Formula ===\n");

console.log("Theorem: sum of phi(d) over all divisors d of n equals n");
console.log("         sum_{d | n} phi(d) = n\n");

function divisors(n: bigint): bigint[] {
  const result: bigint[] = [];
  for (let d = 1n; d * d <= n; d++) {
    if (n % d === 0n) {
      result.push(d);
      if (d !== n / d) {
        result.push(n / d);
      }
    }
  }
  return result.sort((a, b) => a < b ? -1 : 1);
}

const gaussExamples = [6n, 10n, 12n, 20n];

for (const n of gaussExamples) {
  const divs = divisors(n);
  const phiValues = divs.map(d => euler_phi(d));
  const sum = phiValues.reduce((a, b) => a + b, 0n);

  console.log(`n = ${n}:`);
  console.log(`  Divisors: {${divs.join(", ")}}`);
  console.log(`  phi values: [${phiValues.join(", ")}]`);
  console.log(`  Sum: ${phiValues.join(" + ")} = ${sum}`);
  console.log(`  Equals n: ${sum === n}\n`);
}

console.log("Intuition: Each integer from 1 to n belongs to exactly one");
console.log("class based on gcd(k, n) = d, which corresponds to phi(n/d) elements.");

// ============================================================================
// Example 7: Values of phi for consecutive integers
// ============================================================================
console.log("\n=== Example 7: phi(n) for n = 1 to 30 ===\n");

console.log(" n | phi(n) | prime? | factorization");
console.log("---+--------+--------+---------------");
for (let n = 1n; n <= 30n; n++) {
  const phi = euler_phi(n);
  const prime = is_prime(n);
  const factors = n > 1n ? formatFactorization(factor(n)) : "1";
  console.log(`${n.toString().padStart(2)} | ${phi.toString().padStart(6)} | ${prime.toString().padStart(6)} | ${factors}`);
}

console.log(`
Observations:
1. phi(p) = p - 1 for primes
2. phi(2^k) = 2^(k-1) (powers of 2)
3. phi(n) is even for n > 2
4. phi(n) <= n - 1, with equality iff n is prime
5. phi(n) = n/2 iff n = 2^k for some k >= 1
`);

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
Euler's Totient Function phi(n):

Definition:
  phi(n) = |{k : 1 <= k <= n, gcd(k,n) = 1}|

Key Properties:
1. phi(1) = 1
2. phi(p) = p - 1 for prime p
3. phi(p^k) = p^(k-1) * (p - 1) for prime power
4. phi(mn) = phi(m) * phi(n) when gcd(m,n) = 1 (multiplicative)
5. sum_{d|n} phi(d) = n (Gauss's formula)

Formula from factorization n = p1^e1 * ... * pk^ek:
  phi(n) = n * product of (1 - 1/p_i)
         = product of p_i^(e_i - 1) * (p_i - 1)

Group-theoretic interpretation:
  phi(n) = |(Z/nZ)*| (order of the unit group)
`);
