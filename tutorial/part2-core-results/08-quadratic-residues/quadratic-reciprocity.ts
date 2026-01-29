/**
 * # The Law of Quadratic Reciprocity
 *
 * Quadratic reciprocity is one of the most important and beautiful
 * theorems in number theory. It relates the solvability of x^2 = p (mod q)
 * to the solvability of x^2 = q (mod p).
 *
 * ## Main Statement
 *
 * For distinct odd primes p and q:
 *
 *   (p/q)(q/p) = (-1)^((p-1)(q-1)/4)
 *
 * This equals 1 unless both p and q are congruent to 3 (mod 4).
 *
 * ## First Supplement
 *
 *   (-1/p) = (-1)^((p-1)/2) = { 1 if p = 1 (mod 4)
 *                             {-1 if p = 3 (mod 4)
 *
 * ## Second Supplement
 *
 *   (2/p) = (-1)^((p^2-1)/8) = { 1 if p = +/-1 (mod 8)
 *                              {-1 if p = +/-3 (mod 8)
 *
 * @module tutorial/part2/quadratic-reciprocity
 */

import {
  legendre_symbol,
  jacobi_symbol,
  is_prime,
  power_mod,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("THE LAW OF QUADRATIC RECIPROCITY");
console.log("=".repeat(70));

// ============================================================================
// The Main Theorem
// ============================================================================
console.log("\n=== The Main Theorem ===\n");

console.log("For distinct odd primes p and q:");
console.log("");
console.log("  (p/q) * (q/p) = (-1)^((p-1)(q-1)/4)");
console.log("");
console.log("Equivalently:");
console.log("  (p/q) = (q/p)     if p = 1 (mod 4) OR q = 1 (mod 4)");
console.log("  (p/q) = -(q/p)    if p = 3 (mod 4) AND q = 3 (mod 4)");
console.log("");

// ============================================================================
// Example 1: Verifying Quadratic Reciprocity
// ============================================================================
console.log("\n=== Example 1: Verifying the Main Theorem ===\n");

const primePairs: [bigint, bigint][] = [
  [3n, 5n],
  [5n, 7n],
  [7n, 11n],
  [3n, 7n],
  [11n, 13n],
  [3n, 11n],
  [13n, 17n],
  [5n, 11n],
];

console.log("  p |  q | (p/q) | (q/p) | product | (-1)^((p-1)(q-1)/4) | Match");
console.log("----+----+-------+-------+---------+--------------------+------");

for (const [p, q] of primePairs) {
  const legPQ = legendre_symbol(p, q);
  const legQP = legendre_symbol(q, p);
  const product = legPQ * legQP;

  const exponent = ((p - 1n) * (q - 1n)) / 4n;
  const expected = exponent % 2n === 0n ? 1n : -1n;

  const match = product === expected;

  console.log(`${p.toString().padStart(3)} |${q.toString().padStart(3)} |   ${legPQ.toString().padStart(2)}  |   ${legQP.toString().padStart(2)}  |    ${product.toString().padStart(2)}   |         ${expected.toString().padStart(2)}          | ${match ? "Yes" : "No"}`);
}

// ============================================================================
// Example 2: The First Supplement (-1/p)
// ============================================================================
console.log("\n=== Example 2: First Supplement (-1/p) ===\n");

console.log("(-1/p) = (-1)^((p-1)/2) = 1 if p = 1 (mod 4), -1 if p = 3 (mod 4)");
console.log("");
console.log("This tells us when -1 is a quadratic residue mod p.\n");

const testPrimes = [3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n];

console.log("  p  | p mod 4 | (-1/p) | -1 is QR?");
console.log("-----+---------+--------+----------");

for (const p of testPrimes) {
  const leg = legendre_symbol(-1n, p);
  const pMod4 = p % 4n;
  const expected = pMod4 === 1n ? 1n : -1n;
  const isQR = leg === 1n;

  console.log(` ${p.toString().padStart(2)}  |    ${pMod4}    |   ${leg.toString().padStart(2)}   | ${isQR ? "Yes" : "No "}`);
}

// Consequences
console.log("\n--> -1 is a QR mod p iff p = 1 (mod 4)");
console.log("--> Primes p = 1 (mod 4) can be written as sum of two squares");
console.log("    Example: 5 = 1 + 4 = 1^2 + 2^2");
console.log("    Example: 13 = 4 + 9 = 2^2 + 3^2");

// ============================================================================
// Example 3: The Second Supplement (2/p)
// ============================================================================
console.log("\n=== Example 3: Second Supplement (2/p) ===\n");

console.log("(2/p) = (-1)^((p^2-1)/8) =  1 if p = +/-1 (mod 8)");
console.log("                          = -1 if p = +/-3 (mod 8)");
console.log("");
console.log("This tells us when 2 is a quadratic residue mod p.\n");

console.log("  p  | p mod 8 | (2/p) | 2 is QR?");
console.log("-----+---------+-------+----------");

for (const p of testPrimes) {
  const leg = legendre_symbol(2n, p);
  const pMod8 = p % 8n;
  const expected = (pMod8 === 1n || pMod8 === 7n) ? 1n : -1n;
  const isQR = leg === 1n;

  console.log(` ${p.toString().padStart(2)}  |    ${pMod8}    |   ${leg.toString().padStart(2)}  | ${isQR ? "Yes" : "No "}`);
}

// ============================================================================
// Example 4: Using QR to Simplify Computation
// ============================================================================
console.log("\n=== Example 4: Using QR to Simplify Computation ===\n");

console.log("Problem: Is 219 a quadratic residue mod 383?\n");

const a = 219n;
const p = 383n;

console.log(`Direct approach: Compute 219^((383-1)/2) mod 383`);
const directResult = power_mod(a, (p - 1n) / 2n, p);
const directAnswer = directResult === 1n ? 1n : -1n;
console.log(`  Result: ${directResult} --> (219/383) = ${directAnswer}\n`);

console.log("Using quadratic reciprocity:");
console.log(`  (219/383) = (3/383) * (73/383) since 219 = 3 * 73`);

const leg3 = legendre_symbol(3n, p);
const leg73 = legendre_symbol(73n, p);

console.log(`  (3/383): Apply QR with p=3, q=383`);
console.log(`    383 mod 4 = ${p % 4n}, so (3/383) = (383/3) = (${p % 3n}/3) = ${legendre_symbol(p % 3n, 3n)}`);

console.log(`  (73/383): Apply QR with p=73, q=383`);
console.log(`    73 mod 4 = 1, so (73/383) = (383/73)`);
console.log(`    383 mod 73 = ${p % 73n}`);
console.log(`    (383/73) = (${p % 73n}/73)`);

// Continue reduction...
const remainder = p % 73n;
console.log(`    Need (${remainder}/73)...`);

const finalLeg3 = legendre_symbol(3n, p);
const finalLeg73 = legendre_symbol(73n, p);
console.log(`\n  (3/383) = ${finalLeg3}`);
console.log(`  (73/383) = ${finalLeg73}`);
console.log(`  (219/383) = ${finalLeg3} * ${finalLeg73} = ${finalLeg3 * finalLeg73}`);
console.log(`\n  Direct computation agrees: ${directAnswer === finalLeg3 * finalLeg73}`);

// ============================================================================
// Example 5: The Quadratic Reciprocity Algorithm
// ============================================================================
console.log("\n=== Example 5: Efficient Legendre Symbol via QR ===\n");

console.log("The Jacobi symbol algorithm essentially uses quadratic reciprocity");
console.log("repeatedly to reduce the computation, similar to the Euclidean algorithm.\n");

function legendreWithQR(a: bigint, p: bigint): { result: bigint; steps: string[] } {
  const steps: string[] = [];
  let result = 1n;

  steps.push(`Computing (${a}/${p})`);

  // Reduce a mod p
  a = ((a % p) + p) % p;
  steps.push(`  Reduce: (${a}/${p})`);

  while (a !== 0n && a !== 1n) {
    // Factor out 2s
    let twos = 0n;
    while (a % 2n === 0n) {
      a /= 2n;
      twos++;
    }

    if (twos > 0n) {
      // Apply second supplement
      const twoSym = (p % 8n === 1n || p % 8n === 7n) ? 1n : -1n;
      if (twos % 2n === 1n) {
        result *= twoSym;
        steps.push(`  Factor out 2^${twos}: (2/${p}) = ${twoSym}`);
      } else {
        steps.push(`  Factor out 2^${twos}: (2/${p})^${twos} = 1`);
      }
      steps.push(`  Now: (${a}/${p})`);
    }

    if (a === 1n) break;

    // Apply quadratic reciprocity
    const signChange = (a % 4n === 3n && p % 4n === 3n) ? -1n : 1n;
    result *= signChange;

    steps.push(`  Apply QR (sign: ${signChange}): (${a}/${p}) -> (${p % a}/${a})`);

    [a, p] = [p % a, a];
  }

  steps.push(`  Result: ${result}`);
  return { result, steps };
}

const examples: [bigint, bigint][] = [
  [23n, 59n],
  [7n, 31n],
  [11n, 47n],
];

for (const [a, p] of examples) {
  const { result, steps } = legendreWithQR(a, p);
  console.log(steps.join("\n"));
  console.log(`  Library verification: ${legendre_symbol(a, p)}\n`);
}

// ============================================================================
// Example 6: Historical Context
// ============================================================================
console.log("\n=== Example 6: Historical Context ===\n");

console.log("Quadratic reciprocity was conjectured by Euler and Legendre,");
console.log("and first proven by Gauss in 1796 (at age 18!).");
console.log("");
console.log("Gauss called it the 'golden theorem' (theorema aureum) and gave");
console.log("at least 8 different proofs throughout his life.");
console.log("");
console.log("Key milestones:");
console.log("  - 1785: Legendre states the law (incomplete proof)");
console.log("  - 1796: Gauss provides first complete proof");
console.log("  - 1872: Eisenstein gives an elegant proof using roots of unity");
console.log("  - Today: Over 200 different proofs are known!");
console.log("");
console.log("Generalizations:");
console.log("  - Cubic reciprocity (Eisenstein, 1844)");
console.log("  - Quartic reciprocity (Gauss, Jacobi)");
console.log("  - Higher reciprocity laws (Artin, class field theory)");

// ============================================================================
// Example 7: Applications
// ============================================================================
console.log("\n=== Example 7: Cryptographic Applications ===\n");

console.log("1. Efficient computation of Legendre/Jacobi symbols");
console.log("   - Without QR, need modular exponentiation: O(log^3 p)");
console.log("   - With QR, like Euclidean algorithm: O(log^2 p)");
console.log("");
console.log("2. Solovay-Strassen primality test");
console.log("   - Tests (a/n) = a^((n-1)/2) mod n");
console.log("   - Uses efficient Jacobi symbol computation");
console.log("");
console.log("3. Quadratic sieve factoring algorithm");
console.log("   - Uses smooth numbers where (a/p) = 1 for small primes p");
console.log("   - QR helps determine which primes to use");
console.log("");
console.log("4. Elliptic curve cryptography");
console.log("   - Point counting uses character sums");
console.log("   - Quadratic character appears naturally");
console.log("");

// Demonstrate Jacobi vs exponentiation speed
console.log("Speed comparison (conceptual):");
const largePrime = 104729n; // 10000th prime
const testVal = 12345n;

const jacobiResult = jacobi_symbol(testVal, largePrime);
const eulerResult = power_mod(testVal, (largePrime - 1n) / 2n, largePrime);
const eulerNorm = eulerResult === largePrime - 1n ? -1n : eulerResult;

console.log(`  Jacobi(${testVal}, ${largePrime}) = ${jacobiResult}`);
console.log(`  Euler(${testVal}, ${largePrime}) = ${eulerNorm}`);
console.log(`  Both give the same answer, but Jacobi uses fewer operations.`);

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
The Law of Quadratic Reciprocity:

Main Theorem (p, q distinct odd primes):
  (p/q)(q/p) = (-1)^((p-1)(q-1)/4)

Interpretation:
  - If p = 1 (mod 4) OR q = 1 (mod 4): (p/q) = (q/p)
  - If p = 3 (mod 4) AND q = 3 (mod 4): (p/q) = -(q/p)

First Supplement (-1/p):
  (-1/p) = 1 if p = 1 (mod 4)
  (-1/p) = -1 if p = 3 (mod 4)

Second Supplement (2/p):
  (2/p) = 1 if p = +/-1 (mod 8)
  (2/p) = -1 if p = +/-3 (mod 8)

Significance:
  - Reduces Legendre symbol computation from O(log^3 p) to O(log^2 p)
  - Foundation for many number-theoretic algorithms
  - Generalizes to higher reciprocity laws
  - One of the most beautiful theorems in mathematics

Historical Note:
  - "Golden theorem" (Gauss)
  - Over 200 known proofs
  - Gateway to algebraic number theory
`);
