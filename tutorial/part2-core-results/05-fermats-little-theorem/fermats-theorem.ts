/**
 * # Fermat's Little Theorem
 *
 * Fermat's Little Theorem is one of the most important results in number theory,
 * forming the foundation for many cryptographic algorithms including RSA.
 *
 * ## Statement
 *
 * If p is a prime number and a is any integer not divisible by p, then:
 *
 *   a^(p-1) = 1 (mod p)
 *
 * Equivalently, for any integer a:
 *
 *   a^p = a (mod p)
 *
 * ## Proof Intuition
 *
 * Consider the set {1, 2, 3, ..., p-1} of nonzero residues modulo p.
 * When we multiply each element by a (where gcd(a,p) = 1), we get
 * {a, 2a, 3a, ..., (p-1)a} mod p.
 *
 * Key insight: This is just a permutation of {1, 2, ..., p-1}!
 * Why? Because multiplication by a invertible element permutes the group.
 *
 * Therefore:
 *   1 * 2 * 3 * ... * (p-1) = a * 2a * 3a * ... * (p-1)a  (mod p)
 *   (p-1)! = a^(p-1) * (p-1)!  (mod p)
 *
 * Since (p-1)! is coprime to p, we can cancel it:
 *   1 = a^(p-1)  (mod p)
 *
 * @module tutorial/part2/fermat
 */

import {
  is_prime,
  power_mod,
  gcd,
  euler_phi,
  Zmod,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("FERMAT'S LITTLE THEOREM");
console.log("=".repeat(70));

// ============================================================================
// Example 1: Basic verification of Fermat's Little Theorem
// ============================================================================
console.log("\n=== Example 1: Verifying a^(p-1) = 1 (mod p) ===\n");

const p = 7n;
console.log(`Let p = ${p} (a prime number)`);
console.log(`Fermat's Little Theorem: For any a with gcd(a, ${p}) = 1, a^${p - 1n} = 1 (mod ${p})\n`);

// Test for all nonzero residues modulo 7
for (let a = 1n; a < p; a++) {
  const result = power_mod(a, p - 1n, p);
  console.log(`  ${a}^${p - 1n} mod ${p} = ${result}`);
}

// ============================================================================
// Example 2: The alternative form a^p = a (mod p)
// ============================================================================
console.log("\n=== Example 2: Alternative Form: a^p = a (mod p) ===\n");

const p2 = 11n;
console.log(`Let p = ${p2}`);
console.log(`This form works even when a is divisible by p:\n`);

// Test various values including multiples of p
const testValues = [0n, 1n, 5n, 11n, 15n, 22n, -3n];

for (const a of testValues) {
  const aPowerP = power_mod(((a % p2) + p2) % p2, p2, p2);
  const aMod = ((a % p2) + p2) % p2;
  console.log(`  a = ${a.toString().padStart(3)}: ${a}^${p2} mod ${p2} = ${aPowerP}, a mod ${p2} = ${aMod} [${aPowerP === aMod ? "OK" : "FAIL"}]`);
}

// ============================================================================
// Example 3: Understanding why it fails for composite moduli
// ============================================================================
console.log("\n=== Example 3: Why Fermat's Theorem Requires Prime Modulus ===\n");

const n = 15n; // Composite: 15 = 3 * 5
console.log(`Let n = ${n} (composite, not prime)`);
console.log(`Does a^(n-1) = 1 (mod n) hold for all a coprime to n?\n`);

let counterexampleFound = false;
for (let a = 1n; a < n; a++) {
  if (gcd(a, n) === 1n) {
    const result = power_mod(a, n - 1n, n);
    const status = result === 1n ? "OK" : "FAIL";
    if (result !== 1n && !counterexampleFound) {
      console.log(`  Counterexample: ${a}^${n - 1n} mod ${n} = ${result} != 1`);
      counterexampleFound = true;
    }
  }
}

console.log(`\nFor composite n, we need Euler's theorem instead:`);
console.log(`  phi(${n}) = ${euler_phi(n)}`);
console.log(`  For any a with gcd(a, ${n}) = 1: a^phi(${n}) = 1 (mod ${n})`);

// ============================================================================
// Example 4: Using the ring structure
// ============================================================================
console.log("\n=== Example 4: Using the Zmod Ring Structure ===\n");

const Z13 = Zmod(13n);
console.log(`Working in Z/13Z (the ring of integers mod 13):`);
console.log(`This is a field since 13 is prime: ${Z13.is_field()}\n`);

// Demonstrate Fermat's theorem using ring elements
const elem = Z13.__call__(5n);
console.log(`Let a = ${elem.value} in Z/13Z`);
console.log(`Computing a^(p-1) = 5^12:`);
const fermatResult = elem.pow(12n);
console.log(`  5^12 mod 13 = ${fermatResult.value}`);
console.log(`  This equals 1, as Fermat's theorem predicts!`);

// ============================================================================
// Example 5: The structure of the multiplicative group
// ============================================================================
console.log("\n=== Example 5: The Multiplicative Group (Z/pZ)* ===\n");

const p5 = 7n;
const Zp = Zmod(p5);
console.log(`The multiplicative group (Z/${p5}Z)* has ${p5 - 1n} elements.`);
console.log(`By Lagrange's theorem, every element's order divides |G| = ${p5 - 1n}.`);
console.log(`Therefore a^${p5 - 1n} = 1 for all a in (Z/${p5}Z)*.\n`);

// Show the order of each element
console.log("Finding the order of each element:");
for (let a = 1n; a < p5; a++) {
  let order = 1n;
  let current = a;
  while (current !== 1n) {
    current = (current * a) % p5;
    order++;
  }
  const divides = (p5 - 1n) % order === 0n;
  console.log(`  ord(${a}) = ${order} ${divides ? `(divides ${p5 - 1n})` : ""}`);
}

// ============================================================================
// Example 6: Primality testing preview
// ============================================================================
console.log("\n=== Example 6: Fermat Primality Test (Preview) ===\n");

console.log("Fermat's theorem gives us a way to test primality:");
console.log("If n is prime, then a^(n-1) = 1 (mod n) for all a coprime to n.");
console.log("Contrapositive: If a^(n-1) != 1 (mod n), then n is composite.\n");

function fermatTest(n: bigint, witnesses: bigint[]): { isProbablePrime: boolean; witness?: bigint } {
  if (n < 2n) return { isProbablePrime: false };
  if (n === 2n) return { isProbablePrime: true };
  if (n % 2n === 0n) return { isProbablePrime: false };

  for (const a of witnesses) {
    if (gcd(a, n) !== 1n) continue;
    const result = power_mod(a, n - 1n, n);
    if (result !== 1n) {
      return { isProbablePrime: false, witness: a };
    }
  }
  return { isProbablePrime: true };
}

const testNumbers = [7n, 11n, 15n, 17n, 561n, 1009n];
const witnesses = [2n, 3n, 5n, 7n];

console.log("Testing numbers with witnesses [2, 3, 5, 7]:");
for (const num of testNumbers) {
  const actual = is_prime(num);
  const test = fermatTest(num, witnesses);
  const status = actual === test.isProbablePrime ? "CORRECT" : "INCORRECT";
  console.log(`  ${num}: Fermat test = ${test.isProbablePrime}, Actual = ${actual} [${status}]`);
  if (!test.isProbablePrime && test.witness) {
    console.log(`       Witness: ${test.witness}^${num - 1n} mod ${num} != 1`);
  }
}

console.log("\nNote: 561 is a Carmichael number - it passes the Fermat test");
console.log("for all bases coprime to it, but is actually composite (3 * 11 * 17).");

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
Key Points:
1. Fermat's Little Theorem: a^(p-1) = 1 (mod p) for prime p and gcd(a,p) = 1
2. Alternative form: a^p = a (mod p) for any integer a
3. The theorem fails for composite moduli (use Euler's theorem instead)
4. Applications: modular inverse computation, primality testing, RSA
5. The multiplicative group (Z/pZ)* has order p-1, making Fermat a
   consequence of Lagrange's theorem
`);
