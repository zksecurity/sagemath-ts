/**
 * # The Legendre Symbol
 *
 * The Legendre symbol is a notation for determining whether an integer
 * is a quadratic residue modulo a prime.
 *
 * ## Definition
 *
 * For odd prime p and integer a, the Legendre symbol (a/p) is:
 *   - (a/p) = 0  if p | a
 *   - (a/p) = 1  if a is a quadratic residue mod p
 *   - (a/p) = -1 if a is a quadratic non-residue mod p
 *
 * ## Euler's Formula
 *
 *   (a/p) = a^((p-1)/2) (mod p)
 *
 * @module tutorial/part2/legendre-symbol
 */

import {
  legendre_symbol,
  jacobi_symbol,
  power_mod,
  is_prime,
  gcd,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("THE LEGENDRE SYMBOL");
console.log("=".repeat(70));

// ============================================================================
// Example 1: Definition and Basic Computation
// ============================================================================
console.log("\n=== Example 1: Definition and Basic Computation ===\n");

console.log("The Legendre symbol (a/p) for odd prime p:");
console.log("  (a/p) = 0  if p divides a");
console.log("  (a/p) = 1  if a is a quadratic residue mod p");
console.log("  (a/p) = -1 if a is a quadratic non-residue mod p\n");

const p = 13n;
console.log(`Legendre symbols (a/${p}) for a = 0 to ${p - 1n}:`);
console.log("-".repeat(50));

for (let a = 0n; a < p; a++) {
  const leg = legendre_symbol(a, p);
  let status: string;
  if (leg === 0n) {
    status = "divisible by p";
  } else if (leg === 1n) {
    // Find a square root
    let sqrt = -1n;
    for (let x = 0n; x < p; x++) {
      if ((x * x) % p === a) {
        sqrt = x;
        break;
      }
    }
    status = `QR (sqrt = ${sqrt})`;
  } else {
    status = "QNR";
  }
  console.log(`  (${a.toString().padStart(2)}/${p}) = ${leg.toString().padStart(2)}  [${status}]`);
}

// ============================================================================
// Example 2: Euler's Formula
// ============================================================================
console.log("\n=== Example 2: Euler's Formula ===\n");

console.log("Euler's formula: (a/p) = a^((p-1)/2) mod p");
console.log("This provides an efficient way to compute the Legendre symbol.\n");

const p2 = 17n;
const exp = (p2 - 1n) / 2n;

console.log(`For p = ${p2}, exponent = (p-1)/2 = ${exp}:`);
console.log("-".repeat(60));
console.log("  a | a^8 mod 17 | (a/17) | Match");
console.log("----+------------+--------+------");

for (let a = 1n; a <= 16n; a++) {
  const eulerResult = power_mod(a, exp, p2);
  const legResult = legendre_symbol(a, p2);

  // Convert Euler result to match Legendre convention (-1 instead of p-1)
  const eulerNormalized = eulerResult === p2 - 1n ? -1n : eulerResult;

  const match = eulerNormalized === legResult;
  console.log(` ${a.toString().padStart(2)} |     ${eulerResult.toString().padStart(2)}     |   ${legResult.toString().padStart(2)}   | ${match ? "Yes" : "No"}`);
}

// ============================================================================
// Example 3: Properties of the Legendre Symbol
// ============================================================================
console.log("\n=== Example 3: Properties of the Legendre Symbol ===\n");

const p3 = 11n;

console.log(`Properties for p = ${p3}:\n`);

// Property 1: Multiplicativity
console.log("1. Multiplicativity: (ab/p) = (a/p)(b/p)");
console.log("-".repeat(50));

const testPairs = [
  [2n, 3n],
  [3n, 4n],
  [5n, 6n],
  [7n, 9n],
];

for (const [a, b] of testPairs) {
  const ab = (a * b) % p3;
  const legAB = legendre_symbol(ab, p3);
  const legA = legendre_symbol(a, p3);
  const legB = legendre_symbol(b, p3);
  const product = legA * legB;

  console.log(`  (${a}*${b}/${p3}) = (${ab}/${p3}) = ${legAB}`);
  console.log(`  (${a}/${p3}) * (${b}/${p3}) = ${legA} * ${legB} = ${product}`);
  console.log(`  Match: ${legAB === product}\n`);
}

// Property 2: (a/p) = (a mod p / p)
console.log("\n2. Reduction: (a/p) = (a mod p / p)");
console.log("-".repeat(50));

const largeValues = [23n, 45n, 100n, -3n];

for (const a of largeValues) {
  const reduced = ((a % p3) + p3) % p3;
  const legA = legendre_symbol(a, p3);
  const legReduced = legendre_symbol(reduced, p3);

  console.log(`  (${a}/${p3}) = (${reduced}/${p3}) = ${legA} = ${legReduced}`);
}

// Property 3: (-1/p)
console.log("\n\n3. Special value: (-1/p)");
console.log("-".repeat(50));
console.log("  (-1/p) = 1 if p = 1 (mod 4)");
console.log("  (-1/p) = -1 if p = 3 (mod 4)\n");

const testPrimes = [5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n];

for (const pp of testPrimes) {
  const legMinus1 = legendre_symbol(-1n, pp);
  const pMod4 = pp % 4n;
  const expected = pMod4 === 1n ? 1n : -1n;

  console.log(`  p = ${pp.toString().padStart(2)}: (-1/${pp}) = ${legMinus1.toString().padStart(2)}, ` +
              `p mod 4 = ${pMod4}, expected ${expected} [${legMinus1 === expected ? "OK" : "FAIL"}]`);
}

// Property 4: (2/p)
console.log("\n\n4. Special value: (2/p)");
console.log("-".repeat(50));
console.log("  (2/p) = 1 if p = +/-1 (mod 8)");
console.log("  (2/p) = -1 if p = +/-3 (mod 8)\n");

for (const pp of testPrimes) {
  const leg2 = legendre_symbol(2n, pp);
  const pMod8 = pp % 8n;
  const expected = (pMod8 === 1n || pMod8 === 7n) ? 1n : -1n;

  console.log(`  p = ${pp.toString().padStart(2)}: (2/${pp}) = ${leg2.toString().padStart(2)}, ` +
              `p mod 8 = ${pMod8}, expected ${expected} [${leg2 === expected ? "OK" : "FAIL"}]`);
}

// ============================================================================
// Example 4: Computing Legendre Symbols Efficiently
// ============================================================================
console.log("\n=== Example 4: Efficient Computation ===\n");

console.log("Methods for computing (a/p):\n");

console.log("1. Euler's criterion: a^((p-1)/2) mod p");
console.log("   Time: O(log p) modular exponentiations\n");

console.log("2. Using properties + quadratic reciprocity");
console.log("   Often faster for specific values\n");

// Example: Compute (1234567/104729) where 104729 is prime
const largePrime = 104729n; // 10000th prime
const testValue = 1234567n;

console.log(`Computing (${testValue}/${largePrime}):`);

// Method 1: Euler
const euler = power_mod(testValue, (largePrime - 1n) / 2n, largePrime);
const eulerResult = euler === largePrime - 1n ? -1n : euler;

// Method 2: Library
const libResult = legendre_symbol(testValue, largePrime);

console.log(`  Euler's criterion: ${eulerResult}`);
console.log(`  Library result:    ${libResult}`);
console.log(`  Match: ${eulerResult === libResult}`);

// ============================================================================
// Example 5: Legendre Symbol Table
// ============================================================================
console.log("\n=== Example 5: Legendre Symbol Tables ===\n");

const p5 = 7n;

console.log(`Complete Legendre symbol table for p = ${p5}:\n`);

// Header
let header = "  a |";
for (let a = 0n; a < p5; a++) {
  header += a.toString().padStart(3);
}
console.log(header);
console.log("----+" + "-".repeat(Number(p5) * 3));

// Row for (a/p)
let row = "(a/p)|";
for (let a = 0n; a < p5; a++) {
  const leg = legendre_symbol(a, p5);
  row += leg.toString().padStart(3);
}
console.log(row);

// Now show pattern for several primes
console.log("\n\nLegendre symbols for first few primes:");
console.log("-".repeat(70));

const smallPrimes = [3n, 5n, 7n, 11n, 13n];

for (const pp of smallPrimes) {
  let qr = "";
  let qnr = "";

  for (let a = 1n; a < pp; a++) {
    const leg = legendre_symbol(a, pp);
    if (leg === 1n) {
      qr += a.toString() + " ";
    } else {
      qnr += a.toString() + " ";
    }
  }

  console.log(`p = ${pp.toString().padStart(2)}: QR = {${qr.trim()}}, QNR = {${qnr.trim()}}`);
}

// ============================================================================
// Example 6: Applications
// ============================================================================
console.log("\n=== Example 6: Applications ===\n");

console.log("1. Primality Testing (Solovay-Strassen):");
console.log("   If n is prime and gcd(a, n) = 1, then");
console.log("   (a/n) = a^((n-1)/2) mod n");
console.log("   If this fails, n is definitely composite.\n");

// Test some numbers
function solovayStrassenWitness(n: bigint, a: bigint): boolean {
  if (gcd(a, n) !== 1n) return false; // Not a valid witness

  // Solovay-Strassen is defined via the Jacobi symbol: n is the number under
  // test, so it may well be composite. legendre_symbol() rejects composite
  // moduli (as SageMath does), so jacobi_symbol() is the correct call here.
  const jacobi = jacobi_symbol(a, n);
  const euler = power_mod(a, (n - 1n) / 2n, n);

  // Convert jacobi to match euler's range
  const jacobiMod = jacobi === -1n ? n - 1n : (jacobi === 0n ? 0n : 1n);

  return euler === jacobiMod;
}

const testNumbers = [7n, 11n, 15n, 21n, 561n];
const witness = 2n;

console.log(`Testing with witness a = ${witness}:`);
for (const n of testNumbers) {
  if (n % 2n === 0n) continue; // Skip even numbers

  const passes = solovayStrassenWitness(n, witness);
  const isPrime = is_prime(n);
  const status = passes && !isPrime ? "FOOLED" : (passes === isPrime ? "OK" : "Detected");

  console.log(`  n = ${n.toString().padStart(3)}: passes = ${passes.toString().padStart(5)}, ` +
              `prime = ${isPrime.toString().padStart(5)} [${status}]`);
}

console.log("\n\n2. Determining solvability of x^2 = a (mod p):");
console.log("   Solvable iff (a/p) = 1\n");

const p6 = 23n;
console.log(`Which values have square roots mod ${p6}?`);
for (let a = 1n; a <= 10n; a++) {
  const leg = legendre_symbol(a, p6);
  console.log(`  ${a}: (${a}/${p6}) = ${leg.toString().padStart(2)} => ${leg === 1n ? "has square root" : "no square root"}`);
}

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
The Legendre Symbol (a/p):

Definition (p odd prime):
  (a/p) = 0  if p | a
  (a/p) = 1  if a is QR mod p
  (a/p) = -1 if a is QNR mod p

Computation:
  (a/p) = a^((p-1)/2) mod p (Euler's criterion)

Key Properties:
  1. Multiplicativity: (ab/p) = (a/p)(b/p)
  2. Periodicity: (a/p) = (a mod p / p)
  3. (-1/p) = (-1)^((p-1)/2) = 1 if p = 1 mod 4, -1 if p = 3 mod 4
  4. (2/p) = (-1)^((p^2-1)/8) = 1 if p = +/-1 mod 8, -1 if p = +/-3 mod 8

Applications:
  - Determining if x^2 = a (mod p) is solvable
  - Solovay-Strassen primality test
  - Quadratic reciprocity law
  - Various cryptographic protocols

Using sagemath-ts:
  legendre_symbol(a, p) computes (a/p)
`);
