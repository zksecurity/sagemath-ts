/**
 * # Modular Congruence and Equivalence Classes
 *
 * Modular arithmetic is the foundation of modern cryptography. Understanding
 * congruence relations and equivalence classes is essential for working with
 * finite fields, RSA, elliptic curves, and virtually every cryptographic system.
 *
 * ## Mathematical Background
 *
 * We say a is congruent to b modulo n (written a = b (mod n)) if n divides (a - b).
 * Equivalently: a = b (mod n) iff a mod n = b mod n
 *
 * Properties of congruence:
 *   - Reflexive: a = a (mod n)
 *   - Symmetric: if a = b (mod n), then b = a (mod n)
 *   - Transitive: if a = b (mod n) and b = c (mod n), then a = c (mod n)
 *
 * These properties make congruence an equivalence relation, which partitions
 * the integers into equivalence classes.
 *
 * ## Why This Matters for Cryptography
 *
 *   - All arithmetic in RSA is done modulo n = pq
 *   - Elliptic curve operations are done modulo a prime p
 *   - Hash functions produce outputs in a fixed range (equivalence class representatives)
 *   - Finite field arithmetic is fundamentally about equivalence classes
 */

import { gcd, Mod, Zmod, IntegerModRing } from 'sagemath-ts';

// =============================================================================
// Example 1: Understanding Congruence
// =============================================================================
console.log("=== Example 1: Understanding Congruence ===\n");

/**
 * Two integers are congruent modulo n if they have the same remainder
 * when divided by n, or equivalently, if their difference is divisible by n.
 */
function areCongruent(a: bigint, b: bigint, n: bigint): boolean {
  return (a - b) % n === 0n;
}

const modulus = 7n;
console.log(`Working with modulus n = ${modulus}\n`);

const testPairs: [bigint, bigint][] = [
  [17n, 3n],      // 17 - 3 = 14 = 2*7
  [25n, 4n],      // 25 - 4 = 21 = 3*7
  [100n, 2n],     // 100 - 2 = 98 = 14*7
  [15n, 3n],      // 15 - 3 = 12, not divisible by 7
];

console.log("Testing congruences mod 7:");
for (const [a, b] of testPairs) {
  const diff = a - b;
  const congruent = areCongruent(a, b, modulus);
  console.log(`  ${a} = ${b} (mod ${modulus})? ${congruent}`);
  console.log(`    Difference: ${a} - ${b} = ${diff}, divisible by ${modulus}? ${diff % modulus === 0n}\n`);
}

// =============================================================================
// Example 2: Equivalence Classes
// =============================================================================
console.log("=== Example 2: Equivalence Classes ===\n");

/**
 * The equivalence class [a]_n is the set of all integers congruent to a mod n.
 * [a]_n = {..., a-2n, a-n, a, a+n, a+2n, ...}
 *
 * There are exactly n distinct equivalence classes: [0], [1], ..., [n-1]
 */
const n = 5n;
console.log(`The ${n} equivalence classes modulo ${n}:\n`);

for (let i = 0n; i < n; i++) {
  // Show some representatives from each class
  const representatives: bigint[] = [];
  for (let k = -2n; k <= 2n; k++) {
    representatives.push(i + k * n);
  }
  console.log(`  [${i}]_${n} = {..., ${representatives.join(', ')}, ...}`);
}

console.log("\nEvery integer belongs to exactly one equivalence class!");

console.log();

// =============================================================================
// Example 3: Creating the Ring Z/nZ
// =============================================================================
console.log("=== Example 3: The Ring Z/nZ ===\n");

/**
 * Z/nZ is the set of equivalence classes {[0], [1], ..., [n-1]} with
 * arithmetic operations defined on representatives.
 *
 * The Zmod() function creates this ring.
 */
const Z7 = Zmod(7);
console.log(`Created ring: ${Z7.toString()}`);
console.log(`Cardinality: ${Z7.cardinality()}`);
console.log(`Is a field? ${Z7.is_field()} (7 is prime)`);

console.log("\nElements of Z/7Z:");
for (const elem of Z7) {
  console.log(`  ${elem.value}`);
}

// Create a non-prime modulus ring
const Z12 = Zmod(12);
console.log(`\nCreated ring: ${Z12.toString()}`);
console.log(`Is a field? ${Z12.is_field()} (12 is not prime)`);

console.log();

// =============================================================================
// Example 4: Coercion and Representatives
// =============================================================================
console.log("=== Example 4: Coercion and Representatives ===\n");

/**
 * When we create an element of Z/nZ, any integer is coerced to its
 * canonical representative in {0, 1, ..., n-1}.
 */
const Z10 = Zmod(10);

console.log("Coercing various integers to Z/10Z:");
const integers = [3n, 17n, 25n, -3n, -17n, 100n, 0n];

for (const x of integers) {
  const elem = Z10.__call__(x);
  console.log(`  ${x} -> ${elem.value} (in Z/10Z)`);
}

console.log("\nThe canonical representative is always in [0, n-1].");

console.log();

// =============================================================================
// Example 5: Negative Numbers and Modular Reduction
// =============================================================================
console.log("=== Example 5: Handling Negative Numbers ===\n");

/**
 * Different programming languages handle negative modulo differently!
 * In cryptography, we always want the result in [0, n-1].
 *
 * For -3 mod 7: some languages give -3, but we want 4.
 * The sagemath-ts library always gives the positive representative.
 */
function properMod(a: bigint, n: bigint): bigint {
  const result = a % n;
  return result < 0n ? result + n : result;
}

const mod = 7n;
const negatives = [-1n, -3n, -7n, -10n, -100n];

console.log("Proper modular reduction of negative numbers (mod 7):\n");
for (const x of negatives) {
  const jsResult = x % mod;        // JavaScript's built-in (can be negative)
  const properResult = properMod(x, mod);  // Always in [0, n-1]
  const libResult = Mod(x, mod).value;      // Library result

  console.log(`  ${x} mod ${mod}:`);
  console.log(`    JavaScript %: ${jsResult}`);
  console.log(`    Proper: ${properResult}`);
  console.log(`    Library: ${libResult}`);
  console.log();
}

// =============================================================================
// Example 6: Equivalence Class Representatives
// =============================================================================
console.log("=== Example 6: Finding Other Representatives ===\n");

/**
 * Any integer in an equivalence class can represent that class.
 * Sometimes we want different representatives for different purposes.
 */
const Z11 = Zmod(11);
const element = Z11.__call__(4);

console.log(`Element 4 in Z/11Z`);
console.log(`Canonical representative: ${element.value}`);

console.log("\nOther representatives of [4]_11:");
for (let k = -3n; k <= 3n; k++) {
  const rep = 4n + k * 11n;
  console.log(`  4 + ${k}*11 = ${rep}`);
}

// Balanced representatives (closest to 0)
console.log("\nBalanced representatives (in [-n/2, n/2]):");
for (const elem of Z11) {
  const v = elem.value;
  const balanced = v > 5n ? v - 11n : v;
  console.log(`  ${v} -> ${balanced} (balanced)`);
}

console.log();

// =============================================================================
// Example 7: Congruence Preserves Operations
// =============================================================================
console.log("=== Example 7: Congruence Preserves Operations ===\n");

/**
 * If a = a' (mod n) and b = b' (mod n), then:
 *   - a + b = a' + b' (mod n)
 *   - a - b = a' - b' (mod n)
 *   - a * b = a' * b' (mod n)
 *   - a^k = (a')^k (mod n)
 *
 * This is why modular arithmetic works!
 */
const Z17 = Zmod(17);

// Two representations of the same element
const a1 = 5n;
const a2 = 22n;  // 5 + 17 = 22, same class
const b1 = 8n;
const b2 = -9n;  // 8 - 17 = -9, same class

console.log(`a: ${a1} = ${a2} (mod 17)`);
console.log(`b: ${b1} = ${b2} (mod 17)\n`);

console.log("Checking that operations preserve congruence:");

// Addition
const sum1 = (a1 + b1) % 17n;
const sum2 = properMod(a2 + b2, 17n);
console.log(`  Addition: ${a1}+${b1}=${(a1+b1)} = ${sum1}, ${a2}+${b2}=${a2+b2} = ${sum2} (mod 17)`);

// Multiplication
const prod1 = (a1 * b1) % 17n;
const prod2 = properMod(a2 * b2, 17n);
console.log(`  Multiply: ${a1}*${b1}=${a1*b1} = ${prod1}, ${a2}*${b2}=${a2*b2} = ${prod2} (mod 17)`);

// Power
const pow1 = (a1 ** 3n) % 17n;
const pow2 = properMod(a2 ** 3n, 17n);
console.log(`  Power: ${a1}^3=${a1**3n} = ${pow1}, ${a2}^3=${a2**3n} = ${pow2} (mod 17)`);

console.log();

// =============================================================================
// Example 8: Residue Systems
// =============================================================================
console.log("=== Example 8: Complete Residue Systems ===\n");

/**
 * A complete residue system modulo n is a set of n integers, one from
 * each equivalence class. The standard system is {0, 1, ..., n-1}.
 *
 * A reduced residue system contains only integers coprime to n.
 */
const modN = 12n;
console.log(`For n = ${modN}:\n`);

// Complete residue system
console.log(`Complete residue system: {0, 1, 2, ..., ${modN - 1n}}`);

// Reduced residue system (units)
const reduced: bigint[] = [];
for (let i = 1n; i < modN; i++) {
  if (gcd(i, modN) === 1n) {
    reduced.push(i);
  }
}
console.log(`Reduced residue system: {${reduced.join(', ')}}`);
console.log(`Size of reduced system = phi(${modN}) = ${reduced.length}`);

// Alternative complete residue system
console.log(`\nAlternative complete system (balanced): {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6}`);

console.log();

// =============================================================================
// Example 9: Congruence and Divisibility
// =============================================================================
console.log("=== Example 9: Congruence and Divisibility ===\n");

/**
 * a = b (mod n) means n | (a - b), i.e., n divides (a - b).
 * This connects congruence to divisibility.
 */
console.log("Connection between congruence and divisibility:\n");

const examples: [bigint, bigint, bigint][] = [
  [20n, 6n, 7n],
  [100n, 1n, 9n],
  [2024n, 0n, 4n],
];

for (const [a, b, m] of examples) {
  const diff = a - b;
  const quotient = diff / m;
  console.log(`  ${a} = ${b} (mod ${m})`);
  console.log(`    Because ${a} - ${b} = ${diff} = ${m} * ${quotient}`);
  console.log(`    So ${m} divides ${diff}\n`);
}

// =============================================================================
// Exercise: Prove Equivalence Class Properties
// =============================================================================
console.log("=== Exercise: Equivalence Class Properties ===\n");

/**
 * EXERCISE: Verify that congruence is an equivalence relation by
 * testing the three properties for specific examples.
 */

// TODO: Test these properties
function testEquivalenceRelation(n: bigint): void {
  console.log(`Testing equivalence relation properties for mod ${n}:\n`);

  const a = 5n, b = 12n, c = 19n;  // All congruent mod 7: 5 = 12 = 19 (mod 7)

  // 1. Reflexivity: a = a (mod n)
  const reflexive = areCongruent(a, a, n);
  console.log(`  Reflexive: ${a} = ${a} (mod ${n})? ${reflexive}`);

  // 2. Symmetry: if a = b (mod n), then b = a (mod n)
  const abCongruent = areCongruent(a, b, n);
  const baCongruent = areCongruent(b, a, n);
  console.log(`  Symmetric: ${a} = ${b} (mod ${n})? ${abCongruent}`);
  console.log(`             ${b} = ${a} (mod ${n})? ${baCongruent}`);
  console.log(`             Both directions match? ${abCongruent === baCongruent}`);

  // 3. Transitivity: if a = b (mod n) and b = c (mod n), then a = c (mod n)
  const bcCongruent = areCongruent(b, c, n);
  const acCongruent = areCongruent(a, c, n);
  console.log(`  Transitive: ${a}=${b}? ${abCongruent}, ${b}=${c}? ${bcCongruent}, ${a}=${c}? ${acCongruent}`);
}

testEquivalenceRelation(7n);

console.log();

// =============================================================================
// Cryptographic Application: Hash Function Output
// =============================================================================
console.log("=== Cryptographic Application: Hash Outputs as Equivalence Classes ===\n");

/**
 * When a hash function produces a 256-bit output, it's selecting a
 * representative from one of 2^256 equivalence classes.
 *
 * The collision resistance property means it's hard to find two
 * inputs that map to the same equivalence class.
 */
console.log("Hash functions as equivalence class selectors:\n");

// Simulate a simple hash (NOT cryptographically secure!)
function simpleHash(input: bigint): bigint {
  // Just for illustration - compresses to a 4-bit output
  const mod = 16n;
  const mixed = (input * 7n + 3n) % mod;
  return mixed;
}

console.log("Simple 4-bit hash (for illustration):");
console.log("Many inputs map to the same output (equivalence class):\n");

const hashMod = 16n;
const buckets: Map<bigint, bigint[]> = new Map();

for (let i = 0n; i < 50n; i++) {
  const h = simpleHash(i);
  if (!buckets.has(h)) {
    buckets.set(h, []);
  }
  buckets.get(h)!.push(i);
}

for (const [hash, inputs] of buckets.entries()) {
  console.log(`  Hash ${hash}: inputs {${inputs.slice(0, 5).join(', ')}${inputs.length > 5 ? ', ...' : ''}}`);
}

console.log("\nAll inputs in a bucket are 'collisions' - they're in the same equivalence class!");
