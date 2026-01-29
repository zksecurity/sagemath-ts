/**
 * # Fields: Prime Fields Z/pZ
 *
 * A field is a ring where every non-zero element has a multiplicative inverse.
 * Prime fields Z/pZ are the simplest finite fields and are fundamental to
 * cryptography.
 *
 * ## Mathematical Background
 *
 * A FIELD is a set F with operations + and * such that:
 *   1. (F, +) is an abelian group with identity 0
 *   2. (F - {0}, *) is an abelian group with identity 1
 *   3. Multiplication distributes over addition
 *
 * Key properties:
 *   - Every non-zero element is invertible
 *   - No zero divisors (integral domain)
 *   - Division is always defined (except by zero)
 *
 * Z/pZ is a field iff p is prime. We denote it GF(p) or F_p.
 *
 * ## Why This Matters for Cryptography
 *
 *   - Elliptic curves are typically defined over prime fields
 *   - Diffie-Hellman uses the multiplicative group of a prime field
 *   - Many protocols need division, which requires fields
 *   - Polynomial arithmetic over fields (e.g., AES)
 */

import { gcd, is_prime, inverse_mod, power_mod, euler_phi, Mod, factor, Zmod } from 'sagemath-ts';

// =============================================================================
// Example 1: What Makes a Field
// =============================================================================
console.log("=== Example 1: What Makes a Field ===\n");

/**
 * Z/nZ is a field iff n is prime.
 * In a field, every non-zero element has a multiplicative inverse.
 */
const moduli = [4n, 5n, 6n, 7n];

console.log("Checking which Z/nZ are fields:\n");

for (const n of moduli) {
  const ring = Zmod(n);
  const isField = ring.is_field();
  const isPrime = is_prime(n);

  console.log(`Z/${n}Z:`);
  console.log(`  Is n prime? ${isPrime}`);
  console.log(`  Is Z/nZ a field? ${isField}`);

  if (isField) {
    // Show all inverses exist
    const inverses: string[] = [];
    for (let a = 1n; a < n; a++) {
      const inv = inverse_mod(a, n);
      inverses.push(`${a}^-1=${inv}`);
    }
    console.log(`  Inverses: ${inverses.join(', ')}`);
  } else {
    // Show missing inverses
    const noInverse: bigint[] = [];
    for (let a = 1n; a < n; a++) {
      if (gcd(a, n) !== 1n) {
        noInverse.push(a);
      }
    }
    console.log(`  Elements without inverse: {${noInverse.join(', ')}}`);
  }
  console.log();
}

// =============================================================================
// Example 2: The Prime Field GF(p)
// =============================================================================
console.log("=== Example 2: The Prime Field GF(p) ===\n");

/**
 * GF(p) denotes the finite field with p elements.
 * For prime p, GF(p) = Z/pZ.
 */
const p = 7n;
const GF7 = Zmod(p);

console.log(`The field GF(${p}) = Z/${p}Z:\n`);

// Addition table
console.log("Addition table:");
let addHeader = "  + |";
for (let j = 0n; j < p; j++) addHeader += ` ${j}`;
console.log(addHeader);
console.log("----+--------------");
for (let i = 0n; i < p; i++) {
  let row = `  ${i} |`;
  for (let j = 0n; j < p; j++) {
    row += ` ${(i + j) % p}`;
  }
  console.log(row);
}

console.log("\nMultiplication table:");
let mulHeader = "  * |";
for (let j = 0n; j < p; j++) mulHeader += ` ${j}`;
console.log(mulHeader);
console.log("----+--------------");
for (let i = 0n; i < p; i++) {
  let row = `  ${i} |`;
  for (let j = 0n; j < p; j++) {
    row += ` ${(i * j) % p}`;
  }
  console.log(row);
}

console.log();

// =============================================================================
// Example 3: Field Division
// =============================================================================
console.log("=== Example 3: Division in a Field ===\n");

/**
 * In a field, division a/b is defined as a * b^(-1) for b != 0.
 * This is NOT possible in general rings!
 */
const divP = 11n;
console.log(`Division in GF(${divP}):\n`);

const divisions: [bigint, bigint][] = [
  [1n, 3n],
  [5n, 2n],
  [7n, 4n],
  [10n, 10n],
];

for (const [a, b] of divisions) {
  const bInv = inverse_mod(b, divP);
  const quotient = (a * bInv) % divP;

  console.log(`  ${a} / ${b}:`);
  console.log(`    ${b}^(-1) = ${bInv}`);
  console.log(`    ${a} * ${bInv} = ${quotient} mod ${divP}`);
  console.log(`    Verify: ${quotient} * ${b} = ${(quotient * b) % divP}\n`);
}

console.log("Division is always possible in a field (except by 0)!");

console.log();

// =============================================================================
// Example 4: Field Characteristic
// =============================================================================
console.log("=== Example 4: Field Characteristic ===\n");

/**
 * The characteristic of a field is the smallest positive n such that
 * n * 1 = 0, or 0 if no such n exists.
 *
 * For GF(p), the characteristic is p.
 */
const charP = 5n;
console.log(`Characteristic of GF(${charP}):\n`);

console.log("Adding 1 repeatedly:");
let sum = 0n;
for (let i = 1n; i <= charP + 1n; i++) {
  sum = (sum + 1n) % charP;
  console.log(`  ${i} * 1 = ${sum} mod ${charP}`);
  if (sum === 0n) {
    console.log(`\nCharacteristic = ${i}`);
    break;
  }
}

console.log("\nIn GF(p), the characteristic is always p.");
console.log("This means p * a = 0 for all elements a.");

console.log();

// =============================================================================
// Example 5: Multiplicative Group of a Field
// =============================================================================
console.log("=== Example 5: The Multiplicative Group GF(p)* ===\n");

/**
 * The non-zero elements of GF(p) form a multiplicative group GF(p)* of order p-1.
 * This group is always cyclic!
 */
const mulP = 13n;
console.log(`GF(${mulP})* has order ${mulP - 1n}\n`);

// Find a generator (primitive root)
function isPrimitiveRoot(g: bigint, p: bigint): boolean {
  const order = p - 1n;
  const factors = factor(order).filter(([f]) => f > 0n).map(([f]) => f);

  for (const q of factors) {
    if (power_mod(g, order / q, p) === 1n) {
      return false;
    }
  }
  return true;
}

let generator = 2n;
while (!isPrimitiveRoot(generator, mulP)) {
  generator++;
}

console.log(`Generator (primitive root): ${generator}`);

// Show powers
const powers: bigint[] = [];
let current = 1n;
for (let i = 0n; i < mulP - 1n; i++) {
  powers.push(current);
  current = (current * generator) % mulP;
}

console.log(`\nPowers of ${generator}:`);
for (let i = 0n; i < mulP - 1n; i++) {
  console.log(`  ${generator}^${i} = ${powers[Number(i)]}`);
}

console.log("\nGF(p)* is always cyclic of order p-1.");

console.log();

// =============================================================================
// Example 6: Solving Equations in Fields
// =============================================================================
console.log("=== Example 6: Solving Equations in Fields ===\n");

/**
 * In a field, linear equations ax = b always have a unique solution
 * (when a != 0): x = b * a^(-1)
 */
const eqP = 17n;
console.log(`Solving equations in GF(${eqP}):\n`);

const equations: [bigint, bigint][] = [
  [3n, 5n],   // 3x = 5
  [7n, 11n],  // 7x = 11
  [2n, 1n],   // 2x = 1
];

for (const [a, b] of equations) {
  const x = (b * inverse_mod(a, eqP)) % eqP;
  const verify = (a * x) % eqP;

  console.log(`  Solve: ${a}x = ${b} (mod ${eqP})`);
  console.log(`    x = ${b} * ${a}^(-1) = ${b} * ${inverse_mod(a, eqP)} = ${x}`);
  console.log(`    Verify: ${a} * ${x} = ${verify}\n`);
}

console.log("In rings with zero divisors, such equations may have no solution or many solutions!");

console.log();

// =============================================================================
// Example 7: Quadratic Residues
// =============================================================================
console.log("=== Example 7: Quadratic Residues ===\n");

/**
 * An element a in GF(p)* is a quadratic residue if a = x^2 for some x.
 * Exactly (p-1)/2 elements are quadratic residues in GF(p)*.
 */
const qrP = 11n;
console.log(`Quadratic residues in GF(${qrP}):\n`);

// Find all squares
const squares = new Map<bigint, bigint[]>();
for (let x = 1n; x < qrP; x++) {
  const sq = (x * x) % qrP;
  if (!squares.has(sq)) {
    squares.set(sq, []);
  }
  squares.get(sq)!.push(x);
}

const qr = [...squares.keys()].sort((a, b) => Number(a - b));
console.log(`Quadratic residues: {${qr.join(', ')}}`);
console.log(`Count: ${qr.length} = (${qrP}-1)/2 = ${(qrP - 1n) / 2n}\n`);

console.log("Square roots:");
for (const [sq, roots] of [...squares].sort((a, b) => Number(a[0] - b[0]))) {
  console.log(`  sqrt(${sq}) = +/-${roots[0]} = {${roots.join(', ')}}`);
}

// Non-residues
const nonResidues: bigint[] = [];
for (let a = 1n; a < qrP; a++) {
  if (!squares.has(a)) {
    nonResidues.push(a);
  }
}
console.log(`\nQuadratic non-residues: {${nonResidues.join(', ')}}`);

console.log();

// =============================================================================
// Example 8: Frobenius Automorphism
// =============================================================================
console.log("=== Example 8: Frobenius Map ===\n");

/**
 * In GF(p), the Frobenius map is x -> x^p.
 * By Fermat's Little Theorem, x^p = x for all x in GF(p).
 * So Frobenius is the identity on prime fields!
 */
const frobP = 5n;
console.log(`Frobenius map in GF(${frobP}): x -> x^${frobP}\n`);

console.log("Verifying x^p = x for all x:");
for (let x = 0n; x < frobP; x++) {
  const frobenius = power_mod(x, frobP, frobP);
  console.log(`  ${x}^${frobP} = ${frobenius}  (equals ${x}? ${frobenius === x})`);
}

console.log("\nFrobenius is the identity on GF(p), but is non-trivial on extension fields!");

console.log();

// =============================================================================
// Example 9: Field Extensions Preview
// =============================================================================
console.log("=== Example 9: Preview of Field Extensions ===\n");

/**
 * GF(p) is the smallest field of characteristic p.
 * Larger finite fields GF(p^n) are built as extensions.
 */
console.log("Finite fields of characteristic 2:\n");
console.log("  GF(2)   = {0, 1}                     - 2 elements");
console.log("  GF(4)   = GF(2)[x]/(x^2 + x + 1)    - 4 elements");
console.log("  GF(8)   = GF(2)[x]/(x^3 + x + 1)    - 8 elements");
console.log("  GF(16)  = GF(2)[x]/(x^4 + x + 1)    - 16 elements");
console.log("  GF(256) = GF(2)[x]/(x^8 + ...)      - 256 elements (used in AES!)\n");

console.log("Finite fields of characteristic 3:\n");
console.log("  GF(3)   = {0, 1, 2}                  - 3 elements");
console.log("  GF(9)   = GF(3)[x]/(x^2 + 1)        - 9 elements");
console.log("  GF(27)  = GF(3)[x]/(x^3 + ...)      - 27 elements\n");

console.log("For any prime p and positive integer n, there is a unique field GF(p^n)!");

console.log();

// =============================================================================
// Exercise: Verify Field Properties
// =============================================================================
console.log("=== Exercise: Verify GF(11) is a Field ===\n");

/**
 * EXERCISE: Verify all field axioms for GF(11).
 */
const verifyP = 11n;
console.log(`Verifying GF(${verifyP}) is a field:\n`);

// 1. Non-zero elements are invertible
console.log("1. Every non-zero element has an inverse:");
let allInvertible = true;
for (let a = 1n; a < verifyP; a++) {
  const inv = inverse_mod(a, verifyP);
  const product = (a * inv) % verifyP;
  if (product !== 1n) allInvertible = false;
  console.log(`   ${a.toString().padStart(2)}^(-1) = ${inv.toString().padStart(2)}  (product = ${product})`);
}
console.log(`   All invertible? ${allInvertible}\n`);

// 2. No zero divisors
console.log("2. No zero divisors:");
let hasZeroDivisors = false;
for (let a = 1n; a < verifyP && !hasZeroDivisors; a++) {
  for (let b = 1n; b < verifyP; b++) {
    if ((a * b) % verifyP === 0n) {
      hasZeroDivisors = true;
      console.log(`   Found: ${a} * ${b} = 0`);
      break;
    }
  }
}
if (!hasZeroDivisors) {
  console.log("   No zero divisors found!");
}

console.log(`\nGF(${verifyP}) is a field: ${allInvertible && !hasZeroDivisors}`);

console.log();

// =============================================================================
// Cryptographic Application: ECDSA Prime Field
// =============================================================================
console.log("=== Cryptographic Application: ECDSA and Prime Fields ===\n");

console.log("Elliptic Curve Cryptography uses prime fields:\n");

console.log("secp256k1 (Bitcoin's curve):");
console.log("  Field: GF(p) where p = 2^256 - 2^32 - 977");
console.log("  This is a prime field with ~2^256 elements");
console.log("  All field operations (add, multiply, divide) are defined\n");

console.log("P-256 (NIST curve):");
console.log("  Field: GF(p) where p = 2^256 - 2^224 + 2^192 + 2^96 - 1");
console.log("  Also a prime field\n");

console.log("Why prime fields for ECC?");
console.log("  1. Simple, well-understood arithmetic");
console.log("  2. Efficient implementations");
console.log("  3. Division (point subtraction) is always defined");
console.log("  4. Security proofs are cleaner\n");

// Show a mini example
const eccP = 23n;
console.log(`Mini example: Elliptic curve over GF(${eccP})`);
console.log(`  Points on y^2 = x^3 + ax + b require field operations:`);
console.log(`  - Addition and subtraction (mod ${eccP})`);
console.log(`  - Multiplication (mod ${eccP})`);
console.log(`  - Division (multiply by inverse) for slope calculations`);
