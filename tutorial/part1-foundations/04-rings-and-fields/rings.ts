/**
 * # Rings: The Algebraic Structure Z/nZ
 *
 * A ring is an algebraic structure with two operations: addition and multiplication.
 * The ring Z/nZ (integers modulo n) is fundamental to cryptography.
 *
 * ## Mathematical Background
 *
 * A RING (R, +, *) is a set R with two binary operations satisfying:
 *   1. (R, +) is an abelian group (with identity 0)
 *   2. (R, *) is a monoid (associative with identity 1)
 *   3. Multiplication distributes over addition
 *
 * Key concepts in rings:
 *   - ZERO DIVISORS: Elements a, b != 0 where ab = 0
 *   - UNITS: Elements with multiplicative inverses
 *   - INTEGRAL DOMAIN: A ring with no zero divisors
 *   - FIELD: A ring where every non-zero element is a unit
 *
 * ## Why This Matters for Cryptography
 *
 *   - RSA operates in the ring Z/nZ where n = pq
 *   - Zero divisors in Z/nZ reveal the factorization of n!
 *   - Polynomial rings are used in lattice cryptography
 *   - AES uses the ring GF(2^8)
 */

import { gcd, factor, euler_phi, Mod, is_prime, divisors, Zmod, IntegerModRing } from 'sagemath-ts';

// =============================================================================
// Example 1: Ring Z/nZ Basics
// =============================================================================
console.log("=== Example 1: The Ring Z/nZ ===\n");

/**
 * Z/nZ is a ring with both addition and multiplication.
 * The Zmod() function creates this ring.
 */
const Z12 = Zmod(12);

console.log(`The ring Z/12Z:`);
console.log(`  Cardinality: ${Z12.cardinality()}`);
console.log(`  Characteristic: ${Z12.characteristic}`);
console.log(`  Is a field? ${Z12.is_field()}\n`);

// Demonstrate both operations
const a = Z12.__call__(7);
const b = Z12.__call__(5);

console.log("Ring operations:");
console.log(`  ${a.value} + ${b.value} = ${a.add(b).value} (mod 12)`);
console.log(`  ${a.value} - ${b.value} = ${a.sub(b).value} (mod 12)`);
console.log(`  ${a.value} * ${b.value} = ${a.mul(b).value} (mod 12)`);

console.log();

// =============================================================================
// Example 2: Ring Axioms Verification
// =============================================================================
console.log("=== Example 2: Verifying Ring Axioms ===\n");

/**
 * Let's verify that Z/nZ satisfies all ring axioms.
 */
const n = 6n;
console.log(`Verifying ring axioms for Z/${n}Z:\n`);

// 1. (Z/nZ, +) is an abelian group
console.log("1. (Z/nZ, +) is an abelian group:");
console.log("   - Closure: sums stay in {0,1,...,5}");
console.log("   - Associativity: (a+b)+c = a+(b+c)");
console.log("   - Identity: 0");
console.log("   - Inverses: -a = n-a");
console.log("   - Commutativity: a+b = b+a");
console.log("   CHECK!\n");

// 2. Multiplication is associative
console.log("2. Multiplication is associative:");
const x = 2n, y = 3n, z = 4n;
const left = ((x * y) % n * z) % n;
const right = (x * ((y * z) % n)) % n;
console.log(`   (${x} * ${y}) * ${z} = ${left}`);
console.log(`   ${x} * (${y} * ${z}) = ${right}`);
console.log(`   Equal? ${left === right}. CHECK!\n`);

// 3. Multiplicative identity
console.log("3. Multiplicative identity exists:");
console.log("   1 is the identity: 1 * a = a * 1 = a for all a");
console.log("   CHECK!\n");

// 4. Distributive laws
console.log("4. Distributive laws:");
const dist_a = 2n, dist_b = 3n, dist_c = 4n;
const distLeft = (dist_a * ((dist_b + dist_c) % n)) % n;
const distRight = ((dist_a * dist_b) % n + (dist_a * dist_c) % n) % n;
console.log(`   ${dist_a} * (${dist_b} + ${dist_c}) = ${distLeft}`);
console.log(`   ${dist_a}*${dist_b} + ${dist_a}*${dist_c} = ${distRight}`);
console.log(`   Equal? ${distLeft === distRight}. CHECK!`);

console.log();

// =============================================================================
// Example 3: Zero Divisors
// =============================================================================
console.log("=== Example 3: Zero Divisors ===\n");

/**
 * A zero divisor is a non-zero element a where ab = 0 for some non-zero b.
 * Zero divisors exist when n is composite.
 */
const zdMod = 12n;
console.log(`Finding zero divisors in Z/${zdMod}Z:\n`);

const zeroDivisors = new Set<bigint>();

for (let a = 1n; a < zdMod; a++) {
  for (let b = 1n; b < zdMod; b++) {
    if ((a * b) % zdMod === 0n) {
      zeroDivisors.add(a);
      zeroDivisors.add(b);
    }
  }
}

console.log(`Zero divisors: {${[...zeroDivisors].sort((a, b) => Number(a - b)).join(', ')}}`);

// Show specific examples
console.log("\nExamples of zero products:");
for (let a = 1n; a < zdMod; a++) {
  for (let b = a; b < zdMod; b++) {
    if ((a * b) % zdMod === 0n) {
      console.log(`  ${a} * ${b} = ${a * b} = 0 (mod ${zdMod})`);
    }
  }
}

// Connect to factorization
console.log(`\nWhy? ${zdMod} = ${factor(zdMod).map(([p, e]) => e === 1n ? `${p}` : `${p}^${e}`).join(' * ')}`);
console.log("Elements sharing factors with n are zero divisors or zero.");

console.log();

// =============================================================================
// Example 4: Units vs Zero Divisors
// =============================================================================
console.log("=== Example 4: Units vs Zero Divisors ===\n");

/**
 * In Z/nZ, an element is either a unit (invertible) OR a zero divisor.
 * a is a unit iff gcd(a, n) = 1.
 */
const partitionMod = 15n;
console.log(`Partitioning Z/${partitionMod}Z elements:\n`);

const units: bigint[] = [];
const zeroDivs: bigint[] = [];

for (let a = 1n; a < partitionMod; a++) {
  const g = gcd(a, partitionMod);
  if (g === 1n) {
    units.push(a);
  } else {
    zeroDivs.push(a);
  }
}

console.log(`Units (gcd = 1): {${units.join(', ')}}`);
console.log(`  Count: ${units.length} = phi(${partitionMod}) = ${euler_phi(partitionMod)}`);
console.log(`Zero divisors (gcd > 1): {${zeroDivs.join(', ')}}`);
console.log(`  Count: ${zeroDivs.length}`);

// Verify zero divisors
console.log("\nVerifying zero divisor products:");
for (const a of zeroDivs.slice(0, 3)) {
  for (const b of zeroDivs) {
    if ((a * b) % partitionMod === 0n && a <= b) {
      console.log(`  ${a} * ${b} = 0 (mod ${partitionMod})`);
    }
  }
}

console.log();

// =============================================================================
// Example 5: RSA and Zero Divisors
// =============================================================================
console.log("=== Example 5: Zero Divisors Reveal RSA Factorization ===\n");

/**
 * In RSA, if we find a zero divisor in Z/nZ, we've factored n!
 * This is the basis of many factoring attacks.
 */
const rsaN = 143n;  // 11 * 13
console.log(`RSA modulus n = ${rsaN}`);
console.log(`(Secret factorization: ${rsaN} = 11 * 13)\n`);

// Suppose we discover a zero divisor pair
const zd1 = 22n;  // = 2 * 11
const zd2 = 13n;  // = 13

console.log(`Suppose we find that ${zd1} * ${zd2} = ${zd1 * zd2} = 0 (mod ${rsaN})`);
console.log(`This means ${rsaN} divides ${zd1 * zd2}`);

// Extract factors using GCD
const factor1 = gcd(zd1, rsaN);
const factor2 = rsaN / factor1;

console.log(`\nCompute gcd(${zd1}, ${rsaN}) = ${factor1}`);
console.log(`Therefore ${rsaN} = ${factor1} * ${factor2}`);
console.log("\nRSA is BROKEN!");

console.log("\nThis is why RSA security depends on factoring being hard.");

console.log();

// =============================================================================
// Example 6: Multiplication Tables
// =============================================================================
console.log("=== Example 6: Multiplication Table of Z/nZ ===\n");

/**
 * The multiplication table reveals the ring structure.
 */
const tableMod = 6n;
console.log(`Multiplication table for Z/${tableMod}Z:\n`);

let header = "  * |";
for (let j = 0n; j < tableMod; j++) header += ` ${j}`;
console.log(header);
console.log("----+------------");

for (let i = 0n; i < tableMod; i++) {
  let row = `  ${i} |`;
  for (let j = 0n; j < tableMod; j++) {
    const prod = (i * j) % tableMod;
    row += ` ${prod}`;
  }
  console.log(row);
}

console.log("\nObservations:");
console.log("  - Row 0 and column 0 are all zeros");
console.log("  - Rows 1 and 5 contain all elements (these are units)");
console.log("  - Rows 2, 3, 4 are missing some elements");
console.log("  - The zeros in rows 2, 3, 4 (besides column 0) indicate zero divisors");

console.log();

// =============================================================================
// Example 7: Integral Domains
// =============================================================================
console.log("=== Example 7: Integral Domains ===\n");

/**
 * An integral domain is a ring with no zero divisors.
 * Z/nZ is an integral domain iff n is prime.
 */
const testModuli = [4n, 5n, 6n, 7n, 8n, 9n, 10n, 11n];

console.log("Which Z/nZ are integral domains?\n");

for (const m of testModuli) {
  const prime = is_prime(m);
  let hasZeroDivisors = false;

  // Check for zero divisors
  outer: for (let a = 1n; a < m; a++) {
    for (let b = 1n; b < m; b++) {
      if ((a * b) % m === 0n) {
        hasZeroDivisors = true;
        break outer;
      }
    }
  }

  console.log(`  Z/${m.toString().padStart(2)}Z: prime? ${prime.toString().padStart(5)}, zero divisors? ${hasZeroDivisors.toString().padStart(5)}, integral domain? ${!hasZeroDivisors}`);
}

console.log("\nZ/nZ is an integral domain iff n is prime.");

console.log();

// =============================================================================
// Example 8: Ideals in Z/nZ
// =============================================================================
console.log("=== Example 8: Ideals in Z/nZ ===\n");

/**
 * An ideal I of ring R is an additive subgroup that "absorbs" multiplication:
 * for all r in R and a in I, ra is in I.
 *
 * In Z/nZ, ideals are generated by divisors of n.
 */
const idealMod = 12n;
console.log(`Ideals of Z/${idealMod}Z:\n`);

const divisorsOfN = divisors(idealMod);
console.log(`Divisors of ${idealMod}: ${divisorsOfN.join(', ')}\n`);

for (const d of divisorsOfN) {
  // Ideal generated by n/d
  const gen = idealMod / d;
  const ideal: bigint[] = [];

  for (let k = 0n; k < d; k++) {
    ideal.push((k * gen) % idealMod);
  }

  ideal.sort((a, b) => Number(a - b));
  console.log(`  <${gen}> = {${ideal.join(', ')}} (order ${ideal.length})`);
}

console.log("\nIdeals of Z/nZ correspond to divisors of n.");
console.log("This is why Z/nZ has a unique factorization structure.");

console.log();

// =============================================================================
// Example 9: Quotient Rings
// =============================================================================
console.log("=== Example 9: Quotient Rings ===\n");

/**
 * Z/nZ is itself a quotient ring: Z/nZ = Z / (n).
 * We can also form quotients of Z/nZ by its ideals.
 */
const qMod = 12n;
const idealGen = 3n;  // Ideal <3> = {0, 3, 6, 9}

console.log(`Quotient ring (Z/${qMod}Z) / <${idealGen}>:\n`);

// Find cosets
const ideal: bigint[] = [];
for (let k = 0n; (k * idealGen) % qMod !== 0n || k === 0n; k++) {
  ideal.push((k * idealGen) % qMod);
  if (ideal.length > Number(qMod)) break;
}
ideal.sort((a, b) => Number(a - b));

console.log(`Ideal <${idealGen}> = {${ideal.join(', ')}}`);

const cosets: { rep: bigint; elements: bigint[] }[] = [];
const seen = new Set<bigint>();

for (let r = 0n; r < qMod; r++) {
  if (!seen.has(r)) {
    const coset = ideal.map(i => (r + i) % qMod).sort((a, b) => Number(a - b));
    cosets.push({ rep: r, elements: coset });
    for (const e of coset) seen.add(e);
  }
}

console.log(`\nCosets (elements of quotient ring):`);
for (const { rep, elements } of cosets) {
  console.log(`  ${rep} + <${idealGen}> = {${elements.join(', ')}}`);
}

console.log(`\nQuotient has ${cosets.length} elements.`);
console.log(`(Z/${qMod}Z) / <${idealGen}> is isomorphic to Z/${cosets.length}Z`);

console.log();

// =============================================================================
// Exercise: Find All Zero Divisors
// =============================================================================
console.log("=== Exercise: Analyze Zero Divisors ===\n");

/**
 * EXERCISE: Find all zero divisor pairs in Z/30Z and relate them to factorization.
 */
const exMod = 30n;
const factorization = factor(exMod);

console.log(`Analyzing Z/${exMod}Z:`);
console.log(`  ${exMod} = ${factorization.map(([p, e]) => e === 1n ? `${p}` : `${p}^${e}`).join(' * ')}\n`);

// Find all distinct zero divisor pairs
const zdPairs: [bigint, bigint][] = [];
for (let a = 1n; a < exMod; a++) {
  for (let b = a; b < exMod; b++) {
    if ((a * b) % exMod === 0n) {
      zdPairs.push([a, b]);
    }
  }
}

console.log(`Zero divisor pairs (a * b = 0 mod ${exMod}):`);
for (const [a, b] of zdPairs.slice(0, 10)) {
  const gcdA = gcd(a, exMod);
  const gcdB = gcd(b, exMod);
  console.log(`  ${a.toString().padStart(2)} * ${b.toString().padStart(2)} = 0  [gcd(${a},${exMod})=${gcdA}, gcd(${b},${exMod})=${gcdB}]`);
}
if (zdPairs.length > 10) {
  console.log(`  ... (${zdPairs.length - 10} more pairs)`);
}

console.log(`\nTotal zero divisor pairs: ${zdPairs.length}`);
console.log("Zero divisors are elements sharing a prime factor with n.");

console.log();

// =============================================================================
// Cryptographic Application: Ring Structure in RSA
// =============================================================================
console.log("=== Cryptographic Application: RSA Ring Structure ===\n");

console.log("RSA operates in the ring Z/nZ where n = p * q\n");

const p = 11n, q = 13n;
const rsaModulus = p * q;

console.log(`Example: n = ${p} * ${q} = ${rsaModulus}\n`);

console.log("Ring structure of Z/143Z:");
console.log(`  - Total elements: ${rsaModulus}`);
console.log(`  - Units: phi(${rsaModulus}) = ${euler_phi(rsaModulus)}`);
console.log(`  - Zero divisors: ${rsaModulus - 1n - euler_phi(rsaModulus)}`);

console.log(`\nSpecifically:`);
console.log(`  - Multiples of ${p}: 0, ${p}, ${2n*p}, ... (${q} elements)`);
console.log(`  - Multiples of ${q}: 0, ${q}, ${2n*q}, ... (${p} elements)`);
console.log(`  - Their intersection: just {0}`);
console.log(`  - Zero divisors = ${p + q - 2n} elements (excluding 0)`);

console.log(`\nKey insight: Finding ANY zero divisor reveals the factorization!`);
console.log(`This is why RSA security = difficulty of finding zero divisors in Z/nZ.`);
