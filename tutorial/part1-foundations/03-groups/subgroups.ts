/**
 * # Subgroups and Cosets
 *
 * Subgroups are groups within groups. Understanding subgroups and cosets
 * is essential for understanding group structure and for cryptographic
 * applications like small-subgroup attacks.
 *
 * ## Mathematical Background
 *
 * A SUBGROUP H of G is a subset of G that is itself a group under the same operation.
 * Equivalently, H is a subgroup iff:
 *   1. The identity is in H
 *   2. H is closed under the group operation
 *   3. H is closed under taking inverses
 *
 * A LEFT COSET of H in G is a set gH = {gh : h in H} for some g in G.
 * A RIGHT COSET is Hg = {hg : h in H}.
 *
 * Key facts:
 *   - All cosets have the same size |H|
 *   - Cosets partition G: every element is in exactly one coset
 *   - |G| = |H| * [G : H] where [G : H] is the number of cosets (index)
 *
 * ## Why This Matters for Cryptography
 *
 *   - Small-subgroup attacks exploit subgroups of small order
 *   - Safe primes: p = 2q + 1 where q is also prime (few subgroups)
 *   - Cofactor multiplication in ECC to avoid subgroup attacks
 *   - Quotient groups and normal subgroups in advanced protocols
 */

import { gcd, euler_phi, power_mod, factor, divisors, lcm, Mod, Zmod } from 'sagemath-ts';

// =============================================================================
// Example 1: Finding Subgroups
// =============================================================================
console.log("=== Example 1: Finding Subgroups ===\n");

/**
 * Let's find all subgroups of (Z/12Z, +).
 * By Lagrange's theorem, subgroup orders must divide 12.
 */
const n = 12n;
console.log(`Finding subgroups of (Z/${n}Z, +):\n`);

// Find subgroups by finding what each element generates
function findAdditiveSubgroup(g: bigint, n: bigint): bigint[] {
  const subgroup: bigint[] = [];
  let current = 0n;
  do {
    subgroup.push(current);
    current = (current + g) % n;
  } while (current !== 0n);
  return subgroup.sort((a, b) => Number(a - b));
}

const subgroups = new Map<string, bigint[]>();

for (let g = 0n; g < n; g++) {
  const H = findAdditiveSubgroup(g, n);
  const key = H.join(',');
  if (!subgroups.has(key)) {
    subgroups.set(key, H);
    console.log(`  <${g}> = {${H.join(', ')}} (order ${H.length})`);
  }
}

console.log(`\nTotal: ${subgroups.size} distinct subgroups`);
console.log(`Divisors of ${n}: ${divisors(n).join(', ')}`);
console.log("One subgroup for each divisor!");

console.log();

// =============================================================================
// Example 2: Subgroups of Multiplicative Groups
// =============================================================================
console.log("=== Example 2: Subgroups of (Z/pZ)* ===\n");

/**
 * For prime p, (Z/pZ)* is cyclic of order p-1.
 * Subgroups have orders that divide p-1.
 */
const p = 13n;
const groupOrder = p - 1n;

console.log(`(Z/${p}Z)* has order ${groupOrder}`);
console.log(`Divisors of ${groupOrder}: ${divisors(groupOrder).join(', ')}\n`);

function findMultiplicativeSubgroup(g: bigint, p: bigint): bigint[] {
  const subgroup: bigint[] = [];
  let current = 1n;
  do {
    subgroup.push(current);
    current = (current * g) % p;
  } while (current !== 1n);
  return subgroup.sort((a, b) => Number(a - b));
}

const multSubgroups = new Map<string, { generator: bigint; elements: bigint[] }>();

for (let g = 1n; g < p; g++) {
  const H = findMultiplicativeSubgroup(g, p);
  const key = H.join(',');
  if (!multSubgroups.has(key)) {
    multSubgroups.set(key, { generator: g, elements: H });
  }
}

console.log("Subgroups of (Z/13Z)*:");
for (const { generator, elements } of multSubgroups.values()) {
  console.log(`  <${generator}> = {${elements.join(', ')}} (order ${elements.length})`);
}

console.log();

// =============================================================================
// Example 3: Cosets
// =============================================================================
console.log("=== Example 3: Cosets ===\n");

/**
 * A coset gH is the set of all products gh for h in H.
 * Cosets partition the group.
 */
const cosetMod = 12n;
const H = findAdditiveSubgroup(4n, cosetMod);  // H = {0, 4, 8}

console.log(`In (Z/${cosetMod}Z, +):`);
console.log(`  Subgroup H = <4> = {${H.join(', ')}}  (order ${H.length})\n`);

console.log("Cosets of H:");
const cosets = new Map<string, bigint>();

for (let g = 0n; g < cosetMod; g++) {
  // Compute g + H
  const coset = H.map(h => (g + h) % cosetMod).sort((a, b) => Number(a - b));
  const key = coset.join(',');

  if (!cosets.has(key)) {
    cosets.set(key, g);
    console.log(`  ${g} + H = {${coset.join(', ')}}`);
  }
}

console.log(`\nNumber of cosets = ${cosets.size} = ${cosetMod}/${BigInt(H.length)} = |G|/|H|`);
console.log("The cosets partition Z/12Z!");

console.log();

// =============================================================================
// Example 4: Lagrange's Theorem via Cosets
// =============================================================================
console.log("=== Example 4: Lagrange's Theorem via Cosets ===\n");

/**
 * Proof sketch of Lagrange's Theorem:
 * 1. All cosets have the same size as H
 * 2. Cosets partition G (every element is in exactly one coset)
 * 3. Therefore |G| = (number of cosets) * |H|
 * 4. So |H| divides |G|
 */
const lagrangeMod = 15n;
const lagrangeH = findAdditiveSubgroup(5n, lagrangeMod);  // H = {0, 5, 10}

console.log(`Demonstrating Lagrange's Theorem:`);
console.log(`  G = (Z/${lagrangeMod}Z, +), |G| = ${lagrangeMod}`);
console.log(`  H = <5> = {${lagrangeH.join(', ')}}, |H| = ${lagrangeH.length}\n`);

// Find all cosets
const lagrangeCosets: bigint[][] = [];
const seenElements = new Set<bigint>();

for (let g = 0n; g < lagrangeMod; g++) {
  if (!seenElements.has(g)) {
    const coset = lagrangeH.map(h => (g + h) % lagrangeMod).sort((a, b) => Number(a - b));
    lagrangeCosets.push(coset);
    for (const elem of coset) {
      seenElements.add(elem);
    }
    console.log(`  Coset ${lagrangeCosets.length}: ${g} + H = {${coset.join(', ')}}`);
  }
}

console.log(`\nVerification:`);
console.log(`  |G| = ${lagrangeMod}`);
console.log(`  |H| = ${lagrangeH.length}`);
console.log(`  Number of cosets [G:H] = ${lagrangeCosets.length}`);
console.log(`  |G| = |H| * [G:H]? ${lagrangeMod} = ${lagrangeH.length} * ${lagrangeCosets.length} = ${BigInt(lagrangeH.length * lagrangeCosets.length)}`);

console.log();

// =============================================================================
// Example 5: Index of a Subgroup
// =============================================================================
console.log("=== Example 5: Index of a Subgroup ===\n");

/**
 * The index [G:H] is the number of cosets of H in G.
 * [G:H] = |G| / |H|
 */
const indexMod = 24n;
console.log(`Subgroup indices in (Z/${indexMod}Z, +):\n`);

const divisors24 = divisors(indexMod);
for (const d of divisors24) {
  // Find generator for subgroup of order d
  const gen = indexMod / d;  // <gen> has order d
  const subgroup = findAdditiveSubgroup(gen, indexMod);
  const index = indexMod / BigInt(subgroup.length);

  console.log(`  H = <${gen}>, |H| = ${subgroup.length}, [G:H] = ${indexMod}/${subgroup.length} = ${index}`);
}

console.log();

// =============================================================================
// Example 6: Normal Subgroups
// =============================================================================
console.log("=== Example 6: Normal Subgroups ===\n");

/**
 * A subgroup H is NORMAL if gH = Hg for all g in G.
 * In abelian groups, ALL subgroups are normal.
 *
 * For non-abelian groups, this is more interesting.
 */
console.log("In abelian groups (like Z/nZ), all subgroups are normal.\n");

console.log("For non-abelian groups, left and right cosets can differ.");
console.log("Example: In S_3 (symmetric group on 3 elements):");
console.log("  H = {e, (12)} is NOT normal: (13)H != H(13)");
console.log("  H = {e, (123), (132)} IS normal (it's the alternating group A_3)\n");

// Let's verify for (Z/6Z, +) that all subgroups are normal
const verifyMod = 6n;
console.log(`Verifying all subgroups of (Z/${verifyMod}Z, +) are normal:\n`);

const subsOf6 = new Map<string, bigint[]>();
for (let g = 0n; g < verifyMod; g++) {
  const H = findAdditiveSubgroup(g, verifyMod);
  const key = H.join(',');
  if (!subsOf6.has(key)) {
    subsOf6.set(key, H);
  }
}

for (const H of subsOf6.values()) {
  console.log(`  H = {${H.join(', ')}}`);
  let isNormal = true;
  for (let g = 0n; g < verifyMod; g++) {
    const leftCoset = H.map(h => (g + h) % verifyMod).sort((a, b) => Number(a - b));
    const rightCoset = H.map(h => (h + g) % verifyMod).sort((a, b) => Number(a - b));
    if (leftCoset.join(',') !== rightCoset.join(',')) {
      isNormal = false;
      break;
    }
  }
  console.log(`    Normal? ${isNormal}`);
}

console.log();

// =============================================================================
// Example 7: Quotient Groups
// =============================================================================
console.log("=== Example 7: Quotient Groups ===\n");

/**
 * If H is a normal subgroup of G, the quotient group G/H consists of
 * all cosets with operation (aH)(bH) = (ab)H.
 */
const quotientMod = 12n;
const quotientH = findAdditiveSubgroup(3n, quotientMod);  // H = {0, 3, 6, 9}

console.log(`Quotient group (Z/${quotientMod}Z) / H where H = {${quotientH.join(', ')}}:\n`);

// Find coset representatives
const cosetReps: bigint[] = [];
const seenCosets = new Set<string>();

for (let g = 0n; g < quotientMod; g++) {
  const coset = quotientH.map(h => (g + h) % quotientMod).sort((a, b) => Number(a - b));
  const key = coset.join(',');
  if (!seenCosets.has(key)) {
    seenCosets.add(key);
    cosetReps.push(g);
    console.log(`  ${g} + H = {${coset.join(', ')}}`);
  }
}

console.log(`\nQuotient group G/H has ${cosetReps.length} elements: cosets of representatives {${cosetReps.join(', ')}}`);
console.log(`G/H is isomorphic to Z/${cosetReps.length}Z`);

console.log();

// =============================================================================
// Example 8: Subgroup Structure and Cryptography
// =============================================================================
console.log("=== Example 8: Safe Primes and Subgroup Structure ===\n");

/**
 * A "safe prime" is p = 2q + 1 where q is also prime.
 * (Z/pZ)* then has only subgroups of order 1, 2, q, and 2q.
 *
 * This is important for Diffie-Hellman security.
 */
function isSafePrime(p: bigint): boolean {
  const q = (p - 1n) / 2n;
  // Check if q is prime (simple check)
  if (q < 2n) return false;
  for (let i = 2n; i * i <= q; i++) {
    if (q % i === 0n) return false;
  }
  for (let i = 2n; i * i <= p; i++) {
    if (p % i === 0n) return false;
  }
  return true;
}

const candidates = [5n, 7n, 11n, 13n, 17n, 23n, 29n, 31n, 47n, 59n];

console.log("Safe primes (p = 2q + 1 where q is prime):\n");
for (const p of candidates) {
  const q = (p - 1n) / 2n;
  const safe = isSafePrime(p);
  if (safe) {
    console.log(`  p = ${p} is safe: p - 1 = ${p - 1n} = 2 * ${q}`);
    console.log(`    Subgroup orders: 1, 2, ${q}, ${p - 1n}`);
  }
}

console.log("\nSafe primes minimize the number of subgroups, making small-subgroup attacks harder.");

console.log();

// =============================================================================
// Example 9: Small Subgroup Attack
// =============================================================================
console.log("=== Example 9: Small Subgroup Attack ===\n");

/**
 * In Diffie-Hellman, if an attacker can get your public key into a small
 * subgroup, they can determine your private key modulo the subgroup order.
 */
const attackP = 31n;  // Not a safe prime!
const attackOrder = attackP - 1n;  // = 30 = 2 * 3 * 5

console.log(`Demonstration of small-subgroup attack:`);
console.log(`  p = ${attackP}, group order = ${attackOrder} = 2 * 3 * 5\n`);

// Find an element of small order (say, order 5)
function findElementOfOrder(targetOrder: bigint, p: bigint): bigint {
  const groupOrder = p - 1n;
  if (groupOrder % targetOrder !== 0n) throw new Error("Order doesn't divide");

  for (let g = 2n; g < p; g++) {
    if (power_mod(g, groupOrder / targetOrder, p) !== 1n) {
      // g is a generator, use g^(groupOrder/targetOrder)
      return power_mod(g, groupOrder / targetOrder, p);
    }
  }
  throw new Error("Not found");
}

const smallGen = findElementOfOrder(5n, attackP);
console.log(`  Element of order 5: h = ${smallGen}`);
console.log(`  Verify: ${smallGen}^5 = ${power_mod(smallGen, 5n, attackP)} mod ${attackP}\n`);

// Alice's private key
const aliceSecret = 17n;
console.log(`  Alice's private key: a = ${aliceSecret}`);
console.log(`  Attacker sends h instead of legitimate public key`);

// Alice computes h^a
const aliceComputes = power_mod(smallGen, aliceSecret, attackP);
console.log(`  Alice computes: h^a = ${smallGen}^${aliceSecret} = ${aliceComputes}`);

// Attacker can now determine a mod 5 by brute force
console.log(`\n  Attacker's brute force:`);
for (let guess = 0n; guess < 5n; guess++) {
  const test = power_mod(smallGen, guess, attackP);
  if (test === aliceComputes) {
    console.log(`    h^${guess} = ${test} = Alice's result!`);
    console.log(`    Attacker learns: a = ${guess} mod 5`);
    console.log(`    (Actual: ${aliceSecret} mod 5 = ${aliceSecret % 5n})`);
  }
}

console.log("\nDefense: Always use safe primes or validate that received public keys are in the correct subgroup!");

console.log();

// =============================================================================
// Exercise: Find All Subgroups of (Z/20Z)*
// =============================================================================
console.log("=== Exercise: Find Subgroups of (Z/20Z)* ===\n");

/**
 * EXERCISE: Find all subgroups of (Z/20Z)*.
 * Note: This group may not be cyclic!
 */
const exMod = 20n;
const exPhi = euler_phi(exMod);

console.log(`Finding subgroups of (Z/${exMod}Z)*:`);
console.log(`  Group order = phi(${exMod}) = ${exPhi}\n`);

// Find all units
const unitsOf20: bigint[] = [];
for (let a = 1n; a < exMod; a++) {
  if (gcd(a, exMod) === 1n) {
    unitsOf20.push(a);
  }
}
console.log(`  Units: {${unitsOf20.join(', ')}}\n`);

// Find all subgroups
const multSubs = new Map<string, { generator: bigint; elements: bigint[] }>();

for (const g of unitsOf20) {
  const subgroup: bigint[] = [];
  let current = 1n;
  do {
    subgroup.push(current);
    current = (current * g) % exMod;
  } while (current !== 1n);
  subgroup.sort((a, b) => Number(a - b));

  const key = subgroup.join(',');
  if (!multSubs.has(key)) {
    multSubs.set(key, { generator: g, elements: subgroup });
  }
}

console.log("Subgroups:");
for (const { generator, elements } of [...multSubs.values()].sort((a, b) => a.elements.length - b.elements.length)) {
  console.log(`  <${generator}> = {${elements.join(', ')}} (order ${elements.length})`);
}

// Check if cyclic
const maxOrder = Math.max(...[...multSubs.values()].map(s => s.elements.length));
const isCyclic = maxOrder === Number(exPhi);
console.log(`\nMaximum element order: ${maxOrder}`);
console.log(`Group order: ${exPhi}`);
console.log(`Is (Z/${exMod}Z)* cyclic? ${isCyclic}`);
