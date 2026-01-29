/**
 * # The Multiplicative Group (Z/nZ)*
 *
 * The multiplicative group (Z/nZ)* consists of all elements in Z/nZ that
 * have multiplicative inverses. This group is fundamental to RSA and
 * Diffie-Hellman cryptography.
 *
 * ## Mathematical Background
 *
 * (Z/nZ)* = {a in Z/nZ : gcd(a, n) = 1}
 *
 * Key properties:
 *   - |(Z/nZ)*| = phi(n) (Euler's totient function)
 *   - (Z/nZ)* is a group under multiplication
 *   - For prime p: (Z/pZ)* is cyclic of order p-1
 *   - For n = 2, 4, p^k, 2p^k (p odd prime): (Z/nZ)* is cyclic
 *   - For other n: (Z/nZ)* is a direct product of cyclic groups
 *
 * ## Why This Matters for Cryptography
 *
 *   - RSA: Security relies on properties of (Z/nZ)* where n = pq
 *   - Diffie-Hellman: Uses the cyclic group (Z/pZ)* for safe prime p
 *   - ElGamal: Encryption in (Z/pZ)*
 *   - Schnorr signatures: Work in subgroups of (Z/pZ)*
 */

import { gcd, euler_phi, power_mod, factor, divisors, lcm, inverse_mod, Mod, Zmod } from 'sagemath-ts';

// =============================================================================
// Example 1: Building (Z/nZ)*
// =============================================================================
console.log("=== Example 1: Building (Z/nZ)* ===\n");

/**
 * (Z/nZ)* contains exactly those elements coprime to n.
 */
function buildMultiplicativeGroup(n: bigint): bigint[] {
  const units: bigint[] = [];
  for (let a = 1n; a < n; a++) {
    if (gcd(a, n) === 1n) {
      units.push(a);
    }
  }
  return units;
}

const testModuli = [8n, 12n, 15n, 16n];

for (const n of testModuli) {
  const units = buildMultiplicativeGroup(n);
  const phi = euler_phi(n);

  console.log(`(Z/${n}Z)* = {${units.join(', ')}}`);
  console.log(`  Order: ${units.length} = phi(${n}) = ${phi}\n`);
}

// =============================================================================
// Example 2: Group Operation in (Z/nZ)*
// =============================================================================
console.log("=== Example 2: Group Operation ===\n");

/**
 * The group operation is multiplication modulo n.
 * Let's verify the group axioms.
 */
const n = 10n;
const units = buildMultiplicativeGroup(n);

console.log(`Verifying (Z/${n}Z)* = {${units.join(', ')}} is a group:\n`);

// Closure
console.log("1. CLOSURE:");
let closure = true;
for (const a of units) {
  for (const b of units) {
    const product = (a * b) % n;
    if (!units.includes(product)) {
      closure = false;
      console.log(`   FAIL: ${a} * ${b} = ${product} not in group`);
    }
  }
}
if (closure) console.log("   All products are in the group. CHECK!");

// Identity
console.log("\n2. IDENTITY: 1 is the identity.");
for (const a of units) {
  if ((1n * a) % n !== a || (a * 1n) % n !== a) {
    console.log(`   FAIL for ${a}`);
  }
}
console.log("   1 * a = a * 1 = a for all a. CHECK!");

// Inverses
console.log("\n3. INVERSES:");
for (const a of units) {
  const inv = inverse_mod(a, n);
  console.log(`   ${a}^(-1) = ${inv}  (${a} * ${inv} = ${(a * inv) % n})`);
}

// Associativity
console.log("\n4. ASSOCIATIVITY: Multiplication is associative. CHECK!");

console.log();

// =============================================================================
// Example 3: Structure of (Z/nZ)* for Various n
// =============================================================================
console.log("=== Example 3: Structure of (Z/nZ)* ===\n");

/**
 * The structure depends on the factorization of n.
 * By the Chinese Remainder Theorem:
 *   (Z/nZ)* = (Z/p1^a1)* x (Z/p2^a2)* x ...
 */
function analyzeGroupStructure(n: bigint): void {
  const phi = euler_phi(n);
  const factors = factor(n).filter(([p]) => p > 0n);
  const units = buildMultiplicativeGroup(n);

  console.log(`(Z/${n}Z)*:`);
  console.log(`  Factorization: ${n} = ${factors.map(([p, e]) => e === 1n ? `${p}` : `${p}^${e}`).join(' * ')}`);
  console.log(`  Order: phi(${n}) = ${phi}`);

  // Find all element orders
  const orders = new Map<bigint, number>();
  let maxOrder = 0n;

  for (const a of units) {
    const ord = computeOrder(a, n);
    orders.set(ord, (orders.get(ord) || 0) + 1);
    if (ord > maxOrder) maxOrder = ord;
  }

  console.log(`  Element orders: ${[...orders].map(([o, c]) => `ord ${o}: ${c} elem`).join(', ')}`);
  console.log(`  Max order (exponent): ${maxOrder}`);
  console.log(`  Is cyclic? ${maxOrder === phi}`);
  console.log();
}

function computeOrder(a: bigint, n: bigint): bigint {
  let order = 1n;
  let current = a % n;
  while (current !== 1n) {
    current = (current * a) % n;
    order++;
  }
  return order;
}

const analyzeModuli = [5n, 7n, 8n, 9n, 12n, 15n];
for (const m of analyzeModuli) {
  analyzeGroupStructure(m);
}

// =============================================================================
// Example 4: Cyclic vs Non-Cyclic
// =============================================================================
console.log("=== Example 4: When is (Z/nZ)* Cyclic? ===\n");

/**
 * (Z/nZ)* is cyclic iff n = 1, 2, 4, p^k, or 2p^k for odd prime p.
 */
function isCyclic(n: bigint): boolean {
  if (n === 1n || n === 2n || n === 4n) return true;

  const factors = factor(n).filter(([p]) => p > 0n);

  // n = p^k for odd prime p
  if (factors.length === 1) {
    const [prime, _exp] = factors[0]!;
    return prime !== 2n;
  }

  // n = 2 * p^k for odd prime p
  if (factors.length === 2) {
    const [[p1, e1], [p2, _e2]] = factors;
    if (p1 === 2n && e1 === 1n && p2! > 2n) return true;
  }

  return false;
}

console.log("Testing which (Z/nZ)* are cyclic:\n");

for (let n = 2n; n <= 20n; n++) {
  const cyclic = isCyclic(n);
  const phi = euler_phi(n);
  const units = buildMultiplicativeGroup(n);
  const maxOrd = units.reduce((max, a) => {
    const ord = computeOrder(a, n);
    return ord > max ? ord : max;
  }, 0n);

  console.log(`  n=${n.toString().padStart(2)}: phi=${phi.toString().padStart(2)}, max_ord=${maxOrd.toString().padStart(2)}, cyclic=${cyclic}`);
}

console.log();

// =============================================================================
// Example 5: The Klein Four-Group (Z/8Z)*
// =============================================================================
console.log("=== Example 5: (Z/8Z)* - The Klein Four-Group ===\n");

/**
 * (Z/8Z)* = {1, 3, 5, 7} is isomorphic to Z/2Z x Z/2Z, not Z/4Z.
 * Every non-identity element has order 2!
 */
const mod8 = 8n;
const units8 = buildMultiplicativeGroup(mod8);

console.log(`(Z/${mod8}Z)* = {${units8.join(', ')}}\n`);

console.log("Multiplication table:");
let header = "  * |";
for (const u of units8) header += ` ${u}`;
console.log(header);
console.log("----+--------");
for (const a of units8) {
  let row = `  ${a} |`;
  for (const b of units8) {
    row += ` ${(a * b) % mod8}`;
  }
  console.log(row);
}

console.log("\nElement orders:");
for (const a of units8) {
  const ord = computeOrder(a, mod8);
  console.log(`  ord(${a}) = ${ord}`);
}

console.log("\nAll non-identity elements have order 2!");
console.log("This is the Klein four-group V_4 = Z/2Z x Z/2Z.");

console.log();

// =============================================================================
// Example 6: Primitive Roots
// =============================================================================
console.log("=== Example 6: Primitive Roots (Generators) ===\n");

/**
 * A primitive root modulo n is a generator of (Z/nZ)*.
 * These exist only when (Z/nZ)* is cyclic.
 */
function findPrimitiveRoots(n: bigint): bigint[] {
  const phi = euler_phi(n);
  const units = buildMultiplicativeGroup(n);
  const primitiveRoots: bigint[] = [];

  for (const g of units) {
    if (computeOrder(g, n) === phi) {
      primitiveRoots.push(g);
    }
  }

  return primitiveRoots;
}

const primRootModuli = [7n, 9n, 11n, 13n, 14n];

for (const n of primRootModuli) {
  const phi = euler_phi(n);
  const roots = findPrimitiveRoots(n);

  console.log(`(Z/${n}Z)*: phi = ${phi}`);
  if (roots.length > 0) {
    console.log(`  Primitive roots: {${roots.join(', ')}}`);
    console.log(`  Count: ${roots.length} = phi(${phi}) = ${euler_phi(phi)}`);
  } else {
    console.log(`  No primitive roots (group is not cyclic)`);
  }
  console.log();
}

// =============================================================================
// Example 7: Euler's Theorem
// =============================================================================
console.log("=== Example 7: Euler's Theorem ===\n");

/**
 * Euler's Theorem: For gcd(a, n) = 1, a^phi(n) = 1 (mod n)
 * This is because phi(n) is the order of (Z/nZ)*.
 */
const eulerN = 15n;
const eulerPhi = euler_phi(eulerN);

console.log(`Euler's Theorem for n = ${eulerN}:`);
console.log(`  phi(${eulerN}) = ${eulerPhi}\n`);

console.log(`For all a coprime to ${eulerN}, a^${eulerPhi} = 1 mod ${eulerN}:`);
const eulerUnits = buildMultiplicativeGroup(eulerN);
for (const a of eulerUnits) {
  const result = power_mod(a, eulerPhi, eulerN);
  console.log(`  ${a.toString().padStart(2)}^${eulerPhi} = ${result}`);
}

console.log("\nThis is the foundation of RSA!");

console.log();

// =============================================================================
// Example 8: Carmichael's Lambda Function
// =============================================================================
console.log("=== Example 8: Carmichael's Lambda Function ===\n");

/**
 * The Carmichael function lambda(n) is the exponent of (Z/nZ)*:
 * the smallest m such that a^m = 1 for all a in (Z/nZ)*.
 *
 * For cyclic groups: lambda(n) = phi(n)
 * For non-cyclic groups: lambda(n) < phi(n)
 */
function carmichaelLambda(n: bigint): bigint {
  const units = buildMultiplicativeGroup(n);
  const orders = units.map(a => computeOrder(a, n));
  return lcm(orders);
}

console.log("Comparing phi(n) and lambda(n):\n");

for (let n = 2n; n <= 20n; n++) {
  const phi = euler_phi(n);
  const lambda = carmichaelLambda(n);
  const ratio = phi / lambda;

  console.log(`  n=${n.toString().padStart(2)}: phi=${phi.toString().padStart(2)}, lambda=${lambda.toString().padStart(2)}, phi/lambda=${ratio}`);
}

console.log("\nlambda(n) always divides phi(n).");
console.log("For RSA, using lambda(n) instead of phi(n) gives a smaller private exponent d.");

console.log();

// =============================================================================
// Example 9: Subgroup Structure
// =============================================================================
console.log("=== Example 9: Subgroups of (Z/nZ)* ===\n");

/**
 * For prime p, (Z/pZ)* has exactly one subgroup for each divisor of p-1.
 * This is because (Z/pZ)* is cyclic.
 */
const subgroupP = 13n;
const subgroupOrder = subgroupP - 1n;  // = 12

console.log(`Subgroups of (Z/${subgroupP}Z)*:`);
console.log(`  Group order: ${subgroupOrder} = 2^2 * 3`);
console.log(`  Divisors: ${divisors(subgroupOrder).join(', ')}\n`);

// Find a generator
const primitiveRoot = findPrimitiveRoots(subgroupP)[0]!;
console.log(`Primitive root: ${primitiveRoot}\n`);

console.log("Subgroups (generated by g^(n/d) for each divisor d):");
for (const d of divisors(subgroupOrder)) {
  const gen = power_mod(primitiveRoot, subgroupOrder / d, subgroupP);

  // Generate the subgroup
  const subgroup: bigint[] = [];
  let current = 1n;
  do {
    subgroup.push(current);
    current = (current * gen) % subgroupP;
  } while (current !== 1n);

  subgroup.sort((a, b) => Number(a - b));
  console.log(`  Order ${d.toString().padStart(2)}: <${gen}> = {${subgroup.join(', ')}}`);
}

console.log();

// =============================================================================
// Exercise: Compute (Z/nZ)* Structure
// =============================================================================
console.log("=== Exercise: Analyze (Z/24Z)* ===\n");

/**
 * EXERCISE: Analyze the structure of (Z/24Z)*.
 */
const exerciseN = 24n;
const exercisePhi = euler_phi(exerciseN);
const exerciseUnits = buildMultiplicativeGroup(exerciseN);

console.log(`Analyzing (Z/${exerciseN}Z)*:`);
console.log(`  ${exerciseN} = ${factor(exerciseN).filter(([p]) => p > 0n).map(([p, e]) => e === 1n ? `${p}` : `${p}^${e}`).join(' * ')}`);
console.log(`  phi(${exerciseN}) = ${exercisePhi}`);
console.log(`  (Z/${exerciseN}Z)* = {${exerciseUnits.join(', ')}}\n`);

// Order distribution
const orderDist = new Map<bigint, bigint[]>();
for (const a of exerciseUnits) {
  const ord = computeOrder(a, exerciseN);
  if (!orderDist.has(ord)) {
    orderDist.set(ord, []);
  }
  orderDist.get(ord)!.push(a);
}

console.log("Elements by order:");
for (const [ord, elems] of [...orderDist].sort((a, b) => Number(a[0] - b[0]))) {
  console.log(`  Order ${ord}: {${elems.join(', ')}}`);
}

const maxOrder = Math.max(...[...orderDist.keys()].map(Number));
console.log(`\nMaximum order: ${maxOrder}`);
console.log(`Is cyclic? ${BigInt(maxOrder) === exercisePhi}`);
console.log(`Lambda(${exerciseN}) = ${carmichaelLambda(exerciseN)}`);

console.log();

// =============================================================================
// Cryptographic Application: RSA Key Size
// =============================================================================
console.log("=== Cryptographic Application: RSA and (Z/nZ)* ===\n");

console.log("RSA operates in (Z/nZ)* where n = p * q:\n");

const rsaP = 11n;
const rsaQ = 13n;
const rsaN = rsaP * rsaQ;
const rsaPhi = (rsaP - 1n) * (rsaQ - 1n);
const rsaLambda = lcm([rsaP - 1n, rsaQ - 1n]);

console.log(`Example: p = ${rsaP}, q = ${rsaQ}, n = ${rsaN}`);
console.log(`  phi(n) = ${rsaPhi}`);
console.log(`  lambda(n) = ${rsaLambda}\n`);

console.log("Key observations:");
console.log(`  1. |(Z/nZ)*| = phi(n) = ${rsaPhi}`);
console.log(`  2. (Z/nZ)* = (Z/pZ)* x (Z/qZ)* (by CRT)`);
console.log(`  3. For encryption: c = m^e mod n`);
console.log(`  4. For decryption: m = c^d mod n where ed = 1 mod lambda(n)`);

// Compute keys
const rsaE = 7n;
const rsaD_phi = inverse_mod(rsaE, rsaPhi);
const rsaD_lambda = inverse_mod(rsaE, rsaLambda);

console.log(`\nPublic exponent: e = ${rsaE}`);
console.log(`Private exponent (using phi): d = ${rsaD_phi}`);
console.log(`Private exponent (using lambda): d = ${rsaD_lambda}`);

console.log("\nUsing lambda gives a smaller d, which can speed up decryption!");
