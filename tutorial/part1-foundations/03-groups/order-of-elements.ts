/**
 * # Order of Elements and Lagrange's Theorem
 *
 * The order of an element is a fundamental concept that connects individual
 * elements to the structure of the group. Lagrange's theorem provides a
 * powerful constraint on possible orders.
 *
 * ## Mathematical Background
 *
 * The ORDER of an element g in group G is the smallest positive integer n
 * such that g^n = e (identity).
 *   ord(g) = min{n > 0 : g^n = e}
 *
 * LAGRANGE'S THEOREM: If G is a finite group and H is a subgroup of G,
 * then |H| divides |G|.
 *
 * Corollary: The order of any element divides the order of the group.
 * Corollary: g^|G| = e for all g in G.
 *
 * ## Why This Matters for Cryptography
 *
 *   - Fermat's Little Theorem: a^(p-1) = 1 mod p (p prime)
 *   - Euler's Theorem: a^phi(n) = 1 mod n when gcd(a,n) = 1
 *   - RSA: Works because d*e = 1 mod phi(n)
 *   - Subgroup attacks: Exploiting small subgroups in DH
 */

import { gcd, euler_phi, power_mod, factor, divisors, lcm, Mod, Zmod } from 'sagemath-ts';

// =============================================================================
// Example 1: Computing Element Orders
// =============================================================================
console.log("=== Example 1: Computing Element Orders ===\n");

/**
 * The order of g is found by computing g^1, g^2, g^3, ... until we get 1.
 */
function computeOrder(g: bigint, groupOp: (a: bigint, b: bigint) => bigint, identity: bigint, maxOrder: bigint): bigint {
  let current = g;
  let order = 1n;

  while (current !== identity && order <= maxOrder) {
    current = groupOp(current, g);
    order++;
  }

  return current === identity ? order : -1n;  // -1 if order not found
}

const p = 11n;
console.log(`Orders of elements in (Z/${p}Z)*:\n`);

const multMod = (a: bigint, b: bigint) => (a * b) % p;

for (let g = 1n; g < p; g++) {
  const order = computeOrder(g, multMod, 1n, p);
  const divides = (p - 1n) % order === 0n;

  // Show the powers
  const powers: bigint[] = [];
  let current = 1n;
  for (let i = 0n; i <= order; i++) {
    powers.push(current);
    current = (current * g) % p;
  }

  console.log(`  g = ${g.toString().padStart(2)}: ord(g) = ${order.toString().padStart(2)}, divides ${p-1n}? ${divides}`);
  console.log(`         ${g}^k = ${powers.join(', ')}`);
}

console.log();

// =============================================================================
// Example 2: Lagrange's Theorem
// =============================================================================
console.log("=== Example 2: Lagrange's Theorem in Action ===\n");

/**
 * Lagrange's Theorem tells us that element orders must divide the group order.
 * Let's verify this experimentally.
 */
const modulus = 12n;
const groupOrder = euler_phi(modulus);

console.log(`Group: (Z/${modulus}Z)*`);
console.log(`Group order: phi(${modulus}) = ${groupOrder}`);
console.log(`Divisors of ${groupOrder}: ${divisors(groupOrder).join(', ')}\n`);

console.log("Element orders (must be divisors of group order):");
const units: bigint[] = [];
for (let a = 1n; a < modulus; a++) {
  if (gcd(a, modulus) === 1n) {
    units.push(a);
  }
}

const orderCounts: Map<bigint, number> = new Map();
for (const g of units) {
  let order = 1n;
  let current = g;
  while (current !== 1n) {
    current = (current * g) % modulus;
    order++;
  }

  orderCounts.set(order, (orderCounts.get(order) || 0) + 1);
  console.log(`  ord(${g.toString().padStart(2)}) = ${order}`);
}

console.log("\nOrder distribution:");
for (const [order, count] of [...orderCounts].sort((a, b) => Number(a[0] - b[0]))) {
  console.log(`  Order ${order}: ${count} element(s)`);
}

console.log();

// =============================================================================
// Example 3: Efficient Order Computation
// =============================================================================
console.log("=== Example 3: Efficient Order Computation ===\n");

/**
 * For large groups, we don't want to compute all powers naively.
 * We can use the fact that ord(g) divides |G| to check only divisors.
 */
function efficientOrder(g: bigint, n: bigint): bigint {
  if (gcd(g, n) !== 1n) {
    throw new Error("g must be coprime to n");
  }

  const groupOrder = euler_phi(n);
  const divs = divisors(groupOrder);

  // Sort divisors in ascending order
  divs.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));

  for (const d of divs) {
    if (power_mod(g, d, n) === 1n) {
      return d;
    }
  }

  return groupOrder;  // Should always divide
}

const bigMod = 101n;
console.log(`Efficient order computation in (Z/${bigMod}Z)*:\n`);

console.log(`Group order: phi(${bigMod}) = ${euler_phi(bigMod)}`);
console.log(`Divisors to check: ${divisors(euler_phi(bigMod)).join(', ')}\n`);

for (const g of [2n, 3n, 5n, 10n]) {
  const ord = efficientOrder(g, bigMod);
  console.log(`  ord(${g}) = ${ord}`);
}

console.log();

// =============================================================================
// Example 4: Fermat's Little Theorem
// =============================================================================
console.log("=== Example 4: Fermat's Little Theorem ===\n");

/**
 * Fermat's Little Theorem: For prime p and a not divisible by p,
 *   a^(p-1) = 1 (mod p)
 *
 * This is a direct consequence of Lagrange's theorem:
 *   ord(a) divides |G| = p - 1
 *   so a^(p-1) = a^(k * ord(a)) = (a^ord(a))^k = 1^k = 1
 */
const fermatP = 17n;
console.log(`Verifying Fermat's Little Theorem for p = ${fermatP}:\n`);
console.log(`  For all a in (Z/${fermatP}Z)*, we should have a^${fermatP - 1n} = 1\n`);

for (let a = 1n; a < fermatP; a++) {
  const result = power_mod(a, fermatP - 1n, fermatP);
  console.log(`  ${a.toString().padStart(2)}^${fermatP - 1n} = ${result} mod ${fermatP}`);
}

console.log("\nAll results are 1! Fermat's Little Theorem confirmed.");

console.log();

// =============================================================================
// Example 5: Euler's Theorem (Generalization)
// =============================================================================
console.log("=== Example 5: Euler's Theorem ===\n");

/**
 * Euler's Theorem: For gcd(a, n) = 1,
 *   a^phi(n) = 1 (mod n)
 *
 * This generalizes Fermat's Little Theorem to composite moduli.
 */
const eulerN = 15n;
const phi = euler_phi(eulerN);

console.log(`Verifying Euler's Theorem for n = ${eulerN}:`);
console.log(`  phi(${eulerN}) = ${phi}`);
console.log(`  For all a coprime to ${eulerN}, a^${phi} = 1 mod ${eulerN}\n`);

for (let a = 1n; a < eulerN; a++) {
  if (gcd(a, eulerN) === 1n) {
    const result = power_mod(a, phi, eulerN);
    console.log(`  ${a.toString().padStart(2)}^${phi} = ${result} mod ${eulerN}`);
  }
}

console.log("\nEuler's Theorem is the foundation of RSA!");

console.log();

// =============================================================================
// Example 6: Order and Subgroups
// =============================================================================
console.log("=== Example 6: Subgroups Generated by Elements ===\n");

/**
 * Each element g generates a cyclic subgroup <g> of order ord(g).
 * The subgroup structure is determined by element orders.
 */
const subMod = 13n;
console.log(`Subgroups of (Z/${subMod}Z)*:\n`);

// Find all distinct subgroups
const subgroups: Map<string, bigint[]> = new Map();

for (let g = 1n; g < subMod; g++) {
  if (gcd(g, subMod) !== 1n) continue;

  // Generate subgroup
  const subgroup: bigint[] = [];
  let current = 1n;
  do {
    subgroup.push(current);
    current = (current * g) % subMod;
  } while (current !== 1n);

  // Use sorted string representation as key
  const key = subgroup.sort((a, b) => Number(a - b)).join(',');
  if (!subgroups.has(key)) {
    subgroups.set(key, subgroup);
  }
}

console.log("Distinct subgroups:");
for (const [_key, elements] of subgroups) {
  console.log(`  {${elements.join(', ')}}  (order ${elements.length})`);
}

console.log(`\n(Z/${subMod}Z)* is cyclic, so it has one subgroup for each divisor of ${subMod - 1n}.`);

console.log();

// =============================================================================
// Example 7: Product of Element Orders
// =============================================================================
console.log("=== Example 7: LCM and Element Orders ===\n");

/**
 * The LCM of all element orders is called the "exponent" of the group.
 * For cyclic groups, exponent = group order.
 * For non-cyclic groups, exponent < group order.
 */
function groupExponent(n: bigint): bigint {
  const orders: bigint[] = [];

  for (let g = 1n; g < n; g++) {
    if (gcd(g, n) === 1n) {
      const ord = efficientOrder(g, n);
      orders.push(ord);
    }
  }

  return lcm(orders);
}

const testModuli = [7n, 8n, 12n, 15n];

console.log("Comparing group order vs exponent:\n");
for (const n of testModuli) {
  const phi = euler_phi(n);
  const exp = groupExponent(n);
  const isCyclic = phi === exp;

  console.log(`  (Z/${n.toString().padStart(2)}Z)*: order = ${phi.toString().padStart(2)}, exponent = ${exp.toString().padStart(2)}, cyclic? ${isCyclic}`);
}

console.log("\nA group is cyclic iff its exponent equals its order.");

console.log();

// =============================================================================
// Example 8: Order of Powers
// =============================================================================
console.log("=== Example 8: Order of Powers ===\n");

/**
 * Key formula: ord(g^k) = ord(g) / gcd(k, ord(g))
 *
 * This tells us how to find elements of any desired order.
 */
const ordMod = 13n;
const generator = 2n;  // A primitive root mod 13
const genOrder = efficientOrder(generator, ordMod);

console.log(`In (Z/${ordMod}Z)*, let g = ${generator} with ord(g) = ${genOrder}\n`);
console.log(`Formula: ord(g^k) = ord(g) / gcd(k, ord(g))\n`);

for (let k = 1n; k <= genOrder; k++) {
  const gk = power_mod(generator, k, ordMod);
  const expectedOrder = genOrder / gcd(k, genOrder);
  const actualOrder = efficientOrder(gk, ordMod);

  console.log(`  g^${k.toString().padStart(2)} = ${gk.toString().padStart(2)}: expected ord = ${genOrder}/${gcd(k, genOrder)} = ${expectedOrder}, actual = ${actualOrder}`);
}

console.log();

// =============================================================================
// Example 9: Finding Elements of Specific Order
// =============================================================================
console.log("=== Example 9: Finding Elements of Specific Order ===\n");

/**
 * To find an element of order d, find a generator g and compute g^(|G|/d).
 * This only works if d divides |G|.
 */
function findElementOfOrder(targetOrder: bigint, n: bigint): bigint | null {
  const groupOrder = euler_phi(n);

  if (groupOrder % targetOrder !== 0n) {
    return null;  // No such element exists
  }

  // Find a generator
  for (let g = 2n; g < n; g++) {
    if (gcd(g, n) === 1n && efficientOrder(g, n) === groupOrder) {
      // g is a generator, compute g^(groupOrder/targetOrder)
      const exponent = groupOrder / targetOrder;
      return power_mod(g, exponent, n);
    }
  }

  return null;
}

const targetMod = 17n;
console.log(`Finding elements of specific orders in (Z/${targetMod}Z)*:\n`);

const possibleOrders = divisors(euler_phi(targetMod));
for (const d of possibleOrders) {
  const elem = findElementOfOrder(d, targetMod);
  if (elem !== null) {
    const verify = efficientOrder(elem, targetMod);
    console.log(`  Order ${d.toString().padStart(2)}: element ${elem?.toString().padStart(2)}, verified ord = ${verify}`);
  }
}

console.log();

// =============================================================================
// Exercise: Verify Lagrange's Theorem
// =============================================================================
console.log("=== Exercise: Verify Lagrange's Theorem ===\n");

/**
 * EXERCISE: Show that every element order in (Z/17Z)* divides 16.
 */
const verifyMod = 17n;
const verifyGroupOrder = verifyMod - 1n;

console.log(`Verifying Lagrange's Theorem for (Z/${verifyMod}Z)*:\n`);
console.log(`Group order = ${verifyGroupOrder}`);
console.log(`Divisors of ${verifyGroupOrder}: ${divisors(verifyGroupOrder).join(', ')}\n`);

let allDivide = true;
for (let g = 1n; g < verifyMod; g++) {
  const ord = efficientOrder(g, verifyMod);
  const divides = verifyGroupOrder % ord === 0n;
  if (!divides) allDivide = false;
  console.log(`  ord(${g.toString().padStart(2)}) = ${ord.toString().padStart(2)}, divides ${verifyGroupOrder}? ${divides}`);
}

console.log(`\nAll element orders divide group order: ${allDivide}`);

console.log();

// =============================================================================
// Cryptographic Application: RSA via Euler's Theorem
// =============================================================================
console.log("=== Cryptographic Application: Why RSA Works ===\n");

/**
 * RSA works because of Euler's theorem:
 *   m^(ed) = m^(1 + k*phi(n)) = m * (m^phi(n))^k = m * 1^k = m (mod n)
 *
 * The last step uses Euler's theorem: m^phi(n) = 1 when gcd(m, n) = 1.
 */
console.log("RSA correctness proof:\n");

const rsaP = 11n;
const rsaQ = 13n;
const rsaN = rsaP * rsaQ;
const rsaPhi = (rsaP - 1n) * (rsaQ - 1n);
const rsaE = 7n;

// Find d such that e*d = 1 mod phi(n)
let rsaD = 1n;
while ((rsaE * rsaD) % rsaPhi !== 1n) rsaD++;

console.log(`RSA parameters:`);
console.log(`  n = ${rsaP} * ${rsaQ} = ${rsaN}`);
console.log(`  phi(n) = ${rsaPhi}`);
console.log(`  e = ${rsaE}, d = ${rsaD}`);
console.log(`  e * d = ${rsaE * rsaD} = 1 + ${(rsaE * rsaD - 1n) / rsaPhi} * ${rsaPhi}`);

const message = 42n;
console.log(`\nEncrypt/Decrypt message m = ${message}:`);

const encrypted = power_mod(message, rsaE, rsaN);
const decrypted = power_mod(encrypted, rsaD, rsaN);

console.log(`  Encrypt: m^e = ${message}^${rsaE} = ${encrypted} mod ${rsaN}`);
console.log(`  Decrypt: c^d = ${encrypted}^${rsaD} = ${decrypted} mod ${rsaN}`);

console.log(`\nWhy it works:`);
console.log(`  m^(ed) = m^(${rsaE * rsaD}) = m^(1 + ${(rsaE * rsaD - 1n) / rsaPhi} * phi(n))`);
console.log(`        = m * (m^phi(n))^${(rsaE * rsaD - 1n) / rsaPhi}`);
console.log(`        = m * 1^${(rsaE * rsaD - 1n) / rsaPhi}  [by Euler's theorem]`);
console.log(`        = m = ${message}`);
