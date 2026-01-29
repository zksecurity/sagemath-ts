/**
 * # Group Structure of E(F_q)
 *
 * The set of rational points E(F_q) on an elliptic curve forms an
 * abelian group under point addition. Understanding this structure
 * is crucial for cryptographic applications.
 *
 * ## The Structure Theorem
 *
 * For any elliptic curve E over F_q:
 *
 *   E(F_q) is isomorphic to Z/n1 x Z/n2
 *
 * where n2 divides n1, and n1 * n2 = #E(F_q).
 *
 * - If n2 = 1: The group is CYCLIC (one generator suffices)
 * - If n2 > 1: The group needs TWO generators
 *
 * ## Cryptographic Implications
 *
 * - We want large prime-order subgroups for ECDLP hardness
 * - Ideally, #E is prime (fully cyclic)
 * - Small cofactor h = #E / (prime order) is acceptable
 *
 * @module tutorial/elliptic-curves/group-structure
 */

import {
  EllipticCurve,
  abelian_group,
  discrete_log,
} from 'sagemath-ts';
import { FiniteFieldPrime, factor } from 'sagemath-ts';

// ============================================================================
// The Group Structure Theorem
// ============================================================================

console.log("=== The Group Structure Theorem ===\n");

/**
 * E(F_q) = Z/n1 x Z/n2 where:
 * - n2 divides n1
 * - n2 divides gcd(n1, q-1)
 * - n1 * n2 = #E(F_q)
 *
 * Most commonly, n2 = 1 and the group is cyclic.
 */

const F101 = new FiniteFieldPrime(101n);
const E = EllipticCurve(F101, [2n, 3n]);

const cardinality = E.cardinality();
console.log(`Curve: ${E.toString()}`);
console.log(`Cardinality: #E(F_101) = ${cardinality}`);
console.log("");

// Get the group structure
const groupInfo = abelian_group(E);

console.log("Group structure:");
console.log(`  Invariants: [${groupInfo.invariants.join(', ')}]`);
console.log(`  Number of generators: ${groupInfo.generators.length}`);

for (let i = 0; i < groupInfo.generators.length; i++) {
  const gen = groupInfo.generators[i];
  const ord = gen!.order();
  console.log(`  Generator ${i + 1}: ${gen!.toString()}, order = ${ord}`);
}

if (groupInfo.invariants.length === 1) {
  console.log("\nThis group is CYCLIC (one generator suffices)");
} else {
  console.log(`\nThis group is Z/${groupInfo.invariants[0]} x Z/${groupInfo.invariants[1]}`);
}
console.log("");

// ============================================================================
// Cyclic Groups (n2 = 1)
// ============================================================================

console.log("=== Cyclic Groups ===\n");

/**
 * When E(F_q) is cyclic:
 * - Any non-identity point of order #E generates the whole group
 * - Every point P satisfies P = nG for some n (discrete log exists)
 * - This is the ideal case for ECDH and ECDSA
 */

// Find a curve with prime order (definitely cyclic)
console.log("Searching for curves with prime order over F_101...");

for (let a = 1n; a <= 20n; a++) {
  for (let b = 1n; b <= 20n; b++) {
    try {
      const E_test = EllipticCurve(F101, [a, b]);
      const order = E_test.cardinality();

      // Check if prime using trial division
      let isPrime = order > 1n;
      for (let d = 2n; d * d <= order && isPrime; d++) {
        if (order % d === 0n) isPrime = false;
      }

      if (isPrime) {
        console.log(`  Found: y^2 = x^3 + ${a}x + ${b}`);
        console.log(`  Order: ${order} (prime)`);
        console.log(`  Group: Z/${order} (cyclic)`);

        // Any non-identity point generates the group
        const G = E_test.random_point();
        console.log(`  Generator: ${G.toString()}`);
        console.log(`  Order of G: ${G.order()}`);
        console.log("");
        break;
      }
    } catch {
      // Singular
    }
  }
}

// ============================================================================
// Non-Cyclic Groups
// ============================================================================

console.log("=== Non-Cyclic Groups ===\n");

/**
 * When n2 > 1, the group has the form Z/n1 x Z/n2.
 * This happens when the curve has "extra" torsion points.
 *
 * For example, if #E = 12 and n2 = 2:
 *   E(F_q) = Z/6 x Z/2
 *
 * This means there are three points of order 2 (the full 2-torsion is rational).
 */

// Search for a non-cyclic group
console.log("Searching for non-cyclic groups...");

for (const p of [7n, 11n, 13n, 17n, 19n, 23n]) {
  const Fp = new FiniteFieldPrime(p);

  for (let a = 0n; a < p; a++) {
    for (let b = 0n; b < p; b++) {
      try {
        const E_nc = EllipticCurve(Fp, [a, b]);
        const gens = E_nc.gens();

        if (gens.length >= 2) {
          const order = E_nc.cardinality();
          const n1 = gens[0]!.order();
          const n2 = order / n1;

          if (n2 > 1n && order % (n1 * n2) === 0n) {
            console.log(`Found non-cyclic group over F_${p}:`);
            console.log(`  Curve: y^2 = x^3 + ${a}x + ${b}`);
            console.log(`  Order: ${order}`);
            console.log(`  Structure: Z/${n1} x Z/${n2}`);
            console.log(`  Generators:`);
            for (const g of gens) {
              console.log(`    ${g.toString()}, order ${g.order()}`);
            }
            console.log("");
            break;
          }
        }
      } catch {
        // Singular
      }
    }
  }
}

// ============================================================================
// Point Orders and Subgroups
// ============================================================================

console.log("=== Point Orders and Subgroups ===\n");

/**
 * By Lagrange's theorem, the order of any point divides #E(F_q).
 *
 * For cryptography, we want points of LARGE PRIME ORDER.
 * If #E = h * n where n is prime:
 * - h is the "cofactor"
 * - Points of order n form a subgroup <G> where G has order n
 */

console.log(`Analyzing point orders on E: y^2 = x^3 + 2x + 3 over F_101`);
console.log(`Curve order: ${cardinality}`);

// Factor the curve order
const factorization = factor(cardinality);
console.log(`Factorization: ${cardinality} = ${factorization.map(([p, e]) => e > 1 ? `${p}^${e}` : `${p}`).join(' * ')}`);
console.log("");

// Find the distribution of point orders
const orderCounts = new Map<bigint, number>();

// Sample random points and compute their orders
console.log("Sampling point orders:");
const sampleSize = 20;

for (let i = 0; i < sampleSize; i++) {
  const P = E.random_point();
  if (!P.isZero()) {
    const ord = P.order();
    orderCounts.set(ord, (orderCounts.get(ord) || 0) + 1);
  }
}

const sortedOrders = Array.from(orderCounts.entries()).sort((a, b) =>
  a[0] < b[0] ? -1 : 1
);

for (const [ord, count] of sortedOrders) {
  const divides = cardinality % ord === 0n;
  console.log(`  Order ${ord}: ${count} points (divides #E? ${divides})`);
}
console.log("");

// ============================================================================
// Finding Generators
// ============================================================================

console.log("=== Finding Generators ===\n");

/**
 * To find a generator of the prime-order subgroup:
 *
 * 1. Factor #E = h * n where n is the largest prime factor
 * 2. Pick a random point P
 * 3. Compute G = h * P
 * 4. If G != O, then G has order n (or a divisor)
 *
 * Repeat until you get a point of order exactly n.
 */

// Find the largest prime factor
let n = cardinality;
let largestPrime = 1n;

for (let d = 2n; d * d <= n; d++) {
  while (n % d === 0n) {
    largestPrime = d;
    n /= d;
  }
}

if (n > 1n) {
  largestPrime = n;
}

const cofactor = cardinality / largestPrime;

console.log(`Curve order: ${cardinality}`);
console.log(`Largest prime factor (n): ${largestPrime}`);
console.log(`Cofactor (h): ${cofactor}`);
console.log("");

// Find a generator of the prime-order subgroup
console.log("Finding generator of the n-torsion subgroup:");

let G = E.zero();
let attempts = 0;

while (G.isZero() || G.order() !== largestPrime) {
  const P = E.random_point();
  G = P.mul(cofactor); // Clear the cofactor

  attempts++;
  if (attempts > 100) {
    console.log("  Could not find generator (unexpected)");
    break;
  }
}

if (!G.isZero()) {
  console.log(`  Generator G = ${G.toString()}`);
  console.log(`  Order of G: ${G.order()}`);
  console.log(`  Attempts needed: ${attempts}`);
}
console.log("");

// ============================================================================
// Discrete Logarithm in Subgroups
// ============================================================================

console.log("=== Discrete Logarithm ===\n");

/**
 * Given G (generator) and P = nG (public key), finding n is the ECDLP.
 *
 * The baby-step giant-step algorithm solves this in O(sqrt(order)) time.
 *
 * For cryptographic sizes (256 bits), sqrt(2^256) = 2^128, which is
 * still infeasible.
 */

// Demonstrate discrete log on small example
if (!G.isZero()) {
  const secretScalar = 42n % largestPrime;
  const publicPoint = G.mul(secretScalar);

  console.log(`Secret scalar: k = ${secretScalar}`);
  console.log(`Generator: G = ${G.toString()}`);
  console.log(`Public point: P = kG = ${publicPoint.toString()}`);
  console.log("");

  // Solve discrete log (only works for small groups!)
  console.log("Solving discrete log (brute force for demonstration):");

  try {
    const recoveredK = discrete_log(G, publicPoint, largestPrime);
    console.log(`  Recovered k = ${recoveredK}`);
    console.log(`  Verification: ${recoveredK}G = ${G.mul(recoveredK).toString()}`);
    console.log(`  Matches P? ${G.mul(recoveredK).eq(publicPoint)}`);
  } catch (e) {
    console.log(`  Error: ${(e as Error).message}`);
  }
}
console.log("");

// ============================================================================
// Torsion Subgroups
// ============================================================================

console.log("=== Torsion Subgroups ===\n");

/**
 * The n-TORSION subgroup E[n] consists of all points P with nP = O.
 *
 * Over an algebraically closed field:
 *   E[n] is isomorphic to Z/n x Z/n (for n coprime to characteristic)
 *
 * Over F_q, E[n] might be smaller:
 *   E[n](F_q) is a subgroup of E[n]
 *
 * The full n-torsion becomes rational over F_{q^k} where k is the
 * embedding degree.
 */

console.log("2-torsion points (points P with 2P = O):");

// Find 2-torsion: points with y = 0 (for short Weierstrass)
for (let x = 0n; x < 101n; x++) {
  const xElem = F101.__call__(x);
  const rhs = xElem.pow(3).add(F101.__call__(2n).mul(xElem)).add(F101.__call__(3n));

  if (rhs.isZero()) {
    // y = 0, so this is a 2-torsion point
    const P = E.point(x, 0n);
    console.log(`  (${x}, 0): 2P = ${P.double().toString()}`);
  }
}
console.log("  Point at infinity O is also 2-torsion (2O = O)");
console.log("");

// ============================================================================
// Summary: Group Structure and Security
// ============================================================================

console.log("=== Group Structure and Security ===\n");

console.log("For cryptographic use, prefer curves where:");
console.log("");
console.log("1. LARGE PRIME ORDER");
console.log("   - #E is prime, OR");
console.log("   - #E = h * n where n is large prime, h is small (1, 2, 4, 8)");
console.log("");
console.log("2. CYCLIC GROUP");
console.log("   - One generator suffices");
console.log("   - Simplifies protocol design");
console.log("");
console.log("3. NO SPECIAL STRUCTURE");
console.log("   - Avoid anomalous curves (#E = p)");
console.log("   - Avoid supersingular curves (low embedding degree)");
console.log("   - Avoid curves with known DLP weaknesses");
console.log("");

// ============================================================================
// Exercises
// ============================================================================

console.log("=== Exercises ===\n");

console.log("1. Find the group structure of E: y^2 = x^3 + x + 1 over F_13");
console.log("   Is it cyclic? What are the possible point orders?");
console.log("");

console.log("2. For E over F_31 with #E = 30:");
console.log("   - What are the possible group structures?");
console.log("   - Find a point of order 5, 6, 10, 15, 30 (if they exist)");
console.log("");

console.log("3. Implement baby-step giant-step to solve ECDLP more efficiently");
console.log("   than brute force");
console.log("");
