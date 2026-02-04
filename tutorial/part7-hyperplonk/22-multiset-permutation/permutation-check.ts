/**
 * # Permutation Check: Verifying Wiring Constraints
 *
 * The permutation check is essential for verifying that values are
 * correctly "wired" between different positions in a circuit.
 * In PLONK-style systems, this is how we enforce copy constraints.
 *
 * ## The Problem
 *
 * Given polynomials f, g and a permutation sigma: B_n -> B_n, prove:
 *   f(b) = g(sigma(b)) for all b in B_n
 *
 * Equivalently: f is a permutation of g according to sigma.
 *
 * ## Why Permutation Checks Matter
 *
 * - Copy constraints: same value at different positions
 * - Wiring: connecting outputs to inputs in circuits
 * - Memory access: matching reads to writes
 *
 * @module tutorial/part7-hyperplonk/22-multiset-permutation/permutation-check
 */

import { NotImplementedError } from 'sagemath-ts';

console.log("=== Permutation Check: Verifying Wiring ===\n");

// ============================================================================
// Field Arithmetic (simple modular arithmetic)
// ============================================================================

const p = 97n;

const mod = (x: bigint): bigint => ((x % p) + p) % p;
const add = (a: bigint, b: bigint): bigint => mod(a + b);
const sub = (a: bigint, b: bigint): bigint => mod(a - b);
const mul = (a: bigint, b: bigint): bigint => mod(a * b);

function booleanHypercube(n: number): number[][] {
  if (n === 0) return [[]];
  const smaller = booleanHypercube(n - 1);
  return [...smaller.map(v => [0, ...v]), ...smaller.map(v => [1, ...v])];
}

// ============================================================================
// Section 1: What is a Permutation Check?
// ============================================================================

console.log("--- Section 1: Permutation Checks ---\n");

/**
 * PERMUTATION CHECK
 *
 * A permutation sigma is a bijection from B_n to B_n.
 * It shuffles the indices without losing or duplicating any.
 *
 * PERMUTATION CHECK CLAIM:
 * Given f, g, and sigma, prove: f(b) = g(sigma(b)) for all b.
 *
 * This means if we apply sigma to the indices of g, we get f.
 *
 * EXAMPLE:
 * sigma: 00 -> 01, 01 -> 10, 10 -> 00, 11 -> 11
 *
 * f(00) = g(01)
 * f(01) = g(10)
 * f(10) = g(00)
 * f(11) = g(11)
 *
 * WHY THIS MATTERS:
 *
 * In circuits, we often need to enforce that values at different
 * wire positions are equal. For example:
 * - Output of gate i must equal input of gate j
 * - Same variable used in multiple places
 *
 * The permutation sigma encodes which positions must match.
 */

console.log("Permutation Check:");
console.log("  Given: f, g polynomials and permutation sigma: B_n -> B_n");
console.log("  Prove: f(b) = g(sigma(b)) for all b\n");

console.log("Example permutation on B_2:");
console.log("  sigma(00) = 01");
console.log("  sigma(01) = 10");
console.log("  sigma(10) = 00");
console.log("  sigma(11) = 11\n");

// ============================================================================
// Section 2: Encoding Permutations
// ============================================================================

console.log("--- Section 2: Encoding Permutations ---\n");

/**
 * ENCODING PERMUTATIONS AS POLYNOMIALS
 *
 * We need to encode sigma as a polynomial for the protocol.
 *
 * APPROACH 1: Index polynomial
 * Define ID(b) = index of b as an integer (e.g., 00=0, 01=1, 10=2, 11=3)
 * Define sigma_poly(b) = index of sigma(b)
 *
 * Then the constraint becomes:
 *   f(b) = g at position sigma_poly(b)
 *
 * APPROACH 2: Tag-based (HyperPlonk style)
 *
 * Instead of checking f(b) = g(sigma(b)) directly, we check:
 *   {(f(b), ID(b))} = {(g(b), sigma(b))} as multisets
 *
 * The second component "tags" each value with its position.
 * If f(b) = g(sigma(b)), then:
 *   (f(b), b) maps to (g(sigma(b)), sigma(b))
 *
 * SIMPLIFIED APPROACH:
 * Use a random combination:
 *   f(b) + gamma * ID(b) should match g(sigma(b)) + gamma * sigma(b)
 *
 * If sigma is correct, the multisets match.
 */

console.log("Encoding permutations:\n");

console.log("Tag-based approach:");
console.log("  Create pairs (value, position)");
console.log("  Check {(f(b), b)} = {(g(b), sigma(b))} as multisets\n");

console.log("Random combination:");
console.log("  h(b) = f(b) + gamma * ID(b)");
console.log("  k(b) = g(b) + gamma * sigma(b)");
console.log("  Check {h(b)} = {k(b)} as multisets\n");

// ============================================================================
// Section 3: The Grand Product Argument
// ============================================================================

console.log("--- Section 3: Grand Product Argument ---\n");

/**
 * GRAND PRODUCT ARGUMENT FOR PERMUTATION
 *
 * To check permutation with tags, we use the multiset equality reduction:
 *
 *   {f(b) + gamma*ID(b)} = {g(b) + gamma*sigma(b)}
 *
 * <=> for random beta:
 *   prod_b (beta + f(b) + gamma*ID(b)) = prod_b (beta + g(b) + gamma*sigma(b))
 *
 * This is a ProductCheck!
 *
 * COMBINING MULTIPLE COLUMNS:
 *
 * In PLONK, we have multiple wire columns (a, b, c for left/right/output).
 * Each needs its own permutation check.
 *
 * We can batch them using another random combination:
 *
 *   prod_b (beta + a(b) + gamma*ID_a(b)) * (beta + b(b) + gamma*ID_b(b)) * ...
 *   = prod_b (beta + a'(b) + gamma*sigma_a(b)) * ...
 *
 * where a', b', ... are the permuted values.
 *
 * HYPERPLONK APPROACH:
 *
 * HyperPlonk uses a similar structure but leverages multilinear
 * polynomials on the boolean hypercube.
 */

console.log("Grand Product for Permutation:\n");

console.log("  Define h(b) = beta + f(b) + gamma * ID(b)");
console.log("  Define k(b) = beta + g(b) + gamma * sigma(b)");
console.log("  Permutation valid iff prod_b h(b) = prod_b k(b)\n");

console.log("This reduces permutation to ProductCheck!\n");

// ============================================================================
// Section 4: Concrete Example
// ============================================================================

console.log("--- Section 4: Concrete Example ---\n");

interface HypercubeFunction {
  numVars: number;
  values: bigint[];
}

// Permutation: 00<->01, 10<->11 (swap within pairs)
// sigma(00) = 01, sigma(01) = 00, sigma(10) = 11, sigma(11) = 10
const sigma = [1, 0, 3, 2];  // sigma[i] = index of sigma(b_i)

const cube = booleanHypercube(2);
console.log("Permutation sigma:");
for (let i = 0; i < cube.length; i++) {
  console.log(`  sigma(${cube[i]!.join(',')}) = ${cube[sigma[i]!]!.join(',')}`);
}
console.log();

// f values
const fEx: HypercubeFunction = {
  numVars: 2,
  values: [5n, 3n, 9n, 7n]  // f(00)=5, f(01)=3, f(10)=9, f(11)=7
};

// g such that f(b) = g(sigma(b))
// g(01) = f(00) = 5
// g(00) = f(01) = 3
// g(11) = f(10) = 9
// g(10) = f(11) = 7
const gEx: HypercubeFunction = {
  numVars: 2,
  values: [3n, 5n, 7n, 9n]  // g(00)=3, g(01)=5, g(10)=7, g(11)=9
};

console.log("f values:");
for (let i = 0; i < cube.length; i++) {
  console.log(`  f(${cube[i]!.join(',')}) = ${fEx.values[i]!}`);
}

console.log("\ng values:");
for (let i = 0; i < cube.length; i++) {
  console.log(`  g(${cube[i]!.join(',')}) = ${gEx.values[i]!}`);
}

// Verify: f(b) = g(sigma(b))
console.log("\nVerifying f(b) = g(sigma(b)):");
let allMatch = true;
for (let i = 0; i < cube.length; i++) {
  const fVal = fEx.values[i]!;
  const gSigmaVal = gEx.values[sigma[i]!]!;
  const match = fVal === gSigmaVal;
  console.log(`  f(${cube[i]!.join(',')}) = ${fVal}, g(sigma(...)) = ${gSigmaVal} ${match ? 'OK' : 'FAIL'}`);
  allMatch = allMatch && match;
}
console.log(`\nAll constraints satisfied? ${allMatch ? 'YES' : 'NO'}\n`);

// ============================================================================
// Section 5: Running the Grand Product
// ============================================================================

console.log("--- Section 5: Grand Product Check ---\n");

// Random challenges
const beta = 11n;
const gamma = 17n;

console.log(`Random challenges: beta=${beta}, gamma=${gamma}\n`);

// Compute tagged values
// h(b) = beta + f(b) + gamma * ID(b)
// k(b) = beta + g(b) + gamma * sigma(b)

console.log("Tagged values:");
let hProd = 1n;
let kProd = 1n;

for (let i = 0; i < cube.length; i++) {
  const id = BigInt(i);
  const sigmaId = BigInt(sigma[i]!);

  const h = add(add(beta, fEx.values[i]!), mul(gamma, id));
  const k = add(add(beta, gEx.values[i]!), mul(gamma, sigmaId));

  console.log(`  b=${cube[i]!.join(',')}: h = ${h}, k = ${k}`);

  hProd = mul(hProd, h);
  kProd = mul(kProd, k);
}

console.log(`\nprod h = ${hProd}`);
console.log(`prod k = ${kProd}`);
console.log(`Products equal? ${hProd === kProd ? 'YES' : 'NO'}\n`);

// ============================================================================
// Section 6: Detecting Invalid Permutation
// ============================================================================

console.log("--- Section 6: Detecting Invalid Permutation ---\n");

// g' that does NOT satisfy the permutation
const gBad: HypercubeFunction = {
  numVars: 2,
  values: [5n, 5n, 7n, 9n]  // g(00)=5 instead of 3
};

console.log("Bad g values (breaks permutation):");
for (let i = 0; i < cube.length; i++) {
  console.log(`  g(${cube[i]!.join(',')}) = ${gBad.values[i]!}`);
}

console.log("\nVerifying f(b) = g_bad(sigma(b)):");
for (let i = 0; i < cube.length; i++) {
  const fVal = fEx.values[i]!;
  const gSigmaVal = gBad.values[sigma[i]!]!;
  const match = fVal === gSigmaVal;
  console.log(`  f(${cube[i]!.join(',')}) = ${fVal}, g_bad(sigma(...)) = ${gSigmaVal} ${match ? 'OK' : 'FAIL'}`);
}

// Run grand product check
console.log("\nGrand product check on bad g:");
let kBadProd = 1n;
for (let i = 0; i < cube.length; i++) {
  const sigmaId = BigInt(sigma[i]!);
  const k = add(add(beta, gBad.values[i]!), mul(gamma, sigmaId));
  kBadProd = mul(kBadProd, k);
}

console.log(`  prod h = ${hProd}`);
console.log(`  prod k_bad = ${kBadProd}`);
console.log(`  Equal? ${hProd === kBadProd ? 'YES (false positive!)' : 'NO (caught!)'}\n`);

// ============================================================================
// Section 7: HyperPlonk's Wiring Check
// ============================================================================

console.log("--- Section 7: HyperPlonk Wiring ---\n");

/**
 * HYPERPLONK WIRING CHECK
 *
 * In HyperPlonk, wiring constraints enforce that wire values are
 * correctly copied between gates.
 *
 * Setup:
 * - Wire polynomials w_1, ..., w_k over B_mu (one per wire type)
 * - Permutation sigma that specifies which wire values must be equal
 *
 * The permutation sigma is over the index space B_mu x {1, ..., k}.
 * sigma(b, i) = (b', j) means w_i(b) must equal w_j(b').
 *
 * ENCODING:
 *
 * Define ID_i(b) as a unique identifier for position (b, i).
 * For example: ID_i(b) = b * k + i
 *
 * Then the permutation check becomes:
 *   prod_{b,i} (beta + w_i(b) + gamma * ID_i(b))
 *   = prod_{b,i} (beta + w_i(b) + gamma * sigma(b,i))
 *
 * MULTILINEAR ENCODING:
 *
 * The permutation sigma is encoded as a multilinear polynomial.
 * ID_i can also be made multilinear.
 *
 * The entire check reduces to ProductCheck on multilinear polynomials.
 */

console.log("HyperPlonk wiring:");
console.log("  Wire polynomials w_1, ..., w_k");
console.log("  Permutation sigma specifies equalities\n");

console.log("Encoding:");
console.log("  ID(b, i) = unique position identifier");
console.log("  Check: w_i(b) = w_j(sigma(b,i)) for all (b,i)\n");

console.log("Reduces to ProductCheck on:");
console.log("  prod (beta + w_i(b) + gamma * ID(b,i))");
console.log("  = prod (beta + w_i(b) + gamma * sigma(b,i))\n");

// ============================================================================
// Section 8: Implementation Stubs
// ============================================================================

console.log("--- Section 8: Implementation Stubs ---\n");

/**
 * Basic permutation check.
 *
 * Proves f(b) = g(sigma(b)) for all b.
 *
 * @throws NotImplementedError - Full implementation pending
 */
export function permutationCheck(
  f: HypercubeFunction,
  g: HypercubeFunction,
  sigma: number[],
  beta: bigint,
  gamma: bigint,
  getChallenge: () => bigint
): { accepted: boolean } {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: permutationCheck - ' +
    'Full permutation check with ProductCheck integration'
  );
}

/**
 * HyperPlonk wiring check for multiple wire columns.
 *
 * @throws NotImplementedError - Wiring check pending
 */
export function hyperplonkWiringCheck(
  wires: HypercubeFunction[],
  sigma: number[][],
  beta: bigint,
  gamma: bigint,
  getChallenge: () => bigint
): { accepted: boolean } {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: hyperplonkWiringCheck - ' +
    'HyperPlonk wiring constraint check'
  );
}

/**
 * Copy constraint check (simplified permutation for copies).
 *
 * @throws NotImplementedError - Copy constraint check pending
 */
export function copyConstraintCheck(
  values: bigint[],
  equalities: [number, number][],
  getChallenge: () => bigint
): { accepted: boolean } {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: copyConstraintCheck - ' +
    'Copy constraint verification via permutation'
  );
}

console.log("Stub functions:");
console.log("  - permutationCheck: basic permutation");
console.log("  - hyperplonkWiringCheck: multi-wire wiring");
console.log("  - copyConstraintCheck: copy constraints\n");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");

console.log("Permutation Check:");
console.log("  Claim: f(b) = g(sigma(b)) for all b");
console.log("  Method: Tagged multiset equality via ProductCheck\n");

console.log("Key techniques:");
console.log("  1. Tag values with positions: (value, position)");
console.log("  2. Random combination: f(b) + gamma * ID(b)");
console.log("  3. ProductCheck on tagged values\n");

console.log("Applications:");
console.log("  - PLONK copy constraints");
console.log("  - HyperPlonk wiring");
console.log("  - Circuit correctness\n");

console.log("Next: Full HyperPlonk protocol");

// Exports
export { type HypercubeFunction };
