/**
 * # Multiset Equality Check
 *
 * A multiset is a set where elements can appear multiple times.
 * The multiset equality check proves that two multisets contain
 * exactly the same elements with the same multiplicities.
 *
 * ## The Problem
 *
 * Given polynomials f, g: B_n -> F, prove that:
 *   {f(b) : b in B_n} = {g(b) : b in B_n}  (as multisets)
 *
 * This means every value appears the same number of times in both.
 *
 * ## Why Multiset Checks Matter
 *
 * - Lookup arguments: prove values exist in a table
 * - Copy constraints: prove values are copied correctly
 * - Memory consistency: prove read values match written values
 *
 * @module tutorial/part7-hyperplonk/22-multiset-permutation/multiset-check
 */

import { NotImplementedError } from 'sagemath-ts';

console.log("=== Multiset Equality Check ===\n");

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
// Section 1: What is a Multiset?
// ============================================================================

console.log("--- Section 1: Multisets ---\n");

/**
 * MULTISET DEFINITION
 *
 * A multiset (or bag) is like a set, but elements can appear multiple times.
 *
 * Examples:
 *   {1, 2, 2, 3} - multiset with 2 appearing twice
 *   {a, a, b, c, c, c} - a appears twice, b once, c three times
 *
 * Two multisets are EQUAL if they contain exactly the same elements
 * with exactly the same multiplicities.
 *
 *   {1, 2, 2, 3} = {2, 1, 3, 2}  (order doesn't matter)
 *   {1, 2, 2, 3} != {1, 2, 3}    (different multiplicities)
 *
 * POLYNOMIAL ENCODING:
 *
 * We encode a multiset {a_1, ..., a_m} as the polynomial:
 *   P(X) = prod_i (X - a_i)
 *
 * Two multisets are equal iff their polynomials are equal.
 *
 * But comparing polynomials directly is expensive!
 *
 * SCHWARTZ-ZIPPEL TRICK:
 *
 * Instead, evaluate at a random point beta:
 *   prod_i (beta - a_i) = prod_j (beta - b_j)
 *
 * If multisets are different, this fails with high probability.
 */

console.log("Multiset: a set with multiplicities");
console.log("  {1, 2, 2, 3} - 2 appears twice");
console.log("  Two multisets equal iff same elements, same counts\n");

console.log("Polynomial encoding:");
console.log("  Multiset {a_1, ..., a_m} -> P(X) = prod(X - a_i)");
console.log("  Equal multisets -> equal polynomials\n");

console.log("Schwartz-Zippel trick:");
console.log("  Evaluate at random beta");
console.log("  prod(beta - a_i) = prod(beta - b_j) whp iff equal\n");

// ============================================================================
// Section 2: Multiset Check via ProductCheck
// ============================================================================

console.log("--- Section 2: Reduction to ProductCheck ---\n");

/**
 * MULTISET EQUALITY VIA PRODUCTCHECK
 *
 * Given f, g on B_n, we want to check:
 *   {f(b) : b in B_n} = {g(b) : b in B_n}
 *
 * Using the polynomial encoding:
 *   prod_b (X - f(b)) = prod_b (X - g(b))
 *
 * Evaluate at random beta:
 *   prod_b (beta - f(b)) = prod_b (beta - g(b))
 *
 * This is exactly ProductCheck with:
 *   f'(b) = beta - f(b)
 *   g'(b) = beta - g(b)
 *
 * PROTOCOL:
 * 1. Verifier sends random beta
 * 2. Define f'(b) = beta - f(b), g'(b) = beta - g(b)
 * 3. Run ProductCheck: prod_b f'(b) = prod_b g'(b)
 *
 * SOUNDNESS:
 * If multisets differ, the polynomials prod(X - f(b)) and prod(X - g(b))
 * are different (they have different roots). By Schwartz-Zippel,
 * they evaluate to the same value at random beta with probability
 * at most (2^n) / |F|.
 */

console.log("Multiset check reduces to ProductCheck:\n");

console.log("  {f(b)} = {g(b)} as multisets");
console.log("  <=> prod_b (X - f(b)) = prod_b (X - g(b)) as polynomials");
console.log("  <=> prod_b (beta - f(b)) = prod_b (beta - g(b)) whp\n");

console.log("Protocol:");
console.log("  1. Verifier sends random beta");
console.log("  2. Define f'(b) = beta - f(b), g'(b) = beta - g(b)");
console.log("  3. Run ProductCheck on (f', g')\n");

// ============================================================================
// Section 3: Concrete Example
// ============================================================================

console.log("--- Section 3: Concrete Example ---\n");

interface HypercubeFunction {
  numVars: number;
  values: bigint[];
}

// f(00)=3, f(01)=5, f(10)=3, f(11)=7
// Multiset: {3, 5, 3, 7} = {3, 3, 5, 7}
const fEx: HypercubeFunction = {
  numVars: 2,
  values: [3n, 5n, 3n, 7n]
};

// g(00)=7, g(01)=3, g(10)=5, g(11)=3
// Multiset: {7, 3, 5, 3} = {3, 3, 5, 7}
const gEx: HypercubeFunction = {
  numVars: 2,
  values: [7n, 3n, 5n, 3n]
};

const cube = booleanHypercube(2);

console.log("f values (multiset: {3, 5, 3, 7}):");
for (let i = 0; i < cube.length; i++) {
  console.log(`  f(${cube[i]!.join(',')}) = ${fEx.values[i]!}`);
}

console.log("\ng values (multiset: {7, 3, 5, 3}):");
for (let i = 0; i < cube.length; i++) {
  console.log(`  g(${cube[i]!.join(',')}) = ${gEx.values[i]!}`);
}

// Sort to compare multisets
const fSorted = [...fEx.values].sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
const gSorted = [...gEx.values].sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
const multisetEqual = fSorted.every((v, i) => v === gSorted[i]);
console.log(`\nMultisets equal? ${multisetEqual ? 'YES' : 'NO'}\n`);

// Verify with random beta
const beta = 11n;
console.log(`Random challenge beta = ${beta}\n`);

// Compute products
let fProd = 1n;
let gProd = 1n;
for (let i = 0; i < fEx.values.length; i++) {
  const fPrime = sub(beta, fEx.values[i]!);
  const gPrime = sub(beta, gEx.values[i]!);
  console.log(`  beta - f(${cube[i]!.join(',')}) = ${fPrime}, beta - g(...) = ${gPrime}`);
  fProd = mul(fProd, fPrime);
  gProd = mul(gProd, gPrime);
}

console.log(`\nprod(beta - f) = ${fProd}`);
console.log(`prod(beta - g) = ${gProd}`);
console.log(`Products equal? ${fProd === gProd ? 'YES' : 'NO'}\n`);

// ============================================================================
// Section 4: Detecting Non-Equal Multisets
// ============================================================================

console.log("--- Section 4: Detecting Non-Equal Multisets ---\n");

// Different multisets
// f: {3, 5, 3, 7} = {3, 3, 5, 7}
// h: {3, 5, 4, 7} = {3, 4, 5, 7}
const hEx: HypercubeFunction = {
  numVars: 2,
  values: [3n, 5n, 4n, 7n]
};

console.log("Comparing f and h:");
console.log("  f multiset: {3, 5, 3, 7}");
console.log("  h multiset: {3, 5, 4, 7}");
console.log("  Different! (f has two 3s, h has one 3 and one 4)\n");

// Check with beta
let hProd = 1n;
for (let i = 0; i < hEx.values.length; i++) {
  hProd = mul(hProd, sub(beta, hEx.values[i]!));
}

console.log(`prod(beta - f) = ${fProd}`);
console.log(`prod(beta - h) = ${hProd}`);
console.log(`Products equal? ${fProd === hProd ? 'YES (false positive!)' : 'NO (correct)'}\n`);

// ============================================================================
// Section 5: The Multiset Check Protocol
// ============================================================================

console.log("--- Section 5: Full Protocol ---\n");

/**
 * MULTISET CHECK PROTOCOL
 *
 * Public: Polynomials f, g (commitments)
 * Claim: {f(b) : b in B_n} = {g(b) : b in B_n} as multisets
 *
 * PROTOCOL:
 *
 * 1. V -> P: Random beta
 *
 * 2. P computes:
 *    - f'(b) = beta - f(b) for all b
 *    - g'(b) = beta - g(b) for all b
 *
 * 3. P, V run ProductCheck on (f', g')
 *    - Prover computes auxiliary h
 *    - ZeroCheck on constraint
 *
 * 4. V accepts iff ProductCheck accepts
 *
 * EFFICIENCY:
 * - Communication: O(n) field elements (same as ProductCheck)
 * - Prover: O(2^n) for computing f', g', plus ProductCheck
 * - Verifier: O(n) plus oracle queries
 *
 * SOUNDNESS:
 * If multisets differ:
 *   - Polynomials prod(X-f(b)) and prod(X-g(b)) are different
 *   - Degree at most 2^n each
 *   - By Schwartz-Zippel: agree at random beta with prob <= 2^n/|F|
 */

console.log("Multiset Check Protocol:\n");

console.log("  1. Verifier sends random beta");
console.log("  2. Prover computes f'(b) = beta - f(b)");
console.log("  3. Prover computes g'(b) = beta - g(b)");
console.log("  4. Run ProductCheck on (f', g')");
console.log("  5. Accept iff ProductCheck accepts\n");

console.log("Soundness: Error <= 2^n / |F|");
console.log(`  For n=${fEx.numVars}, |F|=${p}: error <= ${2 ** fEx.numVars}/${p} ~ ${(4 / Number(p)).toFixed(4)}\n`);

// ============================================================================
// Section 6: Application - Lookup Arguments
// ============================================================================

console.log("--- Section 6: Application - Lookup Arguments ---\n");

/**
 * LOOKUP ARGUMENTS
 *
 * Lookup arguments prove that a set of values appears in a table.
 *
 * Setup:
 * - Table T with allowed values
 * - Values V that should all be in T
 *
 * Claim: Every element of V appears in T
 *
 * Using multisets:
 * - Define f(b) = V[b] (the values to look up)
 * - Define g(b) = T[b] (the table, possibly with repetition)
 *
 * But wait - V and T might have different sizes!
 *
 * SOLUTION:
 * Extend both to the same size:
 * - Pad V with dummy lookups (arbitrary table values)
 * - Pad T with repetitions of table values
 *
 * Then multiset check proves the multiset relationship.
 *
 * COMPRESSED LOOKUP (PLOOKUP style):
 *
 * More efficient: use randomness to compress table + values.
 *
 * Define sorted polynomials s, t that combine T and V.
 * Use permutation check to verify sorting.
 *
 * This is how modern lookup arguments work!
 */

console.log("Lookup Argument:");
console.log("  Table T = allowed values");
console.log("  Values V = values to look up");
console.log("  Prove: every V[i] appears in T\n");

console.log("Approach: Multiset equality with padding");
console.log("  Extend V with dummy values from T");
console.log("  Extend T with repetitions");
console.log("  Multiset check proves lookup validity\n");

// Example: Lookup
// Table T = {1, 2, 3, 4} (the allowed values)
// Values V = {2, 4, 2, 1} (must all be in T)

console.log("Example:");
console.log("  Table T = {1, 2, 3, 4}");
console.log("  Lookup V = {2, 4, 2, 1}");
console.log("  All of V in T? YES\n");

// With repetition counting:
// V = {1, 2, 2, 4}
// T extended = {1, 2, 2, 3, 3, 4} ... needs more thought

console.log("Note: Real lookup arguments use more sophisticated");
console.log("techniques (PLOOKUP, LogUp, etc.) for efficiency.\n");

// ============================================================================
// Section 7: Implementation Stubs
// ============================================================================

console.log("--- Section 7: Implementation Stubs ---\n");

/**
 * Multiset equality check protocol.
 *
 * @throws NotImplementedError - Full implementation pending
 */
export function multisetCheck(
  f: HypercubeFunction,
  g: HypercubeFunction,
  beta: bigint,
  getChallenge: () => bigint
): { accepted: boolean } {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: multisetCheck - ' +
    'Full multiset check with ProductCheck integration'
  );
}

/**
 * Lookup argument: prove all values in V appear in table T.
 *
 * @throws NotImplementedError - Full lookup argument pending
 */
export function lookupArgument(
  table: bigint[],
  values: bigint[],
  getChallenge: () => bigint
): { accepted: boolean } {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: lookupArgument - ' +
    'Lookup argument using multiset/permutation checks'
  );
}

/**
 * PLOOKUP-style compressed lookup.
 *
 * @throws NotImplementedError - PLOOKUP implementation pending
 */
export function plookup(
  table: bigint[],
  values: bigint[],
  getChallenge: () => bigint
): { accepted: boolean } {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: plookup - ' +
    'PLOOKUP compressed lookup argument'
  );
}

console.log("Stub functions:");
console.log("  - multisetCheck: full multiset equality");
console.log("  - lookupArgument: basic lookup");
console.log("  - plookup: compressed lookup\n");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");

console.log("Multiset Equality Check:");
console.log("  Claim: {f(b)} = {g(b)} as multisets");
console.log("  Method: ProductCheck on (beta - f, beta - g)\n");

console.log("Key insight:");
console.log("  prod_b (X - f(b)) = prod_b (X - g(b)) as polynomials");
console.log("  iff multisets are equal\n");

console.log("Applications:");
console.log("  - Lookup arguments");
console.log("  - Copy constraints");
console.log("  - Memory consistency\n");

console.log("Next: Permutation checks for wiring");

// Exports
export { type HypercubeFunction };
