/**
 * # ProductCheck: Verifying Product Relations
 *
 * ProductCheck is a protocol that allows the prover to demonstrate that
 * the product of polynomial evaluations over the boolean hypercube equals
 * a claimed value (typically 1).
 *
 * ## The Problem
 *
 * Given polynomials f and g, prove that:
 *   prod_{b in {0,1}^n} f(b) = prod_{b in {0,1}^n} g(b)
 *
 * Or equivalently (when g(b) != 0):
 *   prod_{b in {0,1}^n} f(b)/g(b) = 1
 *
 * ## Why ProductCheck Matters
 *
 * ProductCheck is essential for:
 * - Permutation arguments (HyperPlonk wiring)
 * - Lookup arguments
 * - Grand product arguments
 *
 * ## The Approach
 *
 * ProductCheck reduces to ZeroCheck using a "fractional sum" technique.
 *
 * @module tutorial/part7-hyperplonk/21-productcheck/productcheck-intro
 */

import { NotImplementedError } from 'sagemath-ts';

console.log("=== ProductCheck: Verifying Product Relations ===\n");

// ============================================================================
// Field Arithmetic (simple modular arithmetic)
// ============================================================================

const p = 97n;

const mod = (x: bigint): bigint => ((x % p) + p) % p;
const add = (a: bigint, b: bigint): bigint => mod(a + b);
const sub = (a: bigint, b: bigint): bigint => mod(a - b);
const mul = (a: bigint, b: bigint): bigint => mod(a * b);

function powMod(base: bigint, exp: bigint, modulus: bigint): bigint {
  let result = 1n;
  base = ((base % modulus) + modulus) % modulus;
  while (exp > 0n) {
    if (exp % 2n === 1n) result = (result * base) % modulus;
    exp = exp / 2n;
    base = (base * base) % modulus;
  }
  return result;
}

const inv = (a: bigint): bigint => powMod(a, p - 2n, p);

/**
 * Generates all points in {0,1}^n.
 */
function booleanHypercube(n: number): number[][] {
  if (n === 0) return [[]];
  const smaller = booleanHypercube(n - 1);
  return [...smaller.map(v => [0, ...v]), ...smaller.map(v => [1, ...v])];
}

// ============================================================================
// Section 1: The ProductCheck Problem
// ============================================================================

console.log("--- Section 1: The ProductCheck Problem ---\n");

/**
 * PRODUCTCHECK PROBLEM
 *
 * Given: Multilinear polynomials f, g: F^n -> F
 * Claim: prod_{b in B_n} f(b) = prod_{b in B_n} g(b)
 *
 * Or equivalently, if g(b) != 0 for all b:
 * Claim: prod_{b in B_n} f(b)/g(b) = 1
 *
 * WHY THIS IS HARD:
 * - The product has 2^n terms
 * - We can't just check each term individually
 * - Need a clever reduction
 *
 * APPLICATIONS:
 *
 * 1. PERMUTATION CHECK:
 *    If {f(b)} and {g(b)} are permutations of each other,
 *    then for random beta:
 *    prod_b (f(b) + beta) = prod_b (g(b) + beta)
 *
 * 2. MULTISET EQUALITY:
 *    Same as permutation, but with multiplicities
 *
 * 3. LOOKUP ARGUMENTS:
 *    Prove that all values in table appear in lookup
 */

console.log("ProductCheck Problem:");
console.log("  Claim: prod_{b in B_n} f(b) = prod_{b in B_n} g(b)");
console.log("  Or: prod_{b in B_n} f(b)/g(b) = 1\n");

console.log("Applications:");
console.log("  1. Permutation arguments (wiring in PLONK)");
console.log("  2. Multiset equality checks");
console.log("  3. Lookup arguments\n");

// ============================================================================
// Section 2: The Running Product Approach
// ============================================================================

console.log("--- Section 2: The Running Product Approach ---\n");

/**
 * RUNNING PRODUCT TECHNIQUE
 *
 * The key insight is to use an AUXILIARY polynomial that tracks
 * the running product.
 *
 * Define: h(b) = prod_{b' <= b} (f(b') / g(b'))
 *
 * where b' <= b means b' comes before b in lexicographic order.
 *
 * Properties:
 * - h(first point) = f(first point) / g(first point)
 * - h(b) = h(prev(b)) * (f(b) / g(b))
 * - h(last point) = prod_{all b} (f(b) / g(b))
 *
 * The claim prod f/g = 1 is equivalent to:
 *   h(last point) = 1
 *
 * But we also need to verify h is computed correctly!
 *
 * RECURRENCE RELATION:
 *   h(b) = h(prev(b)) * f(b) / g(b)
 *
 * Rearranging:
 *   h(b) * g(b) = h(prev(b)) * f(b)
 *
 * Or:
 *   h(b) * g(b) - h(prev(b)) * f(b) = 0
 *
 * This is a constraint that must hold for all b (except the first).
 */

console.log("Running Product Technique:\n");

console.log("Define h(b) = prod_{b' <= b} (f(b') / g(b'))");
console.log("  - h tracks the running product");
console.log("  - h(last point) = full product\n");

console.log("Recurrence relation:");
console.log("  h(b) * g(b) = h(prev(b)) * f(b)");
console.log("  For the first point: h(first) = f(first) / g(first)\n");

console.log("ProductCheck claim equivalent to:");
console.log("  1. h satisfies recurrence at all points");
console.log("  2. h(last point) = 1\n");

// ============================================================================
// Section 3: Concrete Example
// ============================================================================

console.log("--- Section 3: Concrete Example ---\n");

/**
 * EXAMPLE:
 *
 * Let n = 2, so B_2 = {00, 01, 10, 11}
 *
 * f(00) = 2, f(01) = 3, f(10) = 4, f(11) = 5  -> prod = 120
 * g(00) = 4, g(01) = 5, g(10) = 2, g(11) = 3  -> prod = 120
 *
 * Products are equal!
 *
 * Running product h:
 * h(00) = f(00)/g(00) = 2/4 = 1/2
 * h(01) = h(00) * f(01)/g(01) = (1/2) * (3/5) = 3/10
 * h(10) = h(01) * f(10)/g(10) = (3/10) * (4/2) = 3/5
 * h(11) = h(10) * f(11)/g(11) = (3/5) * (5/3) = 1
 *
 * Final value is 1, so products are equal!
 */

interface HypercubeFunction {
  numVars: number;
  values: bigint[];
}

// f(00)=2, f(01)=3, f(10)=4, f(11)=5
const fFunc: HypercubeFunction = {
  numVars: 2,
  values: [2n, 3n, 4n, 5n]
};

// g(00)=4, g(01)=5, g(10)=2, g(11)=3
const gFunc: HypercubeFunction = {
  numVars: 2,
  values: [4n, 5n, 2n, 3n]
};

console.log("Example on {0,1}^2:\n");

console.log("f values:");
const cube = booleanHypercube(2);
let fProd = 1n;
for (let i = 0; i < cube.length; i++) {
  console.log(`  f(${cube[i]!.join(',')}) = ${fFunc.values[i]!}`);
  fProd = mul(fProd, fFunc.values[i]!);
}
console.log(`  Product = ${fProd}\n`);

console.log("g values:");
let gProd = 1n;
for (let i = 0; i < cube.length; i++) {
  console.log(`  g(${cube[i]!.join(',')}) = ${gFunc.values[i]!}`);
  gProd = mul(gProd, gFunc.values[i]!);
}
console.log(`  Product = ${gProd}\n`);

console.log(`Products equal? ${fProd === gProd ? 'YES' : 'NO'}\n`);

// Compute running product
console.log("Running product h = f/g:");
let runningProd = 1n;
const hValues: bigint[] = [];
for (let i = 0; i < cube.length; i++) {
  const ratio = mul(fFunc.values[i]!, inv(gFunc.values[i]!));
  runningProd = mul(runningProd, ratio);
  hValues.push(runningProd);
  console.log(`  h(${cube[i]!.join(',')}) = ${runningProd}`);
}
console.log(`\nFinal h value: ${runningProd} (should be 1 if products equal)\n`);

// ============================================================================
// Section 4: The GKR-Style ProductCheck
// ============================================================================

console.log("--- Section 4: GKR-Style ProductCheck ---\n");

/**
 * GKR-STYLE PRODUCTCHECK
 *
 * Instead of tracking running product directly, we use a more algebraic approach.
 *
 * KEY OBSERVATION:
 * prod_b (f(b)/g(b)) = 1
 * <=> prod_b f(b) = prod_b g(b)
 *
 * LOGARITHMIC APPROACH (conceptual):
 * If we could take logs: sum_b log(f(b)) = sum_b log(g(b))
 *
 * But we're in a finite field, no logarithms!
 *
 * POLYNOMIAL COMMITMENT APPROACH:
 *
 * 1. Prover computes h(b) for all b (running product)
 *
 * 2. Prover commits to h
 *
 * 3. Verifier checks constraint:
 *    h(b) * g(b) - h(prev(b)) * f(b) = 0 for all b > first
 *
 *    This is a ZeroCheck on a polynomial!
 *
 * 4. Verifier also checks boundary conditions:
 *    - h(first) = f(first) / g(first)
 *    - h(last) = 1
 *
 * COMPLICATION: prev(b) for multilinear indexing
 *
 * In the boolean hypercube with lexicographic order:
 * prev(00...01) = 00...00
 * prev(00...10) = 00...01
 * ...
 * prev(10...00) = 01...11
 *
 * This is essentially binary subtraction by 1!
 */

console.log("GKR-Style approach:");
console.log("  1. Prover commits to auxiliary polynomial h");
console.log("  2. Verify recurrence: h(b)*g(b) = h(prev(b))*f(b)");
console.log("  3. This becomes a ZeroCheck!\n");

console.log("Constraint polynomial:");
console.log("  C(b) = h(b)*g(b) - h(prev(b))*f(b)");
console.log("  Must prove C(b) = 0 for all b (except first)\n");

// ============================================================================
// Section 5: The HyperPlonk ProductCheck
// ============================================================================

console.log("--- Section 5: HyperPlonk ProductCheck ---\n");

/**
 * HYPERPLONK'S PRODUCTCHECK (Simplified)
 *
 * HyperPlonk uses a different formulation based on the paper's approach.
 *
 * DEFINITION:
 * For f: B_n -> F*, define the product:
 *   P_f = prod_{b in B_n} f(b)
 *
 * FRACTIONAL SUM APPROACH:
 *
 * Key insight: prod_b f(b) = prod_b g(b) iff:
 *   sum_b log(f(b)) = sum_b log(g(b))
 *
 * But in finite fields, we use "formal derivatives":
 *   sum_b (f(b)'/f(b)) where f' is derivative
 *
 * Actually, HyperPlonk uses a RUNNING PRODUCT polynomial h.
 *
 * FORMAL PROTOCOL:
 *
 * Claim: prod_b f(b) = prod_b g(b)
 *
 * 1. Prover constructs h where h(first) = 1 and:
 *    h(next(b)) = h(b) * f(b) / g(b)
 *
 * 2. This implies h(after-last) = prod_b f(b)/g(b)
 *
 * 3. Claim is equivalent to h(after-last) = 1
 *
 * 4. Verifier checks:
 *    a) h(first) = 1
 *    b) The recurrence holds at all points (ZeroCheck)
 *    c) The claimed final product is correct
 *
 * ENCODING THE RECURRENCE:
 *
 * For b in B_n, next(b) in B_n is the successor in lex order.
 * The last point wraps around (or we handle boundary separately).
 *
 * The constraint:
 *   h(next(b)) * g(b) - h(b) * f(b) = 0
 *
 * Can be encoded using the "shift" operation on multilinear polynomials.
 */

console.log("HyperPlonk ProductCheck:");
console.log("  Define h with h(first) = 1");
console.log("  Recurrence: h(next(b)) = h(b) * f(b) / g(b)\n");

console.log("Implies: h(after-last) = prod_b (f(b)/g(b))");
console.log("Claim: h(after-last) = 1\n");

console.log("Verification:");
console.log("  1. Check h(first) = 1");
console.log("  2. ZeroCheck: h(next(b))*g(b) - h(b)*f(b) = 0");
console.log("  3. Check final boundary\n");

// ============================================================================
// Section 6: The Shift Polynomial
// ============================================================================

console.log("--- Section 6: The Shift Polynomial ---\n");

/**
 * THE SHIFT OPERATION
 *
 * For multilinear polynomials on B_n, we need to express:
 *   h(next(b))
 * in terms of a polynomial.
 *
 * LEXICOGRAPHIC ORDER ON B_n:
 * For n=2: 00 -> 01 -> 10 -> 11
 * For n=3: 000 -> 001 -> 010 -> 011 -> 100 -> 101 -> 110 -> 111
 *
 * next(b) = b + 1 (as binary number)
 *
 * ENCODING NEXT():
 *
 * For x in B_n, let x = (x_1, ..., x_n) with x_1 most significant.
 * next(x) = x + 1 mod 2^n
 *
 * This can be expressed as:
 *   next(x)_n = 1 - x_n (flip last bit)
 *   next(x)_{n-1} = x_{n-1} XOR x_n (carry)
 *   ...
 *
 * For polynomials, we use:
 *   h_shifted(x) = sum_b h(next(b)) * eq(x, b)
 *                = sum_b h(b) * eq(x, prev(b))
 *
 * This requires careful bookkeeping but is doable!
 *
 * MERGE POLYNOMIAL:
 * An alternative is to use the "merge" operation from HyperPlonk
 * which combines two polynomials into one with an extra variable.
 */

console.log("Shift operation on B_n:");
console.log("  next(b) = b + 1 (binary addition mod 2^n)");
console.log("  For n=2: 00 -> 01 -> 10 -> 11 -> 00\n");

console.log("Polynomial encoding:");
console.log("  h_shifted(x) = sum_b h(next(b)) * eq(x, b)");
console.log("  This 'shifts' the function by one position\n");

// Demonstrate shift
console.log("Shift on our running product h:");
for (let i = 0; i < cube.length; i++) {
  const nextIdx = (i + 1) % cube.length;
  console.log(`  h(${cube[i]!.join(',')}) = ${hValues[i]!}, h(next) = ${hValues[nextIdx]!}`);
}
console.log();

// ============================================================================
// Section 7: Reducing to ZeroCheck
// ============================================================================

console.log("--- Section 7: Reduction to ZeroCheck ---\n");

/**
 * PRODUCTCHECK TO ZEROCHECK REDUCTION
 *
 * Given f, g, and auxiliary h, the ProductCheck constraints are:
 *
 * 1. BOUNDARY: h(first) = 1
 *
 * 2. RECURRENCE: For all b except last:
 *    h(next(b)) * g(b) - h(b) * f(b) = 0
 *
 * 3. FINAL: h(last) * g(last) = 1 * f(last)
 *    (because next(last) wraps to first where h = 1)
 *
 * The recurrence constraint is perfect for ZeroCheck!
 *
 * Define: C(b) = h(next(b)) * g(b) - h(b) * f(b)
 *
 * ProductCheck reduces to:
 *   ZeroCheck(C) where C(b) = 0 for all b in B_n
 *
 * Plus boundary checks:
 *   - Check h(0...0) = 1
 *   - The wrap-around at the end is handled by the constraint
 *
 * SUBTLETY: The polynomial C involves h_shifted, which requires
 * some care in construction.
 */

console.log("Constraint polynomial:");
console.log("  C(b) = h(next(b)) * g(b) - h(b) * f(b)\n");

console.log("ProductCheck = ZeroCheck(C) + boundary checks\n");

// Verify constraint for our example
console.log("Verifying constraint on our example:");
for (let i = 0; i < cube.length; i++) {
  const nextIdx = (i + 1) % cube.length;
  const hNext = hValues[nextIdx]!;
  const hCurr = hValues[i]!;
  const fCurr = fFunc.values[i]!;
  const gCurr = gFunc.values[i]!;

  // C(b) = h(next(b)) * g(b) - h(b) * f(b)
  const constraint = sub(mul(hNext, gCurr), mul(hCurr, fCurr));
  console.log(`  C(${cube[i]!.join(',')}) = h(next)*g - h*f = ${hNext}*${gCurr} - ${hCurr}*${fCurr} = ${constraint}`);
}
console.log();

// Note: For the last point, next wraps to first where h = 1
// So the constraint at last point checks: 1 * g(last) - h(last) * f(last) = 0
// i.e., g(last) = h(last) * f(last)

// ============================================================================
// Section 8: Summary and Preview
// ============================================================================

console.log("--- Section 8: Summary ---\n");

/**
 * PRODUCTCHECK SUMMARY
 *
 * Problem: Prove prod_b f(b) = prod_b g(b)
 *
 * Approach:
 * 1. Construct auxiliary h tracking running product f/g
 * 2. Constraint: h(next(b)) * g(b) = h(b) * f(b)
 * 3. Reduce constraint to ZeroCheck
 * 4. Check boundary: h(first) = 1
 *
 * Complexity:
 * - Prover: O(2^n) to compute h + ZeroCheck prover
 * - Verifier: ZeroCheck verifier + boundary checks
 * - Communication: O(n) field elements
 *
 * Applications:
 * - Permutation checks (HyperPlonk wiring)
 * - Multiset equality (lookup arguments)
 * - Grand products
 */

console.log("ProductCheck Summary:\n");

console.log("Claim: prod_b f(b) = prod_b g(b)");
console.log("Approach: Running product + ZeroCheck\n");

console.log("Steps:");
console.log("  1. Construct h with h(next(b)) = h(b) * f(b)/g(b)");
console.log("  2. Define C(b) = h(next(b))*g(b) - h(b)*f(b)");
console.log("  3. Run ZeroCheck on C");
console.log("  4. Check h(first) = 1\n");

console.log("Properties:");
console.log("  - Sound: if products differ, ZeroCheck fails");
console.log("  - Complete: correct h satisfies all constraints");
console.log("  - Efficient: O(n) communication\n");

console.log("Next: productcheck-protocol.ts - Full implementation");

// Export types for other tutorials
export { type HypercubeFunction };
