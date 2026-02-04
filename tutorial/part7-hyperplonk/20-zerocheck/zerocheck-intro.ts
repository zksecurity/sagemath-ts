/**
 * # ZeroCheck: Proving Polynomial Evaluates to Zero
 *
 * ZeroCheck is a key building block in HyperPlonk that allows the prover
 * to demonstrate that a polynomial evaluates to zero at ALL points on
 * the boolean hypercube {0,1}^n.
 *
 * ## The Problem
 *
 * Given a polynomial f: F^n -> F, prove that:
 *   f(b) = 0 for all b in {0,1}^n
 *
 * This is crucial for verifying circuit satisfaction:
 * - Gate constraints: C(x) = 0 for all wire assignments
 * - Lookup constraints: L(x) = 0 for all table entries
 *
 * ## The Solution: Reduce to SumCheck
 *
 * ZeroCheck reduces to SumCheck using a clever randomization technique.
 * Instead of checking 2^n zeros directly, we check a random combination.
 *
 * @module tutorial/part7-hyperplonk/20-zerocheck/zerocheck-intro
 */

import { NotImplementedError } from 'sagemath-ts';

console.log("=== ZeroCheck: Proving f = 0 on the Boolean Hypercube ===\n");

// ============================================================================
// Field Arithmetic (simple modular arithmetic)
// ============================================================================

const p = 97n;  // Small prime for demonstration

/** Modular reduction */
const mod = (x: bigint): bigint => ((x % p) + p) % p;

/** Field addition */
const add = (a: bigint, b: bigint): bigint => mod(a + b);

/** Field subtraction */
const sub = (a: bigint, b: bigint): bigint => mod(a - b);

/** Field multiplication */
const mul = (a: bigint, b: bigint): bigint => mod(a * b);

/**
 * Generates all points in {0,1}^n.
 */
function booleanHypercube(n: number): number[][] {
  if (n === 0) return [[]];
  const smaller = booleanHypercube(n - 1);
  return [...smaller.map(v => [0, ...v]), ...smaller.map(v => [1, ...v])];
}

// ============================================================================
// Section 1: The ZeroCheck Problem
// ============================================================================

console.log("--- Section 1: The ZeroCheck Problem ---\n");

/**
 * ZEROCHECK PROBLEM STATEMENT
 *
 * Given: Multilinear polynomial f: F^n -> F
 * Claim: f(b) = 0 for all b in {0,1}^n
 *
 * DIRECT APPROACH (NAIVE):
 *   Check f(b) for all 2^n points b in {0,1}^n.
 *   Cost: O(2^n) evaluations - infeasible!
 *
 * POLYNOMIAL IDENTITY:
 *   If f(b) = 0 for all b in B_n, then f must be divisible by
 *   the vanishing polynomial of B_n.
 *
 *   However, working with the vanishing polynomial directly
 *   is expensive (it has 2^n roots).
 *
 * THE INSIGHT:
 *   Use randomization to compress the check into a single SumCheck.
 */

console.log("ZeroCheck Problem:");
console.log("  Given: polynomial f(x_1, ..., x_n)");
console.log("  Prove: f(b) = 0 for ALL b in {0,1}^n\n");

console.log("Why it matters:");
console.log("  - Circuit satisfaction: C(witness) = 0 everywhere");
console.log("  - Gate constraints must hold for all gates");
console.log("  - Lookup constraints must hold for all entries\n");

console.log("Challenge: 2^n points to check - exponential!\n");

// ============================================================================
// Section 2: The Reduction to SumCheck
// ============================================================================

console.log("--- Section 2: Reduction to SumCheck ---\n");

/**
 * ZEROCHECK TO SUMCHECK REDUCTION
 *
 * Key observation: If f(b) = 0 for all b in B_n, then:
 *
 *   sum_{b in B_n} alpha^{|b|} * f(b) = 0
 *
 * for ANY alpha, where |b| = b_1 + b_2 + ... + b_n (Hamming weight).
 *
 * But this is just testing that the sum is zero.
 * We need something stronger that the verifier can check.
 *
 * THE ACTUAL REDUCTION:
 *
 * Let r = (r_1, ..., r_n) be verifier's random challenge.
 * Define:
 *
 *   g(x) = f(x) * eq(x, r)
 *
 * where eq(x, r) = prod_i (x_i * r_i + (1 - x_i)(1 - r_i))
 *
 * Claim: f(b) = 0 for all b in B_n
 *        <=> sum_{b in B_n} g(b) = f(r) * [sum_{b in B_n} eq(b, r)]
 *
 * Wait, that's not quite right either. Let's be more precise.
 *
 * CORRECT REDUCTION (HyperPlonk style):
 *
 * 1. Verifier sends random point r = (r_1, ..., r_n)
 *
 * 2. Define h(x) = eq(x, r) * f(x)
 *
 * 3. If f(b) = 0 for all b in B_n, then h(b) = 0 for all b in B_n
 *
 * 4. Therefore: sum_{b in B_n} h(b) = 0
 *
 * 5. Run SumCheck to verify sum_{b in B_n} h(b) = 0
 *
 * 6. SumCheck reduces to evaluating h(s) = eq(s, r) * f(s)
 *    at a random point s.
 *
 * 7. Verifier computes eq(s, r) and needs f(s) from an oracle.
 */

console.log("Reduction to SumCheck:\n");

console.log("Step 1: Verifier sends random r = (r_1, ..., r_n)");
console.log("Step 2: Define h(x) = eq(x, r) * f(x)");
console.log("Step 3: If f(b) = 0 for all b, then h(b) = 0 for all b");
console.log("Step 4: Therefore sum_{b in B_n} h(b) = 0");
console.log("Step 5: Run SumCheck to verify this sum");
console.log("Step 6: Final check: h(s) = eq(s, r) * f(s) at random s\n");

// ============================================================================
// Section 3: The eq(x, r) Function Revisited
// ============================================================================

console.log("--- Section 3: The eq(x, r) Function ---\n");

/**
 * THE EQUALITY POLYNOMIAL
 *
 * eq(x, r) = prod_{i=1}^{n} (x_i * r_i + (1 - x_i)(1 - r_i))
 *
 * Key properties:
 *
 * 1. For b in {0,1}^n and any r:
 *    eq(b, r) = prod_i [b_i * r_i + (1-b_i)(1-r_i)]
 *
 *    When b_i = 0: factor = (1 - r_i)
 *    When b_i = 1: factor = r_i
 *
 *    So: eq(b, r) = prod_{i: b_i=1} r_i * prod_{i: b_i=0} (1-r_i)
 *
 * 2. sum_{b in B_n} eq(b, r) = 1 for any r
 *    (This is because {eq(b, r)}_b forms a partition of unity)
 *
 * 3. eq(b, b') = delta_{b,b'} for b, b' in {0,1}^n
 *    (Selection property on boolean points)
 *
 * WHY eq MATTERS FOR ZEROCHECK:
 *
 * When we multiply f(x) by eq(x, r), we're creating a "random combination"
 * of the values f(b). If all f(b) = 0, the sum is 0.
 * If any f(b) != 0, the random combination is nonzero with high probability.
 */

/**
 * Computes eq(x, r) = prod_i (x_i * r_i + (1 - x_i)(1 - r_i)).
 */
function eq(x: bigint[], r: bigint[]): bigint {
  if (x.length !== r.length) {
    throw new Error(`Dimension mismatch: ${x.length} vs ${r.length}`);
  }

  let result = 1n;
  for (let i = 0; i < x.length; i++) {
    const xi = x[i]!;
    const ri = r[i]!;
    // x_i * r_i + (1 - x_i)(1 - r_i)
    // = x_i * r_i + 1 - x_i - r_i + x_i * r_i
    // = 2 * x_i * r_i - x_i - r_i + 1
    const factor = add(
      mul(xi, ri),
      mul(sub(1n, xi), sub(1n, ri))
    );
    result = mul(result, factor);
  }

  return result;
}

console.log("eq(x, r) = prod_i (x_i * r_i + (1 - x_i)(1 - r_i))\n");

// Demonstrate eq at boolean points
const r = [3n, 7n];  // Random challenge
console.log(`Random challenge r = (${r.join(', ')})\n`);

console.log("eq(b, r) for b in {0,1}^2:");
let eqSum = 0n;
for (const b of booleanHypercube(2)) {
  const bField = b.map(v => BigInt(v));
  const eqVal = eq(bField, r);
  console.log(`  eq((${b.join(',')}), r) = ${eqVal}`);
  eqSum = add(eqSum, eqVal);
}
console.log(`\nSum of eq values: ${eqSum} (should be 1)\n`);

// ============================================================================
// Section 4: Why the Reduction is Sound
// ============================================================================

console.log("--- Section 4: Soundness of the Reduction ---\n");

/**
 * SOUNDNESS ARGUMENT
 *
 * Claim: If f(b*) != 0 for some b* in B_n, the verifier rejects with
 *        high probability.
 *
 * Proof sketch:
 *
 * 1. The verifier's random r is chosen AFTER the prover commits to f.
 *
 * 2. If f(b*) != 0 for some b*, then:
 *    sum_{b in B_n} eq(b, r) * f(b) = sum_{b in B_n} eq(b, r) * f(b)
 *
 *    This is a random linear combination of the f(b) values,
 *    with coefficients eq(b, r) depending on r.
 *
 * 3. If f(b*) != 0, then the coefficient eq(b*, r) is a polynomial in r
 *    that is not identically zero (it equals 1 when r = b*).
 *
 * 4. By Schwartz-Zippel, the sum is nonzero with probability >= 1 - n/|F|.
 *
 * 5. If the sum is nonzero, SumCheck will catch the cheating prover.
 *
 * TOTAL SOUNDNESS ERROR:
 * ZeroCheck soundness error <= n/|F| (from Schwartz-Zippel)
 * + SumCheck soundness error <= d*n/|F| (for degree d)
 *
 * For multilinear f with 256-bit field: error < 2n/2^256 ~ negligible
 */

console.log("Soundness Analysis:");
console.log("  If f(b*) != 0 for some b* in B_n:");
console.log("    - sum_{b} eq(b, r) * f(b) is a random linear combination");
console.log("    - Coefficient eq(b*, r) is nonzero for most r");
console.log("    - Sum is nonzero with probability >= 1 - n/|F|");
console.log("    - SumCheck catches cheating prover\n");

console.log("Total error: O(n/|F|) - negligible for large fields!\n");

// ============================================================================
// Section 5: Example Walkthrough
// ============================================================================

console.log("--- Section 5: Concrete Example ---\n");

/**
 * EXAMPLE: Verify that f(x,y) = (x - x^2)(y - y^2) is zero on {0,1}^2
 *
 * Note: x - x^2 = x(1-x), which is 0 when x in {0,1}
 * Similarly for y.
 * So f(b) = 0 for all b in {0,1}^2.
 *
 * Let's trace through the ZeroCheck protocol.
 */

// Define f(x,y) = (x - x^2)(y - y^2) = x(1-x) * y(1-y)
function f(x: bigint, y: bigint): bigint {
  const xPart = mul(x, sub(1n, x));  // x(1-x)
  const yPart = mul(y, sub(1n, y));  // y(1-y)
  return mul(xPart, yPart);
}

console.log("Example: f(x,y) = x(1-x) * y(1-y)");
console.log("  This should be zero on all of {0,1}^2\n");

// Verify f is zero on boolean points
console.log("Checking f on {0,1}^2:");
for (const [bx, by] of booleanHypercube(2)) {
  const val = f(BigInt(bx), BigInt(by));
  console.log(`  f(${bx}, ${by}) = ${val}`);
}
console.log();

// ZeroCheck protocol
console.log("=== ZeroCheck Protocol ===\n");

// Step 1: Verifier sends random r
const rChallenge = [3n, 7n];
console.log(`Step 1: Verifier sends r = (${rChallenge.join(', ')})\n`);

// Step 2: Define h(x) = eq(x, r) * f(x)
console.log("Step 2: Define h(x, y) = eq((x,y), r) * f(x, y)\n");

// Step 3: Compute sum_{b in B_2} h(b)
console.log("Step 3: Compute sum_{b in B_2} h(b):");
let hSum = 0n;
for (const b of booleanHypercube(2)) {
  const bField = b.map(v => BigInt(v));
  const eqVal = eq(bField, rChallenge);
  const fVal = f(bField[0]!, bField[1]!);
  const hVal = mul(eqVal, fVal);
  console.log(`  h(${b.join(',')}) = eq * f = ${eqVal} * ${fVal} = ${hVal}`);
  hSum = add(hSum, hVal);
}
console.log(`\n  Sum = ${hSum} (should be 0)\n`);

// Step 4: If sum = 0, run SumCheck (claim is 0)
console.log("Step 4: Run SumCheck with claimed sum = 0");
console.log("  (SumCheck rounds would proceed here...)\n");

// Final evaluation
console.log("Step 5: Final check at SumCheck's random point s");
const s = [5n, 9n];  // Suppose SumCheck ends at s
console.log(`  Suppose SumCheck ends at s = (${s.join(', ')})`);

const eqAtS = eq(s, rChallenge);
const fAtS = f(s[0]!, s[1]!);
const hAtS = mul(eqAtS, fAtS);

console.log(`  eq(s, r) = ${eqAtS}`);
console.log(`  f(s) = ${fAtS}`);
console.log(`  h(s) = eq(s, r) * f(s) = ${hAtS}\n`);

// ============================================================================
// Section 6: ZeroCheck with Non-Zero Polynomial
// ============================================================================

console.log("--- Section 6: Detecting Non-Zero Polynomial ---\n");

/**
 * What happens if f is NOT zero everywhere?
 * Let's try g(x,y) = xy + 1, which is:
 *   g(0,0) = 1 != 0
 *   g(0,1) = 1
 *   g(1,0) = 1
 *   g(1,1) = 2
 */

function g(x: bigint, y: bigint): bigint {
  return add(mul(x, y), 1n);
}

console.log("Cheating example: g(x,y) = xy + 1");
console.log("  g(0,0) = 1, g(0,1) = 1, g(1,0) = 1, g(1,1) = 2");
console.log("  g is NOT zero on {0,1}^2!\n");

// Compute sum with ZeroCheck encoding
console.log("Computing sum_{b} eq(b, r) * g(b):");
let gSum = 0n;
for (const b of booleanHypercube(2)) {
  const bField = b.map(v => BigInt(v));
  const eqVal = eq(bField, rChallenge);
  const gVal = g(bField[0]!, bField[1]!);
  const product = mul(eqVal, gVal);
  console.log(`  eq(${b.join(',')}, r) * g(${b.join(',')}) = ${eqVal} * ${gVal} = ${product}`);
  gSum = add(gSum, product);
}
console.log(`\n  Sum = ${gSum} (NOT zero!)`);
console.log("  => Cheating detected: sum should be 0 if g = 0 everywhere\n");

// ============================================================================
// Section 7: Degree Analysis
// ============================================================================

console.log("--- Section 7: Degree Analysis ---\n");

/**
 * DEGREE OF h(x) = eq(x, r) * f(x)
 *
 * eq(x, r) is multilinear (degree 1 in each variable).
 *
 * If f is multilinear, then h has degree 2 in each variable.
 *
 * For SumCheck:
 * - Round polynomial g_i has degree 2
 * - Prover sends 3 field elements per round (g_i(0), g_i(1), g_i(2))
 *
 * If f has degree d in each variable:
 * - h has degree d+1 in each variable
 * - Prover sends d+2 elements per round
 */

console.log("Degree analysis:");
console.log("  eq(x, r) has degree 1 in each variable");
console.log("  If f has degree d, then h = eq * f has degree d+1\n");

console.log("For multilinear f:");
console.log("  h has degree 2 in each variable");
console.log("  Each round: prover sends 3 field elements");
console.log("  Total: 3n field elements for ZeroCheck\n");

// ============================================================================
// Section 8: Preview of ZeroCheck Protocol
// ============================================================================

console.log("--- Section 8: Full Protocol Structure ---\n");

/**
 * COMPLETE ZEROCHECK PROTOCOL
 *
 * Public: Polynomial f (commitment), number of variables n
 * Claim: f(b) = 0 for all b in {0,1}^n
 *
 * PROTOCOL:
 *
 * 1. V -> P: Random challenge r = (r_1, ..., r_n)
 *
 * 2. P -> V: Run SumCheck on h(x) = eq(x, r) * f(x)
 *            with claimed sum = 0
 *
 *    This produces:
 *    - n round polynomials g_1, ..., g_n
 *    - Verifier challenges s_1, ..., s_n
 *    - Claimed final value h(s)
 *
 * 3. V: Verify SumCheck transcript
 *    - Check each g_i(0) + g_i(1) = previous claim
 *    - Final claim should equal h(s)
 *
 * 4. V: Check final evaluation
 *    - Compute eq(s, r) (can do locally)
 *    - Query oracle for f(s)
 *    - Verify h(s) = eq(s, r) * f(s)
 *
 * ACCEPT iff all checks pass.
 */

console.log("ZeroCheck Protocol:");
console.log("  1. Verifier sends random r");
console.log("  2. Prover runs SumCheck on h(x) = eq(x,r) * f(x) with sum = 0");
console.log("  3. Verifier checks SumCheck transcript");
console.log("  4. Verifier checks h(s) = eq(s,r) * f(s) using oracle\n");

console.log("Properties:");
console.log("  - Completeness: Honest prover always convinces verifier");
console.log("  - Soundness: Cheating prover caught with prob >= 1 - O(n/|F|)");
console.log("  - Communication: O(n) field elements\n");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");

console.log("ZeroCheck Protocol:");
console.log("  Claim: f(b) = 0 for all b in {0,1}^n");
console.log("  Reduction: sum_{b} eq(b, r) * f(b) = 0 iff f = 0 everywhere\n");

console.log("Key function eq(x, r):");
console.log("  eq(x, r) = prod_i (x_i * r_i + (1-x_i)(1-r_i))");
console.log("  Creates random linear combination of f values\n");

console.log("Soundness:");
console.log("  If f(b*) != 0 for some b*, sum is nonzero whp");
console.log("  Error probability: O(n/|F|)\n");

console.log("Next: zerocheck-protocol.ts - Full implementation");

// Export for other tutorials
export { eq, booleanHypercube };
