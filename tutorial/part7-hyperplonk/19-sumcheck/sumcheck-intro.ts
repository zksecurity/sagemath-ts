/**
 * # Introduction to the SumCheck Protocol
 *
 * The SumCheck protocol is one of the most fundamental building blocks in
 * modern zero-knowledge proof systems. It allows a prover to convince a
 * verifier about the sum of a multivariate polynomial over a boolean hypercube.
 *
 * ## Why SumCheck Matters
 *
 * SumCheck is the backbone of:
 * - HyperPlonk (this tutorial series)
 * - GKR protocol
 * - Spartan
 * - Lasso/Jolt
 * - Many lookup arguments
 *
 * ## The Problem Statement
 *
 * Given a multivariate polynomial f(x_1, x_2, ..., x_n) over a field F,
 * the prover wants to convince the verifier of the claimed sum:
 *
 *   H = sum_{(b_1,...,b_n) in {0,1}^n} f(b_1, b_2, ..., b_n)
 *
 * This sum has 2^n terms, which is exponential in n!
 *
 * ## The Key Insight
 *
 * The verifier CANNOT check all 2^n evaluations (that would be expensive).
 * Instead, SumCheck reduces the claim to a SINGLE random evaluation point.
 *
 * After n rounds, the verifier only needs to evaluate f at ONE point:
 *   f(r_1, r_2, ..., r_n)
 *
 * This reduction is what makes SumCheck so powerful.
 *
 * @module tutorial/part7-hyperplonk/19-sumcheck/sumcheck-intro
 */

console.log("=== Introduction to the SumCheck Protocol ===\n");

// ============================================================================
// Field Arithmetic Setup
// ============================================================================

// Simple modular arithmetic for demonstrations
const p = 97n;  // Small prime for readable output
const mod = (x: bigint): bigint => ((x % p) + p) % p;
const add = (a: bigint, b: bigint): bigint => mod(a + b);
const sub = (a: bigint, b: bigint): bigint => mod(a - b);
const mul = (a: bigint, b: bigint): bigint => mod(a * b);

// ============================================================================
// Section 1: The Boolean Hypercube
// ============================================================================

console.log("--- Section 1: The Boolean Hypercube B_n = {0,1}^n ---\n");

/**
 * THE BOOLEAN HYPERCUBE
 *
 * B_n = {0, 1}^n is the set of all n-bit binary strings.
 *
 * For n = 2: B_2 = {(0,0), (0,1), (1,0), (1,1)} - 4 points
 * For n = 3: B_3 = {(0,0,0), (0,0,1), ..., (1,1,1)} - 8 points
 *
 * |B_n| = 2^n (exponential in n)
 *
 * The boolean hypercube is fundamental because:
 * 1. It represents all possible variable assignments
 * 2. Multilinear polynomials are uniquely determined by values on B_n
 * 3. Circuit satisfiability can be encoded as sums over B_n
 */

console.log("Boolean Hypercube B_n = {0, 1}^n");
console.log("  |B_n| = 2^n points\n");

// Let's enumerate B_2 and B_3
function generateHypercube(n: number): number[][] {
  if (n === 0) return [[]];
  const smaller = generateHypercube(n - 1);
  return [...smaller.map(v => [0, ...v]), ...smaller.map(v => [1, ...v])];
}

const B_2 = generateHypercube(2);
console.log("B_2 (4 points):");
B_2.forEach(pt => console.log(`  (${pt.join(', ')})`));

const B_3 = generateHypercube(3);
console.log(`\nB_3 (${B_3.length} points):`);
console.log(`  ${B_3.map(pt => `(${pt.join(',')})`).join(', ')}\n`);

// ============================================================================
// Section 2: Multilinear Polynomials
// ============================================================================

console.log("--- Section 2: Multilinear Polynomials ---\n");

/**
 * MULTILINEAR POLYNOMIALS
 *
 * A polynomial f(x_1, ..., x_n) is MULTILINEAR if each variable
 * appears with degree at most 1 in every term.
 *
 * Examples:
 *   f(x, y) = 3xy + 2x + y + 5  -- multilinear
 *   g(x, y) = x^2 + y           -- NOT multilinear (x has degree 2)
 *   h(x, y, z) = xyz + x + z    -- multilinear
 *
 * UNIQUE EXTENSION PROPERTY:
 *
 * Given ANY function f: B_n -> F (values on the boolean hypercube),
 * there is a UNIQUE multilinear polynomial that agrees with f on B_n.
 *
 * This is because:
 * - A multilinear polynomial in n variables has at most 2^n terms
 * - We have exactly 2^n constraints (one for each point in B_n)
 * - The system is uniquely determined
 *
 * This unique polynomial is called the MULTILINEAR EXTENSION (MLE).
 */

console.log("Multilinear Polynomial: each variable has degree <= 1");
console.log("  Valid:   f(x,y) = 3xy + 2x + y + 5");
console.log("  Invalid: g(x,y) = x^2 + y\n");

console.log("Key Property: UNIQUE MULTILINEAR EXTENSION");
console.log("  Any function f: B_n -> F uniquely extends to a");
console.log("  multilinear polynomial over all of F^n.\n");

// Define a multilinear polynomial f(x, y) = 2xy + 3x + 5y + 7
// We'll work in F_97 for readable output
function evalF(x: bigint, y: bigint): bigint {
  // f(x, y) = 2xy + 3x + 5y + 7 (mod 97)
  return add(add(add(mul(2n, mul(x, y)), mul(3n, x)), mul(5n, y)), 7n);
}

console.log(`Working over F_${p}`);
console.log(`f(x, y) = 2xy + 3x + 5y + 7`);

// Evaluate at points on B_2
console.log("\nEvaluations on B_2:");
for (const [b1, b2] of B_2) {
  const val = evalF(BigInt(b1), BigInt(b2));
  console.log(`  f(${b1}, ${b2}) = ${val}`);
}

// ============================================================================
// Section 3: The SumCheck Problem
// ============================================================================

console.log("\n--- Section 3: The SumCheck Problem ---\n");

/**
 * THE SUMCHECK CLAIM
 *
 * Prover claims: H = sum_{b in B_n} f(b)
 *
 * Naive verification: Check all 2^n evaluations and add them up.
 *   - Time: O(2^n * cost to evaluate f)
 *   - Completely impractical for large n!
 *
 * EXAMPLE:
 * For n = 100 variables, the hypercube has 2^100 points.
 * Even at 10^15 evaluations per second, this would take
 * more than 10^14 years.
 *
 * THE SUMCHECK SOLUTION:
 * Reduce verification to:
 *   - n rounds of interaction
 *   - O(n) total verifier work (excluding final evaluation)
 *   - ONE evaluation of f at a random point
 *
 * The verifier's work is LINEAR in n, not exponential!
 */

console.log("The SumCheck Problem:");
console.log("  Claim: H = sum_{b in B_n} f(b)");
console.log("  Naive check: O(2^n) evaluations - EXPONENTIAL!");
console.log("  SumCheck: O(n) rounds + 1 evaluation - EFFICIENT!\n");

// Compute the sum over B_2 for our example
let sumOverB2 = 0n;
for (const [b1, b2] of B_2) {
  const val = evalF(BigInt(b1), BigInt(b2));
  sumOverB2 = add(sumOverB2, val);
}

console.log(`For our example f(x,y) = 2xy + 3x + 5y + 7:`);
console.log(`  H = sum_{(b1,b2) in B_2} f(b1, b2) = ${sumOverB2}\n`);

// Verify: f(0,0) + f(0,1) + f(1,0) + f(1,1)
// = 7 + (7+5) + (7+3) + (7+3+5+2) = 7 + 12 + 10 + 17 = 46
console.log("Verification:");
console.log("  f(0,0) = 7");
console.log("  f(0,1) = 7 + 5 = 12");
console.log("  f(1,0) = 7 + 3 = 10");
console.log("  f(1,1) = 7 + 3 + 5 + 2 = 17");
console.log("  Sum = 7 + 12 + 10 + 17 = 46\n");

// ============================================================================
// Section 4: The Protocol Intuition
// ============================================================================

console.log("--- Section 4: Protocol Intuition ---\n");

/**
 * SUMCHECK PROTOCOL INTUITION
 *
 * The key idea is to "peel off" one variable at a time.
 *
 * Original claim: H = sum_{x1 in {0,1}} sum_{x2 in {0,1}} ... f(x1, x2, ...)
 *
 * ROUND 1:
 * Define g_1(X_1) = sum_{x2,...,xn in {0,1}} f(X_1, x2, ..., xn)
 *
 * Note: g_1(X_1) is a polynomial in one variable X_1.
 * Its degree in X_1 equals f's degree in x_1.
 *
 * Crucially: H = g_1(0) + g_1(1)
 *
 * The prover sends g_1(X_1) to the verifier.
 * The verifier:
 *   1. Checks that g_1(0) + g_1(1) = H  (the claimed sum)
 *   2. Picks random r_1 and sends it to prover
 *   3. New claim: g_1(r_1) = sum_{x2,...,xn} f(r_1, x2, ..., xn)
 *
 * The verifier has reduced a sum over 2^n points to a sum over 2^{n-1} points!
 *
 * ROUNDS 2 through n:
 * Repeat: prover sends a univariate polynomial, verifier checks and sends challenge.
 *
 * FINAL CHECK:
 * After n rounds, the claim is: g_n(r_n) = f(r_1, r_2, ..., r_n)
 * The verifier evaluates f at this random point (or uses an oracle).
 */

console.log("Protocol Idea: Peel off one variable per round\n");

console.log("Round 1:");
console.log("  g_1(X) = sum_{x2,...,xn in {0,1}} f(X, x2, ..., xn)");
console.log("  Prover sends g_1");
console.log("  Verifier checks: g_1(0) + g_1(1) = H");
console.log("  Verifier sends random r_1\n");

console.log("Round i:");
console.log("  g_i(X) = sum_{x_{i+1},...,xn in {0,1}} f(r_1,...,r_{i-1}, X, x_{i+1},...,xn)");
console.log("  Prover sends g_i");
console.log("  Verifier checks: g_i(0) + g_i(1) = g_{i-1}(r_{i-1})");
console.log("  Verifier sends random r_i\n");

console.log("Final Check:");
console.log("  Verifier checks: g_n(r_n) = f(r_1, r_2, ..., r_n)");
console.log("  (evaluates f at ONE random point)\n");

// ============================================================================
// Section 5: Concrete Example
// ============================================================================

console.log("--- Section 5: Worked Example ---\n");

/**
 * WORKED EXAMPLE WITH f(x,y) = 2xy + 3x + 5y + 7
 *
 * Claim: H = 46 (computed above)
 *
 * ROUND 1:
 * g_1(X) = f(X, 0) + f(X, 1)
 *        = (2*X*0 + 3*X + 5*0 + 7) + (2*X*1 + 3*X + 5*1 + 7)
 *        = (3X + 7) + (2X + 3X + 5 + 7)
 *        = (3X + 7) + (5X + 12)
 *        = 8X + 19
 *
 * Check: g_1(0) + g_1(1) = 19 + 27 = 46 = H  PASS
 *
 * Verifier sends random r_1 (say, r_1 = 2)
 *
 * ROUND 2:
 * Now we need: g_2(Y) = f(r_1, Y) = f(2, Y)
 *            = 2*2*Y + 3*2 + 5*Y + 7
 *            = 4Y + 6 + 5Y + 7
 *            = 9Y + 13
 *
 * Check: g_2(0) + g_2(1) = 13 + 22 = 35 = g_1(2) = 8*2 + 19 = 35  PASS
 *
 * Verifier sends random r_2 (say, r_2 = 3)
 *
 * FINAL CHECK:
 * g_2(r_2) = g_2(3) = 9*3 + 13 = 40
 * f(r_1, r_2) = f(2, 3) = 2*2*3 + 3*2 + 5*3 + 7 = 12 + 6 + 15 + 7 = 40  PASS
 */

console.log("Claim: H = 46 for f(x,y) = 2xy + 3x + 5y + 7\n");

// Round 1: g_1(X) = sum_{y in {0,1}} f(X, y) = f(X,0) + f(X,1)
console.log("ROUND 1:");
console.log("  g_1(X) = f(X, 0) + f(X, 1)");
console.log("         = (3X + 7) + (5X + 12)");
console.log("         = 8X + 19\n");

// g_1(X) = 8X + 19
function evalG1(X: bigint): bigint {
  return add(mul(8n, X), 19n);
}

const g1_0 = evalG1(0n);
const g1_1 = evalG1(1n);
console.log(`  g_1(0) = ${g1_0}`);
console.log(`  g_1(1) = ${g1_1}`);
console.log(`  g_1(0) + g_1(1) = ${add(g1_0, g1_1)} (should equal ${sumOverB2})\n`);

// Verifier sends r_1 = 2
const r_1 = 2n;
console.log(`  Verifier sends r_1 = ${r_1}`);
const g1_r1 = evalG1(r_1);
console.log(`  New claim: g_1(r_1) = ${g1_r1}\n`);

// Round 2: g_2(Y) = f(r_1, Y)
console.log("ROUND 2:");
console.log("  g_2(Y) = f(2, Y) = 4Y + 6 + 5Y + 7 = 9Y + 13\n");

// g_2(Y) = 9Y + 13
function evalG2(Y: bigint): bigint {
  return add(mul(9n, Y), 13n);
}

const g2_0 = evalG2(0n);
const g2_1 = evalG2(1n);
console.log(`  g_2(0) = ${g2_0}`);
console.log(`  g_2(1) = ${g2_1}`);
console.log(`  g_2(0) + g_2(1) = ${add(g2_0, g2_1)} (should equal g_1(${r_1}) = ${g1_r1})\n`);

// Verifier sends r_2 = 3
const r_2 = 3n;
console.log(`  Verifier sends r_2 = ${r_2}\n`);

// Final check
console.log("FINAL CHECK:");
const g2_r2 = evalG2(r_2);
console.log(`  g_2(r_2) = g_2(${r_2}) = ${g2_r2}`);

// f(r_1, r_2) = f(2, 3)
const f_at_r = evalF(r_1, r_2);
console.log(`  f(r_1, r_2) = f(${r_1}, ${r_2}) = ${f_at_r}`);
console.log(`  Match: ${g2_r2 === f_at_r ? 'YES' : 'NO'}\n`);

// ============================================================================
// Section 6: Security Guarantees
// ============================================================================

console.log("--- Section 6: Security Analysis ---\n");

/**
 * SECURITY OF SUMCHECK
 *
 * COMPLETENESS:
 * If the prover is honest and H is correct, the verifier always accepts.
 *
 * SOUNDNESS:
 * If H is incorrect, the verifier rejects with high probability.
 *
 * Why soundness holds:
 * - In each round, the prover commits to a univariate polynomial g_i
 * - The check g_i(0) + g_i(1) = (previous claim) constrains g_i
 * - The verifier's random challenge r_i "tests" g_i at a random point
 * - If the prover cheated (sent wrong g_i), it can only pass if
 *   the true polynomial and fake polynomial agree at r_i
 * - By Schwartz-Zippel lemma, this happens with probability <= d/|F|
 *   where d is the degree
 *
 * SOUNDNESS ERROR:
 * After n rounds, total error probability <= (d_1 + d_2 + ... + d_n) / |F|
 *
 * For multilinear polynomials: each d_i = 1
 * So error <= n / |F|
 *
 * For |F| >> n, this is negligible!
 */

console.log("Soundness Analysis:");
console.log("  If prover cheats at round i:");
console.log("    - Sends fake polynomial g'_i != g_i");
console.log("    - Must have g'_i(r_i) = g_i(r_i) to avoid detection");
console.log("    - Probability: <= deg(g_i) / |F| (Schwartz-Zippel)\n");

console.log("  Total soundness error:");
console.log("    <= (d_1 + d_2 + ... + d_n) / |F|");
console.log("    For multilinear f: <= n / |F|\n");

console.log(`  For our field F_${p} with n = 2:`);
console.log(`    Error <= ${2}/${p} ~ ${(Number(2) / Number(p)).toExponential(2)}\n`);

console.log("  For a 256-bit field with n = 1000:");
console.log("    Error <= 1000 / 2^256 ~ 10^{-74}\n");

// ============================================================================
// Section 7: Why SumCheck is Revolutionary
// ============================================================================

console.log("--- Section 7: Why SumCheck is Revolutionary ---\n");

/**
 * THE POWER OF SUMCHECK
 *
 * 1. REDUCES EXPONENTIAL TO LINEAR
 *    Verifier work: O(n) rounds instead of O(2^n) evaluations
 *
 * 2. COMPOSABILITY
 *    Can be combined with other protocols (polynomial commitments, IOP)
 *
 * 3. MULTILINEAR FRIENDLY
 *    Perfect for constraint systems expressed as multilinear polynomials
 *
 * 4. LINEAR-TIME PROVER (for structured polynomials)
 *    If f has special structure, prover can work in O(n * 2^n) time
 *    rather than O(n^2 * 2^n)
 *
 * 5. MINIMAL INTERACTION
 *    n rounds, each with one polynomial (few field elements) from prover
 *    and one random challenge from verifier
 *
 * APPLICATIONS IN ZK PROOFS:
 * - Spartan: constraint satisfaction via SumCheck
 * - HyperPlonk: gate constraints + wiring via SumCheck
 * - Lasso/Jolt: lookup arguments via SumCheck
 * - GKR: layered circuit computation via SumCheck
 */

console.log("Why SumCheck is Powerful:");
console.log("  1. Verification: O(2^n) -> O(n)");
console.log("  2. Composable with polynomial commitments");
console.log("  3. Perfect for multilinear constraints");
console.log("  4. Can achieve linear-time prover");
console.log("  5. Minimal rounds and communication\n");

console.log("Used in:");
console.log("  - HyperPlonk (this tutorial)");
console.log("  - Spartan");
console.log("  - Lasso / Jolt");
console.log("  - GKR protocol");
console.log("  - Many lookup arguments\n");

// ============================================================================
// Section 8: Preview: From SumCheck to ZeroCheck
// ============================================================================

console.log("--- Section 8: Preview of What's Next ---\n");

/**
 * ZEROCHECK: A SUMCHECK VARIANT
 *
 * ZeroCheck asks: Does f(b) = 0 for all b in B_n?
 *
 * This can be reduced to SumCheck!
 *
 * Key idea: f(b) = 0 for all b in B_n
 *           <=> sum_{b in B_n} eq(x, b) * f(b) = 0 for random x
 *
 * where eq(x, b) is the multilinear extension of the "equality" function.
 *
 * HyperPlonk uses ZeroCheck to verify:
 * - Gate constraints (all gates satisfied)
 * - Wiring constraints (all connections correct)
 *
 * Next: sumcheck-protocol.ts - The full protocol implementation
 */

console.log("Coming next:");
console.log("  - sumcheck-protocol.ts: Full prover/verifier implementation");
console.log("  - sumcheck-optimizations.ts: High-degree polynomial handling");
console.log("  - ZeroCheck: Proving f = 0 everywhere on B_n");
console.log("  - ProductCheck: Verifying product relations");
console.log("  - Full HyperPlonk: Putting it all together\n");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");

console.log("SumCheck Protocol:");
console.log("  Input: Multilinear polynomial f, claimed sum H");
console.log("  Claim: H = sum_{b in B_n} f(b)\n");

console.log("Protocol (n rounds):");
console.log("  1. Prover sends univariate g_i");
console.log("  2. Verifier checks g_i(0) + g_i(1) = previous claim");
console.log("  3. Verifier sends random r_i");
console.log("  4. New claim: g_i(r_i) = ... ");
console.log("  5. After n rounds: verify f(r_1, ..., r_n)\n");

console.log("Properties:");
console.log("  - Perfect completeness");
console.log("  - Soundness error <= n/|F|");
console.log("  - Verifier work: O(n) + one f evaluation");
console.log("  - Communication: O(n * d) field elements\n");

console.log("Next: Implementation of the full protocol");
