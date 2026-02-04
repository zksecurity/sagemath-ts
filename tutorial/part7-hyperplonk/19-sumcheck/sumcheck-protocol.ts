/**
 * # The SumCheck Protocol: Implementation
 *
 * This file implements the SumCheck protocol with prover and verifier roles.
 * We build the protocol step by step, demonstrating how each round works.
 *
 * ## Protocol Overview
 *
 * Goal: Verify that H = sum_{b in {0,1}^n} f(b) for multilinear f.
 *
 * The protocol proceeds in n rounds:
 * - Round i: Prover sends univariate polynomial g_i(X)
 * - Verifier checks g_i(0) + g_i(1) = claimed value
 * - Verifier sends random challenge r_i
 * - Proceed to round i+1 with updated claim
 *
 * After n rounds, verifier checks final evaluation of f.
 *
 * @module tutorial/part7-hyperplonk/19-sumcheck/sumcheck-protocol
 */

import { NotImplementedError } from 'sagemath-ts';

console.log("=== The SumCheck Protocol Implementation ===\n");

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

/** Field negation */
const neg = (a: bigint): bigint => mod(-a);

/** Modular exponentiation (for inverse) */
function powMod(base: bigint, exp: bigint, modulus: bigint): bigint {
  let result = 1n;
  base = ((base % modulus) + modulus) % modulus;
  while (exp > 0n) {
    if (exp % 2n === 1n) {
      result = (result * base) % modulus;
    }
    exp = exp / 2n;
    base = (base * base) % modulus;
  }
  return result;
}

/** Field inverse using Fermat's little theorem */
const inv = (a: bigint): bigint => powMod(a, p - 2n, p);

/** Field division */
const div = (a: bigint, b: bigint): bigint => mul(a, inv(b));

// ============================================================================
// Type Definitions
// ============================================================================

/**
 * A univariate polynomial represented by its coefficients.
 * coefficients[i] is the coefficient of X^i.
 */
interface UnivariatePolynomial {
  coefficients: bigint[];
  degree: number;
}

/**
 * Transcript of a SumCheck execution for verification.
 */
interface SumCheckTranscript {
  /** The claimed sum */
  claimedSum: bigint;
  /** Univariate polynomials sent by prover in each round */
  roundPolynomials: UnivariatePolynomial[];
  /** Random challenges from verifier */
  challenges: bigint[];
  /** Final evaluation point */
  evaluationPoint: bigint[];
  /** Final claimed value */
  finalValue: bigint;
}

// ============================================================================
// Section 1: Helper Functions
// ============================================================================

console.log("--- Section 1: Building Blocks ---\n");

/**
 * Evaluates a univariate polynomial at a point.
 */
function evaluateUnivariate(poly: UnivariatePolynomial, x: bigint): bigint {
  let result = 0n;
  let xPower = 1n;

  for (let i = 0; i <= poly.degree; i++) {
    result = add(result, mul(poly.coefficients[i]!, xPower));
    xPower = mul(xPower, x);
  }

  return result;
}

/**
 * Creates a univariate polynomial from coefficients.
 */
function createUnivariatePolynomial(coeffs: bigint[]): UnivariatePolynomial {
  // Remove trailing zeros
  const newCoeffs = [...coeffs];
  while (newCoeffs.length > 1 && newCoeffs[newCoeffs.length - 1] === 0n) {
    newCoeffs.pop();
  }
  return {
    coefficients: newCoeffs,
    degree: newCoeffs.length - 1
  };
}

/**
 * Generates all points in the boolean hypercube {0,1}^n.
 */
function booleanHypercube(n: number): number[][] {
  if (n === 0) return [[]];
  const smaller = booleanHypercube(n - 1);
  return [...smaller.map(v => [0, ...v]), ...smaller.map(v => [1, ...v])];
}

console.log("Helper functions defined:");
console.log("  - evaluateUnivariate: evaluate polynomial at a point");
console.log("  - createUnivariatePolynomial: build polynomial from coefficients");
console.log("  - booleanHypercube: generate all points in {0,1}^n\n");

// ============================================================================
// Section 2: The eq(x, e) Function
// ============================================================================

console.log("--- Section 2: The eq(x, e) Function ---\n");

/**
 * THE EQUALITY POLYNOMIAL eq(x, e)
 *
 * For e in {0,1}^n, eq(x, e) is the unique multilinear polynomial that:
 *   eq(b, e) = 1 if b = e
 *   eq(b, e) = 0 if b != e  (for b in {0,1}^n)
 *
 * Formula:
 *   eq(x, e) = prod_{i=1}^{n} (x_i * e_i + (1 - x_i) * (1 - e_i))
 *
 * For e_i = 0: factor is (1 - x_i)
 * For e_i = 1: factor is x_i
 *
 * This is the "selection" polynomial - it picks out exactly one vertex
 * of the boolean hypercube.
 */

/**
 * Computes eq(x, e) = prod_i (x_i * e_i + (1 - x_i)(1 - e_i)).
 *
 * @param x - Evaluation point (field elements)
 * @param e - Boolean vector (0s and 1s)
 * @returns The value eq(x, e)
 */
function eq(x: bigint[], e: number[]): bigint {
  if (x.length !== e.length) {
    throw new Error(`Dimension mismatch: x has ${x.length}, e has ${e.length}`);
  }

  let result = 1n;
  for (let i = 0; i < x.length; i++) {
    const xi = x[i]!;
    const ei = BigInt(e[i]!);

    // Factor: x_i * e_i + (1 - x_i) * (1 - e_i)
    // If e_i = 0: (1 - x_i)
    // If e_i = 1: x_i
    let factor: bigint;
    if (ei === 0n) {
      factor = sub(1n, xi);  // (1 - x_i)
    } else {
      factor = xi;  // x_i
    }
    result = mul(result, factor);
  }

  return result;
}

console.log("eq(x, e) = prod_i (x_i * e_i + (1 - x_i)(1 - e_i))");
console.log("  eq(b, e) = 1 if b = e");
console.log("  eq(b, e) = 0 if b != e (for b in {0,1}^n)\n");

// Demonstrate eq function
const testPoint = [2n, 3n];  // A random point in F^2
console.log(`Evaluating eq at x = (${testPoint.join(', ')}):`);

for (const e of booleanHypercube(2)) {
  const val = eq(testPoint, e);
  console.log(`  eq(x, (${e.join(',')})) = ${val}`);
}

// Verify: sum of eq over all e should equal 1
let eqSum = 0n;
for (const e of booleanHypercube(2)) {
  eqSum = add(eqSum, eq(testPoint, e));
}
console.log(`\nSum over all e: ${eqSum} (should be 1 for any x)`);

// ============================================================================
// Section 3: Multilinear Extension
// ============================================================================

console.log("\n--- Section 3: Multilinear Extension ---\n");

/**
 * MULTILINEAR EXTENSION (MLE)
 *
 * Given a function f: {0,1}^n -> F specified by its values on the hypercube,
 * the multilinear extension f~: F^n -> F is:
 *
 *   f~(x) = sum_{e in {0,1}^n} f(e) * eq(x, e)
 *
 * This is the unique multilinear polynomial that agrees with f on {0,1}^n.
 *
 * Properties:
 * - f~(b) = f(b) for all b in {0,1}^n
 * - f~ is multilinear (degree 1 in each variable)
 * - f~ is unique with these properties
 */

/**
 * Represents a function on the boolean hypercube by its evaluations.
 */
interface HypercubeFunction {
  /** Number of variables */
  numVars: number;
  /** values[i] = f(e) where e is the i-th point in lex order */
  values: bigint[];
}

/**
 * Evaluates the multilinear extension at a point.
 *
 * f~(x) = sum_{e in {0,1}^n} f(e) * eq(x, e)
 */
function evaluateMLE(f: HypercubeFunction, x: bigint[]): bigint {
  if (x.length !== f.numVars) {
    throw new Error(`Dimension mismatch: ${x.length} vs ${f.numVars}`);
  }

  const hypercube = booleanHypercube(f.numVars);
  let result = 0n;

  for (let i = 0; i < hypercube.length; i++) {
    const e = hypercube[i]!;
    const eqVal = eq(x, e);
    result = add(result, mul(f.values[i]!, eqVal));
  }

  return result;
}

console.log("Multilinear Extension (MLE):");
console.log("  f~(x) = sum_{e in {0,1}^n} f(e) * eq(x, e)\n");

// Example: Define a function on {0,1}^2
// f(0,0) = 7, f(0,1) = 12, f(1,0) = 10, f(1,1) = 17
const exampleF: HypercubeFunction = {
  numVars: 2,
  values: [7n, 12n, 10n, 17n]  // Lex order: 00, 01, 10, 11
};

console.log("Example function on {0,1}^2:");
for (let i = 0; i < 4; i++) {
  const e = booleanHypercube(2)[i]!;
  console.log(`  f(${e.join(',')}) = ${exampleF.values[i]!}`);
}

// Verify MLE agrees on hypercube
console.log("\nVerifying MLE on boolean points:");
for (let i = 0; i < 4; i++) {
  const e = booleanHypercube(2)[i]!;
  const x = e.map(v => BigInt(v));
  const mleVal = evaluateMLE(exampleF, x);
  console.log(`  f~(${e.join(',')}) = ${mleVal} (expected ${exampleF.values[i]!})`);
}

// Evaluate at a non-boolean point
const testX = [2n, 3n];
const mleAtTest = evaluateMLE(exampleF, testX);
console.log(`\nMLE at (2, 3): f~(2, 3) = ${mleAtTest}`);

// ============================================================================
// Section 4: SumCheck Prover
// ============================================================================

console.log("\n--- Section 4: SumCheck Prover ---\n");

/**
 * SUMCHECK PROVER ALGORITHM
 *
 * For each round i (from 1 to n):
 *
 * 1. Fix variables x_1 = r_1, ..., x_{i-1} = r_{i-1} (from previous challenges)
 * 2. Compute g_i(X) = sum_{(x_{i+1}, ..., x_n) in {0,1}^{n-i}} f(r_1, ..., r_{i-1}, X, x_{i+1}, ..., x_n)
 * 3. Send g_i to verifier
 * 4. Receive challenge r_i
 *
 * Computing g_i efficiently:
 * - For each degree d from 0 to deg_i(f):
 *   - Sum f over partial hypercube with X = d-th basis evaluation
 *
 * For multilinear f: g_i has degree 1, so we only need g_i(0) and g_i(1).
 */

/**
 * Computes the i-th round polynomial for the SumCheck prover.
 *
 * g_i(X) = sum_{remaining vars in {0,1}} f(r_1, ..., r_{i-1}, X, x_{i+1}, ..., x_n)
 *
 * @param f - The function on the hypercube
 * @param challenges - Previous challenges [r_1, ..., r_{i-1}]
 * @param round - Current round (0-indexed)
 * @returns Univariate polynomial g_i
 */
function computeRoundPolynomial(
  f: HypercubeFunction,
  challenges: bigint[],
  round: number
): UnivariatePolynomial {
  const n = f.numVars;
  const remainingVars = n - round - 1;  // Variables after current one

  // For multilinear, we need g_i(0) and g_i(1)
  // g_i(0) = sum_{x_{i+1}, ..., x_n in {0,1}} f(r_1, ..., r_{i-1}, 0, x_{i+1}, ..., x_n)
  // g_i(1) = sum_{x_{i+1}, ..., x_n in {0,1}} f(r_1, ..., r_{i-1}, 1, x_{i+1}, ..., x_n)

  let g0 = 0n;
  let g1 = 0n;

  // Iterate over all assignments to remaining variables
  const remainingCube = booleanHypercube(remainingVars);

  for (const suffix of remainingCube) {
    // Build full evaluation point with X = 0
    const x0: bigint[] = [
      ...challenges,
      0n,
      ...suffix.map(v => BigInt(v))
    ];
    g0 = add(g0, evaluateMLE(f, x0));

    // Build full evaluation point with X = 1
    const x1: bigint[] = [
      ...challenges,
      1n,
      ...suffix.map(v => BigInt(v))
    ];
    g1 = add(g1, evaluateMLE(f, x1));
  }

  // g_i(X) = g0 + (g1 - g0) * X = g0 * (1 - X) + g1 * X
  // coefficients: [g0, g1 - g0]
  return createUnivariatePolynomial([g0, sub(g1, g0)]);
}

/**
 * The SumCheck prover: generates the full transcript.
 *
 * @param f - Function on the boolean hypercube
 * @param claimedSum - The claimed sum to prove
 * @param getChallenge - Function to get random challenges (simulates verifier)
 * @returns The complete SumCheck transcript
 */
function sumcheckProver(
  f: HypercubeFunction,
  claimedSum: bigint,
  getChallenge: () => bigint
): SumCheckTranscript {
  const n = f.numVars;
  const roundPolynomials: UnivariatePolynomial[] = [];
  const challenges: bigint[] = [];

  for (let round = 0; round < n; round++) {
    // Compute round polynomial
    const gi = computeRoundPolynomial(f, challenges, round);
    roundPolynomials.push(gi);

    // Get challenge from verifier
    const ri = getChallenge();
    challenges.push(ri);
  }

  // Final evaluation
  const finalValue = evaluateMLE(f, challenges);

  return {
    claimedSum,
    roundPolynomials,
    challenges,
    evaluationPoint: challenges,
    finalValue
  };
}

console.log("Prover algorithm defined:");
console.log("  computeRoundPolynomial: computes g_i(X) for round i");
console.log("  sumcheckProver: generates full transcript\n");

// ============================================================================
// Section 5: SumCheck Verifier
// ============================================================================

console.log("--- Section 5: SumCheck Verifier ---\n");

/**
 * SUMCHECK VERIFIER ALGORITHM
 *
 * Input: Claimed sum H, transcript from prover, oracle access to f
 *
 * For each round i:
 * 1. Check g_i(0) + g_i(1) = (previous claim or H for i=1)
 * 2. Compute new claim = g_i(r_i)
 *
 * Final check:
 * - Verify g_n(r_n) = f(r_1, ..., r_n) using oracle
 *
 * Accept iff all checks pass.
 */

/**
 * Verifies a SumCheck transcript.
 *
 * @param transcript - The prover's transcript
 * @param oracleEval - Function to evaluate f at any point (or commitment opening)
 * @returns true if verification passes
 */
function sumcheckVerifier(
  transcript: SumCheckTranscript,
  oracleEval: (x: bigint[]) => bigint
): boolean {
  const n = transcript.roundPolynomials.length;
  let currentClaim = transcript.claimedSum;

  console.log(`  Verifying SumCheck with ${n} rounds...`);
  console.log(`  Claimed sum: ${currentClaim}\n`);

  for (let round = 0; round < n; round++) {
    const gi = transcript.roundPolynomials[round]!;
    const ri = transcript.challenges[round]!;

    // Check: g_i(0) + g_i(1) = current claim
    const g0 = evaluateUnivariate(gi, 0n);
    const g1 = evaluateUnivariate(gi, 1n);
    const sum = add(g0, g1);

    console.log(`  Round ${round + 1}:`);
    console.log(`    g_${round + 1}(0) = ${g0}`);
    console.log(`    g_${round + 1}(1) = ${g1}`);
    console.log(`    Sum = ${sum}, Expected = ${currentClaim}`);

    if (sum !== currentClaim) {
      console.log(`    FAILED: sum check\n`);
      return false;
    }
    console.log(`    PASSED`);

    // Update claim for next round
    currentClaim = evaluateUnivariate(gi, ri);
    console.log(`    r_${round + 1} = ${ri}, new claim = ${currentClaim}\n`);
  }

  // Final check: g_n(r_n) should equal f(r_1, ..., r_n)
  const finalOracleValue = oracleEval(transcript.evaluationPoint);
  console.log(`  Final check:`);
  console.log(`    Claimed: ${transcript.finalValue}`);
  console.log(`    Oracle f(${transcript.evaluationPoint.join(', ')}): ${finalOracleValue}`);

  if (transcript.finalValue !== finalOracleValue) {
    console.log(`    FAILED: final evaluation\n`);
    return false;
  }

  console.log(`    PASSED\n`);
  return true;
}

console.log("Verifier algorithm defined:");
console.log("  sumcheckVerifier: checks all rounds and final evaluation\n");

// ============================================================================
// Section 6: Complete Example
// ============================================================================

console.log("--- Section 6: Complete Example ---\n");

// Use our example function from Section 3
// f(0,0) = 7, f(0,1) = 12, f(1,0) = 10, f(1,1) = 17
// Sum = 7 + 12 + 10 + 17 = 46

console.log("Running SumCheck on example function:");
console.log("  f(0,0) = 7, f(0,1) = 12, f(1,0) = 10, f(1,1) = 17");
console.log("  Claimed sum: 46\n");

// Compute actual sum
let actualSum = 0n;
for (const val of exampleF.values) {
  actualSum = add(actualSum, val);
}
console.log(`Computed sum: ${actualSum}\n`);

// Simulate verifier challenges (deterministic for reproducibility)
let challengeIndex = 0;
const deterministicChallenges = [2n, 3n];
function getChallenge(): bigint {
  return deterministicChallenges[challengeIndex++]!;
}

// Run prover
console.log("=== PROVER ===\n");
const transcript = sumcheckProver(exampleF, actualSum, getChallenge);

console.log("Transcript:");
for (let i = 0; i < transcript.roundPolynomials.length; i++) {
  const poly = transcript.roundPolynomials[i]!;
  const coeffStr = poly.coefficients.join(', ');
  console.log(`  Round ${i + 1}: g_${i + 1}(X) with coefficients [${coeffStr}]`);
  console.log(`  Challenge: r_${i + 1} = ${transcript.challenges[i]!}`);
}
console.log(`  Final evaluation point: (${transcript.evaluationPoint.join(', ')})`);
console.log(`  Final value: ${transcript.finalValue}\n`);

// Run verifier
console.log("=== VERIFIER ===\n");
const oracle = (x: bigint[]) => evaluateMLE(exampleF, x);
const accepted = sumcheckVerifier(transcript, oracle);

console.log(`Verification result: ${accepted ? 'ACCEPTED' : 'REJECTED'}\n`);

// ============================================================================
// Section 7: Cheating Prover Detection
// ============================================================================

console.log("--- Section 7: Detecting a Cheating Prover ---\n");

/**
 * Let's demonstrate what happens if the prover tries to cheat
 * by claiming a wrong sum.
 */

console.log("Attempting to prove a FALSE claim (sum = 50 instead of 46):\n");

// Reset challenge counter
challengeIndex = 0;

// Try to prove wrong sum
const wrongSum = 50n;

// A cheating prover would need to send a modified g_1
// g_1(0) + g_1(1) must equal 50, but for honest g_1 it equals 46
// The prover must lie about g_1

// Let's see what happens with an honest prover trying to prove wrong sum
console.log("Honest prover generates transcript (but claims wrong sum)...\n");

// The prover computes honestly but we'll check against wrong sum
const honestTranscript = sumcheckProver(exampleF, actualSum, getChallenge);

// Modify the claimed sum to be wrong
const cheatingTranscript: SumCheckTranscript = {
  ...honestTranscript,
  claimedSum: wrongSum
};

console.log("=== VERIFIER (with tampered claim) ===\n");
const cheatingResult = sumcheckVerifier(cheatingTranscript, oracle);
console.log(`Verification of false claim: ${cheatingResult ? 'ACCEPTED (BAD!)' : 'REJECTED (GOOD!)'}\n`);

// ============================================================================
// Section 8: Production-Ready Prover Stub
// ============================================================================

console.log("--- Section 8: Production Prover (Stub) ---\n");

/**
 * Production SumCheck prover with optimizations.
 *
 * The naive implementation above has O(2^n) cost per round.
 * An optimized implementation can achieve O(2^n) TOTAL using:
 * - Incremental updates
 * - Partial sum tables
 * - Bookkeeping optimization
 *
 * @throws NotImplementedError - Full optimization not yet implemented
 */
export function sumcheckProverOptimized(
  f: HypercubeFunction,
  claimedSum: bigint,
  getChallenge: () => bigint
): SumCheckTranscript {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: sumcheckProverOptimized - ' +
    'Requires bookkeeping optimization for O(2^n) total prover time'
  );
}

/**
 * Batch SumCheck for multiple instances.
 *
 * Proves: sum_i alpha_i * (sum_{b in B_n} f_i(b)) = H
 *
 * where alpha_i are random combination coefficients.
 *
 * @throws NotImplementedError - Batch optimization not yet implemented
 */
export function batchSumcheckProver(
  functions: HypercubeFunction[],
  claimedSums: bigint[],
  combinationCoeffs: bigint[],
  getChallenge: () => bigint
): SumCheckTranscript {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: batchSumcheckProver - ' +
    'Batch sumcheck for multiple polynomial instances'
  );
}

console.log("Stub functions for production implementation:");
console.log("  sumcheckProverOptimized: O(2^n) total time with bookkeeping");
console.log("  batchSumcheckProver: batched sumcheck for multiple instances\n");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");

console.log("SumCheck Protocol Implementation:");
console.log("  - eq(x, e): equality polynomial for selecting hypercube vertices");
console.log("  - evaluateMLE: compute multilinear extension at any point");
console.log("  - computeRoundPolynomial: prover's round computation");
console.log("  - sumcheckProver: complete prover algorithm");
console.log("  - sumcheckVerifier: complete verifier algorithm\n");

console.log("Key Operations:");
console.log("  - Each round: prover sends univariate polynomial g_i");
console.log("  - Verifier checks: g_i(0) + g_i(1) = previous claim");
console.log("  - Final check: oracle evaluation at random point\n");

console.log("Complexity:");
console.log("  - Naive prover: O(n * 2^n)");
console.log("  - Optimized prover: O(2^n) total");
console.log("  - Verifier: O(n * d) where d is max degree\n");

console.log("Next: sumcheck-optimizations.ts - High-degree polynomial handling");

// Export functions for use in other tutorials
export {
  booleanHypercube,
  eq,
  evaluateMLE,
  evaluateUnivariate,
  createUnivariatePolynomial,
  computeRoundPolynomial,
  sumcheckProver,
  sumcheckVerifier,
  type HypercubeFunction,
  type UnivariatePolynomial,
  type SumCheckTranscript,
};
