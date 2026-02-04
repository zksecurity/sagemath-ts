/**
 * # ZeroCheck Protocol: Full Implementation
 *
 * This file implements the complete ZeroCheck protocol, building on
 * the SumCheck implementation from the previous tutorial.
 *
 * ## Protocol Overview
 *
 * ZeroCheck proves that f(b) = 0 for all b in {0,1}^n by:
 * 1. Verifier sends random r
 * 2. Define h(x) = eq(x, r) * f(x)
 * 3. Run SumCheck to prove sum_{b in B_n} h(b) = 0
 * 4. Verify final evaluation h(s) = eq(s, r) * f(s)
 *
 * @module tutorial/part7-hyperplonk/20-zerocheck/zerocheck-protocol
 */

import { NotImplementedError } from 'sagemath-ts';

console.log("=== ZeroCheck Protocol Implementation ===\n");

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

// ============================================================================
// Type Definitions
// ============================================================================

/**
 * A function on the boolean hypercube.
 */
interface HypercubeFunction {
  numVars: number;
  /** values[i] = f(e) where e is the i-th point in lex order */
  values: bigint[];
}

/**
 * Univariate polynomial in evaluation form.
 */
interface EvaluationPolynomial {
  evaluations: bigint[];
  degree: number;
}

/**
 * Transcript of a ZeroCheck execution.
 */
interface ZeroCheckTranscript {
  /** Initial random challenge from verifier */
  initialChallenge: bigint[];
  /** SumCheck round polynomials */
  roundPolynomials: EvaluationPolynomial[];
  /** SumCheck challenges */
  sumcheckChallenges: bigint[];
  /** Final evaluation point from SumCheck */
  evaluationPoint: bigint[];
  /** Final claimed h value */
  finalHValue: bigint;
}

// ============================================================================
// Section 1: Helper Functions
// ============================================================================

console.log("--- Section 1: Building Blocks ---\n");

/**
 * Generates all points in {0,1}^n in lexicographic order.
 */
function booleanHypercube(n: number): number[][] {
  if (n === 0) return [[]];
  const smaller = booleanHypercube(n - 1);
  return [...smaller.map(v => [0, ...v]), ...smaller.map(v => [1, ...v])];
}

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
    const factor = add(mul(xi, ri), mul(sub(1n, xi), sub(1n, ri)));
    result = mul(result, factor);
  }

  return result;
}

/**
 * Evaluates multilinear extension at a point.
 */
function evaluateMLE(f: HypercubeFunction, x: bigint[]): bigint {
  const hypercube = booleanHypercube(f.numVars);
  let result = 0n;

  for (let i = 0; i < hypercube.length; i++) {
    const e = hypercube[i]!;
    const eField = e.map(v => BigInt(v));
    const eqVal = eq(x, eField);
    result = add(result, mul(f.values[i]!, eqVal));
  }

  return result;
}

/**
 * Lagrange interpolation from evaluations at 0, 1, ..., d.
 */
function evaluateFromEvaluations(
  evaluations: bigint[],
  r: bigint
): bigint {
  const d = evaluations.length - 1;
  let result = 0n;

  for (let i = 0; i <= d; i++) {
    let basis = 1n;
    for (let j = 0; j <= d; j++) {
      if (j !== i) {
        const num = sub(r, BigInt(j));
        const den = mod(BigInt(i - j));
        basis = mul(mul(basis, num), inv(den));
      }
    }
    result = add(result, mul(evaluations[i]!, basis));
  }

  return result;
}

console.log("Helper functions defined:");
console.log("  - booleanHypercube: generate {0,1}^n");
console.log("  - eq(x, r): equality polynomial");
console.log("  - evaluateMLE: multilinear extension evaluation");
console.log("  - evaluateFromEvaluations: Lagrange interpolation\n");

// ============================================================================
// Section 2: Computing h(x) = eq(x, r) * f(x)
// ============================================================================

console.log("--- Section 2: The Combined Polynomial h ---\n");

/**
 * COMPUTING h(x) = eq(x, r) * f(x)
 *
 * For ZeroCheck, we need to run SumCheck on h(x) = eq(x, r) * f(x).
 *
 * Properties of h:
 * - If f is multilinear, h has degree 2 in each variable
 * - h(b) = eq(b, r) * f(b) for boolean b
 * - If f(b) = 0 for all b, then h(b) = 0 for all b
 *
 * The prover computes h values on the hypercube:
 *   h(b) = eq(b, r) * f(b) for each b in B_n
 */

/**
 * Creates the h function h(x) = eq(x, r) * f(x).
 * Returns values on the boolean hypercube.
 */
function createH(f: HypercubeFunction, r: bigint[]): HypercubeFunction {
  const hypercube = booleanHypercube(f.numVars);
  const hValues: bigint[] = [];

  for (let i = 0; i < hypercube.length; i++) {
    const b = hypercube[i]!;
    const bField = b.map(v => BigInt(v));
    const eqVal = eq(bField, r);
    const hVal = mul(eqVal, f.values[i]!);
    hValues.push(hVal);
  }

  return {
    numVars: f.numVars,
    values: hValues
  };
}

/**
 * Evaluates h(x) = eq(x, r) * f(x) at any point x (not just boolean).
 */
function evaluateH(f: HypercubeFunction, r: bigint[], x: bigint[]): bigint {
  const fVal = evaluateMLE(f, x);
  const eqVal = eq(x, r);
  return mul(eqVal, fVal);
}

console.log("h(x) = eq(x, r) * f(x)");
console.log("  - Combines verifier's randomness with prover's polynomial");
console.log("  - h(b) = eq(b, r) * f(b) on boolean points");
console.log("  - If f = 0 everywhere, then h = 0 everywhere\n");

// ============================================================================
// Section 3: ZeroCheck to SumCheck Reduction
// ============================================================================

console.log("--- Section 3: Computing Round Polynomials ---\n");

/**
 * ROUND POLYNOMIAL FOR ZEROCHECK
 *
 * In round i, we compute:
 *   g_i(X) = sum_{remaining vars in {0,1}} h(r_1, ..., r_{i-1}, X, ...)
 *
 * Since h = eq * f and both are multilinear, h has degree 2.
 * So g_i has degree 2, and we need g_i(0), g_i(1), g_i(2).
 */

/**
 * Computes the round polynomial for ZeroCheck.
 * Returns evaluations at 0, 1, 2 (since h has degree 2).
 */
function computeZeroCheckRoundPoly(
  f: HypercubeFunction,
  r: bigint[],
  sumcheckChallenges: bigint[],
  round: number
): EvaluationPolynomial {
  const n = f.numVars;
  const remainingVars = n - round - 1;
  const remainingCube = booleanHypercube(remainingVars);

  // h(x) = eq(x, r) * f(x) has degree 2 in each variable
  const degree = 2;
  const evaluations: bigint[] = [];

  for (let j = 0; j <= degree; j++) {
    let sum = 0n;

    for (const suffix of remainingCube) {
      // Build evaluation point: (r_1, ..., r_{i-1}, j, suffix...)
      const x: bigint[] = [
        ...sumcheckChallenges,
        BigInt(j),
        ...suffix.map(v => BigInt(v))
      ];

      // Evaluate h(x) = eq(x, r) * f~(x)
      const hVal = evaluateH(f, r, x);
      sum = add(sum, hVal);
    }

    evaluations.push(sum);
  }

  return { evaluations, degree };
}

console.log("Round polynomial g_i(X) for ZeroCheck:");
console.log("  g_i(X) = sum_{remaining vars} h(challenges, X, ...boolean...)");
console.log("  Degree: 2 (since h = eq * f has degree 2)");
console.log("  Prover sends: g_i(0), g_i(1), g_i(2)\n");

// ============================================================================
// Section 4: ZeroCheck Prover
// ============================================================================

console.log("--- Section 4: ZeroCheck Prover ---\n");

/**
 * ZeroCheck prover algorithm.
 *
 * @param f - Function on the boolean hypercube (should be zero everywhere)
 * @param r - Initial verifier challenge
 * @param getChallenge - Function to get random challenges
 * @returns Complete ZeroCheck transcript
 */
function zeroCheckProver(
  f: HypercubeFunction,
  r: bigint[],
  getChallenge: () => bigint
): ZeroCheckTranscript {
  const n = f.numVars;
  const roundPolynomials: EvaluationPolynomial[] = [];
  const sumcheckChallenges: bigint[] = [];

  // Run SumCheck on h(x) = eq(x, r) * f(x) with claimed sum = 0
  for (let round = 0; round < n; round++) {
    const gi = computeZeroCheckRoundPoly(f, r, sumcheckChallenges, round);
    roundPolynomials.push(gi);

    const challenge = getChallenge();
    sumcheckChallenges.push(challenge);
  }

  // Final evaluation
  const evaluationPoint = sumcheckChallenges;
  const finalHValue = evaluateH(f, r, evaluationPoint);

  return {
    initialChallenge: r,
    roundPolynomials,
    sumcheckChallenges,
    evaluationPoint,
    finalHValue
  };
}

console.log("ZeroCheck Prover:");
console.log("  1. Receive random r from verifier");
console.log("  2. For each round, compute g_i for h(x) = eq(x, r) * f(x)");
console.log("  3. Receive challenge s_i, update partial point");
console.log("  4. Output final value h(s) where s = (s_1, ..., s_n)\n");

// ============================================================================
// Section 5: ZeroCheck Verifier
// ============================================================================

console.log("--- Section 5: ZeroCheck Verifier ---\n");

/**
 * ZeroCheck verifier algorithm.
 *
 * @param transcript - The prover's transcript
 * @param oracleF - Oracle access to f(x)
 * @returns true if verification passes
 */
function zeroCheckVerifier(
  transcript: ZeroCheckTranscript,
  oracleF: (x: bigint[]) => bigint
): boolean {
  const n = transcript.roundPolynomials.length;
  const r = transcript.initialChallenge;

  console.log(`  Verifying ZeroCheck with ${n} rounds...`);
  console.log(`  Initial challenge r = (${r.join(', ')})\n`);

  // The claimed sum for ZeroCheck is 0
  let currentClaim = 0n;

  for (let round = 0; round < n; round++) {
    const gi = transcript.roundPolynomials[round]!;
    const si = transcript.sumcheckChallenges[round]!;

    // Check: g_i(0) + g_i(1) = current claim
    const g0 = gi.evaluations[0]!;
    const g1 = gi.evaluations[1]!;
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

    // Update claim: g_i(s_i) via Lagrange interpolation
    currentClaim = evaluateFromEvaluations(gi.evaluations, si);
    console.log(`    s_${round + 1} = ${si}, new claim = ${currentClaim}\n`);
  }

  // Final check: h(s) = eq(s, r) * f(s)
  const s = transcript.evaluationPoint;
  const eqSR = eq(s, r);
  const fAtS = oracleF(s);
  const expectedH = mul(eqSR, fAtS);

  console.log(`  Final check:`);
  console.log(`    s = (${s.join(', ')})`);
  console.log(`    eq(s, r) = ${eqSR}`);
  console.log(`    f(s) = ${fAtS}`);
  console.log(`    Expected h(s) = eq(s, r) * f(s) = ${expectedH}`);
  console.log(`    Claimed h(s) = ${currentClaim}`);

  if (expectedH !== currentClaim) {
    console.log(`    FAILED: final evaluation\n`);
    return false;
  }

  console.log(`    PASSED\n`);
  return true;
}

console.log("ZeroCheck Verifier:");
console.log("  1. Check each round: g_i(0) + g_i(1) = previous claim");
console.log("  2. Compute g_i(s_i) via Lagrange interpolation");
console.log("  3. Final check: h(s) = eq(s, r) * f(s)\n");

// ============================================================================
// Section 6: Complete Example - Zero Polynomial
// ============================================================================

console.log("--- Section 6: Example with Zero Polynomial ---\n");

// Create a polynomial that IS zero on {0,1}^2
// On {0,1}^2: f(0,0)=0, f(0,1)=0, f(1,0)=0, f(1,1)=0
const zeroF: HypercubeFunction = {
  numVars: 2,
  values: [0n, 0n, 0n, 0n]
};

console.log("Example: f that is zero on {0,1}^2");
console.log("  f(0,0) = 0, f(0,1) = 0, f(1,0) = 0, f(1,1) = 0\n");

// Verifier's initial challenge
const initR = [3n, 7n];
console.log(`Verifier's challenge: r = (${initR.join(', ')})\n`);

// Deterministic challenges for reproducibility
let challengeIdx = 0;
const challenges = [5n, 9n];
function getChallenge(): bigint {
  return challenges[challengeIdx++]!;
}

// Run prover
console.log("=== PROVER ===\n");
const zeroTranscript = zeroCheckProver(zeroF, initR, getChallenge);

console.log("Transcript:");
for (let i = 0; i < zeroTranscript.roundPolynomials.length; i++) {
  const poly = zeroTranscript.roundPolynomials[i]!;
  console.log(`  Round ${i + 1}: g(0)=${poly.evaluations[0]!}, g(1)=${poly.evaluations[1]!}, g(2)=${poly.evaluations[2]!}`);
  console.log(`  Challenge: s_${i + 1} = ${zeroTranscript.sumcheckChallenges[i]!}`);
}
console.log(`  Final evaluation point: (${zeroTranscript.evaluationPoint.join(', ')})`);
console.log(`  Final h value: ${zeroTranscript.finalHValue}\n`);

// Run verifier
console.log("=== VERIFIER ===\n");
const oracle = (x: bigint[]) => evaluateMLE(zeroF, x);
const accepted = zeroCheckVerifier(zeroTranscript, oracle);
console.log(`Verification result: ${accepted ? 'ACCEPTED' : 'REJECTED'}\n`);

// ============================================================================
// Section 7: Example - Non-Zero Polynomial (Cheating Detection)
// ============================================================================

console.log("--- Section 7: Detecting Non-Zero Polynomial ---\n");

// Create a polynomial that is NOT zero everywhere
// f(0,0)=1, f(0,1)=2, f(1,0)=3, f(1,1)=4
const nonZeroF: HypercubeFunction = {
  numVars: 2,
  values: [1n, 2n, 3n, 4n]
};

console.log("Cheating example: f is NOT zero on {0,1}^2");
console.log("  f(0,0) = 1, f(0,1) = 2, f(1,0) = 3, f(1,1) = 4\n");

// Reset challenge counter
challengeIdx = 0;

// Run "cheating" prover (trying to prove f = 0 when it's not)
console.log("=== CHEATING PROVER ===\n");
const cheatTranscript = zeroCheckProver(nonZeroF, initR, getChallenge);

console.log("Transcript (cheating):");
for (let i = 0; i < cheatTranscript.roundPolynomials.length; i++) {
  const poly = cheatTranscript.roundPolynomials[i]!;
  console.log(`  Round ${i + 1}: g(0)=${poly.evaluations[0]!}, g(1)=${poly.evaluations[1]!}, g(2)=${poly.evaluations[2]!}`);
}
console.log(`  Final h value: ${cheatTranscript.finalHValue}\n`);

// Run verifier
console.log("=== VERIFIER ===\n");
const cheatOracle = (x: bigint[]) => evaluateMLE(nonZeroF, x);
const cheatAccepted = zeroCheckVerifier(cheatTranscript, cheatOracle);
console.log(`Verification result: ${cheatAccepted ? 'ACCEPTED (BAD!)' : 'REJECTED (GOOD!)'}\n`);

// ============================================================================
// Section 8: Production Stubs
// ============================================================================

console.log("--- Section 8: Production Implementation Stubs ---\n");

/**
 * Optimized ZeroCheck prover with O(n * 2^n) total work.
 *
 * Uses bookkeeping to avoid recomputing sums from scratch each round.
 * The key insight is that h(x) = eq(x, r) * f(x) can be incrementally
 * updated as we fix variables.
 *
 * @throws NotImplementedError - Optimization not yet implemented
 */
export function zeroCheckProverOptimized(
  f: HypercubeFunction,
  r: bigint[],
  getChallenge: () => bigint
): ZeroCheckTranscript {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: zeroCheckProverOptimized - ' +
    'Requires bookkeeping for O(n * 2^n) prover time'
  );
}

/**
 * Batch ZeroCheck for multiple polynomials.
 *
 * Proves: f_1 = 0 AND f_2 = 0 AND ... AND f_k = 0 on B_n.
 *
 * Uses random combination: g = alpha_1 * f_1 + ... + alpha_k * f_k
 * and proves g = 0 with a single ZeroCheck.
 *
 * @throws NotImplementedError - Batch optimization not yet implemented
 */
export function batchZeroCheck(
  functions: HypercubeFunction[],
  combinationCoeffs: bigint[],
  r: bigint[],
  getChallenge: () => bigint
): ZeroCheckTranscript {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: batchZeroCheck - ' +
    'Batch ZeroCheck for multiple polynomial instances'
  );
}

/**
 * ZeroCheck with polynomial commitment opening.
 *
 * Instead of oracle access to f, the prover provides a commitment
 * and opening proof for f(s).
 *
 * @throws NotImplementedError - Commitment integration not yet implemented
 */
export function zeroCheckWithCommitment(
  fCommitment: unknown,
  r: bigint[],
  getChallenge: () => bigint
): ZeroCheckTranscript {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: zeroCheckWithCommitment - ' +
    'ZeroCheck with polynomial commitment scheme'
  );
}

console.log("Production stubs:");
console.log("  - zeroCheckProverOptimized: O(n * 2^n) prover");
console.log("  - batchZeroCheck: prove multiple f_i = 0");
console.log("  - zeroCheckWithCommitment: integrated with commitments\n");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");

console.log("ZeroCheck Protocol Implementation:");
console.log("  - createH: computes h(x) = eq(x, r) * f(x) on hypercube");
console.log("  - computeZeroCheckRoundPoly: prover's round polynomial");
console.log("  - zeroCheckProver: complete prover algorithm");
console.log("  - zeroCheckVerifier: complete verifier algorithm\n");

console.log("Protocol Properties:");
console.log("  - Proves f(b) = 0 for all b in {0,1}^n");
console.log("  - Reduces to SumCheck with claimed sum = 0");
console.log("  - Communication: O(3n) field elements (degree 2)");
console.log("  - Verifier: O(n) + one f evaluation\n");

console.log("Next: ProductCheck - verifying product relations");

// Export for other tutorials
export {
  eq,
  evaluateMLE,
  createH,
  evaluateH,
  zeroCheckProver,
  zeroCheckVerifier,
  type HypercubeFunction,
  type ZeroCheckTranscript
};
