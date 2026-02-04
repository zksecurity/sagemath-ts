/**
 * # ProductCheck Protocol: Full Implementation
 *
 * This file implements the ProductCheck protocol that verifies:
 *   prod_{b in B_n} f(b) = prod_{b in B_n} g(b)
 *
 * The protocol works by:
 * 1. Constructing an auxiliary polynomial h tracking the running product
 * 2. Reducing to ZeroCheck on the constraint h(next(b))*g(b) = h(b)*f(b)
 *
 * @module tutorial/part7-hyperplonk/21-productcheck/productcheck-protocol
 */

import { NotImplementedError } from 'sagemath-ts';

console.log("=== ProductCheck Protocol Implementation ===\n");

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

interface HypercubeFunction {
  numVars: number;
  values: bigint[];
}

interface ProductCheckTranscript {
  /** The auxiliary polynomial h values on the hypercube */
  hValues: bigint[];
  /** ZeroCheck initial challenge */
  zeroCheckChallenge: bigint[];
  /** ZeroCheck round polynomials */
  roundPolynomials: EvaluationPolynomial[];
  /** SumCheck challenges */
  sumcheckChallenges: bigint[];
  /** Final evaluation point */
  evaluationPoint: bigint[];
  /** Final claimed values */
  finalValues: {
    h: bigint;
    hShifted: bigint;
    f: bigint;
    g: bigint;
  };
}

interface EvaluationPolynomial {
  evaluations: bigint[];
  degree: number;
}

// ============================================================================
// Section 1: Helper Functions
// ============================================================================

console.log("--- Section 1: Building Blocks ---\n");

function booleanHypercube(n: number): number[][] {
  if (n === 0) return [[]];
  const smaller = booleanHypercube(n - 1);
  return [...smaller.map(v => [0, ...v]), ...smaller.map(v => [1, ...v])];
}

function eq(x: bigint[], r: bigint[]): bigint {
  let result = 1n;
  for (let i = 0; i < x.length; i++) {
    const factor = add(mul(x[i]!, r[i]!), mul(sub(1n, x[i]!), sub(1n, r[i]!)));
    result = mul(result, factor);
  }
  return result;
}

function evaluateMLE(f: HypercubeFunction, x: bigint[]): bigint {
  const hypercube = booleanHypercube(f.numVars);
  let result = 0n;
  for (let i = 0; i < hypercube.length; i++) {
    const e = hypercube[i]!.map(v => BigInt(v));
    result = add(result, mul(f.values[i]!, eq(x, e)));
  }
  return result;
}

function evaluateFromEvaluations(evaluations: bigint[], r: bigint): bigint {
  const d = evaluations.length - 1;
  let result = 0n;
  for (let i = 0; i <= d; i++) {
    let basis = 1n;
    for (let j = 0; j <= d; j++) {
      if (j !== i) {
        basis = mul(mul(basis, sub(r, BigInt(j))), inv(mod(BigInt(i - j))));
      }
    }
    result = add(result, mul(evaluations[i]!, basis));
  }
  return result;
}

console.log("Helper functions loaded:");
console.log("  - booleanHypercube, eq, evaluateMLE, evaluateFromEvaluations\n");

// ============================================================================
// Section 2: Computing the Auxiliary Polynomial h
// ============================================================================

console.log("--- Section 2: Computing Auxiliary Polynomial h ---\n");

/**
 * COMPUTING THE RUNNING PRODUCT h
 *
 * h is defined such that:
 *   h(first) = 1
 *   h(next(b)) = h(b) * f(b) / g(b)
 *
 * After all points: h(wrap-around) = prod_b (f(b)/g(b))
 *
 * For ProductCheck to pass, this final value must equal 1.
 */

/**
 * Computes the running product polynomial h.
 *
 * @param f - Numerator polynomial on hypercube
 * @param g - Denominator polynomial on hypercube
 * @returns h values where h[i] = prod_{j<i} (f[j]/g[j])
 */
function computeRunningProduct(
  f: HypercubeFunction,
  g: HypercubeFunction
): HypercubeFunction {
  const n = f.numVars;
  const numPoints = 2 ** n;
  const hValues: bigint[] = [];

  // h(first) = 1
  let runningProd = 1n;
  hValues.push(runningProd);

  // h(next(b)) = h(b) * f(b) / g(b)
  for (let i = 0; i < numPoints - 1; i++) {
    const ratio = mul(f.values[i]!, inv(g.values[i]!));
    runningProd = mul(runningProd, ratio);
    hValues.push(runningProd);
  }

  // The last computed value is prod_{b < last} (f(b)/g(b))
  // The wrap-around product (including last) determines correctness

  return { numVars: n, values: hValues };
}

/**
 * Computes the shifted version of h: h_shifted(b) = h(next(b)).
 */
function computeShiftedH(h: HypercubeFunction): HypercubeFunction {
  const n = h.numVars;
  const numPoints = 2 ** n;
  const shiftedValues: bigint[] = [];

  for (let i = 0; i < numPoints; i++) {
    const nextIdx = (i + 1) % numPoints;
    shiftedValues.push(h.values[nextIdx]!);
  }

  return { numVars: n, values: shiftedValues };
}

console.log("Running product h:");
console.log("  h(first) = 1");
console.log("  h(next(b)) = h(b) * f(b) / g(b)");
console.log("  => h encodes the cumulative product\n");

// Example
const fEx: HypercubeFunction = {
  numVars: 2,
  values: [2n, 3n, 4n, 5n]
};

const gEx: HypercubeFunction = {
  numVars: 2,
  values: [4n, 5n, 2n, 3n]
};

const hEx = computeRunningProduct(fEx, gEx);
const cube = booleanHypercube(2);

console.log("Example:");
for (let i = 0; i < cube.length; i++) {
  console.log(`  h(${cube[i]!.join(',')}) = ${hEx.values[i]!}`);
}

// Compute final product
let finalProd = 1n;
for (let i = 0; i < fEx.values.length; i++) {
  finalProd = mul(finalProd, mul(fEx.values[i]!, inv(gEx.values[i]!)));
}
console.log(`\nFull product prod f/g = ${finalProd}`);
console.log(`Products equal? ${finalProd === 1n ? 'YES' : 'NO'}\n`);

// ============================================================================
// Section 3: The Constraint Polynomial
// ============================================================================

console.log("--- Section 3: Constraint Polynomial ---\n");

/**
 * THE CONSTRAINT POLYNOMIAL C
 *
 * C(b) = h(next(b)) * g(b) - h(b) * f(b)
 *
 * For ProductCheck to be valid:
 *   C(b) = 0 for all b in B_n
 *
 * This is exactly the ZeroCheck condition!
 *
 * Note: We need h_shifted(b) = h(next(b)) in polynomial form.
 */

/**
 * Computes the constraint polynomial C on the hypercube.
 * C(b) = h_shifted(b) * g(b) - h(b) * f(b)
 */
function computeConstraint(
  f: HypercubeFunction,
  g: HypercubeFunction,
  h: HypercubeFunction
): HypercubeFunction {
  const hShifted = computeShiftedH(h);
  const values: bigint[] = [];

  for (let i = 0; i < f.values.length; i++) {
    // C(b) = h_shifted(b) * g(b) - h(b) * f(b)
    const term1 = mul(hShifted.values[i]!, g.values[i]!);
    const term2 = mul(h.values[i]!, f.values[i]!);
    values.push(sub(term1, term2));
  }

  return { numVars: f.numVars, values };
}

const constraintEx = computeConstraint(fEx, gEx, hEx);

console.log("Constraint C(b) = h(next(b))*g(b) - h(b)*f(b):\n");
for (let i = 0; i < cube.length; i++) {
  console.log(`  C(${cube[i]!.join(',')}) = ${constraintEx.values[i]!}`);
}

// Check all zeros
const allZeros = constraintEx.values.every(v => v === 0n);
console.log(`\nAll constraints zero? ${allZeros ? 'YES' : 'NO'}\n`);

// ============================================================================
// Section 4: ProductCheck via ZeroCheck
// ============================================================================

console.log("--- Section 4: ProductCheck Prover ---\n");

/**
 * PRODUCTCHECK = ZEROCHECK(C) + BOUNDARY
 *
 * The ProductCheck prover:
 * 1. Computes h (running product)
 * 2. Computes C (constraint polynomial)
 * 3. Runs ZeroCheck on C
 * 4. Proves boundary conditions
 *
 * We'll implement a simplified version that shows the structure.
 */

/**
 * Computes round polynomial for ZeroCheck on constraint C.
 *
 * The constraint involves h, h_shifted, f, g.
 * When combined with eq(x, r), the degree is 3 in each variable.
 */
function computeProductCheckRoundPoly(
  f: HypercubeFunction,
  g: HypercubeFunction,
  h: HypercubeFunction,
  r: bigint[],
  sumcheckChallenges: bigint[],
  round: number
): EvaluationPolynomial {
  const n = f.numVars;
  const remainingVars = n - round - 1;
  const remainingCube = booleanHypercube(remainingVars);

  // Degree is 3: h, h_shifted, f/g are multilinear, eq is multilinear
  // h * f, h_shifted * g -> degree 2
  // times eq -> degree 3
  const degree = 3;
  const evaluations: bigint[] = [];

  for (let j = 0; j <= degree; j++) {
    let sum = 0n;

    for (const suffix of remainingCube) {
      const x: bigint[] = [
        ...sumcheckChallenges,
        BigInt(j),
        ...suffix.map(v => BigInt(v))
      ];

      // Evaluate C(x) * eq(x, r)
      const constraint = evaluateConstraintMLE(f, g, h, x);
      const eqVal = eq(x, r);
      sum = add(sum, mul(constraint, eqVal));
    }

    evaluations.push(sum);
  }

  return { evaluations, degree };
}

/**
 * Evaluates the constraint C at any point using MLE.
 * C(x) = h_shifted_MLE(x) * g_MLE(x) - h_MLE(x) * f_MLE(x)
 */
function evaluateConstraintMLE(
  f: HypercubeFunction,
  g: HypercubeFunction,
  h: HypercubeFunction,
  x: bigint[]
): bigint {
  const hVal = evaluateMLE(h, x);
  const hShifted = computeShiftedH(h);
  const hShiftedVal = evaluateMLE(hShifted, x);
  const fVal = evaluateMLE(f, x);
  const gVal = evaluateMLE(g, x);

  return sub(mul(hShiftedVal, gVal), mul(hVal, fVal));
}

/**
 * ProductCheck prover.
 */
function productCheckProver(
  f: HypercubeFunction,
  g: HypercubeFunction,
  r: bigint[],
  getChallenge: () => bigint
): ProductCheckTranscript {
  const n = f.numVars;

  // Step 1: Compute auxiliary h
  const h = computeRunningProduct(f, g);

  // Step 2: Run ZeroCheck on constraint
  const roundPolynomials: EvaluationPolynomial[] = [];
  const sumcheckChallenges: bigint[] = [];

  for (let round = 0; round < n; round++) {
    const gi = computeProductCheckRoundPoly(f, g, h, r, sumcheckChallenges, round);
    roundPolynomials.push(gi);
    sumcheckChallenges.push(getChallenge());
  }

  // Final evaluation
  const evaluationPoint = sumcheckChallenges;
  const hShifted = computeShiftedH(h);

  return {
    hValues: h.values,
    zeroCheckChallenge: r,
    roundPolynomials,
    sumcheckChallenges,
    evaluationPoint,
    finalValues: {
      h: evaluateMLE(h, evaluationPoint),
      hShifted: evaluateMLE(hShifted, evaluationPoint),
      f: evaluateMLE(f, evaluationPoint),
      g: evaluateMLE(g, evaluationPoint)
    }
  };
}

console.log("ProductCheck Prover:");
console.log("  1. Compute h (running product)");
console.log("  2. Run ZeroCheck on C = h_shifted*g - h*f");
console.log("  3. Output transcript with h values and ZeroCheck proof\n");

// ============================================================================
// Section 5: ProductCheck Verifier
// ============================================================================

console.log("--- Section 5: ProductCheck Verifier ---\n");

/**
 * ProductCheck verifier.
 */
function productCheckVerifier(
  transcript: ProductCheckTranscript,
  oracleF: (x: bigint[]) => bigint,
  oracleG: (x: bigint[]) => bigint
): boolean {
  const n = transcript.roundPolynomials.length;
  const r = transcript.zeroCheckChallenge;

  console.log(`  Verifying ProductCheck with ${n} variables...`);
  console.log(`  Initial challenge r = (${r.join(', ')})\n`);

  // Check 1: h(first) = 1
  const hFirst = transcript.hValues[0]!;
  console.log(`  Boundary check: h(0...0) = ${hFirst} (should be 1)`);
  if (hFirst !== 1n) {
    console.log(`  FAILED: boundary\n`);
    return false;
  }
  console.log(`  PASSED\n`);

  // Check 2: ZeroCheck on constraint (claimed sum = 0)
  let currentClaim = 0n;

  for (let round = 0; round < n; round++) {
    const gi = transcript.roundPolynomials[round]!;
    const si = transcript.sumcheckChallenges[round]!;

    const g0 = gi.evaluations[0]!;
    const g1 = gi.evaluations[1]!;
    const sum = add(g0, g1);

    console.log(`  Round ${round + 1}:`);
    console.log(`    g(0) + g(1) = ${sum}, expected = ${currentClaim}`);

    if (sum !== currentClaim) {
      console.log(`    FAILED: sum check\n`);
      return false;
    }
    console.log(`    PASSED`);

    currentClaim = evaluateFromEvaluations(gi.evaluations, si);
    console.log(`    s_${round + 1} = ${si}, new claim = ${currentClaim}\n`);
  }

  // Check 3: Final evaluation
  // C(s) * eq(s, r) should equal currentClaim
  const s = transcript.evaluationPoint;
  const eqSR = eq(s, r);
  const fAtS = oracleF(s);
  const gAtS = oracleG(s);

  // We need h(s) and h_shifted(s) - these come from commitment openings
  // For now, we trust the prover's values
  const hAtS = transcript.finalValues.h;
  const hShiftedAtS = transcript.finalValues.hShifted;

  // C(s) = h_shifted(s) * g(s) - h(s) * f(s)
  const constraintAtS = sub(mul(hShiftedAtS, gAtS), mul(hAtS, fAtS));
  const expectedFinal = mul(constraintAtS, eqSR);

  console.log(`  Final check:`);
  console.log(`    h(s) = ${hAtS}`);
  console.log(`    h_shifted(s) = ${hShiftedAtS}`);
  console.log(`    f(s) = ${fAtS}`);
  console.log(`    g(s) = ${gAtS}`);
  console.log(`    C(s) = ${constraintAtS}`);
  console.log(`    eq(s, r) = ${eqSR}`);
  console.log(`    Expected: ${expectedFinal}, Got: ${currentClaim}`);

  if (expectedFinal !== currentClaim) {
    console.log(`    FAILED: final evaluation\n`);
    return false;
  }
  console.log(`    PASSED\n`);

  return true;
}

console.log("ProductCheck Verifier:");
console.log("  1. Check h(first) = 1");
console.log("  2. Verify ZeroCheck rounds");
console.log("  3. Check final: C(s)*eq(s,r) = claimed value\n");

// ============================================================================
// Section 6: Complete Example
// ============================================================================

console.log("--- Section 6: Complete Example ---\n");

console.log("Testing ProductCheck with equal products...\n");

// Deterministic challenges
let challengeIdx = 0;
const challengeValues = [11n, 13n];
function getChallenge(): bigint {
  return challengeValues[challengeIdx++]!;
}

// Initial ZeroCheck challenge
const rInit = [5n, 7n];

// Run prover
console.log("=== PROVER ===\n");
const transcript = productCheckProver(fEx, gEx, rInit, getChallenge);

console.log("h values:");
for (let i = 0; i < cube.length; i++) {
  console.log(`  h(${cube[i]!.join(',')}) = ${transcript.hValues[i]!}`);
}
console.log();

// Run verifier
console.log("=== VERIFIER ===\n");
const oracle = {
  f: (x: bigint[]) => evaluateMLE(fEx, x),
  g: (x: bigint[]) => evaluateMLE(gEx, x)
};

const accepted = productCheckVerifier(transcript, oracle.f, oracle.g);
console.log(`Verification result: ${accepted ? 'ACCEPTED' : 'REJECTED'}\n`);

// ============================================================================
// Section 7: Example with Unequal Products
// ============================================================================

console.log("--- Section 7: Detecting Unequal Products ---\n");

// f and g with different products
const fBad: HypercubeFunction = {
  numVars: 2,
  values: [2n, 3n, 4n, 5n]  // prod = 120
};

const gBad: HypercubeFunction = {
  numVars: 2,
  values: [1n, 1n, 1n, 1n]  // prod = 1
};

console.log("Testing with UNEQUAL products...");
let fBadProd = 1n;
let gBadProd = 1n;
for (let i = 0; i < fBad.values.length; i++) {
  fBadProd = mul(fBadProd, fBad.values[i]!);
  gBadProd = mul(gBadProd, gBad.values[i]!);
}
console.log(`  prod f = ${fBadProd}`);
console.log(`  prod g = ${gBadProd}\n`);

// Reset challenges
challengeIdx = 0;

// Run prover
console.log("=== CHEATING PROVER ===\n");
const badTranscript = productCheckProver(fBad, gBad, rInit, getChallenge);

// Check h values
console.log("h values (should NOT end at 1):");
for (let i = 0; i < cube.length; i++) {
  console.log(`  h(${cube[i]!.join(',')}) = ${badTranscript.hValues[i]!}`);
}
console.log();

// Run verifier
console.log("=== VERIFIER ===\n");
const badOracle = {
  f: (x: bigint[]) => evaluateMLE(fBad, x),
  g: (x: bigint[]) => evaluateMLE(gBad, x)
};

const badAccepted = productCheckVerifier(badTranscript, badOracle.f, badOracle.g);
console.log(`Verification result: ${badAccepted ? 'ACCEPTED (BAD!)' : 'REJECTED (GOOD!)'}\n`);

// ============================================================================
// Section 8: Production Stubs
// ============================================================================

console.log("--- Section 8: Production Stubs ---\n");

/**
 * Optimized ProductCheck prover with O(n * 2^n) time.
 *
 * Uses bookkeeping to avoid recomputation.
 *
 * @throws NotImplementedError - Optimization not yet implemented
 */
export function productCheckProverOptimized(
  f: HypercubeFunction,
  g: HypercubeFunction,
  r: bigint[],
  getChallenge: () => bigint
): ProductCheckTranscript {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: productCheckProverOptimized - ' +
    'Requires bookkeeping for O(n * 2^n) prover time'
  );
}

/**
 * ProductCheck with polynomial commitment openings.
 *
 * Instead of sending h values, prover commits to h and
 * provides opening proofs at the evaluation point.
 *
 * @throws NotImplementedError - Commitment integration not implemented
 */
export function productCheckWithCommitments(
  fCommitment: unknown,
  gCommitment: unknown,
  r: bigint[],
  getChallenge: () => bigint
): ProductCheckTranscript {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: productCheckWithCommitments - ' +
    'ProductCheck with polynomial commitment scheme'
  );
}

/**
 * Batch ProductCheck for multiple (f_i, g_i) pairs.
 *
 * Uses random linear combination to batch multiple checks.
 *
 * @throws NotImplementedError - Batch optimization not implemented
 */
export function batchProductCheck(
  fs: HypercubeFunction[],
  gs: HypercubeFunction[],
  combinationCoeffs: bigint[],
  r: bigint[],
  getChallenge: () => bigint
): ProductCheckTranscript {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: batchProductCheck - ' +
    'Batch ProductCheck for multiple polynomial pairs'
  );
}

console.log("Production stubs defined:");
console.log("  - productCheckProverOptimized");
console.log("  - productCheckWithCommitments");
console.log("  - batchProductCheck\n");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");

console.log("ProductCheck Protocol:");
console.log("  Claim: prod_b f(b) = prod_b g(b)");
console.log("  Method: Running product h + ZeroCheck on constraint\n");

console.log("Implementation:");
console.log("  - computeRunningProduct: builds auxiliary h");
console.log("  - computeConstraint: C = h_shifted*g - h*f");
console.log("  - productCheckProver: runs ZeroCheck on C");
console.log("  - productCheckVerifier: checks all conditions\n");

console.log("Complexity:");
console.log("  - Prover: O(n * 2^n) with optimization");
console.log("  - Verifier: O(n) + oracle queries");
console.log("  - Communication: O(n) field elements\n");

console.log("Next: Multiset and Permutation checks");

// Exports
export {
  computeRunningProduct,
  computeShiftedH,
  computeConstraint,
  productCheckProver,
  productCheckVerifier,
  type HypercubeFunction,
  type ProductCheckTranscript
};
