/**
 * @module sage/zk/sumcheck
 * @description Sumcheck Protocol Implementation
 *
 * The sumcheck protocol allows a prover to convince a verifier that
 * the sum of a multivariate polynomial over the boolean hypercube equals H.
 *
 * Protocol Overview:
 * -----------------
 * Given a multivariate polynomial f(x_1, ..., x_n) and claimed sum H, the protocol
 * proceeds in n rounds:
 *
 * Round i:
 *   1. Prover sends univariate polynomial p_i(X) = sum_{x_{i+1},...,x_n in {0,1}} f(r_1,...,r_{i-1},X,x_{i+1},...,x_n)
 *   2. Verifier checks: p_i(0) + p_i(1) == (previous round's claim)
 *   3. Verifier sends random challenge r_i
 *   4. New claim becomes p_i(r_i)
 *
 * Final Check:
 *   Verifier evaluates f(r_1, ..., r_n) and checks it equals the final claim.
 *
 * Used in: GKR, Lasso, Jolt, and many modern proof systems.
 */

import { ValueError } from '../errors.js';
import type { FiniteFieldElement, FiniteFieldPrime } from '../rings/finite_rings/finite_field_prime.js';
import type { MPolynomial } from '../rings/polynomial/multi_polynomial_element.js';
import type { Polynomial } from '../rings/polynomial/polynomial_element.js';
import { PolynomialRingConstructor } from '../rings/polynomial/polynomial_ring.js';
import { closestPowerOfTwo, intToBinary, multilinearExtension } from './multilinear.js';

/**
 * A sumcheck proof consists of the round polynomials, challenges, and final evaluation.
 */
export interface SumcheckProof {
  /** The univariate polynomials sent by the prover in each round */
  rounds: Polynomial<FiniteFieldElement>[];
  /** The random challenges sent by the verifier in each round */
  challenges: FiniteFieldElement[];
  /** The final evaluation of the polynomial at all challenges */
  finalEvaluation: FiniteFieldElement;
}

/**
 * Result of running the sumcheck protocol (for combined prover-verifier).
 */
export interface SumcheckResult {
  /** Whether the sumcheck verification passed */
  valid: boolean;
  /** The random challenges used */
  challenges: FiniteFieldElement[];
  /** The final sum after all rounds */
  finalSum: FiniteFieldElement;
}

/**
 * Generate the prover's polynomial for one round of sumcheck.
 *
 * This computes p(X) = sum_{x_{i+1},...,x_n in {0,1}} f(r_1,...,r_{i-1},X,x_{i+1},...,x_n)
 * where the challenges r_1,...,r_{i-1} have already been fixed.
 *
 * @param poly - The multivariate polynomial
 * @param challenges - The challenges from previous rounds
 * @param numVars - Total number of variables
 * @param field - The finite field
 * @returns The univariate polynomial for this round
 */
export function sumcheckRoundProver(
  poly: MPolynomial<FiniteFieldElement>,
  challenges: FiniteFieldElement[],
  numVars: number,
  field: FiniteFieldPrime
): Polynomial<FiniteFieldElement> {
  const freeVarIdx = challenges.length;
  const [uniRing] = PolynomialRingConstructor(field, 'x');

  // Substitute challenges for already-bound variables
  let partialPoly = poly;
  for (let i = 0; i < challenges.length; i++) {
    const evalPoint: Record<string, FiniteFieldElement> = {};
    evalPoint[`x${i}`] = challenges[i]!;
    partialPoly = partialPoly.subs(evalPoint);
  }

  // Sum over remaining variables (except the free variable)
  // Layer poly p(x_freeVarIdx) = sum_{x_{freeVarIdx+1}, ..., x_{n-1}} partialPoly
  let p0 = field.__call__(0n); // coefficient of x^0
  let p1 = field.__call__(0n); // coefficient of x^1

  const remainingVars = numVars - freeVarIdx - 1;
  const numPoints = 1 << remainingVars;

  for (let i = 0; i < numPoints; i++) {
    // Build evaluation point for remaining variables
    const bits = intToBinary(i, remainingVars);
    const evalPoint0: Record<string, FiniteFieldElement> = {};
    const evalPoint1: Record<string, FiniteFieldElement> = {};

    evalPoint0[`x${freeVarIdx}`] = field.__call__(0n);
    evalPoint1[`x${freeVarIdx}`] = field.__call__(1n);

    for (let j = 0; j < remainingVars; j++) {
      const varName = `x${freeVarIdx + 1 + j}`;
      evalPoint0[varName] = field.__call__(BigInt(bits[j]!));
      evalPoint1[varName] = field.__call__(BigInt(bits[j]!));
    }

    const val0 = partialPoly.evaluate(evalPoint0) as FiniteFieldElement;
    const val1 = partialPoly.evaluate(evalPoint1) as FiniteFieldElement;

    p0 = p0.add(val0);
    p1 = p1.add(val1.sub(val0));
  }

  return uniRing.__call__([p0, p1]);
}

/**
 * Verify one round of the sumcheck protocol.
 *
 * Checks that p(0) + p(1) equals the expected sum, then computes the
 * new sum as p(challenge).
 *
 * @param expectedSum - The sum that p(0) + p(1) should equal
 * @param layerPoly - The prover's polynomial for this round
 * @param field - The finite field
 * @param challenge - Optional challenge (if not provided, a random one is generated)
 * @returns The challenge used and the new sum
 * @throws {ValueError} If the consistency check fails
 */
export function sumcheckRoundVerifier(
  expectedSum: FiniteFieldElement,
  layerPoly: Polynomial<FiniteFieldElement>,
  field: FiniteFieldPrime,
  challenge?: FiniteFieldElement
): { challenge: FiniteFieldElement; newSum: FiniteFieldElement } {
  // Consistency check: p(0) + p(1) == expectedSum
  const p0 = layerPoly.evaluate(field.__call__(0n));
  const p1 = layerPoly.evaluate(field.__call__(1n));
  const computedSum = p0.add(p1);

  if (!computedSum.eq(expectedSum)) {
    throw new ValueError(
      `Sumcheck failed: p(0) + p(1) = ${computedSum.value} != ${expectedSum.value}`
    );
  }

  // Use provided challenge or generate random one
  const r = challenge ?? field.random_element();

  return {
    challenge: r,
    newSum: layerPoly.evaluate(r),
  };
}

/**
 * Generate a sumcheck proof that the sum of a polynomial over the boolean hypercube
 * equals the claimed value.
 *
 * @param poly - The multivariate polynomial
 * @param claimedSum - The claimed sum
 * @param numVars - The number of variables
 * @param field - The finite field
 * @param challengeGenerator - Optional function to generate challenges (for testing)
 * @returns The sumcheck proof
 *
 * @example
 * ```typescript
 * const F = GF(101n);
 * const values = [F.__call__(9n), F.__call__(2n), F.__call__(5n), F.__call__(4n)];
 * const poly = multilinearExtension(values, F);
 * const sum = values.reduce((a, b) => a.add(b), F.__call__(0n));
 *
 * const proof = sumcheckProve(poly, sum, 2, F);
 * ```
 */
export function sumcheckProve(
  poly: MPolynomial<FiniteFieldElement>,
  claimedSum: FiniteFieldElement,
  numVars: number,
  field: FiniteFieldPrime,
  challengeGenerator?: () => FiniteFieldElement
): SumcheckProof {
  const rounds: Polynomial<FiniteFieldElement>[] = [];
  const challenges: FiniteFieldElement[] = [];
  let currentSum = claimedSum;

  for (let round = 0; round < numVars; round++) {
    // Prover generates layer polynomial
    const layerPoly = sumcheckRoundProver(poly, challenges, numVars, field);
    rounds.push(layerPoly);

    // Verifier checks and produces challenge
    const r = challengeGenerator ? challengeGenerator() : field.random_element();
    const result = sumcheckRoundVerifier(currentSum, layerPoly, field, r);

    challenges.push(result.challenge);
    currentSum = result.newSum;
  }

  // Compute final evaluation
  const evalPoint: Record<string, FiniteFieldElement> = {};
  challenges.forEach((c, i) => {
    evalPoint[`x${i}`] = c;
  });
  const finalEvaluation = poly.evaluate(evalPoint) as FiniteFieldElement;

  return {
    rounds,
    challenges,
    finalEvaluation,
  };
}

/**
 * Verify a sumcheck proof.
 *
 * @param proof - The sumcheck proof to verify
 * @param claimedSum - The claimed sum
 * @param polyEvaluator - A function that evaluates the polynomial at a given point
 * @param field - The finite field
 * @param degreeCheck - Maximum allowed degree for round polynomials (default: 1)
 * @returns True if the proof is valid, false otherwise
 *
 * @example
 * ```typescript
 * const valid = sumcheckVerify(
 *   proof,
 *   claimedSum,
 *   (point) => poly.evaluate(Object.fromEntries(point.map((v, i) => [`x${i}`, v])))
 * );
 * ```
 */
export function sumcheckVerify(
  proof: SumcheckProof,
  claimedSum: FiniteFieldElement,
  polyEvaluator: (point: FiniteFieldElement[]) => FiniteFieldElement,
  field: FiniteFieldPrime,
  degreeCheck: number = 1
): boolean {
  const { rounds, challenges, finalEvaluation } = proof;

  if (rounds.length !== challenges.length) {
    return false;
  }

  let expectedSum = claimedSum;

  // Check each round
  for (let i = 0; i < rounds.length; i++) {
    const layerPoly = rounds[i]!;
    const challenge = challenges[i]!;

    // Degree check
    if (layerPoly.degree() > degreeCheck) {
      return false;
    }

    // Consistency check: p(0) + p(1) == expectedSum
    const p0 = layerPoly.evaluate(field.__call__(0n));
    const p1 = layerPoly.evaluate(field.__call__(1n));
    const computedSum = p0.add(p1);

    if (!computedSum.eq(expectedSum)) {
      return false;
    }

    // Update expected sum for next round
    expectedSum = layerPoly.evaluate(challenge);
  }

  // Final check: polynomial evaluation at challenges should match last round's claim
  const actualFinalEval = polyEvaluator(challenges);

  if (!actualFinalEval.eq(expectedSum)) {
    return false;
  }

  // Also verify the proof's stored final evaluation
  if (!actualFinalEval.eq(finalEvaluation)) {
    return false;
  }

  return true;
}

/**
 * Run the complete sumcheck protocol (combined prover-verifier for testing).
 *
 * This function computes the multilinear extension of the given values,
 * then runs the sumcheck protocol to verify that the sum over the boolean
 * hypercube equals the sum of the input values.
 *
 * @param values - The function values to sum
 * @param field - The finite field
 * @param degreeCheck - Maximum allowed degree for round polynomials (default: 1)
 * @returns The result of the sumcheck protocol
 *
 * @example
 * ```typescript
 * const F = GF(101n);
 * const values = [F.__call__(9n), F.__call__(2n), F.__call__(5n), F.__call__(4n)];
 * const result = sumcheckRun(values, F);
 * console.log(result.valid); // true
 * ```
 */
export function sumcheckRun(
  values: FiniteFieldElement[],
  field: FiniteFieldPrime,
  degreeCheck: number = 1
): SumcheckResult {
  if (values.length === 0) {
    throw new ValueError('values cannot be empty');
  }

  const poly = multilinearExtension(values, field);
  const numVars = closestPowerOfTwo(values.length);

  // Compute expected sum
  let expectedSum = field.__call__(0n);
  for (const v of values) {
    expectedSum = expectedSum.add(v);
  }

  // If we have fewer values than 2^numVars, pad with zeros (implicit)
  // The MLE will evaluate to 0 at the missing points

  const challenges: FiniteFieldElement[] = [];

  for (let round = 0; round < numVars; round++) {
    // Prover generates layer polynomial
    const layerPoly = sumcheckRoundProver(poly, challenges, numVars, field);

    // Degree check
    if (layerPoly.degree() > degreeCheck) {
      return {
        valid: false,
        challenges,
        finalSum: expectedSum,
      };
    }

    // Verifier checks and produces challenge
    try {
      const result = sumcheckRoundVerifier(expectedSum, layerPoly, field);
      challenges.push(result.challenge);
      expectedSum = result.newSum;
    } catch {
      return {
        valid: false,
        challenges,
        finalSum: expectedSum,
      };
    }
  }

  // Final check: verify poly(challenges) == finalSum
  const evalPoint: Record<string, FiniteFieldElement> = {};
  challenges.forEach((c, i) => {
    evalPoint[`x${i}`] = c;
  });
  const finalEval = poly.evaluate(evalPoint) as FiniteFieldElement;

  const valid = finalEval.eq(expectedSum);

  return {
    valid,
    challenges,
    finalSum: expectedSum,
  };
}

/**
 * Utility function to create a polynomial evaluator from an MPolynomial.
 *
 * @param poly - The multivariate polynomial
 * @returns A function that evaluates the polynomial at a given point
 */
export function createPolyEvaluator(
  poly: MPolynomial<FiniteFieldElement>
): (point: FiniteFieldElement[]) => FiniteFieldElement {
  return (point: FiniteFieldElement[]) => {
    const evalPoint: Record<string, FiniteFieldElement> = {};
    point.forEach((v, i) => {
      evalPoint[`x${i}`] = v;
    });
    return poly.evaluate(evalPoint) as FiniteFieldElement;
  };
}
