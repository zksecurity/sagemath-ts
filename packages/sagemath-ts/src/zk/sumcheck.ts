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
import type { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';
import { PolynomialRingConstructor } from '../rings/polynomial/polynomial_ring.js';
import { booleanHypercube, multilinearExtension } from './multilinear.js';

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
 * Return the variable names of a multivariate polynomial's parent ring.
 *
 * The blueprint always addresses variables through `poly.args()`
 * (`sumcheck.sage:90,101,114`); this port addresses them by name because
 * `MPolynomial.subs` / `MPolynomial.evaluate` are keyed by name. Using the
 * ring's actual names (rather than assuming `x0, x1, ...`) is what makes the
 * protocol work on a ring whose generators are called anything else.
 */
function argNames(poly: MPolynomial<FiniteFieldElement>): readonly string[] {
  return poly.parent.names;
}

/**
 * Convert a multivariate polynomial that involves at most the single variable
 * `varIdx` into an element of the univariate ring `uniRing`.
 */
function toUnivariate(
  poly: MPolynomial<FiniteFieldElement>,
  varIdx: number,
  uniRing: PolynomialRing<FiniteFieldElement>
): Polynomial<FiniteFieldElement> {
  const deg = poly.degreeIn(varIdx);
  if (deg < 0) {
    // zero polynomial
    return uniRing.zero();
  }

  const coeffs: FiniteFieldElement[] = [];
  for (let d = 0; d <= deg; d++) {
    const exponent = new Array<number>(poly.parent.ngens_value).fill(0);
    exponent[varIdx] = d;
    coeffs.push(poly.monomial_coefficient(exponent));
  }
  return uniRing.__call__(coeffs);
}

/**
 * Generate the prover's polynomial for one round of sumcheck.
 *
 * This computes p(X) = sum_{x_{i+1},...,x_n in {0,1}} f(r_1,...,r_{i-1},X,x_{i+1},...,x_n)
 * where the challenges r_1,...,r_{i-1} have already been fixed.
 *
 * Port of `sumcheck_round_prover` (`reference/sage_blueprints/sumcheck.sage:82-140`):
 * the round polynomial is built *symbolically*, by substituting the challenges
 * and then summing the partially evaluated polynomial over the boolean
 * hypercube of the trailing variables. It therefore has whatever degree the
 * polynomial actually has in the free variable -- the blueprint explicitly
 * disables the "layer polynomial is linear" assertion because higher degrees
 * are legitimate (GKR).
 *
 * @param poly - The multivariate polynomial
 * @param challenges - The challenges from previous rounds
 * @param numVars - Total number of variables (must equal `poly.args().length`)
 * @param field - The finite field
 * @returns The univariate polynomial for this round
 * @throws {ValueError} If `numVars` disagrees with the polynomial's ring, if
 *   there is no free variable left, or if the round polynomial depends on a
 *   variable other than the free one
 */
export function sumcheckRoundProver(
  poly: MPolynomial<FiniteFieldElement>,
  challenges: FiniteFieldElement[],
  numVars: number,
  field: FiniteFieldPrime
): Polynomial<FiniteFieldElement> {
  const names = argNames(poly);

  // The blueprint derives everything from poly.args(); an inconsistent
  // numVars would silently evaluate missing variables at 0.
  if (numVars !== names.length) {
    throw new ValueError(
      `prover: numVars must equal the number of variables of poly (got ${numVars} vs ${names.length})`
    );
  }

  const freeVarIdx = challenges.length;
  if (freeVarIdx >= numVars) {
    throw new ValueError('prover: no free variable left');
  }

  const [uniRing] = PolynomialRingConstructor(field, 'x');

  // Substitute challenges for the already-bound variables r_1, ..., r_n
  let partialPoly = poly;
  if (challenges.length > 0) {
    const substitution: Record<string, FiniteFieldElement> = {};
    for (let i = 0; i < challenges.length; i++) {
      substitution[names[i]!] = challenges[i]!;
    }
    partialPoly = poly.subs(substitution);
  }

  // Sum over the boolean hypercube of the trailing variables (a_0, ..., a_m)
  const toBeEvaluated = names.slice(freeVarIdx + 1);
  let res = poly.parent.zero();
  for (const evals of booleanHypercube(toBeEvaluated.length)) {
    const substitution: Record<string, FiniteFieldElement> = {};
    for (let j = 0; j < toBeEvaluated.length; j++) {
      substitution[toBeEvaluated[j]!] = field.__call__(BigInt(evals[j]!));
    }
    res = res.add(partialPoly.subs(substitution));
  }

  // Sanity check (blueprint `sumcheck.sage:126-130`): the free variable is the
  // only variable left, and it is the right one.
  //
  // The blueprint spells this as `len(res.variables()) != 1`, which also
  // rejects a *constant* round polynomial. That is too strict: the zero
  // function has zero round polynomials in every round, and running the
  // blueprint on it does not raise its own error but crashes with
  // `AttributeError: 'IntegerMod_int' object has no attribute 'variables'`
  // (verified in SageMath). We therefore only reject a round polynomial that
  // depends on the *wrong* variable.
  //
  // @see Deviation: Sumcheck constant round polynomial
  for (let i = 0; i < numVars; i++) {
    if (i !== freeVarIdx && res.degreeIn(i) > 0) {
      throw new ValueError('prover: Layer polynomial is not built from the correct variable');
    }
  }

  return toUnivariate(res, freeVarIdx, uniRing);
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
  const names = argNames(poly);
  if (numVars !== names.length) {
    throw new ValueError(
      `numVars must equal the number of variables of poly (got ${numVars} vs ${names.length})`
    );
  }

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
    evalPoint[names[i]!] = c;
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
 * The number of rounds is dictated by the *statement*, not by the proof: the
 * blueprint runs `num_rounds = len(poly.args())` rounds
 * (`sumcheck.sage:211`). A verifier that trusts the prover's round count
 * accepts the empty proof for the false claim "sum = f(0, ..., 0)", so
 * `numVars` is a required argument here and a proof with a different number of
 * rounds is rejected outright.
 *
 * @param proof - The sumcheck proof to verify
 * @param claimedSum - The claimed sum
 * @param polyEvaluator - A function that evaluates the polynomial at a given point
 * @param field - The finite field
 * @param numVars - The number of variables of the polynomial, i.e. the number
 *   of rounds the proof must contain
 * @param options.degreeCheck - Maximum allowed degree for round polynomials.
 *   Defaults to no check, matching the blueprint's `degree_checks=None`
 *   (`sumcheck.sage:199`); pass 1 for a multilinear statement.
 * @returns True if the proof is valid, false otherwise
 *
 * @example
 * ```typescript
 * const valid = sumcheckVerify(
 *   proof,
 *   claimedSum,
 *   createPolyEvaluator(poly),
 *   F,
 *   2,
 *   { degreeCheck: 1 }
 * );
 * ```
 */
export function sumcheckVerify(
  proof: SumcheckProof,
  claimedSum: FiniteFieldElement,
  polyEvaluator: (point: FiniteFieldElement[]) => FiniteFieldElement,
  field: FiniteFieldPrime,
  numVars: number,
  options?: { degreeCheck?: number }
): boolean {
  const { rounds, challenges, finalEvaluation } = proof;
  const degreeCheck = options?.degreeCheck;

  if (rounds.length !== challenges.length) {
    return false;
  }

  // The proof must run exactly as many rounds as the polynomial has variables.
  if (rounds.length !== numVars) {
    return false;
  }

  let expectedSum = claimedSum;

  // Check each round
  for (let i = 0; i < rounds.length; i++) {
    const layerPoly = rounds[i]!;
    const challenge = challenges[i]!;

    // Degree check
    if (degreeCheck !== undefined && layerPoly.degree() > degreeCheck) {
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
 * Port of `sumcheck_run` / `sumcheck_run_poly`
 * (`reference/sage_blueprints/sumcheck.sage:191-236`), including the initial
 * sanity check that the polynomial really does sum to the claimed value over
 * the hypercube of *all* of its variables. Note that the multilinear extension
 * of a single value lives in a one-variable ring (`mle.sage:51-54`), so its sum
 * over {0,1} is twice the value and the initial check fails -- exactly as in
 * the blueprint.
 *
 * @param values - The function values to sum
 * @param field - The finite field
 * @param degreeCheck - Maximum allowed degree for round polynomials (default: 1,
 *   which is correct here because the polynomial is a multilinear extension)
 * @returns The result of the sumcheck protocol
 * @throws {ValueError} If `values` is empty, or if the polynomial does not sum
 *   to the expected value over the boolean hypercube
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
  const names = argNames(poly);
  // The round count is the polynomial's variable count (blueprint
  // `sumcheck.sage:211`), not ceil(log2(#values)): for a single value the MLE
  // still lives in a one-variable ring, and running zero rounds would "verify"
  // nothing at all.
  const numVars = names.length;

  // Compute expected sum
  let expectedSum = field.__call__(0n);
  for (const v of values) {
    expectedSum = expectedSum.add(v);
  }

  // If we have fewer values than 2^numVars, pad with zeros (implicit)
  // The MLE will evaluate to 0 at the missing points

  // Sanity check (blueprint `sumcheck.sage:205-208`): sum(poly) == expectedSum
  let hypercubeSum = field.__call__(0n);
  for (const point of booleanHypercube(numVars)) {
    const evalPoint: Record<string, FiniteFieldElement> = {};
    for (let i = 0; i < numVars; i++) {
      evalPoint[names[i]!] = field.__call__(BigInt(point[i]!));
    }
    hypercubeSum = hypercubeSum.add(poly.evaluate(evalPoint) as FiniteFieldElement);
  }
  if (!hypercubeSum.eq(expectedSum)) {
    throw new ValueError('Sumcheck failed: initial check failed');
  }

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
    evalPoint[names[i]!] = c;
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
  const names = argNames(poly);
  return (point: FiniteFieldElement[]) => {
    if (point.length !== names.length) {
      throw new ValueError(
        `point must have ${names.length} coordinates (got ${point.length})`
      );
    }
    const evalPoint: Record<string, FiniteFieldElement> = {};
    point.forEach((v, i) => {
      evalPoint[names[i]!] = v;
    });
    return poly.evaluate(evalPoint) as FiniteFieldElement;
  };
}
