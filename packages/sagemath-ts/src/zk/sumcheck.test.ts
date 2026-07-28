/**
 * Unit tests for the sumcheck protocol
 */
import { describe, expect, test } from 'bun:test';
import { FiniteFieldPrime } from '../rings/finite_rings/finite_field_prime.js';
import { MPolynomialRingConstructor } from '../rings/polynomial/multi_polynomial_ring.js';
import { PolynomialRingConstructor } from '../rings/polynomial/polynomial_ring.js';
import {
  closestPowerOfTwo,
  intToBinary,
  binaryToInt,
  booleanHypercube,
  eqPolynomial,
  multilinearExtension,
  sparseMultilinearExtension,
  isMultilinear,
} from './multilinear.js';
import {
  sumcheckProve,
  sumcheckVerify,
  sumcheckRun,
  sumcheckRoundProver,
  sumcheckRoundVerifier,
  createPolyEvaluator,
} from './sumcheck.js';

// Test field: a large prime for realistic testing
const F = new FiniteFieldPrime(101n);

// Helper to create field elements
const fe = (n: number | bigint) => F.__call__(n);

describe('Multilinear utilities', () => {
  describe('closestPowerOfTwo', () => {
    test('returns 0 for n <= 1', () => {
      expect(closestPowerOfTwo(0)).toBe(0);
      expect(closestPowerOfTwo(1)).toBe(0);
    });

    test('returns correct values for powers of 2', () => {
      expect(closestPowerOfTwo(2)).toBe(1);
      expect(closestPowerOfTwo(4)).toBe(2);
      expect(closestPowerOfTwo(8)).toBe(3);
      expect(closestPowerOfTwo(16)).toBe(4);
    });

    test('returns ceiling for non-powers of 2', () => {
      expect(closestPowerOfTwo(3)).toBe(2);
      expect(closestPowerOfTwo(5)).toBe(3);
      expect(closestPowerOfTwo(6)).toBe(3);
      expect(closestPowerOfTwo(7)).toBe(3);
      expect(closestPowerOfTwo(9)).toBe(4);
      expect(closestPowerOfTwo(100)).toBe(7);
    });
  });

  describe('intToBinary', () => {
    test('converts integers to binary arrays', () => {
      expect(intToBinary(0, 4)).toEqual([0, 0, 0, 0]);
      expect(intToBinary(5, 4)).toEqual([0, 1, 0, 1]);
      expect(intToBinary(7, 3)).toEqual([1, 1, 1]);
      expect(intToBinary(5, 8)).toEqual([0, 0, 0, 0, 0, 1, 0, 1]);
    });

    test('throws for negative values', () => {
      expect(() => intToBinary(-1, 4)).toThrow();
    });
  });

  describe('binaryToInt', () => {
    test('converts binary arrays to integers', () => {
      expect(binaryToInt([0, 0, 0, 0])).toBe(0n);
      expect(binaryToInt([0, 1, 0, 1])).toBe(5n);
      expect(binaryToInt([1, 1, 1])).toBe(7n);
      expect(binaryToInt([])).toBe(0n);
    });

    test('roundtrip with intToBinary', () => {
      for (let i = 0; i < 16; i++) {
        expect(binaryToInt(intToBinary(i, 4))).toBe(BigInt(i));
      }
    });
  });

  describe('booleanHypercube', () => {
    test('generates all points in {0,1}^n', () => {
      expect(booleanHypercube(0)).toEqual([[]]);
      expect(booleanHypercube(1)).toEqual([[0], [1]]);
      expect(booleanHypercube(2)).toEqual([
        [0, 0],
        [0, 1],
        [1, 0],
        [1, 1],
      ]);
      expect(booleanHypercube(3).length).toBe(8);
    });
  });

  describe('eqPolynomial', () => {
    test('evaluates to 1 at matching point', () => {
      const [R] = MPolynomialRingConstructor(F, ['x0', 'x1']);
      const vars = R.gens();

      // eq([0, 1], x) should be 1 at (0, 1)
      const eq01 = eqPolynomial([0, 1], vars, R);
      const result = eq01.evaluate({ x0: fe(0), x1: fe(1) });
      expect(result.eq(1)).toBe(true);
    });

    test('evaluates to 0 at non-matching points', () => {
      const [R] = MPolynomialRingConstructor(F, ['x0', 'x1']);
      const vars = R.gens();

      const eq01 = eqPolynomial([0, 1], vars, R);

      // Should be 0 at all other boolean hypercube points
      expect(eq01.evaluate({ x0: fe(0), x1: fe(0) }).eq(0)).toBe(true);
      expect(eq01.evaluate({ x0: fe(1), x1: fe(0) }).eq(0)).toBe(true);
      expect(eq01.evaluate({ x0: fe(1), x1: fe(1) }).eq(0)).toBe(true);
    });

    test('works with integer base (auto bit decomposition)', () => {
      const [R] = MPolynomialRingConstructor(F, ['x0', 'x1', 'x2']);
      const vars = R.gens();

      // 5 in binary (3 bits) is [1, 0, 1]
      const eq5 = eqPolynomial(5, vars, R);
      expect(eq5.evaluate({ x0: fe(1), x1: fe(0), x2: fe(1) }).eq(1)).toBe(true);
      expect(eq5.evaluate({ x0: fe(0), x1: fe(0), x2: fe(0) }).eq(0)).toBe(true);
    });
  });

  describe('multilinearExtension', () => {
    test('interpolates values on boolean hypercube', () => {
      const values = [fe(1), fe(2), fe(3), fe(4)];
      const poly = multilinearExtension(values, F);

      // Evaluate at all boolean hypercube points
      expect(poly.evaluate({ x0: fe(0), x1: fe(0) }).eq(1)).toBe(true);
      expect(poly.evaluate({ x0: fe(0), x1: fe(1) }).eq(2)).toBe(true);
      expect(poly.evaluate({ x0: fe(1), x1: fe(0) }).eq(3)).toBe(true);
      expect(poly.evaluate({ x0: fe(1), x1: fe(1) }).eq(4)).toBe(true);
    });

    test('handles non-power-of-2 value counts', () => {
      // 3 values requires 2 variables (ceil(log2(3)) = 2)
      const values = [fe(9), fe(2), fe(5)];
      const poly = multilinearExtension(values, F);

      // First 3 values should match
      expect(poly.evaluate({ x0: fe(0), x1: fe(0) }).eq(9)).toBe(true);
      expect(poly.evaluate({ x0: fe(0), x1: fe(1) }).eq(2)).toBe(true);
      expect(poly.evaluate({ x0: fe(1), x1: fe(0) }).eq(5)).toBe(true);
      // 4th point should be 0 (implicit padding)
      expect(poly.evaluate({ x0: fe(1), x1: fe(1) }).eq(0)).toBe(true);
    });

    test('single value produces constant polynomial', () => {
      const values = [fe(42)];
      const poly = multilinearExtension(values, F);

      // Should be constant everywhere
      expect(poly.evaluate({ x0: fe(0) }).eq(42)).toBe(true);
      expect(poly.evaluate({ x0: fe(1) }).eq(42)).toBe(true);
      expect(poly.evaluate({ x0: fe(50) }).eq(42)).toBe(true);
    });

    test('result is multilinear', () => {
      const values = [fe(1), fe(2), fe(3), fe(4), fe(5), fe(6), fe(7), fe(8)];
      const poly = multilinearExtension(values, F);
      expect(isMultilinear(poly)).toBe(true);
    });
  });

  describe('sparseMultilinearExtension', () => {
    test('selector polynomial is 1 at specified indices', () => {
      const [R] = MPolynomialRingConstructor(F, ['x0', 'x1', 'x2']);
      const vars = R.gens();

      const selector = sparseMultilinearExtension([0, 2, 5], vars, R);

      // Should be 1 at indices 0, 2, 5
      expect(selector.evaluate({ x0: fe(0), x1: fe(0), x2: fe(0) }).eq(1)).toBe(true); // 0
      expect(selector.evaluate({ x0: fe(0), x1: fe(1), x2: fe(0) }).eq(1)).toBe(true); // 2
      expect(selector.evaluate({ x0: fe(1), x1: fe(0), x2: fe(1) }).eq(1)).toBe(true); // 5

      // Should be 0 at other indices
      expect(selector.evaluate({ x0: fe(0), x1: fe(0), x2: fe(1) }).eq(0)).toBe(true); // 1
      expect(selector.evaluate({ x0: fe(0), x1: fe(1), x2: fe(1) }).eq(0)).toBe(true); // 3
    });

    test('sparse values interpolation', () => {
      const [R] = MPolynomialRingConstructor(F, ['x0', 'x1']);
      const vars = R.gens();

      const sparseVals = new Map([
        [0, fe(9)],
        [2, fe(5)],
      ]);
      const poly = sparseMultilinearExtension(sparseVals, vars, R);

      expect(poly.evaluate({ x0: fe(0), x1: fe(0) }).eq(9)).toBe(true); // 0
      expect(poly.evaluate({ x0: fe(1), x1: fe(0) }).eq(5)).toBe(true); // 2
      expect(poly.evaluate({ x0: fe(0), x1: fe(1) }).eq(0)).toBe(true); // 1
      expect(poly.evaluate({ x0: fe(1), x1: fe(1) }).eq(0)).toBe(true); // 3
    });
  });
});

describe('Sumcheck Protocol', () => {
  describe('sumcheckRun (combined prover-verifier)', () => {
    test('passes for valid inputs', () => {
      const values = [fe(9), fe(2), fe(5), fe(4)];
      const result = sumcheckRun(values, F);

      expect(result.valid).toBe(true);
      expect(result.challenges.length).toBe(2); // 2 variables
    });

    test('rejects a single value: the MLE lives in a one-variable ring', () => {
      // The MLE of one value is the constant c in a ring with one variable
      // (mle.sage:51-54), so its sum over {0, 1} is 2c, not c. The blueprint
      // raises "Sumcheck failed: initial check failed" here -- verified by
      // running sumcheck_run([F(42)]) in SageMath. The previous expectation
      // (valid === true with zero rounds) came from deriving the round count
      // as ceil(log2(1)) = 0, which "verified" nothing at all.
      expect(() => sumcheckRun([fe(42)], F)).toThrow('Sumcheck failed: initial check failed');
    });

    test('passes for non-power-of-2 value counts', () => {
      const values = [fe(9), fe(2), fe(5), fe(4), fe(5), fe(7)];
      const result = sumcheckRun(values, F);

      expect(result.valid).toBe(true);
      expect(result.challenges.length).toBe(3); // 3 variables
    });

    test('passes for 8 values', () => {
      const values = [fe(1), fe(2), fe(3), fe(4), fe(5), fe(6), fe(7), fe(8)];
      const result = sumcheckRun(values, F);

      expect(result.valid).toBe(true);
      expect(result.challenges.length).toBe(3);
    });

    test('throws for empty values', () => {
      expect(() => sumcheckRun([], F)).toThrow();
    });
  });

  describe('sumcheckProve and sumcheckVerify', () => {
    test('prover generates valid proof', () => {
      const values = [fe(9), fe(2), fe(5), fe(4)];
      const poly = multilinearExtension(values, F);
      const claimedSum = values.reduce((a, b) => a.add(b), fe(0));

      const proof = sumcheckProve(poly, claimedSum, 2, F);

      expect(proof.rounds.length).toBe(2);
      expect(proof.challenges.length).toBe(2);
    });

    test('verifier accepts valid proof', () => {
      const values = [fe(9), fe(2), fe(5), fe(4)];
      const poly = multilinearExtension(values, F);
      const claimedSum = values.reduce((a, b) => a.add(b), fe(0));

      // Use deterministic challenges for testing
      let challengeIdx = 0;
      const challenges = [fe(7), fe(13)];
      const proof = sumcheckProve(poly, claimedSum, 2, F, () => challenges[challengeIdx++]!);

      const evaluator = createPolyEvaluator(poly);
      const valid = sumcheckVerify(proof, claimedSum, evaluator, F, 2, { degreeCheck: 1 });

      expect(valid).toBe(true);
    });

    test('verifier rejects proof with wrong claimed sum', () => {
      const values = [fe(9), fe(2), fe(5), fe(4)];
      const poly = multilinearExtension(values, F);
      const correctSum = values.reduce((a, b) => a.add(b), fe(0));
      const wrongSum = correctSum.add(fe(1));

      let challengeIdx = 0;
      const challenges = [fe(7), fe(13)];
      const proof = sumcheckProve(poly, correctSum, 2, F, () => challenges[challengeIdx++]!);

      const evaluator = createPolyEvaluator(poly);
      const valid = sumcheckVerify(proof, wrongSum, evaluator, F, 2, { degreeCheck: 1 });

      expect(valid).toBe(false);
    });

    test('verifier rejects proof with tampered final evaluation', () => {
      const values = [fe(9), fe(2), fe(5), fe(4)];
      const poly = multilinearExtension(values, F);
      const claimedSum = values.reduce((a, b) => a.add(b), fe(0));

      let challengeIdx = 0;
      const challenges = [fe(7), fe(13)];
      const proof = sumcheckProve(poly, claimedSum, 2, F, () => challenges[challengeIdx++]!);

      const tamperedProof = {
        ...proof,
        finalEvaluation: proof.finalEvaluation.add(fe(1)),
      };

      const evaluator = createPolyEvaluator(poly);
      const valid = sumcheckVerify(tamperedProof, claimedSum, evaluator, F, 2);

      expect(valid).toBe(false);
    });

    test('verifier rejects proof with a tampered round polynomial', () => {
      const values = [fe(9), fe(2), fe(5), fe(4)];
      const poly = multilinearExtension(values, F);
      const claimedSum = values.reduce((a, b) => a.add(b), fe(0));

      let challengeIdx = 0;
      const challenges = [fe(7), fe(13)];
      const proof = sumcheckProve(poly, claimedSum, 2, F, () => challenges[challengeIdx++]!);

      const [uniRing] = PolynomialRingConstructor(F, 'x');
      const evaluator = createPolyEvaluator(poly);

      // Round 0 replaced by a polynomial that still sums correctly but is a
      // different line: the final check must catch it.
      const p0 = proof.rounds[0]!;
      const shifted = p0.add(uniRing.__call__([fe(1), fe(99)])); // + (1 - 2x), sums to 0
      expect(shifted.evaluate(fe(0)).add(shifted.evaluate(fe(1))).eq(claimedSum)).toBe(true);
      expect(
        sumcheckVerify({ ...proof, rounds: [shifted, proof.rounds[1]!] }, claimedSum, evaluator, F, 2)
      ).toBe(false);

      // Round 1 replaced by a polynomial that does not sum to p0(r0).
      const bumped = proof.rounds[1]!.add(uniRing.__call__([fe(1)]));
      expect(
        sumcheckVerify({ ...proof, rounds: [proof.rounds[0]!, bumped] }, claimedSum, evaluator, F, 2)
      ).toBe(false);
    });

    test('verifier rejects a proof with the wrong number of rounds', () => {
      const values = [fe(9), fe(2), fe(5), fe(4)];
      const poly = multilinearExtension(values, F);
      const claimedSum = values.reduce((a, b) => a.add(b), fe(0));
      const evaluator = createPolyEvaluator(poly);

      // The empty-proof forgery: a prover who supplies no rounds at all and
      // claims the sum is f(0, 0). Trusting the proof's round count accepts it.
      const forged = {
        rounds: [],
        challenges: [],
        finalEvaluation: poly.evaluate({ x0: fe(0), x1: fe(0) }),
      };
      expect(forged.finalEvaluation.eq(9)).toBe(true);
      expect(sumcheckVerify(forged, fe(9), evaluator, F, 2)).toBe(false);

      // A truncated honest proof is rejected as well.
      let challengeIdx = 0;
      const challenges = [fe(7), fe(13)];
      const proof = sumcheckProve(poly, claimedSum, 2, F, () => challenges[challengeIdx++]!);
      expect(
        sumcheckVerify(
          { ...proof, rounds: [proof.rounds[0]!], challenges: [proof.challenges[0]!] },
          claimedSum,
          evaluator,
          F,
          2
        )
      ).toBe(false);

      // And so is a padded one.
      expect(
        sumcheckVerify(
          {
            ...proof,
            rounds: [...proof.rounds, proof.rounds[1]!],
            challenges: [...proof.challenges, fe(3)],
          },
          claimedSum,
          evaluator,
          F,
          2
        )
      ).toBe(false);
    });

    test('handles a polynomial of degree 2 in a variable', () => {
      // f = x0^2*x1 + x0 + 1 over GF(101). Verified in SageMath with the
      // blueprint: sumcheck_round_prover(f, []) is x0^2 + 2*x0 + 2 and
      // sum_{x in {0,1}^2} f(x) = 7.
      const [R, x0, x1] = MPolynomialRingConstructor(F, ['x0', 'x1']);
      const f = x0.pow(2n).mul(x1).add(x0).add(R.one());

      let claimedSum = fe(0);
      for (const [a, b] of [
        [0, 0],
        [0, 1],
        [1, 0],
        [1, 1],
      ]) {
        claimedSum = claimedSum.add(f.evaluate({ x0: fe(a!), x1: fe(b!) }));
      }
      expect(claimedSum.eq(7)).toBe(true);

      const round1 = sumcheckRoundProver(f, [], 2, F);
      expect(round1.toString()).toBe('x^2 + 2*x + 2');
      expect(round1.degree()).toBe(2);
      expect(round1.evaluate(fe(0)).add(round1.evaluate(fe(1))).eq(claimedSum)).toBe(true);

      let challengeIdx = 0;
      const challenges = [fe(7), fe(13)];
      const proof = sumcheckProve(f, claimedSum, 2, F, () => challenges[challengeIdx++]!);
      const evaluator = createPolyEvaluator(f);

      expect(sumcheckVerify(proof, claimedSum, evaluator, F, 2)).toBe(true);
      // A multilinear degree bound must reject this (legitimate) proof.
      expect(sumcheckVerify(proof, claimedSum, evaluator, F, 2, { degreeCheck: 1 })).toBe(false);
      expect(sumcheckVerify(proof, claimedSum, evaluator, F, 2, { degreeCheck: 2 })).toBe(true);
    });

    test('works on a ring whose variables are not named x0, x1, ...', () => {
      // Verified in SageMath: with vals = [9, 2, 5, 4] on PolynomialRing(GF(101), ['a','b']),
      // the interpolant is 6*a*b - 4*a - 7*b + 9 and the first round polynomial
      // is -2*a + 11 (i.e. 99*a + 11 mod 101).
      const [Rab, a, b] = MPolynomialRingConstructor(F, ['a', 'b']);
      const values = [fe(9), fe(2), fe(5), fe(4)];
      const poly = multilinearExtension(values, F, Rab, [a, b]);
      expect(poly.toString()).toBe('6*a*b + 97*a + 94*b + 9');

      const claimedSum = values.reduce((x, y) => x.add(y), fe(0));
      expect(claimedSum.eq(20)).toBe(true);

      const round1 = sumcheckRoundProver(poly, [], 2, F);
      expect(round1.toString()).toBe('99*x + 11');

      let challengeIdx = 0;
      const challenges = [fe(7), fe(13)];
      const proof = sumcheckProve(poly, claimedSum, 2, F, () => challenges[challengeIdx++]!);
      const evaluator = createPolyEvaluator(poly);
      expect(sumcheckVerify(proof, claimedSum, evaluator, F, 2, { degreeCheck: 1 })).toBe(true);
      expect(sumcheckVerify(proof, claimedSum.add(fe(1)), evaluator, F, 2)).toBe(false);
    });

    test('rejects a numVars that disagrees with the polynomial', () => {
      const values = [fe(9), fe(2), fe(5), fe(4)];
      const poly = multilinearExtension(values, F);
      const claimedSum = values.reduce((x, y) => x.add(y), fe(0));

      expect(() => sumcheckProve(poly, claimedSum, 3, F)).toThrow(
        'numVars must equal the number of variables of poly'
      );
      expect(() => sumcheckRoundProver(poly, [], 1, F)).toThrow(
        'numVars must equal the number of variables of poly'
      );
    });
  });

  describe('sumcheckRoundProver', () => {
    test('produces the exact round polynomial for an MLE', () => {
      // MLE of [1, 2, 3, 4] over GF(101) is 1 + 2*x0 + x1 (see the
      // "evaluates correctly at non-boolean points" case in multilinear.test.ts),
      // so p(x) = sum_{x1 in {0,1}} f(x, x1) = 2 + 4x + 1 = 4x + 3.
      const values = [fe(1), fe(2), fe(3), fe(4)];
      const poly = multilinearExtension(values, F);

      const layerPoly = sumcheckRoundProver(poly, [], 2, F);

      expect(layerPoly.toString()).toBe('4*x + 3');
      expect(layerPoly.degree()).toBe(1);
    });

    test('subsequent rounds substitute the earlier challenges', () => {
      // f = 1 + 2*x0 + x1; after the challenge r0 = 7 the round polynomial is
      // f(7, x) = 1 + 14 + x = x + 15.
      const values = [fe(1), fe(2), fe(3), fe(4)];
      const poly = multilinearExtension(values, F);

      const layerPoly = sumcheckRoundProver(poly, [fe(7)], 2, F);
      expect(layerPoly.toString()).toBe('x + 15');
    });

    test('throws once every variable has been bound', () => {
      const values = [fe(1), fe(2), fe(3), fe(4)];
      const poly = multilinearExtension(values, F);

      expect(() => sumcheckRoundProver(poly, [fe(7), fe(13)], 2, F)).toThrow(
        'prover: no free variable left'
      );
    });

    test('p(0) + p(1) equals claimed sum in first round', () => {
      const values = [fe(9), fe(2), fe(5), fe(4)];
      const poly = multilinearExtension(values, F);
      const claimedSum = values.reduce((a, b) => a.add(b), fe(0));

      const layerPoly = sumcheckRoundProver(poly, [], 2, F);
      const p0 = layerPoly.evaluate(fe(0));
      const p1 = layerPoly.evaluate(fe(1));
      const computedSum = p0.add(p1);

      expect(computedSum.eq(claimedSum)).toBe(true);
    });
  });

  describe('sumcheckRoundVerifier', () => {
    test('returns challenge and new sum on valid round', () => {
      const values = [fe(9), fe(2), fe(5), fe(4)];
      const poly = multilinearExtension(values, F);
      const claimedSum = values.reduce((a, b) => a.add(b), fe(0));

      const layerPoly = sumcheckRoundProver(poly, [], 2, F);
      const result = sumcheckRoundVerifier(claimedSum, layerPoly, F, fe(7));

      expect(result.challenge.eq(7)).toBe(true);
      expect(result.newSum.eq(layerPoly.evaluate(fe(7)))).toBe(true);
    });

    test('throws on invalid sum', () => {
      const values = [fe(9), fe(2), fe(5), fe(4)];
      const poly = multilinearExtension(values, F);
      const claimedSum = values.reduce((a, b) => a.add(b), fe(0));
      const wrongSum = claimedSum.add(fe(1));

      const layerPoly = sumcheckRoundProver(poly, [], 2, F);

      expect(() => sumcheckRoundVerifier(wrongSum, layerPoly, F)).toThrow();
    });
  });

  describe('createPolyEvaluator', () => {
    test('creates evaluator that works with point arrays', () => {
      const values = [fe(1), fe(2), fe(3), fe(4)];
      const poly = multilinearExtension(values, F);

      const evaluator = createPolyEvaluator(poly);

      expect(evaluator([fe(0), fe(0)]).eq(1)).toBe(true);
      expect(evaluator([fe(0), fe(1)]).eq(2)).toBe(true);
      expect(evaluator([fe(1), fe(0)]).eq(3)).toBe(true);
      expect(evaluator([fe(1), fe(1)]).eq(4)).toBe(true);
    });
  });
});

describe('Sumcheck with larger field', () => {
  // Use a cryptographically-sized prime field
  const Flarge = new FiniteFieldPrime(
    8444461749428370424248824938781546531375899335154063827935233455917409239041n
  );
  const fel = (n: bigint) => Flarge.__call__(n);

  test('works with large field', () => {
    const values = [fel(9n), fel(2n), fel(5n), fel(4n), fel(5n), fel(7n)];
    const result = sumcheckRun(values, Flarge);

    expect(result.valid).toBe(true);
  });

  test('prover-verifier roundtrip with large field', () => {
    const values = [fel(1n), fel(2n), fel(3n), fel(4n)];
    const poly = multilinearExtension(values, Flarge);
    const claimedSum = values.reduce((a, b) => a.add(b), fel(0n));

    // Use fixed challenges for deterministic test
    let challengeIdx = 0;
    const challenges = [fel(123456789n), fel(987654321n)];
    const proof = sumcheckProve(poly, claimedSum, 2, Flarge, () => challenges[challengeIdx++]!);

    const evaluator = createPolyEvaluator(poly);
    const valid = sumcheckVerify(proof, claimedSum, evaluator, Flarge, 2, { degreeCheck: 1 });

    expect(valid).toBe(true);
  });
});

describe('Edge cases', () => {
  test('all zeros', () => {
    // Documented deviation: the blueprint's `len(res.variables()) != 1` check
    // makes a constant round polynomial an error (and in fact crashes on the
    // zero polynomial); we accept it, since the zero function legitimately
    // sums to zero in every round.
    const values = [fe(0), fe(0), fe(0), fe(0)];
    const poly = multilinearExtension(values, F);
    expect(poly.isZero()).toBe(true);

    const round1 = sumcheckRoundProver(poly, [], 2, F);
    expect(round1.isZero()).toBe(true);

    const result = sumcheckRun(values, F);

    expect(result.valid).toBe(true);
    expect(result.challenges.length).toBe(2);
    expect(result.finalSum.eq(0)).toBe(true);
  });

  test('single non-zero value', () => {
    const values = [fe(0), fe(0), fe(42), fe(0)];
    const result = sumcheckRun(values, F);

    expect(result.valid).toBe(true);
  });

  test('values that sum to field modulus (zero)', () => {
    // 50 + 51 = 101 = 0 in F_101
    const values = [fe(50), fe(51), fe(0), fe(0)];
    const result = sumcheckRun(values, F);

    expect(result.valid).toBe(true);
  });
});
