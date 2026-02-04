/**
 * Unit tests for the sumcheck protocol
 */
import { describe, expect, test } from 'bun:test';
import { FiniteFieldPrime } from '../rings/finite_rings/finite_field_prime.js';
import { MPolynomialRingConstructor } from '../rings/polynomial/multi_polynomial_ring.js';
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
      expect(binaryToInt([0, 0, 0, 0])).toBe(0);
      expect(binaryToInt([0, 1, 0, 1])).toBe(5);
      expect(binaryToInt([1, 1, 1])).toBe(7);
      expect(binaryToInt([])).toBe(0);
    });

    test('roundtrip with intToBinary', () => {
      for (let i = 0; i < 16; i++) {
        expect(binaryToInt(intToBinary(i, 4))).toBe(i);
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

    test('passes for single value', () => {
      const values = [fe(42)];
      const result = sumcheckRun(values, F);

      expect(result.valid).toBe(true);
      expect(result.challenges.length).toBe(0); // 0 variables (constant)
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
      const valid = sumcheckVerify(proof, claimedSum, evaluator, F);

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
      const valid = sumcheckVerify(proof, wrongSum, evaluator, F);

      expect(valid).toBe(false);
    });

    test('verifier rejects proof with tampered round polynomial', () => {
      const values = [fe(9), fe(2), fe(5), fe(4)];
      const poly = multilinearExtension(values, F);
      const claimedSum = values.reduce((a, b) => a.add(b), fe(0));

      let challengeIdx = 0;
      const challenges = [fe(7), fe(13)];
      const proof = sumcheckProve(poly, claimedSum, 2, F, () => challenges[challengeIdx++]!);

      // Tamper with the first round polynomial by modifying the final evaluation
      const tamperedProof = {
        ...proof,
        finalEvaluation: proof.finalEvaluation.add(fe(1)),
      };

      const evaluator = createPolyEvaluator(poly);
      const valid = sumcheckVerify(tamperedProof, claimedSum, evaluator, F);

      expect(valid).toBe(false);
    });
  });

  describe('sumcheckRoundProver', () => {
    test('produces degree-1 polynomial for MLEs', () => {
      const values = [fe(1), fe(2), fe(3), fe(4)];
      const poly = multilinearExtension(values, F);

      const layerPoly = sumcheckRoundProver(poly, [], 2, F);

      expect(layerPoly.degree()).toBeLessThanOrEqual(1);
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
    const valid = sumcheckVerify(proof, claimedSum, evaluator, Flarge);

    expect(valid).toBe(true);
  });
});

describe('Edge cases', () => {
  test('all zeros', () => {
    const values = [fe(0), fe(0), fe(0), fe(0)];
    const result = sumcheckRun(values, F);

    expect(result.valid).toBe(true);
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
