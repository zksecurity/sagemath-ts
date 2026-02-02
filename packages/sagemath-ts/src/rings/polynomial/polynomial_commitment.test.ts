/**
 * Unit tests for polynomial commitment operations (KZG, FRI)
 */
import { describe, expect, test } from 'bun:test';
import { type FiniteFieldElement, FiniteFieldPrime } from '../finite_rings/finite_field_prime.js';
import {
  type EvaluationPoint,
  barycentric_weights,
  batch_quotient,
  check_degree_bound,
  compose_with_linear,
  compute_quotient,
  compute_vanishing_polynomial,
  coset_vanishing,
  degree_bound,
  evaluate_basis,
  fri_fold,
  generate_powers,
  lagrange_interpolation,
  linearization,
  multi_evaluate,
  recombine_chunks,
  roots_of_unity_vanishing,
  split_poly,
  verify_quotient_proof,
} from './polynomial_commitment.js';
import { PolynomialRingConstructor } from './polynomial_ring.js';

// Use a prime field for testing (p = 101 for small examples)
const F101 = new FiniteFieldPrime(101n);
const [R, x] = PolynomialRingConstructor(F101, 'x');

// Helper to create field elements
const fe = (n: number | bigint) => F101.__call__(n);

// Helper to create polynomials from coefficient arrays
const poly = (coeffs: number[]) => R.__call__(coeffs.map(fe));

describe('Quotient polynomial operations (KZG)', () => {
  describe('compute_quotient', () => {
    test('basic quotient: (x^2 - 1) / (x - 1) = x + 1', () => {
      // f(X) = X^2 - 1
      const f = x.pow(2).sub(R.one());
      // z = 1, y = f(1) = 0
      const z = fe(1);
      const y = f.evaluate(z);
      expect(y.isZero()).toBe(true);

      // q(X) should be X + 1
      const q = compute_quotient(f, z, y);
      expect(q.degree()).toBe(1);
      expect(q.getCoeff(0).eq(1)).toBe(true);
      expect(q.getCoeff(1).eq(1)).toBe(true);
    });

    test('quotient with non-zero evaluation', () => {
      // f(X) = X^2 + 2X + 3
      const f = poly([3, 2, 1]);
      // z = 2, y = f(2) = 4 + 4 + 3 = 11
      const z = fe(2);
      const y = f.evaluate(z);
      expect(y.eq(11)).toBe(true);

      // q(X) = (f(X) - 11) / (X - 2)
      const q = compute_quotient(f, z, y);

      // Verify: q(X) * (X - 2) + 11 = f(X)
      const divisor = x.sub(R.__call__(z));
      const reconstructed = q.mul(divisor).add(R.__call__(y));
      expect(reconstructed.eq(f)).toBe(true);
    });

    test('quotient degree is one less than original', () => {
      // f(X) = X^5 + X^3 + X + 1
      const f = poly([1, 1, 0, 1, 0, 1]);
      const z = fe(3);
      const y = f.evaluate(z);

      const q = compute_quotient(f, z, y);
      expect(q.degree()).toBe(f.degree() - 1);
    });

    test('throws when f(z) != y', () => {
      const f = x.pow(2).add(R.one()); // f(X) = X^2 + 1
      const z = fe(2);
      const wrong_y = fe(10); // f(2) = 5, not 10

      expect(() => compute_quotient(f, z, wrong_y)).toThrow();
    });

    test('constant polynomial quotient', () => {
      // f(X) = 5 (constant)
      const f = R.__call__(fe(5));
      const z = fe(7);
      const y = f.evaluate(z);
      expect(y.eq(5)).toBe(true);

      // q(X) should be 0 since f(X) - y = 0
      const q = compute_quotient(f, z, y);
      expect(q.isZero()).toBe(true);
    });
  });

  describe('batch_quotient', () => {
    test('batch quotient with two points', () => {
      // f(X) that passes through (1, 3) and (2, 7)
      // First create interpolation: the line through these points
      // Slope = (7-3)/(2-1) = 4, so y = 4x - 1
      // f(X) = 4X - 1 = X^2 + 4X - 1 would have quotient when divided by (X-1)(X-2)

      // Actually let's construct f = I(X) + Z_H(X) * q(X) for known q
      // I(X) through (1, 3) and (2, 7): slope = 4, I(X) = 4X - 1
      const I = poly([100, 4]); // -1 mod 101 = 100
      // Z_H(X) = (X - 1)(X - 2) = X^2 - 3X + 2
      const ZH = poly([2, 98, 1]); // -3 mod 101 = 98
      // Let q(X) = X + 1
      const q = poly([1, 1]);
      // f(X) = I(X) + Z_H(X) * q(X)
      const f = I.add(ZH.mul(q));

      const points: EvaluationPoint<FiniteFieldElement>[] = [
        { x: fe(1), y: fe(3) },
        { x: fe(2), y: fe(7) },
      ];

      const computed_q = batch_quotient(f, points);
      expect(computed_q.eq(q)).toBe(true);
    });

    test('batch quotient matches individual quotients when applicable', () => {
      // For a single point, batch_quotient should match compute_quotient
      const f = poly([1, 2, 3, 4]); // f(X) = 1 + 2X + 3X^2 + 4X^3
      const z = fe(5);
      const y = f.evaluate(z);

      const single_q = compute_quotient(f, z, y);
      const batch_q = batch_quotient(f, [{ x: z, y }]);

      expect(batch_q.eq(single_q)).toBe(true);
    });

    test('throws on empty points array', () => {
      const f = x.pow(2);
      expect(() => batch_quotient(f, [])).toThrow();
    });

    test('throws when evaluations do not match', () => {
      const f = x.pow(2).add(x).add(R.one());
      const points: EvaluationPoint<FiniteFieldElement>[] = [
        { x: fe(1), y: fe(99) }, // Wrong value
        { x: fe(2), y: fe(99) }, // Wrong value
      ];

      expect(() => batch_quotient(f, points)).toThrow();
    });
  });

  describe('verify_quotient_proof', () => {
    test('valid proof passes verification', () => {
      const f = poly([1, 2, 3]); // f(X) = 1 + 2X + 3X^2
      const z = fe(4);
      const y = f.evaluate(z);
      const q = compute_quotient(f, z, y);

      expect(verify_quotient_proof(f, z, y, q)).toBe(true);
    });

    test('invalid quotient fails verification', () => {
      const f = poly([1, 2, 3]);
      const z = fe(4);
      const y = f.evaluate(z);
      const wrong_q = poly([1, 1]); // Wrong quotient

      expect(verify_quotient_proof(f, z, y, wrong_q)).toBe(false);
    });

    test('wrong evaluation point fails verification', () => {
      const f = poly([1, 2, 3]);
      const z = fe(4);
      const y = f.evaluate(z);
      const q = compute_quotient(f, z, y);

      // Try to verify with different z
      const wrong_z = fe(5);
      expect(verify_quotient_proof(f, wrong_z, y, q)).toBe(false);
    });
  });
});

describe('Vanishing polynomial and interpolation', () => {
  describe('compute_vanishing_polynomial', () => {
    test('vanishing polynomial on two points', () => {
      const domain = [fe(1), fe(2)];
      const ZH = compute_vanishing_polynomial(R, domain);

      // Z_H(X) = (X - 1)(X - 2) = X^2 - 3X + 2
      expect(ZH.degree()).toBe(2);
      expect(ZH.evaluate(fe(1)).isZero()).toBe(true);
      expect(ZH.evaluate(fe(2)).isZero()).toBe(true);
      expect(ZH.evaluate(fe(3)).isZero()).toBe(false);
    });

    test('vanishing polynomial on empty domain is 1', () => {
      const ZH = compute_vanishing_polynomial(R, []);
      expect(ZH.eq(R.one())).toBe(true);
    });

    test('vanishing polynomial has correct degree', () => {
      const domain = [fe(1), fe(3), fe(5), fe(7), fe(9)];
      const ZH = compute_vanishing_polynomial(R, domain);
      expect(ZH.degree()).toBe(5);
    });
  });

  describe('lagrange_interpolation', () => {
    test('interpolation through two points', () => {
      // Points: (1, 5) and (3, 11)
      // Line: slope = (11-5)/(3-1) = 3, y = 3x + 2
      const points: EvaluationPoint<FiniteFieldElement>[] = [
        { x: fe(1), y: fe(5) },
        { x: fe(3), y: fe(11) },
      ];

      const I = lagrange_interpolation(R, points);

      expect(I.degree()).toBe(1);
      expect(I.evaluate(fe(1)).eq(5)).toBe(true);
      expect(I.evaluate(fe(3)).eq(11)).toBe(true);
    });

    test('interpolation through three points', () => {
      // Parabola through (0, 1), (1, 0), (2, 3)
      const points: EvaluationPoint<FiniteFieldElement>[] = [
        { x: fe(0), y: fe(1) },
        { x: fe(1), y: fe(0) },
        { x: fe(2), y: fe(3) },
      ];

      const I = lagrange_interpolation(R, points);

      expect(I.degree()).toBe(2);
      expect(I.evaluate(fe(0)).eq(1)).toBe(true);
      expect(I.evaluate(fe(1)).eq(0)).toBe(true);
      expect(I.evaluate(fe(2)).eq(3)).toBe(true);
    });

    test('interpolation single point is constant', () => {
      const points: EvaluationPoint<FiniteFieldElement>[] = [{ x: fe(5), y: fe(7) }];
      const I = lagrange_interpolation(R, points);

      expect(I.isConstant()).toBe(true);
      expect(I.evaluate(fe(5)).eq(7)).toBe(true);
      expect(I.evaluate(fe(0)).eq(7)).toBe(true);
    });

    test('interpolation empty points is zero', () => {
      const I = lagrange_interpolation(R, []);
      expect(I.isZero()).toBe(true);
    });
  });

  describe('roots_of_unity_vanishing', () => {
    test('X^n - 1 for n = 4', () => {
      const ZH = roots_of_unity_vanishing(R, 4);
      expect(ZH.degree()).toBe(4);
      expect(ZH.getCoeff(0).eq(100)).toBe(true); // -1 mod 101
      expect(ZH.getCoeff(4).eq(1)).toBe(true);

      // Check that 1 is a root (1^4 - 1 = 0)
      expect(ZH.evaluate(fe(1)).isZero()).toBe(true);
    });

    test('throws for non-positive n', () => {
      expect(() => roots_of_unity_vanishing(R, 0)).toThrow();
      expect(() => roots_of_unity_vanishing(R, -1)).toThrow();
    });
  });

  describe('coset_vanishing', () => {
    test('coset vanishing X^n - k^n', () => {
      const k = fe(3);
      const ZH = coset_vanishing(R, 4, k);

      // Should be X^4 - 81
      expect(ZH.degree()).toBe(4);
      expect(ZH.getCoeff(4).eq(1)).toBe(true);
      // 3^4 = 81, so constant term is -81 mod 101 = 20
      expect(ZH.getCoeff(0).eq(20)).toBe(true);
    });
  });
});

describe('Linearization', () => {
  describe('linearization', () => {
    test('linear combination of polynomials', () => {
      const f0 = poly([1, 0, 1]); // 1 + X^2
      const f1 = poly([0, 1]); // X
      const f2 = poly([2]); // 2

      const alpha = fe(3);
      const alpha2 = fe(9);
      const alpha3 = fe(27);

      const L = linearization([f0, f1, f2], [alpha, alpha2, alpha3]);

      // L = 3 * (1 + X^2) + 9 * X + 27 * 2
      // = 3 + 3X^2 + 9X + 54
      // = 57 + 9X + 3X^2
      expect(L.getCoeff(0).eq(57)).toBe(true);
      expect(L.getCoeff(1).eq(9)).toBe(true);
      expect(L.getCoeff(2).eq(3)).toBe(true);
    });

    test('throws on mismatched lengths', () => {
      const f0 = poly([1]);
      const f1 = poly([2]);

      expect(() => linearization([f0, f1], [fe(1)])).toThrow();
    });

    test('throws on empty arrays', () => {
      expect(() => linearization([], [])).toThrow();
    });
  });

  describe('generate_powers', () => {
    test('generates powers correctly', () => {
      const alpha = fe(2);
      const powers = generate_powers(alpha, 5);

      expect(powers.length).toBe(5);
      expect(powers[0]!.eq(1)).toBe(true); // 2^0 = 1
      expect(powers[1]!.eq(2)).toBe(true); // 2^1 = 2
      expect(powers[2]!.eq(4)).toBe(true); // 2^2 = 4
      expect(powers[3]!.eq(8)).toBe(true); // 2^3 = 8
      expect(powers[4]!.eq(16)).toBe(true); // 2^4 = 16
    });

    test('empty for n = 0', () => {
      const powers = generate_powers(fe(5), 0);
      expect(powers.length).toBe(0);
    });
  });
});

describe('Lagrange basis evaluation', () => {
  describe('barycentric_weights', () => {
    test('weights for two-point domain', () => {
      const domain = [fe(1), fe(3)];
      const weights = barycentric_weights(domain);

      expect(weights.length).toBe(2);
      // w_0 = 1 / (1 - 3) = 1 / (-2) = -1/2 mod 101
      // w_1 = 1 / (3 - 1) = 1 / 2 = 51 mod 101
      expect(weights[1]!.eq(51)).toBe(true);
    });

    test('throws on duplicate points', () => {
      expect(() => barycentric_weights([fe(1), fe(1)])).toThrow();
    });
  });

  describe('evaluate_basis', () => {
    test('basis evaluates to 1 at domain point', () => {
      const domain = [fe(1), fe(2), fe(3)];
      const basis_at_1 = evaluate_basis(domain, fe(1));

      expect(basis_at_1.length).toBe(3);
      expect(basis_at_1[0]!.eq(1)).toBe(true);
      expect(basis_at_1[1]!.eq(0)).toBe(true);
      expect(basis_at_1[2]!.eq(0)).toBe(true);
    });

    test('basis evaluates correctly at non-domain point', () => {
      const domain = [fe(0), fe(1)];
      const point = fe(2);
      const basis = evaluate_basis(domain, point);

      // L_0(X) = (X - 1) / (0 - 1) = -(X - 1) = 1 - X
      // L_1(X) = (X - 0) / (1 - 0) = X
      // L_0(2) = 1 - 2 = -1 = 100 mod 101
      // L_1(2) = 2
      expect(basis[0]!.eq(100)).toBe(true);
      expect(basis[1]!.eq(2)).toBe(true);
    });

    test('sum of basis evaluations is 1', () => {
      const domain = [fe(5), fe(10), fe(15), fe(20)];
      const point = fe(7);
      const basis = evaluate_basis(domain, point);

      let sum = fe(0);
      for (const b of basis) {
        sum = sum.add(b);
      }
      expect(sum.eq(1)).toBe(true);
    });
  });
});

describe('Degree operations', () => {
  describe('degree_bound', () => {
    test('returns correct degree', () => {
      expect(degree_bound(poly([1, 2, 3]))).toBe(2);
      expect(degree_bound(poly([5]))).toBe(0);
      expect(degree_bound(R.zero())).toBe(-1);
    });
  });

  describe('check_degree_bound', () => {
    test('returns true when within bound', () => {
      const f = poly([1, 2, 3]); // degree 2
      expect(check_degree_bound(f, 3)).toBe(true);
      expect(check_degree_bound(f, 4)).toBe(true);
    });

    test('returns false when exceeding bound', () => {
      const f = poly([1, 2, 3]); // degree 2
      expect(check_degree_bound(f, 2)).toBe(false);
      expect(check_degree_bound(f, 1)).toBe(false);
    });
  });
});

describe('Polynomial splitting (multi-part commitments)', () => {
  describe('split_poly', () => {
    test('splits polynomial into chunks', () => {
      // f(X) = a0 + a1*X + a2*X^2 + a3*X^3 + a4*X^4 + a5*X^5
      const f = poly([1, 2, 3, 4, 5, 6]);
      const chunks = split_poly(f, 2);

      expect(chunks.length).toBe(3);
      // chunks[0] = 1 + 2X
      expect(chunks[0]!.getCoeff(0).eq(1)).toBe(true);
      expect(chunks[0]!.getCoeff(1).eq(2)).toBe(true);
      // chunks[1] = 3 + 4X
      expect(chunks[1]!.getCoeff(0).eq(3)).toBe(true);
      expect(chunks[1]!.getCoeff(1).eq(4)).toBe(true);
      // chunks[2] = 5 + 6X
      expect(chunks[2]!.getCoeff(0).eq(5)).toBe(true);
      expect(chunks[2]!.getCoeff(1).eq(6)).toBe(true);
    });

    test('handles uneven splits', () => {
      const f = poly([1, 2, 3, 4, 5]); // 5 coeffs
      const chunks = split_poly(f, 2);

      expect(chunks.length).toBe(3);
      expect(chunks[2]!.degree()).toBe(0); // Only one coefficient
    });

    test('split of zero polynomial', () => {
      const chunks = split_poly(R.zero(), 3);
      expect(chunks.length).toBe(1);
      expect(chunks[0]!.isZero()).toBe(true);
    });

    test('throws on non-positive chunk size', () => {
      expect(() => split_poly(x.pow(2), 0)).toThrow();
      expect(() => split_poly(x.pow(2), -1)).toThrow();
    });
  });

  describe('recombine_chunks', () => {
    test('recombines chunks correctly', () => {
      const f = poly([1, 2, 3, 4, 5, 6]);
      const chunk_size = 2;
      const chunks = split_poly(f, chunk_size);
      const recombined = recombine_chunks(chunks, chunk_size);

      expect(recombined.eq(f)).toBe(true);
    });

    test('round-trip preserves polynomial', () => {
      // Test with various chunk sizes
      for (const chunk_size of [1, 2, 3, 4, 5]) {
        const f = poly([7, 11, 13, 17, 19, 23, 29]);
        const chunks = split_poly(f, chunk_size);
        const recombined = recombine_chunks(chunks, chunk_size);
        expect(recombined.eq(f)).toBe(true);
      }
    });

    test('throws on empty chunks', () => {
      expect(() => recombine_chunks([], 2)).toThrow();
    });
  });
});

describe('FRI operations', () => {
  describe('compose_with_linear', () => {
    test('f(aX + b) for linear f', () => {
      // f(X) = 2X + 3
      const f = poly([3, 2]);
      // f(5X + 7) = 2(5X + 7) + 3 = 10X + 17
      const result = compose_with_linear(f, fe(5), fe(7));

      expect(result.degree()).toBe(1);
      expect(result.getCoeff(0).eq(17)).toBe(true);
      expect(result.getCoeff(1).eq(10)).toBe(true);
    });

    test('f(aX + b) for quadratic f', () => {
      // f(X) = X^2
      const f = poly([0, 0, 1]);
      // f(2X + 1) = (2X + 1)^2 = 4X^2 + 4X + 1
      const result = compose_with_linear(f, fe(2), fe(1));

      expect(result.degree()).toBe(2);
      expect(result.getCoeff(0).eq(1)).toBe(true);
      expect(result.getCoeff(1).eq(4)).toBe(true);
      expect(result.getCoeff(2).eq(4)).toBe(true);
    });

    test('composition preserves degree', () => {
      const f = poly([1, 2, 3, 4, 5]); // degree 4
      const result = compose_with_linear(f, fe(7), fe(11));
      expect(result.degree()).toBe(4);
    });

    test('zero polynomial stays zero', () => {
      expect(compose_with_linear(R.zero(), fe(2), fe(3)).isZero()).toBe(true);
    });
  });

  describe('fri_fold', () => {
    test('basic FRI folding', () => {
      // f(X) = a0 + a1*X + a2*X^2 + a3*X^3
      // f_even(Y) = a0 + a2*Y
      // f_odd(Y) = a1 + a3*Y
      // fri_fold(f, alpha) = f_even(Y) + alpha * f_odd(Y)
      //                    = (a0 + alpha*a1) + (a2 + alpha*a3)*Y
      const f = poly([1, 2, 3, 4]); // a0=1, a1=2, a2=3, a3=4
      const alpha = fe(5);

      const folded = fri_fold(f, alpha);

      // Expected: (1 + 5*2) + (3 + 5*4)*Y = 11 + 23*Y
      expect(folded.degree()).toBe(1);
      expect(folded.getCoeff(0).eq(11)).toBe(true);
      expect(folded.getCoeff(1).eq(23)).toBe(true);
    });

    test('FRI folding reduces degree by half', () => {
      const f = poly([1, 2, 3, 4, 5, 6, 7, 8]); // degree 7
      const folded = fri_fold(f, fe(3));
      expect(folded.degree()).toBe(3); // (7 + 1) / 2 - 1 = 3
    });

    test('FRI folding odd-degree polynomial', () => {
      const f = poly([1, 2, 3, 4, 5]); // degree 4, 5 coeffs
      const folded = fri_fold(f, fe(2));
      // Even: [1, 3, 5], Odd: [2, 4]
      // folded = (1 + 2*2) + (3 + 2*4)*Y + 5*Y^2 = 5 + 11*Y + 5*Y^2
      expect(folded.degree()).toBe(2);
    });

    test('zero polynomial stays zero', () => {
      expect(fri_fold(R.zero(), fe(10)).isZero()).toBe(true);
    });
  });
});

describe('Multi-evaluation', () => {
  describe('multi_evaluate', () => {
    test('evaluates at multiple points', () => {
      const f = poly([1, 2, 1]); // 1 + 2X + X^2 = (X + 1)^2
      const points = [fe(0), fe(1), fe(2), fe(99)]; // 99 = -2 mod 101

      const evals = multi_evaluate(f, points);

      expect(evals.length).toBe(4);
      expect(evals[0]!.eq(1)).toBe(true); // f(0) = 1
      expect(evals[1]!.eq(4)).toBe(true); // f(1) = 4
      expect(evals[2]!.eq(9)).toBe(true); // f(2) = 9
      expect(evals[3]!.eq(1)).toBe(true); // f(-2) = 1 - 4 + 4 = 1
    });

    test('empty points returns empty array', () => {
      const evals = multi_evaluate(x.pow(2), []);
      expect(evals.length).toBe(0);
    });
  });
});

describe('Integration tests', () => {
  test('KZG opening proof flow', () => {
    // 1. Commit to polynomial (simulated)
    const f = poly([5, 3, 1, 7, 2]); // Some polynomial

    // 2. Prover receives challenge point z
    const z = fe(13);

    // 3. Prover computes evaluation y = f(z)
    const y = f.evaluate(z);

    // 4. Prover computes quotient polynomial (the proof)
    const q = compute_quotient(f, z, y);

    // 5. Verify proof
    expect(verify_quotient_proof(f, z, y, q)).toBe(true);

    // 6. Verify degree bound
    expect(q.degree()).toBe(f.degree() - 1);
  });

  test('batch opening proof flow', () => {
    // Polynomial to commit to
    const f = poly([1, 2, 3, 4, 5]);

    // Multiple evaluation points
    const z0 = fe(7);
    const z1 = fe(11);
    const z2 = fe(13);

    // Compute evaluations
    const points: EvaluationPoint<FiniteFieldElement>[] = [
      { x: z0, y: f.evaluate(z0) },
      { x: z1, y: f.evaluate(z1) },
      { x: z2, y: f.evaluate(z2) },
    ];

    // Compute batch quotient
    const q = batch_quotient(f, points);

    // Verify: f(X) = I(X) + Z_H(X) * q(X)
    const I = lagrange_interpolation(R, points);
    const ZH = compute_vanishing_polynomial(
      R,
      points.map((p) => p.x)
    );
    const reconstructed = I.add(ZH.mul(q));

    expect(reconstructed.eq(f)).toBe(true);
  });

  test('FRI-like protocol simulation', () => {
    // Start with a polynomial of degree 7
    let f = poly([1, 2, 3, 4, 5, 6, 7, 8]);
    expect(f.degree()).toBe(7);

    // Simulate FRI folding rounds
    const challenges = [fe(3), fe(7), fe(11)];

    for (const alpha of challenges) {
      f = fri_fold(f, alpha);
    }

    // After 3 rounds of folding, degree should be 0 (constant)
    expect(f.degree()).toBe(0);
  });

  test('polynomial commitment with splitting', () => {
    // Large polynomial
    const f = poly([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);

    // Split into chunks of size 4
    const chunks = split_poly(f, 4);
    expect(chunks.length).toBe(3);

    // Each chunk can be committed separately
    for (const chunk of chunks) {
      expect(chunk.degree()).toBeLessThan(4);
    }

    // Recombine to verify
    const recombined = recombine_chunks(chunks, 4);
    expect(recombined.eq(f)).toBe(true);
  });

  test('linearization in multi-polynomial protocol', () => {
    // Multiple polynomials in a protocol
    const f_A = poly([1, 0, 1]); // A(X) = 1 + X^2
    const f_B = poly([0, 1, 0]); // B(X) = X
    const f_C = poly([1, 1, 1]); // C(X) = 1 + X + X^2

    // Random challenge
    const alpha = fe(5);
    const powers = generate_powers(alpha, 3);

    // Compute linearization
    const L = linearization([f_A, f_B, f_C], powers);

    // Verify at a random point
    const z = fe(7);
    const L_z = L.evaluate(z);

    // Should equal: 1*A(7) + 5*B(7) + 25*C(7)
    const expected = f_A
      .evaluate(z)
      .add(f_B.evaluate(z).mul(powers[1]!))
      .add(f_C.evaluate(z).mul(powers[2]!));

    expect(L_z.eq(expected)).toBe(true);
  });
});

describe('Edge cases and error handling', () => {
  test('operations with zero polynomial', () => {
    const zero = R.zero();
    const z = fe(5);
    const y = fe(0);

    // Zero polynomial evaluates to 0 everywhere
    expect(zero.evaluate(z).isZero()).toBe(true);

    // Quotient of zero is zero
    const q = compute_quotient(zero, z, y);
    expect(q.isZero()).toBe(true);
  });

  test('operations with constant polynomial', () => {
    const c = R.__call__(fe(42));
    const z = fe(7);
    const y = c.evaluate(z);

    expect(y.eq(42)).toBe(true);

    const q = compute_quotient(c, z, y);
    expect(q.isZero()).toBe(true);
  });

  test('large polynomial operations', () => {
    // Create a larger polynomial
    const coeffs: number[] = [];
    for (let i = 0; i < 100; i++) {
      coeffs.push(i + 1);
    }
    const f = poly(coeffs);

    const z = fe(17);
    const y = f.evaluate(z);

    // Should not throw
    const q = compute_quotient(f, z, y);
    expect(verify_quotient_proof(f, z, y, q)).toBe(true);
  });
});
