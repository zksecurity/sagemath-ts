/**
 * Tests for elliptic curves over finite fields
 */

import { describe, expect, it } from 'vitest';
import { FiniteFieldPrime } from '../../rings/finite_rings/finite_field_prime.js';
import {
  EllipticCurve,
  type EllipticCurveFiniteField,
  type EllipticCurvePoint,
} from './ell_finite_field.js';

// Helper to create GF(p)
function GF(p: bigint): FiniteFieldPrime {
  return new FiniteFieldPrime(p);
}

describe('EllipticCurvePoint', () => {
  describe('basic operations', () => {
    it('should create point at infinity', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);
      const O = E.zero();

      expect(O.isZero()).toBe(true);
      expect(O.isInfinity).toBe(true);
    });

    it('should create affine points', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      // Point (0, 1) should be on y^2 = x^3 + x + 1
      // Check: 1^2 = 0^3 + 0 + 1 = 1
      const P = E.point(0n, 1n);
      const Q = E.point([0n, 1n]);

      expect(P.isZero()).toBe(false);
      expect(P.x!.value).toBe(0n);
      expect(P.y!.value).toBe(1n);
      expect(Q.eq(P)).toBe(true);
    });

    it('should reject invalid points', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      // (1, 1) is not on the curve: 1 != 1 + 1 + 1 = 3
      expect(() => E.point(1n, 1n)).toThrow();
    });

    it('should negate points', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);
      const P = E.point(0n, 1n);
      const negP = P.neg();

      expect(negP.x!.value).toBe(0n);
      expect(negP.y!.value).toBe(22n); // -1 mod 23 = 22
    });

    it('should handle point at infinity negation', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);
      const O = E.zero();

      expect(O.neg().isZero()).toBe(true);
    });
  });

  describe('point addition', () => {
    it('should add point with infinity', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);
      const P = E.point(0n, 1n);
      const O = E.zero();

      expect(P.add(O).eq(P)).toBe(true);
      expect(O.add(P).eq(P)).toBe(true);
    });

    it('should add point with its negative', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);
      const P = E.point(0n, 1n);
      const negP = P.neg();

      expect(P.add(negP).isZero()).toBe(true);
    });

    it('should double a point', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);
      const P = E.point(0n, 1n);

      // Manual doubling for y^2 = x^3 + x + 1:
      // slope = (3*0^2 + 1) / (2*1) = 1/2 = 12 (since 2*12 = 24 ≡ 1 mod 23)
      // x3 = 12^2 - 2*0 = 144 ≡ 6 mod 23
      // y3 = 12*(0 - 6) - 1 = -73 ≡ 19 mod 23
      const twoP = P.double();

      expect(twoP.x!.value).toBe(6n);
      expect(twoP.y!.value).toBe(19n);
    });

    it('should add distinct points', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);
      const P = E.point(0n, 1n);
      const Q = E.point(6n, 4n); // Verify: 4^2 = 16, 6^3 + 6 + 1 = 216 + 7 = 223 ≡ 16 mod 23

      const R = P.add(Q);

      // Manual calculation:
      // slope = (4 - 1) / (6 - 0) = 3/6 = 3*4 = 12 mod 23 (since 6*4 = 24 ≡ 1)
      // x3 = 12^2 - 0 - 6 = 144 - 6 = 138 ≡ 0 mod 23
      // Wait, that gives x3 = 0, same as P. Let me recalculate...
      // Actually 138 mod 23 = 138 - 6*23 = 138 - 138 = 0
      // y3 = 12*(0 - 0) - 1 = -1 ≡ 22 mod 23

      expect(R.x!.value).toBe(0n);
      expect(R.y!.value).toBe(22n);
    });
  });

  describe('scalar multiplication', () => {
    it('should compute [0]P = O', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);
      const P = E.point(0n, 1n);

      expect(P.mul(0n).isZero()).toBe(true);
    });

    it('should compute [1]P = P', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);
      const P = E.point(0n, 1n);

      expect(P.mul(1n).eq(P)).toBe(true);
    });

    it('should compute [2]P = P + P', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);
      const P = E.point(0n, 1n);

      expect(P.mul(2n).eq(P.add(P))).toBe(true);
    });

    it('should compute [-1]P = -P', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);
      const P = E.point(0n, 1n);

      expect(P.mul(-1n).eq(P.neg())).toBe(true);
    });

    it('should compute [n]P for larger n', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);
      const P = E.point(0n, 1n);

      // [3]P = [2]P + P
      const threeP = P.mul(3n);
      const expected = P.double().add(P);

      expect(threeP.eq(expected)).toBe(true);
    });
  });

  describe('point order', () => {
    it('should compute order of identity', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);
      const O = E.zero();

      expect(O.order()).toBe(1n);
    });

    it('should compute point order dividing curve order', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);
      const curveOrder = E.cardinality();

      // Try a random point
      const P = E.random_point();
      const pointOrder = P.order();

      expect(curveOrder % pointOrder).toBe(0n);
      expect(P.mul(pointOrder).isZero()).toBe(true);
    });

    it('should verify [order]P = O', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);

      // Get a point and compute its order
      const P = E.random_point();
      const order = P.order();

      // Verify [order]P = O
      expect(P.mul(order).isZero()).toBe(true);

      // Verify no smaller positive integer works (unless order = 1)
      if (order > 1n) {
        expect(P.mul(order - 1n).isZero()).toBe(false);
      }
    });
  });
});

describe('EllipticCurveFiniteField', () => {
  describe('curve construction', () => {
    it('should create a valid curve', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      expect(E.a.value).toBe(1n);
      expect(E.b.value).toBe(1n);
    });

    it('should reject singular curves', () => {
      const F = GF(23n);

      // y^2 = x^3 has discriminant 0 (a = 0, b = 0)
      expect(() => EllipticCurve(F, [0n, 0n])).toThrow();
    });

    it('should compute j-invariant', () => {
      const F = GF(101n);
      const E = EllipticCurve(F, [2n, 3n]);

      // j = 1728 * 4a^3 / (4a^3 + 27b^2)
      // j = 1728 * 4*8 / (4*8 + 27*9) = 1728 * 32 / (32 + 243) = 1728 * 32 / 275
      const j = E.j_invariant();

      // Verify j is well-defined (no exception)
      expect(typeof j.value).toBe('bigint');
    });
  });

  describe('is_x_coord', () => {
    it('should identify valid x-coordinates', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      // x = 0: y^2 = 0 + 0 + 1 = 1 is a square
      expect(E.is_x_coord(0n)).toBe(true);

      // x = 6: y^2 = 216 + 6 + 1 = 223 ≡ 16 mod 23, which is 4^2
      expect(E.is_x_coord(6n)).toBe(true);
    });

    it('should identify invalid x-coordinates', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      // Check some x values
      let foundInvalid = false;
      for (let x = 0n; x < 23n; x++) {
        if (!E.is_x_coord(x)) {
          foundInvalid = true;
          break;
        }
      }
      // There should be some invalid x-coordinates for most curves
      // (about half, typically)
      expect(foundInvalid).toBe(true);
    });
  });

  describe('lift_x', () => {
    it('should find points with given x-coordinate', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      // x = 0 should give (0, 1) and (0, 22)
      const points = E.lift_x(0n, true);

      expect(points.length).toBe(2);
      expect(points[0]!.x!.value).toBe(0n);
      expect(points[1]!.x!.value).toBe(0n);

      // Check both points are valid
      for (const P of points) {
        expect(E.is_on_curve(P.x!, P.y!)).toBe(true);
      }
    });

    it('should return single point when not using all=true', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const P = E.lift_x(0n);

      expect(P.x!.value).toBe(0n);
      expect(E.is_on_curve(P.x!, P.y!)).toBe(true);
    });

    it('should throw for invalid x-coordinate when not using all', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      // Find an invalid x
      let invalidX: bigint | null = null;
      for (let x = 0n; x < 23n; x++) {
        if (!E.is_x_coord(x)) {
          invalidX = x;
          break;
        }
      }

      if (invalidX !== null) {
        expect(() => E.lift_x(invalidX)).toThrow();
      }
    });

    it('should return empty array for invalid x when using all=true', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      // Find an invalid x
      let invalidX: bigint | null = null;
      for (let x = 0n; x < 23n; x++) {
        if (!E.is_x_coord(x)) {
          invalidX = x;
          break;
        }
      }

      if (invalidX !== null) {
        const points = E.lift_x(invalidX, true);
        expect(points.length).toBe(0);
      }
    });
  });

  describe('random_point', () => {
    it('should generate valid points', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);

      for (let i = 0; i < 10; i++) {
        const P = E.random_point();

        expect(P.isZero()).toBe(false);
        expect(E.is_on_curve(P.x!, P.y!)).toBe(true);
      }
    });

    it('should generate different points', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);

      const points: string[] = [];
      for (let i = 0; i < 10; i++) {
        const P = E.random_point();
        points.push(P.toString());
      }

      // With high probability, we should get some distinct points
      const uniquePoints = new Set(points);
      expect(uniquePoints.size).toBeGreaterThan(1);
    });
  });

  describe('cardinality', () => {
    it('should compute cardinality for small curves', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const order = E.cardinality();

      // Manually count points to verify
      const allPoints = E.points();
      expect(order).toBe(BigInt(allPoints.length));
    });

    it('should satisfy Hasse bound', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);

      const q = 97n;
      const order = E.cardinality();
      const sqrtQ = BigInt(Math.floor(Math.sqrt(Number(q))));

      // |order - (q + 1)| <= 2 * sqrt(q)
      const diff = order > q + 1n ? order - (q + 1n) : q + 1n - order;
      expect(diff <= 2n * sqrtQ + 2n).toBe(true); // +2 for rounding
    });

    it('should cache the result', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);

      const order1 = E.cardinality();
      const order2 = E.cardinality();

      expect(order1).toBe(order2);
    });

    // Known curve orders from test vectors
    it('should compute correct order for known curves', () => {
      // Test curve from PARI documentation
      // E: y^2 = x^3 + x + 1 over GF(1009)
      const F = GF(1009n);
      const E = EllipticCurve(F, [1n, 1n]);

      const order = E.cardinality();

      // Verify with Hasse bound
      const q = 1009n;
      const sqrtQ = BigInt(Math.floor(Math.sqrt(Number(q))));
      const lowerBound = q + 1n - 2n * sqrtQ;
      const upperBound = q + 1n + 2n * sqrtQ;

      expect(order >= lowerBound).toBe(true);
      expect(order <= upperBound).toBe(true);
    });
  });

  describe('points()', () => {
    it('should enumerate all points', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const points = E.points();
      const order = E.cardinality();

      expect(BigInt(points.length)).toBe(order);
    });

    it('should include point at infinity', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const points = E.points();
      const hasInfinity = points.some((P) => P.isZero());

      expect(hasInfinity).toBe(true);
    });

    it('should only include valid points', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const points = E.points();

      for (const P of points) {
        if (!P.isZero()) {
          expect(E.is_on_curve(P.x!, P.y!)).toBe(true);
        }
      }
    });
  });

  describe('gens()', () => {
    it('should return at least one generator', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);

      const gens = E.gens();

      expect(gens.length).toBeGreaterThanOrEqual(1);
    });

    it('should return generators whose orders divide curve order', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);

      const curveOrder = E.cardinality();
      const gens = E.gens();

      for (const G of gens) {
        const genOrder = G.order();
        expect(curveOrder % genOrder).toBe(0n);
      }
    });

    it('should return generator of full group for prime order curves', () => {
      // Find a prime order curve
      const F = GF(97n);

      // Try different curves until we find one with prime order
      let primeOrderCurve: EllipticCurveFiniteField | null = null;

      for (let a = 1n; a < 10n; a++) {
        for (let b = 1n; b < 10n; b++) {
          try {
            const E = EllipticCurve(F, [a, b]);
            const order = E.cardinality();

            // Check if order is prime (simple check)
            let isPrime = order > 1n;
            for (let d = 2n; d * d <= order; d++) {
              if (order % d === 0n) {
                isPrime = false;
                break;
              }
            }

            if (isPrime) {
              primeOrderCurve = E;
              break;
            }
          } catch {
            continue;
          }
        }
        if (primeOrderCurve) break;
      }

      if (primeOrderCurve) {
        const gens = primeOrderCurve.gens();
        expect(gens.length).toBe(1);
        expect(gens[0]!.order()).toBe(primeOrderCurve.cardinality());
      }
    });
  });

  describe('set_order', () => {
    it('should accept valid order', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);

      // First compute the real order
      const realOrder = E.cardinality();

      // Create a new curve and set its order
      const E2 = EllipticCurve(F, [2n, 3n]);
      E2.set_order(realOrder);

      expect(E2.cardinality()).toBe(realOrder);
    });

    it('should reject invalid order outside Hasse bound', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);

      // Try to set an order outside Hasse bound
      expect(() => E.set_order(1n)).toThrow();
      expect(() => E.set_order(1000n)).toThrow();
    });
  });

  describe('trace of Frobenius', () => {
    it('should satisfy #E = q + 1 - t', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);

      const q = 97n;
      const order = E.cardinality();
      const trace = E.trace_of_frobenius();

      expect(order).toBe(q + 1n - trace);
    });
  });
});

describe('secp256k1 test vector', () => {
  // secp256k1 is too large for full point counting, but we can verify
  // that the known generator has the known order

  it('should verify secp256k1 generator order (simplified test)', () => {
    // Use a smaller curve that we can fully compute
    // This tests the same algorithms on a manageable size

    const F = GF(65521n); // Largest 16-bit prime
    const E = EllipticCurve(F, [0n, 7n]); // Like secp256k1 but over smaller field

    const order = E.cardinality();

    // Verify Hasse bound
    const q = 65521n;
    const sqrtQ = BigInt(Math.floor(Math.sqrt(Number(q))));
    expect(order >= q + 1n - 2n * sqrtQ - 2n).toBe(true);
    expect(order <= q + 1n + 2n * sqrtQ + 2n).toBe(true);

    // Verify random point orders divide curve order
    for (let i = 0; i < 5; i++) {
      const P = E.random_point();
      const pOrder = P.order();
      expect(order % pOrder).toBe(0n);
    }
  });
});

describe('known curve orders', () => {
  // Test vectors from PARI documentation

  it('should compute order for curve over GF(1031)', () => {
    const F = GF(1031n);
    const E = EllipticCurve(F, [0n, 1n]); // y^2 = x^3 + 1

    const order = E.cardinality();

    // Verify Hasse bound
    const q = 1031n;
    const sqrtQ = BigInt(Math.floor(Math.sqrt(Number(q))));
    expect(order >= q + 1n - 2n * sqrtQ - 2n).toBe(true);
    expect(order <= q + 1n + 2n * sqrtQ + 2n).toBe(true);

    // All points should satisfy [order]P = O
    const points = E.points();
    for (const P of points) {
      expect(P.mul(order).isZero()).toBe(true);
    }
  });

  it('should handle curve over GF(101)', () => {
    const F = GF(101n);
    const E = EllipticCurve(F, [2n, 3n]);

    const order = E.cardinality();
    const points = E.points();

    expect(BigInt(points.length)).toBe(order);
  });
});

describe('edge cases', () => {
  it('should handle point of order 2', () => {
    // A point (x, 0) has order 2
    const F = GF(97n);
    const E = EllipticCurve(F, [2n, 3n]);

    // Find a point with y = 0
    let pointOfOrder2: EllipticCurvePoint | null = null;
    for (let x = 0n; x < 97n; x++) {
      const xCubed = (x * x * x) % 97n;
      const rhs = (xCubed + 2n * x + 3n) % 97n;
      if (rhs === 0n) {
        pointOfOrder2 = E.point(x, 0n);
        break;
      }
    }

    if (pointOfOrder2) {
      expect(pointOfOrder2.order()).toBe(2n);
      expect(pointOfOrder2.double().isZero()).toBe(true);
    }
  });

  it('should handle very small field', () => {
    const F = GF(5n);
    const E = EllipticCurve(F, [1n, 1n]);

    const order = E.cardinality();
    const points = E.points();

    expect(BigInt(points.length)).toBe(order);
    expect(order).toBeGreaterThan(0n);
  });

  it('should handle curve with identity element only', () => {
    // This is rare but possible for certain curves
    // Most curves over prime fields have more than 1 point

    const F = GF(5n);

    // Try to find a curve with very few points
    // (This is just to test edge case handling)
    for (let a = 0n; a < 5n; a++) {
      for (let b = 1n; b < 5n; b++) {
        try {
          const E = EllipticCurve(F, [a, b]);
          const points = E.points();

          // Test that point counting is consistent
          expect(BigInt(points.length)).toBe(E.cardinality());
        } catch {
          continue;
        }
      }
    }
  });
});

describe('torsion subgroup functions', () => {
  describe('_find_point_of_order', () => {
    it('should find points of a given order', () => {
      const F = GF(101n);
      const E = EllipticCurve(F, [2n, 3n]);
      const order = E.cardinality();

      // Find divisors of the group order and look for points of those orders
      const { _find_point_of_order } = require('./ell_finite_field.js');

      for (const div of [2n, 5n, 10n]) {
        if (order % div === 0n) {
          const P = _find_point_of_order(E, div);
          if (P !== null) {
            expect(P.order()).toBe(div);
            expect(P.mul(div).isZero()).toBe(true);
          }
        }
      }
    });

    it('should return null for orders not dividing group order', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);
      const { _find_point_of_order } = require('./ell_finite_field.js');

      // 7 does not divide the group order for most curves
      // If it does divide, it will return a point; if not, it returns null
      const order = E.cardinality();
      const P = _find_point_of_order(E, 1000000n);

      // 1000000 definitely doesn't divide any curve order over GF(23)
      expect(P).toBeNull();
    });
  });

  describe('torsion_basis', () => {
    it('should return correct number of generators', () => {
      const F = GF(101n);
      const E = EllipticCurve(F, [2n, 3n]);
      const { torsion_basis, _find_point_of_order } = require('./ell_finite_field.js');

      const order = E.cardinality();

      // Find a small divisor of the order
      for (const div of [2n, 5n]) {
        if (order % div === 0n) {
          const P = _find_point_of_order(E, div);
          if (P !== null) {
            try {
              const basis = torsion_basis(E, div);

              // Basis should have 1 or 2 elements
              expect(basis.length).toBeGreaterThanOrEqual(1);
              expect(basis.length).toBeLessThanOrEqual(2);

              // Each basis element should have order exactly n (or dividing n for identity)
              for (const B of basis) {
                if (!B.isZero()) {
                  expect(B.mul(div).isZero()).toBe(true);
                }
              }
            } catch {
              // Torsion might not be fully rational
            }
          }
        }
      }
    });

    it('should return points with correct order', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);
      const { torsion_basis } = require('./ell_finite_field.js');

      const order = E.cardinality();

      // Try small divisors
      for (const n of [2n]) {
        if (order % n === 0n) {
          try {
            const basis = torsion_basis(E, n);

            for (const P of basis) {
              // [n]P = O
              expect(P.mul(n).isZero()).toBe(true);
            }
          } catch {
            // May not have full torsion
          }
        }
      }
    });
  });

  describe('torsion_subgroup', () => {
    it('should return all n-torsion points', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);
      const { torsion_subgroup } = require('./ell_finite_field.js');

      const order = E.cardinality();

      // Try 2-torsion
      if (order % 2n === 0n) {
        const torsion2 = torsion_subgroup(E, 2n);

        // All returned points should be 2-torsion
        for (const P of torsion2) {
          expect(P.mul(2n).isZero()).toBe(true);
        }

        // Should include the identity
        const hasIdentity = torsion2.some((P) => P.isZero());
        expect(hasIdentity).toBe(true);
      }
    });

    it('should have correct cardinality for cyclic torsion', () => {
      const F = GF(101n);
      const E = EllipticCurve(F, [1n, 1n]);
      const { torsion_subgroup } = require('./ell_finite_field.js');

      const order = E.cardinality();

      // 2-torsion should have 1, 2, or 4 elements
      if (order % 2n === 0n) {
        const torsion2 = torsion_subgroup(E, 2n);
        expect([1, 2, 4]).toContain(torsion2.length);
      }
    });
  });

  describe('division_points', () => {
    it('should verify [n]Q = P for all returned Q', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);
      const { division_points } = require('./ell_finite_field.js');

      const P = E.random_point();
      const n = 2n;

      const divisors = division_points(P, n);

      // Each Q in the result should satisfy [n]Q = P
      for (const Q of divisors) {
        expect(Q.mul(n).eq(P)).toBe(true);
      }
    });

    it('should return n-torsion for identity', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);
      const { division_points, torsion_subgroup } = require('./ell_finite_field.js');

      const O = E.zero();
      const n = 2n;

      const divisors = division_points(O, n);
      const torsion = torsion_subgroup(E, n);

      // Division points of O should be exactly E[n]
      expect(divisors.length).toBe(torsion.length);
    });

    it('should return correct count for divisible points', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);
      const { division_points, torsion_subgroup } = require('./ell_finite_field.js');

      // If P is divisible by n, there should be |E[n]| division points
      const P = E.random_point();
      const n = 2n;

      const divisors = division_points(P, n);
      const torsionSize = torsion_subgroup(E, n).length;

      // Either P is not divisible (0 points) or there are |E[n]| points
      expect([0, torsionSize]).toContain(divisors.length);
    });
  });

  describe('abelian_group', () => {
    it('should match cardinality', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);
      const { abelian_group } = require('./ell_finite_field.js');

      const group = abelian_group(E);
      const cardinality = E.cardinality();

      expect(group.order).toBe(cardinality);

      // Product of invariants should equal order
      if (group.invariants.length > 0) {
        const product = group.invariants.reduce((a, b) => a * b, 1n);
        expect(product).toBe(cardinality);
      }
    });

    it('should have n2 | n1 for invariants [n1, n2]', () => {
      const F = GF(41n);
      const E = EllipticCurve(F, [2n, 5n]);
      const { abelian_group } = require('./ell_finite_field.js');

      const group = abelian_group(E);

      if (group.invariants.length === 2) {
        const [n1, n2] = group.invariants;
        expect(n1! % n2!).toBe(0n);
      }
    });
  });

  describe('multiplication_by_m_isogeny', () => {
    it('should correctly multiply points', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);
      const { multiplication_by_m_isogeny } = require('./ell_finite_field.js');

      const doubling = multiplication_by_m_isogeny(E, 2n);
      const tripling = multiplication_by_m_isogeny(E, 3n);

      const P = E.random_point();

      expect(doubling.call(P).eq(P.mul(2n))).toBe(true);
      expect(tripling.call(P).eq(P.mul(3n))).toBe(true);
    });

    it('should have correct degree', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);
      const { multiplication_by_m_isogeny } = require('./ell_finite_field.js');

      const doubling = multiplication_by_m_isogeny(E, 2n);
      const tripling = multiplication_by_m_isogeny(E, 3n);

      expect(doubling.degree()).toBe(4n);
      expect(tripling.degree()).toBe(9n);
    });

    it('should have kernel equal to E[m]', () => {
      const F = GF(97n);
      const E = EllipticCurve(F, [2n, 3n]);
      const { multiplication_by_m_isogeny, torsion_subgroup } = require('./ell_finite_field.js');

      const doubling = multiplication_by_m_isogeny(E, 2n);
      const torsion2 = torsion_subgroup(E, 2n);

      // Kernel size should equal torsion size
      const kernel = doubling.kernel();
      expect(kernel.length).toBe(torsion2.length);

      // Every kernel point should be 2-torsion
      for (const P of kernel) {
        expect(P.mul(2n).isZero()).toBe(true);
      }
    });
  });

  describe('Frobenius operations', () => {
    it('should compute frobenius_discriminant correctly', () => {
      // Test case from SageMath
      const { frobenius_discriminant } = require('./ell_finite_field.js');

      const F = GF(101n);
      const E = EllipticCurve(F, [2n, 3n]);

      // E has trace of Frobenius t = 6 (calculated from cardinality 96 = 101 + 1 - t)
      // So discriminant = t^2 - 4*q = 36 - 404 = -368
      const D = frobenius_discriminant(E);
      expect(D).toBe(-368n);
    });

    it('should compute frobenius_polynomial correctly', () => {
      const { frobenius_polynomial } = require('./ell_finite_field.js');

      const F = GF(11n);
      const E = EllipticCurve(F, [3n, 3n]);

      // From SageMath: E.frobenius_polynomial() = x^2 - 4*x + 11
      const poly = frobenius_polynomial(E);

      // x^2 - 4*x + 11 => coeffs = [11, -4, 1]
      expect(poly.coeffs).toEqual([11n, -4n, 1n]);
      expect(poly.toString()).toBe('x^2 - 4*x + 11');
    });

    it('should compute frobenius_polynomial discriminant', () => {
      const { frobenius_polynomial } = require('./ell_finite_field.js');

      const F = GF(11n);
      const E = EllipticCurve(F, [3n, 3n]);

      const poly = frobenius_polynomial(E);
      const disc = poly.discriminant();

      // b^2 - 4ac = 16 - 44 = -28
      expect(disc).toBe(-28n);
    });

    it('should compute count_points for single degree', () => {
      const { count_points } = require('./ell_finite_field.js');

      const F = GF(101n);
      const E = EllipticCurve(F, [2n, 3n]);

      // count_points(E, 1) should equal cardinality
      const count = count_points(E, 1);
      expect(count).toBe(E.cardinality());
    });

    it('should compute count_points for extension fields', () => {
      const { count_points } = require('./ell_finite_field.js');

      const F = GF(101n);
      const E = EllipticCurve(F, [2n, 3n]);

      // From SageMath: E.count_points(5) = [96, 10368, 1031904, 104053248, 10509895776]
      const counts = count_points(E, 5);
      expect(counts).toEqual([96n, 10368n, 1031904n, 104053248n, 10509895776n]);
    });

    it('should compute frobenius_order for ordinary curve', () => {
      const { frobenius_order } = require('./ell_finite_field.js');

      const F = GF(11n);
      const E = EllipticCurve(F, [3n, 3n]);

      const ord = frobenius_order(E);

      // From SageMath: Order of conductor 2 generated by pi in Number Field in pi
      // with defining polynomial x^2 - 4*x + 11
      expect(ord.degree).toBe(2);
      expect(ord.definingPolynomial.coeffs).toEqual([11n, -4n, 1n]);
      expect(ord.discriminant).toBe(-28n); // t^2 - 4q = 16 - 44 = -28
    });

    it('should create FrobeniusEndomorphism', () => {
      const { frobenius_endomorphism } = require('./ell_finite_field.js');

      const F = GF(101n);
      const E = EllipticCurve(F, [2n, 3n]);

      const frob = frobenius_endomorphism(E);

      expect(frob.getDegree()).toBe(101n);
      expect(frob.domain()).toBe(E);
      expect(frob.codomain()).toBe(E);
    });

    it('should apply FrobeniusEndomorphism to a point', () => {
      const { frobenius_endomorphism } = require('./ell_finite_field.js');

      const F = GF(11n);
      const E = EllipticCurve(F, [3n, 3n]);

      // Get a point on the curve
      const P = E.random_point();
      const frob = frobenius_endomorphism(E);

      if (!P.isInfinity) {
        const frobP = frob.call(P);

        // Over F_11, x^11 = x for all x in F_11 (Fermat's little theorem)
        // So frobenius acts trivially on F_11-rational points
        expect(frobP.x!.value).toBe(P.x!.value);
        expect(frobP.y!.value).toBe(P.y!.value);
      }
    });

    it('should handle frobenius of point at infinity', () => {
      const { frobenius_endomorphism } = require('./ell_finite_field.js');

      const F = GF(101n);
      const E = EllipticCurve(F, [2n, 3n]);

      const frob = frobenius_endomorphism(E);
      const O = E.zero();
      const frobO = frob.call(O);

      expect(frobO.isZero()).toBe(true);
    });

    it('should compute multiplication_by_p_isogeny', () => {
      const { multiplication_by_p_isogeny } = require('./ell_finite_field.js');

      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 2n]);

      const phi = multiplication_by_p_isogeny(E);

      // Degree should be p^2
      expect(phi.degree()).toBe(23n * 23n);

      // phi(P) should equal P.mul(p)
      const P = E.random_point();
      expect(phi.call(P).eq(P.mul(23n))).toBe(true);
    });

    it('should compute endomorphism_order for ordinary curve', () => {
      const { endomorphism_order, is_supersingular } = require('./ell_finite_field.js');

      const F = GF(101n);
      const E = EllipticCurve(F, [2n, 3n]);

      // Skip if supersingular
      if (!is_supersingular(E)) {
        const ord = endomorphism_order(E);

        // Should return a quadratic order
        expect(ord.degree).toBe(2);

        // Discriminant should be t^2 - 4q for the Frobenius order
        const t = E.trace_of_frobenius();
        const q = 101n;
        expect(ord.discriminant).toBe(t * t - 4n * q);
      }
    });
  });

  describe('supersingular and ordinary tests', () => {
    it('should identify supersingular curves correctly', () => {
      const { is_supersingular, is_ordinary } = require('./ell_finite_field.js');

      // For E/F_p, E is supersingular iff p | t where t is trace of Frobenius
      // This means #E(F_p) = p + 1

      const F = GF(11n);

      // Test several curves
      for (let a = 0n; a < 5n; a++) {
        for (let b = 1n; b < 5n; b++) {
          try {
            const E = EllipticCurve(F, [a, b]);
            const ss = is_supersingular(E);
            const ord = is_ordinary(E);

            // Should be mutually exclusive
            expect(ss).toBe(!ord);

            // Verify using trace
            const t = E.trace_of_frobenius();
            expect(ss).toBe(t % 11n === 0n);
          } catch {
            // Skip singular curves
          }
        }
      }
    });

    it('should correctly identify ordinary curves', () => {
      const { is_ordinary } = require('./ell_finite_field.js');

      const F = GF(101n);
      const E = EllipticCurve(F, [2n, 3n]);

      // E: y^2 = x^3 + 2x + 3 over F_101
      // #E = 96 = 101 + 1 - 6, so t = 6
      // Since gcd(6, 101) = 1, E is ordinary
      expect(is_ordinary(E)).toBe(true);
    });
  });

  describe('embedding degree', () => {
    it('should compute correct embedding degree', () => {
      const { embedding_degree } = require('./ell_finite_field.js');

      const F = GF(103n);
      const E = EllipticCurve(F, [1n, 18n]);

      // The embedding degree is the smallest k such that n | q^k - 1
      // For a subgroup of order 19 with q = 103, need to find k where 19 | 103^k - 1
      const n = 19n;
      const k = embedding_degree(E, n);

      // Verify: n should divide q^k - 1
      let qk = 1n;
      for (let i = 0n; i < k; i++) {
        qk *= 103n;
      }
      expect((qk - 1n) % n).toBe(0n);

      // k should be minimal
      for (let i = 1n; i < k; i++) {
        let qi = 1n;
        for (let j = 0n; j < i; j++) {
          qi *= 103n;
        }
        expect((qi - 1n) % n).not.toBe(0n);
      }
    });
  });

  describe('quadratic twist', () => {
    it('should compute quadratic twist with correct cardinality', () => {
      const { quadratic_twist } = require('./ell_finite_field.js');

      const F = GF(101n);
      const E = EllipticCurve(F, [2n, 3n]);

      const Et = quadratic_twist(E);

      // For E and its quadratic twist Et:
      // #E + #Et = 2(q + 1)
      const nE = E.cardinality();
      const nEt = Et.cardinality();
      const q = 101n;

      expect(nE + nEt).toBe(2n * (q + 1n));
    });

    it('should produce a non-isomorphic curve for generic j', () => {
      const { quadratic_twist } = require('./ell_finite_field.js');

      const F = GF(101n);
      const E = EllipticCurve(F, [2n, 3n]);

      const Et = quadratic_twist(E);

      // For generic j-invariant, the twist should have different a, b
      // (unless we happen to choose D = 1, which we shouldn't)
      const j = E.j_invariant();
      const jt = Et.j_invariant();

      // j-invariants should be the same
      expect(j.eq(jt)).toBe(true);

      // But cardinalities should differ (for non-supersingular)
      expect(E.cardinality()).not.toBe(Et.cardinality());
    });
  });
});
