/**
 * Tests for elliptic curve point arithmetic
 *
 * Tests the EllipticCurvePoint class implementation covering:
 * - Point addition
 * - Point negation
 * - Point doubling
 * - Scalar multiplication
 * - Group law associativity
 * - Known test vectors (secp256k1)
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import { EllipticCurve, EllipticCurveFiniteField, EllipticCurvePoint } from './ell_finite_field.js';

describe('EllipticCurvePoint', () => {
  describe('basic point operations', () => {
    // Curve y^2 = x^3 + x + 1 over GF(23)
    const F23 = GF(23n);
    const E = EllipticCurve(F23, [1n, 1n]);

    it('should create point at infinity', () => {
      const O = E.zero();
      expect(O.isZero()).toBe(true);
      expect(O.isInfinity).toBe(true);
    });

    it('should create affine points', () => {
      // (0, 1) is on y^2 = x^3 + x + 1 since 1 = 0 + 0 + 1
      const P = E.point(0n, 1n);
      expect(P.isZero()).toBe(false);
      expect(P.x!.value).toBe(0n);
      expect(P.y!.value).toBe(1n);
    });

    it('should reject points not on the curve', () => {
      expect(() => E.point(0n, 0n)).toThrow();
    });

    it('should verify points on curve', () => {
      // Check that (0, 1) is on the curve
      expect(E.is_on_curve(F23.__call__(0n), F23.__call__(1n))).toBe(true);
      expect(E.is_on_curve(F23.__call__(0n), F23.__call__(22n))).toBe(true); // -1 = 22 mod 23
      expect(E.is_on_curve(F23.__call__(0n), F23.__call__(0n))).toBe(false);
    });
  });

  describe('point negation', () => {
    const F23 = GF(23n);
    const E = EllipticCurve(F23, [1n, 1n]);

    it('should negate the point at infinity to itself', () => {
      const O = E.zero();
      expect(O.neg().eq(O)).toBe(true);
    });

    it('should negate affine points correctly', () => {
      // -P = (x, -y) for short Weierstrass
      const P = E.point(0n, 1n);
      const negP = P.neg();

      expect(negP.x!.value).toBe(0n);
      expect(negP.y!.value).toBe(22n); // -1 mod 23 = 22
    });

    it('should satisfy -(-P) = P', () => {
      const P = E.point(0n, 1n);
      expect(P.neg().neg().eq(P)).toBe(true);
    });
  });

  describe('point addition', () => {
    const F23 = GF(23n);
    const E = EllipticCurve(F23, [1n, 1n]);

    it('should satisfy P + O = P', () => {
      const P = E.point(0n, 1n);
      const O = E.zero();

      expect(P.add(O).eq(P)).toBe(true);
      expect(O.add(P).eq(P)).toBe(true);
    });

    it('should satisfy P + (-P) = O', () => {
      const P = E.point(0n, 1n);
      const negP = P.neg();

      expect(P.add(negP).isZero()).toBe(true);
    });

    it('should add distinct points correctly', () => {
      // Find another point on the curve
      const P = E.point(0n, 1n);

      // Double P first to get a different point
      const Q = P.double();

      // P + Q should be a new point
      const R = P.add(Q);
      expect(R.isZero()).toBe(false);

      // Verify result is on the curve
      expect(E.is_on_curve(R.x!, R.y!)).toBe(true);
    });
  });

  describe('point doubling', () => {
    const F23 = GF(23n);
    const E = EllipticCurve(F23, [1n, 1n]);

    it('should double the point at infinity to itself', () => {
      const O = E.zero();
      expect(O.double().eq(O)).toBe(true);
    });

    it('should compute 2P correctly', () => {
      const P = E.point(0n, 1n);
      const twoP = P.double();

      expect(twoP.isZero()).toBe(false);
      expect(E.is_on_curve(twoP.x!, twoP.y!)).toBe(true);
    });

    it('should satisfy P + P = 2P', () => {
      const P = E.point(0n, 1n);
      expect(P.add(P).eq(P.double())).toBe(true);
    });
  });

  describe('scalar multiplication', () => {
    const F23 = GF(23n);
    const E = EllipticCurve(F23, [1n, 1n]);

    it('should compute 0*P = O', () => {
      const P = E.point(0n, 1n);
      expect(P.mul(0n).isZero()).toBe(true);
    });

    it('should compute 1*P = P', () => {
      const P = E.point(0n, 1n);
      expect(P.mul(1n).eq(P)).toBe(true);
    });

    it('should compute 2*P = P + P', () => {
      const P = E.point(0n, 1n);
      expect(P.mul(2n).eq(P.add(P))).toBe(true);
    });

    it('should compute 3*P = P + P + P', () => {
      const P = E.point(0n, 1n);
      expect(P.mul(3n).eq(P.add(P).add(P))).toBe(true);
    });

    it('should handle negative scalars: (-n)*P = n*(-P)', () => {
      const P = E.point(0n, 1n);
      const negP = P.neg();

      expect(P.mul(-3n).eq(negP.mul(3n))).toBe(true);
    });

    it('should satisfy n*O = O', () => {
      const O = E.zero();
      expect(O.mul(100n).isZero()).toBe(true);
    });
  });

  describe('subtraction', () => {
    const F23 = GF(23n);
    const E = EllipticCurve(F23, [1n, 1n]);

    it('should compute P - P = O', () => {
      const P = E.point(0n, 1n);
      expect(P.sub(P).isZero()).toBe(true);
    });

    it('should compute P - O = P', () => {
      const P = E.point(0n, 1n);
      const O = E.zero();
      expect(P.sub(O).eq(P)).toBe(true);
    });

    it('should satisfy P - Q = P + (-Q)', () => {
      const P = E.point(0n, 1n);
      const Q = P.double();

      expect(P.sub(Q).eq(P.add(Q.neg()))).toBe(true);
    });
  });

  describe('group law associativity', () => {
    // Use a small curve to test all points
    const F7 = GF(7n);
    // y^2 = x^3 + 3 over GF(7)
    // Points: O, (1, 2), (1, 5), (2, 3), (2, 4), (4, 1), (4, 6), (5, 2), (5, 5), (6, 0)
    const E = EllipticCurve(F7, [0n, 3n]);

    it('should satisfy (P + Q) + R = P + (Q + R)', () => {
      // Get some points on the curve
      const points = E.points()
        .filter((p) => !p.isZero())
        .slice(0, 3);

      if (points.length >= 3) {
        const [P, Q, R] = points;

        const left = P!.add(Q!).add(R!);
        const right = P!.add(Q!.add(R!));

        expect(left.eq(right)).toBe(true);
      }
    });

    it('should satisfy commutativity: P + Q = Q + P', () => {
      const points = E.points()
        .filter((p) => !p.isZero())
        .slice(0, 2);

      if (points.length >= 2) {
        const [P, Q] = points;
        expect(P!.add(Q!).eq(Q!.add(P!))).toBe(true);
      }
    });
  });

  describe('point order', () => {
    const F23 = GF(23n);
    const E = EllipticCurve(F23, [1n, 1n]);

    it('should have order 1 for point at infinity', () => {
      const O = E.zero();
      expect(O.order()).toBe(1n);
    });

    it('should compute finite order for non-identity points', () => {
      const P = E.point(0n, 1n);
      const order = P.order();

      // Order should be positive
      expect(order).toBeGreaterThan(1n);

      // nP should be O
      expect(P.mul(order).isZero()).toBe(true);
    });

    it('should satisfy nP = O for point order n', () => {
      const points = E.points()
        .filter((p) => !p.isZero())
        .slice(0, 5);

      for (const P of points) {
        const n = P.order();
        expect(P.mul(n).isZero()).toBe(true);
      }
    });
  });

  describe('secp256k1 test vectors', () => {
    // secp256k1 parameters
    // p = 2^256 - 2^32 - 977
    const p = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2fn;
    const Fp = GF(p);

    // Curve: y^2 = x^3 + 7
    const secp256k1 = EllipticCurve(Fp, [0n, 7n]);

    // Generator point G
    const Gx = 0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798n;
    const Gy = 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8n;

    // Group order n
    const n = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141n;

    it('should verify generator point is on curve', () => {
      expect(secp256k1.is_on_curve(Fp.__call__(Gx), Fp.__call__(Gy))).toBe(true);
    });

    it('should create generator point', () => {
      const G = secp256k1.point(Gx, Gy);
      expect(G.x!.value).toBe(Gx);
      expect(G.y!.value).toBe(Gy);
    });

    it('should compute 2G correctly', () => {
      const G = secp256k1.point(Gx, Gy);
      const twoG = G.double();

      // Known value of 2G for secp256k1
      const expectedX = 0xc6047f9441ed7d6d3045406e95c07cd85c778e4b8cef3ca7abac09b95c709ee5n;
      const expectedY = 0x1ae168fea63dc339a3c58419466ceaeef7f632653266d0e1236431a950cfe52an;

      expect(twoG.x!.value).toBe(expectedX);
      expect(twoG.y!.value).toBe(expectedY);
    });

    it('should satisfy nG = O (group order)', () => {
      const G = secp256k1.point(Gx, Gy);
      secp256k1.set_order(n, false); // Set known order

      // nG should be point at infinity
      const result = G.mul(n);
      expect(result.isZero()).toBe(true);
    });

    it('should satisfy (n-1)G + G = O', () => {
      const G = secp256k1.point(Gx, Gy);
      const nMinus1G = G.mul(n - 1n);

      expect(nMinus1G.add(G).isZero()).toBe(true);
    });

    it('should verify scalar multiplication test vector', () => {
      const G = secp256k1.point(Gx, Gy);

      // Test vector: 7*G
      const expected7Gx = 0x5cbdf0646e5db4eaa398f365f2ea7a0e3d419b7e0330e39ce92bddedcac4f9bcn;
      const expected7Gy = 0x6aebca40ba255960a3178d6d861a54dba813d0b813fde7b5a5082628087264dan;

      const sevenG = G.mul(7n);
      expect(sevenG.x!.value).toBe(expected7Gx);
      expect(sevenG.y!.value).toBe(expected7Gy);
    });
  });

  describe('curve over small field - exhaustive tests', () => {
    // Use a small curve to exhaustively test
    // y^2 = x^3 + 2x + 3 over GF(97)
    const F97 = GF(97n);
    const E = EllipticCurve(F97, [2n, 3n]);

    it('should enumerate all points', () => {
      const points = E.points();

      // Count should include point at infinity
      expect(points.length).toBeGreaterThan(0);

      // First point should be infinity
      expect(points[0]!.isZero()).toBe(true);

      // All other points should be on the curve
      for (const P of points) {
        if (!P.isZero()) {
          expect(E.is_on_curve(P.x!, P.y!)).toBe(true);
        }
      }
    });

    it('should compute curve cardinality', () => {
      const order = E.cardinality();
      const points = E.points();

      expect(order).toBe(BigInt(points.length));
    });

    it('should verify all points have order dividing curve order', () => {
      const curveOrder = E.cardinality();
      const points = E.points()
        .filter((p) => !p.isZero())
        .slice(0, 10); // Test first 10

      for (const P of points) {
        const result = P.mul(curveOrder);
        expect(result.isZero()).toBe(true);
      }
    });
  });

  describe('edge cases', () => {
    const F23 = GF(23n);
    const E = EllipticCurve(F23, [1n, 1n]);

    it('should handle point with y=0 (order 2 points)', () => {
      // Find a point with y=0 if it exists
      const points = E.points();
      const yZeroPoints = points.filter((p) => !p.isZero() && p.y!.isZero());

      for (const P of yZeroPoints) {
        // 2P should be O for a point with y=0
        expect(P.double().isZero()).toBe(true);
        expect(P.add(P).isZero()).toBe(true);
      }
    });

    it('should correctly compute large scalar multiples', () => {
      const P = E.point(0n, 1n);
      const order = P.order();

      // (order + 1)*P should equal P
      expect(P.mul(order + 1n).eq(P)).toBe(true);

      // (2*order)*P should be O
      expect(P.mul(2n * order).isZero()).toBe(true);
    });
  });
});

/**
 * Note: The functions height(), archimedean_local_height(), non_archimedean_local_height(),
 * point_has_good_reduction(), and reduction() in ell_point.ts are designed for elliptic
 * curves over number fields (using the generic EllipticCurvePoint interface).
 *
 * For curves over finite fields (using EllipticCurveFiniteField), these functions
 * return trivial values (0 for heights, true for good reduction) because:
 * - Over finite fields, all points are torsion (height = 0)
 * - Finite fields have no archimedean places
 * - Points are already in their "reduced" form
 *
 * The tests below verify the main functionality using EllipticCurveFiniteField from
 * ell_finite_field.ts. The height functions can be tested when number field support
 * is implemented.
 */

describe('EllipticCurveFiniteField', () => {
  describe('curve creation', () => {
    const F23 = GF(23n);

    it('should create a non-singular curve', () => {
      const E = EllipticCurve(F23, [1n, 1n]);
      expect(E).toBeDefined();
    });

    it('should reject singular curves', () => {
      // Curve y^2 = x^3 (a=0, b=0) is singular
      expect(() => EllipticCurve(F23, [0n, 0n])).toThrow();
    });

    it('should compute discriminant', () => {
      const E = EllipticCurve(F23, [1n, 1n]);
      const disc = E.discriminant();

      // Discriminant should be non-zero for valid curve
      expect(disc.isZero()).toBe(false);
    });

    it('should compute j-invariant', () => {
      const E = EllipticCurve(F23, [1n, 1n]);
      const j = E.j_invariant();

      // j-invariant should be defined
      expect(j).toBeDefined();
    });
  });

  describe('lift_x', () => {
    const F23 = GF(23n);
    const E = EllipticCurve(F23, [1n, 1n]);

    it('should find points with given x-coordinate', () => {
      const P = E.lift_x(0n);

      expect(P.x!.value).toBe(0n);
      expect(E.is_on_curve(P.x!, P.y!)).toBe(true);
    });

    it('should find all points with given x-coordinate', () => {
      const points = E.lift_x(0n, true);

      // Should return 2 points (y and -y) unless y=0
      expect(points.length).toBeGreaterThanOrEqual(1);
      expect(points.length).toBeLessThanOrEqual(2);

      for (const P of points) {
        expect(P.x!.value).toBe(0n);
        expect(E.is_on_curve(P.x!, P.y!)).toBe(true);
      }
    });

    it('should throw for invalid x-coordinate', () => {
      // Find an x that has no point
      let invalidX: bigint | null = null;
      for (let x = 0n; x < 23n; x++) {
        if (!E.is_x_coord(x)) {
          invalidX = x;
          break;
        }
      }

      if (invalidX !== null) {
        expect(() => E.lift_x(invalidX!)).toThrow();
        expect(E.lift_x(invalidX!, true)).toEqual([]);
      }
    });
  });

  describe('random_point', () => {
    const F23 = GF(23n);
    const E = EllipticCurve(F23, [1n, 1n]);

    it('should return a point on the curve', () => {
      const P = E.random_point();

      expect(P.isZero()).toBe(false);
      expect(E.is_on_curve(P.x!, P.y!)).toBe(true);
    });
  });

  describe('trace of Frobenius', () => {
    const F23 = GF(23n);
    const E = EllipticCurve(F23, [1n, 1n]);

    it('should satisfy Hasse bound', () => {
      const q = 23n;
      const t = E.trace_of_frobenius();

      // |t| <= 2*sqrt(q)
      const sqrtQ = BigInt(Math.floor(Math.sqrt(Number(q))));
      expect(t <= 2n * sqrtQ + 1n).toBe(true);
      expect(t >= -(2n * sqrtQ + 1n)).toBe(true);
    });

    it('should satisfy #E = q + 1 - t', () => {
      const q = 23n;
      const t = E.trace_of_frobenius();
      const order = E.cardinality();

      expect(order).toBe(q + 1n - t);
    });
  });
});
