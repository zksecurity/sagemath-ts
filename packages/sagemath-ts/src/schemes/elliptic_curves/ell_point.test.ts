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
import { EllipticCurve as EllipticCurveGeneric_ } from './constructor.js';
import { EllipticCurve, EllipticCurveFiniteField, EllipticCurvePoint } from './ell_finite_field.js';
import {
  type EllipticCurvePoint as EllipticCurvePointGeneric,
  ate_pairing as generic_ate_pairing,
  division_points as generic_division_points,
  is_divisible_by as generic_is_divisible_by,
  point_log as generic_point_log,
  tate_pairing as generic_tate_pairing,
  weil_pairing as generic_weil_pairing,
} from './ell_point.js';

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

  describe('has_order', () => {
    const F23 = GF(23n);
    const E = EllipticCurve(F23, [1n, 1n]);

    it('should return true for exact order', () => {
      const P = E.point(0n, 1n);
      const order = P.order();
      expect(P.has_order(order)).toBe(true);
    });

    it('should return false for multiples of order', () => {
      const P = E.point(0n, 1n);
      const order = P.order();
      // 2*order satisfies (2*order)*P = O, but is not the exact order
      expect(P.has_order(2n * order)).toBe(false);
    });

    it('should return false for divisors of order (when order is composite)', () => {
      // Find a point with composite order
      const points = E.points().filter((p) => !p.isZero());
      for (const P of points) {
        const order = P.order();
        // If order > 1 and has a proper divisor
        for (let d = 2n; d * d <= order; d++) {
          if (order % d === 0n) {
            // d is a proper divisor of order
            // d*P should NOT be O (so P does not have order d)
            expect(P.has_order(d)).toBe(false);
            break;
          }
        }
      }
    });

    it('should correctly verify order for point at infinity', () => {
      const O = E.zero();
      expect(O.has_order(1n)).toBe(true);
      expect(O.has_order(2n)).toBe(false);
    });

    it('should return false for invalid orders', () => {
      const P = E.point(0n, 1n);
      expect(P.has_order(0n)).toBe(false);
      expect(P.has_order(-1n)).toBe(false);
    });

    it('should distinguish between order and its multiples for all points', () => {
      // This is the key test - the bug was that has_order only checked n*P = O
      // which would pass for ANY multiple of the actual order
      const curveOrder = E.cardinality();
      const points = E.points()
        .filter((p) => !p.isZero())
        .slice(0, 10);

      for (const P of points) {
        const actualOrder = P.order();
        // has_order should be true ONLY for the actual order
        expect(P.has_order(actualOrder)).toBe(true);

        // For curve order (which is always a multiple of point order),
        // has_order should be false unless it equals the actual order
        if (curveOrder !== actualOrder) {
          expect(P.has_order(curveOrder)).toBe(false);
        }
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

// ============================================================================
// ell_point.ts entry points (audit M96): division_points, is_divisible_by,
// point_log and the exported pairings.
// ============================================================================

describe('ell_point entry points', () => {
  describe('division_points', () => {
    it('matches brute force on every curve and point over GF(11) and GF(13)', () => {
      // Sage computes division_points from the distinct roots of the division
      // polynomial; the old implementation returned every solution twice.
      const key = (P: EllipticCurvePointGeneric) => (P.is_zero() ? 'O' : `${P.x()},${P.y()}`);

      let cases = 0;
      for (const p of [11n, 13n]) {
        const K = GF(p);
        for (let a = 0n; a < p; a++) {
          for (let b = 0n; b < p; b++) {
            let E: ReturnType<typeof EllipticCurveGeneric_>;
            try {
              E = EllipticCurveGeneric_(K, [a, b]);
            } catch {
              continue;
            }
            const pts = E.torsion_points();
            for (const P of pts) {
              for (const m of [2n, 3n, 4n]) {
                const expected = pts
                  .filter((Q) => Q.mul(m).eq(P))
                  .map(key)
                  .sort();
                const got = generic_division_points(P, m).map(key);
                expect([...got].sort()).toEqual(expected);
                cases++;
              }
            }
          }
        }
      }
      expect(cases).toBeGreaterThan(10000);
    });

    it('is sorted by (Z, X, Y) like Sage', () => {
      const K = GF(11n);
      const E = EllipticCurveGeneric_(K, [0n, 1n]);
      const pts = generic_division_points(E.zero(), 2n);
      // The identity comes first, then the affine points by increasing x.
      expect(pts[0]!.is_zero()).toBe(true);
      for (let i = 2; i < pts.length; i++) {
        const prev = pts[i - 1]!;
        const cur = pts[i]!;
        const px = (prev.x() as unknown as { value: bigint }).value;
        const cx = (cur.x() as unknown as { value: bigint }).value;
        const py = (prev.y() as unknown as { value: bigint }).value;
        const cy = (cur.y() as unknown as { value: bigint }).value;
        expect(px < cx || (px === cx && py <= cy)).toBe(true);
      }
    });

    it('works in characteristic 2', () => {
      // y^2 + y = x^3 over GF(2): (0:1:1) doubles to (0:0:1).
      const F2 = GF(2n);
      const E = EllipticCurveGeneric_(F2, [0n, 0n, 1n, 0n, 0n]);
      const P = E.point([F2.__call__(0n), F2.__call__(0n)]);
      const pts = generic_division_points(P, 2n);
      expect(pts.length).toBe(1);
      expect(pts[0]!.mul(2n).eq(P)).toBe(true);
    });

    it('returns [m*P] for m = 1 and m = -1', () => {
      const K = GF(11n);
      const E = EllipticCurveGeneric_(K, [0n, 1n]);
      const P = E.point([K.__call__(0n), K.__call__(1n)]);
      expect(generic_division_points(P, 1n)[0]!.eq(P)).toBe(true);
      expect(generic_division_points(P, -1n)[0]!.eq(P.neg())).toBe(true);
    });
  });

  describe('is_divisible_by', () => {
    it("matches Sage's doctest over GF(101) for E = [23, 34]", () => {
      // sage: E = EllipticCurve(GF(101), [23, 34])
      // sage: len([T for T in E.points() if T.is_divisible_by(2)])  ->  53
      // sage: len([T for T in E.points() if T.is_divisible_by(3)])  ->  106
      const K = GF(101n);
      const E = EllipticCurveGeneric_(K, [23n, 34n]);
      const pts = E.torsion_points();
      expect(pts.length).toBe(106);
      expect(pts.filter((T) => generic_is_divisible_by(T, 2n)).length).toBe(53);
      expect(pts.filter((T) => generic_is_divisible_by(T, 3n)).length).toBe(106);
    });

    it('agrees with division_points on every curve over GF(11)', () => {
      const p = 11n;
      const K = GF(p);
      for (let a = 0n; a < p; a++) {
        for (let b = 0n; b < p; b++) {
          let E: ReturnType<typeof EllipticCurveGeneric_>;
          try {
            E = EllipticCurveGeneric_(K, [a, b]);
          } catch {
            continue;
          }
          for (const P of E.torsion_points()) {
            for (const m of [2n, 3n, 4n]) {
              expect(generic_is_divisible_by(P, m)).toBe(generic_division_points(P, m).length > 0);
            }
          }
        }
      }
    });
  });

  describe('point_log', () => {
    it('recovers every exponent on curves over GF(11) and GF(23)', () => {
      for (const p of [11n, 23n]) {
        const K = GF(p);
        for (let a = 0n; a < p; a += 3n) {
          for (let b = 0n; b < p; b += 5n) {
            let E: ReturnType<typeof EllipticCurveGeneric_>;
            try {
              E = EllipticCurveGeneric_(K, [a, b]);
            } catch {
              continue;
            }
            for (const P of E.torsion_points()) {
              if (P.is_zero()) continue;
              const ord = P.order();
              for (let x = 0n; x < ord; x++) {
                expect(generic_point_log(P.mul(x), P)).toBe(x);
              }
            }
          }
        }
      }
    });

    it('raises when the target is not in the subgroup', () => {
      const K = GF(11n);
      const E = EllipticCurveGeneric_(K, [1n, 9n]); // Z/4 x Z/2
      const pts = E.torsion_points().filter((P) => !P.is_zero());
      const P = pts.find((Q) => Q.order() === 4n)!;
      const Q = pts.find((R) => R.order() === 2n && !R.eq(P.mul(2n)))!;
      expect(() => generic_point_log(Q, P)).toThrow(
        'ECDLog problem has no solution (non-trivial Weil pairing)'
      );
    });
  });

  describe('pairings', () => {
    // E: y^2 = x^3 + 7*x over GF(13); #E = 18, and E[3] is fully rational.
    const K = GF(13n);
    const E = EllipticCurveGeneric_(K, [7n, 0n]);
    const n = 3n;
    const t = -4n;
    const tors = E.torsion_points().filter((P) => P.mul(n).is_zero());
    const P = tors.find((X) => !X.is_zero())!;
    const Q = tors.find((X) => !X.is_zero() && !generic_weil_pairing(P, X, n).eq(K.one()))!;

    it('weil_pairing is a primitive n-th root of unity on an independent pair', () => {
      const w = generic_weil_pairing(P, Q, n);
      expect(w.pow(n).eq(K.one())).toBe(true);
      expect(w.eq(K.one())).toBe(false);
      // Alternating: e(Q,P) = e(P,Q)^-1
      expect(generic_weil_pairing(Q, P, n).eq(w.inv())).toBe(true);
      // e(P,P) = 1
      expect(generic_weil_pairing(P, P, n).eq(K.one())).toBe(true);
    });

    it('weil_pairing rejects points that are not n-torsion', () => {
      const R = E.torsion_points().find((X) => !X.mul(n).is_zero())!;
      expect(() => generic_weil_pairing(P, R, n)).toThrow('points must both be n-torsion');
    });

    it('tate_pairing rejects a non n-torsion P and a bad q', () => {
      expect(() => generic_tate_pairing(P, Q, 4n, 1, 13n)).toThrow('The point P must be n-torsion');
      expect(() => generic_tate_pairing(P, Q, n, 1, 5n)).toThrow(
        'n does not divide (q^k - 1) for the supplied value of q'
      );
    });

    it('ate_pairing equals tate_pairing^M (Sage identity, negative trace)', () => {
      // c = (k*q^(k-1)) mod n; T = t-1; N = gcd(T^k-1, q^k-1); s = N/n;
      // L = (T^k-1)/N; M = (L*s*c^-1) mod n
      // ate(P,Q,n,k,t) == tate(Q,P,n,k)^M     [ell_point.py:2681-2683]
      const q = 13n;
      const T = t - 1n;
      const gcdBig = (x: bigint, y: bigint): bigint => {
        let a = x < 0n ? -x : x;
        let b = y < 0n ? -y : y;
        while (b) [a, b] = [b, a % b];
        return a;
      };
      const N = gcdBig(T - 1n, q - 1n);
      const s = N / n;
      const L = (T - 1n) / N;
      const c = 1n % n;
      let cinv = 1n;
      while ((c * cinv) % n !== 1n % n) cinv++;
      const M = ((((L % n) * (s % n)) % n) * cinv) % n;

      const ate = generic_ate_pairing(P, Q, n, 1, t, q);
      const tate = generic_tate_pairing(Q, P, n, 1, q);
      expect(ate.eq(tate.pow(M))).toBe(true);
    });

    it('ate_pairing rejects a point of the wrong order', () => {
      const R = E.torsion_points().find((X) => !X.mul(n).is_zero() && !X.is_zero())!;
      expect(() => generic_ate_pairing(R, Q, n, 1, t, 13n)).toThrow('is not of order n=3');
    });
  });
});
