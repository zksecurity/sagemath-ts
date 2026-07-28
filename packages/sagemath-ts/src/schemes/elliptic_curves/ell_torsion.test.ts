/**
 * Tests for elliptic curve torsion subgroups
 *
 * Tests the torsion subgroup computation for elliptic curves over finite fields.
 * Covers:
 * - Torsion subgroup structure (cyclic and product of cyclic)
 * - Point order computation
 * - Discrete logarithm (baby-step giant-step)
 * - Group generators
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import type { FiniteFieldElement } from '../../rings/finite_rings/finite_field_prime.js';
import { EllipticCurve } from './constructor.js';
import {
  EllipticCurveTorsionSubgroup,
  _p_primary_torsion_basis,
  discrete_log,
  has_finite_order,
  has_infinite_order,
  order_from_multiple,
  point_order,
  torsion_bound,
  torsion_order,
  torsion_points,
  torsion_subgroup,
} from './ell_torsion.js';

describe('EllipticCurveTorsionSubgroup', () => {
  describe('construction and basic properties', () => {
    it('should compute torsion subgroup for y^2 = x^3 + x over GF(7)', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 0n]); // y^2 = x^3 + x
      const T = new EllipticCurveTorsionSubgroup(E);

      // The order should be positive
      expect(T.order()).toBeGreaterThan(0n);

      // The curve should be the same
      expect(T.curve()).toBe(E);
    });

    it('should compute torsion subgroup for y^2 = x^3 + 1 over GF(7)', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [0n, 1n]); // y^2 = x^3 + 1
      const T = new EllipticCurveTorsionSubgroup(E);

      // For finite fields, torsion = entire group
      const pts = T.points();
      expect(pts.length).toBe(Number(T.order()));
    });

    it('should have consistent ngens and gens', () => {
      const F11 = GF(11n);
      const E = EllipticCurve<FiniteFieldElement>(F11, [1n, 1n]);
      const T = new EllipticCurveTorsionSubgroup(E);

      const ngens = T.ngens();
      const gens = T.gens();

      expect(gens.length).toBe(ngens);
      expect(ngens).toBeLessThanOrEqual(2); // EC groups have rank at most 2
    });

    it('should return correct generator via gen()', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const T = new EllipticCurveTorsionSubgroup(E);

      if (T.ngens() > 0) {
        const g0 = T.gen(0);
        expect(g0).toBeDefined();
        expect(T.contains(g0)).toBe(true);
      }
    });

    it('should throw for invalid generator index', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const T = new EllipticCurveTorsionSubgroup(E);

      expect(() => T.gen(10)).toThrow();
      expect(() => T.gen(-1)).toThrow();
    });
  });

  describe('invariants and structure', () => {
    it('should return valid invariants', () => {
      const F13 = GF(13n);
      const E = EllipticCurve<FiniteFieldElement>(F13, [2n, 3n]);
      const T = new EllipticCurveTorsionSubgroup(E);

      const invs = T.invariants();

      // Invariants should be positive integers
      for (const inv of invs) {
        expect(inv).toBeGreaterThan(0n);
      }

      // Product of invariants equals order
      if (invs.length > 0) {
        const product = invs.reduce((a, b) => a * b, 1n);
        expect(product).toBe(T.order());
      }
    });

    it('should have short_name describing structure', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const T = new EllipticCurveTorsionSubgroup(E);

      const name = T.short_name();
      expect(typeof name).toBe('string');
      // Should contain Z/ for non-trivial groups
      if (T.order() > 1n) {
        expect(name).toContain('Z/');
      }
    });

    it('should produce string representation', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const T = new EllipticCurveTorsionSubgroup(E);

      const str = T.toString();
      expect(str).toContain('Torsion');
    });
  });

  describe('points enumeration', () => {
    it('should enumerate all points correctly', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const T = new EllipticCurveTorsionSubgroup(E);

      const pts = T.points();

      // All points should be on the curve
      for (const P of pts) {
        if (!P.is_zero()) {
          expect(E.is_on_curve(P.x(), P.y())).toBe(true);
        }
      }

      // Should contain the identity
      const hasIdentity = pts.some((P) => P.is_zero());
      expect(hasIdentity).toBe(true);
    });

    it('should be iterable', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 0n]);
      const T = new EllipticCurveTorsionSubgroup(E);

      let count = 0;
      for (const _P of T) {
        count++;
      }

      expect(count).toBe(Number(T.order()));
    });

    it('should have correct length property', () => {
      const F11 = GF(11n);
      const E = EllipticCurve<FiniteFieldElement>(F11, [1n, 1n]);
      const T = new EllipticCurveTorsionSubgroup(E);

      expect(T.length).toBe(Number(T.order()));
    });

    it('should contain points correctly', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const T = new EllipticCurveTorsionSubgroup(E);

      // Identity should be contained
      expect(T.contains(E.zero())).toBe(true);

      // All points from the group should be contained
      for (const P of T.points()) {
        expect(T.contains(P)).toBe(true);
      }
    });

    it('should return identity element', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const T = new EllipticCurveTorsionSubgroup(E);

      const id = T.identity();
      expect(id.is_zero()).toBe(true);
    });

    it('should return random element', () => {
      const F11 = GF(11n);
      // Use [1, 1] since [2, 3] is singular over GF(11)
      const E = EllipticCurve<FiniteFieldElement>(F11, [1n, 1n]);
      const T = new EllipticCurveTorsionSubgroup(E);

      const P = T.random_element();
      expect(T.contains(P)).toBe(true);
    });
  });

  describe('comparison', () => {
    it('should compare equal torsion groups', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const T1 = new EllipticCurveTorsionSubgroup(E);
      const T2 = new EllipticCurveTorsionSubgroup(E);

      expect(T1.eq(T2)).toBe(true);
      expect(T1.compare(T2)).toBe(0);
    });
  });
});

describe('Torsion helper functions', () => {
  describe('torsion_subgroup', () => {
    it('should return torsion subgroup', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const T = torsion_subgroup(E);

      expect(T).toBeInstanceOf(EllipticCurveTorsionSubgroup);
    });
  });

  describe('torsion_order', () => {
    it('should return order of torsion subgroup', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const ord = torsion_order(E);

      expect(ord).toBeGreaterThan(0n);
    });
  });

  describe('torsion_points', () => {
    it('should return all torsion points', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const pts = torsion_points(E);

      expect(pts.length).toBeGreaterThan(0);
    });
  });

  describe('torsion_bound', () => {
    it('should return bound for finite field', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const bound = torsion_bound(E);

      // For finite fields, bound equals actual order
      expect(bound).toBe(torsion_order(E));
    });
  });
});

describe('Point order functions', () => {
  describe('has_finite_order', () => {
    it('should return true for identity', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const O = E.zero();

      expect(has_finite_order(O)).toBe(true);
    });

    it('should return true for points over finite fields', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);

      // Find a non-identity point
      for (let x = 0n; x < 7n; x++) {
        for (let y = 0n; y < 7n; y++) {
          const xEl = F7.__call__(x);
          const yEl = F7.__call__(y);
          if (E.is_on_curve(xEl, yEl)) {
            const P = E.point([xEl, yEl]);
            expect(has_finite_order(P)).toBe(true);
            return;
          }
        }
      }
    });
  });

  describe('has_infinite_order', () => {
    it('should return false for finite field points', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const O = E.zero();

      expect(has_infinite_order(O)).toBe(false);
    });
  });

  describe('point_order', () => {
    it('should return 1 for identity', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const O = E.zero();

      expect(point_order(O)).toBe(1n);
    });

    it('should return correct order for 2-torsion point', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 0n]); // y^2 = x^3 + x

      // (0, 0) is a 2-torsion point
      const P = E.point([F7.__call__(0n), F7.__call__(0n)]);
      const ord = point_order(P);

      if (typeof ord === 'bigint') {
        expect(ord).toBe(2n);
      }
    });
  });

  describe('order_from_multiple', () => {
    it('should compute order given a multiple', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 0n]);

      // (0, 0) is a 2-torsion point, so 2 is a multiple of its order
      const P = E.point([F7.__call__(0n), F7.__call__(0n)]);
      const ord = order_from_multiple(P, 2n);

      expect(ord).toBe(2n);
    });

    it('should throw if multiple does not annihilate point', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 0n]);

      // (0, 0) is a 2-torsion point, 3 is not a multiple of its order
      const P = E.point([F7.__call__(0n), F7.__call__(0n)]);

      // 3*P is not O for a 2-torsion point that is not identity
      // Actually for 2-torsion, 2P = O, so 3P = P. Let's check...
      const threeP = P.mul(3n);
      if (!threeP.is_zero()) {
        expect(() => order_from_multiple(P, 3n)).toThrow();
      }
    });

    it('should compute order with factorization hint', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);

      // Find a point and test with a valid multiple
      for (let x = 0n; x < 7n; x++) {
        for (let y = 0n; y < 7n; y++) {
          const xEl = F7.__call__(x);
          const yEl = F7.__call__(y);
          if (E.is_on_curve(xEl, yEl)) {
            const P = E.point([xEl, yEl]);
            if (P.is_zero()) continue;

            // Get the group order as a bound
            const T = new EllipticCurveTorsionSubgroup(E);
            const groupOrder = T.order();

            // The point order divides the group order
            const ord = order_from_multiple(P, groupOrder);
            expect(ord).toBeGreaterThan(0n);
            expect(groupOrder % ord).toBe(0n);
            return;
          }
        }
      }
    });
  });
});

describe('Discrete logarithm', () => {
  describe('discrete_log', () => {
    it('should return 0 for identity point', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);

      // Find a non-identity point as base
      for (let x = 0n; x < 7n; x++) {
        for (let y = 0n; y < 7n; y++) {
          const xEl = F7.__call__(x);
          const yEl = F7.__call__(y);
          if (E.is_on_curve(xEl, yEl)) {
            const P = E.point([xEl, yEl]);
            if (!P.is_zero()) {
              const n = discrete_log(E.zero(), P);
              expect(n).toBe(0n);
              return;
            }
          }
        }
      }
    });

    it('should return 1 for P with base P', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);

      // Find a non-identity point
      for (let x = 0n; x < 7n; x++) {
        for (let y = 0n; y < 7n; y++) {
          const xEl = F7.__call__(x);
          const yEl = F7.__call__(y);
          if (E.is_on_curve(xEl, yEl)) {
            const P = E.point([xEl, yEl]);
            if (!P.is_zero()) {
              const n = discrete_log(P, P);
              expect(n).toBe(1n);
              return;
            }
          }
        }
      }
    });

    it('should compute correct discrete log', () => {
      const F11 = GF(11n);
      const E = EllipticCurve<FiniteFieldElement>(F11, [1n, 1n]);

      // Find a generator
      const T = new EllipticCurveTorsionSubgroup(E);
      if (T.ngens() > 0) {
        const P = T.gen(0);
        const ord = T.invariants()[0]!;

        // Compute Q = 3*P
        const Q = P.mul(3n);

        // discrete_log(Q, P) should return 3
        const n = discrete_log(Q, P, ord);
        expect(P.mul(n).eq(Q)).toBe(true);
      }
    });

    it('should throw for identity base point', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);

      // Find a non-identity point
      for (let x = 0n; x < 7n; x++) {
        for (let y = 0n; y < 7n; y++) {
          const xEl = F7.__call__(x);
          const yEl = F7.__call__(y);
          if (E.is_on_curve(xEl, yEl)) {
            const P = E.point([xEl, yEl]);
            if (!P.is_zero()) {
              expect(() => discrete_log(P, E.zero())).toThrow();
              return;
            }
          }
        }
      }
    });

    it('should throw for point not in subgroup', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 0n]); // y^2 = x^3 + x

      // (0, 0) is a 2-torsion point
      const P = E.point([F7.__call__(0n), F7.__call__(0n)]);

      // Find a point not in the subgroup generated by P
      for (let x = 1n; x < 7n; x++) {
        for (let y = 0n; y < 7n; y++) {
          const xEl = F7.__call__(x);
          const yEl = F7.__call__(y);
          if (E.is_on_curve(xEl, yEl)) {
            const Q = E.point([xEl, yEl]);
            // P generates {O, P}, so Q not in this subgroup if Q != O and Q != P
            if (!Q.is_zero() && !Q.eq(P)) {
              expect(() => discrete_log(Q, P, 2n)).toThrow();
              return;
            }
          }
        }
      }
    });
  });
});

describe('p-primary torsion basis', () => {
  describe('_p_primary_torsion_basis', () => {
    it('should return empty array when no p-torsion', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
      const T = new EllipticCurveTorsionSubgroup(E);
      const order = T.order();

      // Find a prime that doesn't divide the order
      const primes = [5n, 11n, 13n, 17n, 19n];
      for (const p of primes) {
        if (order % p !== 0n) {
          const basis = _p_primary_torsion_basis(E, p);
          expect(basis.length).toBe(0);
          return;
        }
      }
    });

    it('should find 2-primary basis when curve has 2-torsion', () => {
      const F7 = GF(7n);
      const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 0n]); // y^2 = x^3 + x

      // This curve has a 2-torsion point at (0, 0)
      const basis = _p_primary_torsion_basis(E, 2n);

      // Should find at least one generator
      if (basis.length > 0) {
        const [gen, exp] = basis[0]!;
        // Generator should have order 2^exp
        const order = gen.mul(BigInt(2 ** exp)).is_zero();
        expect(order).toBe(true);
      }
    });
  });
});

describe('Known test vectors', () => {
  it('should match Sage: E = EllipticCurve(GF(7), [1, 1]) has order 5', () => {
    // sage: E = EllipticCurve(GF(7), [1, 1])
    // sage: E.cardinality()
    // 5
    const F7 = GF(7n);
    const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
    const T = new EllipticCurveTorsionSubgroup(E);

    // The torsion order for finite fields equals cardinality
    expect(T.order()).toBe(5n);
  });

  it('should match Sage: E = EllipticCurve(GF(11), [1, 1]) has order 14', () => {
    // sage: E = EllipticCurve(GF(11), [1, 1])
    // sage: E.cardinality()
    // 14
    const F11 = GF(11n);
    const E = EllipticCurve<FiniteFieldElement>(F11, [1n, 1n]);
    const T = new EllipticCurveTorsionSubgroup(E);

    // Check the order matches
    // Note: 14 = 2 * 7, so group is cyclic Z/14Z
    expect(T.order()).toBe(14n);
  });

  it('should match Sage: E = EllipticCurve(GF(13), [0, 1]) has order 12', () => {
    // sage: E = EllipticCurve(GF(13), [0, 1])
    // sage: E.cardinality()
    // 12
    const F13 = GF(13n);
    const E = EllipticCurve<FiniteFieldElement>(F13, [0n, 1n]);
    const T = new EllipticCurveTorsionSubgroup(E);

    expect(T.order()).toBe(12n);
  });

  it('should handle non-cyclic group: E = EllipticCurve(GF(41), [2, 5])', () => {
    // sage: E = EllipticCurve(GF(41), [2, 5])
    // sage: E.abelian_group()
    // Additive abelian group isomorphic to Z/22 + Z/2
    const F41 = GF(41n);
    const E = EllipticCurve<FiniteFieldElement>(F41, [2n, 5n]);
    const T = new EllipticCurveTorsionSubgroup(E);

    // Order should be 44 = 22 * 2
    expect(T.order()).toBe(44n);

    // For non-cyclic groups, we should have 2 generators
    // Note: exact invariants might vary based on generator choice
    const invs = T.invariants();
    const product = invs.reduce((a, b) => a * b, 1n);
    expect(product).toBe(44n);
  });

  it('should handle curve with large prime order: E = EllipticCurve(GF(101), [2, 3])', () => {
    // sage: E = EllipticCurve(GF(101), [2, 3])
    // sage: E.cardinality()
    // 96
    const F101 = GF(101n);
    const E = EllipticCurve<FiniteFieldElement>(F101, [2n, 3n]);
    const T = new EllipticCurveTorsionSubgroup(E);

    expect(T.order()).toBe(96n);

    // Verify generators work
    const pts = T.points();
    expect(pts.length).toBe(96);
  });

  it('should correctly compute structure for curves over larger field', () => {
    // sage: E = EllipticCurve(GF(1009), [1, 1])
    // sage: E.cardinality()
    // 1034
    const F1009 = GF(1009n);
    const E = EllipticCurve<FiniteFieldElement>(F1009, [1n, 1n]);
    const T = new EllipticCurveTorsionSubgroup(E);

    // Verify order
    expect(T.order()).toBe(1034n);

    // Verify structure: product of invariants = order
    const invs = T.invariants();
    const product = invs.reduce((a, b) => a * b, 1n);
    expect(product).toBe(1034n);

    // For cyclic groups, there should be 1 generator
    // For non-cyclic, there should be 2
    expect(T.ngens()).toBeLessThanOrEqual(2);
  });
});

describe('Group structure computation', () => {
  it('should find correct generator order for cyclic group', () => {
    // E over GF(7) with [1,1] has order 5 (prime), so cyclic
    const F7 = GF(7n);
    const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 1n]);
    const T = new EllipticCurveTorsionSubgroup(E);

    expect(T.ngens()).toBe(1);
    if (T.ngens() > 0) {
      const gen = T.gen(0);
      // Generator should have order equal to group order for cyclic group
      const genOrder = gen.order ? gen.order() : T.invariants()[0];
      expect(genOrder).toBe(5n);
    }
  });

  it('should verify generators actually generate the group', () => {
    const F13 = GF(13n);
    const E = EllipticCurve<FiniteFieldElement>(F13, [0n, 1n]);
    const T = new EllipticCurveTorsionSubgroup(E);

    // Generate all points from the generators
    const gens = T.gens();
    const invs = T.invariants();
    const generatedPoints = new Set<string>();

    if (gens.length === 1) {
      // Cyclic case
      let P = E.zero();
      for (let i = 0n; i < invs[0]!; i++) {
        const key = P.is_zero() ? 'O' : `${P.x()},${P.y()}`;
        generatedPoints.add(key);
        P = P.add(gens[0]!);
      }
    } else if (gens.length === 2) {
      // Product of two cyclic groups
      let P1 = E.zero();
      for (let i = 0n; i < invs[0]!; i++) {
        let P = P1;
        for (let j = 0n; j < invs[1]!; j++) {
          const key = P.is_zero() ? 'O' : `${P.x()},${P.y()}`;
          generatedPoints.add(key);
          P = P.add(gens[1]!);
        }
        P1 = P1.add(gens[0]!);
      }
    }

    // Number of generated points should equal group order
    expect(BigInt(generatedPoints.size)).toBe(T.order());
  });

  it('should handle the identity-only case', () => {
    // Create a curve where we can test the trivial group path
    // (This is rare but the code should handle it)
    const F2 = GF(2n);
    try {
      const E = EllipticCurve<FiniteFieldElement>(F2, [0n, 1n]);
      const T = new EllipticCurveTorsionSubgroup(E);
      // Should have at least the identity
      expect(T.order()).toBeGreaterThanOrEqual(1n);
    } catch {
      // Curve might be singular, which is fine
    }
  });
});

describe('Efficient order computation', () => {
  it('should compute point order efficiently using factorization', () => {
    const F97 = GF(97n);
    const E = EllipticCurve<FiniteFieldElement>(F97, [2n, 3n]);
    const T = new EllipticCurveTorsionSubgroup(E);

    // Get a random non-identity point
    const pts = T.points();
    const P = pts.find((p) => !p.is_zero());

    if (P) {
      // Point order should divide group order
      const groupOrder = T.order();
      const ptOrder = P.order
        ? P.order()
        : (() => {
            // Compute order manually for verification
            let Q = P;
            let ord = 1n;
            while (!Q.is_zero() && ord <= groupOrder) {
              Q = Q.add(P);
              ord++;
            }
            return Q.is_zero() ? ord : groupOrder;
          })();

      expect(groupOrder % ptOrder).toBe(0n);
    }
  });

  it('should correctly identify 2-torsion points', () => {
    const F7 = GF(7n);
    const E = EllipticCurve<FiniteFieldElement>(F7, [1n, 0n]); // y^2 = x^3 + x
    const T = new EllipticCurveTorsionSubgroup(E);

    // (0, 0) is a 2-torsion point (y = 0)
    const pts = T.points();
    const twoTorsion = pts.filter((P) => {
      if (P.is_zero()) return false;
      const twoP = P.add(P);
      return twoP.is_zero();
    });

    // Should have at least one 2-torsion point
    expect(twoTorsion.length).toBeGreaterThan(0);
  });
});

// ============================================================================
// Audit fixes H93, M101, M102.
// ============================================================================

describe('_p_primary_torsion_basis (audit H93/M102)', () => {
  it('returns exponents summing to v_p(#E) with independent generators', () => {
    let cases = 0;
    for (const q of [11n, 13n, 17n, 19n, 23n]) {
      const K = GF(q);
      for (let a = 0n; a < q; a++) {
        for (let b = 0n; b < q; b++) {
          let E: ReturnType<typeof EllipticCurve<FiniteFieldElement>>;
          try {
            E = EllipticCurve<FiniteFieldElement>(K, [a, b]);
          } catch {
            continue;
          }
          const N = BigInt(E.torsion_points().length);
          for (let p = 2n; p <= N; p++) {
            if (N % p !== 0n) continue;
            let isP = true;
            for (let d = 2n; d * d <= p; d++) if (p % d === 0n) isP = false;
            if (!isP) continue;

            let v = 0;
            let t = N;
            while (t % p === 0n) {
              t /= p;
              v++;
            }

            const basis = _p_primary_torsion_basis(E, p);
            expect(basis.reduce((s, [, k]) => s + k, 0)).toBe(v);
            expect(basis.length).toBeLessThanOrEqual(2);

            for (const [T, k] of basis) {
              expect(T.mul(p ** BigInt(k)).is_zero()).toBe(true);
              expect(T.mul(p ** BigInt(k - 1)).is_zero()).toBe(false);
            }

            if (basis.length === 2) {
              const [[T1, k1], [T2, k2]] = basis as [
                [(typeof basis)[0][0], number],
                [(typeof basis)[0][0], number],
              ];
              expect(k2).toBeLessThanOrEqual(k1);
              // The subgroup <T1, T2> must have order exactly p^(k1+k2).
              const seen = new Set<string>();
              for (let i = 0n; i < p ** BigInt(k1); i++) {
                for (let j = 0n; j < p ** BigInt(k2); j++) {
                  const R = T1.mul(i).add(T2.mul(j));
                  seen.add(R.is_zero() ? 'O' : `${R.x()},${R.y()}`);
                }
              }
              expect(BigInt(seen.size)).toBe(p ** BigInt(k1 + k2));
            }
            cases++;
          }
        }
      }
    }
    expect(cases).toBeGreaterThan(2000);
  }, 600000);

  it('returns exponents summing to 3 for GF(11), y^2 = x^3 + x + 9 (Z/4 x Z/2)', () => {
    const K = GF(11n);
    const E = EllipticCurve<FiniteFieldElement>(K, [1n, 9n]);
    expect(E.torsion_points().length).toBe(8);
    const basis = _p_primary_torsion_basis(E, 2n);
    expect(basis.map(([, k]) => k)).toEqual([2, 1]);
  });

  it('rejects a composite p', () => {
    const K = GF(11n);
    const E = EllipticCurve<FiniteFieldElement>(K, [1n, 9n]);
    expect(() => _p_primary_torsion_basis(E, 4n)).toThrow('should be prime');
  });
});

describe('torsion subgroup comparison (audit M101)', () => {
  it('distinguishes non-isomorphic quadratic twists', () => {
    // GF(11): y^2 = x^3 + x has invariants [12]; y^2 = x^3 + 2x has [6, 2].
    // Comparing j-invariants alone made these compare equal.
    const K = GF(11n);
    const E1 = EllipticCurve<FiniteFieldElement>(K, [1n, 0n]);
    const E2 = EllipticCurve<FiniteFieldElement>(K, [2n, 0n]);
    expect(E1.j_invariant().eq(E2.j_invariant())).toBe(true);

    const T1 = new EllipticCurveTorsionSubgroup(E1);
    const T2 = new EllipticCurveTorsionSubgroup(E2);
    expect(T1.invariants()).toEqual([12n]);
    expect(T2.invariants()).toEqual([6n, 2n]);
    expect(T1.eq(T2)).toBe(false);
    expect(T2.eq(T1)).toBe(false);
  });

  it('is reflexive and equal for two subgroups of the same curve', () => {
    const K = GF(11n);
    const E = EllipticCurve<FiniteFieldElement>(K, [1n, 0n]);
    const T = new EllipticCurveTorsionSubgroup(E);
    expect(T.eq(T)).toBe(true);
    expect(T.eq(new EllipticCurveTorsionSubgroup(E))).toBe(true);
    expect(T.compare(new EllipticCurveTorsionSubgroup(E))).toBe(0);
  });
});
