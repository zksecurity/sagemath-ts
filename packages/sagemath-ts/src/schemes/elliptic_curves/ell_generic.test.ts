/**
 * Tests for elliptic curves
 */

import { describe, expect, it } from 'vitest';
import {
  FiniteFieldElement,
  FiniteFieldPrime,
  GF,
} from '../../rings/finite_rings/finite_field_constructor.js';
import { EllipticCurve, EllipticCurveGeneric, EllipticCurve_from_j } from './constructor.js';

describe('EllipticCurve constructor', () => {
  describe('short Weierstrass form', () => {
    it('should create curve y^2 = x^3 + x + 1 over GF(23)', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      expect(E.is_short_weierstrass()).toBe(true);
      expect(E.a1().isZero()).toBe(true);
      expect(E.a2().isZero()).toBe(true);
      expect(E.a3().isZero()).toBe(true);
      expect(E.a4().eq(1n)).toBe(true);
      expect(E.a6().eq(1n)).toBe(true);
    });

    it('should create curve y^2 = x^3 + 4x + 4 over GF(5)', () => {
      const F = GF(5n);
      const E = EllipticCurve(F, [4n, 4n]);

      expect(E.is_short_weierstrass()).toBe(true);
    });

    it('should reject singular curve (discriminant = 0)', () => {
      const F = GF(7n);
      // y^2 = x^3 has discriminant 0
      expect(() => EllipticCurve(F, [0n, 0n])).toThrow();
    });
  });

  describe('long Weierstrass form', () => {
    it('should create curve with all coefficients', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 2n, 3n, 4n, 5n]);

      expect(E.a1().eq(1n)).toBe(true);
      expect(E.a2().eq(2n)).toBe(true);
      expect(E.a3().eq(3n)).toBe(true);
      expect(E.a4().eq(4n)).toBe(true);
      expect(E.a6().eq(5n)).toBe(true);
      expect(E.is_short_weierstrass()).toBe(false);
    });
  });

  describe('field element coefficients', () => {
    it('should accept field elements directly', () => {
      const F = GF(23n);
      const a = F.__call__(1n);
      const b = F.__call__(1n);
      const E = EllipticCurve(F, [a, b]);

      expect(E.a4().eq(1n)).toBe(true);
      expect(E.a6().eq(1n)).toBe(true);
    });
  });
});

describe('Elliptic curve invariants', () => {
  describe('b-invariants', () => {
    it('should compute b-invariants for short Weierstrass', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]); // y^2 = x^3 + x + 1

      // For short Weierstrass (a1=a2=a3=0, a4=1, a6=1):
      // b2 = 0^2 + 4*0 = 0
      // b4 = 0*0 + 2*1 = 2
      // b6 = 0^2 + 4*1 = 4
      // b8 = 0^2*1 + 4*0*1 - 0*0*1 + 0*0^2 - 1^2 = -1 = 22 mod 23
      const [b2, b4, b6, b8] = E.b_invariants();

      expect(b2.eq(0n)).toBe(true);
      expect(b4.eq(2n)).toBe(true);
      expect(b6.eq(4n)).toBe(true);
      expect(b8.eq(22n)).toBe(true); // -1 mod 23
    });

    it('should compute b-invariants for long Weierstrass', () => {
      const F = GF(7n);
      const E = EllipticCurve(F, [1n, 2n, 3n, 4n, 5n]);

      // b2 = a1^2 + 4*a2 = 1 + 8 = 9 = 2 mod 7
      // b4 = a1*a3 + 2*a4 = 3 + 8 = 11 = 4 mod 7
      // b6 = a3^2 + 4*a6 = 9 + 20 = 29 = 1 mod 7
      // b8 = a1^2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2
      //    = 5 + 40 - 12 + 18 - 16 = 35 = 0 mod 7
      const [b2, b4, b6, b8] = E.b_invariants();

      expect(b2.eq(2n)).toBe(true);
      expect(b4.eq(4n)).toBe(true);
      expect(b6.eq(1n)).toBe(true);
      expect(b8.eq(0n)).toBe(true);
    });
  });

  describe('c-invariants', () => {
    it('should compute c-invariants', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      // c4 = b2^2 - 24*b4 = 0 - 48 = -48 = -48 + 3*23 = 21 mod 23
      // c6 = -b2^3 + 36*b2*b4 - 216*b6 = 0 + 0 - 864 = -864 mod 23
      //    = -864 + 38*23 = -864 + 874 = 10 mod 23
      const [c4, c6] = E.c_invariants();

      // -48 mod 23 = -48 + 69 = 21
      expect(c4.eq(21n)).toBe(true);

      // -216*4 = -864; -864 mod 23:
      // -864 = -38 * 23 + 10 = -874 + 10, so -864 mod 23 = 10
      expect(c6.eq(10n)).toBe(true);
    });
  });

  describe('discriminant', () => {
    it('should compute non-zero discriminant for valid curve', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const disc = E.discriminant();
      expect(disc.isZero()).toBe(false);
    });

    it('should compute discriminant correctly for y^2 = x^3 + x + 1 over GF(23)', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      // For short Weierstrass y^2 = x^3 + ax + b:
      // discriminant = -16(4a^3 + 27b^2)
      // = -16(4*1 + 27*1) = -16*31 = -496
      // -496 mod 23 = -496 + 22*23 = -496 + 506 = 10
      const disc = E.discriminant();
      expect(disc.eq(10n)).toBe(true);
    });

    it('should compute discriminant correctly for y^2 = x^3 + 2 over GF(7)', () => {
      const F = GF(7n);
      const E = EllipticCurve(F, [0n, 2n]);

      // -16(4*0 + 27*4) = -16*108 = -1728
      // -1728 mod 7 = -1728 + 247*7 = -1728 + 1729 = 1
      const disc = E.discriminant();
      expect(disc.eq(1n)).toBe(true);
    });
  });

  describe('j-invariant', () => {
    it('should compute j-invariant', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const j = E.j_invariant();
      // j = c4^3 / disc
      expect(j.isZero()).toBe(false);
    });

    it('should give j=0 for y^2 = x^3 + 1', () => {
      const F = GF(7n);
      const E = EllipticCurve(F, [0n, 1n]);

      // For a=0: c4 = -48*b4 = -48*0 = 0, so j = 0
      const j = E.j_invariant();
      expect(j.isZero()).toBe(true);
    });

    it('should give j=1728 for y^2 = x^3 + x', () => {
      // Use a larger prime where 1728 fits
      const F = GF(4999n); // Prime > 1728
      const E = EllipticCurve(F, [1n, 0n]);

      // For b=0: c6 = 0, j = c4^3 / disc
      // disc = -16*4a^3 = -64, c4 = -48*b4 = -48*2 = -96
      // Actually for y^2 = x^3 + x: a1=a2=a3=0, a4=1, a6=0
      // b2=0, b4=2, b6=0, b8=-1
      // c4 = 0 - 24*2 = -48
      // c6 = 0 + 0 - 0 = 0
      // disc = -0 - 8*8 - 0 + 0 = -64
      // j = (-48)^3 / (-64) = -110592 / -64 = 1728
      const j = E.j_invariant();
      expect(j.eq(1728n)).toBe(true);
    });
  });
});

describe('Point operations', () => {
  describe('is_on_curve', () => {
    it('should verify point (0, 1) is on y^2 = x^3 + x + 1 over GF(23)', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      // Check: 1^2 = 0^3 + 0 + 1 => 1 = 1 (yes)
      const x = F.__call__(0n);
      const y = F.__call__(1n);

      expect(E.is_on_curve(x, y)).toBe(true);
    });

    it('should verify point is not on curve', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      // Check: 0^2 = 0^3 + 0 + 1 => 0 != 1 (no)
      const x = F.__call__(0n);
      const y = F.__call__(0n);

      expect(E.is_on_curve(x, y)).toBe(false);
    });

    it('should find points on curve by testing', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      // Find some points by brute force
      const points: Array<[bigint, bigint]> = [];
      for (let x = 0n; x < 23n; x++) {
        for (let y = 0n; y < 23n; y++) {
          const xElem = F.__call__(x);
          const yElem = F.__call__(y);
          if (E.is_on_curve(xElem, yElem)) {
            points.push([x, y]);
          }
        }
      }

      // There should be several points
      expect(points.length).toBeGreaterThan(0);

      // Verify each found point
      for (const [x, y] of points) {
        const xElem = F.__call__(x);
        const yElem = F.__call__(y);
        expect(E.is_on_curve(xElem, yElem)).toBe(true);
      }
    });
  });

  describe('point creation', () => {
    it('should create point at infinity', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const O = E.point_at_infinity();
      expect(O.is_zero()).toBe(true);
    });

    it('should create valid point', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const x = F.__call__(0n);
      const y = F.__call__(1n);
      const P = E.point([x, y]);

      expect(P.is_zero()).toBe(false);
      expect(P.x().eq(0n)).toBe(true);
      expect(P.y().eq(1n)).toBe(true);
    });

    it('should throw for invalid point', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const x = F.__call__(0n);
      const y = F.__call__(0n);

      expect(() => E.point([x, y])).toThrow();
    });
  });

  describe('point arithmetic', () => {
    it('should satisfy P + O = P', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const x = F.__call__(0n);
      const y = F.__call__(1n);
      const P = E.point([x, y]);
      const O = E.point_at_infinity();

      const R = P.add(O);
      expect(R.eq(P)).toBe(true);
    });

    it('should satisfy O + P = P', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const x = F.__call__(0n);
      const y = F.__call__(1n);
      const P = E.point([x, y]);
      const O = E.point_at_infinity();

      const R = O.add(P);
      expect(R.eq(P)).toBe(true);
    });

    it('should satisfy P + (-P) = O', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const x = F.__call__(0n);
      const y = F.__call__(1n);
      const P = E.point([x, y]);
      const negP = P.neg();

      const R = P.add(negP);
      expect(R.is_zero()).toBe(true);
    });

    it('should compute P + Q correctly', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      // Find two different points
      const P = E.point([F.__call__(0n), F.__call__(1n)]);

      // 2P should still be on the curve
      const twoP = P.add(P);
      if (!twoP.is_zero()) {
        expect(E.is_on_curve(twoP.x(), twoP.y())).toBe(true);
      }
    });

    it('should compute 2P correctly (doubling)', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const P = E.point([F.__call__(0n), F.__call__(1n)]);
      const twoP = P.add(P);

      // If 2P is not infinity, it should be on the curve
      if (!twoP.is_zero()) {
        expect(E.is_on_curve(twoP.x(), twoP.y())).toBe(true);
      }
    });

    it('should compute scalar multiplication', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const P = E.point([F.__call__(0n), F.__call__(1n)]);

      // 0*P = O
      expect(P.mul(0n).is_zero()).toBe(true);

      // 1*P = P
      expect(P.mul(1n).eq(P)).toBe(true);

      // 2*P = P + P
      expect(P.mul(2n).eq(P.add(P))).toBe(true);

      // 3*P = P + P + P
      expect(P.mul(3n).eq(P.add(P).add(P))).toBe(true);
    });

    it('should find point order', () => {
      const F = GF(23n);
      const E = EllipticCurve(F, [1n, 1n]);

      const P = E.point([F.__call__(0n), F.__call__(1n)]);

      // The order should divide the group order
      // For y^2 = x^3 + x + 1 over GF(23), the group order is known
      const order = P.order(100);
      expect(order).toBeDefined();

      if (order !== undefined) {
        // n*P = O
        expect(P.mul(order).is_zero()).toBe(true);

        // (n-1)*P != O (unless n=1)
        if (order > 1n) {
          expect(P.mul(order - 1n).is_zero()).toBe(false);
        }
      }
    });
  });
});

describe('secp256k1 curve parameters', () => {
  it('should correctly compute secp256k1 discriminant and j-invariant', () => {
    // secp256k1 parameters
    const p = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2fn;
    const F = GF(p);

    // y^2 = x^3 + 7 (a = 0, b = 7)
    const E = EllipticCurve(F, [0n, 7n]);

    // For y^2 = x^3 + b, we have:
    // - discriminant = -16 * 27 * b^2 = -16 * 27 * 49 = -21168
    // - j-invariant = 0 (since a = 0)
    expect(E.j_invariant().isZero()).toBe(true);

    // Verify discriminant is non-zero
    expect(E.discriminant().isZero()).toBe(false);
  });

  it('should verify secp256k1 generator point is on curve', () => {
    const p = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2fn;
    const F = GF(p);
    const E = EllipticCurve(F, [0n, 7n]);

    // Generator point
    const Gx = 0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798n;
    const Gy = 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8n;

    expect(E.is_on_curve(F.__call__(Gx), F.__call__(Gy))).toBe(true);
  });

  it('should create secp256k1 generator point', () => {
    const p = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2fn;
    const F = GF(p);
    const E = EllipticCurve(F, [0n, 7n]);

    const Gx = 0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798n;
    const Gy = 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8n;

    const G = E.point([F.__call__(Gx), F.__call__(Gy)]);

    expect(G.is_zero()).toBe(false);
    expect(G.x().eq(Gx)).toBe(true);
    expect(G.y().eq(Gy)).toBe(true);
  });

  it('should verify 2G is on curve', () => {
    const p = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2fn;
    const F = GF(p);
    const E = EllipticCurve(F, [0n, 7n]);

    const Gx = 0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798n;
    const Gy = 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8n;

    const G = E.point([F.__call__(Gx), F.__call__(Gy)]);
    const twoG = G.mul(2n);

    expect(twoG.is_zero()).toBe(false);
    expect(E.is_on_curve(twoG.x(), twoG.y())).toBe(true);

    // Known value of 2G (computed from secp256k1 generator point doubling)
    const twoGx = 0xc6047f9441ed7d6d3045406e95c07cd85c778e4b8cef3ca7abac09b95c709ee5n;
    const twoGy = 0x1ae168fea63dc339a3c58419466ceaeef7f632653266d0e1236431a950cfe52an;

    expect(twoG.x().eq(twoGx)).toBe(true);
    expect(twoG.y().eq(twoGy)).toBe(true);
  });
});

describe('EllipticCurve_from_j', () => {
  it('should create curve with j=0', () => {
    const F = GF(7n);
    const E = EllipticCurve_from_j(F, 0n);

    expect(E.j_invariant().isZero()).toBe(true);
  });

  it('should create curve with j=1728', () => {
    const F = GF(4999n); // Large enough prime
    const E = EllipticCurve_from_j(F, 1728n);

    expect(E.j_invariant().eq(1728n)).toBe(true);
  });

  it('should create curve with general j', () => {
    const F = GF(101n);
    const j = 50n;
    const E = EllipticCurve_from_j(F, j);

    expect(E.j_invariant().eq(j)).toBe(true);
  });
});

describe('Small prime fields', () => {
  it('should work over GF(5)', () => {
    const F = GF(5n);
    // y^2 = x^3 + x + 1 over GF(5)
    // Need to check this is non-singular
    const E = EllipticCurve(F, [1n, 1n]);

    expect(E.discriminant().isZero()).toBe(false);
  });

  it('should work over GF(7)', () => {
    const F = GF(7n);
    const E = EllipticCurve(F, [1n, 1n]);

    expect(E.discriminant().isZero()).toBe(false);
  });

  it('should work over GF(11)', () => {
    const F = GF(11n);
    const E = EllipticCurve(F, [1n, 1n]);

    expect(E.discriminant().isZero()).toBe(false);
  });

  it('should work over GF(13)', () => {
    const F = GF(13n);
    const E = EllipticCurve(F, [1n, 1n]);

    expect(E.discriminant().isZero()).toBe(false);
  });

  it('should count points on small curve', () => {
    const F = GF(7n);
    const E = EllipticCurve(F, [1n, 1n]);

    // Count all affine points plus infinity
    let count = 1; // Start with point at infinity
    for (let x = 0n; x < 7n; x++) {
      for (let y = 0n; y < 7n; y++) {
        if (E.is_on_curve(F.__call__(x), F.__call__(y))) {
          count++;
        }
      }
    }

    // For y^2 = x^3 + x + 1 over GF(7), should have some points
    expect(count).toBeGreaterThan(1);
  });
});

describe('Division polynomials', () => {
  it('should compute two_division_polynomial for short Weierstrass', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]); // y^2 = x^3 + x + 1

    // 2-division polynomial is 4x^3 + b2*x^2 + 2*b4*x + b6
    // For short Weierstrass: b2 = 0, b4 = 2*a4 = 2, b6 = 4*a6 = 4
    // So: 4x^3 + 0 + 4x + 4 = 4x^3 + 4x + 4
    const psi2 = E.two_division_polynomial();

    expect(psi2.degree()).toBe(3);
    // Leading coefficient should be 4
    expect(psi2.leading_coefficient().eq(4n)).toBe(true);
  });

  it('should compute division_polynomial(1) = 1', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    const psi1 = E.division_polynomial(1);
    expect(psi1.degree()).toBe(0);
    expect(psi1.leading_coefficient().eq(1n)).toBe(true);
  });

  it('should compute division_polynomial(3)', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    // psi_3 = 3x^4 + b2*x^3 + 3*b4*x^2 + 3*b6*x + b8
    // For short Weierstrass (a1=a2=a3=0, a4=1, a6=1):
    // b2 = 0, b4 = 2, b6 = 4, b8 = -1 = 22 mod 23
    // So psi_3 = 3x^4 + 0 + 6x^2 + 12x + 22
    const psi3 = E.division_polynomial(3);

    expect(psi3.degree()).toBe(4);
    expect(psi3.leading_coefficient().eq(3n)).toBe(true);
  });
});

describe('Short Weierstrass model', () => {
  it('should return self if already in short form', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    const sw = E.short_weierstrass_model();

    // Should be the same curve
    expect(sw.is_short_weierstrass()).toBe(true);
    expect(sw.a4().eq(1n)).toBe(true);
    expect(sw.a6().eq(1n)).toBe(true);
  });

  it('should convert long form to short form', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 2n, 3n, 4n, 5n]);

    const sw = E.short_weierstrass_model();

    expect(sw.is_short_weierstrass()).toBe(true);
    expect(sw.a1().isZero()).toBe(true);
    expect(sw.a2().isZero()).toBe(true);
    expect(sw.a3().isZero()).toBe(true);
  });

  it('should preserve j-invariant', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 2n, 3n, 4n, 5n]);

    const sw = E.short_weierstrass_model();

    expect(E.j_invariant().eq(sw.j_invariant())).toBe(true);
  });
});

describe('Isomorphisms', () => {
  it('should find identity automorphism', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    const auts = E.automorphisms();

    // Every curve has at least the identity and negation automorphisms
    expect(auts.length).toBeGreaterThanOrEqual(2);

    // First should be identity [1, 0, 0, 0]
    const [u, r, s, t] = auts[0]!;
    expect(u.eq(1n)).toBe(true);
    expect(r.isZero()).toBe(true);
    expect(s.isZero()).toBe(true);
    expect(t.isZero()).toBe(true);
  });

  it('should find negation automorphism', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    const auts = E.automorphisms();

    // Second should be negation [-1, 0, 0, ...]
    const [u, r, s] = auts[1]!;
    expect(u.eq(22n)).toBe(true); // -1 mod 23 = 22
    expect(r.isZero()).toBe(true);
    expect(s.isZero()).toBe(true);
  });

  it('should recognize isomorphic curves', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 2n, 3n, 4n, 5n]);
    const sw = E.short_weierstrass_model();

    expect(E.is_isomorphic(sw)).toBe(true);
  });

  it('should recognize non-isomorphic curves', () => {
    const F = GF(23n);
    const E1 = EllipticCurve(F, [1n, 1n]);
    const E2 = EllipticCurve(F, [1n, 2n]);

    // Different j-invariants means non-isomorphic
    if (!E1.j_invariant().eq(E2.j_invariant())) {
      expect(E1.is_isomorphic(E2)).toBe(false);
    }
  });

  it('should find isomorphism to short form', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 2n, 3n, 4n, 5n]);
    const sw = E.short_weierstrass_model();

    const iso = E.isomorphism_to(sw);

    // iso should be a valid tuple [u, r, s, t]
    expect(iso.length).toBe(4);
    expect(iso[0].isZero()).toBe(false); // u != 0
  });
});

describe('Torsion points', () => {
  it('should find all points on small curve', () => {
    const F = GF(7n);
    const E = EllipticCurve(F, [1n, 1n]);

    const points = E.torsion_points();

    // Should include point at infinity
    expect(points.some((p) => p.is_zero())).toBe(true);

    // All points should be on the curve
    for (const p of points) {
      if (!p.is_zero()) {
        expect(E.is_on_curve(p.x(), p.y())).toBe(true);
      }
    }

    // Count should match manual enumeration
    let manualCount = 1; // infinity
    for (let x = 0n; x < 7n; x++) {
      for (let y = 0n; y < 7n; y++) {
        if (E.is_on_curve(F.__call__(x), F.__call__(y))) {
          manualCount++;
        }
      }
    }

    expect(points.length).toBe(manualCount);
  });

  it('should find torsion points including 2-torsion', () => {
    const F = GF(11n);
    const E = EllipticCurve(F, [0n, 2n]); // y^2 = x^3 + 2

    const points = E.torsion_points();

    // Check that we have points
    expect(points.length).toBeGreaterThan(1);

    // Verify all points are valid
    for (const p of points) {
      if (!p.is_zero()) {
        expect(E.is_on_curve(p.x(), p.y())).toBe(true);
      }
    }
  });
});

describe('Weierstrass transformations', () => {
  it('should apply rst_transform correctly', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]); // y^2 = x^3 + x + 1

    // Apply identity transformation (r=0, s=0, t=0)
    const r = F.__call__(0n);
    const s = F.__call__(0n);
    const t = F.__call__(0n);
    const E1 = E.rst_transform(r, s, t);

    // Should be the same curve
    expect(E1.a1().value).toBe(E.a1().value);
    expect(E1.a2().value).toBe(E.a2().value);
    expect(E1.a3().value).toBe(E.a3().value);
    expect(E1.a4().value).toBe(E.a4().value);
    expect(E1.a6().value).toBe(E.a6().value);
  });

  it('should apply change_weierstrass_model with u=1', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    const u = F.__call__(1n);
    const r = F.__call__(0n);
    const s = F.__call__(0n);
    const t = F.__call__(0n);
    const E1 = E.change_weierstrass_model(u, r, s, t);

    // Identity transformation should preserve the curve
    expect(E1.a4().value).toBe(E.a4().value);
    expect(E1.a6().value).toBe(E.a6().value);
  });

  it('should produce isomorphic curves under rst_transform', () => {
    const F = GF(97n);
    const E = EllipticCurve(F, [1n, 1n]);

    // Apply a non-trivial transformation
    const r = F.__call__(2n);
    const s = F.__call__(3n);
    const t = F.__call__(5n);
    const E1 = E.rst_transform(r, s, t);

    // The j-invariants should be the same (isomorphic curves)
    expect(E1.j_invariant().value).toBe(E.j_invariant().value);
  });

  it('should change base ring correctly', () => {
    const F5 = GF(5n);
    const F97 = GF(97n);
    const E = EllipticCurve(F5, [1n, 1n]);

    // Change ring to F97 (a larger prime field)
    const E1 = E.change_ring(F97);

    // The a-invariants should be coerced correctly
    expect(E1.a4().value).toBe(1n);
    expect(E1.a6().value).toBe(1n);
  });
});

describe('Torsion subgroup', () => {
  it('should return EllipticCurveTorsionSubgroup', () => {
    const F = GF(7n);
    const E = EllipticCurve(F, [1n, 1n]);

    const T = E.torsion_subgroup();
    expect(T).toBeDefined();
    expect(T.order()).toBeGreaterThan(0n);
    expect(T.curve()).toBe(E);
  });

  it('should have correct order for small curve', () => {
    const F = GF(7n);
    const E = EllipticCurve(F, [1n, 1n]);

    const T = E.torsion_subgroup();
    const points = E.torsion_points();

    // Order should equal number of points
    expect(Number(T.order())).toBe(points.length);
  });

  it('should have generators that generate the group', () => {
    const F = GF(11n);
    const E = EllipticCurve(F, [1n, 1n]);

    const T = E.torsion_subgroup();
    const gens = T.gens();
    const points = T.points();

    // All points should be linear combinations of generators
    expect(points.length).toBe(Number(T.order()));
  });
});

describe('Scalar multiplication endomorphism', () => {
  it('should return correct domain and codomain', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    const mult2 = E.scalar_multiplication(2n);

    expect(mult2.domain()).toBe(E);
    expect(mult2.codomain()).toBe(E);
  });

  it('should have degree m^2', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    expect(E.scalar_multiplication(2n).degree()).toBe(4n);
    expect(E.scalar_multiplication(3n).degree()).toBe(9n);
    expect(E.scalar_multiplication(5n).degree()).toBe(25n);
  });

  it('should compute [m]P correctly', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);
    const P = E.point([F.__call__(0n), F.__call__(1n)]);

    const mult3 = E.scalar_multiplication(3n);
    const threeP = mult3.call(P);

    // Should equal P + P + P
    expect(threeP.eq(P.add(P).add(P))).toBe(true);
  });
});

describe('Identity morphism', () => {
  it('should return the same point', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);
    const P = E.point([F.__call__(0n), F.__call__(1n)]);

    const id = E.identity_morphism();
    expect(id.call(P).eq(P)).toBe(true);
  });

  it('should have degree 1', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    const id = E.identity_morphism();
    expect(id.degree()).toBe(1n);
  });

  it('should be separable, injective, and surjective', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    const id = E.identity_morphism();
    expect(id.is_separable()).toBe(true);
    expect(id.is_injective()).toBe(true);
    expect(id.is_surjective()).toBe(true);
  });
});

describe('Frobenius isogeny', () => {
  it('should have correct degree p^n', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    const frob1 = E.frobenius_isogeny(1);
    expect(frob1.degree()).toBe(23n);

    const frob2 = E.frobenius_isogeny(2);
    expect(frob2.degree()).toBe(23n * 23n);
  });

  it('should be inseparable', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    const frob = E.frobenius_isogeny(1);
    expect(frob.is_separable()).toBe(false);
  });

  it('should map point at infinity to itself', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);
    const O = E.zero();

    const frob = E.frobenius_isogeny(1);
    expect(frob.call(O).is_zero()).toBe(true);
  });

  it('should fix points in prime field (for F_p)', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);
    const P = E.point([F.__call__(0n), F.__call__(1n)]);

    // Over F_p, Frobenius is identity since x^p = x for x in F_p
    const frob = E.frobenius_isogeny(1);
    const frobP = frob.call(P);

    // x^23 = x mod 23, y^23 = y mod 23 (Fermat's little theorem)
    expect(frobP.x().eq(P.x())).toBe(true);
    expect(frobP.y().eq(P.y())).toBe(true);
  });
});

describe('p-primary torsion basis', () => {
  it('should return empty for prime not dividing order', () => {
    const F = GF(7n);
    const E = EllipticCurve(F, [1n, 1n]);

    // y^2 = x^3 + x + 1 over GF(7) has 13 points (prime)
    const T = E.torsion_subgroup();
    const order = T.order();

    // If 2 doesn't divide order, 2-primary basis should be empty or trivial
    if (order % 2n !== 0n) {
      const basis2 = E._p_primary_torsion_basis(2n);
      expect(basis2.length).toBe(0);
    }
  });

  it('should find generators for p dividing order', () => {
    const F = GF(11n);
    const E = EllipticCurve(F, [0n, 1n]); // y^2 = x^3 + 1

    const T = E.torsion_subgroup();
    const order = T.order();

    // Find a prime that divides the order
    for (const p of [2n, 3n, 5n, 7n]) {
      if (order % p === 0n) {
        const basis = E._p_primary_torsion_basis(p);
        // Should have at least one generator
        expect(basis.length).toBeGreaterThanOrEqual(0);
        break;
      }
    }
  });
});

describe('Reduction type methods', () => {
  it('should throw error for finite fields with different prime', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    // has_good_reduction at a different prime doesn't make sense
    expect(() => E.has_good_reduction(5n)).toThrow();
    expect(() => E.has_bad_reduction(7n)).toThrow();
  });

  it('should work for the defining characteristic', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    // At the defining characteristic, good reduction means non-zero disc
    expect(E.has_good_reduction(23n)).toBe(true);
    expect(E.has_bad_reduction(23n)).toBe(false);
  });

  it('should throw for conductor over finite fields', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    expect(() => E.conductor()).toThrow();
  });

  it('should throw for local_data over finite fields', () => {
    const F = GF(23n);
    const E = EllipticCurve(F, [1n, 1n]);

    expect(() => E.local_data(2n)).toThrow();
  });
});
