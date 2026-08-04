/**
 * @module parigp-ts/elliptic/points.test
 * @description Tests for elliptic curve point operations
 */

import { describe, expect, it } from 'bun:test';
import {
  type EllipticPoint,
  FpE_isoncurve,
  FpE_to_FpJ,
  FpJ_is_inf,
  FpJ_to_FpE,
  type JacobianPoint,
  type ShortWeierstrassCurve,
  ell_is_inf,
  ellinf,
  ellinf_FpJ,
  ellordinate,
  mkpoint,
  random_FpE,
} from './points.js';

describe('ellinf', () => {
  it('returns a point at infinity', () => {
    const inf = ellinf();
    expect(inf.isInfinity).toBe(true);
  });

  it('is recognized by ell_is_inf', () => {
    const inf = ellinf();
    expect(ell_is_inf(inf)).toBe(true);
  });
});

describe('mkpoint', () => {
  it('creates a finite point', () => {
    const P = mkpoint(3n, 5n);
    expect(P.isInfinity).toBe(false);
    if (!P.isInfinity) {
      expect(P.x).toBe(3n);
      expect(P.y).toBe(5n);
    }
  });

  it('is not recognized as infinity', () => {
    const P = mkpoint(0n, 0n);
    expect(ell_is_inf(P)).toBe(false);
  });
});

describe('ellordinate', () => {
  // Curve: y^2 = x^3 + x + 1 over F_23
  const E1: ShortWeierstrassCurve = { a4: 1n, a6: 1n, p: 23n };

  it('finds y-coordinates for x=0 on y^2 = x^3 + x + 1 mod 23', () => {
    // x = 0: y^2 = 0 + 0 + 1 = 1, so y = 1 or y = 22
    const ys = ellordinate(E1, 0n);
    expect(ys).toEqual([1n, 22n]);
  });

  it('finds y-coordinates for x=1 on y^2 = x^3 + x + 1 mod 23', () => {
    // x = 1: y^2 = 1 + 1 + 1 = 3
    // Need to check if 3 is a QR mod 23
    // 3^11 mod 23 = 3^11 mod 23
    // 3^2 = 9, 3^4 = 81 mod 23 = 12, 3^8 = 144 mod 23 = 6
    // 3^11 = 3^8 * 3^2 * 3 = 6 * 9 * 3 = 162 mod 23 = 1
    // So 3 is a QR mod 23
    const ys = ellordinate(E1, 1n);
    expect(ys.length).toBe(2);
    // Verify both are valid square roots of 3 mod 23
    for (const y of ys) {
      expect((y * y) % 23n).toBe(3n);
    }
  });

  it('returns empty array when no y exists', () => {
    // Find an x where x^3 + x + 1 is not a QR mod 23
    // Try x = 2: 8 + 2 + 1 = 11
    // 11^11 mod 23: 11^2 = 121 mod 23 = 6, 11^4 = 36 mod 23 = 13
    // 11^8 = 169 mod 23 = 8, 11^11 = 8 * 6 * 11 = 528 mod 23 = 22 = -1
    // So 11 is NOT a QR mod 23
    const ys = ellordinate(E1, 2n);
    expect(ys).toEqual([]);
  });

  it('handles y = 0 case', () => {
    // Curve: y^2 = x^3 over F_23 (a4 = 0, a6 = 0)
    const E: ShortWeierstrassCurve = { a4: 0n, a6: 0n, p: 23n };
    const ys = ellordinate(E, 0n);
    expect(ys).toEqual([0n]);
  });

  // Test with secp256k1 parameters (small test)
  it('works with larger primes', () => {
    // Curve: y^2 = x^3 + 7 over a small prime
    const E: ShortWeierstrassCurve = { a4: 0n, a6: 7n, p: 1009n };

    // x = 1: y^2 = 1 + 7 = 8
    const ys = ellordinate(E, 1n);
    // Verify results if any
    for (const y of ys) {
      expect((y * y) % 1009n).toBe(8n);
    }
  });

  it('handles negative x values by reducing mod p', () => {
    const ys1 = ellordinate(E1, -1n);
    const ys2 = ellordinate(E1, 22n); // -1 mod 23 = 22
    expect(ys1).toEqual(ys2);
  });

  /**
   * PARI returns [(d-b)/2, (d-b)/2 - d] with d = Fp_sqrt(b^2+4a) — it does
   * NOT sort the two ordinates.  Oracle (PARI 2.x through Sage 10.3):
   *   E = ellinit([1,1],23);
   *   E.ellordinate(1) -> [Mod(16,23), Mod(7,23)]
   *   E.ellordinate(3) -> [Mod(13,23), Mod(10,23)]
   *   E.ellordinate(0) -> [Mod(1,23),  Mod(22,23)]
   */
  it('returns PARI order (descending when d/2 > d/2 - d)', () => {
    expect(ellordinate(E1, 1n)).toEqual([16n, 7n]);
    expect(ellordinate(E1, 3n)).toEqual([13n, 10n]);
    expect(ellordinate(E1, 7n)).toEqual([12n, 11n]);
    // ascending case, for contrast
    expect(ellordinate(E1, 0n)).toEqual([1n, 22n]);
  });
});

describe('random_FpE', () => {
  // Curve: y^2 = x^3 + x + 1 over F_23
  const E1: ShortWeierstrassCurve = { a4: 1n, a6: 1n, p: 23n };

  it('generates a point on the curve', () => {
    const P = random_FpE(E1);
    expect(ell_is_inf(P)).toBe(false);
    expect(FpE_isoncurve(E1, P)).toBe(true);
  });

  it('generates different points on repeated calls', () => {
    const points: EllipticPoint[] = [];
    const maxTries = 100;
    let foundDifferent = false;

    for (let i = 0; i < maxTries; i++) {
      const P = random_FpE(E1);
      expect(FpE_isoncurve(E1, P)).toBe(true);

      // Check if we got a different point
      if (points.length > 0) {
        const last = points[points.length - 1]!;
        if (!ell_is_inf(P) && !ell_is_inf(last)) {
          if (P.x !== last.x || P.y !== last.y) {
            foundDifferent = true;
          }
        }
      }
      points.push(P);
    }

    // It's statistically almost certain we'll get different points
    expect(foundDifferent).toBe(true);
  });

  it('works with y^2 = x^3 + 7 (secp256k1-like)', () => {
    const E: ShortWeierstrassCurve = { a4: 0n, a6: 7n, p: 1009n };

    for (let i = 0; i < 10; i++) {
      const P = random_FpE(E);
      expect(FpE_isoncurve(E, P)).toBe(true);
    }
  });

  it('handles curves with y = 0 points correctly', () => {
    // Curve: y^2 = x^3 over F_23
    // This has singular point at origin, but random_FpE should avoid it
    // Actually y^2 = x^3 has a cusp singularity at (0,0)
    // For testing, use a non-singular curve
    const E: ShortWeierstrassCurve = { a4: 1n, a6: 0n, p: 23n };

    for (let i = 0; i < 10; i++) {
      const P = random_FpE(E);
      expect(FpE_isoncurve(E, P)).toBe(true);
    }
  });

  /**
   * PARI's guard `(!signe(rhs) && !signe(3x^2+a4))` rejects the singular
   * point of a singular curve.  y^2 = x^3 over F_11 has a cusp at (0,0);
   * random_FpE must never return it.
   */
  it('never returns the singular point of a singular curve', () => {
    const E: ShortWeierstrassCurve = { a4: 0n, a6: 0n, p: 11n };

    for (let i = 0; i < 3000; i++) {
      const P = random_FpE(E);
      expect(ell_is_inf(P)).toBe(false);
      if (!P.isInfinity) {
        expect(P.x === 0n && P.y === 0n).toBe(false);
      }
    }
  });

  /**
   * PARI returns Fp_sqrt(rhs, p), the canonical (smallest) square root; the
   * sign is not randomized.
   */
  it('returns the canonical (smallest) square root, like PARI', () => {
    const E: ShortWeierstrassCurve = { a4: 1n, a6: 1n, p: 23n };
    for (let i = 0; i < 200; i++) {
      const P = random_FpE(E);
      expect(P.isInfinity).toBe(false);
      if (!P.isInfinity) {
        expect(P.y <= 23n - P.y).toBe(true);
      }
    }
  });
});

describe('FpE_to_FpJ', () => {
  it('converts point at infinity', () => {
    const inf = ellinf();
    const J = FpE_to_FpJ(inf);
    expect(J.X).toBe(1n);
    expect(J.Y).toBe(1n);
    expect(J.Z).toBe(0n);
  });

  it('converts finite point', () => {
    const P = mkpoint(3n, 5n);
    const J = FpE_to_FpJ(P);
    expect(J.X).toBe(3n);
    expect(J.Y).toBe(5n);
    expect(J.Z).toBe(1n);
  });
});

describe('FpJ_to_FpE', () => {
  const p = 23n;

  it('converts Jacobian infinity to affine infinity', () => {
    const J: JacobianPoint = { X: 1n, Y: 1n, Z: 0n };
    const P = FpJ_to_FpE(J, p);
    expect(ell_is_inf(P)).toBe(true);
  });

  it('converts Jacobian point with Z=1 to affine', () => {
    const J: JacobianPoint = { X: 3n, Y: 5n, Z: 1n };
    const P = FpJ_to_FpE(J, p);
    expect(ell_is_inf(P)).toBe(false);
    if (!P.isInfinity) {
      expect(P.x).toBe(3n);
      expect(P.y).toBe(5n);
    }
  });

  it('correctly scales Jacobian coordinates', () => {
    // Jacobian (X, Y, Z) = affine (X/Z^2, Y/Z^3)
    // Let's use Z = 2, so affine x = X/4, y = Y/8
    // To get affine (3, 5), we need X = 3*4 = 12, Y = 5*8 = 40 mod 23 = 17
    const J: JacobianPoint = { X: 12n, Y: 17n, Z: 2n };
    const P = FpJ_to_FpE(J, p);
    expect(ell_is_inf(P)).toBe(false);
    if (!P.isInfinity) {
      expect(P.x).toBe(3n);
      expect(P.y).toBe(5n);
    }
  });

  it('roundtrips correctly', () => {
    const original = mkpoint(7n, 11n);
    const J = FpE_to_FpJ(original);
    const back = FpJ_to_FpE(J, p);
    expect(ell_is_inf(back)).toBe(false);
    if (!back.isInfinity) {
      expect(back.x).toBe(7n);
      expect(back.y).toBe(11n);
    }
  });

  it('roundtrips infinity correctly', () => {
    const inf = ellinf();
    const J = FpE_to_FpJ(inf);
    const back = FpJ_to_FpE(J, p);
    expect(ell_is_inf(back)).toBe(true);
  });
});

describe('ellinf_FpJ', () => {
  it('returns Jacobian point at infinity', () => {
    const J = ellinf_FpJ();
    expect(J.X).toBe(1n);
    expect(J.Y).toBe(1n);
    expect(J.Z).toBe(0n);
  });

  it('is recognized by FpJ_is_inf', () => {
    const J = ellinf_FpJ();
    expect(FpJ_is_inf(J)).toBe(true);
  });
});

describe('FpJ_is_inf', () => {
  it('returns true for Z = 0', () => {
    const J: JacobianPoint = { X: 5n, Y: 7n, Z: 0n };
    expect(FpJ_is_inf(J)).toBe(true);
  });

  it('returns false for Z != 0', () => {
    const J: JacobianPoint = { X: 5n, Y: 7n, Z: 1n };
    expect(FpJ_is_inf(J)).toBe(false);
  });
});

describe('FpE_isoncurve', () => {
  // Curve: y^2 = x^3 + x + 1 over F_23
  const E: ShortWeierstrassCurve = { a4: 1n, a6: 1n, p: 23n };

  it('returns true for point at infinity', () => {
    expect(FpE_isoncurve(E, ellinf())).toBe(true);
  });

  it('returns true for valid points', () => {
    // (0, 1) is on curve: 1 = 0 + 0 + 1
    expect(FpE_isoncurve(E, mkpoint(0n, 1n))).toBe(true);
    // (0, 22) is also on curve: 22^2 mod 23 = 484 mod 23 = 1
    expect(FpE_isoncurve(E, mkpoint(0n, 22n))).toBe(true);
  });

  it('returns false for invalid points', () => {
    // (0, 2) is not on curve: 2^2 = 4 != 1
    expect(FpE_isoncurve(E, mkpoint(0n, 2n))).toBe(false);
    // (1, 1) is not on curve: 1 != 1 + 1 + 1 = 3
    expect(FpE_isoncurve(E, mkpoint(1n, 1n))).toBe(false);
  });

  it('validates points from ellordinate', () => {
    const ys = ellordinate(E, 0n);
    for (const y of ys) {
      expect(FpE_isoncurve(E, mkpoint(0n, y))).toBe(true);
    }
  });

  it('validates points from random_FpE', () => {
    for (let i = 0; i < 10; i++) {
      const P = random_FpE(E);
      expect(FpE_isoncurve(E, P)).toBe(true);
    }
  });
});

describe('integration tests', () => {
  it('ellordinate and FpE_isoncurve agree', () => {
    const E: ShortWeierstrassCurve = { a4: 1n, a6: 1n, p: 23n };

    for (let x = 0n; x < 23n; x++) {
      const ys = ellordinate(E, x);
      for (const y of ys) {
        expect(FpE_isoncurve(E, mkpoint(x, y))).toBe(true);
      }
    }
  });

  it('coordinate conversions preserve curve membership', () => {
    const E: ShortWeierstrassCurve = { a4: 1n, a6: 1n, p: 23n };

    for (let i = 0; i < 10; i++) {
      const P = random_FpE(E);
      const J = FpE_to_FpJ(P);
      const P2 = FpJ_to_FpE(J, E.p);

      expect(FpE_isoncurve(E, P2)).toBe(true);
      if (!ell_is_inf(P) && !ell_is_inf(P2)) {
        expect(P2.x).toBe(P.x);
        expect(P2.y).toBe(P.y);
      }
    }
  });
});
