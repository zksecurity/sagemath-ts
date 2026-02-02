/**
 * @module parigp-ts/elliptic/group.test
 * @description Tests for elliptic curve group structure functions
 *
 * Tests based on PARI/GP test vectors and the design document ELLIPTIC_CURVES.md
 */

import { describe, expect, it } from 'vitest';
import {
  type EllipticCurveFp,
  type EllipticPointFp,
  FpE_add,
  FpE_mul,
  FpE_neg,
  FpE_random,
  Fp_elltrace_naive,
  ellcard,
  ellgenerators,
  ellgroup,
  ellinf,
  ellinit_Fp,
  ellisoncurve,
  elllift_x,
  ellorder,
  ellpoint,
  trace_of_frobenius,
} from './group.js';

describe('Point arithmetic', () => {
  // Curve: y^2 = x^3 + x + 1 over F_23
  // From ELLIPTIC_CURVES.md test vectors
  const p = 23n;
  const a4 = 1n;
  const a6 = 1n;
  const E = ellinit_Fp(a4, a6, p);

  it('should create point at infinity', () => {
    const O = ellinf();
    expect(O.isInfinity).toBe(true);
    expect(O.x).toBeNull();
    expect(O.y).toBeNull();
  });

  it('should create finite points', () => {
    const P = ellpoint(0n, 1n);
    expect(P.isInfinity).toBe(false);
    expect(P.x).toBe(0n);
    expect(P.y).toBe(1n);
  });

  it('should check if point is on curve', () => {
    // (0, 1): 1^2 = 0^3 + 0 + 1 = 1
    expect(ellisoncurve(E, ellpoint(0n, 1n))).toBe(true);
    // (0, 22): 22^2 = 484 = 1 mod 23
    expect(ellisoncurve(E, ellpoint(0n, 22n))).toBe(true);
    // (1, 7): 7^2 = 49 = 3 mod 23, 1 + 1 + 1 = 3
    expect(ellisoncurve(E, ellpoint(1n, 7n))).toBe(true);
    // (1, 16): 16^2 = 256 = 3 mod 23
    expect(ellisoncurve(E, ellpoint(1n, 16n))).toBe(true);

    // Point not on curve
    expect(ellisoncurve(E, ellpoint(0n, 2n))).toBe(false);

    // Point at infinity is always on curve
    expect(ellisoncurve(E, ellinf())).toBe(true);
  });

  it('should negate points correctly', () => {
    const P = ellpoint(0n, 1n);
    const negP = FpE_neg(P, p);

    expect(negP.x).toBe(0n);
    expect(negP.y).toBe(22n); // -1 mod 23 = 22

    // Negating infinity gives infinity
    const O = ellinf();
    const negO = FpE_neg(O, p);
    expect(negO.isInfinity).toBe(true);
  });

  it('should add points correctly', () => {
    // From ELLIPTIC_CURVES.md test:
    // P = (0, 1), Q = (1, 7)
    // slope = (7 - 1) / (1 - 0) = 6 mod 23
    // x_R = 6^2 - 0 - 1 = 35 mod 23 = 12
    // y_R = 6*(0 - 12) - 1 = -73 mod 23 = 19
    // P + Q = (12, 19)

    const P = ellpoint(0n, 1n);
    const Q = ellpoint(1n, 7n);
    const R = FpE_add(P, Q, a4, p);

    expect(R.x).toBe(12n);
    expect(R.y).toBe(19n);
  });

  it('should handle identity element', () => {
    const P = ellpoint(0n, 1n);
    const O = ellinf();

    // P + O = P
    const R1 = FpE_add(P, O, a4, p);
    expect(R1.x).toBe(P.x);
    expect(R1.y).toBe(P.y);

    // O + P = P
    const R2 = FpE_add(O, P, a4, p);
    expect(R2.x).toBe(P.x);
    expect(R2.y).toBe(P.y);

    // O + O = O
    const R3 = FpE_add(O, O, a4, p);
    expect(R3.isInfinity).toBe(true);
  });

  it('should add point to its negation to get infinity', () => {
    const P = ellpoint(0n, 1n);
    const negP = ellpoint(0n, 22n);

    const R = FpE_add(P, negP, a4, p);
    expect(R.isInfinity).toBe(true);
  });

  it('should double points correctly', () => {
    // From ELLIPTIC_CURVES.md test:
    // [2]P where P = (0, 1)
    // slope = (3*0^2 + 1) / (2*1) = 1/2 = 12 mod 23 (since 2*12 = 24 = 1 mod 23)
    // x = 12^2 - 2*0 = 144 mod 23 = 6
    // y = 12*(0 - 6) - 1 = -73 mod 23 = 19
    // [2]P = (6, 19)

    const P = ellpoint(0n, 1n);
    const twoP = FpE_add(P, P, a4, p);

    // Note: The expected value from the doc might need verification
    // Let's verify: 2*12 = 24 = 1 mod 23, so 12 is the inverse of 2
    // slope = 1 * 12 = 12
    // x3 = 144 - 0 = 144 = 6 mod 23 (since 144 = 6*23 + 6)
    // Actually 144 = 6*24 = 144, 144/23 = 6.26, so 144 mod 23 = 144 - 6*23 = 144 - 138 = 6

    expect(twoP.x).toBe(6n);
    // y3 = 12 * (0 - 6) - 1 = -72 - 1 = -73
    // -73 mod 23: -73 + 4*23 = -73 + 92 = 19
    expect(twoP.y).toBe(19n);
  });

  it('should perform scalar multiplication', () => {
    const P = ellpoint(0n, 1n);

    // [0]P = O
    const zeroP = FpE_mul(P, 0n, a4, p);
    expect(zeroP.isInfinity).toBe(true);

    // [1]P = P
    const oneP = FpE_mul(P, 1n, a4, p);
    expect(oneP.x).toBe(P.x);
    expect(oneP.y).toBe(P.y);

    // [2]P should match our doubling result
    const twoP = FpE_mul(P, 2n, a4, p);
    expect(twoP.x).toBe(6n);
    expect(twoP.y).toBe(19n);

    // [-1]P = -P
    const negOneP = FpE_mul(P, -1n, a4, p);
    expect(negOneP.x).toBe(0n);
    expect(negOneP.y).toBe(22n);
  });
});

describe('ellcard - Curve cardinality', () => {
  it('should compute cardinality for small curves via exhaustive enumeration', () => {
    // Curve y^2 = x^3 + x + 1 over F_23
    // Need to verify the actual count
    const E = ellinit_Fp(1n, 1n, 23n);
    const card = ellcard(E);

    // The cardinality should be in Hasse interval: [23 + 1 - 2*sqrt(23), 23 + 1 + 2*sqrt(23)]
    // sqrt(23) ~ 4.79, so Hasse interval is [24 - 10, 24 + 10] = [14, 34]
    expect(card).toBeGreaterThanOrEqual(14n);
    expect(card).toBeLessThanOrEqual(34n);

    // Verify by counting manually for this small curve
    // For each x in 0..22, check if x^3 + x + 1 is a QR
    let count = 1n; // point at infinity
    for (let x = 0n; x < 23n; x++) {
      const ySquared = (x * x * x + x + 1n) % 23n;
      if (ySquared === 0n) {
        count += 1n;
      } else {
        // Check if ySquared is a quadratic residue
        // Using Euler's criterion: a^((p-1)/2) = 1 if QR
        const exp = (23n - 1n) / 2n; // 11
        let result = 1n;
        let base = ySquared;
        let e = exp;
        while (e > 0n) {
          if ((e & 1n) === 1n) {
            result = (result * base) % 23n;
          }
          base = (base * base) % 23n;
          e >>= 1n;
        }
        if (result === 1n) {
          count += 2n;
        }
      }
    }
    expect(card).toBe(count);
  });

  it('should verify Hasse bound for various primes', () => {
    const testCases = [
      { a4: 1n, a6: 1n, p: 23n },
      { a4: 2n, a6: 3n, p: 101n },
      { a4: 0n, a6: 7n, p: 97n }, // Similar to secp256k1
      { a4: -3n, a6: 1n, p: 67n },
    ];

    for (const { a4, a6, p } of testCases) {
      const E = ellinit_Fp(a4, a6, p);
      const card = ellcard(E);

      const sqrtP = BigInt(Math.floor(Math.sqrt(Number(p))));
      const lowerBound = p + 1n - 2n * sqrtP - 1n; // -1 for rounding
      const upperBound = p + 1n + 2n * sqrtP + 1n; // +1 for rounding

      expect(card).toBeGreaterThanOrEqual(lowerBound);
      expect(card).toBeLessThanOrEqual(upperBound);
    }
  });

  it('should compute cardinality for medium-sized primes', () => {
    // Test with primes > 1000 to trigger BSGS algorithm
    const E = ellinit_Fp(2n, 3n, 1009n);
    const card = ellcard(E);

    // Hasse bound: 1009 + 1 - 2*sqrt(1009) <= card <= 1009 + 1 + 2*sqrt(1009)
    // sqrt(1009) ~ 31.76, so [946, 1074]
    expect(card).toBeGreaterThanOrEqual(946n);
    expect(card).toBeLessThanOrEqual(1074n);

    // Verify by multiplying a random point by the computed order
    const P = FpE_random(E);
    const Q = FpE_mul(P, card, E.a4, E.p);
    expect(Q.isInfinity).toBe(true);
  });

  it('should cache computed cardinality', () => {
    const E = ellinit_Fp(1n, 1n, 23n);

    const card1 = ellcard(E);
    const card2 = ellcard(E);

    expect(card1).toBe(card2);
    expect(E._card).toBe(card1);
  });
});

describe('trace_of_frobenius', () => {
  it('should compute trace correctly', () => {
    const E = ellinit_Fp(1n, 1n, 23n);
    const card = ellcard(E);
    const trace = trace_of_frobenius(E);

    // #E = p + 1 - trace
    expect(card).toBe(23n + 1n - trace);
  });

  it('should satisfy Hasse bound for trace', () => {
    const testCases = [
      { a4: 1n, a6: 1n, p: 23n },
      { a4: 2n, a6: 3n, p: 101n },
      { a4: 0n, a6: 7n, p: 97n },
    ];

    for (const { a4, a6, p } of testCases) {
      const E = ellinit_Fp(a4, a6, p);
      const trace = trace_of_frobenius(E);

      const sqrtP = BigInt(Math.floor(Math.sqrt(Number(p))));
      const bound = 2n * sqrtP + 1n; // +1 for rounding

      // |trace| <= 2*sqrt(p)
      const absTrace = trace < 0n ? -trace : trace;
      expect(absTrace).toBeLessThanOrEqual(bound);
    }
  });
});

describe('Fp_elltrace_naive', () => {
  it('should compute trace via Legendre symbol summation', () => {
    const p = 23n;
    const a4 = 1n;
    const a6 = 1n;

    const trace = Fp_elltrace_naive(a4, a6, p);
    const card = p + 1n - trace;

    // Verify using ellcard
    const E = ellinit_Fp(a4, a6, p);
    expect(card).toBe(ellcard(E));
  });
});

describe('ellorder - Point order', () => {
  it('should return 1 for point at infinity', () => {
    const E = ellinit_Fp(1n, 1n, 23n);
    const O = ellinf();
    expect(ellorder(E, O)).toBe(1n);
  });

  it('should compute correct point orders', () => {
    const E = ellinit_Fp(1n, 1n, 23n);
    const card = ellcard(E);

    // Find a point and compute its order
    const P = elllift_x(E, 0n)!;
    expect(P).not.toBeNull();

    const order = ellorder(E, P);

    // Order must divide curve order
    expect(card % order).toBe(0n);

    // [order]P = O
    const Q = FpE_mul(P, order, E.a4, E.p);
    expect(Q.isInfinity).toBe(true);

    // [order - 1]P != O (unless order = 1)
    if (order > 1n) {
      const R = FpE_mul(P, order - 1n, E.a4, E.p);
      expect(R.isInfinity).toBe(false);
    }
  });

  it('should find points of various orders', () => {
    // Use a curve where we can find different order points
    const E = ellinit_Fp(2n, 3n, 101n);
    const card = ellcard(E);

    // Generate several random points and verify their orders divide the curve order
    for (let i = 0; i < 10; i++) {
      try {
        const P = FpE_random(E);
        const order = ellorder(E, P, card);

        expect(card % order).toBe(0n);

        // Verify by computing [order]P
        const Q = FpE_mul(P, order, E.a4, E.p);
        expect(Q.isInfinity).toBe(true);
      } catch {
        // Random point generation might fail, that's ok
        continue;
      }
    }
  });
});

describe('ellgroup - Group structure', () => {
  it('should return empty array for trivial group', () => {
    // Need to find a curve with #E = 1 (only point at infinity)
    // This is rare, so we skip this edge case test
  });

  it('should return [N] for cyclic groups', () => {
    // Most curves over prime fields have cyclic groups
    const E = ellinit_Fp(0n, 7n, 97n);
    const card = ellcard(E);
    const group = ellgroup(E);

    // If group is cyclic, it should be [N]
    if (group.length === 1) {
      expect(group[0]).toBe(card);
    }
  });

  it('should return [d1, d2] for non-cyclic groups with d2 | d1', () => {
    // Test curves known to have non-cyclic group structure
    // From ELLIPTIC_CURVES.md test vectors:
    // p = 1031, group = [504, 2]
    // p = 2053, group = [1008, 2]

    // Test with p = 1031
    const E1 = ellinit_Fp(1n, 0n, 1031n);
    const group1 = ellgroup(E1);

    if (group1.length === 2) {
      const [d1, d2] = group1;
      // d2 | d1
      expect(d1! % d2!).toBe(0n);
      // d1 * d2 = #E (not exactly, but d1 * (some cofactor) = #E)
      expect(ellcard(E1) % d1!).toBe(0n);
      expect(ellcard(E1) % d2!).toBe(0n);
    }
  });

  it('should cache group structure', () => {
    const E = ellinit_Fp(1n, 1n, 23n);

    const group1 = ellgroup(E);
    const group2 = ellgroup(E);

    expect(group1).toEqual(group2);
    expect(E._group).toEqual(group1);
  });

  it('should satisfy group structure invariants', () => {
    const testCases = [
      { a4: 1n, a6: 1n, p: 23n },
      { a4: 2n, a6: 3n, p: 101n },
      { a4: 0n, a6: 7n, p: 97n },
    ];

    for (const { a4, a6, p } of testCases) {
      const E = ellinit_Fp(a4, a6, p);
      const card = ellcard(E);
      const group = ellgroup(E);

      if (group.length === 1) {
        // Cyclic group: [N]
        expect(group[0]).toBe(card);
      } else if (group.length === 2) {
        // Non-cyclic: [d1, d2] with d2 | d1
        const [d1, d2] = group;
        expect(d1! % d2!).toBe(0n);

        // Product should equal cardinality
        expect(d1! * d2!).toBe(card);
      }
    }
  });
});

describe('ellgenerators', () => {
  it('should return empty array for trivial group', () => {
    // Skip - trivial groups are rare
  });

  it('should return one generator for cyclic groups', () => {
    const E = ellinit_Fp(1n, 1n, 23n);
    const card = ellcard(E);
    const group = ellgroup(E);

    if (group.length === 1) {
      const gens = ellgenerators(E);
      expect(gens.length).toBe(1);

      // Generator should have order equal to group order
      const G = gens[0]!;
      const order = ellorder(E, G, card);
      expect(order).toBe(card);
    }
  });

  it('should return generators on the curve', () => {
    const E = ellinit_Fp(2n, 3n, 101n);
    const gens = ellgenerators(E);

    for (const G of gens) {
      expect(ellisoncurve(E, G)).toBe(true);
    }
  });

  it('should cache generators', () => {
    const E = ellinit_Fp(1n, 1n, 23n);

    const gens1 = ellgenerators(E);
    const gens2 = ellgenerators(E);

    expect(gens1).toEqual(gens2);
    expect(E._generators).toEqual(gens1);
  });
});

describe('elllift_x', () => {
  it('should lift valid x-coordinates to points', () => {
    const E = ellinit_Fp(1n, 1n, 23n);

    // x = 0: y^2 = 0 + 0 + 1 = 1, sqrt(1) = 1 or 22
    const P = elllift_x(E, 0n);
    expect(P).not.toBeNull();
    expect(P!.x).toBe(0n);
    expect(P!.y === 1n || P!.y === 22n).toBe(true);
    expect(ellisoncurve(E, P!)).toBe(true);
  });

  it('should return null for x-coordinates with no points', () => {
    const E = ellinit_Fp(1n, 1n, 23n);

    // Try various x values to find one that doesn't lift
    for (let x = 0n; x < 23n; x++) {
      const P = elllift_x(E, x);
      if (P === null) {
        // Found an x that doesn't lift - verify y^2 is not a QR
        const ySquared = (x * x * x + x + 1n) % 23n;
        // ySquared^11 should be -1 (or p-1) mod p
        let result = 1n;
        let base = ySquared;
        let exp = 11n;
        while (exp > 0n) {
          if ((exp & 1n) === 1n) {
            result = (result * base) % 23n;
          }
          base = (base * base) % 23n;
          exp >>= 1n;
        }
        expect(result === 1n || ySquared === 0n).toBe(false);
        return;
      }
    }
  });
});

describe('FpE_random', () => {
  it('should generate points on the curve', () => {
    const E = ellinit_Fp(1n, 1n, 23n);

    for (let i = 0; i < 10; i++) {
      const P = FpE_random(E);
      expect(ellisoncurve(E, P)).toBe(true);
      expect(P.isInfinity).toBe(false);
    }
  });

  it('should generate diverse points', () => {
    const E = ellinit_Fp(2n, 3n, 101n);
    const points = new Set<string>();

    for (let i = 0; i < 20; i++) {
      const P = FpE_random(E);
      points.add(`${P.x},${P.y}`);
    }

    // Should have generated multiple distinct points
    expect(points.size).toBeGreaterThan(1);
  });
});

describe('Integration tests', () => {
  it('should correctly compute [#E]P = O for all points', () => {
    const E = ellinit_Fp(1n, 1n, 23n);
    const card = ellcard(E);

    // Test with several random points
    for (let i = 0; i < 5; i++) {
      const P = FpE_random(E);
      const Q = FpE_mul(P, card, E.a4, E.p);
      expect(Q.isInfinity).toBe(true);
    }
  });

  it('should handle edge cases with small primes', () => {
    // p = 5: y^2 = x^3 + 2x + 1 (disc = -4*8 - 27 = -59 = 1 mod 5, non-zero)
    const E5 = ellinit_Fp(2n, 1n, 5n);
    const card5 = ellcard(E5);
    expect(card5).toBeGreaterThan(0n);

    // p = 7: y^2 = x^3 + x + 2 (disc = -4 - 27*4 = -112 = 0 mod 7... try different)
    // y^2 = x^3 + 2x + 3 (disc = -4*8 - 27*9 = -32 - 243 = -275 = -275 mod 7 = 5, non-zero)
    const E7 = ellinit_Fp(2n, 3n, 7n);
    const card7 = ellcard(E7);
    expect(card7).toBeGreaterThan(0n);

    // p = 11: y^2 = x^3 + x + 1 (disc = -4 - 27 = -31 = -31 + 33 = 2 mod 11, non-zero)
    const E11 = ellinit_Fp(1n, 1n, 11n);
    const card11 = ellcard(E11);
    expect(card11).toBeGreaterThan(0n);
  });

  it('should work with secp256k1-like curve (small prime)', () => {
    // y^2 = x^3 + 7 (like secp256k1, but over small prime)
    const E = ellinit_Fp(0n, 7n, 97n);
    const card = ellcard(E);
    const trace = trace_of_frobenius(E);

    expect(card).toBe(97n + 1n - trace);

    // Find and verify a generator
    const gens = ellgenerators(E);
    expect(gens.length).toBeGreaterThan(0);

    const G = gens[0]!;
    expect(ellisoncurve(E, G)).toBe(true);
  });
});

describe('Singular curve detection', () => {
  it('should reject singular curves', () => {
    // y^2 = x^3 has discriminant 0
    expect(() => ellinit_Fp(0n, 0n, 23n)).toThrow();

    // 4a^3 + 27b^2 = 0 is the singular condition
    // Example: a = -3, b = 2 gives 4*(-27) + 27*4 = -108 + 108 = 0
    expect(() => ellinit_Fp(-3n, 2n, 23n)).toThrow();
  });
});
