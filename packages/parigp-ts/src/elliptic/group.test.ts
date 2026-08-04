/**
 * @module parigp-ts/elliptic/group.test
 * @description Tests for elliptic curve group structure functions
 *
 * Tests based on PARI/GP test vectors and the design document ELLIPTIC_CURVES.md
 */

import { describe, expect, it } from 'bun:test';
import { kronecker } from '../ff.js';
import {
  type EllipticCurveFp,
  type EllipticPointFp,
  FpE_add,
  FpE_mul,
  FpE_neg,
  FpE_random,
  Fp_ellcard_CM,
  Fp_ellcard_Shanks,
  Fp_ellj_get_CM,
  Fp_ellj_nodiv,
  Fp_elltrace_naive,
  ec_ap_cm,
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
  /**
   * Oracle: PARI/GP (via Sage 10.3), `ellinit([a4,a6],p).ellgroup()`.
   */
  it('should return [N] for cyclic groups', () => {
    const E = ellinit_Fp(0n, 7n, 97n);
    expect(ellcard(E)).toBe(79n);
    expect(ellgroup(E)).toEqual([79n]);

    const E2 = ellinit_Fp(1n, 0n, 1031n);
    expect(ellcard(E2)).toBe(1032n);
    expect(ellgroup(E2)).toEqual([1032n]);

    const E3 = ellinit_Fp(2n, 5n, 5003n);
    expect(ellcard(E3)).toBe(4928n);
    expect(ellgroup(E3)).toEqual([4928n]);
  });

  /**
   * Oracle: PARI/GP (via Sage 10.3)
   *   ellinit([1,0],2053).ellgroup() -> [1010, 2]
   *   ellinit([0,7],3001).ellgroup() -> [724, 4]
   *   ellinit([3,4],10007).ellgroup() -> [5036, 2]
   *   ellinit([2,3],1000003).ellgroup() -> [499854, 2]
   *   ellinit([5,0],1000003).ellgroup() -> [500002, 2]
   */
  it('should return [d1, d2] for non-cyclic groups with d2 | d1', () => {
    const E1 = ellinit_Fp(1n, 0n, 2053n);
    expect(ellcard(E1)).toBe(2020n);
    expect(ellgroup(E1)).toEqual([1010n, 2n]);

    const E2 = ellinit_Fp(0n, 7n, 3001n);
    expect(ellcard(E2)).toBe(2896n);
    expect(ellgroup(E2)).toEqual([724n, 4n]);

    const E3 = ellinit_Fp(3n, 4n, 10007n);
    expect(ellcard(E3)).toBe(10072n);
    expect(ellgroup(E3)).toEqual([5036n, 2n]);

    const E4 = ellinit_Fp(2n, 3n, 1000003n);
    expect(ellcard(E4)).toBe(999708n);
    expect(ellgroup(E4)).toEqual([499854n, 2n]);

    // j = 1728 (a6 = 0): ap_j1728 path of Fp_ellcard_Shanks
    const E5 = ellinit_Fp(5n, 0n, 1000003n);
    expect(ellcard(E5)).toBe(1000004n);
    expect(ellgroup(E5)).toEqual([500002n, 2n]);
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

      expect(group.length === 1 || group.length === 2).toBe(true);
      if (group.length === 1) {
        // Cyclic group: [N]
        expect(group[0]).toBe(card);
      } else {
        // Non-cyclic: [d1, d2] with d2 | d1 and d1*d2 = N
        const [d1, d2] = group;
        expect(d1! % d2!).toBe(0n);
        expect(d1! * d2!).toBe(card);
      }
    }
  });
});

describe('ellgenerators', () => {
  it('should return one generator for cyclic groups', () => {
    const E = ellinit_Fp(1n, 1n, 23n);
    const card = ellcard(E);
    expect(ellgroup(E)).toEqual([card]);

    const gens = ellgenerators(E);
    expect(gens.length).toBe(1);
    expect(ellorder(E, gens[0]!, card)).toBe(card);
  });

  it('should return two independent generators for non-cyclic groups', () => {
    // PARI: ellinit([1,0],2053).ellgroup() == [1010, 2]
    const E = ellinit_Fp(1n, 0n, 2053n);
    const [d1, d2] = ellgroup(E) as [bigint, bigint];
    expect([d1, d2]).toEqual([1010n, 2n]);

    const gens = ellgenerators(E);
    expect(gens.length).toBe(2);
    expect(ellorder(E, gens[0]!)).toBe(d1);

    // <G1, G2> must be the whole group: enumerate i*G1 + j*G2
    const seen = new Set<string>();
    for (let i = 0n; i < d1; i++) {
      const P = FpE_mul(gens[0]!, i, E.a4, E.p);
      for (let j = 0n; j < d2; j++) {
        const Q = FpE_add(P, FpE_mul(gens[1]!, j, E.a4, E.p), E.a4, E.p);
        seen.add(Q.isInfinity ? 'O' : `${Q.x},${Q.y}`);
      }
    }
    expect(BigInt(seen.size)).toBe(ellcard(E));
  });

  it('should return generators on the curve', () => {
    const E = ellinit_Fp(2n, 3n, 101n);
    const gens = ellgenerators(E);

    expect(gens.length).toBeGreaterThan(0);
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

  /**
   * PARI's `gen_ellgens` tests independence with the Weil pairing order, not
   * by enumerating the multiples of a point.  Regression for the curves that
   * used to make `ellgenerators` throw.
   */
  it('finds generators for every curve over GF(11)..GF(23)', () => {
    for (const p of [11n, 13n, 17n, 19n, 23n]) {
      for (let a4 = 0n; a4 < p; a4++) {
        for (let a6 = 0n; a6 < p; a6++) {
          const disc = mod(-16n * (4n * a4 * a4 * a4 + 27n * a6 * a6), p);
          if (disc === 0n) continue;
          const E = ellinit_Fp(a4, a6, p);
          const D = ellgroup(E);
          const gens = ellgenerators(E);
          expect(gens.length).toBe(D.length);
          for (const G of gens) expect(ellisoncurve(E, G)).toBe(true);
          if (gens.length >= 1) expect(ellorder(E, gens[0]!)).toBe(D[0]!);
        }
      }
    }
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

// ============================================================================
// Brute-force oracle sweeps
//
// These enumerate every affine point of the curve and derive #E and the
// elementary divisors independently of the implementation, then require
// *exact* equality with ellcard/ellgroup.  Replacing `ellgroup` by
// `[ellcard(E)]` or perturbing `ellcard` by one must make these fail.
// ============================================================================

/** a mod p, in [0, p). */
function mod(a: bigint, p: bigint): bigint {
  const r = a % p;
  return r < 0n ? r + p : r;
}

/** All affine points of y^2 = x^3 + a4 x + a6 over F_p. */
function affinePoints(a4: bigint, a6: bigint, p: bigint): [bigint, bigint][] {
  const pts: [bigint, bigint][] = [];
  for (let x = 0n; x < p; x++) {
    const rhs = mod(((x * x) % p) * x + a4 * x + a6, p);
    for (let y = 0n; y < p; y++) {
      if ((y * y) % p === rhs) pts.push([x, y]);
    }
  }
  return pts;
}

/** Order of an affine point, by repeated addition. */
function bruteOrder(P: [bigint, bigint], a4: bigint, p: bigint): bigint {
  const base: EllipticPointFp = { x: P[0], y: P[1], isInfinity: false };
  let Q: EllipticPointFp = base;
  let n = 1n;
  while (!Q.isInfinity) {
    Q = FpE_add(Q, base, a4, p);
    n++;
  }
  return n;
}

/** [] / [d1] / [d1, d2] from an exhaustive enumeration. */
function bruteGroup(a4: bigint, a6: bigint, p: bigint): { N: bigint; D: bigint[] } {
  const pts = affinePoints(a4, a6, p);
  const N = BigInt(pts.length) + 1n;
  if (N === 1n) return { N, D: [] };
  let d1 = 1n;
  for (const P of pts) {
    const o = bruteOrder(P, a4, p);
    if (o > d1) d1 = o;
  }
  const d2 = N / d1;
  return { N, D: d2 === 1n ? [d1] : [d1, d2] };
}

const smallPrimes = (lo: number, hi: number): bigint[] => {
  const out: bigint[] = [];
  for (let n = lo; n <= hi; n++) {
    let isP = n > 1;
    for (let d = 2; d * d <= n; d++) {
      if (n % d === 0) {
        isP = false;
        break;
      }
    }
    if (isP) out.push(BigInt(n));
  }
  return out;
};

describe('ellcard / ellgroup vs brute-force oracle', () => {
  it('matches exactly for every curve with 100 < p < 130, a4,a6 in [0,5)', () => {
    let checked = 0;
    for (const p of smallPrimes(101, 130)) {
      for (let a4 = 0n; a4 < 5n; a4++) {
        for (let a6 = 0n; a6 < 5n; a6++) {
          if (mod(-16n * (4n * a4 * a4 * a4 + 27n * a6 * a6), p) === 0n) continue;
          const { N, D } = bruteGroup(a4, a6, p);
          // 4 independent runs: ellgroup is randomized
          for (let run = 0; run < 4; run++) {
            const E = ellinit_Fp(a4, a6, p);
            expect(ellcard(E)).toBe(N);
            expect(ellgroup(E)).toEqual(D);
          }
          checked++;
        }
      }
    }
    expect(checked).toBeGreaterThan(100);
  });

  it('Fp_ellcard_Shanks matches point counting for 1000 < p < 1100', () => {
    let checked = 0;
    for (const p of smallPrimes(1000, 1100)) {
      for (let a4 = 0n; a4 < 8n; a4++) {
        for (let a6 = 0n; a6 < 8n; a6++) {
          if (mod(-16n * (4n * a4 * a4 * a4 + 27n * a6 * a6), p) === 0n) continue;
          // count points directly (cheaper than enumerating y as well)
          let N = 1n;
          for (let x = 0n; x < p; x++) {
            const rhs = mod(((x * x) % p) * x + a4 * x + a6, p);
            if (rhs === 0n) N += 1n;
            else {
              let r = 1n;
              let b = rhs;
              let e = (p - 1n) / 2n;
              while (e > 0n) {
                if (e & 1n) r = (r * b) % p;
                b = (b * b) % p;
                e >>= 1n;
              }
              if (r === 1n) N += 2n;
            }
          }
          // BSGS branch (ellcard() itself uses naive enumeration below 2048)
          expect(Fp_ellcard_Shanks(a4, a6, p)).toBe(N);
          checked++;
        }
      }
    }
    expect(checked).toBeGreaterThan(900);
  });

  /**
   * Regression for the three counterexamples reported in AUDIT-2026-07 C17.
   * Oracle: exhaustive point count (and PARI).
   */
  it('gets the C17 counterexamples right', () => {
    expect(Fp_ellcard_Shanks(0n, 3n, 1063n)).toBe(1129n);
    expect(Fp_ellcard_Shanks(0n, 1n, 1069n)).toBe(1008n);
    expect(Fp_ellcard_Shanks(2n, 0n, 1013n)).toBe(1058n);
  });

  /**
   * Oracle: PARI/GP (via Sage 10.3), `ellinit([a4,a6],p).ellcard()`.
   * These primes are far beyond exhaustive enumeration, so they exercise the
   * Shanks/Mestre BSGS path end to end.
   */
  it('matches PARI ellcard for large primes', () => {
    const cases: [bigint, bigint, bigint, bigint][] = [
      [1000003n, 2n, 3n, 999708n],
      [1000003n, 0n, 7n, 999007n],
      [1000003n, 5n, 0n, 1000004n],
      [1000033n, 1n, 1n, 1001287n],
      [999983n, 3n, 4n, 999072n],
      [1000000007n, 0n, 7n, 1000000008n],
      [1000000007n, 2n, 3n, 1000004178n],
      [1000000009n, 1n, 1n, 999970651n],
      [2147483647n, 0n, 7n, 2147444533n],
    ];
    for (const [p, a4, a6, want] of cases) {
      expect(ellcard(ellinit_Fp(a4, a6, p))).toBe(want);
    }
  });
});

describe('ellorder delegates factorization to Z_factor', () => {
  /**
   * The old local trial-division `factor()` had no primality short-circuit:
   * a prime bound of ~10^18 took ~16 s.  `Z_factor` returns immediately.
   */
  it('handles a 10^18 prime bound quickly', () => {
    const E = ellinit_Fp(2n, 3n, 1000003n);
    const P = FpE_random(E);
    const t0 = Date.now();
    const n = ellorder(E, P, 1000000000000000003n);
    expect(Date.now() - t0).toBeLessThan(1000);
    expect(n).toBeGreaterThan(0n);
  });
});

describe('elllift_x uses PARI Fp_sqrt normalization', () => {
  /**
   * PARI's Fp_sqrt returns the smallest of the two roots (arith1.c:1277).
   * y^2 = x^3 + 3 over F_11 at x = 1 gives y^2 = 4, so PARI returns 2 (the
   * old private Tonelli-Shanks in this file returned 9).
   */
  it('returns the smallest root', () => {
    const E = ellinit_Fp(0n, 3n, 11n);
    const P = elllift_x(E, 1n);
    expect(P).not.toBeNull();
    expect(P!.y).toBe(2n);
  });
});

describe('ellgenerators terminates for split-prime group structures', () => {
  /**
   * E/F_43: y^2 = x^3 + 7x + 8 has group [12, 3] with N0 = 2^2*3^2, and
   * `gen_ellgroup` frequently settles the primes 2 and 3 in different
   * iterations.  With PARI 2.18-dev's `*pm = m` the resulting exponent can be
   * 4, and `gen_ellgens` then loops forever (the pairing of two 4-torsion
   * points can never have order 3).  Using `g1` keeps d2 | m | d1.
   */
  it('never fails over 500 fresh runs', () => {
    for (let i = 0; i < 500; i++) {
      const E = ellinit_Fp(7n, 8n, 43n);
      expect(ellgroup(E)).toEqual([12n, 3n]);
      const gens = ellgenerators(E);
      expect(gens.length).toBe(2);
      expect(E._m! % 3n).toBe(0n);
      expect(12n % E._m!).toBe(0n);
    }
  });
});

describe('ellcard / ellgroup vs PARI on random large curves', () => {
  /**
   * Oracle: PARI/GP 2.15.4 (Sage 10.3),
   *   E = ellinit([a4,a6],p); [E.ellcard(), E.ellgroup()]
   * for pseudo-random (a4, a6) over primes from 10^5 to 2^32.
   */
  it('matches exactly', () => {
    const cases: [bigint, bigint, bigint, bigint, bigint[]][] = [
      [100003n, 18462n, 62415n, 100280n, [50140n, 2n]],
      [1000003n, 355500n, 801085n, 1001482n, [1001482n]],
      [15485863n, 4863367n, 1890876n, 15487877n, [15487877n]],
      [32452843n, 11082502n, 5056359n, 32451460n, [16225730n, 2n]],
      [1000000007n, 239810037n, 543121245n, 1000047980n, [1000047980n]],
      [2147483647n, 87844563n, 1974046288n, 2147461164n, [715820388n, 3n]],
      [4294967291n, 2686620090n, 2731624997n, 4295053965n, [4295053965n]],
      [4294967291n, 3276826518n, 1724795536n, 4294944022n, [4294944022n]],
    ];
    for (const [p, a4, a6, card, group] of cases) {
      const E = ellinit_Fp(a4, a6, p);
      expect(ellcard(E)).toBe(card);
      expect(ellgroup(E)).toEqual(group);
    }
  });
});

// ============================================================================
// CM by a principal order (FpE.c:1280-1421)
// ============================================================================

/** #E(Fp) by naive enumeration of the Legendre symbols. */
function bruteCard(a4: bigint, a6: bigint, p: bigint): bigint {
  const md = (a: bigint) => ((a % p) + p) % p;
  let s = 0n;
  for (let x = 0n; x < p; x++) s += BigInt(kronecker(md(x * x * x + a4 * x + a6), p));
  return p + 1n + s;
}

/** PARI FpE.c:1387-1402 - Fp_ellj_to_a4a6 */
function ellj_to_a4a6(j: bigint, p: bigint): [bigint, bigint] {
  const md = (a: bigint) => ((a % p) + p) % p;
  j = md(j);
  if (j === 0n) return [0n, 1n];
  if (j === md(1728n)) return [1n, 0n];
  const k = md(1728n - j);
  const kj = md(k * j);
  return [md(3n * kj), md(2n * kj * k)];
}

/** The thirteen class-number-one CM discriminants (FpE.c:646-661). */
const CM_TABLE: [number, bigint][] = [
  [-3, 0n],
  [-4, 1728n],
  [-7, -3375n],
  [-8, 8000n],
  [-11, -32768n],
  [-12, 54000n],
  [-16, 287496n],
  [-19, -884736n],
  [-27, -12288000n],
  [-28, 16581375n],
  [-43, -884736000n],
  [-67, -147197952000n],
  [-163, -262537412640768000n],
];

describe('Fp_ellcard_CM', () => {
  it('detects each of the thirteen CM discriminants', () => {
    // p = 1 (mod 4*163*...) so that all thirteen j-invariants stay distinct
    const p = 916169n;
    const seen = new Set<number>();
    for (const [CM, J] of CM_TABLE) {
      const [a4, a6] = ellj_to_a4a6(J, p);
      const [jn, jd] = Fp_ellj_nodiv(a4, a6, p);
      expect(Fp_ellj_get_CM(jn, jd, p)).toBe(CM);
      seen.add(CM);
    }
    expect(seen.size).toBe(13);
  });

  it('returns null for a curve without CM by a class-number-one order', () => {
    // j = 5 is not in the table modulo this prime
    const p = 1000003n;
    const [a4, a6] = ellj_to_a4a6(5n, p);
    const [jn, jd] = Fp_ellj_nodiv(a4, a6, p);
    expect(Fp_ellj_get_CM(jn, jd, p)).toBe(0);
    expect(Fp_ellcard_CM(a4, a6, p)).toBeNull();
  });

  it('agrees with naive point counting on every CM curve and its twists', () => {
    const md = (a: bigint, p: bigint) => ((a % p) + p) % p;
    const primes = [101n, 103n, 107n, 109n, 113n, 127n, 131n, 137n, 139n, 149n, 151n, 157n];
    let checked = 0;
    const seen = new Set<number>();
    for (const p of primes) {
      for (const [, J] of CM_TABLE) {
        const [a4b, a6b] = ellj_to_a4a6(J, p);
        for (let u = 1n; u < p; u++) {
          const a4 = md(a4b * u * u, p);
          const a6 = md(a6b * u * u * u, p);
          if (md(4n * a4 * a4 * a4 + 27n * a6 * a6, p) === 0n) continue;
          const [jn, jd] = Fp_ellj_nodiv(a4, a6, p);
          seen.add(Fp_ellj_get_CM(jn, jd, p));
          const got = Fp_ellcard_CM(a4, a6, p);
          expect(got).not.toBeNull();
          expect(got).toBe(bruteCard(a4, a6, p));
          checked++;
        }
      }
    }
    expect(checked).toBeGreaterThan(10000);
    // every discriminant of the table was exercised
    expect([...seen].sort((a, b) => b - a)).toEqual(CM_TABLE.map(([d]) => d));
  }, 60000);

  it('agrees with Shanks/Mestre on 32-bit primes', () => {
    const md = (a: bigint, p: bigint) => ((a % p) + p) % p;
    const primes = [2147483647n, 4294967291n, 3221225473n, 2971215073n];
    let checked = 0;
    for (const p of primes) {
      for (const [, J] of CM_TABLE) {
        const [a4b, a6b] = ellj_to_a4a6(J, p);
        for (const u of [1n, 2n, 3n, 5n, 7n, 11n]) {
          const a4 = md(a4b * u * u, p);
          const a6 = md(a6b * u * u * u, p);
          if (md(4n * a4 * a4 * a4 + 27n * a6 * a6, p) === 0n) continue;
          const cm = Fp_ellcard_CM(a4, a6, p);
          expect(cm).not.toBeNull();
          expect(cm).toBe(Fp_ellcard_Shanks(a4, a6, p));
          const t = p + 1n - cm!;
          expect(t * t <= 4n * p).toBe(true);
          checked++;
        }
      }
    }
    expect(checked).toBeGreaterThan(250);
  }, 60000);

  it('reproduces the published orders of the SECG Koblitz curves (j = 0)', () => {
    // #E = n * h with h = 1 for all four; values from SEC 2 v2, section 2.
    const cases: [bigint, bigint, bigint][] = [
      // secp160k1: p = 2^160 - 2^32 - 21389, y^2 = x^3 + 7
      [
        0xfffffffffffffffffffffffffffffffeffffac73n,
        7n,
        0x0100000000000000000001b8fa16dfab9aca16b6b3n,
      ],
      // secp192k1: p = 2^192 - 2^32 - 4553, y^2 = x^3 + 3
      [
        0xfffffffffffffffffffffffffffffffffffffffeffffee37n,
        3n,
        0xfffffffffffffffffffffffe26f2fc170f69466a74defd8dn,
      ],
      // secp224k1: p = 2^224 - 2^32 - 6803, y^2 = x^3 + 5
      [
        0xfffffffffffffffffffffffffffffffffffffffffffffffeffffe56dn,
        5n,
        0x010000000000000000000000000001dce8d2ec6184caf0a971769fb1f7n,
      ],
      // secp256k1: p = 2^256 - 2^32 - 977, y^2 = x^3 + 7
      [
        0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2fn,
        7n,
        0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141n,
      ],
    ];
    for (const [p, b, card] of cases) {
      expect(Fp_ellcard_CM(0n, b, p)).toBe(card);
      expect(ellcard(ellinit_Fp(0n, b, p))).toBe(card);
    }
  });

  it('ec_ap_cm agrees with Fp_ellcard_CM and flips sign on the twist', () => {
    const md = (a: bigint, p: bigint) => ((a % p) + p) % p;
    for (const p of [916169n, 2147483647n]) {
      let u = 2n;
      while (kronecker(u, p) !== -1) u++;
      for (const [CM, J] of CM_TABLE) {
        const [a4, a6] = ellj_to_a4a6(J, p);
        const ap = ec_ap_cm(CM, a4, a6, p);
        expect(ap).not.toBeNull();
        expect(Fp_ellcard_CM(a4, a6, p)).toBe(p + 1n - ap!);
        // the quadratic twist has trace -a_p (#E + #E' = 2p + 2)
        const t4 = md(a4 * u * u, p);
        const t6 = md(a6 * u * u * u, p);
        expect(ec_ap_cm(CM, t4, t6, p)).toBe(-ap!);
      }
    }
  });

  it('[#E]P = O for every CM curve (independent of the counting method)', () => {
    const md = (a: bigint, p: bigint) => ((a % p) + p) % p;
    for (const p of [1000003n, 2147483647n, 4294967291n]) {
      for (const [, J] of CM_TABLE) {
        const [a4b, a6b] = ellj_to_a4a6(J, p);
        for (const u of [1n, 3n, 5n]) {
          const a4 = md(a4b * u * u, p);
          const a6 = md(a6b * u * u * u, p);
          if (md(4n * a4 * a4 * a4 + 27n * a6 * a6, p) === 0n) continue;
          const N = Fp_ellcard_CM(a4, a6, p);
          expect(N).not.toBeNull();
          const t = p + 1n - N!;
          expect(t * t <= 4n * p).toBe(true);
          const E = ellinit_Fp(a4, a6, p);
          for (let k = 0; k < 4; k++) {
            const P = FpE_random(E);
            expect(FpE_mul(P, N!, a4, p).isInfinity).toBe(true);
          }
        }
      }
    }
  });

  it('handles the supersingular branches (a_p = 0)', () => {
    // p = 3 mod 4: y^2 = x^3 + a4 x is supersingular, ap_j1728 returns 0
    for (const p of [10007n, 100003n, 1000003n]) {
      expect(p % 4n).toBe(3n);
      expect(Fp_ellcard_CM(1n, 0n, p)).toBe(p + 1n);
    }
    // p = 2 mod 3: y^2 = x^3 + a6 is supersingular, ap_j0 returns 0
    for (const p of [10007n, 100019n, 1000037n]) {
      expect(p % 3n).toBe(2n);
      expect(Fp_ellcard_CM(0n, 1n, p)).toBe(p + 1n);
    }
  });
});
