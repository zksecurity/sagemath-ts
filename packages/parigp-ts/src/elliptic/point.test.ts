/**
 * Tests for Elliptic Curve Point Operations
 *
 * Test vectors from ELLIPTIC_CURVES.md section 5 and PARI/GP test files.
 */

import { describe, expect, test } from 'bun:test';
import {
  type EllipticPoint,
  FpE_to_FpJ,
  FpJ_add,
  FpJ_dbl,
  FpJ_is_inf,
  FpJ_neg,
  FpJ_to_FpE,
  type JacobianPoint,
  type ShortWeierstrassCurve,
  ell_is_inf,
  elladd,
  ellinf,
  ellinf_FpJ,
  ellisoncurve,
  ellmul,
  ellneg,
  ellsub,
  mkpoint,
} from './point.js';

// =============================================================================
// Test Curves
// =============================================================================

// Small test curve: y^2 = x^3 + x + 1 over F_23
// From ELLIPTIC_CURVES.md section 5.4
const smallCurve: ShortWeierstrassCurve = { a4: 1n, a6: 1n, p: 23n };

// Curve y^2 = x^3 over F_1009 (a4 = 0, a6 = 0)
// From ELLIPTIC_CURVES.md section 5.1
const zeroCurve: ShortWeierstrassCurve = { a4: 0n, a6: 0n, p: 1009n };

// secp256k1 (Bitcoin curve)
// From ELLIPTIC_CURVES.md section 5.3
const secp256k1: ShortWeierstrassCurve = {
  a4: 0n,
  a6: 7n,
  p: 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2fn,
};

// secp256k1 generator point
const secp256k1_G: EllipticPoint = mkpoint(
  0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798n,
  0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8n
);

// secp256k1 order
const secp256k1_n = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141n;

// =============================================================================
// Point at Infinity Tests
// =============================================================================

describe('Point at Infinity', () => {
  test('ellinf returns infinity point', () => {
    const inf = ellinf();
    expect(ell_is_inf(inf)).toBe(true);
  });

  test('ell_is_inf correctly identifies infinity', () => {
    expect(ell_is_inf(ellinf())).toBe(true);
    expect(ell_is_inf(mkpoint(0n, 0n))).toBe(false);
    expect(ell_is_inf(mkpoint(1n, 2n))).toBe(false);
  });

  test('ellinf_FpJ returns (1:1:0)', () => {
    const inf = ellinf_FpJ();
    expect(inf.X).toBe(1n);
    expect(inf.Y).toBe(1n);
    expect(inf.Z).toBe(0n);
  });

  test('FpJ_is_inf correctly identifies Jacobian infinity', () => {
    expect(FpJ_is_inf(ellinf_FpJ())).toBe(true);
    expect(FpJ_is_inf({ X: 1n, Y: 1n, Z: 1n })).toBe(false);
  });
});

// =============================================================================
// Jacobian Coordinate Tests
// =============================================================================

describe('Jacobian Coordinates', () => {
  test('FpE_to_FpJ converts infinity correctly', () => {
    const inf = FpE_to_FpJ(ellinf());
    expect(FpJ_is_inf(inf)).toBe(true);
  });

  test('FpE_to_FpJ converts finite point correctly', () => {
    const P = mkpoint(3n, 5n);
    const J = FpE_to_FpJ(P);
    expect(J.X).toBe(3n);
    expect(J.Y).toBe(5n);
    expect(J.Z).toBe(1n);
  });

  test('FpJ_to_FpE converts infinity correctly', () => {
    const E = FpJ_to_FpE(ellinf_FpJ(), 23n);
    expect(ell_is_inf(E)).toBe(true);
  });

  test('FpJ_to_FpE converts (X:Y:1) correctly', () => {
    const J: JacobianPoint = { X: 3n, Y: 5n, Z: 1n };
    const E = FpJ_to_FpE(J, 23n);
    expect(ell_is_inf(E)).toBe(false);
    if (!ell_is_inf(E)) {
      expect(E.x).toBe(3n);
      expect(E.y).toBe(5n);
    }
  });

  test('FpJ_to_FpE handles non-trivial Z', () => {
    // (6 : 20 : 2) mod 23
    // x = 6 / 2^2 = 6 / 4 = 6 * 6 = 36 mod 23 = 13
    // y = 20 / 2^3 = 20 / 8 = 20 * 3 = 60 mod 23 = 14
    const J: JacobianPoint = { X: 6n, Y: 20n, Z: 2n };
    const E = FpJ_to_FpE(J, 23n);
    expect(ell_is_inf(E)).toBe(false);
    if (!ell_is_inf(E)) {
      // Verify: 2^-1 mod 23 = 12, 4^-1 mod 23 = 6, 8^-1 mod 23 = 3
      expect(E.x).toBe((6n * 6n) % 23n); // 13n
      expect(E.y).toBe((20n * 3n) % 23n); // 14n
    }
  });

  test('FpJ_neg negates Y coordinate', () => {
    const J: JacobianPoint = { X: 3n, Y: 5n, Z: 1n };
    const neg = FpJ_neg(J, 23n);
    expect(neg.X).toBe(3n);
    expect(neg.Y).toBe(18n); // 23 - 5 = 18
    expect(neg.Z).toBe(1n);
  });

  test('round-trip conversion preserves point', () => {
    const P = mkpoint(0n, 1n); // Point on smallCurve
    const J = FpE_to_FpJ(P);
    const P2 = FpJ_to_FpE(J, 23n);
    expect(ell_is_inf(P2)).toBe(false);
    if (!ell_is_inf(P2)) {
      expect(P2.x).toBe(P.x);
      expect(P2.y).toBe(P.y);
    }
  });
});

// =============================================================================
// ellisoncurve Tests
// =============================================================================

describe('ellisoncurve', () => {
  test('point at infinity is always on curve', () => {
    expect(ellisoncurve(smallCurve, ellinf())).toBe(true);
    expect(ellisoncurve(secp256k1, ellinf())).toBe(true);
  });

  test('(0, 0) is on y^2 = x^3', () => {
    expect(ellisoncurve(zeroCurve, mkpoint(0n, 0n))).toBe(true);
  });

  test('(0, 1) is on y^2 = x^3 + x + 1 over F_23', () => {
    // 1^2 = 0^3 + 0 + 1 = 1
    expect(ellisoncurve(smallCurve, mkpoint(0n, 1n))).toBe(true);
  });

  test('(0, 22) is on y^2 = x^3 + x + 1 over F_23', () => {
    // 22^2 mod 23 = 484 mod 23 = 1
    // 0^3 + 0 + 1 = 1
    expect(ellisoncurve(smallCurve, mkpoint(0n, 22n))).toBe(true);
  });

  test('(1, 7) is on y^2 = x^3 + x + 1 over F_23', () => {
    // 7^2 mod 23 = 49 mod 23 = 3
    // 1^3 + 1 + 1 = 3
    expect(ellisoncurve(smallCurve, mkpoint(1n, 7n))).toBe(true);
  });

  test('(1, 16) is on y^2 = x^3 + x + 1 over F_23', () => {
    // 16^2 mod 23 = 256 mod 23 = 3
    // 1^3 + 1 + 1 = 3
    expect(ellisoncurve(smallCurve, mkpoint(1n, 16n))).toBe(true);
  });

  test('secp256k1 generator is on curve', () => {
    expect(ellisoncurve(secp256k1, secp256k1_G)).toBe(true);
  });

  test('invalid point is not on curve', () => {
    expect(ellisoncurve(smallCurve, mkpoint(0n, 2n))).toBe(false);
    expect(ellisoncurve(smallCurve, mkpoint(1n, 1n))).toBe(false);
  });
});

// =============================================================================
// ellneg Tests
// =============================================================================

describe('ellneg', () => {
  test('negation of infinity is infinity', () => {
    expect(ell_is_inf(ellneg(smallCurve, ellinf()))).toBe(true);
  });

  test('negation of (0, 1) is (0, 22) on F_23', () => {
    const neg = ellneg(smallCurve, mkpoint(0n, 1n));
    expect(ell_is_inf(neg)).toBe(false);
    if (!ell_is_inf(neg)) {
      expect(neg.x).toBe(0n);
      expect(neg.y).toBe(22n);
    }
  });

  test('negation of (1, 7) is (1, 16) on F_23', () => {
    const neg = ellneg(smallCurve, mkpoint(1n, 7n));
    expect(ell_is_inf(neg)).toBe(false);
    if (!ell_is_inf(neg)) {
      expect(neg.x).toBe(1n);
      expect(neg.y).toBe(16n);
    }
  });

  test('double negation returns original', () => {
    const P = mkpoint(1n, 7n);
    const neg = ellneg(smallCurve, P);
    const doubleNeg = ellneg(smallCurve, neg);
    expect(ell_is_inf(doubleNeg)).toBe(false);
    if (!ell_is_inf(doubleNeg)) {
      expect(doubleNeg.x).toBe(P.x);
      expect(doubleNeg.y).toBe(P.y);
    }
  });

  test('P + (-P) = O', () => {
    const P = mkpoint(0n, 1n);
    const negP = ellneg(smallCurve, P);
    const sum = elladd(smallCurve, P, negP);
    expect(ell_is_inf(sum)).toBe(true);
  });
});

// =============================================================================
// elladd Tests
// =============================================================================

describe('elladd', () => {
  test('P + O = P', () => {
    const P = mkpoint(0n, 1n);
    const sum = elladd(smallCurve, P, ellinf());
    expect(ell_is_inf(sum)).toBe(false);
    if (!ell_is_inf(sum)) {
      expect(sum.x).toBe(P.x);
      expect(sum.y).toBe(P.y);
    }
  });

  test('O + P = P', () => {
    const P = mkpoint(0n, 1n);
    const sum = elladd(smallCurve, ellinf(), P);
    expect(ell_is_inf(sum)).toBe(false);
    if (!ell_is_inf(sum)) {
      expect(sum.x).toBe(P.x);
      expect(sum.y).toBe(P.y);
    }
  });

  test('O + O = O', () => {
    const sum = elladd(smallCurve, ellinf(), ellinf());
    expect(ell_is_inf(sum)).toBe(true);
  });

  test('P + Q distinct points', () => {
    // From ELLIPTIC_CURVES.md section 5.4:
    // P = (0, 1), Q = (1, 7)
    // slope = (7 - 1) / (1 - 0) = 6 mod 23
    // x_R = 6^2 - 0 - 1 = 35 mod 23 = 12
    // y_R = 6*(0 - 12) - 1 = -73 mod 23 = 19 (since -73 + 4*23 = 19)
    const P = mkpoint(0n, 1n);
    const Q = mkpoint(1n, 7n);
    const R = elladd(smallCurve, P, Q);
    expect(ell_is_inf(R)).toBe(false);
    if (!ell_is_inf(R)) {
      expect(R.x).toBe(12n);
      expect(R.y).toBe(19n);
    }
    expect(ellisoncurve(smallCurve, R)).toBe(true);
  });

  test('P + P (doubling)', () => {
    // From ELLIPTIC_CURVES.md section 5.4:
    // [2]P where P = (0, 1)
    // slope = (3*0^2 + 1) / (2*1) = 1/2 = 12 mod 23 (since 2*12 = 24 = 1 mod 23)
    // x = 12^2 - 2*0 = 144 mod 23 = 6
    // y = 12*(0 - 6) - 1 = -73 mod 23 = 19
    const P = mkpoint(0n, 1n);
    const R = elladd(smallCurve, P, P);
    expect(ell_is_inf(R)).toBe(false);
    if (!ell_is_inf(R)) {
      expect(R.x).toBe(6n);
      expect(R.y).toBe(19n);
    }
    expect(ellisoncurve(smallCurve, R)).toBe(true);
  });

  test('P + (-P) = O', () => {
    const P = mkpoint(1n, 7n);
    const negP = mkpoint(1n, 16n);
    const sum = elladd(smallCurve, P, negP);
    expect(ell_is_inf(sum)).toBe(true);
  });

  test('result is on curve', () => {
    const P = mkpoint(0n, 1n);
    const Q = mkpoint(1n, 7n);
    const R = elladd(smallCurve, P, Q);
    expect(ellisoncurve(smallCurve, R)).toBe(true);
  });
});

// =============================================================================
// ellsub Tests
// =============================================================================

describe('ellsub', () => {
  test('P - O = P', () => {
    const P = mkpoint(0n, 1n);
    const diff = ellsub(smallCurve, P, ellinf());
    expect(ell_is_inf(diff)).toBe(false);
    if (!ell_is_inf(diff)) {
      expect(diff.x).toBe(P.x);
      expect(diff.y).toBe(P.y);
    }
  });

  test('O - P = -P', () => {
    const P = mkpoint(0n, 1n);
    const diff = ellsub(smallCurve, ellinf(), P);
    expect(ell_is_inf(diff)).toBe(false);
    if (!ell_is_inf(diff)) {
      expect(diff.x).toBe(0n);
      expect(diff.y).toBe(22n);
    }
  });

  test('P - P = O', () => {
    const P = mkpoint(0n, 1n);
    const diff = ellsub(smallCurve, P, P);
    expect(ell_is_inf(diff)).toBe(true);
  });

  test('P - Q = P + (-Q)', () => {
    const P = mkpoint(0n, 1n);
    const Q = mkpoint(1n, 7n);
    const diff = ellsub(smallCurve, P, Q);
    const negQ = ellneg(smallCurve, Q);
    const sum = elladd(smallCurve, P, negQ);
    expect(ell_is_inf(diff)).toBe(ell_is_inf(sum));
    if (!ell_is_inf(diff) && !ell_is_inf(sum)) {
      expect(diff.x).toBe(sum.x);
      expect(diff.y).toBe(sum.y);
    }
  });
});

// =============================================================================
// ellmul Tests
// =============================================================================

describe('ellmul', () => {
  test('[0]P = O', () => {
    const P = mkpoint(0n, 1n);
    const R = ellmul(smallCurve, P, 0n);
    expect(ell_is_inf(R)).toBe(true);
  });

  test('[1]P = P', () => {
    const P = mkpoint(0n, 1n);
    const R = ellmul(smallCurve, P, 1n);
    expect(ell_is_inf(R)).toBe(false);
    if (!ell_is_inf(R)) {
      expect(R.x).toBe(P.x);
      expect(R.y).toBe(P.y);
    }
  });

  test('[-1]P = -P', () => {
    const P = mkpoint(0n, 1n);
    const R = ellmul(smallCurve, P, -1n);
    expect(ell_is_inf(R)).toBe(false);
    if (!ell_is_inf(R)) {
      expect(R.x).toBe(0n);
      expect(R.y).toBe(22n);
    }
  });

  test('[2]P = P + P', () => {
    const P = mkpoint(0n, 1n);
    const R1 = ellmul(smallCurve, P, 2n);
    const R2 = elladd(smallCurve, P, P);
    expect(ell_is_inf(R1)).toBe(ell_is_inf(R2));
    if (!ell_is_inf(R1) && !ell_is_inf(R2)) {
      expect(R1.x).toBe(R2.x);
      expect(R1.y).toBe(R2.y);
      expect(R1.x).toBe(6n);
      expect(R1.y).toBe(19n);
    }
  });

  test('[3]P = [2]P + P', () => {
    const P = mkpoint(0n, 1n);
    const R1 = ellmul(smallCurve, P, 3n);
    const twoP = ellmul(smallCurve, P, 2n);
    const R2 = elladd(smallCurve, twoP, P);
    expect(ell_is_inf(R1)).toBe(ell_is_inf(R2));
    if (!ell_is_inf(R1) && !ell_is_inf(R2)) {
      expect(R1.x).toBe(R2.x);
      expect(R1.y).toBe(R2.y);
    }
  });

  test('[n]O = O', () => {
    expect(ell_is_inf(ellmul(smallCurve, ellinf(), 5n))).toBe(true);
    expect(ell_is_inf(ellmul(smallCurve, ellinf(), 100n))).toBe(true);
  });

  test('result is on curve', () => {
    const P = mkpoint(0n, 1n);
    for (let n = 1n; n <= 10n; n++) {
      const R = ellmul(smallCurve, P, n);
      if (!ell_is_inf(R)) {
        expect(ellisoncurve(smallCurve, R)).toBe(true);
      }
    }
  });

  test('scalar multiplication consistency', () => {
    const P = mkpoint(1n, 7n);
    // [5]P = [3]P + [2]P
    const fiveP = ellmul(smallCurve, P, 5n);
    const threeP = ellmul(smallCurve, P, 3n);
    const twoP = ellmul(smallCurve, P, 2n);
    const sum = elladd(smallCurve, threeP, twoP);
    expect(ell_is_inf(fiveP)).toBe(ell_is_inf(sum));
    if (!ell_is_inf(fiveP) && !ell_is_inf(sum)) {
      expect(fiveP.x).toBe(sum.x);
      expect(fiveP.y).toBe(sum.y);
    }
  });

  test('negative scalar multiplication', () => {
    const P = mkpoint(1n, 7n);
    // [-3]P = -([3]P)
    const negThreeP = ellmul(smallCurve, P, -3n);
    const threeP = ellmul(smallCurve, P, 3n);
    const negated = ellneg(smallCurve, threeP);
    expect(ell_is_inf(negThreeP)).toBe(ell_is_inf(negated));
    if (!ell_is_inf(negThreeP) && !ell_is_inf(negated)) {
      expect(negThreeP.x).toBe(negated.x);
      expect(negThreeP.y).toBe(negated.y);
    }
  });
});

// =============================================================================
// Group Law Tests
// =============================================================================

describe('Group Law Properties', () => {
  const P = mkpoint(0n, 1n);
  const Q = mkpoint(1n, 7n);
  const R = mkpoint(1n, 16n);

  test('identity: P + O = O + P = P', () => {
    const sum1 = elladd(smallCurve, P, ellinf());
    const sum2 = elladd(smallCurve, ellinf(), P);
    expect(ell_is_inf(sum1)).toBe(false);
    expect(ell_is_inf(sum2)).toBe(false);
    if (!ell_is_inf(sum1) && !ell_is_inf(sum2)) {
      expect(sum1.x).toBe(P.x);
      expect(sum1.y).toBe(P.y);
      expect(sum2.x).toBe(P.x);
      expect(sum2.y).toBe(P.y);
    }
  });

  test('inverse: P + (-P) = O', () => {
    const negP = ellneg(smallCurve, P);
    expect(ell_is_inf(elladd(smallCurve, P, negP))).toBe(true);
  });

  test('commutativity: P + Q = Q + P', () => {
    const PQ = elladd(smallCurve, P, Q);
    const QP = elladd(smallCurve, Q, P);
    expect(ell_is_inf(PQ)).toBe(ell_is_inf(QP));
    if (!ell_is_inf(PQ) && !ell_is_inf(QP)) {
      expect(PQ.x).toBe(QP.x);
      expect(PQ.y).toBe(QP.y);
    }
  });

  test('associativity: (P + Q) + R = P + (Q + R)', () => {
    const left = elladd(smallCurve, elladd(smallCurve, P, Q), R);
    const right = elladd(smallCurve, P, elladd(smallCurve, Q, R));
    expect(ell_is_inf(left)).toBe(ell_is_inf(right));
    if (!ell_is_inf(left) && !ell_is_inf(right)) {
      expect(left.x).toBe(right.x);
      expect(left.y).toBe(right.y);
    }
  });

  test('scalar distribution: [m+n]P = [m]P + [n]P', () => {
    const m = 3n;
    const n = 5n;
    const left = ellmul(smallCurve, P, m + n);
    const right = elladd(smallCurve, ellmul(smallCurve, P, m), ellmul(smallCurve, P, n));
    expect(ell_is_inf(left)).toBe(ell_is_inf(right));
    if (!ell_is_inf(left) && !ell_is_inf(right)) {
      expect(left.x).toBe(right.x);
      expect(left.y).toBe(right.y);
    }
  });

  test('scalar multiplication: [mn]P = [m]([n]P)', () => {
    const m = 3n;
    const n = 4n;
    const left = ellmul(smallCurve, P, m * n);
    const right = ellmul(smallCurve, ellmul(smallCurve, P, n), m);
    expect(ell_is_inf(left)).toBe(ell_is_inf(right));
    if (!ell_is_inf(left) && !ell_is_inf(right)) {
      expect(left.x).toBe(right.x);
      expect(left.y).toBe(right.y);
    }
  });
});

// =============================================================================
// Jacobian Doubling and Addition Tests
// =============================================================================

describe('Jacobian Operations', () => {
  test('FpJ_dbl of infinity returns infinity', () => {
    const inf = ellinf_FpJ();
    const result = FpJ_dbl(inf, 1n, 23n);
    expect(FpJ_is_inf(result)).toBe(true);
  });

  test('FpJ_dbl matches affine doubling', () => {
    const P = mkpoint(0n, 1n);
    const PJ = FpE_to_FpJ(P);
    const doubledJ = FpJ_dbl(PJ, 1n, 23n);
    const doubledAffine = FpJ_to_FpE(doubledJ, 23n);
    const expected = elladd(smallCurve, P, P);
    expect(ell_is_inf(doubledAffine)).toBe(ell_is_inf(expected));
    if (!ell_is_inf(doubledAffine) && !ell_is_inf(expected)) {
      expect(doubledAffine.x).toBe(expected.x);
      expect(doubledAffine.y).toBe(expected.y);
    }
  });

  test('FpJ_add of P + O returns P', () => {
    const P: JacobianPoint = { X: 3n, Y: 5n, Z: 1n };
    const inf = ellinf_FpJ();
    const result = FpJ_add(P, inf, 1n, 23n);
    expect(result.X).toBe(P.X);
    expect(result.Y).toBe(P.Y);
    expect(result.Z).toBe(P.Z);
  });

  test('FpJ_add of O + P returns P', () => {
    const P: JacobianPoint = { X: 3n, Y: 5n, Z: 1n };
    const inf = ellinf_FpJ();
    const result = FpJ_add(inf, P, 1n, 23n);
    expect(result.X).toBe(P.X);
    expect(result.Y).toBe(P.Y);
    expect(result.Z).toBe(P.Z);
  });

  test('FpJ_add matches affine addition', () => {
    const P = mkpoint(0n, 1n);
    const Q = mkpoint(1n, 7n);
    const PJ = FpE_to_FpJ(P);
    const QJ = FpE_to_FpJ(Q);
    const sumJ = FpJ_add(PJ, QJ, 1n, 23n);
    const sumAffine = FpJ_to_FpE(sumJ, 23n);
    const expected = elladd(smallCurve, P, Q);
    expect(ell_is_inf(sumAffine)).toBe(ell_is_inf(expected));
    if (!ell_is_inf(sumAffine) && !ell_is_inf(expected)) {
      expect(sumAffine.x).toBe(expected.x);
      expect(sumAffine.y).toBe(expected.y);
    }
  });

  test('FpJ_add handles same point (doubling)', () => {
    const P = mkpoint(0n, 1n);
    const PJ = FpE_to_FpJ(P);
    const sumJ = FpJ_add(PJ, PJ, 1n, 23n);
    const sumAffine = FpJ_to_FpE(sumJ, 23n);
    const expected = elladd(smallCurve, P, P);
    expect(ell_is_inf(sumAffine)).toBe(ell_is_inf(expected));
    if (!ell_is_inf(sumAffine) && !ell_is_inf(expected)) {
      expect(sumAffine.x).toBe(expected.x);
      expect(sumAffine.y).toBe(expected.y);
    }
  });

  test('FpJ_add handles inverse points', () => {
    const P = mkpoint(1n, 7n);
    const negP = mkpoint(1n, 16n);
    const PJ = FpE_to_FpJ(P);
    const negPJ = FpE_to_FpJ(negP);
    const sumJ = FpJ_add(PJ, negPJ, 1n, 23n);
    expect(FpJ_is_inf(sumJ)).toBe(true);
  });
});

// =============================================================================
// secp256k1 Tests
// =============================================================================

describe('secp256k1 Operations', () => {
  test('generator is on curve', () => {
    expect(ellisoncurve(secp256k1, secp256k1_G)).toBe(true);
  });

  test('[2]G is on curve', () => {
    const twoG = ellmul(secp256k1, secp256k1_G, 2n);
    expect(ellisoncurve(secp256k1, twoG)).toBe(true);
  });

  test('[3]G = [2]G + G', () => {
    const threeG = ellmul(secp256k1, secp256k1_G, 3n);
    const twoG = ellmul(secp256k1, secp256k1_G, 2n);
    const sum = elladd(secp256k1, twoG, secp256k1_G);
    expect(ell_is_inf(threeG)).toBe(ell_is_inf(sum));
    if (!ell_is_inf(threeG) && !ell_is_inf(sum)) {
      expect(threeG.x).toBe(sum.x);
      expect(threeG.y).toBe(sum.y);
    }
  });

  test('[n]G = O where n is the group order', () => {
    // This test may be slow for large n, so we just verify the concept
    // by checking that [n-1]G + G = O, which is equivalent
    const nMinus1_G = ellmul(secp256k1, secp256k1_G, secp256k1_n - 1n);
    expect(ellisoncurve(secp256k1, nMinus1_G)).toBe(true);

    // [n-1]G = -G
    const negG = ellneg(secp256k1, secp256k1_G);
    expect(ell_is_inf(nMinus1_G)).toBe(ell_is_inf(negG));
    if (!ell_is_inf(nMinus1_G) && !ell_is_inf(negG)) {
      expect(nMinus1_G.x).toBe(negG.x);
      expect(nMinus1_G.y).toBe(negG.y);
    }

    // Therefore [n]G = [n-1]G + G = -G + G = O
    const nG = elladd(secp256k1, nMinus1_G, secp256k1_G);
    expect(ell_is_inf(nG)).toBe(true);
  });

  test('known scalar multiplication value', () => {
    // [2]G on secp256k1
    const twoG = ellmul(secp256k1, secp256k1_G, 2n);
    expect(ell_is_inf(twoG)).toBe(false);
    if (!ell_is_inf(twoG)) {
      // Known values for 2*G on secp256k1 (verified by manual computation)
      expect(twoG.x).toBe(0xc6047f9441ed7d6d3045406e95c07cd85c778e4b8cef3ca7abac09b95c709ee5n);
      expect(twoG.y).toBe(0x1ae168fea63dc339a3c58419466ceaeef7f632653266d0e1236431a950cfe52an);
    }
  });
});

// =============================================================================
// Edge Cases
// =============================================================================

describe('Edge Cases', () => {
  test('point with y = 0 doubles to infinity', () => {
    // Find a point with y = 0 on some curve
    // y^2 = x^3 over F_7: (0, 0) has y = 0
    const curve: ShortWeierstrassCurve = { a4: 0n, a6: 0n, p: 7n };
    const P = mkpoint(0n, 0n);
    expect(ellisoncurve(curve, P)).toBe(true);
    const doubled = elladd(curve, P, P);
    expect(ell_is_inf(doubled)).toBe(true);
  });

  test('large scalar multiplication', () => {
    const P = mkpoint(0n, 1n);
    const largeN = 1000000n;
    const R = ellmul(smallCurve, P, largeN);
    // Just verify it completes and result is on curve
    if (!ell_is_inf(R)) {
      expect(ellisoncurve(smallCurve, R)).toBe(true);
    }
  });

  test('handles modular arithmetic edge cases', () => {
    // Find a valid point with x near p-1
    // For x = 17: RHS = 17^3 + 17 + 1 mod 23 = 4913 + 18 mod 23
    // 4913 mod 23 = 4913 - 213*23 = 4913 - 4899 = 14
    // RHS = 14 + 18 mod 23 = 32 mod 23 = 9 = 3^2
    // So (17, 3) is on the curve
    const P2 = mkpoint(17n, 3n);
    expect(ellisoncurve(smallCurve, P2)).toBe(true);

    const twoP2 = ellmul(smallCurve, P2, 2n);
    if (!ell_is_inf(twoP2)) {
      expect(ellisoncurve(smallCurve, twoP2)).toBe(true);
    }
  });
});
