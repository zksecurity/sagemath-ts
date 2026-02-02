/**
 * Tests for p-adic L-series of elliptic curves
 *
 * Reference: sage/schemes/elliptic_curves/padic_lseries.py
 *
 * NOTE: Many tests are marked as skipped because the full implementation
 * requires modular symbols and p-adic arithmetic which are not yet
 * available in sagemath-ts.
 */

import { describe, expect, it, test } from 'bun:test';
import { NotImplementedError, ValueError } from '../../errors.js';
import type { EllipticCurveGeneric } from './ell_generic.js';
import {
  RationalElement,
  RationalRing,
  pAdicLseries,
  pAdicLseriesOrdinary,
  pAdicLseriesSupersingular,
  rational,
} from './padic_lseries.js';
import type { FieldElement } from './types.js';

interface PAdicLseriesTestAccess {
  _e_bounds: (n: number, prec: number) => [number, number];
  _prec_bounds: (n: number, prec: number) => [number, number];
  _get_series_from_cache: (n: number, prec: number, D: bigint, eta: number) => unknown;
}

// Mock elliptic curve for testing
const mockCurve = {
  a_invariants: () => [
    { isZero: () => true, toString: () => '0' },
    { isZero: () => true, toString: () => '0' },
    { isZero: () => true, toString: () => '1' },
    { isZero: () => false, toString: () => '-1' },
    { isZero: () => true, toString: () => '0' },
  ],
  ainvs: () => [0n, 0n, 1n, -1n, 0n],
  is_isomorphic: () => false,
  j_invariant: () => ({ isZero: () => false, eq: () => false }),
  toString: () => 'Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field',
} as unknown as EllipticCurveGeneric<FieldElement>;

describe('pAdicLseries', () => {
  describe('constructor', () => {
    it('should create a p-adic L-series with valid prime', () => {
      const L = new pAdicLseries(mockCurve, 5n);
      expect(L.prime()).toBe(5n);
      expect(L.elliptic_curve()).toBe(mockCurve);
    });

    it('should accept number for prime', () => {
      const L = new pAdicLseries(mockCurve, 7);
      expect(L.prime()).toBe(7n);
    });

    it('should reject non-prime values', () => {
      expect(() => new pAdicLseries(mockCurve, 4n)).toThrow(ValueError);
      expect(() => new pAdicLseries(mockCurve, 1n)).toThrow(ValueError);
      expect(() => new pAdicLseries(mockCurve, 0n)).toThrow(ValueError);
    });

    it('should accept different implementations', () => {
      const L1 = new pAdicLseries(mockCurve, 5n, { implementation: 'eclib' });
      const L2 = new pAdicLseries(mockCurve, 5n, { implementation: 'sage' });
      const L3 = new pAdicLseries(mockCurve, 5n, { implementation: 'num' });
      expect(L1.prime()).toBe(5n);
      expect(L2.prime()).toBe(5n);
      expect(L3.prime()).toBe(5n);
    });

    it('should reject invalid implementations', () => {
      const badOptions = { implementation: 'invalid' as unknown as 'eclib' | 'sage' | 'num' };
      expect(() => new pAdicLseries(mockCurve, 5n, badOptions)).toThrow(ValueError);
    });
  });

  describe('elliptic_curve', () => {
    it('should return the associated elliptic curve', () => {
      const L = new pAdicLseries(mockCurve, 5n);
      expect(L.elliptic_curve()).toBe(mockCurve);
    });
  });

  describe('prime', () => {
    it('should return the prime p', () => {
      const L3 = new pAdicLseries(mockCurve, 3n);
      const L5 = new pAdicLseries(mockCurve, 5n);
      const L7 = new pAdicLseries(mockCurve, 7n);

      expect(L3.prime()).toBe(3n);
      expect(L5.prime()).toBe(5n);
      expect(L7.prime()).toBe(7n);
    });
  });

  describe('toString', () => {
    it('should include the prime and curve', () => {
      const L = new pAdicLseries(mockCurve, 5n);
      const str = L.toString();
      expect(str).toContain('5');
      expect(str).toContain('adic L-series');
    });

    it('should indicate when not normalized', () => {
      const L = new pAdicLseries(mockCurve, 5n, { normalize: 'none' });
      expect(L.toString()).toContain('not normalized');
    });
  });

  describe('modular_symbol', () => {
    it('should throw NotImplementedError', () => {
      const L = new pAdicLseries(mockCurve, 5n);
      expect(() => L.modular_symbol(rational(0n, 1n))).toThrow(NotImplementedError);
    });
  });

  describe('measure', () => {
    it('should throw NotImplementedError', () => {
      const L = new pAdicLseries(mockCurve, 5n);
      expect(() => L.measure(1n, 2, 10)).toThrow(NotImplementedError);
    });
  });

  describe('alpha', () => {
    it('should throw NotImplementedError', () => {
      const L = new pAdicLseries(mockCurve, 5n);
      expect(() => L.alpha(10)).toThrow(NotImplementedError);
    });
  });

  describe('teichmuller', () => {
    it('should return Teichmuller lifts for residues mod p', () => {
      const L = new pAdicLseries(mockCurve, 5n);
      const lifts = L.teichmuller(10);
      // Should have p elements: [0, teich(1), teich(2), ..., teich(p-1)]
      expect(lifts.length).toBe(5);
      expect(lifts[0]).toBe(0n);
      // teich(1) = 1 for any p
      expect(lifts[1]).toBe(1n);
      // The other values should be non-zero
      expect(lifts[2]).not.toBe(0n);
      expect(lifts[3]).not.toBe(0n);
      expect(lifts[4]).not.toBe(0n);
    });

    it('should return p elements for prime p', () => {
      const L3 = new pAdicLseries(mockCurve, 3n);
      expect(L3.teichmuller(10).length).toBe(3);

      const L7 = new pAdicLseries(mockCurve, 7n);
      expect(L7.teichmuller(10).length).toBe(7);
    });
  });

  describe('series', () => {
    it('should throw NotImplementedError', () => {
      const L = new pAdicLseries(mockCurve, 5n);
      expect(() => L.series(2, 1, 5, 0)).toThrow(NotImplementedError);
    });
  });

  describe('is_ordinary / is_supersingular', () => {
    it('should throw NotImplementedError for base class', () => {
      const L = new pAdicLseries(mockCurve, 5n);
      expect(() => L.is_ordinary()).toThrow(NotImplementedError);
    });
  });

  describe('frobenius', () => {
    it('should throw NotImplementedError', () => {
      const L = new pAdicLseries(mockCurve, 5n);
      expect(() => L.frobenius(20)).toThrow(NotImplementedError);
    });
  });

  describe('bernardi_sigma_function', () => {
    it('should throw NotImplementedError', () => {
      const L = new pAdicLseries(mockCurve, 5n);
      expect(() => L.bernardi_sigma_function(20)).toThrow(NotImplementedError);
    });
  });
});

describe('pAdicLseriesOrdinary', () => {
  describe('is_ordinary / is_supersingular', () => {
    it('should return true for is_ordinary', () => {
      const L = new pAdicLseriesOrdinary(mockCurve, 5n);
      expect(L.is_ordinary()).toBe(true);
    });

    it('should return false for is_supersingular', () => {
      const L = new pAdicLseriesOrdinary(mockCurve, 5n);
      expect(L.is_supersingular()).toBe(false);
    });
  });

  describe('series', () => {
    it('should validate n parameter', () => {
      const L = new pAdicLseriesOrdinary(mockCurve, 5n);
      expect(() => L.series(0)).toThrow(ValueError);
    });

    it('should require n >= 2 for p=2', () => {
      const L = new pAdicLseriesOrdinary(mockCurve, 2n);
      expect(() => L.series(1)).toThrow(ValueError);
    });

    it('should validate prec parameter', () => {
      const L = new pAdicLseriesOrdinary(mockCurve, 5n);
      expect(() => L.series(2, 1, 0)).toThrow(ValueError);
    });

    it('should reject quadratic twists with non-zero eta', () => {
      const L = new pAdicLseriesOrdinary(mockCurve, 5n);
      expect(() => L.series(2, -3, 5, 1)).toThrow(NotImplementedError);
    });

    it('should accept valid fundamental discriminants', () => {
      const L5 = new pAdicLseriesOrdinary(mockCurve, 5n);
      const L3 = new pAdicLseriesOrdinary(mockCurve, 3n);
      const L7 = new pAdicLseriesOrdinary(mockCurve, 7n);
      // D = 1 is always valid (trivial twist)
      expect(() => L5.series(2, 1, 5, 0)).toThrow(NotImplementedError);
      // D = -3 is a fundamental discriminant (-3 = 1 mod 4, squarefree), coprime to 5
      expect(() => L5.series(2, -3, 5, 0)).toThrow(NotImplementedError);
      // D = -4 is a fundamental discriminant (-4/4 = -1 = 3 mod 4), coprime to 5
      expect(() => L5.series(2, -4, 5, 0)).toThrow(NotImplementedError);
      // D = 5 is a fundamental discriminant (5 = 1 mod 4, squarefree), coprime to 3
      expect(() => L3.series(2, 5, 5, 0)).toThrow(NotImplementedError);
      // D = 8 is a fundamental discriminant (8/4 = 2 = 2 mod 4), coprime to 3
      expect(() => L3.series(2, 8, 5, 0)).toThrow(NotImplementedError);
      // D = 12 is a fundamental discriminant (12/4 = 3 = 3 mod 4), coprime to 5
      expect(() => L5.series(2, 12, 5, 0)).toThrow(NotImplementedError);
    });

    it('should reject non-fundamental discriminants', () => {
      const L = new pAdicLseriesOrdinary(mockCurve, 5n);
      // D = 2 is not fundamental (2 = 2 mod 4)
      expect(() => L.series(2, 2, 5, 0)).toThrow(ValueError);
      // D = 3 is not fundamental (3 = 3 mod 4)
      expect(() => L.series(2, 3, 5, 0)).toThrow(ValueError);
      // D = -5 is not fundamental (-5 = 3 mod 4)
      expect(() => L.series(2, -5, 3, 0)).toThrow(ValueError);
      // D = 12 with non-squarefree D/4 is not fundamental
      // But 12/4 = 3 which is squarefree and 3 mod 4, so 12 IS fundamental
      // Try D = 16: 16/4 = 4 which is not squarefree
      expect(() => L.series(2, 16, 5, 0)).toThrow(ValueError);
      // D = 9 is not squarefree
      expect(() => L.series(2, 9, 5, 0)).toThrow(ValueError);
    });

    it('should reject twist not coprime to p', () => {
      const L = new pAdicLseriesOrdinary(mockCurve, 5n);
      // D = -20 would need to be coprime to p=5, but gcd(20, 5) = 5
      // First check that -20 is a fundamental discriminant
      // -20 = 0 mod 4, -20/4 = -5, -5 mod 4 = 3, and -5 is squarefree
      // So -20 is a fundamental discriminant but not coprime to 5
      expect(() => L.series(2, -20, 5, 0)).toThrow(ValueError);
    });
  });
});

describe('pAdicLseriesSupersingular', () => {
  describe('is_ordinary / is_supersingular', () => {
    it('should return false for is_ordinary', () => {
      const L = new pAdicLseriesSupersingular(mockCurve, 5n);
      expect(L.is_ordinary()).toBe(false);
    });

    it('should return true for is_supersingular', () => {
      const L = new pAdicLseriesSupersingular(mockCurve, 5n);
      expect(L.is_supersingular()).toBe(true);
    });
  });

  describe('series', () => {
    it('should validate n parameter', () => {
      const L = new pAdicLseriesSupersingular(mockCurve, 5n);
      expect(() => L.series(0)).toThrow(ValueError);
    });

    it('should require n >= 2 for p=2', () => {
      const L = new pAdicLseriesSupersingular(mockCurve, 2n);
      expect(() => L.series(1)).toThrow(ValueError);
    });

    it('should accept valid fundamental discriminants', () => {
      const L = new pAdicLseriesSupersingular(mockCurve, 5n);
      // D = -4 is fundamental
      expect(() => L.series(3, -4, 5, 0)).toThrow(NotImplementedError);
      // D = -7 is fundamental (-7 = 1 mod 4, squarefree)
      expect(() => L.series(3, -7, 5, 0)).toThrow(NotImplementedError);
    });

    it('should reject non-fundamental discriminants', () => {
      const L = new pAdicLseriesSupersingular(mockCurve, 5n);
      // D = 6 is not fundamental (6 = 2 mod 4)
      expect(() => L.series(3, 6, 5, 0)).toThrow(ValueError);
    });
  });

  describe('Dp_valued_series', () => {
    it('should throw NotImplementedError', () => {
      const L = new pAdicLseriesSupersingular(mockCurve, 5n);
      expect(() => L.Dp_valued_series(3, 1, 5)).toThrow(NotImplementedError);
    });
  });

  describe('Dp_valued_height', () => {
    it('should throw NotImplementedError', () => {
      const L = new pAdicLseriesSupersingular(mockCurve, 5n);
      expect(() => L.Dp_valued_height(20)).toThrow(NotImplementedError);
    });
  });

  describe('Dp_valued_regulator', () => {
    it('should throw NotImplementedError', () => {
      const L = new pAdicLseriesSupersingular(mockCurve, 5n);
      expect(() => L.Dp_valued_regulator(20)).toThrow(NotImplementedError);
    });
  });
});

describe('rational helper', () => {
  it('should create a rational number', () => {
    const r = rational(3n, 4n);
    expect(r.num).toBe(3n);
    expect(r.den).toBe(4n);
  });

  it('should reduce fractions', () => {
    const r = rational(6n, 8n);
    expect(r.num).toBe(3n);
    expect(r.den).toBe(4n);
  });

  it('should handle negative numerator', () => {
    const r = rational(-3n, 4n);
    expect(r.num).toBe(-3n);
    expect(r.den).toBe(4n);
  });

  it('should normalize negative denominator', () => {
    const r = rational(3n, -4n);
    expect(r.num).toBe(-3n);
    expect(r.den).toBe(4n);
  });

  it('should reject zero denominator', () => {
    expect(() => rational(1n, 0n)).toThrow(ValueError);
  });

  it('should handle zero numerator', () => {
    const r = rational(0n, 5n);
    expect(r.num).toBe(0n);
    expect(r.den).toBe(1n);
  });
});

describe('RationalElement', () => {
  it('should create from bigints', () => {
    const r = new RationalElement(3n, 4n);
    expect(r.num).toBe(3n);
    expect(r.den).toBe(4n);
  });

  it('should reduce fractions automatically', () => {
    const r = new RationalElement(6n, 8n);
    expect(r.num).toBe(3n);
    expect(r.den).toBe(4n);
  });

  it('should normalize negative denominators', () => {
    const r = new RationalElement(3n, -4n);
    expect(r.num).toBe(-3n);
    expect(r.den).toBe(4n);
  });

  it('should add correctly', () => {
    const a = new RationalElement(1n, 2n);
    const b = new RationalElement(1n, 3n);
    const sum = a.add(b);
    expect(sum.num).toBe(5n);
    expect(sum.den).toBe(6n);
  });

  it('should subtract correctly', () => {
    const a = new RationalElement(1n, 2n);
    const b = new RationalElement(1n, 3n);
    const diff = a.sub(b);
    expect(diff.num).toBe(1n);
    expect(diff.den).toBe(6n);
  });

  it('should multiply correctly', () => {
    const a = new RationalElement(2n, 3n);
    const b = new RationalElement(3n, 4n);
    const prod = a.mul(b);
    expect(prod.num).toBe(1n);
    expect(prod.den).toBe(2n);
  });

  it('should divide correctly', () => {
    const a = new RationalElement(1n, 2n);
    const b = new RationalElement(3n, 4n);
    const quot = a.div(b);
    expect(quot.num).toBe(2n);
    expect(quot.den).toBe(3n);
  });

  it('should negate correctly', () => {
    const a = new RationalElement(3n, 4n);
    const neg = a.neg();
    expect(neg.num).toBe(-3n);
    expect(neg.den).toBe(4n);
  });

  it('should test equality correctly', () => {
    const a = new RationalElement(1n, 2n);
    const b = new RationalElement(2n, 4n);
    expect(a.eq(b)).toBe(true);
    expect(a.eq(1)).toBe(false);
    expect(a.eq(0)).toBe(false);
    expect(new RationalElement(1n, 1n).eq(1)).toBe(true);
  });

  it('should test zero correctly', () => {
    expect(new RationalElement(0n, 1n).isZero()).toBe(true);
    expect(new RationalElement(1n, 2n).isZero()).toBe(false);
  });

  it('should test one correctly', () => {
    expect(new RationalElement(1n, 1n).isOne()).toBe(true);
    expect(new RationalElement(3n, 3n).isOne()).toBe(true);
    expect(new RationalElement(1n, 2n).isOne()).toBe(false);
  });

  it('should invert correctly', () => {
    const a = new RationalElement(3n, 4n);
    const inv = a.inv();
    expect(inv.num).toBe(4n);
    expect(inv.den).toBe(3n);
  });

  it('should throw on division by zero', () => {
    const a = new RationalElement(1n, 2n);
    const zero = new RationalElement(0n, 1n);
    expect(() => a.div(zero)).toThrow();
  });

  it('should convert to string correctly', () => {
    expect(new RationalElement(3n, 4n).toString()).toBe('3/4');
    expect(new RationalElement(5n, 1n).toString()).toBe('5');
    expect(new RationalElement(-2n, 3n).toString()).toBe('-2/3');
  });
});

describe('RationalRing', () => {
  const QQ = new RationalRing();

  it('should provide zero', () => {
    const zero = QQ.zero();
    expect(zero.isZero()).toBe(true);
    expect(zero.num).toBe(0n);
    expect(zero.den).toBe(1n);
  });

  it('should provide one', () => {
    const one = QQ.one();
    expect(one.isOne()).toBe(true);
    expect(one.num).toBe(1n);
    expect(one.den).toBe(1n);
  });

  it('should coerce bigint', () => {
    const r = QQ.__call__(5n);
    expect(r.num).toBe(5n);
    expect(r.den).toBe(1n);
  });

  it('should coerce number', () => {
    const r = QQ.__call__(5);
    expect(r.num).toBe(5n);
    expect(r.den).toBe(1n);
  });

  it('should coerce rational object', () => {
    const r = QQ.__call__({ num: 3n, den: 4n });
    expect(r.num).toBe(3n);
    expect(r.den).toBe(4n);
  });

  it('should be a field', () => {
    expect(QQ.is_field()).toBe(true);
  });

  it('should have characteristic zero', () => {
    expect(QQ.characteristic()).toBe(0n);
  });
});

describe('pAdicLseries._e_bounds', () => {
  it('should return correct bounds for small cases', () => {
    const L = new pAdicLseries(mockCurve, 2n);
    // Access protected method through the class
    const bounds = (L as unknown as PAdicLseriesTestAccess)._e_bounds(1, 10);
    // bounds[0] should be Infinity
    expect(bounds[0]).toBe(Number.POSITIVE_INFINITY);
    // bounds should be decreasing
    for (let i = 1; i < bounds.length - 1; i++) {
      expect(bounds[i]).toBeGreaterThanOrEqual(bounds[i + 1]);
    }
  });

  it('should compute bounds for p=5', () => {
    const L = new pAdicLseries(mockCurve, 5n);
    const bounds = (L as unknown as PAdicLseriesTestAccess)._e_bounds(2, 10);
    // First entry is Infinity
    expect(bounds[0]).toBe(Number.POSITIVE_INFINITY);
    // Should have prec entries
    expect(bounds.length).toBe(10);
  });
});

describe('pAdicLseriesOrdinary._prec_bounds', () => {
  it('should return bounds with c_bound correction', () => {
    const L = new pAdicLseriesOrdinary(mockCurve, 5n);
    // Access protected methods
    const bounds = (L as unknown as PAdicLseriesTestAccess)._prec_bounds(3, 10);
    expect(bounds[0]).toBe(Number.POSITIVE_INFINITY);
    expect(bounds.length).toBe(10);
  });
});

describe('pAdicLseriesSupersingular._prec_bounds', () => {
  it('should return alpha-adic bounds', () => {
    const L = new pAdicLseriesSupersingular(mockCurve, 5n);
    const bounds = (L as unknown as PAdicLseriesTestAccess)._prec_bounds(3, 5);
    expect(bounds[0]).toBe(Number.POSITIVE_INFINITY);
    expect(bounds.length).toBe(5);
  });
});

describe('pAdicLseries caching', () => {
  it('should cache and retrieve series', () => {
    const L = new pAdicLseries(mockCurve, 5n);
    // Test the cache methods
    const cached = (L as unknown as PAdicLseriesTestAccess)._get_series_from_cache(2, 5, 1n, 0);
    expect(cached).toBe(null);
    // Note: We can't easily test setting the cache without a proper series
  });
});
