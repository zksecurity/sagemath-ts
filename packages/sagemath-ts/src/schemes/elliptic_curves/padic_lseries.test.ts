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
import { QQ } from '../../rings/rational_field.js';
import { EllipticCurve } from './constructor.js';
import type { EllipticCurveGeneric } from './ell_generic.js';
import {
  RationalElement,
  RationalRing,
  type pAdicEisensteinQuadraticElement,
  pAdicEisensteinQuadraticExtension,
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

/**
 * Mock curve exposing conductor() and ap(), used to exercise the reduction-type
 * branches of the constructor and of alpha().
 */
function makeMockCurve(conductor: bigint, aps: Map<bigint, bigint>) {
  return {
    a_invariants: () => [
      { isZero: () => true, toString: () => '0' },
      { isZero: () => true, toString: () => '0' },
      { isZero: () => false, toString: () => '1' },
      { isZero: () => false, toString: () => '-1' },
      { isZero: () => true, toString: () => '0' },
    ],
    ainvs: () => [0n, 0n, 1n, -1n, 0n],
    conductor: () => conductor,
    ap: (p: bigint) => {
      const v = aps.get(p);
      if (v === undefined) throw new Error(`no a_p for p=${p}`);
      return v;
    },
    is_isomorphic: () => false,
    j_invariant: () => ({ isZero: () => false, eq: () => false }),
    toString: () => `Elliptic Curve of conductor ${conductor}`,
  } as unknown as EllipticCurveGeneric<FieldElement>;
}

// Curve 11a1: y^2 + y = x^3 - x^2 - 10x - 20, conductor 11, a_5 = 1, a_11 = 1
const curve11a1 = makeMockCurve(
  11n,
  new Map([
    [5n, 1n],
    [11n, 1n],
  ])
);

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

    // padic_lseries.py:512-514: multiplicative reduction (p | N) gives K(a_p)
    // directly, before the ordinary/supersingular split.  Curve 11a1 has
    // N = 11 and a_11 = 1, so alpha(10) is 1 (not a unit root of x^2 - x + 11).
    it('returns a_p when p divides the conductor', () => {
      const L = new pAdicLseries(curve11a1, 11n);
      expect((L.alpha(10) as unknown as { lift(): bigint }).lift()).toBe(1n);
    });

    // Good ordinary reduction is unchanged: 11a1 at p = 5 has a_5 = 1, so
    // alpha is the unit root of x^2 - x + 5 in Z_5.
    it('returns the unit root for good ordinary reduction', () => {
      const L = new pAdicLseries(curve11a1, 5n);
      const a = L.alpha(5);
      const val = (a as unknown as { lift(): bigint }).lift();
      const p5 = 5n ** 5n;
      // alpha^2 - a_p*alpha + p == 0 mod p^5 and alpha is a unit
      expect((((val * val - 1n * val + 5n) % p5) + p5) % p5).toBe(0n);
      expect(val % 5n).not.toBe(0n);
    });
  });

  describe('semi-stability check', () => {
    // padic_lseries.py:182-183
    it('rejects primes with p^2 dividing the conductor', () => {
      const badCurve = makeMockCurve(49n, new Map([[7n, 0n]]));
      expect(() => new pAdicLseries(badCurve, 7n)).toThrow(NotImplementedError);
    });

    it('accepts primes of semi-stable reduction', () => {
      expect(() => new pAdicLseries(curve11a1, 11n)).not.toThrow();
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

    // padic_lseries.py:868 reduces eta mod (p-1) (mod 2 when p = 2) *before*
    // the quadratic-twist compatibility check, so eta = p-1 is the trivial
    // Teichmueller component and a twist is accepted.
    it('reduces eta modulo p-1 before the quadratic-twist check', () => {
      const L = new pAdicLseriesOrdinary(mockCurve, 5n);
      // eta = 4 == 0 (mod 4): must NOT be rejected as a non-zero component
      expect(() => L.series(2, -3, 5, 4)).toThrow(/requires modular symbols and p-adic arithmetic/);
      // eta = 5 == 1 (mod 4): still a non-trivial component, so rejected
      expect(() => L.series(2, -3, 5, 5)).toThrow(
        /quadratic twists only implemented for the 0th Teichmueller component/
      );
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
  const INF = Number.POSITIVE_INFINITY;

  // Sage doctests (padic_lseries.py:_e_bounds) for E = 11a1, p = 2.
  it('matches SageMath for p=2', () => {
    const L = new pAdicLseries(mockCurve, 2n);
    const acc = L as unknown as PAdicLseriesTestAccess;
    expect(acc._e_bounds(1, 10)).toEqual([INF, 1, 0, 0, 0, 0, 0, 0, 0, 0]);
    expect(acc._e_bounds(2, 10)).toEqual([INF, 2, 1, 1, 0, 0, 0, 0, 0, 0]);
    expect(acc._e_bounds(3, 10)).toEqual([INF, 3, 2, 2, 1, 1, 1, 1, 0, 0]);
    expect(acc._e_bounds(4, 10)).toEqual([INF, 4, 3, 3, 2, 2, 2, 2, 1, 1]);
  });

  it('should compute bounds for p=5', () => {
    const L = new pAdicLseries(mockCurve, 5n);
    const bounds = (L as unknown as PAdicLseriesTestAccess)._e_bounds(2, 10);
    // Derived from Sage's Lp._prec_bounds(3,10) == [+Infinity,1,1,1,1,0,0,0,0,0]
    // for 11a1 at p = 5, where _c_bound() == 1.
    expect(bounds).toEqual([INF, 2, 2, 2, 2, 1, 1, 1, 1, 1]);
  });
});

describe('pAdicLseriesOrdinary._prec_bounds', () => {
  // Sage's _prec_bounds subtracts _c_bound() from the e-bounds.  _c_bound()
  // needs E.galois_representation() and modular-symbol denominators, which are
  // not ported, so it now raises instead of silently returning 0 (which
  // over-reported the precision: Sage gives [+Infinity,1,1,1,1,0,...] for 11a1
  // at p = 5, the c = 0 version gave [+Infinity,2,2,2,2,1,...]).
  it('propagates the NotImplementedError from _c_bound', () => {
    const L = new pAdicLseriesOrdinary(mockCurve, 5n);
    expect(() => (L as unknown as PAdicLseriesTestAccess)._prec_bounds(3, 10)).toThrow(
      NotImplementedError
    );
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

// ---------------------------------------------------------------------------
// bernardi_sigma_function (padic_lseries.py:1613-1641)
// ---------------------------------------------------------------------------

/** Coefficients c_0 .. c_{n-1} of a power series, as strings. */
function coeffStrings(f: { list(): Array<{ toString(): string }> }, n: number): string[] {
  const list = f.list();
  const out: string[] = [];
  for (let i = 0; i < n; i++) out.push(list[i] === undefined ? '0' : String(list[i]));
  return out;
}

const curveOverQQ = (ainvs: bigint[]): EllipticCurveGeneric<FieldElement> =>
  EllipticCurve(QQ as never, ainvs as never) as unknown as EllipticCurveGeneric<FieldElement>;

describe('pAdicLseries.bernardi_sigma_function', () => {
  // padic_lseries.py:1621-1624
  //   sage: E = EllipticCurve('14a')
  //   sage: L = E.padic_lseries(5)
  //   sage: L.bernardi_sigma_function(prec=5)
  //   z + 1/24*z^3 + 29/384*z^5 - 8399/322560*z^7 - 291743/92897280*z^9 + O(z^10)
  it("matches SageMath's doctest for curve 14a", () => {
    const L = new pAdicLseries(curveOverQQ([1n, 0n, 1n, 4n, -6n]), 5n);
    const sigma = L.bernardi_sigma_function(5);
    expect(sigma.prec()).toBe(10);
    expect(coeffStrings(sigma, 10)).toEqual([
      '0',
      '1',
      '0',
      '1/24',
      '0',
      '29/384',
      '0',
      '-8399/322560',
      '0',
      '-291743/92897280',
    ]);
    expect(sigma.toString()).toBe(
      'z + 1/24*z^3 + 29/384*z^5 - 8399/322560*z^7 - 291743/92897280*z^9 + O(z^10)'
    );
  });

  // Independent oracle (no doctest involved): for a short Weierstrass curve
  // y^2 = x^3 + A*x + B the formal x in the coordinate z = log(t) is the
  // Weierstrass p-function of the lattice with g2 = -4A, g3 = -4B, and
  // sigma = z*exp(h) is built so that h'' = 1/z^2 - x(z) = -sum_n c_n z^{2n}
  // where c_1 = g2/20, c_2 = g3/28 and
  // c_n = 3/((2n+3)(n-2)) * sum_{m=1}^{n-2} c_m*c_{n-1-m}.
  it("satisfies the Weierstrass p-function recursion for h''", () => {
    type R = { n: bigint; d: bigint };
    const g = (a: bigint, b: bigint): bigint => (b ? g(b, a % b) : a < 0n ? -a : a);
    const mk = (n: bigint, d: bigint): R => {
      if (d < 0n) {
        n = -n;
        d = -d;
      }
      const k = g(n < 0n ? -n : n, d) || 1n;
      return { n: n / k, d: d / k };
    };
    const add = (a: R, b: R) => mk(a.n * b.d + b.n * a.d, a.d * b.d);
    const mul = (a: R, b: R) => mk(a.n * b.n, a.d * b.d);
    const div = (a: R, b: R) => mk(a.n * b.d, a.d * b.n);
    const str = (a: R) => (a.d === 1n ? `${a.n}` : `${a.n}/${a.d}`);

    for (const [A, B] of [
      [-1n, 0n],
      [0n, 1n],
      [-43n, 166n],
      [2n, -7n],
    ] as const) {
      const L = new pAdicLseries(curveOverQQ([0n, 0n, 0n, A, B]), 5n);
      const sigma = L.bernardi_sigma_function(8);
      const h = sigma._shiftRight(1).log();
      const hpp = h.derivative().derivative();

      const c: R[] = [mk(0n, 1n)];
      c[1] = div(mk(-4n * A, 1n), mk(20n, 1n));
      c[2] = div(mk(-4n * B, 1n), mk(28n, 1n));
      for (let n = 3; n <= 10; n++) {
        let s = mk(0n, 1n);
        for (let m = 1; m <= n - 2; m++) s = add(s, mul(c[m]!, c[n - 1 - m]!));
        c[n] = mul(mk(3n, BigInt((2 * n + 3) * (n - 2))), s);
      }
      const got = coeffStrings(hpp, hpp.prec());
      for (let k = 0; k < hpp.prec(); k++) {
        const want = k === 0 || k % 2 === 1 ? '0' : str(mul(mk(-1n, 1n), c[k / 2]!));
        expect(`z^${k}: ${got[k]}`).toBe(`z^${k}: ${want}`);
      }
      expect(hpp.prec()).toBeGreaterThan(8);
    }
  });

  it('is an odd series starting with z, for every curve tried', () => {
    for (const ainvs of [
      [1n, 0n, 1n, 4n, -6n], // 14a
      [0n, 0n, 1n, -1n, 0n], // 37a
      [1n, -1n, 1n, 0n, 0n], // 53a
      [0n, 1n, 1n, 0n, 0n], // 43a
      [3n, 2n, -4n, -2n, 5n],
    ]) {
      const sigma = new pAdicLseries(curveOverQQ(ainvs), 5n).bernardi_sigma_function(6);
      const cs = coeffStrings(sigma, 11);
      expect(cs[0]).toBe('0');
      expect(cs[1]).toBe('1');
      for (let k = 2; k < 11; k += 2) expect(cs[k]).toBe('0');
    }
  });

  it('rejects non-positive precision', () => {
    const L = new pAdicLseries(curveOverQQ([1n, 0n, 1n, 4n, -6n]), 5n);
    expect(() => L.bernardi_sigma_function(0)).toThrow(ValueError);
  });

  it('still refuses a curve with no formal group', () => {
    const L = new pAdicLseries(mockCurve, 5n);
    expect(() => L.bernardi_sigma_function(5)).toThrow(NotImplementedError);
  });
});

// ---------------------------------------------------------------------------
// alpha at a supersingular prime (padic_lseries.py:513-518)
// ---------------------------------------------------------------------------

describe('pAdicLseries.alpha at a supersingular prime', () => {
  // Curve 37a has a_3 = -3, so 3 is a prime of supersingular reduction.
  // padic_lseries.py:476-480
  //   sage: L = E.padic_lseries(3)
  //   sage: alpha = L.alpha(10); alpha
  //   alpha + O(alpha^21)
  //   sage: alpha^2 - E.ap(3)*alpha + 3
  //   O(alpha^22)
  const curve37a = makeMockCurve(37n, new Map([[3n, -3n]]));

  it("matches SageMath's doctest for 37a at p = 3", () => {
    const L = new pAdicLseries(curve37a, 3n);
    const alpha = L.alpha(10) as pAdicEisensteinQuadraticElement;
    expect(alpha.toString()).toBe('alpha + O(alpha^21)');
    expect(alpha.valuation()).toBe(1);
    expect(alpha.precision_relative()).toBe(20);

    const A = alpha.parent();
    const check = alpha.mul(alpha).sub(A.__call__(-3n).mul(alpha)).add(A.__call__(3n));
    expect(check.toString()).toBe('O(alpha^22)');
    expect(check.is_zero()).toBe(true);
  });

  it('describes the extension the way SageMath builds it', () => {
    const A = (
      new pAdicLseries(curve37a, 3n).alpha(10) as pAdicEisensteinQuadraticElement
    ).parent();
    expect(A.degree()).toBe(2);
    expect(A.e()).toBe(2);
    expect(A.f()).toBe(1);
    expect(A.precision_cap()).toBe(20);
    expect(A.defining_polynomial()).toEqual([3n, 3n, 1n]); // x^2 + 3x + 3
  });

  it('caches alpha per precision', () => {
    const L = new pAdicLseries(curve37a, 3n);
    expect(L.alpha(10)).toBe(L.alpha(10));
  });
});

describe('pAdicEisensteinQuadraticExtension', () => {
  // Z -> A is a ring homomorphism, alpha is a root of the defining polynomial,
  // v(alpha) = 1, v(p) = 2, and the alpha-adic expansion reconstructs its element.
  it('is a ring, for several Eisenstein polynomials', () => {
    for (const [p, ap] of [
      [3n, -3n],
      [3n, 0n],
      [3n, 3n],
      [5n, 0n],
      [5n, 5n],
      [7n, -7n],
      [2n, 0n],
      [2n, 2n],
    ] as const) {
      const A = new pAdicEisensteinQuadraticExtension(p, ap, 8);
      const al = A.gen();
      expect(al.mul(al).sub(A.__call__(ap).mul(al)).add(A.__call__(p)).is_zero()).toBe(true);
      expect(al.valuation()).toBe(1);
      expect(A.__call__(p).valuation()).toBe(2);
      for (let a = -12; a <= 12; a++) {
        for (let b = -12; b <= 12; b += 5) {
          expect(
            A.__call__(BigInt(a))
              .add(A.__call__(BigInt(b)))
              .eq(A.__call__(BigInt(a + b)))
          ).toBe(true);
          expect(
            A.__call__(BigInt(a))
              .mul(A.__call__(BigInt(b)))
              .eq(A.__call__(BigInt(a * b)))
          ).toBe(true);
        }
      }
      for (let a = 1; a <= 20; a++) {
        const xs = [
          A.__call__(BigInt(a)),
          A.__call__(BigInt(a)).mul(al),
          al.pow(2).add(A.__call__(BigInt(a))),
          al.pow(3).add(A.__call__(BigInt(a))),
          // negative valuation: exercises the denominator path of expansion()
          A.__call__(BigInt(a)).inv(),
          A.__call__(BigInt(a)).mul(al).inv(),
        ];
        for (const x of xs) {
          if (x.is_zero()) continue;
          // reconstruct from the alpha-adic expansion
          let acc = A.zero();
          for (const [k, d] of x.expansion()) acc = acc.add(A.__call__(d).mul(al.pow(k)));
          expect(acc.sub(x).is_zero()).toBe(true);
          // inverse
          expect(x.mul(x.inv()).eq(A.one())).toBe(true);
        }
      }
    }
  });

  it('rejects a non-Eisenstein defining polynomial', () => {
    expect(() => new pAdicEisensteinQuadraticExtension(5n, 3n, 10)).toThrow(ValueError);
  });
});
