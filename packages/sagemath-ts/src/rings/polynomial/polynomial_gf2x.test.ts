/**
 * Unit tests for GF2X polynomial module
 */
import { describe, expect, test } from 'bun:test';
import { set_random_seed } from '../../misc/randstate.js';
import {
  GF2X,
  GF2X_BuildIrred_list,
  GF2X_BuildSparseIrred_list,
  GF2X_Ring,
  buildIrred,
  buildRandomIrred,
  buildSparseIrred,
  distinctDegreeFactorization,
  equalDegreeFactorization,
  factor,
  squareFreeDecomp,
} from './polynomial_gf2x.js';

describe('GF2X basic operations', () => {
  test('zero and one', () => {
    expect(GF2X.zero().isZero()).toBe(true);
    expect(GF2X.zero().bits).toBe(0n);
    expect(GF2X.one().isOne()).toBe(true);
    expect(GF2X.one().bits).toBe(1n);
    expect(GF2X.x().isX()).toBe(true);
    expect(GF2X.x().bits).toBe(2n);
  });

  test('fromBigInt', () => {
    // x^3 + x + 1 = binary 1011 = 11
    const p = GF2X.fromBigInt(11n);
    expect(p.degree()).toBe(3);
    expect(p.getCoeff(0)).toBe(1);
    expect(p.getCoeff(1)).toBe(1);
    expect(p.getCoeff(2)).toBe(0);
    expect(p.getCoeff(3)).toBe(1);
  });

  test('fromCoeffs', () => {
    const p = GF2X.fromCoeffs([1, 1, 0, 1]); // x^3 + x + 1
    expect(p.bits).toBe(11n);
    expect(p.degree()).toBe(3);
  });

  test('fromHex', () => {
    const p = GF2X.fromHex('b'); // 11 in hex = 1011 in binary
    expect(p.bits).toBe(11n);

    const p2 = GF2X.fromHex('0xFF');
    expect(p2.bits).toBe(255n);
  });

  test('monomial', () => {
    expect(GF2X.monomial(0).bits).toBe(1n);
    expect(GF2X.monomial(1).bits).toBe(2n);
    expect(GF2X.monomial(5).bits).toBe(32n);
    expect(GF2X.monomial(10).degree()).toBe(10);
  });

  test('degree', () => {
    expect(GF2X.zero().degree()).toBe(-1);
    expect(GF2X.one().degree()).toBe(0);
    expect(GF2X.x().degree()).toBe(1);
    expect(GF2X.fromBigInt(11n).degree()).toBe(3);
    expect(GF2X.fromBigInt(255n).degree()).toBe(7);
  });

  test('getCoeff and setCoeff', () => {
    const p = GF2X.fromBigInt(11n); // x^3 + x + 1
    expect(p.getCoeff(-1)).toBe(0); // out of bounds
    expect(p.getCoeff(0)).toBe(1);
    expect(p.getCoeff(1)).toBe(1);
    expect(p.getCoeff(2)).toBe(0);
    expect(p.getCoeff(3)).toBe(1);
    expect(p.getCoeff(100)).toBe(0); // out of bounds

    const p2 = p.setCoeff(2, 1);
    expect(p2.getCoeff(2)).toBe(1);
    expect(p2.bits).toBe(15n);

    const p3 = p2.setCoeff(0, 0);
    expect(p3.getCoeff(0)).toBe(0);
    expect(p3.bits).toBe(14n);
  });

  test('weight', () => {
    expect(GF2X.zero().weight()).toBe(0);
    expect(GF2X.one().weight()).toBe(1);
    expect(GF2X.fromBigInt(11n).weight()).toBe(3); // 1011 has 3 ones
    expect(GF2X.fromBigInt(255n).weight()).toBe(8); // all 8 bits set
  });

  test('toString', () => {
    expect(GF2X.zero().toString()).toBe('0');
    expect(GF2X.one().toString()).toBe('1');
    expect(GF2X.x().toString()).toBe('x');
    expect(GF2X.fromBigInt(11n).toString()).toBe('x^3 + x + 1');
    expect(GF2X.fromBigInt(7n).toString()).toBe('x^2 + x + 1');
    expect(GF2X.fromBigInt(3n).toString()).toBe('x + 1');
    expect(GF2X.fromBigInt(5n).toString()).toBe('x^2 + 1');
  });

  test('toCoeffs', () => {
    expect(GF2X.zero().toCoeffs()).toEqual([0]);
    expect(GF2X.one().toCoeffs()).toEqual([1]);
    expect(GF2X.x().toCoeffs()).toEqual([0, 1]);
    expect(GF2X.fromBigInt(11n).toCoeffs()).toEqual([1, 1, 0, 1]);
  });
});

describe('GF2X arithmetic', () => {
  test('addition (XOR)', () => {
    const p1 = GF2X.fromBigInt(11n); // x^3 + x + 1
    const p2 = GF2X.fromBigInt(7n); // x^2 + x + 1

    const sum = p1.add(p2);
    // (x^3 + x + 1) + (x^2 + x + 1) = x^3 + x^2 (x and 1 cancel)
    expect(sum.bits).toBe(12n); // 1100 in binary

    // In GF(2), a + a = 0
    expect(p1.add(p1).isZero()).toBe(true);
  });

  test('subtraction equals addition in GF(2)', () => {
    const p1 = GF2X.fromBigInt(11n);
    const p2 = GF2X.fromBigInt(7n);

    expect(p1.sub(p2).eq(p1.add(p2))).toBe(true);
  });

  test('negation is identity in GF(2)', () => {
    const p = GF2X.fromBigInt(11n);
    expect(p.neg().eq(p)).toBe(true);
  });

  test('multiplication', () => {
    // (x + 1) * (x + 1) = x^2 + 2x + 1 = x^2 + 1 (mod 2)
    const p = GF2X.fromBigInt(3n); // x + 1
    const p2 = p.mul(p);
    expect(p2.bits).toBe(5n); // x^2 + 1

    // x * x = x^2
    const x = GF2X.x();
    expect(x.mul(x).bits).toBe(4n);

    // (x^2 + 1) * x = x^3 + x
    const q = GF2X.fromBigInt(5n);
    expect(q.mul(x).bits).toBe(10n); // 1010 = x^3 + x

    // Multiply by zero
    expect(GF2X.fromBigInt(11n).mul(GF2X.zero()).isZero()).toBe(true);

    // Multiply by one
    const r = GF2X.fromBigInt(11n);
    expect(r.mul(GF2X.one()).eq(r)).toBe(true);
  });

  test('squaring', () => {
    // x^2 = x^2
    const x = GF2X.x();
    expect(x.sqr().bits).toBe(4n);

    // (x + 1)^2 = x^2 + 2x + 1 = x^2 + 1 (mod 2)
    const p = GF2X.fromBigInt(3n);
    expect(p.sqr().bits).toBe(5n);

    // (x^2 + x + 1)^2 = x^4 + x^2 + 1 (coefficients spread out)
    const q = GF2X.fromBigInt(7n);
    expect(q.sqr().bits).toBe(21n); // 10101

    // sqr() should equal mul(self)
    const r = GF2X.fromBigInt(13n);
    expect(r.sqr().eq(r.mul(r))).toBe(true);
  });

  test('shift operations', () => {
    const p = GF2X.fromBigInt(11n); // x^3 + x + 1

    // Left shift = multiply by x^n
    expect(p.leftShift(1).bits).toBe(22n);
    expect(p.leftShift(2).bits).toBe(44n);
    expect(p.mulByX().eq(p.leftShift(1))).toBe(true);

    // Right shift = divide by x^n (floor)
    expect(p.rightShift(1).bits).toBe(5n); // x^2 + 1
    expect(p.rightShift(2).bits).toBe(2n); // x

    // Negative shifts reverse direction
    expect(p.leftShift(-1).eq(p.rightShift(1))).toBe(true);
  });

  test('division', () => {
    // x^3 + x + 1 divided by x + 1
    const dividend = GF2X.fromBigInt(11n); // x^3 + x + 1
    const divisor = GF2X.fromBigInt(3n); // x + 1

    const [q, r] = dividend.divRem(divisor);
    // Verify: dividend = quotient * divisor + remainder
    expect(q.mul(divisor).add(r).eq(dividend)).toBe(true);

    // Division by smaller degree polynomial
    expect(GF2X.fromBigInt(3n).divRem(GF2X.fromBigInt(11n))).toEqual([
      GF2X.zero(),
      GF2X.fromBigInt(3n),
    ]);

    // Division by one
    expect(dividend.div(GF2X.one()).eq(dividend)).toBe(true);

    // Exact division
    const a = GF2X.fromBigInt(7n); // x^2 + x + 1
    const b = GF2X.fromBigInt(3n); // x + 1
    const c = a.mul(b); // product
    expect(c.div(b).eq(a)).toBe(true);
    expect(c.mod(b).isZero()).toBe(true);
  });

  test('division by zero throws', () => {
    expect(() => GF2X.fromBigInt(11n).divRem(GF2X.zero())).toThrow();
  });

  test('pow', () => {
    const p = GF2X.fromBigInt(3n); // x + 1

    expect(p.pow(0).isOne()).toBe(true);
    expect(p.pow(1).eq(p)).toBe(true);
    expect(p.pow(2).bits).toBe(5n); // x^2 + 1
    expect(p.pow(3).eq(p.mul(p).mul(p))).toBe(true);

    expect(GF2X.zero().pow(0).isOne()).toBe(true);
    expect(GF2X.zero().pow(5).isZero()).toBe(true);

    // Negative exponent throws
    expect(() => p.pow(-1)).toThrow();
  });

  test('powMod', () => {
    const base = GF2X.fromBigInt(7n); // x^2 + x + 1
    const mod = GF2X.fromBigInt(11n); // x^3 + x + 1 (irreducible)

    expect(base.powMod(0, mod).isOne()).toBe(true);
    expect(base.powMod(1, mod).eq(base.mod(mod))).toBe(true);

    // In GF(2^3) with modulus x^3 + x + 1, the group has order 7
    // So base^7 = base (for non-zero base)
    // Actually base^(2^3 - 1) = base^7 should be 1 for units
    const result = base.powMod(7, mod);
    expect(result.isOne()).toBe(true);
  });

  test('trunc', () => {
    const p = GF2X.fromBigInt(255n); // 8 bits set

    expect(p.trunc(4).bits).toBe(15n); // lower 4 bits
    expect(p.trunc(1).bits).toBe(1n); // just constant
    expect(p.trunc(8).bits).toBe(255n); // all bits
    expect(p.trunc(0).isZero()).toBe(true);
  });
});

describe('GF2X GCD and extended GCD', () => {
  test('gcd', () => {
    const a = GF2X.fromBigInt(15n); // x^3 + x^2 + x + 1 = (x+1)(x^2+1) = (x+1)^3
    const b = GF2X.fromBigInt(3n); // x + 1

    const g = a.gcd(b);
    expect(a.isDivisibleBy(g)).toBe(true);
    expect(b.isDivisibleBy(g)).toBe(true);

    // gcd(0, a) = a
    expect(GF2X.zero().gcd(a).eq(a)).toBe(true);

    // gcd(a, a) = a
    expect(a.gcd(a).eq(a)).toBe(true);

    // Coprime polynomials
    const f = GF2X.fromBigInt(11n); // x^3 + x + 1 (irreducible)
    const h = GF2X.fromBigInt(3n); // x + 1
    expect(f.gcd(h).isOne()).toBe(true);
  });

  test('xgcd', () => {
    const a = GF2X.fromBigInt(11n); // x^3 + x + 1
    const b = GF2X.fromBigInt(7n); // x^2 + x + 1

    const [g, s, t] = a.xgcd(b);

    // Verify: s*a + t*b = g
    expect(s.mul(a).add(t.mul(b)).eq(g)).toBe(true);
  });

  test('invMod', () => {
    const mod = GF2X.fromBigInt(11n); // x^3 + x + 1 (irreducible)
    const a = GF2X.fromBigInt(7n); // x^2 + x + 1

    const inv = a.invMod(mod);
    // a * inv = 1 (mod mod)
    expect(a.mul(inv).mod(mod).isOne()).toBe(true);
  });

  test('invMod throws for non-invertible', () => {
    const mod = GF2X.fromBigInt(6n); // x^2 + x (not irreducible)
    const a = GF2X.fromBigInt(2n); // x, shares factor with mod

    expect(() => a.invMod(mod)).toThrow();
  });
});

describe('GF2X irreducibility testing', () => {
  test('is_irreducible for small polynomials', () => {
    // Degree 1: x + 1 is the only irreducible
    expect(GF2X.fromBigInt(3n).is_irreducible()).toBe(true); // x + 1
    expect(GF2X.fromBigInt(2n).is_irreducible()).toBe(true); // x

    // Degree 2: x^2 + x + 1 is the only irreducible (x^2 + 1 = (x+1)^2)
    expect(GF2X.fromBigInt(7n).is_irreducible()).toBe(true); // x^2 + x + 1
    expect(GF2X.fromBigInt(5n).is_irreducible()).toBe(false); // x^2 + 1 = (x+1)^2

    // Degree 3: x^3 + x + 1 and x^3 + x^2 + 1 are the irreducibles
    expect(GF2X.fromBigInt(11n).is_irreducible()).toBe(true); // x^3 + x + 1
    expect(GF2X.fromBigInt(13n).is_irreducible()).toBe(true); // x^3 + x^2 + 1
    expect(GF2X.fromBigInt(15n).is_irreducible()).toBe(false); // x^3 + x^2 + x + 1 = (x+1)^3
  });

  test('is_irreducible edge cases', () => {
    expect(GF2X.zero().is_irreducible()).toBe(false);
    expect(GF2X.one().is_irreducible()).toBe(false);

    // Polynomials without constant term are divisible by x
    expect(GF2X.fromBigInt(6n).is_irreducible()).toBe(false); // x^2 + x = x(x+1)
    expect(GF2X.fromBigInt(10n).is_irreducible()).toBe(false); // x^3 + x = x(x^2+1)
  });

  test('degree 4 irreducibles', () => {
    // Degree 4 irreducibles: x^4 + x + 1, x^4 + x^3 + 1, x^4 + x^3 + x^2 + x + 1
    expect(GF2X.fromBigInt(19n).is_irreducible()).toBe(true); // x^4 + x + 1
    expect(GF2X.fromBigInt(25n).is_irreducible()).toBe(true); // x^4 + x^3 + 1
    expect(GF2X.fromBigInt(31n).is_irreducible()).toBe(true); // x^4 + x^3 + x^2 + x + 1

    // Not irreducible
    expect(GF2X.fromBigInt(21n).is_irreducible()).toBe(false); // x^4 + x^2 + 1 = (x^2+x+1)^2
  });
});

describe('GF2X derivative and reverse', () => {
  test('derivative', () => {
    // Derivative of x^n is n*x^(n-1), which in GF(2) is:
    // - x^(n-1) if n is odd
    // - 0 if n is even

    expect(GF2X.one().derivative().isZero()).toBe(true);
    expect(GF2X.x().derivative().isOne()).toBe(true); // d(x)/dx = 1

    // d(x^2)/dx = 2x = 0 in GF(2)
    expect(GF2X.fromBigInt(4n).derivative().isZero()).toBe(true);

    // d(x^3)/dx = 3x^2 = x^2 in GF(2)
    expect(GF2X.fromBigInt(8n).derivative().bits).toBe(4n);

    // d(x^3 + x + 1)/dx = x^2 + 1
    expect(GF2X.fromBigInt(11n).derivative().bits).toBe(5n);

    // d(x^4 + x^3 + x^2 + x + 1)/dx = x^2 + 1
    expect(GF2X.fromBigInt(31n).derivative().bits).toBe(5n);
  });

  test('reverse', () => {
    // Reverse of x^3 + x + 1 (11n) is x^3 + x^2 + 1 (13n)
    expect(GF2X.fromBigInt(11n).reverse().bits).toBe(13n);

    // Reverse of x^2 + 1 is 1 + x^2 (same)
    expect(GF2X.fromBigInt(5n).reverse().bits).toBe(5n);

    // Reverse of x + 1 is 1 + x (same)
    expect(GF2X.fromBigInt(3n).reverse().bits).toBe(3n);

    // Reverse of x^3 is 1
    expect(GF2X.fromBigInt(8n).reverse().bits).toBe(1n);

    // Reverse with explicit hi parameter
    // x + 1 (3n) padded to degree 3 and reversed: x^3 + x^2 (12n)
    expect(GF2X.fromBigInt(3n).reverse(3).bits).toBe(12n);
  });
});

describe('buildIrred functions', () => {
  test('buildIrred', () => {
    // Should return the lexicographically smallest irreducible.
    // NTL special-cases n = 1 and returns x, not x + 1
    // (`GF2XFactoring.cpp:481-484`: `if (n == 1) { SetX(f); return; }`),
    // and Sage's GF2X_BuildIrred_list(1) is therefore [0, 1].
    expect(buildIrred(1).bits).toBe(2n); // x
    expect(buildSparseIrred(1).bits).toBe(2n); // x
    expect(buildRandomIrred(1).bits).toBe(2n); // x
    expect(GF2X_BuildIrred_list(1)).toEqual([0, 1]);
    expect(GF2X_BuildSparseIrred_list(1)).toEqual([0, 1]);
    expect(buildIrred(2).bits).toBe(7n); // x^2 + x + 1
    expect(buildIrred(3).bits).toBe(11n); // x^3 + x + 1
    expect(buildIrred(4).bits).toBe(19n); // x^4 + x + 1

    // Verify all are irreducible
    for (let n = 1; n <= 8; n++) {
      const irr = buildIrred(n);
      expect(irr.is_irreducible()).toBe(true);
      expect(irr.degree()).toBe(n);
    }
  });

  test('buildSparseIrred', () => {
    // Should prefer trinomials, then pentanomials
    for (let n = 1; n <= 8; n++) {
      const irr = buildSparseIrred(n);
      expect(irr.is_irreducible()).toBe(true);
      expect(irr.degree()).toBe(n);
    }
  });

  test('buildRandomIrred', () => {
    set_random_seed(12345n);

    for (let n = 1; n <= 6; n++) {
      const irr = buildRandomIrred(n);
      expect(irr.is_irreducible()).toBe(true);
      expect(irr.degree()).toBe(n);
    }
  });

  test('GF2X_BuildIrred_list matches SageMath', () => {
    // From SageMath:
    // sage: from sage.rings.polynomial.polynomial_gf2x import GF2X_BuildIrred_list
    // sage: GF2X_BuildIrred_list(2)
    // [1, 1, 1]
    expect(GF2X_BuildIrred_list(2)).toEqual([1, 1, 1]);

    // sage: GF2X_BuildIrred_list(3)
    // [1, 1, 0, 1]
    expect(GF2X_BuildIrred_list(3)).toEqual([1, 1, 0, 1]);

    // sage: GF2X_BuildIrred_list(4)
    // [1, 1, 0, 0, 1]
    expect(GF2X_BuildIrred_list(4)).toEqual([1, 1, 0, 0, 1]);
  });

  test('GF2X_BuildSparseIrred_list matches for small degrees', () => {
    // For small degrees, sparse irred should equal canonical irred
    for (let n = 1; n <= 4; n++) {
      expect(GF2X_BuildSparseIrred_list(n)).toEqual(GF2X_BuildIrred_list(n));
    }
  });
});

describe('GF2X factorization', () => {
  test('squareFreeDecomp', () => {
    // (x + 1)^2 = x^2 + 1
    const p1 = GF2X.fromBigInt(5n);
    const decomp1 = squareFreeDecomp(p1);
    expect(decomp1.length).toBe(1);
    expect(decomp1[0]![0].bits).toBe(3n); // x + 1
    expect(decomp1[0]![1]).toBe(2);

    // Square-free polynomial
    const p2 = GF2X.fromBigInt(11n); // x^3 + x + 1 (irreducible)
    const decomp2 = squareFreeDecomp(p2);
    expect(decomp2.length).toBe(1);
    expect(decomp2[0]![0].bits).toBe(11n);
    expect(decomp2[0]![1]).toBe(1);
  });

  test('distinctDegreeFactorization', () => {
    // x^3 + x^2 + 1 is irreducible of degree 3
    const f = GF2X.fromBigInt(13n);
    const ddf = distinctDegreeFactorization(f);
    expect(ddf.length).toBe(1);
    expect(ddf[0]![1]).toBe(3); // degree 3

    // Product of two irreducibles of degree 2 and 3
    const irr2 = GF2X.fromBigInt(7n); // x^2 + x + 1
    const irr3 = GF2X.fromBigInt(11n); // x^3 + x + 1
    const prod = irr2.mul(irr3);
    const ddf2 = distinctDegreeFactorization(prod);
    expect(ddf2.length).toBe(2);
  });

  test('factor simple cases', () => {
    // Factor x^2 + 1 = (x + 1)^2
    const f1 = GF2X.fromBigInt(5n);
    const factors1 = factor(f1);
    expect(factors1.length).toBe(1);
    expect(factors1[0]![0].bits).toBe(3n); // x + 1
    expect(factors1[0]![1]).toBe(2);

    // Factor x^3 + x + 1 (irreducible)
    const f2 = GF2X.fromBigInt(11n);
    const factors2 = factor(f2);
    expect(factors2.length).toBe(1);
    expect(factors2[0]![0].bits).toBe(11n);
    expect(factors2[0]![1]).toBe(1);

    // Factor x^4 + x^2 + 1 = (x^2 + x + 1)^2
    const f3 = GF2X.fromBigInt(21n);
    const factors3 = factor(f3);
    expect(factors3.length).toBe(1);
    expect(factors3[0]![0].bits).toBe(7n); // x^2 + x + 1
    expect(factors3[0]![1]).toBe(2);
  });

  test('factor product of distinct irreducibles', () => {
    set_random_seed(42n);

    const irr2 = GF2X.fromBigInt(7n); // x^2 + x + 1
    const irr3 = GF2X.fromBigInt(11n); // x^3 + x + 1
    const prod = irr2.mul(irr3);

    const factors = factor(prod);
    expect(factors.length).toBe(2);

    // Reconstruct and verify
    let reconstructed = GF2X.one();
    for (const [f, m] of factors) {
      reconstructed = reconstructed.mul(f.pow(m));
    }
    expect(reconstructed.eq(prod)).toBe(true);
  });

  test('factor with multiplicities', () => {
    set_random_seed(123n);

    // (x + 1)^2 * (x^2 + x + 1)
    const p1 = GF2X.fromBigInt(3n); // x + 1
    const p2 = GF2X.fromBigInt(7n); // x^2 + x + 1
    const prod = p1.pow(2).mul(p2);

    const factors = factor(prod);

    // Verify reconstruction
    let reconstructed = GF2X.one();
    for (const [f, m] of factors) {
      reconstructed = reconstructed.mul(f.pow(m));
    }
    expect(reconstructed.eq(prod)).toBe(true);
  });
});

describe('GF2X_Ring', () => {
  test('ring operations', () => {
    const R = GF2X_Ring;

    expect(R.zero().isZero()).toBe(true);
    expect(R.one().isOne()).toBe(true);
    expect(R.gen().isX()).toBe(true);
    expect(R.is_field()).toBe(false);
    expect(R.is_finite()).toBe(false);
  });

  test('__call__ coercion', () => {
    const R = GF2X_Ring;

    expect(R.__call__(11n).bits).toBe(11n);
    expect(R.__call__(11).bits).toBe(11n);
    expect(R.__call__([1, 1, 0, 1]).bits).toBe(11n);

    const p = GF2X.fromBigInt(7n);
    expect(R.__call__(p)).toBe(p);
  });

  test('random_element', () => {
    set_random_seed(999n);
    const R = GF2X_Ring;

    // Sage's `random_element(d)` returns a polynomial of degree *exactly* d
    // (`polynomial_ring.py:1369`: `R.random_element(6).degree() == 6`).
    for (let i = 0; i < 20; i++) {
      expect(R.random_element(10).degree()).toBe(10);
    }
    // Degree -1 is the zero polynomial.
    expect(R.random_element(-1).isZero()).toBe(true);
    // A (min, max) range gives a degree in that range.
    for (let i = 0; i < 20; i++) {
      const d = R.random_element([0, 4]).degree();
      expect(d).toBeGreaterThanOrEqual(0);
      expect(d).toBeLessThanOrEqual(4);
    }
    // Monic polynomials are never zero.
    for (let i = 0; i < 20; i++) {
      const r = R.random_element([-1, 3], true);
      expect(r.isZero()).toBe(false);
    }
    expect(() => R.random_element(-5)).toThrow('must be at least -1');
    expect(() => R.random_element([5, 4])).toThrow(
      'minimum degree must be less or equal than maximum degree'
    );
  });
});

describe('GF2X edge cases and large polynomials', () => {
  test('operations with large polynomials', () => {
    const large1 = GF2X.fromBigInt((1n << 100n) | (1n << 50n) | 1n);
    const large2 = GF2X.fromBigInt((1n << 80n) | (1n << 40n) | 1n);

    // Addition
    const sum = large1.add(large2);
    expect(sum.add(large2).eq(large1)).toBe(true);

    // Multiplication
    const prod = large1.mul(large2);
    expect(prod.degree()).toBe(180);

    // Division
    const [q, r] = prod.divRem(large2);
    expect(q.mul(large2).add(r).eq(prod)).toBe(true);
  });

  test('irreducibility of AES polynomial', () => {
    // AES uses x^8 + x^4 + x^3 + x + 1 which is irreducible over GF(2)
    const aes = GF2X.fromBigInt(0x11bn); // 100011011 in binary
    expect(aes.is_irreducible()).toBe(true);
    expect(aes.degree()).toBe(8);
  });

  test('equality with bigint', () => {
    const p = GF2X.fromBigInt(11n);
    expect(p.eq(11n)).toBe(true);
    expect(p.eq(12n)).toBe(false);
  });
});
