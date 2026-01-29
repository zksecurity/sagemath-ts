/**
 * Unit tests for polynomial module
 */
import { describe, expect, test } from 'bun:test';
import { GF2, GF2Element } from '../finite_rings/gf2.js';
import { Polynomial, PolynomialRing, PolynomialRingConstructor, QuotientRing } from './index.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';

/**
 * Helper function to compute GCD of two bigints.
 */
function bigintGcd(a: bigint, b: bigint): bigint {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) {
    const t = b;
    b = a % b;
    a = t;
  }
  return a;
}

/**
 * A simple wrapper around bigint that implements RingElement.
 * Used for testing polynomials over integers.
 */
class IntegerElement implements RingElement {
  readonly value: bigint;

  constructor(value: bigint | number) {
    this.value = typeof value === 'number' ? BigInt(value) : value;
  }

  add(other: IntegerElement): IntegerElement {
    return new IntegerElement(this.value + other.value);
  }

  sub(other: IntegerElement): IntegerElement {
    return new IntegerElement(this.value - other.value);
  }

  mul(other: IntegerElement): IntegerElement {
    return new IntegerElement(this.value * other.value);
  }

  div(other: IntegerElement): IntegerElement {
    return new IntegerElement(this.value / other.value);
  }

  neg(): IntegerElement {
    return new IntegerElement(-this.value);
  }

  eq(other: IntegerElement | number): boolean {
    if (typeof other === 'number') {
      return this.value === BigInt(other);
    }
    return this.value === other.value;
  }

  isZero(): boolean {
    return this.value === 0n;
  }

  gcd(other: IntegerElement): IntegerElement {
    return new IntegerElement(bigintGcd(this.value, other.value));
  }

  toString(): string {
    return this.value.toString();
  }
}

/**
 * Integer ring that produces IntegerElement objects.
 */
const IntegerRing: CoefficientRing<IntegerElement> = {
  zero(): IntegerElement {
    return new IntegerElement(0n);
  },
  one(): IntegerElement {
    return new IntegerElement(1n);
  },
  __call__(x: unknown): IntegerElement {
    if (x instanceof IntegerElement) {
      return x;
    }
    if (typeof x === 'bigint') {
      return new IntegerElement(x);
    }
    if (typeof x === 'number') {
      return new IntegerElement(BigInt(x));
    }
    throw new Error(`Cannot convert ${typeof x} to IntegerElement`);
  },
  is_field(): boolean {
    return false;
  },
};

describe('Polynomial over integers', () => {
  const [R, x] = PolynomialRingConstructor(IntegerRing, 'x');

  test('zero and one', () => {
    expect(R.zero().isZero()).toBe(true);
    expect(R.one().isConstant()).toBe(true);
    expect(R.one().leading_coefficient().eq(1)).toBe(true);
  });

  test('generator', () => {
    expect(x.degree()).toBe(1);
    expect(x.getCoeff(0).isZero()).toBe(true);
    expect(x.getCoeff(1).eq(1)).toBe(true);
  });

  test('addition', () => {
    // (x + 1) + (x^2 + 2x) = x^2 + 3x + 1
    const p1 = x.add(R.one());
    const p2 = x.pow(2).add(x.scalar_mul(new IntegerElement(2)));
    const sum = p1.add(p2);

    expect(sum.degree()).toBe(2);
    expect(sum.getCoeff(0).eq(1)).toBe(true);
    expect(sum.getCoeff(1).eq(3)).toBe(true);
    expect(sum.getCoeff(2).eq(1)).toBe(true);
  });

  test('subtraction', () => {
    // (x^2 + 2x + 1) - (x + 1) = x^2 + x
    const p1 = x
      .pow(2)
      .add(x.scalar_mul(new IntegerElement(2)))
      .add(R.one());
    const p2 = x.add(R.one());
    const diff = p1.sub(p2);

    expect(diff.degree()).toBe(2);
    expect(diff.getCoeff(0).isZero()).toBe(true);
    expect(diff.getCoeff(1).eq(1)).toBe(true);
    expect(diff.getCoeff(2).eq(1)).toBe(true);
  });

  test('multiplication', () => {
    // (x + 1) * (x - 1) = x^2 - 1
    const p1 = x.add(R.one());
    const minusOne = R.__call__(new IntegerElement(-1));
    const p2 = x.add(minusOne);
    const prod = p1.mul(p2);

    expect(prod.degree()).toBe(2);
    expect(prod.getCoeff(0).eq(-1)).toBe(true);
    expect(prod.getCoeff(1).isZero()).toBe(true);
    expect(prod.getCoeff(2).eq(1)).toBe(true);
  });

  test('negation', () => {
    const p = x.add(R.one());
    const neg = p.neg();

    expect(neg.getCoeff(0).eq(-1)).toBe(true);
    expect(neg.getCoeff(1).eq(-1)).toBe(true);
  });

  test('power', () => {
    // (x + 1)^3 = x^3 + 3x^2 + 3x + 1
    const p = x.add(R.one());
    const cubed = p.pow(3);

    expect(cubed.degree()).toBe(3);
    expect(cubed.getCoeff(0).eq(1)).toBe(true);
    expect(cubed.getCoeff(1).eq(3)).toBe(true);
    expect(cubed.getCoeff(2).eq(3)).toBe(true);
    expect(cubed.getCoeff(3).eq(1)).toBe(true);
  });

  test('degree', () => {
    expect(R.zero().degree()).toBe(-1);
    expect(R.one().degree()).toBe(0);
    expect(x.degree()).toBe(1);
    expect(x.pow(5).degree()).toBe(5);
  });

  test('leading_coefficient', () => {
    const p = x
      .pow(3)
      .add(x.scalar_mul(new IntegerElement(2)))
      .add(R.one());
    expect(p.leading_coefficient().eq(1)).toBe(true);

    const p2 = x.scalar_mul(new IntegerElement(5));
    expect(p2.leading_coefficient().eq(5)).toBe(true);
  });

  test('evaluate', () => {
    // p(x) = x^2 + 2x + 1
    const p = x
      .pow(2)
      .add(x.scalar_mul(new IntegerElement(2)))
      .add(R.one());

    // p(0) = 1
    expect(p.evaluate(new IntegerElement(0)).eq(1)).toBe(true);

    // p(1) = 1 + 2 + 1 = 4
    expect(p.evaluate(new IntegerElement(1)).eq(4)).toBe(true);

    // p(2) = 4 + 4 + 1 = 9
    expect(p.evaluate(new IntegerElement(2)).eq(9)).toBe(true);

    // p(-1) = 1 - 2 + 1 = 0
    expect(p.evaluate(new IntegerElement(-1)).isZero()).toBe(true);
  });

  test('is_monic', () => {
    expect(x.is_monic()).toBe(true);
    expect(x.pow(3).add(x).add(R.one()).is_monic()).toBe(true);

    const p = x.scalar_mul(new IntegerElement(2));
    expect(p.is_monic()).toBe(false);

    expect(R.zero().is_monic()).toBe(false);
  });

  test('equality', () => {
    const p1 = x.pow(2).add(x).add(R.one());
    const p2 = x.pow(2).add(x).add(R.one());
    const p3 = x.pow(2).add(x);

    expect(p1.eq(p2)).toBe(true);
    expect(p1.eq(p3)).toBe(false);
    expect(R.zero().eq(R.zero())).toBe(true);
  });

  test('toString', () => {
    const p = x.pow(2).add(x).add(R.one());
    const str = p.toString();
    expect(str).toContain('x^2');
    expect(str).toContain('x');
  });

  test('content', () => {
    // f = 6x^2 + 4x + 2
    const f = x
      .pow(2)
      .scalar_mul(new IntegerElement(6))
      .add(x.scalar_mul(new IntegerElement(4)))
      .add(R.__call__(new IntegerElement(2)));

    // content = gcd(6, 4, 2) = 2
    expect(f.content().eq(2)).toBe(true);
  });

  test('content of monic polynomial', () => {
    // f = x^2 + 2x + 1
    const f = x.pow(2).add(x.scalar_mul(new IntegerElement(2))).add(R.one());

    // content = gcd(1, 2, 1) = 1
    expect(f.content().eq(1)).toBe(true);
  });

  test('content of zero polynomial', () => {
    expect(R.zero().content().eq(0)).toBe(true);
  });

  test('primitive_part', () => {
    // f = 6x^2 + 4x + 2
    const f = x
      .pow(2)
      .scalar_mul(new IntegerElement(6))
      .add(x.scalar_mul(new IntegerElement(4)))
      .add(R.__call__(new IntegerElement(2)));

    // primitive_part = 3x^2 + 2x + 1
    const pp = f.primitive_part();
    expect(pp.getCoeff(0).eq(1)).toBe(true);
    expect(pp.getCoeff(1).eq(2)).toBe(true);
    expect(pp.getCoeff(2).eq(3)).toBe(true);
  });

  test('primitive_part of monic polynomial', () => {
    // f = x^2 + 2x + 1, content = 1
    const f = x.pow(2).add(x.scalar_mul(new IntegerElement(2))).add(R.one());

    // primitive_part should be the same polynomial
    expect(f.primitive_part().eq(f)).toBe(true);
  });

  test('shift positive', () => {
    // f = x^2 + 2x + 4
    const f = x
      .pow(2)
      .add(x.scalar_mul(new IntegerElement(2)))
      .add(R.__call__(new IntegerElement(4)));

    // shift(2) -> x^4 + 2x^3 + 4x^2
    const shifted = f.shift(2);
    expect(shifted.degree()).toBe(4);
    expect(shifted.getCoeff(0).eq(0)).toBe(true);
    expect(shifted.getCoeff(1).eq(0)).toBe(true);
    expect(shifted.getCoeff(2).eq(4)).toBe(true);
    expect(shifted.getCoeff(3).eq(2)).toBe(true);
    expect(shifted.getCoeff(4).eq(1)).toBe(true);
  });

  test('shift negative', () => {
    // f = x^2 + 2x + 4
    const f = x
      .pow(2)
      .add(x.scalar_mul(new IntegerElement(2)))
      .add(R.__call__(new IntegerElement(4)));

    // shift(-1) -> x + 2
    const shifted = f.shift(-1);
    expect(shifted.degree()).toBe(1);
    expect(shifted.getCoeff(0).eq(2)).toBe(true);
    expect(shifted.getCoeff(1).eq(1)).toBe(true);
  });
});

describe('Polynomial over GF(2)', () => {
  const [R, x] = PolynomialRingConstructor(GF2, 'x');

  test('zero and one', () => {
    expect(R.zero().isZero()).toBe(true);
    expect(R.one().isConstant()).toBe(true);
  });

  test('addition (same as XOR)', () => {
    // In GF(2), 1 + 1 = 0
    const p1 = x.add(R.one()); // x + 1
    const p2 = x.add(R.one()); // x + 1
    const sum = p1.add(p2); // (x + 1) + (x + 1) = 0 in GF(2)

    expect(sum.isZero()).toBe(true);
  });

  test('subtraction equals addition', () => {
    // In GF(2), subtraction is the same as addition
    const p1 = x.pow(2).add(x);
    const p2 = x.add(R.one());
    const diff = p1.sub(p2);
    const sum = p1.add(p2);

    expect(diff.eq(sum)).toBe(true);
  });

  test('multiplication', () => {
    // (x + 1) * (x + 1) = x^2 + 2x + 1 = x^2 + 1 in GF(2)
    const p = x.add(R.one());
    const prod = p.mul(p);

    expect(prod.degree()).toBe(2);
    expect(prod.getCoeff(0).eq(1)).toBe(true);
    expect(prod.getCoeff(1).isZero()).toBe(true); // 2x = 0 in GF(2)
    expect(prod.getCoeff(2).eq(1)).toBe(true);
  });

  test('power', () => {
    // (x + 1)^2 = x^2 + 1 in GF(2) (Freshman's dream)
    const p = x.add(R.one());
    const squared = p.pow(2);

    expect(squared.degree()).toBe(2);
    expect(squared.getCoeff(0).eq(1)).toBe(true);
    expect(squared.getCoeff(1).isZero()).toBe(true);
    expect(squared.getCoeff(2).eq(1)).toBe(true);
  });

  test('evaluate', () => {
    // p(x) = x^2 + x + 1
    const p = x.pow(2).add(x).add(R.one());

    // p(0) = 0 + 0 + 1 = 1
    expect(p.evaluate(GF2.zero()).eq(1)).toBe(true);

    // p(1) = 1 + 1 + 1 = 1 in GF(2)
    expect(p.evaluate(GF2.one()).eq(1)).toBe(true);
  });

  test('quo_rem basic', () => {
    // (x^2 + 1) / (x + 1) = (x + 1) remainder 0 in GF(2)
    // because (x + 1)^2 = x^2 + 1 in GF(2)
    const dividend = x.pow(2).add(R.one());
    const divisor = x.add(R.one());
    const [q, r] = dividend.quo_rem(divisor);

    expect(q.eq(divisor)).toBe(true);
    expect(r.isZero()).toBe(true);
  });

  test('quo_rem with remainder', () => {
    // x^3 / (x^2 + x + 1) = x remainder (x^2 + x)
    // x * (x^2 + x + 1) = x^3 + x^2 + x
    // x^3 - (x^3 + x^2 + x) = x^2 + x (in GF(2), - is +)
    const dividend = x.pow(3);
    const divisor = x.pow(2).add(x).add(R.one());
    const [q, r] = dividend.quo_rem(divisor);

    // Verify: q * divisor + r = dividend
    const reconstructed = q.mul(divisor).add(r);
    expect(reconstructed.eq(dividend)).toBe(true);
    expect(r.degree()).toBeLessThan(divisor.degree());
  });

  test('mod', () => {
    // x^3 mod (x^2 + x + 1)
    const dividend = x.pow(3);
    const modulus = x.pow(2).add(x).add(R.one());
    const r = dividend.mod(modulus);

    expect(r.degree()).toBeLessThan(modulus.degree());
  });
});

describe('QuotientRing (GF(2^n) construction)', () => {
  // Create GF(4) = GF(2)[x] / (x^2 + x + 1)
  const [R, x] = PolynomialRingConstructor(GF2, 'a');
  const modulus = x.pow(2).add(x).add(R.one()); // x^2 + x + 1 is irreducible over GF(2)
  const GF4 = new QuotientRing(R, modulus);

  test('creation', () => {
    expect(GF4.degree).toBe(2);
  });

  test('zero and one', () => {
    expect(GF4.zero().isZero()).toBe(true);
    expect(GF4.one().eq(1)).toBe(true);
  });

  test('generator', () => {
    const a = GF4.gen();
    expect(a.lift.degree()).toBe(1);
  });

  test('addition', () => {
    const a = GF4.gen();
    const one = GF4.one();

    // a + 1
    const sum = a.add(one);
    expect(sum.lift.getCoeff(0).eq(1)).toBe(true);
    expect(sum.lift.getCoeff(1).eq(1)).toBe(true);

    // a + a = 0 in characteristic 2
    const zero = a.add(a);
    expect(zero.isZero()).toBe(true);
  });

  test('multiplication with reduction', () => {
    const a = GF4.gen();

    // a^2 = a + 1 (since a^2 + a + 1 = 0)
    const aSquared = a.mul(a);

    // Should be reduced modulo x^2 + x + 1
    expect(aSquared.lift.degree()).toBeLessThan(2);

    // a^2 + a + 1 = 0, so a^2 = -a - 1 = a + 1 (in char 2)
    const expected = a.add(GF4.one());
    expect(aSquared.eq(expected)).toBe(true);
  });

  test('inverse', () => {
    const a = GF4.gen();

    // Find inverse of a
    const aInv = a.inv();

    // a * a^(-1) = 1
    const product = a.mul(aInv);
    expect(product.eq(1)).toBe(true);
  });

  test('division', () => {
    const a = GF4.gen();
    const one = GF4.one();
    const aPlus1 = a.add(one);

    // (a + 1) / a should give some element
    const quotient = aPlus1.div(a);

    // Verify: quotient * a = a + 1
    const reconstructed = quotient.mul(a);
    expect(reconstructed.eq(aPlus1)).toBe(true);
  });

  test('power', () => {
    const a = GF4.gen();

    // a^0 = 1
    expect(a.pow(0).eq(1)).toBe(true);

    // a^1 = a
    expect(a.pow(1).eq(a)).toBe(true);

    // a^3 = 1 (since GF(4)* has order 3)
    const aCubed = a.pow(3);
    expect(aCubed.eq(1)).toBe(true);
  });

  test('cardinality', () => {
    // GF(4) has 4 elements
    expect(GF4.cardinality()).toBe(4n);
  });

  test('iteration over elements', () => {
    const elements = [...GF4];
    expect(elements.length).toBe(4);

    // Check all elements are distinct
    for (let i = 0; i < elements.length; i++) {
      for (let j = i + 1; j < elements.length; j++) {
        expect(elements[i]!.eq(elements[j]!)).toBe(false);
      }
    }
  });

  test('field properties', () => {
    const elements = [...GF4];

    // Every non-zero element should have an inverse
    for (const elem of elements) {
      if (!elem.isZero()) {
        const inv = elem.inv();
        expect(elem.mul(inv).eq(1)).toBe(true);
      }
    }
  });
});

describe('QuotientRing - GF(8) construction', () => {
  // Create GF(8) = GF(2)[x] / (x^3 + x + 1)
  const [R, x] = PolynomialRingConstructor(GF2, 'b');
  const modulus = x.pow(3).add(x).add(R.one()); // x^3 + x + 1 is irreducible over GF(2)
  const GF8 = new QuotientRing(R, modulus);

  test('cardinality', () => {
    expect(GF8.cardinality()).toBe(8n);
  });

  test('generator order', () => {
    const b = GF8.gen();

    // b^7 = 1 (since GF(8)* has order 7)
    const bToThe7 = b.pow(7);
    expect(bToThe7.eq(1)).toBe(true);

    // b^i != 1 for 1 <= i < 7
    for (let i = 1; i < 7; i++) {
      expect(b.pow(i).eq(1)).toBe(false);
    }
  });

  test('multiplication table consistency', () => {
    const b = GF8.gen();
    const elements = [...GF8];

    // Check associativity for a few random elements
    const a1 = b;
    const a2 = b.add(GF8.one());
    const a3 = b.pow(2);

    // (a1 * a2) * a3 = a1 * (a2 * a3)
    const left = a1.mul(a2).mul(a3);
    const right = a1.mul(a2.mul(a3));
    expect(left.eq(right)).toBe(true);
  });
});

describe('Polynomial edge cases', () => {
  const [R, x] = PolynomialRingConstructor(GF2, 'x');

  test('zero polynomial operations', () => {
    const zero = R.zero();
    const one = R.one();

    expect(zero.add(one).eq(one)).toBe(true);
    expect(zero.mul(x).isZero()).toBe(true);
    expect(zero.pow(0).eq(one)).toBe(true);
    expect(zero.neg().isZero()).toBe(true);
  });

  test('constant polynomial operations', () => {
    const one = R.one();

    expect(one.mul(x).eq(x)).toBe(true);
    expect(one.pow(100).eq(one)).toBe(true);
    expect(one.add(one).isZero()).toBe(true); // 1 + 1 = 0 in GF(2)
  });

  test('fromTerms', () => {
    // Create x^3 + x using fromTerms
    const p = R.fromTerms([
      [GF2.one(), 3],
      [GF2.one(), 1],
    ]);

    expect(p.degree()).toBe(3);
    expect(p.getCoeff(3).eq(1)).toBe(true);
    expect(p.getCoeff(2).isZero()).toBe(true);
    expect(p.getCoeff(1).eq(1)).toBe(true);
    expect(p.getCoeff(0).isZero()).toBe(true);
  });

  test('monomial', () => {
    const x5 = R.monomial(5);

    expect(x5.degree()).toBe(5);
    expect(x5.leading_coefficient().eq(1)).toBe(true);
    expect(x5.getCoeff(0).isZero()).toBe(true);
    expect(x5.getCoeff(4).isZero()).toBe(true);
    expect(x5.getCoeff(5).eq(1)).toBe(true);
  });

  test('scalar_mul', () => {
    const p = x.pow(2).add(x).add(R.one());

    // Multiply by 0
    expect(p.scalar_mul(GF2.zero()).isZero()).toBe(true);

    // Multiply by 1
    expect(p.scalar_mul(GF2.one()).eq(p)).toBe(true);
  });
});

describe('Polynomial shift, truncate, reverse', () => {
  const [R, x] = PolynomialRingConstructor(GF2, 'x');

  test('shift positive', () => {
    // p = x^2 + 1 (in GF(2))
    const p = x.pow(2).add(R.one());

    // shift(2) -> x^4 + x^2
    const shifted = p.shift(2);
    expect(shifted.degree()).toBe(4);
    expect(shifted.getCoeff(0).isZero()).toBe(true);
    expect(shifted.getCoeff(1).isZero()).toBe(true);
    expect(shifted.getCoeff(2).eq(1)).toBe(true);
    expect(shifted.getCoeff(3).isZero()).toBe(true);
    expect(shifted.getCoeff(4).eq(1)).toBe(true);
  });

  test('shift negative', () => {
    // p = x^3 + x^2 + x + 1 (coeffs: 1, 1, 1, 1)
    const p = x.pow(3).add(x.pow(2)).add(x).add(R.one());

    // shift(-1) -> x^2 + x + 1 (drop lowest coefficient)
    const shifted = p.shift(-1);
    expect(shifted.degree()).toBe(2);
    expect(shifted.getCoeff(0).eq(1)).toBe(true);
    expect(shifted.getCoeff(1).eq(1)).toBe(true);
    expect(shifted.getCoeff(2).eq(1)).toBe(true);
  });

  test('shift zero', () => {
    const p = x.pow(2).add(R.one());
    expect(p.shift(0).eq(p)).toBe(true);
  });

  test('shift zero polynomial', () => {
    const zero = R.zero();
    expect(zero.shift(5).isZero()).toBe(true);
    expect(zero.shift(-5).isZero()).toBe(true);
  });

  test('shift negative too large', () => {
    const p = x.pow(2).add(R.one()); // degree 2
    // shift(-5) should result in zero (drop more than we have)
    expect(p.shift(-5).isZero()).toBe(true);
  });

  test('truncate', () => {
    // p = x^3 + x^2 + x + 1
    const p = x.pow(3).add(x.pow(2)).add(x).add(R.one());

    // truncate(2) -> x + 1 (keep only terms of degree < 2)
    const truncated = p.truncate(2);
    expect(truncated.degree()).toBe(1);
    expect(truncated.getCoeff(0).eq(1)).toBe(true);
    expect(truncated.getCoeff(1).eq(1)).toBe(true);
  });

  test('truncate to zero', () => {
    const p = x.pow(2).add(R.one());
    expect(p.truncate(0).isZero()).toBe(true);
  });

  test('truncate larger than degree', () => {
    const p = x.pow(2).add(R.one());
    // truncate(10) keeps all terms (degree < 10)
    expect(p.truncate(10).eq(p)).toBe(true);
  });

  test('reverse basic', () => {
    // p = x^2 + 1 (coeffs: [1, 0, 1])
    const p = x.pow(2).add(R.one());

    // reverse -> 1 + 0*x + x^2 = x^2 + 1 (same polynomial!)
    const reversed = p.reverse();
    expect(reversed.eq(p)).toBe(true);
  });

  test('reverse asymmetric', () => {
    // p = x^2 + x (coeffs: [0, 1, 1])
    const p = x.pow(2).add(x);

    // reverse -> coeffs become [1, 1, 0] = x + 1 (after removing trailing zero)
    const reversed = p.reverse();
    expect(reversed.degree()).toBe(1);
    expect(reversed.getCoeff(0).eq(1)).toBe(true);
    expect(reversed.getCoeff(1).eq(1)).toBe(true);
  });

  test('reverse with degree argument', () => {
    // p = x^2 + 1 (coeffs: [1, 0, 1], degree 2)
    const p = x.pow(2).add(R.one());

    // reverse(degree=4) -> pad to degree 4, then reverse
    // [1, 0, 1, 0, 0] -> [0, 0, 1, 0, 1] = x^4 + x^2
    const reversed = p.reverse(4);
    expect(reversed.degree()).toBe(4);
    expect(reversed.getCoeff(0).isZero()).toBe(true);
    expect(reversed.getCoeff(2).eq(1)).toBe(true);
    expect(reversed.getCoeff(4).eq(1)).toBe(true);
  });
});

describe('Polynomial resultant and discriminant over GF(2)', () => {
  const [R, x] = PolynomialRingConstructor(GF2, 'x');

  test('resultant of coprime polynomials', () => {
    // f = x + 1, g = x (no common roots in algebraic closure)
    const f = x.add(R.one());
    const g = x;

    // Resultant should be non-zero
    const res = f.resultant(g);
    expect(res.eq(1)).toBe(true); // f(0) = 1
  });

  test('resultant of polynomials with common root', () => {
    // f = x^2 + x = x(x+1), g = x + 1
    // They share root x = 1 (in GF(2), 1+1=0)
    const f = x.pow(2).add(x);
    const g = x.add(R.one());

    // Resultant should be zero
    const res = f.resultant(g);
    expect(res.isZero()).toBe(true);
  });

  test('resultant with zero polynomial', () => {
    const f = x.add(R.one());
    const zero = R.zero();

    expect(f.resultant(zero).isZero()).toBe(true);
    expect(zero.resultant(f).isZero()).toBe(true);
  });

  test('discriminant of x^2 + x + 1', () => {
    // x^2 + x + 1 is irreducible over GF(2), discriminant should be non-zero
    const f = x.pow(2).add(x).add(R.one());
    const disc = f.discriminant();
    // In GF(2): disc = (-1)^1 * Res(f, f') / a_n
    // f' = 2x + 1 = 1 (in GF(2))
    // So disc = Res(x^2+x+1, 1) = 1
    expect(disc.eq(1)).toBe(true);
  });
});
