/**
 * Tests for quotient_ring.ts.
 *
 * Expected values are SageMath's
 * (`polynomial_quotient_ring.py:767` / `:1006` doctests):
 *
 *     sage: R.<z> = ZZ[]; R.quo(z^2 - 2).is_field()
 *     False
 *     sage: R.quo(1).cardinality()
 *     1
 *     sage: R.quo(z^3 - 2).cardinality()
 *     +Infinity
 */
import { describe, expect, test } from 'bun:test';
import { ZeroDivisionError } from '../../errors.js';
import { FiniteFieldPrime } from '../finite_rings/finite_field_prime.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';
import { PolynomialRingConstructor } from './polynomial_ring.js';
import { QuotientRing } from './quotient_ring.js';

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

class IntegerElement implements RingElement {
  readonly value: bigint;
  constructor(value: bigint | number) {
    this.value = typeof value === 'number' ? BigInt(value) : value;
  }
  add(o: IntegerElement): IntegerElement {
    return new IntegerElement(this.value + o.value);
  }
  sub(o: IntegerElement): IntegerElement {
    return new IntegerElement(this.value - o.value);
  }
  mul(o: IntegerElement): IntegerElement {
    return new IntegerElement(this.value * o.value);
  }
  div(o: IntegerElement): IntegerElement {
    return new IntegerElement(this.value / o.value);
  }
  neg(): IntegerElement {
    return new IntegerElement(-this.value);
  }
  eq(o: IntegerElement | number): boolean {
    return typeof o === 'number' ? this.value === BigInt(o) : this.value === o.value;
  }
  isZero(): boolean {
    return this.value === 0n;
  }
  gcd(o: IntegerElement): IntegerElement {
    return new IntegerElement(bigintGcd(this.value, o.value));
  }
  toString(): string {
    return this.value.toString();
  }
}

const ZZ: CoefficientRing<IntegerElement> = {
  zero: () => new IntegerElement(0n),
  one: () => new IntegerElement(1n),
  __call__(v: unknown): IntegerElement {
    if (v instanceof IntegerElement) return v;
    if (typeof v === 'bigint') return new IntegerElement(v);
    if (typeof v === 'number') return new IntegerElement(BigInt(v));
    throw new Error('cannot convert to IntegerElement');
  },
  is_field: () => false,
  toString: () => 'Integer Ring',
};

const F5 = new FiniteFieldPrime(5n);
const [R5, x] = PolynomialRingConstructor(F5, 'x');
const [RZ, z] = PolynomialRingConstructor(ZZ, 'z');

describe('is_field (H13)', () => {
  test('GF(5)[x]/(x^2-1) is not a field', () => {
    const Q = new QuotientRing(R5, x.pow(2).sub(R5.one()));
    expect(Q.modulus.is_irreducible()).toBe(false);
    expect(Q.is_field()).toBe(false);
  });

  test('GF(5)[x]/(x^2+2) is a field with 25 elements', () => {
    const Q = new QuotientRing(R5, x.pow(2).add(R5.__call__(F5.__call__(2))));
    expect(Q.is_field()).toBe(true);
    expect(Q.cardinality()).toBe(25n);
  });

  test('ZZ[z]/(z^2-2) is not a field (base ring is not a field)', () => {
    const Q = new QuotientRing(RZ, z.pow(2).sub(RZ.__call__(new IntegerElement(2))));
    expect(Q.is_field()).toBe(false);
  });
});

describe('cardinality and is_finite (L16)', () => {
  test('degree-0 modulus is allowed and gives the zero ring', () => {
    const Q = new QuotientRing(RZ, RZ.one());
    expect(Q.cardinality()).toBe(1n);
    expect(Q.is_finite()).toBe(true);
    // Everything reduces to zero
    expect(Q.__call__(z).isZero()).toBe(true);
    expect(Q.one().isZero()).toBe(true);
  });

  test('infinite base ring gives Infinity instead of throwing', () => {
    const Q = new QuotientRing(RZ, z.pow(3).sub(RZ.__call__(new IntegerElement(2))));
    expect(Q.cardinality()).toBe(Number.POSITIVE_INFINITY);
    expect(Q.order()).toBe(Number.POSITIVE_INFINITY);
    expect(Q.is_finite()).toBe(false);
  });

  test('zero modulus is still rejected', () => {
    expect(() => new QuotientRing(RZ, RZ.zero())).toThrow('modulus cannot be zero');
  });
});

describe('inverse of a non-unit (L18)', () => {
  test('raises ZeroDivisionError with Sage message', () => {
    const Q = new QuotientRing(R5, x.pow(2).sub(R5.one()));
    const zeroDivisor = Q.__call__(x.sub(R5.one()));
    expect(() => zeroDivisor.inv()).toThrow(ZeroDivisionError);
    expect(() => zeroDivisor.inv()).toThrow('of quotient polynomial ring not invertible');
  });

  test('units are still invertible', () => {
    const Q = new QuotientRing(R5, x.pow(2).add(R5.__call__(F5.__call__(2))));
    const g = Q.gen();
    expect(g.mul(g.inv()).isOne()).toBe(true);
  });

  test('inverting zero raises ZeroDivisionError', () => {
    const Q = new QuotientRing(R5, x.pow(2).add(R5.__call__(F5.__call__(2))));
    expect(() => Q.zero().inv()).toThrow(ZeroDivisionError);
  });
});
