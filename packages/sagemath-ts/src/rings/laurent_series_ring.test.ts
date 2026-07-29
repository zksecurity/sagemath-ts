/**
 * Unit tests for the Laurent series ring module.
 *
 * Every value assertion below is copied verbatim from a SageMath doctest in
 * `reference/sage/src/sage/rings/laurent_series_ring_element.pyx` (or
 * `laurent_series_ring.py`); the source line is quoted next to each block.
 */
import { describe, expect, test } from 'bun:test';
import { EllipticCurve } from '../schemes/elliptic_curves/constructor.js';
import { GF as GFring } from './finite_rings/finite_field_constructor.js';
import { LaurentSeriesElement, LaurentSeriesRing } from './laurent_series_ring.js';
import type { CoefficientRing, RingElement } from './power_series_ring.js';
import { Rational } from './rational.js';
import { QQ as QQring } from './rational_field.js';

// ---------------------------------------------------------------------------
// Coefficient rings used by the doctests: QQ, ZZ and GF(p).
// ---------------------------------------------------------------------------

class RationalElement implements RingElement {
  readonly value: Rational;

  constructor(value: Rational | bigint | number | string) {
    if (value instanceof Rational) {
      this.value = value;
    } else if (typeof value === 'bigint') {
      this.value = new Rational(value, 1n);
    } else if (typeof value === 'number') {
      this.value = new Rational(BigInt(value), 1n);
    } else {
      this.value = Rational.fromString(value);
    }
  }

  add(other: RationalElement): RationalElement {
    return new RationalElement(this.value.add(other.value));
  }
  sub(other: RationalElement): RationalElement {
    return new RationalElement(this.value.sub(other.value));
  }
  mul(other: RationalElement): RationalElement {
    return new RationalElement(this.value.mul(other.value));
  }
  div(other: RationalElement): RationalElement {
    return new RationalElement(this.value.div(other.value));
  }
  neg(): RationalElement {
    return new RationalElement(this.value.neg());
  }
  inv(): RationalElement {
    return new RationalElement(this.value.inv());
  }
  eq(other: RationalElement | number | bigint): boolean {
    if (typeof other === 'number' || typeof other === 'bigint') {
      return this.value.eq(BigInt(other));
    }
    return this.value.eq(other.value);
  }
  isZero(): boolean {
    return this.value.eq(0n);
  }
  isOne(): boolean {
    return this.value.eq(1n);
  }
  isUnit(): boolean {
    return !this.isZero();
  }
  is_square(): boolean {
    const n = this.value.numerator;
    const d = this.value.denominator;
    if (n < 0n) return false;
    return isPerfectSquare(n) && isPerfectSquare(d);
  }
  sqrt(): RationalElement {
    const n = isqrt(this.value.numerator);
    const d = isqrt(this.value.denominator);
    if (n * n !== this.value.numerator || d * d !== this.value.denominator) {
      throw new Error(`unable to take the square root of ${this}`);
    }
    return new RationalElement(new Rational(n, d));
  }
  toString(): string {
    return this.value.toString();
  }
}

function isqrt(n: bigint): bigint {
  if (n < 0n) throw new Error('negative');
  if (n < 2n) return n;
  let x = n;
  let y = (x + 1n) / 2n;
  while (y < x) {
    x = y;
    y = (x + n / x) / 2n;
  }
  return x;
}

function isPerfectSquare(n: bigint): boolean {
  if (n < 0n) return false;
  const r = isqrt(n);
  return r * r === n;
}

const QQ: CoefficientRing<RationalElement> = {
  zero: () => new RationalElement(0n),
  one: () => new RationalElement(1n),
  __call__(x: unknown): RationalElement {
    if (x instanceof RationalElement) return x;
    if (x instanceof Rational) return new RationalElement(x);
    if (typeof x === 'bigint') return new RationalElement(x);
    if (typeof x === 'number') return new RationalElement(BigInt(x));
    if (typeof x === 'string') return new RationalElement(x);
    throw new Error(`Cannot convert ${typeof x} to RationalElement`);
  },
  is_field: () => true,
  characteristic: () => 0n,
};

/** The integers; division that is not exact throws, as in SageMath. */
class IntegerElement implements RingElement {
  readonly value: bigint;

  constructor(value: bigint | number) {
    this.value = typeof value === 'bigint' ? value : BigInt(value);
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
    if (other.value === 0n) throw new Error('division by zero');
    if (this.value % other.value !== 0n) {
      // SageMath lands in QQ and then fails to coerce back into ZZ with a
      // TypeError; the port's ring throws instead.
      throw new TypeError(`${this.value}/${other.value} is not an integer`);
    }
    return new IntegerElement(this.value / other.value);
  }
  neg(): IntegerElement {
    return new IntegerElement(-this.value);
  }
  eq(other: IntegerElement | number | bigint): boolean {
    if (typeof other === 'number') return this.value === BigInt(other);
    if (typeof other === 'bigint') return this.value === other;
    return this.value === other.value;
  }
  isZero(): boolean {
    return this.value === 0n;
  }
  isOne(): boolean {
    return this.value === 1n;
  }
  isUnit(): boolean {
    return this.value === 1n || this.value === -1n;
  }
  is_square(): boolean {
    return isPerfectSquare(this.value);
  }
  sqrt(): IntegerElement {
    if (!isPerfectSquare(this.value)) {
      throw new Error(`unable to take the square root of ${this}`);
    }
    return new IntegerElement(isqrt(this.value));
  }
  toString(): string {
    return this.value.toString();
  }
}

const ZZ: CoefficientRing<IntegerElement> = {
  zero: () => new IntegerElement(0n),
  one: () => new IntegerElement(1n),
  __call__(x: unknown): IntegerElement {
    if (x instanceof IntegerElement) return x;
    if (typeof x === 'bigint') return new IntegerElement(x);
    if (typeof x === 'number') return new IntegerElement(BigInt(x));
    throw new Error(`Cannot convert ${typeof x} to IntegerElement`);
  },
  is_field: () => false,
  characteristic: () => 0n,
};

/** GF(p) for a small prime p. */
class GFElement implements RingElement {
  readonly value: bigint;
  readonly p: bigint;

  constructor(value: bigint, p: bigint) {
    this.value = ((value % p) + p) % p;
    this.p = p;
  }

  add(other: GFElement): GFElement {
    return new GFElement(this.value + other.value, this.p);
  }
  sub(other: GFElement): GFElement {
    return new GFElement(this.value - other.value, this.p);
  }
  mul(other: GFElement): GFElement {
    return new GFElement(this.value * other.value, this.p);
  }
  div(other: GFElement): GFElement {
    return this.mul(other.inv());
  }
  neg(): GFElement {
    return new GFElement(-this.value, this.p);
  }
  inv(): GFElement {
    if (this.value === 0n) throw new Error('division by zero');
    let r = 1n;
    let b = this.value;
    let e = this.p - 2n;
    while (e > 0n) {
      if (e & 1n) r = (r * b) % this.p;
      b = (b * b) % this.p;
      e >>= 1n;
    }
    return new GFElement(r, this.p);
  }
  eq(other: GFElement | number | bigint): boolean {
    if (typeof other === 'number')
      return this.value === ((BigInt(other) % this.p) + this.p) % this.p;
    if (typeof other === 'bigint') return this.value === ((other % this.p) + this.p) % this.p;
    return this.value === other.value;
  }
  isZero(): boolean {
    return this.value === 0n;
  }
  isOne(): boolean {
    return this.value === 1n;
  }
  isUnit(): boolean {
    return this.value !== 0n;
  }
  toString(): string {
    return this.value.toString();
  }
}

function GF(p: bigint): CoefficientRing<GFElement> {
  return {
    zero: () => new GFElement(0n, p),
    one: () => new GFElement(1n, p),
    __call__(x: unknown): GFElement {
      if (x instanceof GFElement) return x;
      if (typeof x === 'bigint') return new GFElement(x, p);
      if (typeof x === 'number') return new GFElement(BigInt(x), p);
      throw new Error(`Cannot convert ${typeof x} to GF(${p}) element`);
    },
    is_field: () => true,
    characteristic: () => p,
  };
}

// Convenience: build a Laurent series from a coefficient list and a shift.
function L<T extends RingElement>(
  R: LaurentSeriesRing<T>,
  coeffs: (number | bigint | string)[],
  n = 0
): LaurentSeriesElement<T> {
  return new LaurentSeriesElement<T>(R, coeffs, n);
}

// ---------------------------------------------------------------------------

describe('LaurentSeriesRing (sage/rings/laurent_series_ring.py)', () => {
  test('repr, gens and basic accessors', () => {
    // sage: LaurentSeriesRing(QQ, 'q')
    // Laurent Series Ring in q over Rational Field
    const R = new LaurentSeriesRing(QQ, 'q');
    expect(R.variable_name()).toBe('q');
    expect(R.ngens()).toBe(1);
    expect(R.gen().toString()).toBe('q');
    expect(R.default_prec()).toBe(20);
    expect(R.is_field()).toBe(true);
    expect(R.is_exact()).toBe(false);
    expect(R.characteristic()).toBe(0n);
    expect(R.uniformizer().toString()).toBe('q');
    // sage: R.<x> = LaurentSeriesRing(GF(17)); R.characteristic() -> 17
    expect(new LaurentSeriesRing(GF(17n), 'x').characteristic()).toBe(17n);
    // sage: LaurentSeriesRing(ZZ,'t').is_field() -> False
    expect(new LaurentSeriesRing(ZZ, 't').is_field()).toBe(false);
    // sage: R.<x> = LaurentSeriesRing(ZZ); R.residue_field()
    // TypeError: the base ring is not a field
    expect(() => new LaurentSeriesRing(ZZ, 'x').residue_field()).toThrow(
      'the base ring is not a field'
    );
    expect(() => new LaurentSeriesRing(ZZ, 'x').uniformizer()).toThrow(
      'the base ring is not a field'
    );
  });

  test('construction from a list and a shift', () => {
    // sage: R.<q> = LaurentSeriesRing(ZZ)
    // sage: R([1,2,3]) -> 1 + 2*q + 3*q^2
    // sage: R([1,2,3], -5) -> q^-5 + 2*q^-4 + 3*q^-3
    const R = new LaurentSeriesRing(ZZ, 'q');
    expect(R.__call__([1, 2, 3]).toString()).toBe('1 + 2*q + 3*q^2');
    expect(R.__call__([1, 2, 3], -5).toString()).toBe('q^-5 + 2*q^-4 + 3*q^-3');
    // sage: P.<x> = LaurentSeriesRing(QQ); P({-3: 1}) -> x^-3
    const P = new LaurentSeriesRing(QQ, 'x');
    expect(P.__call__(new Map([[-3, 1]])).toString()).toBe('x^-3');
  });

  test('1/(1-t+O(t^10)) over GF(7)', () => {
    // sage: R.<t> = LaurentSeriesRing(GF(7), 't')
    // sage: f = 1/(1-t+O(t^10)); f
    // 1 + t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + t^8 + t^9 + O(t^10)
    // sage: f[2] -> 1
    const R = new LaurentSeriesRing(GF(7n), 't');
    const f = R.one().div(L(R, [1, -1]).add_bigoh(10));
    expect(f.toString()).toBe('1 + t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + t^8 + t^9 + O(t^10)');
    expect(f.__getitem__(2).toString()).toBe('1');
  });
});

describe('LaurentSeries predicates (laurent_series_ring_element.pyx)', () => {
  test('is_unit / inverse', () => {
    // sage: R.<t> = LaurentSeriesRing(QQ)
    // sage: (2 + t).is_unit() -> True
    // sage: f = 2 + t^2 + O(t^10); f.is_unit() -> True
    // sage: 1/f -> 1/2 - 1/4*t^2 + 1/8*t^4 - 1/16*t^6 + 1/32*t^8 + O(t^10)
    // sage: R(0).is_unit() -> False
    const R = new LaurentSeriesRing(QQ, 't');
    expect(L(R, [2, 1]).is_unit()).toBe(true);
    const f = L(R, [2, 0, 1]).add_bigoh(10);
    expect(f.is_unit()).toBe(true);
    expect(R.one().div(f).toString()).toBe(
      '1/2 - 1/4*t^2 + 1/8*t^4 - 1/16*t^6 + 1/32*t^8 + O(t^10)'
    );
    expect(R.__call__(0).is_unit()).toBe(false);

    // sage: R.<s> = LaurentSeriesRing(ZZ)
    // sage: f = 2 + s^2 + O(s^10); f.is_unit() -> False
    // sage: 1/f -> ValueError: constant term 2 is not a unit
    const S = new LaurentSeriesRing(ZZ, 's');
    const g = L(S, [2, 0, 1]).add_bigoh(10);
    expect(g.is_unit()).toBe(false);
    expect(() => S.one().div(g)).toThrow('constant term 2 is not a unit');
  });

  test('is_zero and bool', () => {
    // sage: x = Frac(QQ[['x']]).0
    // sage: f = 1/x + x + x^2 + 3*x^4 + O(x^7); f.is_zero() -> 0
    // sage: z = 0*f; z.is_zero() -> 1
    const R = new LaurentSeriesRing(QQ, 'x');
    const f = L(R, [1, 0, 1, 1, 0, 3], -1).add_bigoh(7);
    expect(f.is_zero()).toBe(false);
    expect(R.__call__(0).mul(f).is_zero()).toBe(true);

    // sage: bool(t), bool(1/t), bool(2+t), bool(O(t^3)), bool(O(t^-3)) ...
    const t = R.gen();
    expect(t.bool()).toBe(true);
    expect(t.inverse().bool()).toBe(true);
    expect(R.zero().add_bigoh(3).bool()).toBe(false);
    expect(R.zero().add_bigoh(-3).bool()).toBe(false);
    expect(R.zero().bool()).toBe(false);
  });

  test('is_monomial', () => {
    // sage: k.<z> = LaurentSeriesRing(QQ, 'z')
    // sage: (30*z).is_monomial() -> False
    // sage: k(1).is_monomial() -> True
    // sage: (z+1).is_monomial() -> False
    // sage: (z^-2909).is_monomial() -> True
    // sage: (3*z^-2909).is_monomial() -> False
    const k = new LaurentSeriesRing(QQ, 'z');
    const z = k.gen();
    expect(z.scalar_mul(QQ.__call__(30)).is_monomial()).toBe(false);
    expect(k.__call__(1).is_monomial()).toBe(true);
    expect(z.add(k.__call__(1)).is_monomial()).toBe(false);
    expect(z.pow(-2909).is_monomial()).toBe(true);
    expect(z.pow(-2909).scalar_mul(QQ.__call__(3)).is_monomial()).toBe(false);
  });
});

describe('LaurentSeries repr (laurent_series_ring_element.pyx:325)', () => {
  const R = new LaurentSeriesRing(QQ, 't');

  test('2 + 2/3*t^3', () => {
    // sage: (2 + (2/3)*t^3).__repr__() -> '2 + 2/3*t^3'
    expect(L(R, [2, 0, 0, '2/3']).toString()).toBe('2 + 2/3*t^3');
  });

  test('zero and big-oh only', () => {
    expect(R.zero().toString()).toBe('0');
    expect(R.zero().add_bigoh(3).toString()).toBe('O(t^3)');
    // The zero branch of Sage's _repr_ always prints ``O(t^prec)`` (the
    // ``O(1)``/``O(t)`` spellings only apply to the big-oh appended to a
    // nonzero series).
    expect(R.zero().add_bigoh(1).toString()).toBe('O(t^1)');
    expect(R.zero().add_bigoh(0).toString()).toBe('O(t^0)');
  });
});

describe('LaurentSeries coefficients (laurent_series_ring_element.pyx)', () => {
  const R = new LaurentSeriesRing(QQ, 't');

  test('__getitem__', () => {
    // sage: f = -5/t^(10) + t + t^2 - 10/3*t^3
    // sage: f[-10] -> -5 ; f[1] -> 1 ; f[3] -> -10/3 ; f[-9] -> 0
    const f = L(R, [-5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, '-10/3'], -10);
    expect(f.toString()).toBe('-5*t^-10 + t + t^2 - 10/3*t^3');
    expect(f.__getitem__(-10).toString()).toBe('-5');
    expect(f.__getitem__(1).toString()).toBe('1');
    expect(f.__getitem__(3).toString()).toBe('-10/3');
    expect(f.__getitem__(-9).toString()).toBe('0');
  });

  test('list, coefficients, exponents, residue', () => {
    // sage: f = -5/t^(2) + t + t^2 - 10/3*t^3
    // sage: f.list() -> [-5, 0, 0, 1, 1, -10/3]
    // sage: f.coefficients() -> [-5, 1, 1, -10/3]
    // sage: f.exponents() -> [-2, 1, 2, 3]
    const f = L(R, [-5, 0, 0, 1, 1, '-10/3'], -2);
    expect(f.toString()).toBe('-5*t^-2 + t + t^2 - 10/3*t^3');
    expect(f.list().map(String)).toEqual(['-5', '0', '0', '1', '1', '-10/3']);
    expect(f.coefficients().map(String)).toEqual(['-5', '1', '1', '-10/3']);
    expect(f.exponents()).toEqual([-2, 1, 2, 3]);

    // sage: t = LaurentSeriesRing(ZZ,'t').gen()
    // sage: f = 1/t**2 + 2/t + 3 + 4*t; f.residue() -> 2
    // sage: f = t + t**2; f.residue() -> 0
    const Rz = new LaurentSeriesRing(ZZ, 't');
    expect(L(Rz, [1, 2, 3, 4], -2).residue().toString()).toBe('2');
    expect(L(Rz, [1, 1], 1).residue().toString()).toBe('0');
  });

  test('lift_to_precision', () => {
    // sage: A.<t> = LaurentSeriesRing(GF(5))
    // sage: x = t^(-1) + t^2 + O(t^5)
    // sage: x.lift_to_precision(10) -> t^-1 + t^2 + O(t^10)
    // sage: x.lift_to_precision() -> t^-1 + t^2
    const A = new LaurentSeriesRing(GF(5n), 't');
    const x = L(A, [1, 0, 0, 1], -1).add_bigoh(5);
    expect(x.toString()).toBe('t^-1 + t^2 + O(t^5)');
    expect(x.lift_to_precision(10).toString()).toBe('t^-1 + t^2 + O(t^10)');
    expect(x.lift_to_precision().toString()).toBe('t^-1 + t^2');
  });
});

describe('LaurentSeries arithmetic (laurent_series_ring_element.pyx)', () => {
  const R = new LaurentSeriesRing(QQ, 't');
  const t = R.gen();

  test('_add_', () => {
    // sage: t + t -> 2*t
    expect(t.add(t).toString()).toBe('2*t');
    // sage: f = 1/t + t^2 + t^3 - 17/3 * t^4 + O(t^5)
    // sage: g = 1/(1-t + O(t^7)); g
    // 1 + t + t^2 + t^3 + t^4 + t^5 + t^6 + O(t^7)
    // sage: f + g -> t^-1 + 1 + t + 2*t^2 + 2*t^3 - 14/3*t^4 + O(t^5)
    const f = L(R, [1, 0, 0, 1, 1, '-17/3'], -1).add_bigoh(5);
    const g = R.one().div(L(R, [1, -1]).add_bigoh(7));
    expect(g.toString()).toBe('1 + t + t^2 + t^3 + t^4 + t^5 + t^6 + O(t^7)');
    expect(f.add(g).toString()).toBe('t^-1 + 1 + t + 2*t^2 + 2*t^3 - 14/3*t^4 + O(t^5)');
    // sage: f + 0 -> t^-1 + t^2 + t^3 - 17/3*t^4 + O(t^5)
    expect(f.add(R.__call__(0)).toString()).toBe('t^-1 + t^2 + t^3 - 17/3*t^4 + O(t^5)');
    expect(R.__call__(0).add(f).toString()).toBe('t^-1 + t^2 + t^3 - 17/3*t^4 + O(t^5)');
    // sage: R(0) + R(0) -> 0
    expect(R.__call__(0).add(R.__call__(0)).toString()).toBe('0');
    // sage: (t^3 + O(t^10)) + (t^-3 + O(t^9)) -> t^-3 + t^3 + O(t^9)
    expect(t.pow(3).add_bigoh(10).add(t.pow(-3).add_bigoh(9)).toString()).toBe(
      't^-3 + t^3 + O(t^9)'
    );
  });

  test('_sub_', () => {
    // sage: t - t -> 0
    expect(t.sub(t).toString()).toBe('0');
    // sage: t^5 + 2 * t^-5 -> 2*t^-5 + t^5
    expect(
      t
        .pow(5)
        .add(t.pow(-5).scalar_mul(QQ.__call__(2)))
        .toString()
    ).toBe('2*t^-5 + t^5');
  });

  test('__neg__', () => {
    // sage: -(1+t^5) -> -1 - t^5
    expect(L(R, [1, 0, 0, 0, 0, 1]).neg().toString()).toBe('-1 - t^5');
    // sage: -(1/(1+t+O(t^5))) -> -1 + t - t^2 + t^3 - t^4 + O(t^5)
    expect(
      R.one()
        .div(L(R, [1, 1]).add_bigoh(5))
        .neg()
        .toString()
    ).toBe('-1 + t - t^2 + t^3 - t^4 + O(t^5)');
  });

  test('_mul_', () => {
    // sage: x = Frac(QQ[['x']]).0
    // sage: f = 1/x^3 + x + x^2 + 3*x^4 + O(x^7)
    // sage: g = 1 - x + x^2 - x^4 + O(x^8)
    // sage: f*g -> x^-3 - x^-2 + x^-1 + 4*x^4 + O(x^5)
    const Rx = new LaurentSeriesRing(QQ, 'x');
    const f = L(Rx, [1, 0, 0, 0, 1, 1, 0, 3], -3).add_bigoh(7);
    const g = L(Rx, [1, -1, 1, 0, -1]).add_bigoh(8);
    expect(f.mul(g).toString()).toBe('x^-3 - x^-2 + x^-1 + 4*x^4 + O(x^5)');
  });

  test('_div_', () => {
    // sage: x = Frac(QQ[['x']]).0
    // sage: f = x + x^2 + 3*x^4 + O(x^7)
    // sage: g = 1/x^7 - x + x^2 - x^4 + O(x^8)
    // sage: f/x -> 1 + x + 3*x^3 + O(x^6)
    // sage: f/g -> x^8 + x^9 + 3*x^11 + O(x^14)
    const Rx = new LaurentSeriesRing(QQ, 'x');
    const x = Rx.gen();
    const f = L(Rx, [1, 1, 0, 3], 1).add_bigoh(7);
    const g = L(Rx, [1, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1], -7).add_bigoh(8);
    expect(f.div(x).toString()).toBe('1 + x + 3*x^3 + O(x^6)');
    expect(f.div(g).toString()).toBe('x^8 + x^9 + 3*x^11 + O(x^14)');
    expect(() => f.div(Rx.zero())).toThrow();
  });

  test('__pow__ with integer and rational exponents', () => {
    // sage: x = Frac(QQ[['x']]).0
    // sage: f = x + x^2 + 3*x^4 + O(x^7)
    // sage: g = 1/x^10 - x + x^2 - x^4 + O(x^8)
    // sage: f^7 -> x^7 + 7*x^8 + 21*x^9 + 56*x^10 + 161*x^11 + 336*x^12 + O(x^13)
    // sage: g^7 -> x^-70 - 7*x^-59 + 7*x^-58 - 7*x^-56 + O(x^-52)
    // sage: g^(1/2) -> x^-5 - 1/2*x^6 + 1/2*x^7 - 1/2*x^9 + O(x^13)
    // sage: g^(1/5) -> x^-2 - 1/5*x^9 + 1/5*x^10 - 1/5*x^12 + O(x^16)
    // sage: g^(2/5) -> x^-4 - 2/5*x^7 + 2/5*x^8 - 2/5*x^10 + O(x^14)
    // sage: h = x^2 + 2*x^4 + x^6; h^(1/2) -> x + x^3
    const Rx = new LaurentSeriesRing(QQ, 'x');
    const f = L(Rx, [1, 1, 0, 3], 1).add_bigoh(7);
    const g = L(Rx, [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1], -10).add_bigoh(8);
    expect(f.pow(7).toString()).toBe(
      'x^7 + 7*x^8 + 21*x^9 + 56*x^10 + 161*x^11 + 336*x^12 + O(x^13)'
    );
    expect(g.pow(7).toString()).toBe('x^-70 - 7*x^-59 + 7*x^-58 - 7*x^-56 + O(x^-52)');
    expect(g.pow(new Rational(1n, 2n)).toString()).toBe(
      'x^-5 - 1/2*x^6 + 1/2*x^7 - 1/2*x^9 + O(x^13)'
    );
    expect(g.pow(new Rational(1n, 5n)).toString()).toBe(
      'x^-2 - 1/5*x^9 + 1/5*x^10 - 1/5*x^12 + O(x^16)'
    );
    expect(g.pow(new Rational(2n, 5n)).toString()).toBe(
      'x^-4 - 2/5*x^7 + 2/5*x^8 - 2/5*x^10 + O(x^14)'
    );
    const h = L(Rx, [1, 0, 2, 0, 1], 2);
    expect(h.pow(new Rational(1n, 2n)).toString()).toBe('x + x^3');
  });
});

describe('LaurentSeries precision and shifting', () => {
  const R = new LaurentSeriesRing(QQ, 't');
  const t = R.gen();

  test('add_bigoh', () => {
    // sage: f = t^2 + t^3 + O(t^10); f.add_bigoh(5) -> t^2 + t^3 + O(t^5)
    const f = L(R, [1, 1], 2).add_bigoh(10);
    expect(f.add_bigoh(5).toString()).toBe('t^2 + t^3 + O(t^5)');
    // sage: (t^(-2)).add_bigoh(-1) -> t^-2 + O(t^-1)
    // sage: (t^(-2)).add_bigoh(-2) -> O(t^-2)
    // sage: (t^(-2)).add_bigoh(-3) -> O(t^-3)
    expect(t.pow(-2).add_bigoh(-1).toString()).toBe('t^-2 + O(t^-1)');
    expect(t.pow(-2).add_bigoh(-2).toString()).toBe('O(t^-2)');
    expect(t.pow(-2).add_bigoh(-3).toString()).toBe('O(t^-3)');
  });

  test('O', () => {
    // sage: f = t^-5 + t^-4 + t^3 + O(t^10)
    // sage: f.O(-4) -> t^-5 + O(t^-4)
    // sage: f.O(15) -> t^-5 + t^-4 + t^3 + O(t^10)
    const f = L(R, [1, 1, 0, 0, 0, 0, 0, 0, 1], -5).add_bigoh(10);
    expect(f.O(-4).toString()).toBe('t^-5 + O(t^-4)');
    expect(f.O(15).toString()).toBe('t^-5 + t^-4 + t^3 + O(t^10)');
  });

  test('degree', () => {
    // sage: g = x^2 - x^4 + O(x^8); g.degree() -> 4
    // sage: g = -10/x^5 + x^2 - x^4 + O(x^8); g.degree() -> 4
    // sage: (x^-2 + O(x^0)).degree() -> -2
    const Rx = new LaurentSeriesRing(QQ, 'x');
    expect(L(Rx, [1, 0, -1], 2).add_bigoh(8).degree()).toBe(4);
    expect(L(Rx, [-10, 0, 0, 0, 0, 0, 0, 1, 0, -1], -5).add_bigoh(8).degree()).toBe(4);
    expect(Rx.gen().pow(-2).add_bigoh(0).degree()).toBe(-2);
  });

  test('valuation, prec, precision_relative', () => {
    // sage: f = 1/x + x^2 + 3*x^4 + O(x^7); f.valuation() -> -1
    // sage: g = 1 - x + x^2 - x^4 + O(x^8); g.valuation() -> 0
    // sage: h = f - f; h -> O(x^7); h.valuation() -> +Infinity
    // sage: R(0).valuation() -> +Infinity
    const Rx = new LaurentSeriesRing(QQ, 'x');
    const f = L(Rx, [1, 0, 0, 1, 0, 3], -1).add_bigoh(7);
    const g = L(Rx, [1, -1, 1, 0, -1]).add_bigoh(8);
    expect(f.valuation()).toBe(-1);
    expect(g.valuation()).toBe(0);
    const h = f.sub(f);
    expect(h.toString()).toBe('O(x^7)');
    expect(h.valuation()).toBe(Number.POSITIVE_INFINITY);
    expect(Rx.__call__(0).valuation()).toBe(Number.POSITIVE_INFINITY);

    // sage: f = x^2 + 3*x^4 + O(x^7); f.prec() -> 7
    // sage: g = 1/x^10 - x + x^2 - x^4 + O(x^8); g.prec() -> 8
    expect(L(Rx, [1, 0, 3], 2).add_bigoh(7).prec()).toBe(7);
    const g2 = L(Rx, [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1], -10).add_bigoh(8);
    expect(g2.prec()).toBe(8);

    // sage: (t^2 + O(t^3)).precision_relative() -> 1
    // sage: (1 - t^2 + O(t^100)).precision_relative() -> 100
    // sage: O(t^4).precision_relative() -> 0
    expect(L(R, [1], 2).add_bigoh(3).precision_relative()).toBe(1);
    expect(L(R, [1, 0, -1]).add_bigoh(100).precision_relative()).toBe(100);
    expect(R.zero().add_bigoh(4).precision_relative()).toBe(0);
    expect(L(R, [1], 2).add_bigoh(3).precision_absolute()).toBe(3);
  });

  test('common_prec and common_valuation', () => {
    // sage: f = t^(-1) + t + t^2 + O(t^3); g = t + t^3 + t^4 + O(t^4)
    // sage: f.common_prec(g) -> 3 ; f.common_valuation(g) -> -1
    const f = L(R, [1, 0, 1, 1], -1).add_bigoh(3);
    const g = L(R, [1, 0, 1, 1], 1).add_bigoh(4);
    expect(f.common_prec(g)).toBe(3);
    expect(g.common_prec(f)).toBe(3);
    expect(f.common_valuation(g)).toBe(-1);
    expect(g.common_valuation(f)).toBe(-1);

    // sage: f = t + t^2 + O(t^3); g = t^(-3) + t^2
    // sage: f.common_prec(g) -> 3 ; f.common_valuation(g) -> -3
    const f2 = L(R, [1, 1], 1).add_bigoh(3);
    const g2 = L(R, [1, 0, 0, 0, 0, 1], -3);
    expect(f2.common_prec(g2)).toBe(3);
    expect(f2.common_valuation(g2)).toBe(-3);

    // sage: f = t + t^2; g = t^2
    // sage: f.common_prec(g) -> +Infinity ; f.common_valuation(g) -> 1
    const f3 = L(R, [1, 1], 1);
    const g3 = L(R, [1], 2);
    expect(f3.common_prec(g3)).toBe(Number.POSITIVE_INFINITY);
    expect(f3.common_valuation(g3)).toBe(1);

    // sage: f = t^(-3) + O(t^(-2)); g = t^(-5) + O(t^(-1))
    // sage: f.common_prec(g) -> -2 ; f.common_valuation(g) -> -5
    const f4 = L(R, [1], -3).add_bigoh(-2);
    const g4 = L(R, [1], -5).add_bigoh(-1);
    expect(f4.common_prec(g4)).toBe(-2);
    expect(f4.common_valuation(g4)).toBe(-5);

    // sage: f = O(t^2); g = O(t^5)
    // sage: f.common_prec(g) -> 2 ; f.common_valuation(g) -> +Infinity
    expect(R.zero().add_bigoh(2).common_prec(R.zero().add_bigoh(5))).toBe(2);
    expect(R.zero().add_bigoh(2).common_valuation(R.zero().add_bigoh(5))).toBe(
      Number.POSITIVE_INFINITY
    );
  });

  test('equality compares up to the common precision', () => {
    // sage: f = x^(-1) + 1 + x + O(x^2); g = x^(-1) + 1 + O(x)
    // sage: f == g -> True
    // sage: g = x^(-1) + 2 + O(x); f == g -> False
    const Rx = new LaurentSeriesRing(QQ, 'x');
    const f = L(Rx, [1, 1, 1], -1).add_bigoh(2);
    const g = L(Rx, [1, 1], -1).add_bigoh(1);
    expect(f.eq(g)).toBe(true);
    const g2 = L(Rx, [1, 2], -1).add_bigoh(1);
    expect(f.eq(g2)).toBe(false);
    // sage: f = x^(-2) + 1 + x + O(x^2); g = x^(-1) + 2 + O(x); f == g -> False
    const f2 = L(Rx, [1, 0, 1, 1], -2).add_bigoh(2);
    expect(f2.eq(g2)).toBe(false);
  });

  test('shift, << and >>', () => {
    // sage: f = (t+t^-1)^4; f -> t^-4 + 4*t^-2 + 6 + 4*t^2 + t^4
    // sage: f.shift(10) -> t^6 + 4*t^8 + 6*t^10 + 4*t^12 + t^14
    // sage: f >> 10 -> t^-14 + 4*t^-12 + 6*t^-10 + 4*t^-8 + t^-6
    // sage: t << 4 -> t^5
    // sage: t + O(t^3) >> 4 -> t^-3 + O(t^-1)
    const f = t.add(t.pow(-1)).pow(4);
    expect(f.toString()).toBe('t^-4 + 4*t^-2 + 6 + 4*t^2 + t^4');
    expect(f.shift(10).toString()).toBe('t^6 + 4*t^8 + 6*t^10 + 4*t^12 + t^14');
    expect(f.rshift(10).toString()).toBe('t^-14 + 4*t^-12 + 6*t^-10 + 4*t^-8 + t^-6');
    expect(t.lshift(4).toString()).toBe('t^5');
    expect(t.add_bigoh(3).rshift(4).toString()).toBe('t^-3 + O(t^-1)');
  });

  test('truncate, truncate_laurentseries, truncate_neg', () => {
    // sage: A.<x> = LaurentSeriesRing(ZZ); f = 1/(1-x)
    // sage: f -> 1 + x + ... + x^19 + O(x^20)
    // sage: f.truncate(10) -> 1 + x + x^2 + ... + x^9
    // sage: f.truncate_laurentseries(10) -> 1 + ... + x^9 + O(x^10)
    const A = new LaurentSeriesRing(ZZ, 'x');
    const f = A.one().div(L(A, [1, -1]));
    expect(f.toString()).toBe(
      '1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + x^10 + x^11 + x^12 + x^13 + x^14 + x^15 + x^16 + x^17 + x^18 + x^19 + O(x^20)'
    );
    expect(f.truncate(10).toString()).toBe('1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9');
    expect(f.truncate_laurentseries(10).toString()).toBe(
      '1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + O(x^10)'
    );

    // sage: A.<t> = LaurentSeriesRing(ZZ); f = 1/(1-t)
    // sage: f.truncate_neg(15) -> t^15 + t^16 + t^17 + t^18 + t^19 + O(t^20)
    const At = new LaurentSeriesRing(ZZ, 't');
    const ft = At.one().div(L(At, [1, -1]));
    expect(ft.truncate_neg(15).toString()).toBe('t^15 + t^16 + t^17 + t^18 + t^19 + O(t^20)');
    // sage: S.<t> = LaurentSeriesRing(QQ)
    // sage: (t+t^2).truncate_neg(-1) -> t + t^2
    // sage: (t+t^2).truncate_neg(-2) -> t + t^2
    const s = L(R, [1, 1], 1);
    expect(s.truncate_neg(-1).toString()).toBe('t + t^2');
    expect(s.truncate_neg(-2).toString()).toBe('t + t^2');
  });
});

describe('LaurentSeries verschiebung (laurent_series_ring_element.pyx:366)', () => {
  const R = new LaurentSeriesRing(QQ, 'x');
  const x = R.gen();

  test('V on exact and inexact series', () => {
    // sage: f = -1/x + 1 + 2*x^2 + 5*x^5
    // sage: f.V(2) -> -x^-2 + 1 + 2*x^4 + 5*x^10
    // sage: f.V(-1) -> 5*x^-5 + 2*x^-2 + 1 - x
    // sage: h = f.add_bigoh(7); h.V(2) -> -x^-2 + 1 + 2*x^4 + 5*x^10 + O(x^14)
    // sage: h.V(-2) -> ValueError: For finite precision only positive arguments allowed
    const f = L(R, [-1, 1, 0, 2, 0, 0, 5], -1);
    expect(f.toString()).toBe('-x^-1 + 1 + 2*x^2 + 5*x^5');
    expect(f.V(2).toString()).toBe('-x^-2 + 1 + 2*x^4 + 5*x^10');
    expect(f.V(-1).toString()).toBe('5*x^-5 + 2*x^-2 + 1 - x');
    const h = f.add_bigoh(7);
    expect(h.V(2).toString()).toBe('-x^-2 + 1 + 2*x^4 + 5*x^10 + O(x^14)');
    expect(() => h.V(-2)).toThrow('For finite precision only positive arguments allowed');
    expect(() => f.V(0)).toThrow('n must be nonzero');

    // sage: f = x; f.V(3) -> x^3 ; f.V(-3) -> x^-3
    expect(x.V(3).toString()).toBe('x^3');
    expect(x.V(-3).toString()).toBe('x^-3');
    // sage: g = 2*x^(-1) + 3 + 5*x; g.V(-1) -> 5*x^-1 + 3 + 2*x
    const g = L(R, [2, 3, 5], -1);
    expect(g.V(-1).toString()).toBe('5*x^-1 + 3 + 2*x');
  });
});

describe('LaurentSeries calculus (laurent_series_ring_element.pyx)', () => {
  test('derivative', () => {
    // sage: R.<x> = LaurentSeriesRing(QQ)
    // sage: g = 1/x^10 - x + x^2 - x^4 + O(x^8)
    // sage: g.derivative() -> -10*x^-11 - 1 + 2*x - 4*x^3 + O(x^7)
    const R = new LaurentSeriesRing(QQ, 'x');
    const g = L(R, [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1], -10).add_bigoh(8);
    expect(g.derivative().toString()).toBe('-10*x^-11 - 1 + 2*x - 4*x^3 + O(x^7)');
    // sage: f = x^2 + 3*x^4 + O(x^7); f._derivative() -> 2*x + 12*x^3 + O(x^6)
    const f = L(R, [1, 0, 3], 2).add_bigoh(7);
    expect(f.derivative().toString()).toBe('2*x + 12*x^3 + O(x^6)');
  });

  test('integral', () => {
    // sage: t = LaurentSeriesRing(ZZ, 't').0
    // sage: f = 2*t^-3 + 3*t^2 + O(t^4); f.integral() -> -t^-2 + t^3 + O(t^5)
    const Rz = new LaurentSeriesRing(ZZ, 't');
    const f = L(Rz, [2, 0, 0, 0, 0, 3], -3).add_bigoh(4);
    expect(f.toString()).toBe('2*t^-3 + 3*t^2 + O(t^4)');
    expect(f.integral().toString()).toBe('-t^-2 + t^3 + O(t^5)');

    // sage: f = t^3; f.integral()
    // ArithmeticError: Coefficients of integral cannot be coerced into the base ring
    expect(() => L(Rz, [1], 3).integral()).toThrow(
      'Coefficients of integral cannot be coerced into the base ring'
    );

    // sage: t = Frac(QQ[['t']]).0
    // sage: f = -1/t^3 - 31/t + O(t^3); f.integral()
    // ArithmeticError: The integral of is not a Laurent series, since t^-1 has
    // nonzero coefficient.
    const Rq = new LaurentSeriesRing(QQ, 't');
    const f2 = L(Rq, [-1, 0, -31], -3).add_bigoh(3);
    expect(() => f2.integral()).toThrow(
      'The integral of is not a Laurent series, since t^-1 has nonzero coefficient.'
    );

    // sage: A.<t> = QQ[[]]; f = -2*t^(-4) + O(t^8)
    // sage: f.integral() -> 2/3*t^-3 + O(t^9)
    // sage: f.integral().derivative() == f -> True
    const f3 = L(Rq, [-2], -4).add_bigoh(8);
    expect(f3.integral().toString()).toBe('2/3*t^-3 + O(t^9)');
    expect(f3.integral().derivative().eq(f3)).toBe(true);
  });

  test('nth_root', () => {
    // sage: R.<x> = LaurentSeriesRing(QQ)
    // sage: (x^-2 + 1 + x).nth_root(2)
    // x^-1 + 1/2*x + 1/2*x^2 - ... - 19437/65536*x^18 + O(x^19)
    // sage: (x^-2 + 1 + x).nth_root(2)**2 -> x^-2 + 1 + x + O(x^18)
    const R = new LaurentSeriesRing(QQ, 'x');
    const f = L(R, [1, 0, 1, 1], -2);
    const r = f.nth_root(2);
    expect(r.toString().startsWith('x^-1 + 1/2*x + 1/2*x^2 -')).toBe(true);
    expect(r.toString().endsWith('- 19437/65536*x^18 + O(x^19)')).toBe(true);
    expect(r.mul(r).toString()).toBe('x^-2 + 1 + x + O(x^18)');
    expect(() => f.nth_root(0)).toThrow('n must be positive');
    expect(() => R.gen().nth_root(2)).toThrow('valuation must be divisible by n');
  });

  test('reverse', () => {
    // sage: R.<x> = Frac(QQ[['x']])
    // sage: f = 2*x + 3*x^2 - x^4 + O(x^5); g = f.reverse()
    // sage: g -> 1/2*x - 3/8*x^2 + 9/16*x^3 - 131/128*x^4 + O(x^5)
    // sage: f(g) -> x + O(x^5) ; g(f) -> x + O(x^5)
    const R = new LaurentSeriesRing(QQ, 'x');
    const f = L(R, [2, 3, 0, -1], 1).add_bigoh(5);
    const g = f.reverse();
    expect(g.toString()).toBe('1/2*x - 3/8*x^2 + 9/16*x^3 - 131/128*x^4 + O(x^5)');
    expect(f.__call__(g).toString()).toBe('x + O(x^5)');
    expect(g.__call__(f).toString()).toBe('x + O(x^5)');

    // sage: A.<t> = LaurentSeriesRing(ZZ)
    // sage: a = t - t^2 - 2*t^4 + t^5 + O(t^6); b = a.reverse()
    // sage: b -> t + t^2 + 2*t^3 + 7*t^4 + 25*t^5 + O(t^6)
    const A = new LaurentSeriesRing(ZZ, 't');
    const a = L(A, [1, -1, 0, -2, 1], 1).add_bigoh(6);
    expect(a.reverse().toString()).toBe('t + t^2 + 2*t^3 + 7*t^4 + 25*t^5 + O(t^6)');

    // sage: R.<x> = LaurentSeriesRing(QQ, default_prec=20)
    // sage: (x - x^2).reverse()  # Catalan numbers
    expect(L(R, [1, -1], 1).reverse().toString()).toBe(
      'x + x^2 + 2*x^3 + 5*x^4 + 14*x^5 + 42*x^6 + 132*x^7 + 429*x^8 + 1430*x^9 + 4862*x^10 + 16796*x^11 + 58786*x^12 + 208012*x^13 + 742900*x^14 + 2674440*x^15 + 9694845*x^16 + 35357670*x^17 + 129644790*x^18 + 477638700*x^19 + O(x^20)'
    );
    // sage: (x - x^2).reverse(precision=3) -> x + x^2 + O(x^3)
    expect(L(R, [1, -1], 1).reverse(3).toString()).toBe('x + x^2 + O(x^3)');

    // sage: f = 1 + 2*x + 3*x^2 - x^4 + O(x^5); f.reverse()
    // ValueError: Series must have valuation one for reversion.
    expect(() => L(R, [1, 2, 3, 0, -1]).add_bigoh(5).reverse()).toThrow(
      'Series must have valuation one for reversion.'
    );
  });

  test('is_square', () => {
    // sage: R.<x> = LaurentSeriesRing(QQ)
    // sage: (x^2).is_square() -> True
    // sage: (x^3).is_square() -> False
    // sage: (4/x^2 + 4/x + 1).is_square(root=True) -> (True, 2*x^-1 + 1)
    const R = new LaurentSeriesRing(QQ, 'x');
    const x = R.gen();
    expect(x.pow(2).is_square()).toBe(true);
    expect(x.pow(3).is_square()).toBe(false);
    const f = L(R, [4, 4, 1], -2);
    const [ok, root] = f.is_square(true);
    expect(ok).toBe(true);
    expect(root!.toString()).toBe('2*x^-1 + 1');

    // sage: R.<t> = LaurentSeriesRing(ZZ)
    // sage: (t^-4).is_square() -> True
    // sage: (2*t^-4).is_square() -> False
    const Rz = new LaurentSeriesRing(ZZ, 't');
    expect(Rz.gen().pow(-4).is_square()).toBe(true);
    expect(Rz.gen().pow(-4).scalar_mul(ZZ.__call__(2)).is_square()).toBe(false);
  });
});

describe('LaurentSeries conversions (laurent_series_ring_element.pyx)', () => {
  test('power_series', () => {
    // sage: R.<t> = LaurentSeriesRing(ZZ)
    // sage: f = 1/(1-t+O(t^10)); g = f.power_series(); g
    // 1 + t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + t^8 + t^9 + O(t^10)
    const R = new LaurentSeriesRing(ZZ, 't');
    const f = R.one().div(L(R, [1, -1]).add_bigoh(10));
    expect(f.power_series().toString()).toBe(
      '1 + t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + t^8 + t^9 + O(t^10)'
    );
    // sage: f = 3/t^2 + t^2 + t^3 + O(t^10); f.power_series()
    // TypeError: self is not a power series
    const f2 = L(R, [3, 0, 0, 0, 1, 1], -2).add_bigoh(10);
    expect(() => f2.power_series()).toThrow('self is not a power series');

    // sage: S.<x> = PowerSeriesRing(QQ); L = Frac(S); s = L(O(x^2))
    // sage: (s*x^(-1)).power_series() -> O(x^1)
    // sage: (s*x^(-2)).power_series() -> O(x^0)
    // sage: (s*x^(-3)).power_series() -> TypeError: self is not a power series
    const Lq = new LaurentSeriesRing(QQ, 'x');
    const s = Lq.zero().add_bigoh(2);
    expect(s.mul(Lq.gen().pow(-1)).power_series().toString()).toBe('O(x^1)');
    expect(s.mul(Lq.gen().pow(-2)).power_series().toString()).toBe('O(x^0)');
    expect(() => s.mul(Lq.gen().pow(-3)).power_series()).toThrow('self is not a power series');
  });

  test('inverse', () => {
    // sage: R.<t> = LaurentSeriesRing(ZZ)
    // sage: t.inverse() -> t^-1
    // sage: (1-t).inverse() -> 1 + t + t^2 + ...
    const R = new LaurentSeriesRing(ZZ, 't');
    expect(R.gen().inverse().toString()).toBe('t^-1');
    expect(L(R, [1, -1]).inverse().toString()).toBe(
      '1 + t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + t^8 + t^9 + t^10 + t^11 + t^12 + t^13 + t^14 + t^15 + t^16 + t^17 + t^18 + t^19 + O(t^20)'
    );
  });

  test('__call__ substitution', () => {
    // sage: R.<t> = LaurentSeriesRing(...)
    // sage: f = t^-2 + t^2 + O(t^8)
    // sage: f(t^3) -> t^-6 + t^6 + O(t^24)     [pyx:2020, with x=y=1]
    // sage: f(t + O(t^5)) -> t^-2 + O(t^2)
    // sage: f(t^-2) -> ValueError: Can only substitute elements of positive valuation
    const R = new LaurentSeriesRing(QQ, 't');
    const t = R.gen();
    const f = L(R, [1, 0, 0, 0, 1], -2).add_bigoh(8);
    expect(f.__call__(t.pow(3)).toString()).toBe('t^-6 + t^6 + O(t^24)');
    expect(f.__call__(t.add_bigoh(5)).toString()).toBe('t^-2 + O(t^2)');
    expect(() => f.__call__(t.pow(-2))).toThrow(
      'Can only substitute elements of positive valuation'
    );
  });
});

// ---------------------------------------------------------------------------
// Algebraic identities on pseudo-random series (deterministic seed).
// ---------------------------------------------------------------------------

describe('LaurentSeries algebraic identities', () => {
  // A tiny deterministic LCG so that the tests are reproducible.
  let state = 20260728n;
  const nextInt = (bound: number): number => {
    state = (state * 6364136223846793005n + 1442695040888963407n) % (1n << 64n);
    return Number((state >> 33n) % BigInt(bound));
  };

  const R = new LaurentSeriesRing(QQ, 't', 15);

  const randomSeries = (): LaurentSeriesElement<RationalElement> => {
    const n = nextInt(9) - 4; // valuation shift in [-4, 4]
    const len = 6 + nextInt(4);
    const coeffs: string[] = [];
    for (let i = 0; i < len; i++) {
      const num = nextInt(11) - 5;
      const den = 1 + nextInt(4);
      coeffs.push(`${num}/${den}`);
    }
    if (coeffs.every((c) => c.startsWith('0/'))) {
      coeffs[0] = '1';
    }
    return L(R, coeffs, n).add_bigoh(n + 10);
  };

  test('a * a^-1 == 1', () => {
    for (let i = 0; i < 25; i++) {
      const a = randomSeries();
      const one = a.mul(a.inverse());
      expect(one.eq(R.one())).toBe(true);
    }
  });

  test('valuation is additive', () => {
    for (let i = 0; i < 25; i++) {
      const a = randomSeries();
      const b = randomSeries();
      expect(a.mul(b).valuation()).toBe(a.valuation() + b.valuation());
    }
  });

  test('(a*b)/b == a', () => {
    for (let i = 0; i < 25; i++) {
      const a = randomSeries();
      const b = randomSeries();
      expect(a.mul(b).div(b).eq(a)).toBe(true);
    }
  });

  test('(a+b)-b == a and distributivity', () => {
    for (let i = 0; i < 25; i++) {
      const a = randomSeries();
      const b = randomSeries();
      const c = randomSeries();
      expect(a.add(b).sub(b).eq(a)).toBe(true);
      expect(a.mul(b.add(c)).eq(a.mul(b).add(a.mul(c)))).toBe(true);
    }
  });

  test('derivative of a product obeys the Leibniz rule', () => {
    for (let i = 0; i < 25; i++) {
      const a = randomSeries();
      const b = randomSeries();
      const lhs = a.mul(b).derivative();
      const rhs = a.derivative().mul(b).add(a.mul(b.derivative()));
      expect(lhs.eq(rhs)).toBe(true);
    }
  });

  test('a^2 has a square root equal to +-a', () => {
    for (let i = 0; i < 10; i++) {
      const a = randomSeries();
      const sq = a.mul(a);
      const [ok, root] = sq.is_square(true);
      expect(ok).toBe(true);
      expect(root!.mul(root!).eq(sq)).toBe(true);
    }
  });
});

// ---------------------------------------------------------------------------
// End-to-end: the elliptic formal group produces Laurent series, and the
// Weierstrass relation between them is checked with Laurent arithmetic.
// ---------------------------------------------------------------------------

describe('Laurent series produced by the elliptic formal group', () => {
  test('x(10) and y(10) print exactly as in SageMath', () => {
    // sage: EllipticCurve([0, 0, 1, -1, 0]).formal_group().x(10)
    // t^-2 - t + t^2 - t^4 + 2*t^5 - t^6 - 2*t^7 + 6*t^8 - 6*t^9 + O(t^10)
    // sage: EllipticCurve([0, 0, 1, -1, 0]).formal_group().y(10)
    // -t^-3 + 1 - t + t^3 - 2*t^4 + t^5 + 2*t^6 - 6*t^7 + 6*t^8 + 3*t^9 + O(t^10)
    // biome-ignore lint/suspicious/noExplicitAny: the curve interfaces are structural
    const E = EllipticCurve(QQring as any, [0n, 0n, 1n, -1n, 0n] as any) as any;
    const fg = E.formal_group();
    expect(fg.x(10).toString()).toBe(
      't^-2 - t + t^2 - t^4 + 2*t^5 - t^6 - 2*t^7 + 6*t^8 - 6*t^9 + O(t^10)'
    );
    expect(fg.y(10).toString()).toBe(
      '-t^-3 + 1 - t + t^3 - 2*t^4 + t^5 + 2*t^6 - 6*t^7 + 6*t^8 + 3*t^9 + O(t^10)'
    );
  });

  test('the formal point (x(t), y(t)) satisfies the Weierstrass equation', () => {
    // y^2 + a1*x*y + a3*y == x^3 + a2*x^2 + a4*x + a6, computed with Laurent
    // series arithmetic; and t = -x/y is the local parameter.
    // biome-ignore lint/suspicious/noExplicitAny: the curve interfaces are structural
    const E = EllipticCurve(GFring(7n) as any, [3n, 4n] as any) as any;
    const fg = E.formal_group();
    const x = fg.x(7);
    const y = fg.y(7);
    const Lr = x.parent();
    const a = E.ainvs();
    const c = (v: unknown) => Lr.__call__(v);
    const lhs = y.mul(y).add(c(a[0]).mul(x).mul(y)).add(c(a[2]).mul(y));
    const rhs = x.pow(3).add(c(a[1]).mul(x).mul(x)).add(c(a[3]).mul(x)).add(c(a[4]));
    expect(lhs.sub(rhs).is_zero()).toBe(true);
    expect(x.neg().div(y).toString().startsWith('t + O(')).toBe(true);
  });
});
