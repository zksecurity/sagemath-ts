/**
 * Unit tests for power series ring module
 */
import { describe, expect, test } from 'bun:test';
import { EllipticCurve } from '../schemes/elliptic_curves/constructor.js';
import { GF } from './finite_rings/finite_field_constructor.js';
import {
  type CoefficientRing,
  MPowerSeries,
  MPowerSeriesRing,
  PowerSeriesElement,
  PowerSeriesRing,
  type RingElement,
} from './power_series_ring.js';
import { Rational } from './rational.js';

/**
 * Adapter for Rational to implement RingElement interface.
 */
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

  toString(): string {
    return this.value.toString();
  }
}

/**
 * Rational field as coefficient ring.
 */
const QQ: CoefficientRing<RationalElement> = {
  zero(): RationalElement {
    return new RationalElement(Rational.zero());
  },
  one(): RationalElement {
    return new RationalElement(Rational.one());
  },
  __call__(x: unknown): RationalElement {
    if (x instanceof RationalElement) {
      return x;
    }
    if (x instanceof Rational) {
      return new RationalElement(x);
    }
    if (typeof x === 'bigint') {
      return new RationalElement(new Rational(x, 1n));
    }
    if (typeof x === 'number') {
      return new RationalElement(new Rational(BigInt(x), 1n));
    }
    if (typeof x === 'string') {
      return new RationalElement(Rational.fromString(x));
    }
    throw new Error(`Cannot convert ${typeof x} to RationalElement`);
  },
  is_field(): boolean {
    return true;
  },
  characteristic(): bigint {
    return 0n;
  },
};

describe('PowerSeriesRing', () => {
  test('constructor and basic properties', () => {
    const R = new PowerSeriesRing(QQ, 'x', 10);

    expect(R.variable_name()).toBe('x');
    expect(R.default_prec()).toBe(10);
    expect(R.ngens()).toBe(1);
    expect(R.characteristic()).toBe(0n);
  });

  test('zero and one', () => {
    const R = new PowerSeriesRing(QQ, 'x', 10);

    const zero = R.zero();
    const one = R.one();

    expect(zero.is_zero()).toBe(true);
    expect(one.is_one()).toBe(true);
    expect(one.is_zero()).toBe(false);
    expect(zero.is_one()).toBe(false);
  });

  test('generator', () => {
    const R = new PowerSeriesRing(QQ, 'x', 10);
    const x = R.gen();

    expect(x.valuation()).toBe(1);
    expect(x.degree()).toBe(1);
    expect(x.__getitem__(0).isZero()).toBe(true);
    expect(x.__getitem__(1).eq(1)).toBe(true);
  });

  test('coercion from list', () => {
    const R = new PowerSeriesRing(QQ, 'x', 10);

    // 1 + 2x + 3x^2
    const f = R.__call__([1, 2, 3]);

    expect(f.degree()).toBe(2);
    expect(f.__getitem__(0).eq(1)).toBe(true);
    expect(f.__getitem__(1).eq(2)).toBe(true);
    expect(f.__getitem__(2).eq(3)).toBe(true);
  });

  test('coercion from scalar', () => {
    const R = new PowerSeriesRing(QQ, 'x', 10);

    const f = R.__call__(5);

    expect(f.degree()).toBe(0);
    expect(f.__getitem__(0).eq(5)).toBe(true);
    expect(f.valuation()).toBe(0);
  });
});

describe('PowerSeriesElement arithmetic', () => {
  const R = new PowerSeriesRing(QQ, 'x', 20);
  const x = R.gen();

  test('addition', () => {
    // (1 + x) + (2 + 3x) = 3 + 4x
    const f = R.__call__([1, 1]);
    const g = R.__call__([2, 3]);
    const sum = f.add(g);

    expect(sum.__getitem__(0).eq(3)).toBe(true);
    expect(sum.__getitem__(1).eq(4)).toBe(true);
  });

  test('subtraction', () => {
    // (1 + x) - (2 + 3x) = -1 - 2x
    const f = R.__call__([1, 1]);
    const g = R.__call__([2, 3]);
    const diff = f.sub(g);

    expect(diff.__getitem__(0).eq(-1)).toBe(true);
    expect(diff.__getitem__(1).eq(-2)).toBe(true);
  });

  test('multiplication', () => {
    // (1 + x) * (1 + x) = 1 + 2x + x^2
    const f = R.__call__([1, 1]);
    const prod = f.mul(f);

    expect(prod.__getitem__(0).eq(1)).toBe(true);
    expect(prod.__getitem__(1).eq(2)).toBe(true);
    expect(prod.__getitem__(2).eq(1)).toBe(true);
  });

  test('negation', () => {
    const f = R.__call__([1, 2, 3]);
    const neg = f.neg();

    expect(neg.__getitem__(0).eq(-1)).toBe(true);
    expect(neg.__getitem__(1).eq(-2)).toBe(true);
    expect(neg.__getitem__(2).eq(-3)).toBe(true);
  });

  test('power', () => {
    // (1 + x)^3 = 1 + 3x + 3x^2 + x^3
    const f = R.__call__([1, 1]);
    const cube = f.pow(3);

    expect(cube.__getitem__(0).eq(1)).toBe(true);
    expect(cube.__getitem__(1).eq(3)).toBe(true);
    expect(cube.__getitem__(2).eq(3)).toBe(true);
    expect(cube.__getitem__(3).eq(1)).toBe(true);
  });

  test('inverse', () => {
    // 1/(1-x) = 1 + x + x^2 + x^3 + ...
    const f = R.__call__([1, -1]).add_bigoh(10);
    const inv = f.inv();

    // Check first few coefficients are all 1
    for (let i = 0; i < 5; i++) {
      expect(inv.__getitem__(i).eq(1)).toBe(true);
    }
  });

  test('division', () => {
    // (1 + x)/(1 - x) with precision 5
    const numer = R.__call__([1, 1]).add_bigoh(5);
    const denom = R.__call__([1, -1]).add_bigoh(5);
    const quot = numer.div(denom);

    // (1+x)/(1-x) = (1+x)(1+x+x^2+...) = 1 + 2x + 2x^2 + 2x^3 + ...
    expect(quot.__getitem__(0).eq(1)).toBe(true);
    expect(quot.__getitem__(1).eq(2)).toBe(true);
    expect(quot.__getitem__(2).eq(2)).toBe(true);
    expect(quot.__getitem__(3).eq(2)).toBe(true);
  });
});

describe('PowerSeriesElement truncation and precision', () => {
  const R = new PowerSeriesRing(QQ, 'x', 20);

  test('add_bigoh', () => {
    const f = R.__call__([1, 2, 3, 4, 5]);
    const g = f.add_bigoh(3);

    expect(g.prec()).toBe(3);
    expect(g.__getitem__(0).eq(1)).toBe(true);
    expect(g.__getitem__(1).eq(2)).toBe(true);
    expect(g.__getitem__(2).eq(3)).toBe(true);
  });

  test('truncate', () => {
    const f = R.__call__([1, 2, 3, 4, 5]).add_bigoh(10);
    const g = f.truncate(3);

    expect(g.prec()).toBe(Number.POSITIVE_INFINITY); // truncate returns a polynomial
    expect(g.degree()).toBe(2);
  });

  test('valuation', () => {
    expect(R.__call__([0, 0, 1, 2]).valuation()).toBe(2);
    expect(R.__call__([1, 2, 3]).valuation()).toBe(0);
    expect(R.zero().valuation()).toBe(Number.POSITIVE_INFINITY);
  });

  test('degree', () => {
    expect(R.__call__([1, 2, 3]).degree()).toBe(2);
    expect(R.__call__([1]).degree()).toBe(0);
    expect(R.zero().degree()).toBe(-1);
  });
});

describe('PowerSeriesElement derivative and integral', () => {
  const R = new PowerSeriesRing(QQ, 'x', 20);

  test('derivative', () => {
    // d/dx (1 + 2x + 3x^2) = 2 + 6x
    const f = R.__call__([1, 2, 3]);
    const df = f.derivative();

    expect(df.__getitem__(0).eq(2)).toBe(true);
    expect(df.__getitem__(1).eq(6)).toBe(true);
  });

  test('integral', () => {
    // integral of (1 + 2x) = x + x^2
    const f = R.__call__([1, 2]);
    const intf = f.integral();

    expect(intf.__getitem__(0).isZero()).toBe(true);
    expect(intf.__getitem__(1).eq(1)).toBe(true);
    expect(intf.__getitem__(2).eq(1)).toBe(true);
  });

  test('derivative of integral is identity', () => {
    const f = R.__call__([1, 2, 3, 4]);
    const result = f.integral().derivative();

    expect(result.__getitem__(0).eq(1)).toBe(true);
    expect(result.__getitem__(1).eq(2)).toBe(true);
    expect(result.__getitem__(2).eq(3)).toBe(true);
    expect(result.__getitem__(3).eq(4)).toBe(true);
  });
});

describe('PowerSeriesElement exp and log', () => {
  const R = new PowerSeriesRing(QQ, 'x', 10);
  const x = R.gen();

  test('exp of x', () => {
    // exp(x) = 1 + x + x^2/2 + x^3/6 + ...
    const expX = x.add_bigoh(5).exp();

    expect(expX.__getitem__(0).eq(1)).toBe(true);
    expect(expX.__getitem__(1).eq(1)).toBe(true);
    // x^2 coefficient is 1/2
    expect(expX.__getitem__(2).value.eq(new Rational(1n, 2n))).toBe(true);
    // x^3 coefficient is 1/6
    expect(expX.__getitem__(3).value.eq(new Rational(1n, 6n))).toBe(true);
    // x^4 coefficient is 1/24
    expect(expX.__getitem__(4).value.eq(new Rational(1n, 24n))).toBe(true);
  });

  test('log of 1+x', () => {
    // log(1+x) = x - x^2/2 + x^3/3 - x^4/4 + ...
    const f = R.__call__([1, 1]).add_bigoh(5);
    const logF = f.log();

    expect(logF.__getitem__(0).isZero()).toBe(true);
    expect(logF.__getitem__(1).eq(1)).toBe(true);
    // x^2 coefficient is -1/2
    expect(logF.__getitem__(2).value.eq(new Rational(-1n, 2n))).toBe(true);
    // x^3 coefficient is 1/3
    expect(logF.__getitem__(3).value.eq(new Rational(1n, 3n))).toBe(true);
  });

  test('exp(log(1+x)) = 1+x', () => {
    const f = R.__call__([1, 1]).add_bigoh(6);
    const result = f.log().exp();

    // Should recover 1 + x
    expect(result.__getitem__(0).eq(1)).toBe(true);
    expect(result.__getitem__(1).eq(1)).toBe(true);
    // Higher coefficients should be zero (or very small)
    for (let i = 2; i < 5; i++) {
      expect(result.__getitem__(i).isZero()).toBe(true);
    }
  });
});

describe('PowerSeriesElement sqrt', () => {
  const R = new PowerSeriesRing(QQ, 'x', 10);

  test('sqrt of 1+x', () => {
    // sqrt(1+x) = 1 + x/2 - x^2/8 + x^3/16 - ...
    const f = R.__call__([1, 1]).add_bigoh(5);
    const sqrtF = f.sqrt();

    expect(sqrtF.__getitem__(0).eq(1)).toBe(true);
    // x coefficient is 1/2
    expect(sqrtF.__getitem__(1).value.eq(new Rational(1n, 2n))).toBe(true);
    // x^2 coefficient is -1/8
    expect(sqrtF.__getitem__(2).value.eq(new Rational(-1n, 8n))).toBe(true);
  });

  test('sqrt squared is identity', () => {
    const f = R.__call__([1, 1]).add_bigoh(6);
    const sqrtF = f.sqrt();
    const result = sqrtF.mul(sqrtF);

    expect(result.__getitem__(0).eq(1)).toBe(true);
    expect(result.__getitem__(1).eq(1)).toBe(true);
    // Higher coefficients should be zero
    for (let i = 2; i < 5; i++) {
      expect(result.__getitem__(i).isZero()).toBe(true);
    }
  });
});

describe('PowerSeriesElement composition', () => {
  const R = new PowerSeriesRing(QQ, 'x', 10);
  const x = R.gen();

  test('composition f(g) where g = x', () => {
    // f(x) = 1 + x + x^2
    const f = R.__call__([1, 1, 1]).add_bigoh(5);
    // f(x) should equal f
    const result = f.__call__(x.add_bigoh(5));

    expect(result.__getitem__(0).eq(1)).toBe(true);
    expect(result.__getitem__(1).eq(1)).toBe(true);
    expect(result.__getitem__(2).eq(1)).toBe(true);
  });

  test('composition f(g) where g = 2x', () => {
    // f(x) = 1 + x
    const f = R.__call__([1, 1]).add_bigoh(5);
    // g(x) = 2x
    const g = R.__call__([0, 2]).add_bigoh(5);
    // f(g(x)) = 1 + 2x
    const result = f.__call__(g);

    expect(result.__getitem__(0).eq(1)).toBe(true);
    expect(result.__getitem__(1).eq(2)).toBe(true);
  });
});

describe('PowerSeriesElement reversion', () => {
  const R = new PowerSeriesRing(QQ, 'x', 10);

  test('reversion of x - x^2', () => {
    // f(x) = x - x^2
    // Catalan numbers: the reversion is x + x^2 + 2x^3 + 5x^4 + 14x^5 + ...
    const f = R.__call__([0, 1, -1]).add_bigoh(6);
    const g = f.reversion();

    expect(g.__getitem__(1).eq(1)).toBe(true);
    expect(g.__getitem__(2).eq(1)).toBe(true);
    expect(g.__getitem__(3).eq(2)).toBe(true);
    expect(g.__getitem__(4).eq(5)).toBe(true);
  });

  test('reversion of 2x', () => {
    // f(x) = 2x, f^{-1}(x) = x/2
    const f = R.__call__([0, 2]).add_bigoh(5);
    const g = f.reversion();

    expect(g.__getitem__(1).value.eq(new Rational(1n, 2n))).toBe(true);
  });
});

describe('PowerSeriesElement toString', () => {
  const R = new PowerSeriesRing(QQ, 't', 10);

  test('string representation of polynomial', () => {
    const f = R.__call__([1, 2, 3]);
    expect(f.toString()).toContain('1');
    expect(f.toString()).toContain('t');
    expect(f.toString()).toContain('t^2');
  });

  test('string representation with precision', () => {
    const f = R.__call__([1, 2]).add_bigoh(5);
    expect(f.toString()).toContain('O(t^5)');
  });

  test('string representation of zero', () => {
    expect(R.zero().toString()).toBe('0');
  });

  test('string representation of O(t^n)', () => {
    const f = R.__call__([]).add_bigoh(3);
    expect(f.toString()).toBe('O(t^3)');
  });
});

describe('LaurentSeriesRing', () => {
  const R = new PowerSeriesRing(QQ, 'x', 10);

  test('laurent_series_ring', () => {
    const L = R.laurent_series_ring();

    expect(L.variable_name()).toBe('x');
    // Use toString() comparison since laurent_series_ring creates a new reference
    expect(L.power_series_ring().toString()).toBe(R.toString());
  });

  test('coercion from power series', () => {
    const L = R.laurent_series_ring();
    const f = R.__call__([1, 2, 3]);
    const g = L.__call__(f);

    expect(g.valuation()).toBe(0);
  });
});

// ---------------------------------------------------------------------------
// Value assertions ported from SageMath doctests.
//
// The log tests above use only 1+x, the single input where the (previously
// wrong) recurrence coincides with the truth; pade and nth_root were never
// exercised at all.
// ---------------------------------------------------------------------------

describe('PowerSeriesElement precision (SageMath doctests)', () => {
  const R = new PowerSeriesRing(QQ, 't', 20);

  test('valuation of O(t^r) is r, only the exact zero is infinite', () => {
    // sage: O(t^7).valuation() -> 7 ; R(0).valuation() -> +Infinity
    expect(R.zero().add_bigoh(7).valuation()).toBe(7);
    expect(R.zero().valuation()).toBe(Number.POSITIVE_INFINITY);
    expect(R.__call__([5, 0, 0, 0, 0, 0, 0, 0, -1]).add_bigoh(11).valuation()).toBe(0);
    expect(R.__call__([0, 0, 0, 0, 0, 0, 0, 0, -1]).add_bigoh(11).valuation()).toBe(8);
  });

  test('precision_absolute and precision_relative', () => {
    // sage: (t^2 + O(t^3)).precision_absolute() -> 3, precision_relative() -> 1
    const f = R.__call__([0, 0, 1]).add_bigoh(3);
    expect(f.precision_absolute()).toBe(3);
    expect(f.precision_relative()).toBe(1);
    expect(R.zero().add_bigoh(4).precision_relative()).toBe(0);
    expect(R.__call__([1, 0, -1]).add_bigoh(100).precision_relative()).toBe(100);
  });
});

describe('PowerSeriesElement log (SageMath doctests)', () => {
  const R10 = new PowerSeriesRing(QQ, 't', 10);
  const t = R10.gen();

  test('log(1+t) to default precision', () => {
    expect(R10.__call__([1, 1]).add_bigoh(10).log().toString()).toBe(
      't - 1/2*t^2 + 1/3*t^3 - 1/4*t^4 + 1/5*t^5 - 1/6*t^6 + 1/7*t^7 - 1/8*t^8 + 1/9*t^9 + O(t^10)'
    );
  });

  test('t.exp().log() and (1+t).log().exp()', () => {
    expect(t.add_bigoh(10).exp().log().toString()).toBe('t + O(t^10)');
    expect(R10.__call__([1, 1]).log().exp().toString()).toBe('1 + t + O(t^10)');
  });

  test('log of a dense series', () => {
    // log(1 + x + x^2) = log((1-x^3)/(1-x))
    const R = new PowerSeriesRing(QQ, 'x', 20);
    expect(R.__call__([1, 1, 1]).log(6).toString()).toBe(
      'x + 1/2*x^2 - 2/3*x^3 + 1/4*x^4 + 1/5*x^5 + O(x^6)'
    );
  });

  test('log/exp round trip on a dense series', () => {
    const R = new PowerSeriesRing(QQ, 'x', 20);
    const f = R.__call__([1, 3, -2, 5, 7, -1, 4, 2]);
    expect(f.log(8).exp(8).toString()).toBe(f.add_bigoh(8).toString());
  });

  test('log requires constant term 1', () => {
    expect(() => R10.__call__([-1, 1]).add_bigoh(10).log()).toThrow(
      'constant term of power series is not 1'
    );
  });
});

describe('PowerSeriesElement nth_root (SageMath doctests)', () => {
  const R = new PowerSeriesRing(QQ, 'x', 20);

  test('(1+x).nth_root(5)', () => {
    expect(R.__call__([1, 1]).add_bigoh(5).nth_root(5).toString()).toBe(
      '1 + 1/5*x - 2/25*x^2 + 6/125*x^3 - 21/625*x^4 + O(x^5)'
    );
  });

  test('exact roots show infinite precision', () => {
    // sage: ((1+x)^5).nth_root(5) -> 1 + x
    expect(R.__call__([1, 1]).pow(5).nth_root(5).toString()).toBe('1 + x');
  });

  test('precision on O(x^r)', () => {
    expect(R.zero().add_bigoh(4).nth_root(2).toString()).toBe('O(x^2)');
    expect(R.zero().add_bigoh(4).nth_root(3).toString()).toBe('O(x^1)');
    expect(R.zero().add_bigoh(4).nth_root(4).toString()).toBe('O(x^1)');
  });

  test('precision on higher valuation series', () => {
    // sage: (x^5+x^6+O(x^7)).nth_root(5) -> x + 1/5*x^2 + O(x^3)
    expect(R.__call__([0, 0, 0, 0, 0, 1, 1]).add_bigoh(7).nth_root(5).toString()).toBe(
      'x + 1/5*x^2 + O(x^3)'
    );
  });

  test('consistent with taking log and exponential', () => {
    // sage: p = (1 + 2*x - x^4)**200 ; p.nth_root(1000, prec) == (p.log()/1000).exp()
    const R40 = new PowerSeriesRing(QQ, 'x', 40);
    const p = R40.__call__([1, 2, 0, 0, -1]).add_bigoh(30).pow(200);
    const p1 = p.nth_root(1000, 30);
    const p2 = p.log(30).scalar_div(new RationalElement(1000n)).exp(30);
    expect(p1.prec()).toBe(30);
    expect(p2.prec()).toBe(30);
    expect(p1.toString()).toBe(p2.toString());
    expect(p1.pow(1000).add_bigoh(30).toString()).toBe(p.toString());
  });

  test('cube root round trip (n >= 3 used to raise a TypeError)', () => {
    const f = R.__call__([1, 4, 5, 6]).add_bigoh(10);
    expect(f.nth_root(3).pow(3).add_bigoh(10).toString()).toBe(f.toString());
  });

  test('nth_root(2) agrees with sqrt', () => {
    const f = R.__call__([1, 1]).add_bigoh(8);
    expect(f.nth_root(2).toString()).toBe(f.sqrt(8).toString());
  });

  test('errors', () => {
    expect(() => R.__call__([0, 1, 1]).add_bigoh(8).nth_root(2)).toThrow(
      'power series valuation is not a multiple of 2'
    );
    expect(() => R.__call__([1, 1]).nth_root(-3)).toThrow('n (=-3) must be positive');
    expect(() => R.__call__([1, 1]).nth_root(0)).toThrow('n (=0) must be positive');
  });
});

describe('PowerSeriesElement pade (SageMath doctests)', () => {
  const Rz = new PowerSeriesRing(QQ, 'z', 20);
  const z = Rz.gen();
  const expz = z.add_bigoh(12).exp(12);

  test('exp(z).pade(4, 0)', () => {
    expect(expz.pade(4, 0).toString()).toBe('1/24*z^4 + 1/6*z^3 + 1/2*z^2 + z + 1');
  });

  test('exp(z).pade(1, 1)', () => {
    expect(expz.pade(1, 1).toString()).toBe('(-z - 2)/(z - 2)');
  });

  test('exp(z).pade(3, 3)', () => {
    expect(expz.pade(3, 3).toString()).toBe(
      '(-z^3 - 12*z^2 - 60*z - 120)/(z^3 - 12*z^2 + 60*z - 120)'
    );
  });

  test('log(1-z).pade(4, 4)', () => {
    expect(Rz.__call__([1, -1]).log(12).pade(4, 4).toString()).toBe(
      '(25/6*z^4 - 130/3*z^3 + 105*z^2 - 70*z)/(z^4 - 20*z^3 + 90*z^2 - 140*z + 70)'
    );
  });

  test('sqrt(1+z).pade(3, 2)', () => {
    expect(Rz.__call__([1, 1]).sqrt(12).pade(3, 2).toString()).toBe(
      '(1/6*z^3 + 3*z^2 + 8*z + 16/3)/(z^2 + 16/3*z + 16/3)'
    );
  });

  test('exp(2*z).pade(3, 3)', () => {
    expect(Rz.__call__([0, 2]).add_bigoh(12).exp(12).pade(3, 3).toString()).toBe(
      '(-z^3 - 6*z^2 - 15*z - 15)/(z^3 - 6*z^2 + 15*z - 15)'
    );
  });

  test('trac 21212 and correct precision', () => {
    const Rx = new PowerSeriesRing(QQ, 'x', 20);
    expect(Rx.__call__([1, 1]).add_bigoh(100).pade(2, 2).toString()).toBe('x + 1');
    expect(Rx.__call__([1, 1]).add_bigoh(2).pade(0, 1).toString()).toBe('-1/(x - 1)');
  });

  test('too low precision raises', () => {
    const Rx = new PowerSeriesRing(QQ, 'x', 20);
    expect(() => Rx.__call__([0, 1]).add_bigoh(6).pade(4, 4)).toThrow(
      'the precision of the series is not large enough'
    );
  });

  test('f - Q/P = O(z^(m+n+1))', () => {
    const approx = expz.pade(3, 3);
    expect(expz.sub(approx.power_series(12)).add_bigoh(7).is_zero()).toBe(true);
    expect(approx.numerator().degree()).toBeLessThanOrEqual(3);
    expect(approx.denominator().degree()).toBeLessThanOrEqual(3);
  });
});

// ---------------------------------------------------------------------------
// Multivariate power series.
//
// Value assertions are SageMath doctests from
// `sage/rings/multi_power_series_ring_element.py`, copied verbatim.
// ---------------------------------------------------------------------------

/** The integers, as a coefficient ring (SageMath's ZZ). */
class IntElement implements RingElement {
  readonly value: bigint;
  constructor(value: bigint | number) {
    this.value = typeof value === 'bigint' ? value : BigInt(value);
  }
  add(o: IntElement): IntElement {
    return new IntElement(this.value + o.value);
  }
  sub(o: IntElement): IntElement {
    return new IntElement(this.value - o.value);
  }
  mul(o: IntElement): IntElement {
    return new IntElement(this.value * o.value);
  }
  div(o: IntElement): IntElement {
    if (o.value === 0n || this.value % o.value !== 0n) {
      throw new TypeError(`${this.value}/${o.value} is not an integer`);
    }
    return new IntElement(this.value / o.value);
  }
  neg(): IntElement {
    return new IntElement(-this.value);
  }
  eq(o: IntElement | number | bigint): boolean {
    if (typeof o === 'number') return this.value === BigInt(o);
    if (typeof o === 'bigint') return this.value === o;
    return this.value === o.value;
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
  toString(): string {
    return this.value.toString();
  }
}

const ZZring: CoefficientRing<IntElement> = {
  zero: () => new IntElement(0n),
  one: () => new IntElement(1n),
  __call__(x: unknown): IntElement {
    if (x instanceof IntElement) return x;
    if (typeof x === 'bigint') return new IntElement(x);
    if (typeof x === 'number') return new IntElement(BigInt(x));
    throw new Error('cannot convert to IntElement');
  },
  is_field: () => false,
  characteristic: () => 0n,
  toString: () => 'Integer Ring',
};

/** Build a multivariate series from [exponents, coefficient] pairs. */
function mps<T extends RingElement>(
  R: MPowerSeriesRing<T>,
  terms: [number[], number | bigint | string][],
  prec = Number.POSITIVE_INFINITY
): MPowerSeries<T> {
  return new MPowerSeries<T>(
    R,
    terms.map(([e, c]) => [e, R.base_ring().__call__(typeof c === 'number' ? BigInt(c) : c)]),
    prec
  );
}

describe('MPowerSeriesRing (sage/rings/multi_power_series_ring.py)', () => {
  test('repr of the ring and of elements', () => {
    // sage: R.<s,t> = PowerSeriesRing(ZZ); R
    // Multivariate Power Series Ring in s, t over Integer Ring
    const R = new MPowerSeriesRing(ZZring, 's,t');
    expect(R.toString()).toBe('Multivariate Power Series Ring in s, t over Integer Ring');
    expect(R.ngens()).toBe(2);
    expect(R.default_prec()).toBe(10);
    // sage: f = 1 + t + s + s*t + R.O(3); f
    // 1 + s + t + s*t + O(s, t)^3
    const [s, t] = R.gens() as [MPowerSeries<IntElement>, MPowerSeries<IntElement>];
    const f = R.one().add(t).add(s).add(s.mul(t)).add(R.O(3));
    expect(f.toString()).toBe('1 + s + t + s*t + O(s, t)^3');

    // sage: S.<s,t> = PowerSeriesRing(ZZ); f = s + 4*t + 3*s*t
    // sage: f.add_bigoh(4) -> s + 4*t + 3*s*t + O(s, t)^4
    const g = mps(R, [
      [[1, 0], 1],
      [[0, 1], 4],
      [[1, 1], 3],
    ]);
    expect(g.add_bigoh(4).toString()).toBe('s + 4*t + 3*s*t + O(s, t)^4');
    // sage: g = 1 + s + t - s*t + S.O(5); g -> 1 + s + t - s*t + O(s, t)^5
    const h = mps(
      R,
      [
        [[0, 0], 1],
        [[1, 0], 1],
        [[0, 1], 1],
        [[1, 1], -1],
      ],
      5
    );
    expect(h.toString()).toBe('1 + s + t - s*t + O(s, t)^5');
  });

  test('repr with three variables over QQ', () => {
    // sage: B.<s,t,v> = PowerSeriesRing(QQ)
    // sage: e = 1 + s - s*t + t*v/2 - 2*s*t*v/8 + B.O(4)
    // sage: e._repr_() -> '1 + s - s*t + 1/2*t*v - 1/4*s*t*v + O(s, t, v)^4'
    const B = new MPowerSeriesRing(QQ, 's,t,v');
    const e = mps(
      B,
      [
        [[0, 0, 0], '1'],
        [[1, 0, 0], '1'],
        [[1, 1, 0], '-1'],
        [[0, 1, 1], '1/2'],
        [[1, 1, 1], '-1/4'],
      ],
      4
    );
    expect(e.toString()).toBe('1 + s - s*t + 1/2*t*v - 1/4*s*t*v + O(s, t, v)^4');
  });
});

describe('MPowerSeries arithmetic (multi_power_series_ring_element.py)', () => {
  const R = new MPowerSeriesRing(ZZring, 'a,b,c');
  // sage: f0 = -a^3*b*c^2 + a^2*b^2*c^4 - 12*a^3*b^3*c^3 + R.O(10)
  const f0 = mps(
    R,
    [
      [[3, 1, 2], -1],
      [[2, 2, 4], 1],
      [[3, 3, 3], -12],
    ],
    10
  );
  // sage: f1 = -6*b*c^3 - 4*a^2*b*c^2 + a^6*b^2*c - 2*a^3*b^3*c^3 + R.O(10)
  const f1 = mps(
    R,
    [
      [[0, 1, 3], -6],
      [[2, 1, 2], -4],
      [[6, 2, 1], 1],
      [[3, 3, 3], -2],
    ],
    10
  );

  test('_add_', () => {
    // sage: g = f0 + f1; g
    // -6*b*c^3 - 4*a^2*b*c^2 - a^3*b*c^2 + a^2*b^2*c^4 + a^6*b^2*c
    //  - 14*a^3*b^3*c^3 + O(a, b, c)^10
    expect(f0.add(f1).toString()).toBe(
      '-6*b*c^3 - 4*a^2*b*c^2 - a^3*b*c^2 + a^2*b^2*c^4 + a^6*b^2*c - 14*a^3*b^3*c^3 + O(a, b, c)^10'
    );
  });

  test('_sub_', () => {
    // sage: g = f0 - f1; g
    // 6*b*c^3 + 4*a^2*b*c^2 - a^3*b*c^2 + a^2*b^2*c^4 - a^6*b^2*c
    //  - 10*a^3*b^3*c^3 + O(a, b, c)^10
    expect(f0.sub(f1).toString()).toBe(
      '6*b*c^3 + 4*a^2*b*c^2 - a^3*b*c^2 + a^2*b^2*c^4 - a^6*b^2*c - 10*a^3*b^3*c^3 + O(a, b, c)^10'
    );
  });

  test('_mul_', () => {
    // sage: g = f0*f1; g
    // 6*a^3*b^2*c^5 + 4*a^5*b^2*c^4 - 6*a^2*b^3*c^7 - 4*a^4*b^3*c^6
    //  + 72*a^3*b^4*c^6 + O(a, b, c)^14
    expect(f0.mul(f1).toString()).toBe(
      '6*a^3*b^2*c^5 + 4*a^5*b^2*c^4 - 6*a^2*b^3*c^7 - 4*a^4*b^3*c^6 + 72*a^3*b^4*c^6 + O(a, b, c)^14'
    );
  });

  test('_lmul_', () => {
    // sage: g = 3*f0; g
    // -3*a^3*b*c^2 + 3*a^2*b^2*c^4 - 36*a^3*b^3*c^3 + O(a, b, c)^10
    expect(f0.scalar_mul(new IntElement(3n)).toString()).toBe(
      '-3*a^3*b*c^2 + 3*a^2*b^2*c^4 - 36*a^3*b^3*c^3 + O(a, b, c)^10'
    );
  });

  test('__invert__', () => {
    // sage: R.<a,b,c> = PowerSeriesRing(ZZ)
    // sage: f = 1 + a + b - a*b - b*c - a*c + R.O(4)
    // sage: ~f
    // 1 - a - b + a^2 + 3*a*b + a*c + b^2 + b*c - a^3 - 5*a^2*b
    //  - 2*a^2*c - 5*a*b^2 - 4*a*b*c - b^3 - 2*b^2*c + O(a, b, c)^4
    const f = mps(
      R,
      [
        [[0, 0, 0], 1],
        [[1, 0, 0], 1],
        [[0, 1, 0], 1],
        [[1, 1, 0], -1],
        [[0, 1, 1], -1],
        [[1, 0, 1], -1],
      ],
      4
    );
    expect(f.inv().toString()).toBe(
      '1 - a - b + a^2 + 3*a*b + a*c + b^2 + b*c - a^3 - 5*a^2*b - 2*a^2*c - 5*a*b^2 - 4*a*b*c - b^3 - 2*b^2*c + O(a, b, c)^4'
    );
    // sage: g = 1/f; g == ~f -> True ; f*g == 1
    expect(R.one().div(f).eq(f.inv())).toBe(true);
    expect(f.mul(f.inv()).eq(R.one())).toBe(true);
    // Non-unit constant term is not implemented (as in SageMath).
    expect(() => R.gen(0).inv()).toThrow(
      'Multiplicative inverse of multivariate power series currently implemented only if constant coefficient is a unit.'
    );
  });

  test('valuation, degree, is_unit, prec', () => {
    // sage: R.<a,b> = PowerSeriesRing(GF(4949717))
    // sage: f = a^2 + a*b + a^3 + R.O(9); f.valuation() -> 2
    // sage: g = 1 + a + a^3; g.valuation() -> 0
    // sage: R.zero().valuation() -> +Infinity
    const R2 = new MPowerSeriesRing(ZZring, 'a,b');
    const f = mps(
      R2,
      [
        [[2, 0], 1],
        [[1, 1], 1],
        [[3, 0], 1],
      ],
      9
    );
    expect(f.valuation()).toBe(2);
    const g = mps(R2, [
      [[0, 0], 1],
      [[1, 0], 1],
      [[3, 0], 1],
    ]);
    expect(g.valuation()).toBe(0);
    expect(R2.zero().valuation()).toBe(Number.POSITIVE_INFINITY);

    // sage: B.<x,y> = PowerSeriesRing(QQ); r = 1 - x*y + x^2
    // sage: r.add_bigoh(4) -> 1 + x^2 - x*y + O(x, y)^4
    // sage: r.add_bigoh(2) -> 1 + O(x, y)^2
    // sage: r.add_bigoh(4).degree() -> 2
    const B = new MPowerSeriesRing(QQ, 'x,y');
    const r = mps(B, [
      [[0, 0], '1'],
      [[1, 1], '-1'],
      [[2, 0], '1'],
    ]);
    expect(r.add_bigoh(4).toString()).toBe('1 + x^2 - x*y + O(x, y)^4');
    expect(r.add_bigoh(2).toString()).toBe('1 + O(x, y)^2');
    expect(r.add_bigoh(4).degree()).toBe(2);
    // sage: r -> 1 + x^2 - x*y  (add_bigoh does not change self)
    expect(r.toString()).toBe('1 + x^2 - x*y');
    // sage: f.truncate().prec() -> +Infinity
    expect(f.truncate().prec()).toBe(Number.POSITIVE_INFINITY);

    // sage: R.<a,b> = PowerSeriesRing(ZZ)
    // sage: f = 2 + a^2 + a*b + a^3 + R.O(9); f.is_unit() -> False
    // sage: (O(a,b)^0).is_unit() -> False
    const f2 = mps(
      R2,
      [
        [[0, 0], 2],
        [[2, 0], 1],
        [[1, 1], 1],
        [[3, 0], 1],
      ],
      9
    );
    expect(f2.is_unit()).toBe(false);
    expect(R2.O(0).is_unit()).toBe(false);
  });

  test('__getitem__', () => {
    // sage: R.<x,y> = QQ[[]]
    // sage: ((x+y)^3)[2,1] -> 3
    // sage: f = 1/(1 + x + y); f[2,5] -> -21
    // sage: f[0,30] -> IndexError
    const R2 = new MPowerSeriesRing(QQ, 'x,y');
    const [x, y] = R2.gens() as [MPowerSeries<RationalElement>, MPowerSeries<RationalElement>];
    expect(x.add(y).pow(3).__getitem__([2, 1]).toString()).toBe('3');
    const f = R2.one().div(R2.one().add(x).add(y));
    expect(f.__getitem__([2, 5]).toString()).toBe('-21');
    expect(() => f.__getitem__([0, 30])).toThrow(
      'Cannot return the coefficients of terms of total degree greater than or equal to precision of self.'
    );
  });

  test('__call__ substitution', () => {
    // sage: R.<s,t> = PowerSeriesRing(ZZ)
    // sage: f = s^2 + s*t + s^3 + s^2*t + 3*s^4 + 3*s^3*t + R.O(5)
    // sage: f(t,s) -> s*t + t^2 + s*t^2 + t^3 + 3*s*t^3 + 3*t^4 + O(s, t)^5
    // sage: f(t,0) -> t^2 + t^3 + 3*t^4 + O(s, t)^5
    // sage: f(t,2) -> TypeError
    // sage: f.truncate()(t,2) -> 2*t + 3*t^2 + 7*t^3 + 3*t^4
    const R2 = new MPowerSeriesRing(ZZring, 's,t');
    const [s, t] = R2.gens() as [MPowerSeries<IntElement>, MPowerSeries<IntElement>];
    const f = mps(
      R2,
      [
        [[2, 0], 1],
        [[1, 1], 1],
        [[3, 0], 1],
        [[2, 1], 1],
        [[4, 0], 3],
        [[3, 1], 3],
      ],
      5
    );
    expect(f.__call__(t, s).toString()).toBe(
      's*t + t^2 + s*t^2 + t^3 + 3*s*t^3 + 3*t^4 + O(s, t)^5'
    );
    expect(f.__call__(t, R2.zero()).toString()).toBe('t^2 + t^3 + 3*t^4 + O(s, t)^5');
    expect(() => f.__call__(t, R2.__call__(2))).toThrow(
      'Substitution defined only for elements of positive valuation, unless self has infinite precision.'
    );
    expect(f.truncate().__call__(t, R2.__call__(2)).toString()).toBe('2*t + 3*t^2 + 7*t^3 + 3*t^4');
  });
});

describe('MPowerSeries: the formal group law in three variables', () => {
  // sage: e = EllipticCurve(GF(7), [3, 4]); ehat = e.formal()
  // sage: F = ehat.group_law(7); F
  // t1 + t2 + t1^4*t2 + 2*t1^3*t2^2 + 2*t1^2*t2^3 + t1*t2^4 + O(t1, t2)^7
  // sage: R.<x,y,z> = GF(7)[[]]
  // sage: F(x, ehat.inverse()(x)) -> 0 + O(x, y, z)^7
  // sage: F(x, y) == F(y, x) -> True
  // sage: F(x, F(y, z)) == F(F(x, y), z) -> True
  const K = GF(7n) as unknown as CoefficientRing<RingElement>;
  const prec = 7;
  const E = EllipticCurve(GF(7n) as never, [3n, 4n] as never) as never as {
    formal_group(): {
      group_law(prec: number): { terms: Map<string, RingElement> };
      inverse(prec: number): { list(): RingElement[] };
    };
  };
  const ehat = E.formal_group();
  const law = ehat.group_law(prec);

  const R2 = new MPowerSeriesRing<RingElement>(K, 't1,t2');
  const R3 = new MPowerSeriesRing<RingElement>(K, 'x,y,z');

  /** The group law as an element of GF(7)[[t1,t2]]. */
  const F = new MPowerSeries<RingElement>(
    R2,
    [...law.terms].map(([key, c]) => [key.split(',').map(Number), c] as [number[], RingElement]),
    prec
  );

  const [x, y, z] = R3.gens();

  test('the group law prints as in the SageMath doctest', () => {
    expect(F.toString()).toBe(
      't1 + t2 + t1^4*t2 + 2*t1^3*t2^2 + 2*t1^2*t2^3 + t1*t2^4 + O(t1, t2)^7'
    );
  });

  test('F(x, y) == F(y, x)', () => {
    expect(F.__call__(x!, y!).eq(F.__call__(y!, x!))).toBe(true);
  });

  test('F(x, F(y, z)) == F(F(x, y), z)', () => {
    const lhs = F.__call__(x!, F.__call__(y!, z!));
    const rhs = F.__call__(F.__call__(x!, y!), z!);
    expect(lhs.eq(rhs)).toBe(true);
    expect(lhs.prec()).toBe(7);
    // The identity is not vacuous: both sides really are t1 + t2 + ... in 3 vars.
    expect(lhs.__getitem__([1, 0, 0]).toString()).toBe('1');
    expect(lhs.__getitem__([0, 1, 0]).toString()).toBe('1');
    expect(lhs.__getitem__([0, 0, 1]).toString()).toBe('1');
  });

  test('F(x, i(x)) == 0 + O(x, y, z)^7', () => {
    const inv = ehat.inverse(prec).list();
    const ix = new MPowerSeries<RingElement>(
      R3,
      inv.map((c, k) => [[k, 0, 0], c] as [number[], RingElement]),
      prec
    );
    expect(F.__call__(x!, ix).toString()).toBe('0 + O(x, y, z)^7');
  });
});
