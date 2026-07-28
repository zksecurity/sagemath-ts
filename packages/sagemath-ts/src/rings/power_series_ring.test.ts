/**
 * Unit tests for power series ring module
 */
import { describe, expect, test } from 'bun:test';
import {
  type CoefficientRing,
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
