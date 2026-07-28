/**
 * Tests for polynomial_element.ts over ZZ, QQ and finite fields.
 *
 * The expected values are SageMath's:
 *
 *     sage: R.<x> = ZZ[]
 *     sage: (x^3+x+1).resultant(x^3-x-1)
 *     -8
 *     sage: (x^3+x+1).discriminant()
 *     -31
 *     sage: (2*x^3+x+1).discriminant()
 *     -116
 *     sage: (x^2-4*x+1).discriminant()
 *     12
 *     sage: (x^4+x+1).discriminant()
 *     229
 */
import { describe, expect, test } from 'bun:test';
import { ArithmeticError, NotImplementedError } from '../../errors.js';
import { FiniteFieldPrime } from '../finite_rings/finite_field_prime.js';
import { GF2 } from '../finite_rings/gf2.js';
import { type CoefficientRing, Polynomial, type RingElement } from './polynomial_element.js';
import { PolynomialRingConstructor } from './polynomial_ring.js';

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

/** ZZ element (integral domain, truncating division). */
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
  __call__(x: unknown): IntegerElement {
    if (x instanceof IntegerElement) return x;
    if (typeof x === 'bigint') return new IntegerElement(x);
    if (typeof x === 'number') return new IntegerElement(BigInt(x));
    throw new Error(`cannot convert ${typeof x} to IntegerElement`);
  },
  is_field: () => false,
  toString: () => 'Integer Ring',
};

/** QQ element. */
class RationalElement implements RingElement {
  readonly numerator: bigint;
  readonly denominator: bigint;
  constructor(n: bigint | number, d: bigint | number = 1n) {
    let num = typeof n === 'number' ? BigInt(n) : n;
    let den = typeof d === 'number' ? BigInt(d) : d;
    if (den === 0n) throw new Error('zero denominator');
    if (den < 0n) {
      num = -num;
      den = -den;
    }
    const g = bigintGcd(num, den) || 1n;
    this.numerator = num / g;
    this.denominator = den / g;
  }
  add(o: RationalElement): RationalElement {
    return new RationalElement(
      this.numerator * o.denominator + o.numerator * this.denominator,
      this.denominator * o.denominator
    );
  }
  sub(o: RationalElement): RationalElement {
    return new RationalElement(
      this.numerator * o.denominator - o.numerator * this.denominator,
      this.denominator * o.denominator
    );
  }
  mul(o: RationalElement): RationalElement {
    return new RationalElement(this.numerator * o.numerator, this.denominator * o.denominator);
  }
  div(o: RationalElement): RationalElement {
    return new RationalElement(this.numerator * o.denominator, this.denominator * o.numerator);
  }
  inv(): RationalElement {
    return new RationalElement(this.denominator, this.numerator);
  }
  neg(): RationalElement {
    return new RationalElement(-this.numerator, this.denominator);
  }
  eq(o: RationalElement | number): boolean {
    if (typeof o === 'number') return this.denominator === 1n && this.numerator === BigInt(o);
    return this.numerator === o.numerator && this.denominator === o.denominator;
  }
  isZero(): boolean {
    return this.numerator === 0n;
  }
  toString(): string {
    return this.denominator === 1n
      ? this.numerator.toString()
      : `${this.numerator}/${this.denominator}`;
  }
}

const QQ: CoefficientRing<RationalElement> = {
  zero: () => new RationalElement(0n),
  one: () => new RationalElement(1n),
  __call__(x: unknown): RationalElement {
    if (x instanceof RationalElement) return x;
    if (typeof x === 'bigint' || typeof x === 'number') return new RationalElement(x as bigint);
    if (x && typeof x === 'object' && 'numer' in x && 'denom' in x) {
      const r = x as { numer: bigint; denom: bigint };
      return new RationalElement(r.numer, r.denom);
    }
    throw new Error(`cannot convert ${typeof x} to RationalElement`);
  },
  is_field: () => true,
  toString: () => 'Rational Field',
};

const [R, x] = PolynomialRingConstructor(ZZ, 'x');
const zz = (n: bigint | number) => new IntegerElement(n);
const zc = (n: bigint | number) => R.__call__(zz(n));

describe('resultant and discriminant over ZZ (H9)', () => {
  test("Sage's doctest values", () => {
    const f = x.pow(3).add(x).add(R.one()); // x^3 + x + 1
    const g = x.pow(3).sub(x).sub(R.one()); // x^3 - x - 1

    expect(f.resultant(g).eq(-8)).toBe(true);
    expect(f.discriminant().eq(-31)).toBe(true);

    // 2*x^3 + x + 1 -> -116  (non-monic: the resultant must be divided exactly)
    const h = x.pow(3).scalar_mul(zz(2)).add(x).add(R.one());
    expect(h.discriminant().eq(-116)).toBe(true);

    // x^2 - 4x + 1 -> 12
    expect(
      x
        .pow(2)
        .sub(x.scalar_mul(zz(4)))
        .add(R.one())
        .discriminant()
        .eq(12)
    ).toBe(true);

    // x^4 + x + 1 -> 229
    expect(x.pow(4).add(x).add(R.one()).discriminant().eq(229)).toBe(true);
  });

  test('resultant is symmetric up to sign and vanishes on common factors', () => {
    const f = x.pow(2).sub(R.one()); // (x-1)(x+1)
    const g = x
      .pow(2)
      .add(x.scalar_mul(zz(2)))
      .add(R.one()); // (x+1)^2
    expect(f.resultant(g).isZero()).toBe(true);

    // res(f, g) = (-1)^(deg f * deg g) res(g, f)
    const a = x
      .pow(3)
      .add(x.scalar_mul(zz(2)))
      .add(zc(5));
    const b = x.pow(2).sub(x).add(zc(7));
    expect(a.resultant(b).eq(b.resultant(a))).toBe(true); // 3*2 even
  });

  test('discriminant over QQ', () => {
    const [RQ, y] = PolynomialRingConstructor(QQ, 'y');
    // sage: R.<y> = QQ[]; (y^3+y+1).discriminant() == -31
    const f = y.pow(3).add(y).add(RQ.one());
    expect(f.discriminant().eq(-31)).toBe(true);
  });
});

describe('quo_rem exactness (H10)', () => {
  test('non-exact division over ZZ raises ArithmeticError', () => {
    const a = x.pow(2).sub(R.one());
    const b = x
      .pow(2)
      .scalar_mul(zz(2))
      .add(x.scalar_mul(zz(4)))
      .add(zc(2));
    expect(() => a.quo_rem(b)).toThrow(ArithmeticError);
    expect(() => a.quo_rem(b)).toThrow('division non exact');
  });

  test('exact division over ZZ still works', () => {
    const f = x.pow(2).sub(R.one());
    const [q, r] = f.quo_rem(x.sub(R.one()));
    expect(q.eq(x.add(R.one()))).toBe(true);
    expect(r.isZero()).toBe(true);
  });

  test('division by zero raises ZeroDivisionError with Sage message', () => {
    expect(() => x.quo_rem(R.zero())).toThrow('division by zero polynomial');
  });
});

describe('gcd / xgcd (H10, H11, L12, M14)', () => {
  test('gcd over ZZ terminates and matches Sage', () => {
    const a = x.pow(2).sub(R.one()); // x^2-1
    const b = x
      .pow(2)
      .add(x.scalar_mul(zz(2)))
      .add(R.one()); // x^2+2x+1
    expect(a.gcd(b).toString()).toBe(x.add(R.one()).toString());

    // sage: (2*x+2).gcd(4*x+4) == 2*x + 2  (content is part of the ZZ[x] gcd)
    const c = x.scalar_mul(zz(2)).add(zc(2));
    const d = x.scalar_mul(zz(4)).add(zc(4));
    expect(c.gcd(d).eq(c)).toBe(true);

    // sage: R.zero().gcd(-2*x-2) == 2*x + 2  (positive leading coefficient)
    expect(R.zero().gcd(c.neg()).eq(c)).toBe(true);
  });

  test('gcd over GF(p)', () => {
    const F7 = new FiniteFieldPrime(7n);
    const [R7, t] = PolynomialRingConstructor(F7, 't');
    const f = t.sub(R7.one()).mul(t.sub(R7.__call__(F7.__call__(3))));
    const g = t.sub(R7.one()).mul(t.add(R7.one()));
    const gcd = f.gcd(g);
    expect(gcd.eq(t.sub(R7.one()))).toBe(true);
    expect(gcd.is_monic()).toBe(true);
    // sage: (2*t^2).gcd(2*t) == t
    expect(
      t
        .pow(2)
        .scalar_mul(F7.__call__(2))
        .gcd(t.scalar_mul(F7.__call__(2)))
        .eq(t)
    ).toBe(true);
  });

  test('xgcd satisfies the Bezout identity and Sage zero conventions', () => {
    const [RQ, y] = PolynomialRingConstructor(QQ, 'y');
    // sage: F = (x^2+2)*x^3; G = (x^2+2)*(x-3)
    // sage: QQ._xgcd_univariate_polynomial(F, G) == (x^2 + 2, 1/27, -1/27*x^2 - 1/9*x - 1/3)
    const F = y
      .pow(2)
      .add(RQ.__call__(new RationalElement(2n)))
      .mul(y.pow(3));
    const G = y
      .pow(2)
      .add(RQ.__call__(new RationalElement(2n)))
      .mul(y.sub(RQ.__call__(new RationalElement(3n))));
    const [g, u, v] = F.xgcd(G);
    expect(g.toString()).toBe('y^2 + 2');
    expect(u.toString()).toBe('1/27');
    expect(v.toString()).toBe('(-1/27)*y^2 + (-1/9)*y + -1/3');
    expect(u.mul(F).add(v.mul(G)).eq(g)).toBe(true);

    // sage: zero.xgcd(zero) == (0, 0, 0)
    const [g0, u0, v0] = RQ.zero().xgcd(RQ.zero());
    expect(g0.isZero() && u0.isZero() && v0.isZero()).toBe(true);

    // sage: zero.xgcd(2*x) == (x, 0, 1/2) and (2*x).xgcd(zero) == (x, 1/2, 0)
    const twoY = y.scalar_mul(new RationalElement(2n));
    const [g1, u1, v1] = RQ.zero().xgcd(twoY);
    expect(g1.eq(y) && u1.isZero() && v1.toString() === '1/2').toBe(true);
    const [g2, u2, v2] = twoY.xgcd(RQ.zero());
    expect(g2.eq(y) && u2.toString() === '1/2' && v2.isZero()).toBe(true);
  });

  test('gcd/xgcd over a ring with no gcd implementation raise NotImplementedError', () => {
    const noGcdRing: CoefficientRing<IntegerElement> = {
      ...ZZ,
      toString: () => 'Some Ring',
    };
    const [Rn, t] = PolynomialRingConstructor(noGcdRing, 't');
    expect(() => t.gcd(t.add(Rn.one()))).toThrow(NotImplementedError);
    expect(() => t.xgcd(t.add(Rn.one()))).toThrow(NotImplementedError);
  });
});

describe('compose (M14)', () => {
  test('f(g) over ZZ', () => {
    // sage: f = x^2+1; g = x+1; f(g) == x^2 + 2*x + 2
    const f = x.pow(2).add(R.one());
    const g = x.add(R.one());
    const fg = f.compose(g);
    expect(fg.getCoeff(0).eq(2)).toBe(true);
    expect(fg.getCoeff(1).eq(2)).toBe(true);
    expect(fg.getCoeff(2).eq(1)).toBe(true);
  });

  test('compose agrees with evaluation', () => {
    const f = x
      .pow(3)
      .sub(x.scalar_mul(zz(2)))
      .add(zc(5));
    const g = x.pow(2).add(x);
    const fg = f.compose(g);
    for (let a = -3n; a <= 3n; a++) {
      const ga = g.evaluate(zz(a));
      expect(fg.evaluate(zz(a)).eq(f.evaluate(ga))).toBe(true);
    }
  });
});

describe('content and primitive_part (L11)', () => {
  test('content carries the sign of the leading coefficient', () => {
    // sage: (-2*x^2-4).content() == -2 ; R(-1).content() == -1
    const f = x.pow(2).scalar_mul(zz(-2)).sub(zc(4));
    expect(f.content().eq(-2)).toBe(true);
    expect(zc(-1).content().eq(-1)).toBe(true);
    expect(zc(0).content().eq(0)).toBe(true);
    // sage: (2*x^2-4*x^4+14*x^7).content() == 2
    const g = x
      .pow(2)
      .scalar_mul(zz(2))
      .sub(x.pow(4).scalar_mul(zz(4)))
      .add(x.pow(7).scalar_mul(zz(14)));
    expect(g.content().eq(2)).toBe(true);
  });

  test('primitive_part has a positive leading coefficient', () => {
    const f = x.pow(2).scalar_mul(zz(-2)).sub(zc(4));
    const pp = f.primitive_part();
    expect(pp.leading_coefficient().eq(1)).toBe(true);
    expect(pp.getCoeff(0).eq(2)).toBe(true);
    expect(pp.scalar_mul(f.content()).eq(f)).toBe(true);
  });
});

describe('derivative (L15)', () => {
  test('values', () => {
    const f = x
      .pow(3)
      .add(x.scalar_mul(zz(2)))
      .add(R.one());
    const d = f.derivative();
    expect(d.getCoeff(0).eq(2)).toBe(true);
    expect(d.getCoeff(1).eq(0)).toBe(true);
    expect(d.getCoeff(2).eq(3)).toBe(true);
  });

  test('high degree derivative is not quadratic time', () => {
    const coeffs: IntegerElement[] = [];
    for (let i = 0; i <= 3000; i++) coeffs.push(zz(1));
    const big = new Polynomial(coeffs, R);
    const start = Date.now();
    const d = big.derivative();
    expect(Date.now() - start).toBeLessThan(1000);
    expect(d.degree()).toBe(2999);
    expect(d.getCoeff(2999).eq(3000)).toBe(true);
  });
});

describe('is_irreducible (C4, L14)', () => {
  test('constants defer to the base ring', () => {
    // sage: R(5).is_irreducible() is True ; R(4) and R(1) are False
    expect(zc(5).is_irreducible()).toBe(true);
    expect(zc(-5).is_irreducible()).toBe(true);
    expect(zc(4).is_irreducible()).toBe(false);
    expect(zc(1).is_irreducible()).toBe(false);
    expect(zc(0).is_irreducible()).toBe(false);

    // over a field every nonzero constant is a unit
    const F7 = new FiniteFieldPrime(7n);
    const [R7] = PolynomialRingConstructor(F7, 't');
    expect(R7.__call__(F7.__call__(3)).is_irreducible()).toBe(false);
  });

  test('ZZ doctests', () => {
    expect(x.pow(3).add(R.one()).is_irreducible()).toBe(false);
    expect(x.pow(2).sub(R.one()).is_irreducible()).toBe(false);
    expect(x.pow(3).add(zc(2)).is_irreducible()).toBe(true);
    expect(x.scalar_mul(zz(2)).is_irreducible()).toBe(false); // not primitive
  });

  test('GF(2): reducible polynomials whose factor degree does not divide n', () => {
    const [R2, t] = PolynomialRingConstructor(GF2, 't');
    // (t^2+t+1)(t^3+t+1) = t^5+t^4+1
    expect(t.pow(5).add(t.pow(4)).add(R2.one()).is_irreducible()).toBe(false);
    expect(t.pow(5).add(t).add(R2.one()).is_irreducible()).toBe(false);
    expect(t.pow(5).add(t.pow(2)).add(R2.one()).is_irreducible()).toBe(true);
  });

  test('GF(2) exhaustive sweep of degree <= 6 against trial division', () => {
    const [R2, t] = PolynomialRingConstructor(GF2, 't');
    const toPoly = (bits: number, deg: number) => {
      let p = R2.zero();
      for (let i = 0; i <= deg; i++) if ((bits >> i) & 1) p = p.add(t.pow(i));
      return p;
    };
    for (let deg = 2; deg <= 6; deg++) {
      for (let low = 0; low < 1 << deg; low++) {
        const p = toPoly((1 << deg) | low, deg);
        let expected = true;
        for (let d = 1; 2 * d <= deg && expected; d++) {
          for (let l = 0; l < 1 << d; l++) {
            if (p.mod(toPoly((1 << d) | l, d)).isZero()) {
              expected = false;
              break;
            }
          }
        }
        expect(p.is_irreducible()).toBe(expected);
      }
    }
  });

  test('GF(5): (y^2+2)(y^3+y+1) is reducible', () => {
    const F5 = new FiniteFieldPrime(5n);
    const [R5, y] = PolynomialRingConstructor(F5, 'y');
    const f = y
      .pow(2)
      .add(R5.__call__(F5.__call__(2)))
      .mul(y.pow(3).add(y).add(R5.one()));
    expect(f.is_irreducible()).toBe(false);
  });
});

describe('factor keeps the unit (M12)', () => {
  test('over GF(5)', () => {
    const F5 = new FiniteFieldPrime(5n);
    const [R5, y] = PolynomialRingConstructor(F5, 'y');
    // sage: (2*y^2+3).factor() == (2) * (y + 1) * (y + 4)
    const f = y
      .pow(2)
      .scalar_mul(F5.__call__(2))
      .add(R5.__call__(F5.__call__(3)));
    const factors = f.factor();
    let product = R5.one();
    for (const [p, e] of factors) product = product.mul(p.pow(e));
    expect(product.eq(f)).toBe(true);
    expect(factors.some(([p]) => p.degree() === 0 && p.getCoeff(0).eq(2))).toBe(true);
  });

  test('over ZZ', () => {
    // sage: (-x^2+4).factor() == (-1) * (x - 2) * (x + 2)
    const f = x.pow(2).neg().add(zc(4));
    const factors = f.factor();
    let product = R.one();
    for (const [p, e] of factors) product = product.mul(p.pow(e));
    expect(product.eq(f)).toBe(true);
    expect(factors.some(([p]) => p.degree() === 0 && p.getCoeff(0).eq(-1))).toBe(true);
  });

  test('multiplicities over ZZ are not duplicated (H11)', () => {
    // sage: ((x-1)^2*(x-2)).factor() == (x - 2) * (x - 1)^2
    const f = x
      .sub(R.one())
      .pow(2)
      .mul(x.sub(zc(2)));
    const factors = f.factor();
    expect(factors.length).toBe(2);
    const byString = new Map(factors.map(([p, e]) => [p.toString(), e]));
    expect(byString.get('x + -1')).toBe(2);
    expect(byString.get('x + -2')).toBe(1);
  });
});

describe('roots over ZZ (M13)', () => {
  test('large constant term does not trial divide', () => {
    const start = Date.now();
    const roots = x.sub(zc(10n ** 17n)).roots();
    expect(Date.now() - start).toBeLessThan(2000);
    expect(roots.length).toBe(1);
    expect(roots[0]![0].eq(zz(10n ** 17n))).toBe(true);
  });

  test('multiplicities', () => {
    const f = x
      .sub(zc(2))
      .pow(3)
      .mul(x.add(zc(3)));
    const roots = new Map(f.roots().map(([r, m]) => [r.toString(), m]));
    expect(roots.get('2')).toBe(3);
    expect(roots.get('-3')).toBe(1);
  });
});

describe('pseudo_quo_rem (L13)', () => {
  test("Sage's doctest over ZZ", () => {
    // sage: p = x^4 + 6*x^3 + x^2 - x + 2; q = 2*x^2 - 3*x - 1
    // sage: p.pseudo_quo_rem(q) == (4*x^2 + 30*x + 51, 175*x + 67)
    const p = x
      .pow(4)
      .add(x.pow(3).scalar_mul(zz(6)))
      .add(x.pow(2))
      .sub(x)
      .add(zc(2));
    const g = x
      .pow(2)
      .scalar_mul(zz(2))
      .sub(x.scalar_mul(zz(3)))
      .sub(R.one());

    const [quo, rem] = p.pseudo_quo_rem(g);
    expect(quo.toString()).toBe('4*x^2 + 30*x + 51');
    expect(rem.toString()).toBe('175*x + 67');

    // 2^(4-2+1) * p == quo*q + rem
    expect(p.scalar_mul(zz(8)).eq(quo.mul(g).add(rem))).toBe(true);
  });

  test('the identity l^(m-n+1) f = q g + r holds in general', () => {
    const samples: Array<[typeof x, typeof x]> = [
      [x.pow(4).add(x.pow(2)).sub(zc(3)), x.pow(2).scalar_mul(zz(5)).add(x).add(R.one())],
      [x.pow(5).sub(R.one()), x.pow(3).scalar_mul(zz(-2)).add(x)],
      [x.pow(3).scalar_mul(zz(7)), x.pow(2).scalar_mul(zz(3)).sub(zc(1))],
      [x.pow(2).add(x), x.pow(4).add(R.one())],
    ];
    for (const [f, g] of samples) {
      const [quo, rem] = f.pseudo_quo_rem(g);
      const e = Math.max(0, f.degree() - g.degree() + 1);
      let scale = new IntegerElement(1n);
      for (let i = 0; i < e; i++) scale = scale.mul(g.leading_coefficient());
      expect(f.scalar_mul(scale).eq(quo.mul(g).add(rem))).toBe(true);
      expect(rem.degree()).toBeLessThan(g.degree());
    }
  });

  test('pseudo-division by zero raises ZeroDivisionError', () => {
    expect(() => x.pseudo_quo_rem(R.zero())).toThrow('Pseudo-division by zero is not possible');
  });
});

/**
 * Factoring over QQ divides each integer factor by its leading coefficient.
 * That quotient is formed with the coefficient ring's own division, so a
 * rational coefficient ring whose `__call__` only understands integers -- like
 * the one `Rational.minpoly` builds -- must work too.  (It used to call
 * `base_ring.__call__({numer, denom})`, which QQ rejected outright.)
 */
describe('factor over QQ through a minimal coefficient ring', () => {
  interface RationalCoefficientRing extends CoefficientRing<RationalElement> {
    toString(): string;
  }

  /** Accepts only Rational/bigint/number, exactly like rational.pyx's _qqRing. */
  const StrictQQ: RationalCoefficientRing = {
    zero: () => new RationalElement(0n),
    one: () => new RationalElement(1n),
    __call__(x: unknown): RationalElement {
      if (x instanceof RationalElement) return x;
      if (typeof x === 'bigint') return new RationalElement(x);
      if (typeof x === 'number') return new RationalElement(BigInt(x));
      throw new Error(`cannot coerce ${typeof x} to RationalElement`);
    },
    is_field: () => true,
    toString: () => 'Rational Field',
  };

  const [RQ, y] = PolynomialRingConstructor(StrictQQ, 'y');
  const q = (n: bigint, d: bigint = 1n) => new RationalElement(n, d);

  test('6*y^2 + y - 2 = 6 * (y - 1/2) * (y + 2/3)', () => {
    // sage: R.<y> = QQ[]; (6*y^2 + y - 2).factor()
    // (6) * (y - 1/2) * (y + 2/3)
    const f = y
      .pow(2)
      .scalar_mul(q(6n))
      .add(y)
      .sub(RQ.__call__(q(2n)));
    const factors = f.factor();
    const nonUnit = factors.filter(([g]) => g.degree() > 0);
    expect(nonUnit.map(([g]) => g.toString()).sort()).toEqual(['y + -1/2', 'y + 2/3']);
    expect(nonUnit.every(([g]) => g.is_monic())).toBe(true);

    let product = RQ.one();
    for (const [g, e] of factors) product = product.mul(g.pow(e));
    expect(product.eq(f)).toBe(true);
  });

  test('y^4 - 1 = (y - 1)(y + 1)(y^2 + 1)', () => {
    const f = y.pow(4).sub(RQ.one());
    const factors = f.factor();
    expect(factors.map(([g, e]) => `(${g})^${e}`).sort()).toEqual(
      ['(y + -1)^1', '(y + 1)^1', '(y^2 + 1)^1'].sort()
    );
    let product = RQ.one();
    for (const [g, e] of factors) product = product.mul(g.pow(e));
    expect(product.eq(f)).toBe(true);
  });
});
