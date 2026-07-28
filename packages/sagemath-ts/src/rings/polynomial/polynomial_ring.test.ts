/**
 * Tests for polynomial_ring.ts: Lagrange/Neville interpolation, divided
 * differences and cyclotomic polynomials.
 *
 * Expected values are SageMath's (`polynomial_ring.py:2295` doctests):
 *
 *     sage: R = PolynomialRing(QQ, 'x')
 *     sage: R.lagrange_polynomial([(0,1), (2,2), (3,-2), (-4,9)], algorithm='neville')
 *     [9, -11/7*x + 19/7, -17/42*x^2 - 83/42*x + 53/7, -23/84*x^3 - 11/84*x^2 + 13/7*x + 1]
 *     sage: points = [(1,-3), (2,0), (3,15), (4,48), (5,105), (6,192)]
 *     sage: R.divided_difference(points)
 *     [-3, 3, 6, 1, 0, 0]
 */
import { describe, expect, test } from 'bun:test';
import { ArithmeticError, ZeroDivisionError } from '../../errors.js';
import { FiniteFieldPrime } from '../finite_rings/finite_field_prime.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';
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

class RationalElement implements RingElement {
  readonly numerator: bigint;
  readonly denominator: bigint;
  constructor(n: bigint | number, d: bigint | number = 1n) {
    let num = typeof n === 'number' ? BigInt(n) : n;
    let den = typeof d === 'number' ? BigInt(d) : d;
    if (den === 0n) throw new ZeroDivisionError('zero denominator');
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
    if (o.isZero()) throw new ZeroDivisionError('rational division by zero');
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
  __call__(v: unknown): RationalElement {
    if (v instanceof RationalElement) return v;
    if (typeof v === 'bigint' || typeof v === 'number') return new RationalElement(v as bigint);
    throw new Error('cannot convert to rational');
  },
  is_field: () => true,
  toString: () => 'Rational Field',
};

const q = (n: bigint | number, d: bigint | number = 1n) => new RationalElement(n, d);
const [R, x] = PolynomialRingConstructor(QQ, 'x');

describe("lagrange_polynomial with algorithm='neville' returns the whole row (M17)", () => {
  const points: Array<[RationalElement, RationalElement]> = [
    [q(0), q(1)],
    [q(2), q(2)],
    [q(3), q(-2)],
    [q(-4), q(9)],
  ];

  test("Sage's doctest row", () => {
    const row = R.lagrange_polynomial(points, 'neville');
    expect(Array.isArray(row)).toBe(true);
    expect(row.length).toBe(4);
    expect(row[0]!.toString()).toBe('9');
    expect(row[1]!.toString()).toBe('(-11/7)*x + 19/7');
    expect(row[2]!.toString()).toBe('(-17/42)*x^2 + (-83/42)*x + 53/7');
    expect(row[3]!.toString()).toBe('(-23/84)*x^3 + (-11/84)*x^2 + 13/7*x + 1');
  });

  test('the last entry equals the divided-difference answer', () => {
    const row = R.lagrange_polynomial(points, 'neville');
    const p = R.lagrange_polynomial(points, 'divided_difference');
    expect(row[row.length - 1]!.eq(p)).toBe(true);
    // sage: R.lagrange_polynomial([...]) == -23/84*x^3 - 11/84*x^2 + 13/7*x + 1
    expect(p.toString()).toBe('(-23/84)*x^3 + (-11/84)*x^2 + 13/7*x + 1');
  });

  test('previous_row extends an earlier computation', () => {
    const firstTwo = R.lagrange_polynomial(points.slice(0, 2), 'neville');
    const full = R.lagrange_polynomial(points, 'neville', firstTwo);
    expect(full.length).toBe(4);
    expect(full[3]!.toString()).toBe('(-23/84)*x^3 + (-11/84)*x^2 + 13/7*x + 1');
  });

  test('empty input returns an empty row', () => {
    expect(R.lagrange_polynomial([], 'neville')).toEqual([]);
    expect(R.lagrange_polynomial([]).isZero()).toBe(true);
  });

  test("algorithm='pari' agrees with the default", () => {
    const data: Array<[RationalElement, RationalElement]> = [
      [q(0), q(1)],
      [q(2), q(5)],
      [q(3), q(10)],
    ];
    // sage: R.lagrange_polynomial(data, algorithm='pari') == x^2 + 1
    expect(R.lagrange_polynomial(data, 'pari').toString()).toBe('x^2 + 1');
  });

  test('unknown algorithm raises ValueError with Sage message', () => {
    expect(() =>
      // biome-ignore lint/suspicious/noExplicitAny: deliberately invalid argument
      R.lagrange_polynomial(points, 'abc' as any)
    ).toThrow("algorithm can be 'divided_difference', 'neville' or 'pari'");
  });

  test('positive characteristic (Sage issue 9787)', () => {
    const F101 = new FiniteFieldPrime(101n);
    const [R101] = PolynomialRingConstructor(F101, 'x');
    const pts: Array<[ReturnType<typeof F101.__call__>, ReturnType<typeof F101.__call__>]> = [
      [F101.__call__(1), F101.__call__(0)],
      [F101.__call__(2), F101.__call__(0)],
      [F101.__call__(3), F101.__call__(0)],
    ];
    expect(R101.lagrange_polynomial(pts).isZero()).toBe(true);
  });
});

describe('divided_difference (L17)', () => {
  const points: Array<[RationalElement, RationalElement]> = [
    [q(1), q(-3)],
    [q(2), q(0)],
    [q(3), q(15)],
    [q(4), q(48)],
    [q(5), q(105)],
    [q(6), q(192)],
  ];

  test('diagonal', () => {
    // sage: R.divided_difference(points) == [-3, 3, 6, 1, 0, 0]
    expect(R.divided_difference(points).map((c) => c.toString())).toEqual([
      '-3',
      '3',
      '6',
      '1',
      '0',
      '0',
    ]);
  });

  test('full_table', () => {
    // sage: R.divided_difference(points, full_table=True)
    const table = R.divided_difference(points, true).map((row) => row.map((c) => c.toString()));
    expect(table).toEqual([
      ['-3'],
      ['0', '3'],
      ['15', '15', '6'],
      ['48', '33', '9', '1'],
      ['105', '57', '12', '1', '0'],
      ['192', '87', '15', '1', '0', '0'],
    ]);
  });
});

describe('duplicate x values (L20)', () => {
  test('surface as a ZeroDivisionError from the base ring, as in Sage', () => {
    expect(() =>
      R.lagrange_polynomial([
        [q(1), q(2)],
        [q(1), q(3)],
      ])
    ).toThrow(ZeroDivisionError);
  });
});

describe('cyclotomic_polynomial (L19)', () => {
  test('values', () => {
    // sage: ZZ['x'].cyclotomic_polynomial(8) == x^4 + 1
    expect(R.cyclotomic_polynomial(8).toString()).toBe('x^4 + 1');
    // sage: ZZ['x'].cyclotomic_polynomial(12) == x^4 - x^2 + 1
    expect(R.cyclotomic_polynomial(12).toString()).toBe('x^4 + (-1)*x^2 + 1');
    // sage: ZZ['x'].cyclotomic_polynomial(1) == x - 1
    expect(R.cyclotomic_polynomial(1).toString()).toBe('x + -1');
  });

  test('n <= 0 raises ArithmeticError with Sage message', () => {
    // sage: R.cyclotomic_polynomial(0) -> ArithmeticError: n=0 must be positive
    expect(() => R.cyclotomic_polynomial(0)).toThrow(ArithmeticError);
    expect(() => R.cyclotomic_polynomial(0)).toThrow('n=0 must be positive');
    expect(() => R.cyclotomic_polynomial(-1)).toThrow('n=-1 must be positive');
  });

  test('non-integral n raises TypeError', () => {
    expect(() => R.cyclotomic_polynomial(1.5)).toThrow(TypeError);
  });
});
