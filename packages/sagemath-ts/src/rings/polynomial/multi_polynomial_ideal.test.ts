/**
 * Tests for multivariate polynomial ideals pinned to SageMath's own values.
 *
 * Every expected value below is taken from an upstream doctest; the file
 * reference is given next to each one. The pre-existing `groebner.test.ts`
 * only asserts that the *input generators* reduce to zero, which any superset
 * of the generators satisfies (audit M20) -- these tests pin the actual basis
 * and check the S-pair criterion.
 */

import { describe, expect, test } from 'bun:test';
import { FiniteFieldPrime } from '../finite_rings/finite_field_prime.js';
import type { Rational } from '../rational.js';
import { RationalField } from '../rational_field.js';
import {
  MPolynomialIdeal,
  groebner_basis,
  ideal,
  reduce,
  sPolynomial,
} from './multi_polynomial_ideal.js';
import { MPolynomialRing } from './multi_polynomial_ring.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';

const QQ = RationalField.getInstance();

/** Rational n/d as an element of QQ. */
const q = (n: number, d = 1): Rational => QQ.__call__(n).div(QQ.__call__(d));

/**
 * Katsura-3 over QQ in lex order, exactly as `sage.rings.ideal.Katsura(P, 3)`
 * prints it (`sage/rings/ideal.py:1854`):
 *   Ideal (x + 2*y + 2*z - 1, x^2 + 2*y^2 + 2*z^2 - x, 2*x*y + 2*y*z - y)
 */
function katsura3() {
  const R = new MPolynomialRing(QQ, ['a', 'b', 'c'], 'lex');
  const [a, b, c] = R.gens() as [
    ReturnType<typeof R.gen>,
    ReturnType<typeof R.gen>,
    ReturnType<typeof R.gen>,
  ];
  const two = QQ.__call__(2);
  const f1 = a.add(b.scalarMul(two)).add(c.scalarMul(two)).sub(R.one());
  const f2 = a.pow(2).add(b.pow(2).scalarMul(two)).add(c.pow(2).scalarMul(two)).sub(a);
  const f3 = a.mul(b).scalarMul(two).add(b.mul(c).scalarMul(two)).sub(b);
  return { R, a, b, c, gens: [f1, f2, f3] };
}

describe('groebner_basis: SageMath doctest values', () => {
  test('Katsura-3 lex basis matches SageMath exactly', () => {
    // sage/rings/polynomial/multi_polynomial_ideal.py:4613
    //   sage: P.<a,b,c> = PolynomialRing(QQbar, 3, order='lex')
    //   sage: I = sage.rings.ideal.Katsura(P,3)
    //   sage: I.groebner_basis()
    //   [a + (-60)*c^3 + 158/7*c^2 + 8/7*c - 1,
    //    b + 30*c^3 + (-79/7)*c^2 + 3/7*c,
    //    c^4 + (-10/21)*c^3 + 1/84*c^2 + 1/84*c]
    const { R, a, b, c, gens } = katsura3();
    const gb = groebner_basis(gens);

    const e1 = a
      .sub(c.pow(3).scalarMul(q(60)))
      .add(c.pow(2).scalarMul(q(158, 7)))
      .add(c.scalarMul(q(8, 7)))
      .sub(R.one());
    const e2 = b
      .add(c.pow(3).scalarMul(q(30)))
      .sub(c.pow(2).scalarMul(q(79, 7)))
      .add(c.scalarMul(q(3, 7)));
    const e3 = c
      .pow(4)
      .sub(c.pow(3).scalarMul(q(10, 21)))
      .add(c.pow(2).scalarMul(q(1, 84)))
      .add(c.scalarMul(q(1, 84)));

    expect(gb.length).toBe(3);
    expect(gb.map((g) => g.toString()).sort()).toEqual(
      [e1, e2, e3].map((g) => g.toString()).sort()
    );
  });

  test('Katsura-3 basis satisfies Buchberger criterion (every S-pair reduces to 0)', () => {
    const { gens } = katsura3();
    const gb = groebner_basis(gens);
    for (let i = 0; i < gb.length; i++) {
      for (let j = i + 1; j < gb.length; j++) {
        expect(reduce(sPolynomial(gb[i]!, gb[j]!), gb).isZero()).toBe(true);
      }
    }
  });

  test('exhausting maxIterations throws rather than returning a non-Gröbner basis', () => {
    // Audit M19: with a cap of 3 the old code returned a set for which an
    // ideal member did not reduce to zero.
    const { gens } = katsura3();
    expect(() => groebner_basis(gens, { maxIterations: 3 })).toThrow(/exhausted maxIterations/);
  });
});

describe('reduce over a non-field base ring', () => {
  /** Minimal ZZ used to reproduce audit H14. */
  class IntegerElement implements RingElement {
    constructor(readonly value: bigint) {}
    add(o: IntegerElement) {
      return new IntegerElement(this.value + o.value) as this;
    }
    sub(o: IntegerElement) {
      return new IntegerElement(this.value - o.value) as this;
    }
    mul(o: IntegerElement) {
      return new IntegerElement(this.value * o.value) as this;
    }
    div(o: IntegerElement) {
      return new IntegerElement(this.value / o.value) as this;
    }
    neg() {
      return new IntegerElement(-this.value) as this;
    }
    eq(o: IntegerElement | number) {
      return typeof o === 'number' ? this.value === BigInt(o) : this.value === o.value;
    }
    isZero() {
      return this.value === 0n;
    }
    toString() {
      return this.value.toString();
    }
  }

  const makeZZ = (withIsField: boolean): CoefficientRing<IntegerElement> => {
    const ring: CoefficientRing<IntegerElement> = {
      zero: () => new IntegerElement(0n),
      one: () => new IntegerElement(1n),
      __call__: (x: unknown) =>
        x instanceof IntegerElement ? x : new IntegerElement(BigInt(x as number | bigint)),
    };
    if (withIsField) {
      (ring as { is_field?: () => boolean }).is_field = () => false;
    }
    return ring;
  };

  test("reduce(x^2, [2x+1]) over ZZ raises SageMath's TypeError instead of hanging", () => {
    // sage/rings/polynomial/multi_polynomial_element.py:2488
    //   if not k.is_field(): raise TypeError("Can only reduce polynomials over fields.")
    const ZZ = makeZZ(true);
    const R = new MPolynomialRing(ZZ, ['x', 'y']);
    const x = R.gen(0);
    const g = x.scalarMul(new IntegerElement(2n)).add(R.one());

    const started = Date.now();
    expect(() => reduce(x.pow(2), [g])).toThrow('Can only reduce polynomials over fields.');
    expect(Date.now() - started).toBeLessThan(1000);
  });

  test('groebner_basis over ZZ raises the same TypeError', () => {
    const ZZ = makeZZ(true);
    const R = new MPolynomialRing(ZZ, ['x', 'y']);
    const x = R.gen(0);
    expect(() =>
      groebner_basis([x.pow(2), x.scalarMul(new IntegerElement(2n)).add(R.one())])
    ).toThrow('Can only reduce polynomials over fields.');
  });

  test('a base ring without is_field() still terminates via the exact-division guard', () => {
    const ZZ = makeZZ(false);
    const R = new MPolynomialRing(ZZ, ['x', 'y']);
    const x = R.gen(0);
    const started = Date.now();
    expect(() => reduce(x.pow(2), [x.scalarMul(new IntegerElement(2n)).add(R.one())])).toThrow(
      'Can only reduce polynomials over fields.'
    );
    expect(Date.now() - started).toBeLessThan(1000);
  });
});

describe('MPolynomialIdeal: zero ideal and dimension', () => {
  test('the zero ideal has Gröbner basis [0]', () => {
    // sage/rings/polynomial/multi_polynomial_ideal.py:4586
    //   sage: P.ideal([]).groebner_basis()
    //   [0]
    //   sage: P.ideal([0]).groebner_basis()
    //   [0]
    const R = new MPolynomialRing(QQ, ['x', 'y']);
    for (const I of [ideal(R, []), ideal(R, [R.zero()])]) {
      const gb = I.groebner_basis();
      expect(gb.length).toBe(1);
      expect(gb[0]!.isZero()).toBe(true);
    }
  });

  test('the empty generator list is accepted when the ring is given', () => {
    const R = new MPolynomialRing(QQ, ['x', 'y']);
    expect(() => ideal(R, [])).not.toThrow();
    // Without a ring the parent cannot be determined, so this still throws.
    expect(() => ideal([])).toThrow();
  });

  test('dimension of the zero ideal is the number of variables (an integer)', () => {
    // sage/rings/polynomial/multi_polynomial_ideal.py:1114
    //   sage: R.<x,y> = PolynomialRing(GF(2147483659^2), order='lex')
    //   sage: I = R.ideal(0); I.dimension()
    //   2
    const R = new MPolynomialRing(QQ, ['x', 'y'], 'lex');
    const dim = ideal(R, []).dimension();
    expect(dim).toBe(2);
    expect(Number.isFinite(dim)).toBe(true);
  });

  test('dimension of ideal(x^2 - y, x^3) in three variables is 1', () => {
    // sage/rings/polynomial/multi_polynomial_ideal.py:1080
    const F = new FiniteFieldPrime(32003n);
    const P = new MPolynomialRing(F, ['x', 'y', 'z'], 'degrevlex');
    const [x, y] = P.gens();
    expect(ideal(P, [x!.pow(2).sub(y!), x!.pow(3)]).dimension()).toBe(1);
  });

  test('SageMath toy-implementation dimension doctests', () => {
    // sage/rings/polynomial/multi_polynomial_ideal.py:1100-1117
    const R = new MPolynomialRing(QQ, ['x', 'y'], 'lex');
    const [x, y] = R.gens() as [ReturnType<typeof R.gen>, ReturnType<typeof R.gen>];

    // I = R.ideal([x*y, x*y + 1]) -> -1 (the total ring)
    expect(ideal(R, [x.mul(y), x.mul(y).add(R.one())]).dimension()).toBe(-1);
    // I = ideal([x*(x*y+1), y*(x*y+1)]) -> 1
    const u = x.mul(y).add(R.one());
    expect(ideal(R, [x.mul(u), y.mul(u)]).dimension()).toBe(1);
    // I = R.ideal([x^3*y, x*y^2]) -> 1
    expect(ideal(R, [x.pow(3).mul(y), x.mul(y.pow(2))]).dimension()).toBe(1);
  });

  test('principal ideals follow Theorem 3.5.1 of [Ger2008]', () => {
    // multi_polynomial_ideal.py:1133-1141
    const R = new MPolynomialRing(QQ, ['x', 'y'], 'lex');
    expect(ideal(R, [R.gen(0)]).dimension()).toBe(1); // nonzero non-unit -> n - 1
    expect(ideal(R, [R.one()]).dimension()).toBe(-1); // unit -> -1
    expect(ideal(R, [R.zero()]).dimension()).toBe(2); // zero -> n
  });

  test('dimension over a non-field base ring raises NotImplementedError', () => {
    // multi_polynomial_ideal.py:1131-1132
    const ZZ: CoefficientRing<never> & { is_field: () => boolean } = {
      zero: () => {
        throw new Error('unused');
      },
      one: () => {
        throw new Error('unused');
      },
      __call__: () => {
        throw new Error('unused');
      },
      is_field: () => false,
    };
    const R = new MPolynomialRing(ZZ as unknown as CoefficientRing<RingElement>, ['x']);
    expect(() => new MPolynomialIdeal(R, []).dimension()).toThrow('implemented only over fields');
  });
});
