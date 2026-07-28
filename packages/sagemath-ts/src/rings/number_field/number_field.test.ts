/**
 * Tests for number field functionality
 */

import { describe, expect, it } from 'bun:test';
import { Rational } from '../rational.js';
import {
  CyclotomicField,
  NumberField,
  NumberFieldConstructor,
  NumberFieldElement,
  QuadraticField,
  RationalPolynomial,
} from './number_field.js';
import { EquationOrder } from './order.js';

describe('RationalPolynomial', () => {
  describe('construction', () => {
    it('should create a polynomial from rationals', () => {
      const p = new RationalPolynomial([new Rational(1n), new Rational(2n), new Rational(3n)]);
      expect(p.degree()).toBe(2);
      expect(p.getCoeff(0).eq(new Rational(1n))).toBe(true);
      expect(p.getCoeff(1).eq(new Rational(2n))).toBe(true);
      expect(p.getCoeff(2).eq(new Rational(3n))).toBe(true);
    });

    it('should create from bigints', () => {
      const p = RationalPolynomial.fromBigInts([1n, 2n, 3n]);
      expect(p.degree()).toBe(2);
    });

    it('should remove trailing zeros', () => {
      const p = new RationalPolynomial([new Rational(1n), new Rational(2n), Rational.zero()]);
      expect(p.degree()).toBe(1);
    });

    it('should handle zero polynomial', () => {
      const p = RationalPolynomial.zero();
      expect(p.isZero()).toBe(true);
      expect(p.degree()).toBe(-1);
    });
  });

  describe('arithmetic', () => {
    it('should add polynomials', () => {
      const p1 = RationalPolynomial.fromBigInts([1n, 2n]);
      const p2 = RationalPolynomial.fromBigInts([3n, 4n, 5n]);
      const sum = p1.add(p2);
      expect(sum.getCoeff(0).eq(new Rational(4n))).toBe(true);
      expect(sum.getCoeff(1).eq(new Rational(6n))).toBe(true);
      expect(sum.getCoeff(2).eq(new Rational(5n))).toBe(true);
    });

    it('should subtract polynomials', () => {
      const p1 = RationalPolynomial.fromBigInts([5n, 6n]);
      const p2 = RationalPolynomial.fromBigInts([1n, 2n]);
      const diff = p1.sub(p2);
      expect(diff.getCoeff(0).eq(new Rational(4n))).toBe(true);
      expect(diff.getCoeff(1).eq(new Rational(4n))).toBe(true);
    });

    it('should multiply polynomials', () => {
      // (1 + x) * (1 - x) = 1 - x^2
      const p1 = RationalPolynomial.fromBigInts([1n, 1n]);
      const p2 = RationalPolynomial.fromBigInts([1n, -1n]);
      const prod = p1.mul(p2);
      expect(prod.getCoeff(0).eq(new Rational(1n))).toBe(true);
      expect(prod.getCoeff(1).isZero()).toBe(true);
      expect(prod.getCoeff(2).eq(new Rational(-1n))).toBe(true);
    });

    it('should divide with remainder', () => {
      // (x^2 + 2x + 1) / (x + 1) = (x + 1) with remainder 0
      const p1 = RationalPolynomial.fromBigInts([1n, 2n, 1n]);
      const p2 = RationalPolynomial.fromBigInts([1n, 1n]);
      const [q, r] = p1.divmod(p2);
      expect(q.getCoeff(0).eq(new Rational(1n))).toBe(true);
      expect(q.getCoeff(1).eq(new Rational(1n))).toBe(true);
      expect(r.isZero()).toBe(true);
    });

    it('should compute GCD', () => {
      const p1 = RationalPolynomial.fromBigInts([1n, 2n, 1n]); // (x+1)^2
      const p2 = RationalPolynomial.fromBigInts([1n, 1n]); // x+1
      const g = p1.gcd(p2);
      expect(g.degree()).toBe(1);
      expect(g.isMonic()).toBe(true);
    });
  });

  describe('evaluation', () => {
    it('should evaluate at a point', () => {
      // x^2 + 2x + 1 at x = 2 should be 9
      const p = RationalPolynomial.fromBigInts([1n, 2n, 1n]);
      const result = p.evaluate(new Rational(2n));
      expect(result.eq(new Rational(9n))).toBe(true);
    });
  });

  describe('derivative', () => {
    it('should compute derivative', () => {
      // d/dx (x^3 + 2x^2 + x) = 3x^2 + 4x + 1
      const p = RationalPolynomial.fromBigInts([0n, 1n, 2n, 1n]);
      const dp = p.derivative();
      expect(dp.getCoeff(0).eq(new Rational(1n))).toBe(true);
      expect(dp.getCoeff(1).eq(new Rational(4n))).toBe(true);
      expect(dp.getCoeff(2).eq(new Rational(3n))).toBe(true);
    });
  });

  describe('discriminant', () => {
    it('should compute discriminant of quadratic', () => {
      // x^2 - 2 has discriminant 8
      const p = RationalPolynomial.fromBigInts([-2n, 0n, 1n]);
      const disc = p.discriminant();
      expect(disc.eq(new Rational(8n))).toBe(true);
    });
  });
});

describe('NumberField', () => {
  describe('construction', () => {
    it('should create a quadratic field', () => {
      // x^2 - 2
      const poly = RationalPolynomial.fromBigInts([-2n, 0n, 1n]);
      const K = new NumberField(poly, 'a');
      expect(K.degree()).toBe(2);
    });

    it('should create a cubic field', () => {
      // x^3 - 2
      const poly = RationalPolynomial.fromBigInts([-2n, 0n, 0n, 1n]);
      const K = new NumberField(poly, 'a');
      expect(K.degree()).toBe(3);
    });

    it('should store monic polynomial', () => {
      // 2x^2 - 4 should become x^2 - 2
      const poly = RationalPolynomial.fromBigInts([-4n, 0n, 2n]);
      const K = new NumberField(poly, 'a');
      expect(K.polynomial().isMonic()).toBe(true);
    });
  });

  describe('basic operations', () => {
    it('should return generator', () => {
      const poly = RationalPolynomial.fromBigInts([-2n, 0n, 1n]);
      const K = new NumberField(poly, 'a');
      const a = K.gen();
      expect(a.__getitem__(0).isZero()).toBe(true);
      expect(a.__getitem__(1).eq(Rational.one())).toBe(true);
    });

    it('should return zero and one', () => {
      const poly = RationalPolynomial.fromBigInts([-2n, 0n, 1n]);
      const K = new NumberField(poly, 'a');
      expect(K.zero().is_zero()).toBe(true);
      expect(K.one().is_one()).toBe(true);
    });
  });

  describe('signature', () => {
    it('should compute signature of real quadratic field', () => {
      // x^2 - 2 (real quadratic)
      const poly = RationalPolynomial.fromBigInts([-2n, 0n, 1n]);
      const K = new NumberField(poly, 'a');
      const [r1, r2] = K.signature();
      expect(r1).toBe(2);
      expect(r2).toBe(0);
    });

    it('should compute signature of imaginary quadratic field', () => {
      // x^2 + 1 (imaginary quadratic)
      const poly = RationalPolynomial.fromBigInts([1n, 0n, 1n]);
      const K = new NumberField(poly, 'a');
      const [r1, r2] = K.signature();
      expect(r1).toBe(0);
      expect(r2).toBe(1);
    });
  });

  describe('discriminant', () => {
    it('should compute discriminant of quadratic field', () => {
      // x^2 - 2 has polynomial discriminant 8
      const poly = RationalPolynomial.fromBigInts([-2n, 0n, 1n]);
      const K = new NumberField(poly, 'a');
      expect(K.discriminant()).toBe(8n);
    });
  });
});

describe('NumberFieldElement', () => {
  const poly = RationalPolynomial.fromBigInts([-2n, 0n, 1n]); // x^2 - 2
  const K = new NumberField(poly, 'a');

  describe('arithmetic', () => {
    it('should add elements', () => {
      const a = K.gen();
      const one = K.one();
      const sum = a.add(one); // a + 1
      expect(sum.__getitem__(0).eq(Rational.one())).toBe(true);
      expect(sum.__getitem__(1).eq(Rational.one())).toBe(true);
    });

    it('should subtract elements', () => {
      const a = K.gen();
      const one = K.one();
      const diff = a.sub(one); // a - 1
      expect(diff.__getitem__(0).eq(new Rational(-1n))).toBe(true);
      expect(diff.__getitem__(1).eq(Rational.one())).toBe(true);
    });

    it('should multiply elements', () => {
      const a = K.gen();
      const a2 = a.mul(a); // a^2 = 2
      expect(a2.__getitem__(0).eq(new Rational(2n))).toBe(true);
      expect(a2.__getitem__(1).isZero()).toBe(true);
    });

    it('should compute powers', () => {
      const a = K.gen();
      const a3 = a.pow(3n); // a^3 = 2a
      expect(a3.__getitem__(0).isZero()).toBe(true);
      expect(a3.__getitem__(1).eq(new Rational(2n))).toBe(true);
    });

    it('should invert elements', () => {
      const a = K.gen();
      const aInv = a.inv(); // 1/a = a/2
      const product = a.mul(aInv);
      expect(product.is_one()).toBe(true);
    });

    it('should divide elements', () => {
      const a = K.gen();
      const two = K.__call__(2n);
      const result = a.div(two); // a/2
      expect(result.__getitem__(0).isZero()).toBe(true);
      expect(result.__getitem__(1).eq(new Rational(1n, 2n))).toBe(true);
    });
  });

  describe('norm and trace', () => {
    it('should compute trace', () => {
      const a = K.gen(); // sqrt(2)
      // trace(sqrt(2)) = sqrt(2) + (-sqrt(2)) = 0
      const tr = a.trace();
      expect(tr.isZero()).toBe(true);
    });

    it('should compute trace of sum', () => {
      const a = K.gen();
      const elem = a.add(K.one()); // sqrt(2) + 1
      // trace(sqrt(2) + 1) = (sqrt(2)+1) + (-sqrt(2)+1) = 2
      const tr = elem.trace();
      expect(tr.eq(new Rational(2n))).toBe(true);
    });

    it('should compute norm', () => {
      const a = K.gen(); // sqrt(2)
      // norm(sqrt(2)) = sqrt(2) * (-sqrt(2)) = -2
      const nm = a.norm();
      expect(nm.eq(new Rational(-2n))).toBe(true);
    });

    it('should compute norm of sum', () => {
      const a = K.gen();
      const elem = a.add(K.one()); // sqrt(2) + 1
      // norm(sqrt(2) + 1) = (sqrt(2)+1)(-sqrt(2)+1) = 1-2 = -1
      const nm = elem.norm();
      expect(nm.eq(new Rational(-1n))).toBe(true);
    });
  });

  describe('characteristic polynomial', () => {
    it('should compute charpoly of generator', () => {
      const a = K.gen();
      const cp = a.charpoly();
      // charpoly of sqrt(2) is x^2 - 2
      expect(cp.degree()).toBe(2);
      expect(cp.getCoeff(0).eq(new Rational(-2n))).toBe(true);
      expect(cp.getCoeff(1).isZero()).toBe(true);
      expect(cp.getCoeff(2).eq(Rational.one())).toBe(true);
    });

    it('should compute charpoly of constant', () => {
      const two = K.__call__(2n);
      const cp = two.charpoly();
      // charpoly of 2 is (x-2)^2 = x^2 - 4x + 4
      expect(cp.degree()).toBe(2);
      expect(cp.getCoeff(0).eq(new Rational(4n))).toBe(true);
      expect(cp.getCoeff(1).eq(new Rational(-4n))).toBe(true);
    });
  });

  describe('minimal polynomial', () => {
    it('should compute minpoly of generator', () => {
      const a = K.gen();
      const mp = a.minpoly();
      // minpoly of sqrt(2) is x^2 - 2
      expect(mp.degree()).toBe(2);
    });

    it('should compute minpoly of rational', () => {
      const two = K.__call__(2n);
      const mp = two.minpoly();
      // minpoly of 2 is x - 2
      expect(mp.degree()).toBe(1);
      expect(mp.getCoeff(0).eq(new Rational(-2n))).toBe(true);
    });
  });

  describe('is_integral', () => {
    it('should recognize integral elements', () => {
      const a = K.gen();
      expect(a.is_integral()).toBe(true);

      const aPlusOne = a.add(K.one());
      expect(aPlusOne.is_integral()).toBe(true);
    });

    it('should recognize non-integral elements', () => {
      const a = K.gen();
      const half = K.__call__(new Rational(1n, 2n));
      expect(half.is_integral()).toBe(false);
    });
  });

  describe('is_unit', () => {
    it('should recognize units', () => {
      const one = K.one();
      expect(one.is_unit()).toBe(true);

      const minusOne = K.__call__(-1n);
      expect(minusOne.is_unit()).toBe(true);
    });

    it('should recognize sqrt(2)+1 as unit (norm = -1)', () => {
      const a = K.gen();
      const unit = a.add(K.one()); // sqrt(2) + 1
      expect(unit.is_unit()).toBe(true);
    });
  });
});

describe('QuadraticField', () => {
  it('should create Q(sqrt(2))', () => {
    const K = QuadraticField.create(2n);
    expect(K.degree()).toBe(2);
    expect(K.d).toBe(2n);
  });

  it('should create Q(sqrt(-1))', () => {
    const K = QuadraticField.create(-1n);
    expect(K.degree()).toBe(2);
    expect(K.d).toBe(-1n);
  });

  it('should compute discriminant correctly', () => {
    // d = 5 (1 mod 4): disc = 5
    const K1 = QuadraticField.create(5n);
    expect(K1.discriminant()).toBe(5n);

    // d = 2 (not 1 mod 4): disc = 8
    const K2 = QuadraticField.create(2n);
    expect(K2.discriminant()).toBe(8n);

    // d = -1 (not 1 mod 4): disc = -4
    const K3 = QuadraticField.create(-1n);
    expect(K3.discriminant()).toBe(-4n);
  });

  it('should identify real vs imaginary', () => {
    const realK = QuadraticField.create(5n);
    expect(realK.is_totally_real()).toBe(true);

    const imagK = QuadraticField.create(-5n);
    expect(imagK.is_totally_real()).toBe(false);
    expect(imagK.is_totally_imaginary()).toBe(true);
  });

  it('should use x^2 - D verbatim, not the squarefree part', () => {
    // Sage: QuadraticField(8).defining_polynomial() is x^2 - 8 and
    // QuadraticField(8).gen()^2 == 8.
    const K = QuadraticField.create(8n);
    expect(K.defining_polynomial().toString()).toBe('x^2 - 8');
    expect(K.gen().pow(2n).eq(K.__call__(8n))).toBe(true);
    expect(K.D).toBe(8n);
    // The squarefree part is exposed separately: Q(sqrt(12)) = Q(sqrt(3)).
    expect(QuadraticField.create(12n).d).toBe(3n);
    // disc(Q(sqrt(8))) = disc(Q(sqrt(2))) = 8
    expect(K.discriminant()).toBe(8n);
  });

  it('should reject perfect squares', () => {
    // Sage: QuadraticField(9) raises ValueError("D must not be a perfect square.")
    expect(() => QuadraticField.create(9n)).toThrow('D must not be a perfect square.');
    expect(() => QuadraticField.create(1n)).toThrow('D must not be a perfect square.');
    expect(() => QuadraticField.create(0n)).toThrow('D must not be a perfect square.');
  });
});

describe('CyclotomicField', () => {
  it('should create Q(zeta_3)', () => {
    const K = CyclotomicField.create(3n);
    // phi(3) = 2
    expect(K.degree()).toBe(2);
    expect(K.n).toBe(3n);
  });

  it('should create Q(zeta_4)', () => {
    const K = CyclotomicField.create(4n);
    // phi(4) = 2
    expect(K.degree()).toBe(2);
  });

  it('should create Q(zeta_5)', () => {
    const K = CyclotomicField.create(5n);
    // phi(5) = 4
    expect(K.degree()).toBe(4);
  });

  it('should create Q(zeta_6)', () => {
    const K = CyclotomicField.create(6n);
    // phi(6) = 2
    expect(K.degree()).toBe(2);
  });

  it('should return primitive roots', () => {
    const K = CyclotomicField.create(4n);
    const zeta = K.zeta();
    // zeta_4 should satisfy zeta^4 = 1 and zeta^2 = -1
    const zeta2 = zeta.pow(2n);
    const minusOne = K.__call__(-1n);
    expect(zeta2.eq(minusOne)).toBe(true);
  });
});

describe('NumberFieldConstructor', () => {
  it('should create field from bigint array', () => {
    const K = NumberFieldConstructor([-2n, 0n, 1n], 'a');
    expect(K.degree()).toBe(2);
  });

  it('should create field from RationalPolynomial', () => {
    const poly = RationalPolynomial.fromBigInts([-2n, 0n, 1n]);
    const K = NumberFieldConstructor(poly, 'a');
    expect(K.degree()).toBe(2);
  });
});

describe('Order', () => {
  describe('discriminant', () => {
    it('should compute discriminant for Q(sqrt(2)) equation order', () => {
      const K = QuadraticField.create(2n);
      const O = K.ring_of_integers();
      const disc = O.discriminant();
      // For Q(sqrt(2)), d=2 is 2 mod 4, so disc = 4*2 = 8
      // The equation order Z[sqrt(2)] has discriminant 8, which equals the field discriminant
      expect(disc).toBe(8n);
    });

    it('should compute discriminant for Q(sqrt(-1)) equation order', () => {
      const K = QuadraticField.create(-1n);
      const O = K.ring_of_integers();
      const disc = O.discriminant();
      // For Q(i), d=-1 is 3 mod 4, so disc = 4*(-1) = -4
      expect(disc).toBe(-4n);
    });
  });

  describe('index_in_maximal_order', () => {
    it('should return 1 for equation order when it equals maximal order', () => {
      // For Q(sqrt(2)), d=2 is not 1 mod 4, so equation order = maximal order
      const K = QuadraticField.create(2n);
      const O = EquationOrder(K);
      const index = O.index_in_maximal_order();
      expect(index).toBe(1n);
    });

    it('should return 2 for equation order when d ≡ 1 (mod 4)', () => {
      // For Q(sqrt(5)), d=5 ≡ 1 (mod 4), so Z[sqrt(5)] has index 2 in O_K.
      // Equation order disc = 4*5 = 20, maximal order disc = 5, index = 2.
      // NOTE: K.ring_of_integers() is the *maximal* order (Sage: index 1), so
      // the equation order has to be built explicitly.
      const K = QuadraticField.create(5n);
      const O = EquationOrder(K);
      const index = O.index_in_maximal_order();
      expect(index).toBe(2n);
      expect(K.ring_of_integers().index_in_maximal_order()).toBe(1n);
    });
  });

  describe('conductor', () => {
    it('should return 1 for maximal order', () => {
      // Use d=2 where equation order = maximal order
      const K = QuadraticField.create(2n);
      const O = K.ring_of_integers();
      const f = O.conductor();
      expect(f).toBe(1n);
    });

    it('should return 2 for equation order when d ≡ 1 (mod 4)', () => {
      const K = QuadraticField.create(5n);
      expect(EquationOrder(K).conductor()).toBe(2n);
      // Sage: QuadraticField(5).ring_of_integers().conductor() == 1
      expect(K.ring_of_integers().conductor()).toBe(1n);
    });
  });

  describe('is_maximal', () => {
    it('should return true when d ≢ 1 (mod 4)', () => {
      // For Q(sqrt(2)), equation order = maximal order
      const K = QuadraticField.create(2n);
      const O = K.ring_of_integers();
      expect(O.is_maximal()).toBe(true);
    });

    it('should return false for the equation order when d ≡ 1 (mod 4)', () => {
      // For Q(sqrt(5)), d ≡ 1 (mod 4), so Z[sqrt(5)] != O_K.
      // Sage: QuadraticField(5).ring_of_integers().is_maximal() is True and its
      // basis is [1/2*a + 1/2, a]; only the equation order is non-maximal.
      const K = QuadraticField.create(5n);
      expect(EquationOrder(K).is_maximal()).toBe(false);
      expect(K.ring_of_integers().is_maximal()).toBe(true);
    });
  });

  describe('class_number', () => {
    it('should be 1 for Q(sqrt(2))', () => {
      const K = QuadraticField.create(2n);
      const h = K.class_number();
      expect(h).toBe(1n);
    });

    it('should be 1 for Q(i)', () => {
      const K = QuadraticField.create(-1n);
      const h = K.class_number();
      expect(h).toBe(1n);
    });

    it('should be 1 for Q(sqrt(-3))', () => {
      const K = QuadraticField.create(-3n);
      const h = K.class_number();
      expect(h).toBe(1n);
    });

    it('should be 3 for Q(sqrt(-23))', () => {
      const K = QuadraticField.create(-23n);
      const h = K.class_number();
      expect(h).toBe(3n);
    });

    it('should be 2 for Q(sqrt(-15))', () => {
      const K = QuadraticField.create(-15n);
      const h = K.class_number();
      expect(h).toBe(2n);
    });
  });
});

describe('UnitGroup', () => {
  describe('torsion_order', () => {
    it('should be 4 for Q(i)', () => {
      const K = QuadraticField.create(-1n);
      const U = K.unit_group();
      expect(U.torsion_order()).toBe(4n);
    });

    it('should be 6 for Q(sqrt(-3))', () => {
      const K = QuadraticField.create(-3n);
      const U = K.unit_group();
      expect(U.torsion_order()).toBe(6n);
    });

    it('should be 2 for other imaginary quadratic fields', () => {
      const K = QuadraticField.create(-7n);
      const U = K.unit_group();
      expect(U.torsion_order()).toBe(2n);
    });

    it('should be 2 for real quadratic fields', () => {
      const K = QuadraticField.create(5n);
      const U = K.unit_group();
      expect(U.torsion_order()).toBe(2n);
    });
  });

  describe('rank', () => {
    it('should be 0 for imaginary quadratic fields', () => {
      const K = QuadraticField.create(-5n);
      const U = K.unit_group();
      expect(U.rank()).toBe(0);
    });

    it('should be 1 for real quadratic fields', () => {
      const K = QuadraticField.create(5n);
      const U = K.unit_group();
      expect(U.rank()).toBe(1);
    });
  });

  describe('order', () => {
    it('should be finite for imaginary quadratic fields', () => {
      const K = QuadraticField.create(-7n);
      const U = K.unit_group();
      expect(U.order()).toBe(2n);
    });

    it('should be infinite for real quadratic fields', () => {
      const K = QuadraticField.create(5n);
      const U = K.unit_group();
      expect(U.order()).toBe('infinity');
    });
  });

  describe('roots_of_unity', () => {
    it('should return 4 roots for Q(i)', () => {
      const K = QuadraticField.create(-1n);
      const U = K.unit_group();
      const roots = U.roots_of_unity();
      expect(roots.length).toBe(4);
    });

    it('should return 6 roots for Q(sqrt(-3))', () => {
      const K = QuadraticField.create(-3n);
      const U = K.unit_group();
      const roots = U.roots_of_unity();
      expect(roots.length).toBe(6);
    });
  });
});

// ---------------------------------------------------------------------------
// Regressions for the 2026-07 audit (H16, H17, H18, H22, H23, H24, M27, M31, L36)
// ---------------------------------------------------------------------------

describe('field discriminant (audit H16)', () => {
  it('returns disc(O_K), not the polynomial discriminant', () => {
    // Sage: K.<t> = NumberField(x^3 + x^2 - 2*x + 8); K.disc() == -503
    //       K.disc([1, t, t^2]) == -2012  (the polynomial discriminant)
    const K = NumberFieldConstructor([8n, -2n, 1n, 1n], 'a');
    expect(K.discriminant()).toBe(-503n);
    expect(K.defining_polynomial().discriminant().numerator).toBe(-2012n);
  });

  it('handles non-integral defining polynomials', () => {
    // Sage: NumberField(x^2 - 1/2, 'a').discriminant() == 8
    const K = new NumberField(
      new RationalPolynomial([new Rational(-1n, 2n), Rational.zero(), Rational.one()]),
      'a'
    );
    expect(K.discriminant()).toBe(8n);
  });

  it('matches Sage on quadratic and cyclotomic fields', () => {
    expect(NumberFieldConstructor([-5n, 0n, 1n], 'a').discriminant()).toBe(5n);
    expect(QuadraticField.create(2n).discriminant()).toBe(8n);
    expect(QuadraticField.create(5n).discriminant()).toBe(5n);
    expect(QuadraticField.create(-1n).discriminant()).toBe(-4n);
    expect(QuadraticField.create(-30n).discriminant()).toBe(-120n);
    expect(CyclotomicField.create(5n).discriminant()).toBe(125n);
    expect(CyclotomicField.create(7n).discriminant()).toBe(-16807n);
  });

  it('computes Sage integral bases', () => {
    // Sage: NumberField(x^3 + x^2 - 2*x + 8).integral_basis()
    //       == [1, 1/2*a^2 + 1/2*a, a^2]
    const K = NumberFieldConstructor([8n, -2n, 1n, 1n], 'a');
    expect(K.integral_basis().map((b) => b.toString())).toEqual(['1', '1/2*a^2 + 1/2*a', 'a^2']);
    // Sage: NumberField(x^5 + 10*x + 1).integral_basis() == [1, a, a^2, a^3, a^4]
    const K5 = NumberFieldConstructor([1n, 10n, 0n, 0n, 0n, 1n], 'a');
    expect(K5.integral_basis().map((b) => b.toString())).toEqual(['1', 'a', 'a^2', 'a^3', 'a^4']);
  });
});

describe('class number and class group (audit H17, H18)', () => {
  it('keys the class number on the field discriminant', () => {
    // All three fields have class number 1 (they are Heegner discriminants);
    // the port used to look the *polynomial* discriminant up in a table.
    expect(NumberFieldConstructor([11n, 0n, 1n], 'a').class_number()).toBe(1n);
    expect(NumberFieldConstructor([9n, 0n, 1n], 'a').class_number()).toBe(1n);
    expect(NumberFieldConstructor([27n, 0n, 1n], 'a').class_number()).toBe(1n);
    expect(NumberFieldConstructor([23n, 0n, 1n], 'a').class_number()).toBe(3n);
  });

  it('returns the true abelian structure', () => {
    // disc = -120 = 8 * 5 * (-3): genus theory gives 2-rank 2, so Cl = C2 x C2
    const G = QuadraticField.create(-30n).class_group();
    expect(G.invariants()).toEqual([2n, 2n]);
    expect(G.order()).toBe(4n);
    expect(G.is_cyclic()).toBe(false);
  });

  it('matches the Sage doctest for x^2 + 20072 (C38 x C2)', () => {
    const G = NumberFieldConstructor([20072n, 0n, 1n], 'a').class_group();
    expect(G.order()).toBe(76n);
    expect(G.invariants()).toEqual([38n, 2n]);
    expect(G.is_cyclic()).toBe(false);
  });

  it('matches the Sage doctest for x^2 + 23 (C3)', () => {
    const G = NumberFieldConstructor([23n, 0n, 1n], 'a').class_group();
    expect(G.order()).toBe(3n);
    expect(G.invariants()).toEqual([3n]);
    expect(G.is_cyclic()).toBe(true);
  });
});

describe('irreducibility of the defining polynomial (audit H22)', () => {
  it('rejects reducible defining polynomials with Sage’s message', () => {
    expect(() => NumberFieldConstructor([-1n, 0n, 1n], 'a')).toThrow(
      'defining polynomial (x^2 - 1) must be irreducible'
    );
    expect(() => NumberFieldConstructor([-4n, 0n, 1n], 'a')).toThrow(
      'defining polynomial (x^2 - 4) must be irreducible'
    );
    expect(() => NumberFieldConstructor([1n, 0n, 2n, 0n, 1n], 'a')).toThrow('must be irreducible');
  });

  it('detects reducible polynomials without rational roots', () => {
    // (x^2 + 1)^2 and (x^2 + 1)(x^2 + 2) have no rational roots
    expect(RationalPolynomial.fromBigInts([1n, 0n, 2n, 0n, 1n]).isIrreducible()).toBe(false);
    expect(RationalPolynomial.fromBigInts([2n, 0n, 3n, 0n, 1n]).isIrreducible()).toBe(false);
    expect(RationalPolynomial.fromBigInts([4n, 0n, 0n, 0n, 1n]).isIrreducible()).toBe(false);
    expect(RationalPolynomial.fromBigInts([1n, 0n, 0n, 0n, 1n]).isIrreducible()).toBe(true);
    expect(RationalPolynomial.fromBigInts([1n, 1n, 1n, 1n, 1n, 1n, 1n]).isIrreducible()).toBe(true);
  });
});

describe('automorphisms (audit H24)', () => {
  it('finds all three automorphisms of the cyclic cubic x^3 - 3x + 1', () => {
    // Sage: K.<a> = NumberField(x^3 - 3*x + 1); K.automorphisms() has 3 elements
    //       a |--> a, a |--> a^2 - 2, a |--> -a^2 - a + 2
    const K = NumberFieldConstructor([1n, -3n, 0n, 1n], 'a');
    const auts = K.automorphisms();
    expect(auts.length).toBe(3);
    expect(auts.map((s) => s.im_gens()[0]!.toString()).sort()).toEqual([
      '-a^2 - a + 2',
      'a',
      'a^2 - 2',
    ]);
    expect(K.is_galois()).toBe(true);
    expect(K.galois_closure()).toBe(K);
  });

  it('returns only the identity for the non-Galois x^3 - 2', () => {
    const K = NumberFieldConstructor([-2n, 0n, 0n, 1n], 'a');
    expect(K.automorphisms().length).toBe(1);
    expect(K.is_galois()).toBe(false);
  });

  it('matches Sage for x^6 - x^4 - 2x^2 + 1', () => {
    // Sage doctest: len(NumberField(x^6 - x^4 - 2*x^2 + 1, 'a').automorphisms()) == 2
    const K = NumberFieldConstructor([1n, 0n, -2n, 0n, -1n, 0n, 1n], 'a');
    expect(K.automorphisms().length).toBe(2);
  });

  it('sends the generator to genuine roots of the defining polynomial', () => {
    const K = NumberFieldConstructor([1n, -3n, 0n, 1n], 'a');
    for (const s of K.automorphisms()) {
      const beta = s.im_gens()[0]!;
      let acc = K.zero();
      let pw = K.one();
      for (let i = 0; i <= 3; i++) {
        acc = acc.add(pw.scalarMul(K.defining_polynomial().getCoeff(i)));
        pw = pw.mul(beta);
      }
      expect(acc.is_zero()).toBe(true);
    }
  });

  it('models cyclotomic automorphisms as zeta -> zeta^k', () => {
    const K = CyclotomicField.create(8n);
    const auts = K.automorphisms();
    expect(auts.length).toBe(4);
    for (const s of auts) {
      const img = s.im_gens()[0]!;
      expect(img.pow(8n).is_one()).toBe(true);
      expect(img.pow(4n).is_one()).toBe(false);
    }
  });
});

describe('is_unit on a field vs an order (audit M27)', () => {
  it('treats every nonzero element of the field as a unit', () => {
    // Sage: K.<a> = NumberField(x^2 - x - 1); K(13).is_unit() is True,
    //       OK(13).is_unit() is False, OK(a).is_unit() is True
    const K = NumberFieldConstructor([-1n, -1n, 1n], 'a');
    expect(K.__call__(13n).is_unit()).toBe(true);
    expect(K.zero().is_unit()).toBe(false);
    expect(K.__call__(13n).is_integral_unit()).toBe(false);
    expect(K.gen().is_integral_unit()).toBe(true);
  });
});

describe('ring_of_integers is the maximal order (audit M31)', () => {
  it('returns O_K, not Z[alpha]', () => {
    const K = QuadraticField.create(5n);
    const O = K.ring_of_integers();
    // Sage: QuadraticField(5).maximal_order().basis() == [1/2*a + 1/2, a]
    expect(O.basis().map((b) => b.toString())).toEqual(['1/2*a + 1/2', 'a']);
    expect(O.is_maximal()).toBe(true);
    expect(O.discriminant()).toBe(5n);
    expect(O.index_in_maximal_order()).toBe(1n);
  });

  it('still exposes the (non-maximal) equation order', () => {
    const K = QuadraticField.create(5n);
    expect(EquationOrder(K).discriminant()).toBe(20n);
    expect(EquationOrder(K).is_maximal()).toBe(false);
  });
});

describe('prime decomposition (audit L36)', () => {
  it('primes_above returns the list and prime_above a single ideal', () => {
    const K = QuadraticField.create(-1n);
    expect(K.primes_above(5n).length).toBe(2); // split
    expect(K.primes_above(3n).length).toBe(1); // inert
    expect(K.primes_above(2n).length).toBe(1); // ramified
    expect(K.decomposition(2n)[0]![1]).toBe(2n);
    const P = K.prime_above(5n);
    expect(P.norm().eq(new Rational(5n))).toBe(true);
    expect(K.prime_above(5n, { degree: 1n }).residue_class_degree()).toBe(1n);
  });

  it('produces a factorisation of p*O_K', () => {
    for (const K of [
      QuadraticField.create(-1n),
      QuadraticField.create(5n),
      NumberFieldConstructor([-2n, 0n, 0n, 1n], 'a'),
    ]) {
      for (const p of [2n, 3n, 5n, 7n, 31n]) {
        const dec = K.decomposition(p);
        const prod = dec.reduce((acc, [Q, e]) => acc.mul(Q.pow(e)), K.ideal(1n));
        expect(prod.eq(K.ideal(p))).toBe(true);
        // sum e_i f_i = n
        const efSum = dec.reduce((acc, [Q, e]) => acc + e * Q.residue_class_degree(), 0n);
        expect(efSum).toBe(BigInt(K.degree()));
      }
    }
  });

  it('handles primes dividing the index of the equation order', () => {
    // Z[sqrt(5)] has index 2 in O_K, so Dedekind-Kummer needs another generator.
    const K = QuadraticField.create(5n);
    const dec = K.decomposition(2n);
    expect(dec.length).toBe(1);
    expect(dec[0]![1]).toBe(1n); // unramified
    expect(dec[0]![0].residue_class_degree()).toBe(2n); // inert
  });
});
