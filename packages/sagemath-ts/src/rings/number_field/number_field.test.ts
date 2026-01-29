/**
 * Tests for number field functionality
 */

import { describe, expect, it } from 'bun:test';
import {
  NumberField,
  NumberFieldElement,
  QuadraticField,
  CyclotomicField,
  RationalPolynomial,
  NumberFieldConstructor,
} from './number_field.js';
import { Rational } from '../rational.js';

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

  it('should reduce to squarefree part', () => {
    // sqrt(12) = 2*sqrt(3), so Q(sqrt(12)) = Q(sqrt(3))
    const K = QuadraticField.create(12n);
    expect(K.d).toBe(3n);
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
      const O = K.ring_of_integers();
      const index = O.index_in_maximal_order();
      expect(index).toBe(1n);
    });

    it('should return 2 for equation order when d ≡ 1 (mod 4)', () => {
      // For Q(sqrt(5)), d=5 ≡ 1 (mod 4), so equation order has index 2 in maximal order
      // Equation order disc = 4*5 = 20, maximal order disc = 5
      // Index = sqrt(20/5) = 2
      const K = QuadraticField.create(5n);
      const O = K.ring_of_integers();
      const index = O.index_in_maximal_order();
      expect(index).toBe(2n);
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
      const O = K.ring_of_integers();
      const f = O.conductor();
      expect(f).toBe(2n);
    });
  });

  describe('is_maximal', () => {
    it('should return true when d ≢ 1 (mod 4)', () => {
      // For Q(sqrt(2)), equation order = maximal order
      const K = QuadraticField.create(2n);
      const O = K.ring_of_integers();
      expect(O.is_maximal()).toBe(true);
    });

    it('should return false when d ≡ 1 (mod 4)', () => {
      // For Q(sqrt(5)), d ≡ 1 (mod 4), so equation order ≠ maximal order
      const K = QuadraticField.create(5n);
      const O = K.ring_of_integers();
      expect(O.is_maximal()).toBe(false);
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
