/**
 * Tests for Order functionality
 */

import { describe, expect, it } from 'bun:test';
import { Order, AbsoluteOrder, RelativeOrder, EquationOrder } from './order.js';
import { NumberField, QuadraticField, CyclotomicField, RationalPolynomial } from './number_field.js';
import { Rational } from '../rational.js';

describe('Order', () => {
  describe('construction', () => {
    it('should create order from quadratic field', () => {
      const K = QuadraticField.create(2n);
      const O = new Order(K);
      expect(O.rank()).toBe(2);
    });

    it('should create order from cyclotomic field', () => {
      const K = CyclotomicField.create(3n);
      const O = new Order(K);
      expect(O.rank()).toBe(2);
    });

    it('should create order from cubic field', () => {
      const poly = RationalPolynomial.fromBigInts([-2n, 0n, 0n, 1n]); // x^3 - 2
      const K = new NumberField(poly, 'a');
      const O = new Order(K);
      expect(O.rank()).toBe(3);
    });
  });

  describe('basis', () => {
    it('should return power basis by default', () => {
      const K = QuadraticField.create(2n);
      const O = new Order(K);
      const basis = O.basis();
      expect(basis.length).toBe(2);

      // First element should be 1
      expect(basis[0]!.is_one()).toBe(true);

      // Second element should be the generator
      const gen = K.gen();
      expect(basis[1]!.eq(gen)).toBe(true);
    });

    it('should return power basis method', () => {
      const K = QuadraticField.create(5n);
      const O = new Order(K);
      const pb = O.power_basis();
      expect(pb.length).toBe(2);
    });
  });

  describe('discriminant', () => {
    it('should compute discriminant for Q(sqrt(2))', () => {
      const K = QuadraticField.create(2n);
      const O = new Order(K);
      const disc = O.discriminant();
      // Discriminant of Z[sqrt(2)] is 8
      expect(disc).toBe(8n);
    });

    it('should compute discriminant for Q(sqrt(5))', () => {
      const K = QuadraticField.create(5n);
      const O = new Order(K);
      const disc = O.discriminant();
      // Discriminant of Z[sqrt(5)] (power basis) is 20
      // But field discriminant is 5 (since 5 = 1 mod 4)
      // The equation order discriminant is n^2 * field_disc for some index n
      expect(disc).toBe(20n);
    });

    it('should compute discriminant for Q(sqrt(-1))', () => {
      const K = QuadraticField.create(-1n);
      const O = new Order(K);
      const disc = O.discriminant();
      // Discriminant of Z[i] is -4
      expect(disc).toBe(-4n);
    });

    it('should compute discriminant for Q(sqrt(-3))', () => {
      const K = QuadraticField.create(-3n);
      const O = new Order(K);
      const disc = O.discriminant();
      // Discriminant of Z[sqrt(-3)] is -12
      // Field discriminant is -3 (since -3 = 1 mod 4)
      expect(disc).toBe(-12n);
    });
  });

  describe('is_maximal', () => {
    it('should identify maximal order for Q(sqrt(2))', () => {
      const K = QuadraticField.create(2n);
      const O = new Order(K);
      // Z[sqrt(2)] is maximal because disc(Z[sqrt(2)]) = 8 = disc(Q(sqrt(2)))
      expect(O.is_maximal()).toBe(true);
    });

    it('should identify maximal order for Q(sqrt(-1))', () => {
      const K = QuadraticField.create(-1n);
      const O = new Order(K);
      // Z[i] is maximal
      expect(O.is_maximal()).toBe(true);
    });
  });

  describe('contains', () => {
    it('should check if element is in order', () => {
      const K = QuadraticField.create(2n);
      const O = new Order(K);
      const a = K.gen();

      // sqrt(2) is in Z[sqrt(2)]
      expect(O.contains(a)).toBe(true);

      // 1 is in Z[sqrt(2)]
      expect(O.contains(K.one())).toBe(true);

      // 2 + 3*sqrt(2) is in Z[sqrt(2)]
      const elem = K.__call__(2n).add(a.scalarMul(new Rational(3n)));
      expect(O.contains(elem)).toBe(true);
    });

    it('should reject non-integral elements', () => {
      const K = QuadraticField.create(2n);
      const O = new Order(K);

      // 1/2 is not in Z[sqrt(2)]
      const half = K.__call__(new Rational(1n, 2n));
      expect(O.contains(half)).toBe(false);

      // sqrt(2)/2 is not in Z[sqrt(2)]
      const a = K.gen();
      const halfA = a.scalarMul(new Rational(1n, 2n));
      expect(O.contains(halfA)).toBe(false);
    });
  });

  describe('coordinates', () => {
    it('should compute coordinates of elements', () => {
      const K = QuadraticField.create(2n);
      const O = new Order(K);
      const a = K.gen();

      // Coordinates of 1 + 2*sqrt(2) should be [1, 2]
      const elem = K.one().add(a.scalarMul(new Rational(2n)));
      const coords = O.coordinates(elem);

      expect(coords).not.toBeNull();
      expect(coords![0]!.eq(new Rational(1n))).toBe(true);
      expect(coords![1]!.eq(new Rational(2n))).toBe(true);
    });

    it('should return fractional coordinates for non-integral elements', () => {
      const K = QuadraticField.create(2n);
      const O = new Order(K);

      // 1/2 has coordinates [1/2, 0]
      const half = K.__call__(new Rational(1n, 2n));
      const coords = O.coordinates(half);

      expect(coords).not.toBeNull();
      expect(coords![0]!.eq(new Rational(1n, 2n))).toBe(true);
      expect(coords![1]!.isZero()).toBe(true);
    });
  });

  describe('random_element', () => {
    it('should return integral elements', () => {
      const K = QuadraticField.create(5n);
      const O = new Order(K);

      for (let i = 0; i < 5; i++) {
        const elem = O.random_element(5n);
        expect(O.contains(elem)).toBe(true);
      }
    });
  });

  describe('characteristic', () => {
    it('should return 0', () => {
      const K = QuadraticField.create(2n);
      const O = new Order(K);
      expect(O.characteristic()).toBe(0n);
    });
  });

  describe('degree', () => {
    it('should return same as rank', () => {
      const K = QuadraticField.create(2n);
      const O = new Order(K);
      expect(O.degree()).toBe(O.rank());
      expect(O.absolute_degree()).toBe(O.rank());
    });
  });
});

describe('AbsoluteOrder', () => {
  it('should always be maximal', () => {
    const K = QuadraticField.create(7n);
    const O = new AbsoluteOrder(K);
    expect(O.is_maximal()).toBe(true);
  });
});

describe('EquationOrder', () => {
  it('should create equation order with power basis', () => {
    const K = QuadraticField.create(3n);
    const O = EquationOrder(K);

    expect(O.rank()).toBe(2);

    const basis = O.basis();
    expect(basis[0]!.is_one()).toBe(true);
    expect(basis[1]!.eq(K.gen())).toBe(true);
  });

  it('should have same discriminant as polynomial discriminant', () => {
    const K = QuadraticField.create(2n);
    const O = EquationOrder(K);

    // For x^2 - 2, disc = 8
    expect(O.discriminant()).toBe(8n);
  });
});

describe('Order string representation', () => {
  it('should show maximal order', () => {
    const K = QuadraticField.create(2n);
    const O = new Order(K);
    expect(O.toString()).toContain('Maximal Order');
  });
});
