/**
 * Tests for binary quadratic forms
 */

import { describe, expect, test } from 'bun:test';
import { BinaryQF, BinaryQF_reduced_representatives, class_number } from './binary_qf.js';

describe('BinaryQF', () => {
  describe('construction', () => {
    test('from three arguments', () => {
      const Q = new BinaryQF(1n, 2n, 3n);
      expect(Q.a).toBe(1n);
      expect(Q.b).toBe(2n);
      expect(Q.c).toBe(3n);
    });

    test('from tuple', () => {
      const Q = new BinaryQF([1n, 2n, 3n]);
      expect(Q.a).toBe(1n);
      expect(Q.b).toBe(2n);
      expect(Q.c).toBe(3n);
    });

    test('accepts bigint literals', () => {
      const Q = new BinaryQF(1n, 2n, 3n);
      expect(Q.a).toBe(1n);
      expect(Q.b).toBe(2n);
      expect(Q.c).toBe(3n);
    });
  });

  describe('principal form', () => {
    test('D = -4', () => {
      const Q = BinaryQF.principal(-4n);
      expect(Q.toTuple()).toEqual([1n, 0n, 1n]);
      expect(Q.discriminant()).toBe(-4n);
    });

    test('D = -3', () => {
      const Q = BinaryQF.principal(-3n);
      expect(Q.toTuple()).toEqual([1n, 1n, 1n]);
      expect(Q.discriminant()).toBe(-3n);
    });

    test('D = 8', () => {
      const Q = BinaryQF.principal(8n);
      expect(Q.toTuple()).toEqual([1n, 0n, -2n]);
      expect(Q.discriminant()).toBe(8n);
    });

    test('D = 5', () => {
      const Q = BinaryQF.principal(5n);
      expect(Q.toTuple()).toEqual([1n, 1n, -1n]);
      expect(Q.discriminant()).toBe(5n);
    });
  });

  describe('discriminant', () => {
    test('x^2 + 2xy + 3y^2', () => {
      const Q = new BinaryQF(1n, 2n, 3n);
      expect(Q.discriminant()).toBe(-8n);
    });

    test('x^2 - y^2', () => {
      const Q = new BinaryQF(1n, 0n, -1n);
      expect(Q.discriminant()).toBe(4n);
    });
  });

  describe('evaluate', () => {
    test('basic evaluation', () => {
      const Q = new BinaryQF(2n, 3n, 4n);
      expect(Q.evaluate(1n, 2n)).toBe(24n);
    });

    test('x^2 + y^2 at (1, 1)', () => {
      const Q = new BinaryQF(1n, 0n, 1n);
      expect(Q.evaluate(1n, 1n)).toBe(2n);
    });
  });

  describe('is_primitive', () => {
    test('primitive form', () => {
      const Q = new BinaryQF(1n, 1n, 1n);
      expect(Q.is_primitive()).toBe(true);
    });

    test('non-primitive form', () => {
      const Q = new BinaryQF(6n, 3n, 9n);
      expect(Q.is_primitive()).toBe(false);
    });
  });

  describe('is_reduced', () => {
    test('positive definite reduced', () => {
      const Q = new BinaryQF(2n, 1n, 3n);
      expect(Q.is_reduced()).toBe(true);
    });

    test('positive definite not reduced', () => {
      const Q = new BinaryQF(1n, 2n, 3n);
      expect(Q.is_reduced()).toBe(false);
    });

    test('x^2 + y^2 is reduced', () => {
      const Q = new BinaryQF(1n, 0n, 1n);
      expect(Q.is_reduced()).toBe(true);
    });
  });

  describe('reduced_form', () => {
    test('x^2 + 2xy + 3y^2 reduces to x^2 + 2y^2', () => {
      const Q = new BinaryQF(1n, 2n, 3n);
      const R = Q.reduced_form();
      expect(R.toTuple()).toEqual([1n, 0n, 2n]);
    });

    test('33x^2 + 11xy + 5y^2', () => {
      const Q = new BinaryQF(33n, 11n, 5n);
      const R = Q.reduced_form();
      expect(R.is_reduced()).toBe(true);
      expect(R.discriminant()).toBe(Q.discriminant());
    });

    test('already reduced form returns itself', () => {
      const Q = new BinaryQF(1n, 0n, 1n);
      const R = Q.reduced_form();
      expect(R.equals(Q)).toBe(true);
    });

    test('with transformation matrix', () => {
      const Q = new BinaryQF(1n, 2n, 3n);
      const [R, M] = Q.reduced_form({ transformation: true });
      // Verify Q * M = R
      const transformed = Q.matrix_action_right(M);
      expect(transformed.equals(R)).toBe(true);
    });
  });

  describe('inverse', () => {
    test('inverse negates b coefficient', () => {
      const Q = new BinaryQF(1n, 2n, 3n);
      const inv = Q.inverse();
      expect(inv.toTuple()).toEqual([1n, -2n, 3n]);
    });
  });

  describe('compose', () => {
    test('composition of same form', () => {
      const Q = new BinaryQF(2n, 1n, 3n);
      const composed = Q.compose(Q);
      expect(composed.discriminant()).toBe(Q.discriminant());
    });

    test('composition with identity', () => {
      const Q = new BinaryQF(2n, 1n, 3n);
      const D = Q.discriminant();
      const identity = BinaryQF.principal(D);
      const composed = Q.compose(identity);
      expect(composed.reduced_form().equals(Q.reduced_form())).toBe(true);
    });
  });

  describe('is_equivalent', () => {
    test('same form is equivalent', () => {
      const Q = new BinaryQF(1n, 2n, 3n);
      expect(Q.is_equivalent(Q)).toBe(true);
    });

    test('form and its reduction are equivalent', () => {
      const Q = new BinaryQF(33n, 11n, 5n);
      const R = Q.reduced_form();
      expect(Q.is_equivalent(R)).toBe(true);
    });

    test('different discriminants are not equivalent', () => {
      const Q1 = new BinaryQF(1n, 0n, 1n);
      const Q2 = new BinaryQF(1n, 0n, 2n);
      expect(Q1.is_equivalent(Q2)).toBe(false);
    });
  });

  describe('toString', () => {
    test('x^2 + 2xy + 3y^2', () => {
      const Q = new BinaryQF(1n, 2n, 3n);
      expect(Q.toString()).toBe('x^2 + 2*x*y + 3*y^2');
    });

    test('-x^2 + y^2', () => {
      const Q = new BinaryQF(-1n, 0n, 1n);
      expect(Q.toString()).toBe('-x^2 + y^2');
    });

    test('zero form', () => {
      const Q = new BinaryQF(0n, 0n, 0n);
      expect(Q.toString()).toBe('0');
    });
  });
});

describe('BinaryQF_reduced_representatives', () => {
  test('D = -4 has one class', () => {
    const reps = BinaryQF_reduced_representatives(-4n);
    expect(reps.length).toBe(1);
    expect(reps[0]!.toTuple()).toEqual([1n, 0n, 1n]);
  });

  test('D = -23 has three classes', () => {
    const reps = BinaryQF_reduced_representatives(-23n);
    expect(reps.length).toBe(3);
    // Verify all have same discriminant
    for (const Q of reps) {
      expect(Q.discriminant()).toBe(-23n);
      expect(Q.is_reduced()).toBe(true);
    }
  });

  test('D = -163 has one class (Heegner number)', () => {
    const reps = BinaryQF_reduced_representatives(-163n);
    expect(reps.length).toBe(1);
  });
});

describe('class_number', () => {
  test('h(-4) = 1', () => {
    expect(class_number(-4n)).toBe(1n);
  });

  test('h(-23) = 3', () => {
    expect(class_number(-23n)).toBe(3n);
  });

  test('h(-163) = 1 (Heegner number)', () => {
    expect(class_number(-163n)).toBe(1n);
  });

  test('h(-67) = 1 (Heegner number)', () => {
    expect(class_number(-67n)).toBe(1n);
  });

  test('h(-11) = 1', () => {
    expect(class_number(-11n)).toBe(1n);
  });
});
