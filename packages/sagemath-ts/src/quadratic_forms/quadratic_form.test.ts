/**
 * Tests for `sage/quadratic_forms/quadratic_form` and
 * `sage/quadratic_forms/quadratic_form__local_field_invariants`.
 *
 * Every expected value below is either copied verbatim from a SageMath
 * doctest in `reference/sage/src/sage/quadratic_forms/` (the file:line is
 * quoted) or was produced by running the installed `sage` binary.
 */

import { describe, expect, test } from 'bun:test';
import { NotImplementedError } from '../errors.js';
import { ZZ } from '../rings/integer_ring.js';
import { Rational } from '../rings/rational.js';
import { QQ } from '../rings/rational_field.js';
import {
  DiagonalQuadraticForm,
  QuadraticForm,
  quadratic_form_from_invariants,
} from './quadratic_form.js';
import {
  hasse_invariant,
  is_padic_square,
  qfgaussred,
  signature_vector,
} from './quadratic_form__local_field_invariants.js';

const R = (a: bigint, b: bigint) => new Rational(a, b);

describe('QuadraticForm construction and representation', () => {
  test('_repr_ (quadratic_form.py:714)', () => {
    const Q = new QuadraticForm(ZZ, 2n, [1n, 3n, 5n]);
    expect(Q.toString()).toBe(
      'Quadratic form in 2 variables over Integer Ring with coefficients: \n[ 1 3 ]\n[ * 5 ]'
    );
  });

  test('dim / base_ring / coefficients (quadratic_form.py:1425,1444,1456)', () => {
    const Q = new QuadraticForm(ZZ, 2n, [1n, 2n, 3n]);
    expect(Q.dim()).toBe(2n);
    expect(Q.base_ring()).toBe(ZZ);
    expect(Q.coefficients().map((c) => c.toString())).toEqual(['1', '2', '3']);
  });

  test('QuadraticForm(ZZ, 4, range(10)) (quadratic_form.py:504)', () => {
    const Q = new QuadraticForm(ZZ, 4n, [0n, 1n, 2n, 3n, 4n, 5n, 6n, 7n, 8n, 9n]);
    expect(Q.dim()).toBe(4n);
    expect(Q.toString()).toBe(
      'Quadratic form in 4 variables over Integer Ring with coefficients: \n' +
        '[ 0 1 2 3 ]\n[ * 4 5 6 ]\n[ * * 7 8 ]\n[ * * * 9 ]'
    );
  });

  test('negative size raises (quadratic_form.py:530)', () => {
    expect(() => new QuadraticForm(ZZ, -1n)).toThrow(
      'the size must be a nonnegative integer, not -1'
    );
  });

  test('wrong number of entries raises', () => {
    expect(() => new QuadraticForm(ZZ, 2n, [1n, 2n])).toThrow(/must be a list of size/);
  });

  test('__getitem__ (quadratic_form.py:765)', () => {
    const Q = new QuadraticForm(ZZ, 3n, [1n, 2n, 3n, 4n, 5n, 6n]);
    // matrix(ZZ, 3, 3, [Q[i,j] for i in range(3) for j in range(3)])
    // == [1 2 3; 2 4 5; 3 5 6]
    const rows = [0, 1, 2].map((i) => [0, 1, 2].map((j) => Q.get(i, j).numerator));
    expect(rows).toEqual([
      [1n, 2n, 3n],
      [2n, 4n, 5n],
      [3n, 5n, 6n],
    ]);
  });

  test('__setitem__ (quadratic_form.py:788)', () => {
    const Q = new QuadraticForm(ZZ, 3n, [1n, 2n, 3n, 4n, 5n, 6n]);
    Q.set(2, 1, 17n);
    expect(Q.toString()).toBe(
      'Quadratic form in 3 variables over Integer Ring with coefficients: \n' +
        '[ 1 2 3 ]\n[ * 4 17 ]\n[ * * 6 ]'
    );
  });

  test('__eq__ (quadratic_form.py:838)', () => {
    const Q = new QuadraticForm(ZZ, 2n, [1n, 4n, 10n]);
    const Q1 = new QuadraticForm(QQ, 2n, [1n, 4n, 10n]);
    const Q2 = new QuadraticForm(ZZ, 2n, [1n, 4n, -10n]);
    expect(Q.equals(Q)).toBe(true);
    expect(Q.equals(Q1)).toBe(false);
    expect(Q.equals(Q2)).toBe(false);
    expect(Q1.equals(Q2)).toBe(false);
  });

  test('DiagonalQuadraticForm (quadratic_form.py:1727)', () => {
    const Q = DiagonalQuadraticForm(ZZ, [1n, 3n, 5n, 7n]);
    expect(Q.toString()).toBe(
      'Quadratic form in 4 variables over Integer Ring with coefficients: \n' +
        '[ 1 0 0 0 ]\n[ * 3 0 0 ]\n[ * * 5 0 ]\n[ * * * 7 ]'
    );
  });

  test('from a symmetric matrix with even diagonal', () => {
    // QuadraticForm(ZZ, 2, [2,3,5]).matrix() == [4 3; 3 10]
    const Q = new QuadraticForm(ZZ, 2n, [2n, 3n, 5n]);
    const A = Q.matrix();
    expect(A.rows().map((r) => r.map((x) => x.numerator))).toEqual([
      [4n, 3n],
      [3n, 10n],
    ]);
    const Q2 = new QuadraticForm(ZZ, A);
    expect(Q2.equals(Q)).toBe(true);
    // Odd diagonal is rejected over ZZ (quadratic_form.py:1067)
    A.set(0, 0, R(1n, 1n));
    expect(QuadraticForm._is_even_symmetric_matrix_(A, ZZ)).toBe(false);
  });
});

describe('QuadraticForm sums', () => {
  test('__add__ is the direct sum (quadratic_form.py:864)', () => {
    const Q = new QuadraticForm(ZZ, 2n, [1n, 4n, 10n]);
    const Q2 = new QuadraticForm(ZZ, 2n, [1n, 4n, -10n]);
    expect(Q.add(Q2).toString()).toBe(
      'Quadratic form in 4 variables over Integer Ring with coefficients: \n' +
        '[ 1 4 0 0 ]\n[ * 10 0 0 ]\n[ * * 1 4 ]\n[ * * * -10 ]'
    );
  });

  test('sum_by_coefficients_with (quadratic_form.py:901)', () => {
    const Q = new QuadraticForm(ZZ, 2n, [1n, 4n, 10n]);
    const Q2 = new QuadraticForm(ZZ, 2n, [1n, 4n, -10n]);
    expect(Q.sum_by_coefficients_with(Q2).toString()).toBe(
      'Quadratic form in 2 variables over Integer Ring with coefficients: \n[ 2 8 ]\n[ * 0 ]'
    );
  });

  test('different base rings cannot be added', () => {
    const Q = new QuadraticForm(ZZ, 2n, [1n, 4n, 10n]);
    const Q2 = new QuadraticForm(QQ, 2n, [1n, 4n, 10n]);
    expect(() => Q.add(Q2)).toThrow(/do not have the same base rings/);
  });
});

describe('QuadraticForm evaluation', () => {
  test('on vectors (quadratic_form.py:954)', () => {
    const Q = new QuadraticForm(QQ, 3n, [0n, 1n, 2n, 3n, 4n, 5n]);
    expect(Q.__call__([1n, 2n, 3n]).toString()).toBe('89');
    expect(Q.__call__([1n, 0n, 0n]).toString()).toBe('0');
    expect(Q.__call__([1n, 1n, 1n]).toString()).toBe('15');
  });

  test('on matrices (quadratic_form.py:984,1001)', () => {
    const Q = new QuadraticForm(QQ, 2n, [1n, 2n, 3n]);
    expect(
      (
        Q.__call__([
          [-1n, 0n],
          [0n, 1n],
        ]) as QuadraticForm
      ).toString()
    ).toBe(
      'Quadratic form in 2 variables over Rational Field with coefficients: \n[ 1 -2 ]\n[ * 3 ]'
    );

    const Q2 = new QuadraticForm(ZZ, 2n, [1n, 0n, 1n]);
    expect(
      (
        Q2.__call__([
          [1n, 1n],
          [0n, 1n],
        ]) as QuadraticForm
      ).toString()
    ).toBe('Quadratic form in 2 variables over Integer Ring with coefficients: \n[ 1 2 ]\n[ * 2 ]');

    const Q3 = DiagonalQuadraticForm(ZZ, [1n, 1n, 1n]);
    expect((Q3.__call__([[1n], [2n], [3n]]) as QuadraticForm).toString()).toBe(
      'Quadratic form in 1 variables over Integer Ring with coefficients: \n[ 14 ]'
    );
    expect(Q3.__call__([1n, 2n, 3n]).toString()).toBe('14');
  });

  test('wrong vector length raises (quadratic_form.py:1050)', () => {
    const Q = DiagonalQuadraticForm(ZZ, [1n, 7n, 7n]);
    expect(() => Q.__call__([1n, 2n])).toThrow('your vector needs to have length 3');
  });

  test('bilinear_map (quadratic_form.py:1666)', () => {
    const Q = new QuadraticForm(ZZ, 3n, [1n, 4n, 0n, 1n, 4n, 1n]);
    expect(Q.bilinear_map([1n, 2n, 0n], [0n, 1n, 1n]).toString()).toBe('8');
    const Q2 = new QuadraticForm(QQ, 2n, [R(1n, 2n), 2n, 1n]);
    expect(Q2.bilinear_map([1n, 1n], [R(1n, 2n), 2n]).toString()).toBe('19/4');
    const Q3 = DiagonalQuadraticForm(ZZ, [1n, 7n, 7n]);
    expect(() => Q3.bilinear_map([1n, 2n], [1n, 1n, 1n])).toThrow('vectors must have length 3');
  });
});

describe('QuadraticForm matrices, determinants and level', () => {
  test('Hessian_matrix (quadratic_form.py:1138)', () => {
    const Q = new QuadraticForm(QQ, 2n, [1n, 2n, 3n]);
    expect(
      Q.Hessian_matrix()
        .rows()
        .map((r) => r.map((x) => x.toString()))
    ).toEqual([
      ['2', '2'],
      ['2', '6'],
    ]);
    const Q2 = new QuadraticForm(ZZ, 3n, [0n, 1n, 2n, 3n, 4n, 5n]);
    expect(
      Q2.matrix()
        .rows()
        .map((r) => r.map((x) => x.toString()))
    ).toEqual([
      ['0', '1', '2'],
      ['1', '6', '4'],
      ['2', '4', '10'],
    ]);
  });

  test('Gram_matrix / Gram_matrix_rational (quadratic_form.py:1164,1188)', () => {
    const Q = DiagonalQuadraticForm(ZZ, [1n, 3n, 5n, 7n]);
    const expected = [
      ['1', '0', '0', '0'],
      ['0', '3', '0', '0'],
      ['0', '0', '5', '0'],
      ['0', '0', '0', '7'],
    ];
    expect(
      Q.Gram_matrix()
        .rows()
        .map((r) => r.map((x) => x.toString()))
    ).toEqual(expected);
    expect(
      Q.Gram_matrix_rational()
        .rows()
        .map((r) => r.map((x) => x.toString()))
    ).toEqual(expected);
  });

  test('has_integral_Gram_matrix (quadratic_form.py:1225)', () => {
    expect(new QuadraticForm(ZZ, 2n, [7n, 8n, 9n]).has_integral_Gram_matrix()).toBe(true);
    expect(new QuadraticForm(ZZ, 2n, [4n, 5n, 6n]).has_integral_Gram_matrix()).toBe(false);
    expect(() => new QuadraticForm(ZZ, 2n, [4n, 5n, 6n]).Gram_matrix()).toThrow(
      'this form does not have an integral Gram matrix'
    );
  });

  test('det and Gram_det (quadratic_form.py:1469,1497)', () => {
    const Q = new QuadraticForm(ZZ, 2n, [1n, 2n, 3n]);
    expect(Q.det().toString()).toBe('8');
    expect(Q.Gram_det().toString()).toBe('2');
    // sage: QuadraticForm(ZZ,3,[1,2,3,4,5,6]).det() -> 82
    expect(new QuadraticForm(ZZ, 3n, [1n, 2n, 3n, 4n, 5n, 6n]).det().toString()).toBe('82');
    expect(DiagonalQuadraticForm(ZZ, [1n, 3n, 5n, 7n]).Gram_det().toString()).toBe('105');
  });

  test('level (quadratic_form.py:1558)', () => {
    expect(new QuadraticForm(ZZ, 2n, [1n, 2n, 3n]).level()).toBe(8n);
    expect(DiagonalQuadraticForm(ZZ, [1n, 3n, 5n, 7n]).level()).toBe(420n);
    // sage: QuadraticForm(ZZ,3,[1,2,3,4,5,6]).level() -> 164
    expect(new QuadraticForm(ZZ, 3n, [1n, 2n, 3n, 4n, 5n, 6n]).level()).toBe(164n);
  });

  test('gcd / primitive / adjoint_primitive (quadratic_form.py:1256,1385,1408)', () => {
    expect(new QuadraticForm(ZZ, 4n, [1n, 3n, 5n, 7n, 9n, 11n, 13n, 15n, 17n, 19n]).gcd()).toBe(1n);
    expect(new QuadraticForm(ZZ, 4n, [0n, 2n, 4n, 6n, 8n, 10n, 12n, 14n, 16n, 18n]).gcd()).toBe(2n);
    expect(new QuadraticForm(ZZ, 2n, [2n, 3n, 4n]).is_primitive()).toBe(true);
    expect(new QuadraticForm(ZZ, 2n, [2n, 4n, 8n]).is_primitive()).toBe(false);
    expect(new QuadraticForm(ZZ, 2n, [2n, 4n, 8n]).primitive().toString()).toBe(
      'Quadratic form in 2 variables over Integer Ring with coefficients: \n[ 1 2 ]\n[ * 4 ]'
    );
    expect(new QuadraticForm(ZZ, 2n, [1n, 2n, 3n]).adjoint_primitive().toString()).toBe(
      'Quadratic form in 2 variables over Integer Ring with coefficients: \n[ 3 -2 ]\n[ * 1 ]'
    );
  });

  test('change_ring to QQ', () => {
    const Q = DiagonalQuadraticForm(ZZ, [1n, 1n]);
    const Q1 = Q.change_ring(QQ);
    expect(Q1.base_ring()).toBe(QQ);
    expect(Q1.toString()).toBe(
      'Quadratic form in 2 variables over Rational Field with coefficients: \n[ 1 0 ]\n[ * 1 ]'
    );
    expect(() => new QuadraticForm(QQ, 1n, [R(1n, 2n)]).change_ring(ZZ)).toThrow(
      /no canonical coercion/
    );
  });
});

describe('rational_diagonal_form', () => {
  test('QuadraticForm(ZZ, 2, [0,1,-1]) (local_field_invariants.py:56)', () => {
    const Q = new QuadraticForm(ZZ, 2n, [0n, 1n, -1n]);
    expect(Q.rational_diagonal_form().toString()).toBe(
      'Quadratic form in 2 variables over Rational Field with coefficients: \n[ 1/4 0 ]\n[ * -1 ]'
    );
  });

  test('a diagonal form stays diagonal (local_field_invariants.py:67)', () => {
    const Q = DiagonalQuadraticForm(ZZ, [1n, 3n, 5n, 7n]);
    expect(Q.rational_diagonal_form().toString()).toBe(
      'Quadratic form in 4 variables over Rational Field with coefficients: \n' +
        '[ 1 0 0 0 ]\n[ * 3 0 0 ]\n[ * * 5 0 ]\n[ * * * 7 ]'
    );
  });

  test('QuadraticForm(ZZ, 4, range(10)) with transformation (local_field_invariants.py:77)', () => {
    const Q = new QuadraticForm(ZZ, 4n, [0n, 1n, 2n, 3n, 4n, 5n, 6n, 7n, 8n, 9n]);
    const [D, T] = Q.rational_diagonal_form(true);
    expect(D.toString()).toBe(
      'Quadratic form in 4 variables over Rational Field with coefficients: \n' +
        '[ -1/16 0 0 0 ]\n[ * 4 0 0 ]\n[ * * 13 0 ]\n[ * * * 563/52 ]'
    );
    expect(T.rows().map((r) => r.map((x) => x.toString()))).toEqual([
      ['1', '0', '11', '149/26'],
      ['-1/8', '1', '-2', '-10/13'],
      ['0', '0', '1', '-29/26'],
      ['0', '0', '0', '1'],
    ]);
    // T.transpose() * Q.matrix() * T == D.matrix()
    expect(T.transpose().mul(Q.matrix()).mul(T).eq(D.matrix())).toBe(true);
  });

  test('QuadraticForm(ZZ, 4, [1,1,0,0,1,0,0,1,0,18]) (local_field_invariants.py:105)', () => {
    const Q = new QuadraticForm(ZZ, 4n, [1n, 1n, 0n, 0n, 1n, 0n, 0n, 1n, 0n, 18n]);
    const [D, T] = Q.rational_diagonal_form(true);
    expect(D.toString()).toBe(
      'Quadratic form in 4 variables over Rational Field with coefficients: \n' +
        '[ 1 0 0 0 ]\n[ * 3/4 0 0 ]\n[ * * 1 0 ]\n[ * * * 18 ]'
    );
    expect(T.rows().map((r) => r.map((x) => x.toString()))).toEqual([
      ['1', '-1/2', '0', '0'],
      ['0', '1', '0', '0'],
      ['0', '0', '1', '0'],
      ['0', '0', '0', '1'],
    ]);
  });

  test('the singular case PARI cannot do (local_field_invariants.py:127)', () => {
    const Q = new QuadraticForm(QQ, 2n, [R(1n, 2n), 1n, R(1n, 2n)]);
    expect(Q.rational_diagonal_form().toString()).toBe(
      'Quadratic form in 2 variables over Rational Field with coefficients: \n[ 1/2 0 ]\n[ * 0 ]'
    );
  });

  test('changing the returned form does not corrupt the cache (local_field_invariants.py:150)', () => {
    const Q1 = new QuadraticForm(ZZ, 4n, [1n, 1n, 0n, 0n, 1n, 0n, 0n, 1n, 0n, 18n]);
    const D = Q1.rational_diagonal_form();
    D.set(0, 0, 13n);
    expect(Q1.rational_diagonal_form().toString()).toBe(
      'Quadratic form in 4 variables over Rational Field with coefficients: \n' +
        '[ 1 0 0 0 ]\n[ * 3/4 0 0 ]\n[ * * 1 0 ]\n[ * * * 18 ]'
    );
  });

  test('qfgaussred matches PARI on QuadraticForm(ZZ, 4, range(10))', () => {
    // sage: QuadraticForm(ZZ, 4, range(10)).__pari__().qfgaussred()
    // [-1/8, 0, -11, -18; 1/8, 8, 5/8, 3/4; 0, 0, 26, 29/26; 0, 0, 0, 563/26]
    const Q = new QuadraticForm(ZZ, 4n, [0n, 1n, 2n, 3n, 4n, 5n, 6n, 7n, 8n, 9n]);
    expect(
      qfgaussred(Q.matrix())
        .rows()
        .map((r) => r.map((x) => x.toString()))
    ).toEqual([
      ['-1/8', '0', '-11', '-18'],
      ['1/8', '8', '5/8', '3/4'],
      ['0', '0', '26', '29/26'],
      ['0', '0', '0', '563/26'],
    ]);
  });

  test('T^t A T is diagonal for a random-ish spread of forms', () => {
    const cases: bigint[][] = [
      [1n, 2n, 3n],
      [0n, 1n, 0n],
      [0n, 0n, 1n, 0n, 1n, 0n],
      [1n, 2n, 3n, 4n, 5n, 6n],
      [2n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 5n],
    ];
    for (const entries of cases) {
      const n = (Math.sqrt(8 * entries.length + 1) - 1) / 2;
      const Q = new QuadraticForm(ZZ, BigInt(n), entries);
      const [D, T] = Q.rational_diagonal_form(true);
      expect(T.transpose().mul(Q.matrix()).mul(T).eq(D.matrix())).toBe(true);
    }
  });
});

describe('signature', () => {
  test('signature_vector (local_field_invariants.py:297)', () => {
    expect(DiagonalQuadraticForm(ZZ, [1n, 0n, 0n, -4n]).signature_vector()).toEqual([1n, 1n, 2n]);
    expect(DiagonalQuadraticForm(ZZ, [1n, 2n, -3n, -4n]).signature_vector()).toEqual([2n, 2n, 0n]);
    const Q = new QuadraticForm(ZZ, 4n, [0n, 1n, 2n, 3n, 4n, 5n, 6n, 7n, 8n, 9n]);
    expect(Q.signature_vector()).toEqual([3n, 1n, 0n]);
    expect(signature_vector(Q)).toEqual([3n, 1n, 0n]);
  });

  test('signature (local_field_invariants.py:346)', () => {
    expect(DiagonalQuadraticForm(ZZ, [1n, 0n, 0n, -4n, 3n, 11n, 3n]).signature()).toBe(3n);
    expect(DiagonalQuadraticForm(ZZ, [1n, 2n, -3n, -4n]).signature()).toBe(0n);
    const Q = new QuadraticForm(ZZ, 4n, [0n, 1n, 2n, 3n, 4n, 5n, 6n, 7n, 8n, 9n]);
    expect(Q.signature()).toBe(2n);
  });
});

describe('Hasse invariants', () => {
  const primes = [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n];

  test('QuadraticForm(ZZ, 2, [1,2,3]) (local_field_invariants.py:381)', () => {
    const Q = new QuadraticForm(ZZ, 2n, [1n, 2n, 3n]);
    expect(Q.rational_diagonal_form().toString()).toBe(
      'Quadratic form in 2 variables over Rational Field with coefficients: \n[ 1 0 ]\n[ * 2 ]'
    );
    expect(primes.map((p) => Q.hasse_invariant(p))).toEqual(primes.map(() => 1n));
    expect(primes.map((p) => Q.hasse_invariant__OMeara(p))).toEqual(primes.map(() => 1n));
  });

  test('DiagonalQuadraticForm(ZZ, [1,-1]) (local_field_invariants.py:508)', () => {
    const Q = DiagonalQuadraticForm(ZZ, [1n, -1n]);
    expect(primes.map((p) => Q.hasse_invariant(p))).toEqual([1n, 1n, 1n, 1n, 1n, 1n, 1n, 1n]);
    expect(primes.map((p) => Q.hasse_invariant__OMeara(p))).toEqual([
      -1n,
      1n,
      1n,
      1n,
      1n,
      1n,
      1n,
      1n,
    ]);
  });

  test('DiagonalQuadraticForm(ZZ, [1,-1,5]) (local_field_invariants.py:427)', () => {
    const Q = DiagonalQuadraticForm(ZZ, [1n, -1n, 5n]);
    expect(primes.map((p) => Q.hasse_invariant(p))).toEqual([1n, 1n, 1n, 1n, 1n, 1n, 1n, 1n]);
    expect(primes.map((p) => Q.hasse_invariant__OMeara(p))).toEqual([
      -1n,
      1n,
      1n,
      1n,
      1n,
      1n,
      1n,
      1n,
    ]);
  });

  test('DiagonalQuadraticForm(ZZ, [1,-1,-1]) (local_field_invariants.py:521)', () => {
    const Q = DiagonalQuadraticForm(ZZ, [1n, -1n, -1n]);
    expect(primes.map((p) => Q.hasse_invariant(p))).toEqual([-1n, 1n, 1n, 1n, 1n, 1n, 1n, 1n]);
    expect(primes.map((p) => Q.hasse_invariant__OMeara(p))).toEqual([
      -1n,
      1n,
      1n,
      1n,
      1n,
      1n,
      1n,
      1n,
    ]);
    expect(hasse_invariant(Q, 2n)).toBe(-1n);
  });

  test('Hilbert reciprocity: the product over all places is 1', () => {
    // Verified against sage:
    //   prod(Q.hasse_invariant(p) for p in [-1] + prime_divisors(2*Q.det())) == 1
    const cases: bigint[][] = [
      [1n, 3n, 5n],
      [1n, -1n, 5n],
      [2n, 3n, -7n],
      [1n, 2n, 3n, 4n],
      [5n, -3n, 2n, 11n],
    ];
    for (const diag of cases) {
      const Q = DiagonalQuadraticForm(QQ, diag);
      const places: bigint[] = [-1n];
      const det2 = Q.det().mul(2n).numerator;
      for (const p of [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n]) {
        if (det2 % p === 0n) places.push(p);
      }
      const prod = places.reduce((acc, p) => acc * Q.hasse_invariant(p), 1n);
      expect(prod).toBe(1n);
    }
  });
});

describe('is_padic_square', () => {
  test('sage/rings/rational.pyx:1751 doctests', () => {
    expect(is_padic_square(R(2n, 1n), 7n)).toBe(true);
    expect(is_padic_square(R(98n, 1n), 7n)).toBe(true);
    expect(is_padic_square(R(2n, 1n), 5n)).toBe(false);
    expect(is_padic_square(R(5n, 7n), 2n)).toBe(false);
    expect(is_padic_square(R(0n, 1n), 3n)).toBe(true);
  });
});

describe('hyperbolicity and (an)isotropy', () => {
  test('is_hyperbolic on [1,1] (local_field_invariants.py:576)', () => {
    const Q = DiagonalQuadraticForm(ZZ, [1n, 1n]);
    expect(Q.is_hyperbolic(-1n)).toBe(false);
    expect(Q.is_hyperbolic(2n)).toBe(false);
    expect(Q.is_hyperbolic(3n)).toBe(false);
    expect(Q.is_hyperbolic(5n)).toBe(true);
    expect(Q.is_hyperbolic(7n)).toBe(false);
    expect(Q.is_hyperbolic(13n)).toBe(true);
  });

  test('the hyperbolic plane really is hyperbolic everywhere', () => {
    // sage: [DiagonalQuadraticForm(ZZ,[1,-1]).is_hyperbolic(p) for p in [-1,2,3,5,7]]
    // [True, True, True, True, True]
    const H = DiagonalQuadraticForm(ZZ, [1n, -1n]);
    for (const p of [-1n, 2n, 3n, 5n, 7n]) {
      expect(H.is_hyperbolic(p)).toBe(true);
    }
    const H2 = DiagonalQuadraticForm(ZZ, [1n, -1n, 1n, -1n]);
    for (const p of [-1n, 2n, 3n, 5n, 7n]) {
      expect(H2.is_hyperbolic(p)).toBe(true);
    }
  });

  test('is_anisotropic / is_isotropic (local_field_invariants.py:619,685)', () => {
    const Q = DiagonalQuadraticForm(ZZ, [1n, 1n]);
    expect(Q.is_anisotropic(2n)).toBe(true);
    expect(Q.is_anisotropic(3n)).toBe(true);
    expect(Q.is_anisotropic(5n)).toBe(false);
    expect(Q.is_isotropic(5n)).toBe(true);

    const Q2 = DiagonalQuadraticForm(ZZ, [1n, -1n]);
    expect(Q2.is_anisotropic(2n)).toBe(false);
    expect(Q2.is_anisotropic(3n)).toBe(false);
    expect(Q2.is_anisotropic(5n)).toBe(false);
  });

  test('every form in >= 5 variables is isotropic over Q_p', () => {
    const Q = DiagonalQuadraticForm(ZZ, [1n, 1n, 1n, 1n, 1n]);
    for (const p of [2n, 3n, 5n, 7n]) {
      expect(Q.is_anisotropic(p)).toBe(false);
    }
  });

  test('anisotropic_primes (local_field_invariants.py:729)', () => {
    expect(DiagonalQuadraticForm(ZZ, [1n, 1n, 1n]).anisotropic_primes()).toEqual([2n, -1n]);
    expect(DiagonalQuadraticForm(ZZ, [1n, 1n, 1n, 1n]).anisotropic_primes()).toEqual([2n, -1n]);
    expect(DiagonalQuadraticForm(ZZ, [1n, 1n, 1n, 1n, 1n]).anisotropic_primes()).toEqual([-1n]);
    // sage: DiagonalQuadraticForm(QQ,[1,3,5]).anisotropic_primes() -> [5, -1]
    expect(DiagonalQuadraticForm(QQ, [1n, 3n, 5n]).anisotropic_primes()).toEqual([5n, -1n]);
    // sage: DiagonalQuadraticForm(QQ,[1,-1,5]).anisotropic_primes() -> []
    expect(DiagonalQuadraticForm(QQ, [1n, -1n, 5n]).anisotropic_primes()).toEqual([]);
  });

  test('anisotropic_primes with a non-integral 2*det', () => {
    // sage: Q = QuadraticForm(QQ, 1, [3/8]); Q.det(), 2*Q.det(), Q.anisotropic_primes()
    // (3/4, 3/2, [2, 3, -1])
    const Q = new QuadraticForm(QQ, 1n, [R(3n, 8n)]);
    expect(Q.det().toString()).toBe('3/4');
    expect(Q.anisotropic_primes()).toEqual([2n, 3n, -1n]);

    // sage: Q2 = QuadraticForm(QQ, 2, [1/3, 1/5, 2/7])
    // Q2.det() == 179/525, Q2.anisotropic_primes() == [3, 7, 179, -1],
    // Q2.signature() == 2, [Q2.hasse_invariant(p) for p in [-1,2,3,5,7]]
    //   == [1, 1, -1, 1, -1]
    const Q2 = new QuadraticForm(QQ, 2n, [R(1n, 3n), R(1n, 5n), R(2n, 7n)]);
    expect(Q2.det().toString()).toBe('179/525');
    expect(Q2.anisotropic_primes()).toEqual([3n, 7n, 179n, -1n]);
    expect(Q2.signature()).toBe(2n);
    expect([-1n, 2n, 3n, 5n, 7n].map((p) => Q2.hasse_invariant(p))).toEqual([1n, 1n, -1n, 1n, -1n]);
  });
});

describe('definiteness', () => {
  test('compute_definiteness (local_field_invariants.py:758)', () => {
    const Q = DiagonalQuadraticForm(ZZ, [1n, 1n, 1n, 1n, 1n]);
    expect(Q.is_positive_definite()).toBe(true);
    expect(Q.is_negative_definite()).toBe(false);
    expect(Q.is_indefinite()).toBe(false);
    expect(Q.is_definite()).toBe(true);

    const Z0 = DiagonalQuadraticForm(ZZ, []);
    expect(Z0.is_positive_definite()).toBe(true);
    expect(Z0.is_negative_definite()).toBe(true);
    expect(Z0.is_indefinite()).toBe(false);
    expect(Z0.is_definite()).toBe(true);

    const Dg = DiagonalQuadraticForm(ZZ, [1n, 0n, -1n]);
    expect(Dg.is_positive_definite()).toBe(false);
    expect(Dg.is_negative_definite()).toBe(false);
    expect(Dg.is_indefinite()).toBe(false);
    expect(Dg.is_definite()).toBe(false);
  });

  test('is_positive_definite / is_negative_definite (local_field_invariants.py:924,957)', () => {
    expect(DiagonalQuadraticForm(ZZ, [1n, 3n, 5n]).is_positive_definite()).toBe(true);
    expect(DiagonalQuadraticForm(ZZ, [1n, -3n, 5n]).is_positive_definite()).toBe(false);
    expect(DiagonalQuadraticForm(ZZ, [-1n, -3n, -5n]).is_negative_definite()).toBe(true);
    expect(DiagonalQuadraticForm(ZZ, [1n, -3n, 5n]).is_negative_definite()).toBe(false);
    expect(DiagonalQuadraticForm(ZZ, [-1n, -3n, -5n]).is_indefinite()).toBe(false);
    expect(DiagonalQuadraticForm(ZZ, [1n, -3n, 5n]).is_indefinite()).toBe(true);
    expect(DiagonalQuadraticForm(ZZ, [-1n, -3n, -5n]).is_definite()).toBe(true);
    expect(DiagonalQuadraticForm(ZZ, [1n, -3n, 5n]).is_definite()).toBe(false);
  });

  test('compute_definiteness_string_by_determinants (local_field_invariants.py:829)', () => {
    expect(
      DiagonalQuadraticForm(ZZ, [1n, 1n, 1n, 1n, 1n]).compute_definiteness_string_by_determinants()
    ).toBe('pos_def');
    expect(DiagonalQuadraticForm(ZZ, []).compute_definiteness_string_by_determinants()).toBe(
      'zero'
    );
    expect(
      DiagonalQuadraticForm(ZZ, [1n, 0n, -1n]).compute_definiteness_string_by_determinants()
    ).toBe('degenerate');
    expect(DiagonalQuadraticForm(ZZ, [1n, -1n]).compute_definiteness_string_by_determinants()).toBe(
      'indefinite'
    );
    expect(
      DiagonalQuadraticForm(ZZ, [-1n, -1n]).compute_definiteness_string_by_determinants()
    ).toBe('neg_def');
  });

  test('is_anisotropic(-1) agrees with is_definite', () => {
    for (const diag of [
      [1n, 1n, 1n],
      [1n, -1n, 1n],
      [-1n, -2n, -3n],
    ]) {
      const Q = DiagonalQuadraticForm(ZZ, diag);
      expect(Q.is_anisotropic(-1n)).toBe(Q.is_definite());
    }
  });
});

describe('is_rationally_isometric (Hasse-Minkowski)', () => {
  test('doctest (quadratic_form__equivalence_testing.py)', () => {
    const V = DiagonalQuadraticForm(QQ, [1n, 1n, 2n]);
    const W = DiagonalQuadraticForm(QQ, [2n, 2n, 2n]);
    expect(V.is_rationally_isometric(W)).toBe(true);
  });

  test('more cases verified with sage', () => {
    const D = (l: bigint[]) => DiagonalQuadraticForm(QQ, l);
    expect(D([1n, 1n, 1n]).is_rationally_isometric(D([1n, 1n, 3n]))).toBe(false);
    expect(D([1n, 1n, 1n]).is_rationally_isometric(D([1n, 1n, 4n]))).toBe(true);
    expect(D([1n, -1n]).is_rationally_isometric(D([2n, -2n]))).toBe(true);
    expect(D([1n, 2n, 3n]).is_rationally_isometric(D([1n, 2n, 3n]))).toBe(true);
  });

  test('a form is isometric to its rational diagonalisation', () => {
    const Q = new QuadraticForm(QQ, 3n, [1n, 2n, 3n, 4n, 5n, 6n]);
    expect(Q.is_rationally_isometric(Q.rational_diagonal_form())).toBe(true);
  });

  test('over ZZ SageMath has no implementation; we say so', () => {
    const V = DiagonalQuadraticForm(ZZ, [1n, 1n, 2n]);
    const W = DiagonalQuadraticForm(ZZ, [2n, 2n, 2n]);
    expect(() => V.is_rationally_isometric(W)).toThrow(NotImplementedError);
  });

  test('a degenerate form is rejected', () => {
    const V = DiagonalQuadraticForm(QQ, [1n, 0n]);
    expect(() => V.is_rationally_isometric(V)).toThrow('this only tests regular forms');
  });
});

describe('theta series (PARI qfrep)', () => {
  test('DiagonalQuadraticForm(ZZ, [1,3,5,7]) (quadratic_form__theta.py:theta_series)', () => {
    const Q = DiagonalQuadraticForm(ZZ, [1n, 3n, 5n, 7n]);
    expect(Q.theta_series(10n)).toEqual([1n, 2n, 0n, 2n, 6n, 2n, 4n, 6n, 8n, 14n]);
    expect(Q.theta_series(25n)).toEqual([
      1n,
      2n,
      0n,
      2n,
      6n,
      2n,
      4n,
      6n,
      8n,
      14n,
      4n,
      12n,
      18n,
      12n,
      12n,
      8n,
      34n,
      12n,
      8n,
      32n,
      10n,
      28n,
      0n,
      16n,
      44n,
    ]);
  });

  test('representation_number_list (quadratic_form__ternary_Tornaria.py)', () => {
    const Q = DiagonalQuadraticForm(ZZ, [1n, 1n, 1n, 1n, 1n, 1n, 1n, 1n]);
    expect(Q.representation_number_list(10n)).toEqual([
      1n,
      16n,
      112n,
      448n,
      1136n,
      2016n,
      3136n,
      5504n,
      9328n,
      12112n,
    ]);
  });

  test('sum of four squares: r_4(n) = 8 * sum of divisors not divisible by 4', () => {
    // quadratic_form__theta.py:theta_by_pari doctest, with Prec = 100.
    const Q = DiagonalQuadraticForm(ZZ, [1n, 1n, 1n, 1n]);
    const computed = Q.theta_by_pari(100n);
    const exact: bigint[] = [1n];
    for (let i = 1n; i < 100n; i++) {
      let s = 0n;
      for (let d = 1n; d <= i; d++) {
        if (i % d === 0n && d % 4n !== 0n) s += d;
      }
      exact.push(8n * s);
    }
    expect(computed).toEqual(exact);
  });
});

describe('honest stubs', () => {
  test('quadratic_form_from_invariants', () => {
    expect(() => quadratic_form_from_invariants(QQ, 2n, -15n, [3n, 5n], 1n)).toThrow(
      /hilbert_symbol_negative_at_S/
    );
  });

  test('unimplemented QuadraticForm methods name what is missing', () => {
    const Q = DiagonalQuadraticForm(ZZ, [1n, 1n, 1n]);
    expect(() => Q.local_normal_form(2n)).toThrow(NotImplementedError);
    expect(() => Q.jordan_blocks_by_scale_and_unimodular(2n)).toThrow(/local_normal_form/);
    expect(() => Q.local_genus_symbol(2n)).toThrow(/genera.genus/);
    expect(() => Q.global_genus_symbol()).toThrow(/genera.genus/);
    expect(() => Q.CS_genus_symbol_list()).toThrow(/genera.genus/);
    expect(() => Q.is_globally_equivalent_to(Q)).toThrow(/qfisom/);
    expect(() => Q.is_locally_equivalent_to(Q)).toThrow(/genus symbols/);
    expect(() => Q.solve()).toThrow(/qfsolve/);
    expect(() => Q.automorphisms()).toThrow(/qfauto/);
    expect(() => Q.number_of_automorphisms()).toThrow(/qfauto/);
    expect(() => Q.siegel_product(1n)).toThrow(/local density/);
    expect(() => Q.local_density(2n, 1n)).toThrow(/count_local_2/);
    expect(() => Q.conway_mass()).toThrow(/mass/);
    expect(() => Q.polynomial()).toThrow(/multivariate polynomial/);
    expect(() => QuadraticForm.from_polynomial(null)).toThrow(/multivariate polynomial/);
  });
});
