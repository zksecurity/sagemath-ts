/**
 * Tests for `sage/quadratic_forms/ternary_qf` (and the `ternary.pyx` helpers).
 *
 * Every expected value below is either copied verbatim from a SageMath
 * doctest in `reference/sage/src/sage/quadratic_forms/` (the file:line is
 * quoted) or was produced by running the installed `sage` binary.
 */

import { describe, expect, test } from 'bun:test';
import { NotImplementedError } from '../errors.js';
import { Rational } from '../rings/rational.js';
import { QuadraticForm } from './quadratic_form.js';
import {
  TernaryQF,
  _basic_lemma,
  _basic_lemma_vec,
  _find_zeros_mod_p_2,
  _find_zeros_mod_p_odd,
  _reduced_ternary_form_eisenstein_without_matrix,
  evaluate,
  extend,
  find_a_ternary_qf_by_level_disc,
  find_all_ternary_qf_by_level_disc,
  primitivize,
  red_mfact,
} from './ternary_qf.js';

const R = (a: bigint, b: bigint) => new Rational(a, b);

/** Reduce each component of a 3-tuple into `[0, p)`. */
function normalizeModP(p: bigint) {
  return (v: [bigint, bigint, bigint]): [bigint, bigint, bigint] =>
    v.map((x) => ((x % p) + p) % p) as [bigint, bigint, bigint];
}

/** Sort a list of 3-tuples lexicographically (so that set comparisons work). */
function sortTuples(v: [bigint, bigint, bigint][]): string[] {
  return v.map((t) => t.join(',')).sort();
}

describe('TernaryQF basics', () => {
  test('_repr_ (ternary_qf.py:169)', () => {
    const Q = new TernaryQF([1n, 2n, 3n, 4n, 5n, 6n]);
    expect(Q.toString()).toBe(
      'Ternary quadratic form with integer coefficients:\n[1 2 3]\n[4 5 6]'
    );
    expect(new TernaryQF([0n, 0n, 0n, 0n, 0n, 0n]).toString()).toBe(
      'Ternary quadratic form with integer coefficients:\n[0 0 0]\n[0 0 0]'
    );
  });

  test('wrong number of coefficients (ternary_qf.py:95)', () => {
    expect(() => new TernaryQF([1n, 2n, 3n])).toThrow(
      'Ternary quadratic form must be given by a list of six coefficients'
    );
  });

  test('coefficients / coefficient (ternary_qf.py:101,130)', () => {
    const Q = new TernaryQF([1n, 2n, 3n, 4n, 5n, 6n]);
    expect(Q.coefficients()).toEqual([1n, 2n, 3n, 4n, 5n, 6n]);
    expect(Q.coefficient(2)).toBe(3n);
    expect(Q.coefficient(5)).toBe(6n);
  });

  test('__eq__ (ternary_qf.py:567)', () => {
    const Q = new TernaryQF([1n, 2n, 3n, 1n, 2n, 3n]);
    expect(Q.equals(Q)).toBe(true);
    expect(Q.equals(new TernaryQF([1n, 2n, 3n, 1n, 2n, 2n]))).toBe(false);
  });

  test('__neg__ (ternary_qf.py:398)', () => {
    const Q = new TernaryQF([1n, 1n, 2n, -2n, 0n, -1n]);
    expect(Q.neg().toString()).toBe(
      'Ternary quadratic form with integer coefficients:\n[-1 -1 -2]\n[2 0 1]'
    );
    const Z = new TernaryQF([0n, 0n, 0n, 0n, 0n, 0n]);
    expect(Z.equals(Z.neg())).toBe(true);
  });

  test('matrix and disc (ternary_qf.py:266,298)', () => {
    const Q = new TernaryQF([1n, 1n, 2n, 0n, -1n, 4n]);
    expect(
      Q.matrix()
        .rows()
        .map((r) => r.map((x) => x.numerator))
    ).toEqual([
      [2n, 4n, -1n],
      [4n, 2n, 0n],
      [-1n, 0n, 4n],
    ]);
    expect(Q.disc()).toBe(-25n);
    expect(Q.__call__([1n, 2n, 3n])).toBe(28n);
  });

  test('quadratic_form (ternary_qf.py:247)', () => {
    const Q = new TernaryQF([1n, 2n, 3n, 1n, 1n, 1n]);
    const QF1 = Q.quadratic_form();
    expect(QF1.toString()).toBe(
      'Quadratic form in 3 variables over Integer Ring with coefficients: \n' +
        '[ 1 1 1 ]\n[ * 2 1 ]\n[ * * 3 ]'
    );
    expect(QF1.equals(new QuadraticForm(QF1.base_ring(), 3n, [1n, 1n, 1n, 2n, 1n, 3n]))).toBe(true);
  });
});

describe('TernaryQF evaluation', () => {
  test('on a 3x3 matrix (ternary_qf.py:190)', () => {
    const Q = new TernaryQF([1n, 2n, 3n, 4n, 5n, 6n]);
    const A = [
      [1n, -7n, 1n],
      [0n, -2n, 1n],
      [0n, -1n, 0n],
    ];
    expect((Q.__call__(A) as TernaryQF).toString()).toBe(
      'Ternary quadratic form with integer coefficients:\n[1 187 9]\n[-85 8 -31]'
    );
  });

  test('on a 3x2 matrix gives a QuadraticForm (ternary_qf.py:216)', () => {
    const Q = new TernaryQF([1n, 1n, 1n, -1n, -2n, -3n]);
    expect(Q.__call__([1n, 1n, 1n])).toBe(-3n);
    const M = [
      [358n, 6n],
      [2n, 0n],
      [0n, 4n],
    ];
    expect((Q.__call__(M) as QuadraticForm).toString()).toBe(
      'Quadratic form in 2 variables over Integer Ring with coefficients: \n[ 126020 1388 ]\n[ * 4 ]'
    );
  });

  test('on a 3x3 matrix, second doctest (ternary_qf.py:224)', () => {
    const Q = new TernaryQF([1n, 1n, 1n, -1n, -2n, -3n]);
    const M = [
      [1n, 3n, 0n],
      [-1n, 4n, 2n],
      [1n, -1n, -1n],
    ];
    expect((Q.__call__(M) as TernaryQF).toString()).toBe(
      'Ternary quadratic form with integer coefficients:\n[5 0 7]\n[12 -13 -16]'
    );
  });

  test('evaluate helper (ternary.pyx:491)', () => {
    expect(evaluate(1n, 2n, 3n, -1n, 0n, 0n, [1n, -1n, 19n])).toBe(1105n);
  });
});

describe('TernaryQF definiteness', () => {
  test('is_definite (ternary_qf.py:315)', () => {
    const Q = new TernaryQF([10n, 10n, 1n, -1n, 2n, 3n]);
    expect(Q.is_definite()).toBe(true);
    expect(Q.neg().is_definite()).toBe(true);
    expect(new TernaryQF([1n, 1n, 2n, -3n, 0n, -1n]).is_definite()).toBe(false);
  });

  test('is_positive_definite (ternary_qf.py:341)', () => {
    const Q = new TernaryQF([10n, 10n, 1n, -1n, 2n, 3n]);
    expect(Q.is_positive_definite()).toBe(true);
    expect(Q.neg().is_positive_definite()).toBe(false);
    expect(new TernaryQF([1n, 1n, 0n, 0n, 0n, 0n]).is_positive_definite()).toBe(false);
    expect(new TernaryQF([1n, 1n, 1n, -1n, -2n, -3n]).is_positive_definite()).toBe(false);
  });

  test('is_negative_definite (ternary_qf.py:372)', () => {
    expect(new TernaryQF([-8n, -9n, -10n, 1n, 9n, -3n]).is_negative_definite()).toBe(true);
    const Q = new TernaryQF([-4n, -1n, 6n, -5n, 1n, -5n]);
    expect(Q.__call__([0n, 0n, 1n])).toBe(6n);
    expect(Q.is_negative_definite()).toBe(false);
  });

  test('definiteness agrees with the general QuadraticForm code', () => {
    for (const c of [
      [10n, 10n, 1n, -1n, 2n, 3n],
      [-8n, -9n, -10n, 1n, 9n, -3n],
      [1n, 1n, 2n, -3n, 0n, -1n],
      [1n, 1n, 1n, 0n, 0n, 0n],
    ]) {
      const T = new TernaryQF(c);
      const Q = T.quadratic_form();
      expect(T.is_positive_definite()).toBe(Q.is_positive_definite());
      expect(T.is_negative_definite()).toBe(Q.is_negative_definite());
    }
  });
});

describe('TernaryQF content, adjoint, level', () => {
  test('content / is_primitive / primitive (ternary_qf.py:610,418,440)', () => {
    expect(new TernaryQF([1n, 2n, 3n, 4n, 5n, 6n]).is_primitive()).toBe(true);
    expect(new TernaryQF([1n, 2n, 3n, 4n, 5n, 6n]).content()).toBe(1n);
    const Q = new TernaryQF([10n, 10n, 10n, 5n, 5n, 5n]);
    expect(Q.content()).toBe(5n);
    expect(Q.is_primitive()).toBe(false);
    expect(Q.primitive().toString()).toBe(
      'Ternary quadratic form with integer coefficients:\n[2 2 2]\n[1 1 1]'
    );
    expect(new TernaryQF([1n, 1n, 2n, 0n, 0n, 0n]).content()).toBe(1n);
    expect(new TernaryQF([2n, 4n, 6n, 0n, 0n, 0n]).content()).toBe(2n);
    expect(
      (new TernaryQF([2n, 4n, 6n, 0n, 0n, 0n]).scale_by_factor(100n) as TernaryQF).content()
    ).toBe(200n);
  });

  test('scale_by_factor (ternary_qf.py:465)', () => {
    const Q = new TernaryQF([2n, 2n, 4n, 0n, -2n, 8n]);
    expect((Q.scale_by_factor(5n) as TernaryQF).toString()).toBe(
      'Ternary quadratic form with integer coefficients:\n[10 10 20]\n[0 -10 40]'
    );
    expect((Q.scale_by_factor(R(1n, 2n)) as TernaryQF).toString()).toBe(
      'Ternary quadratic form with integer coefficients:\n[1 1 2]\n[0 -1 4]'
    );
    expect((Q.scale_by_factor(R(1n, 3n)) as QuadraticForm).toString()).toBe(
      'Quadratic form in 3 variables over Rational Field with coefficients: \n' +
        '[ 2/3 8/3 -2/3 ]\n[ * 2/3 0 ]\n[ * * 4/3 ]'
    );
  });

  test('adjoint (ternary_qf.py:584)', () => {
    const Q = new TernaryQF([1n, 1n, 17n, 0n, 0n, 1n]);
    expect(Q.adjoint().toString()).toBe(
      'Ternary quadratic form with integer coefficients:\n[68 68 3]\n[0 0 -68]'
    );
    // Q.adjoint().matrix() == 2*Q.matrix().adjoint_classical()
    const A = Q.adjoint().matrix();
    const H = Q.matrix();
    // classical adjoint of a 3x3 matrix: det * H^{-1}; check A == 2 * adj(H)
    const det = 2n * Q.disc();
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        // (adj H)_{ij} = cofactor_{ji}
        const rows = [0, 1, 2].filter((x) => x !== j);
        const cols = [0, 1, 2].filter((x) => x !== i);
        const m = H.get(rows[0]!, cols[0]!)
          .mul(H.get(rows[1]!, cols[1]!))
          .sub(H.get(rows[0]!, cols[1]!).mul(H.get(rows[1]!, cols[0]!)));
        const sign = (i + j) % 2 === 0 ? 1n : -1n;
        expect(A.get(i, j).eq(m.mul(sign).mul(2n))).toBe(true);
      }
    }
    expect(det).not.toBe(0n);
  });

  test('reciprocal (ternary_qf.py:506)', () => {
    const Q = new TernaryQF([2n, 2n, 14n, 0n, 0n, 0n]);
    expect(Q.reciprocal().toString()).toBe(
      'Ternary quadratic form with integer coefficients:\n[14 14 2]\n[0 0 0]'
    );
    expect(Q.content()).toBe(2n);
    expect(Q.reciprocal().content()).toBe(2n);
    expect(Q.adjoint().content()).toBe(16n);
  });

  test('reciprocal_reduced (ternary_qf.py:529)', () => {
    const Q = new TernaryQF([1n, 1n, 3n, 0n, -1n, 0n]);
    const Qrr = Q.reciprocal_reduced();
    expect(Qrr.toString()).toBe(
      'Ternary quadratic form with integer coefficients:\n[4 11 12]\n[0 -4 0]'
    );
    expect(Q.is_eisenstein_reduced()).toBe(true);
    expect(Q.reciprocal().reduced_form_eisenstein(false).equals(Qrr)).toBe(true);
  });

  test('divisor / omega / delta / level (ternary_qf.py:548,627,642,662)', () => {
    expect(new TernaryQF([1n, 1n, 17n, 0n, 0n, 0n]).divisor()).toBe(4n);
    const Q = new TernaryQF([4n, 11n, 12n, 0n, -4n, 0n]);
    expect(Q.omega()).toBe(176n);
    expect(Q.primitive().adjoint().content()).toBe(176n);
    const Q2 = new TernaryQF([1n, 2n, 2n, -1n, 0n, -1n]);
    expect(Q2.delta()).toBe(208n);
    expect(Q2.adjoint().omega()).toBe(208n);
    const Q3 = new TernaryQF([1n, -1n, 1n, 0n, 0n, 0n]);
    expect(Q3.delta()).toBe(4n);
    expect(Q3.omega()).toBe(4n);
    expect(Q2.level()).toBe(52n);
    expect((4n * Q2.disc()) / Q2.divisor()).toBe(52n);
  });
});

describe('Eisenstein reduction', () => {
  test('is_eisenstein_reduced (ternary_qf.py:676)', () => {
    expect(new TernaryQF([1n, 1n, 1n, 0n, 0n, 0n]).is_eisenstein_reduced()).toBe(true);
    expect(new TernaryQF([34n, 14n, 44n, 12n, 25n, -22n]).is_eisenstein_reduced()).toBe(false);
  });

  test('red_mfact (ternary.pyx:22)', () => {
    expect(red_mfact(0n, 3n)).toBe(0n);
    expect(red_mfact(-5n, 100n)).toBe(9n);
  });

  test('reduced_form_eisenstein with matrix (ternary_qf.py:748)', () => {
    const Q = new TernaryQF([293n, 315n, 756n, 908n, 929n, 522n]);
    const [Qr, m] = Q.reduced_form_eisenstein();
    expect(Qr.toString()).toBe(
      'Ternary quadratic form with integer coefficients:\n[1 2 2]\n[-1 0 -1]'
    );
    expect(Qr.is_eisenstein_reduced()).toBe(true);
    expect(m).toEqual([
      [-54n, 137n, -38n],
      [-23n, 58n, -16n],
      [47n, -119n, 33n],
    ]);
    // m.det() == 1
    const det =
      m[0]![0]! * (m[1]![1]! * m[2]![2]! - m[1]![2]! * m[2]![1]!) -
      m[0]![1]! * (m[1]![0]! * m[2]![2]! - m[1]![2]! * m[2]![0]!) +
      m[0]![2]! * (m[1]![0]! * m[2]![1]! - m[1]![1]! * m[2]![0]!);
    expect(det).toBe(1n);
    expect((Q.__call__(m) as TernaryQF).equals(Qr)).toBe(true);
  });

  test('reduced_form_eisenstein without matrix (ternary_qf.py:770)', () => {
    expect(
      new TernaryQF([12n, 36n, 3n, 14n, -7n, -19n]).reduced_form_eisenstein(false).toString()
    ).toBe('Ternary quadratic form with integer coefficients:\n[3 8 20]\n[3 2 1]');
    expect(
      _reduced_ternary_form_eisenstein_without_matrix(293n, 315n, 756n, 908n, 929n, 522n)
    ).toEqual([1n, 2n, 2n, -1n, 0n, -1n]);
  });

  test('reduction is an equivalence: Q(M) == Qr and Qr is reduced', () => {
    const cases: bigint[][] = [
      [293n, 315n, 756n, 908n, 929n, 522n],
      [12n, 36n, 3n, 14n, -7n, -19n],
      [1n, 1n, 3n, 0n, -1n, 0n],
      [4n, 7n, 8n, -4n, -1n, -3n],
      [26n, 42n, 53n, -36n, -17n, -3n],
      [9n, 12n, 30n, -26n, -28n, 20n],
    ];
    for (const c of cases) {
      const Q = new TernaryQF(c);
      const [Qr, M] = Q.reduced_form_eisenstein();
      expect(Qr.is_eisenstein_reduced()).toBe(true);
      expect((Q.__call__(M) as TernaryQF).equals(Qr)).toBe(true);
      expect(Qr.disc()).toBe(Q.disc());
    }
  });
});

describe('ternary.pyx helpers', () => {
  test('extend (ternary.pyx:933)', () => {
    expect(extend([6n, 4n, 12n])).toEqual([6n, 1n, 0n, 4n, 1n, 0n, 12n, 0n, 1n]);
    expect(extend([-12n, 20n, 30n])).toEqual([-12n, 1n, 0n, 20n, -2n, 1n, 30n, 0n, -7n]);
    expect(extend([0n, 0n, 1n])).toEqual([0n, 1n, 0n, 0n, 0n, 1n, 1n, 0n, 0n]);
    expect(extend([0n, 0n, -2n])).toEqual([0n, 0n, 1n, 0n, 1n, 0n, -2n, 0n, 0n]);
  });

  test('primitivize (ternary.pyx:470)', () => {
    expect(primitivize(12n, 13n, 14n, 5n)).toEqual([3n, 2n, 1n]);
    expect(primitivize(12n, 13n, 15n, 5n)).toEqual([4n, 1n, 0n]);
  });

  test('_basic_lemma / _basic_lemma_vec (ternary.pyx:1124,1089)', () => {
    expect(_basic_lemma(5n, 2n, 3n, -1n, 0n, 0n, 5n)).toBe(2n);
    expect(_basic_lemma_vec(5n, 2n, 3n, -1n, 0n, 0n, 5n)).toEqual([0n, 1n, 0n]);
    expect(new TernaryQF([3n, 3n, 3n, -2n, 0n, -1n]).basic_lemma(3n)).toBe(4n);
  });

  test('_find_zeros_mod_p_odd from a fixed zero (ternary.pyx:575)', () => {
    // sage: _find_zeros_mod_p_odd(1, 2, 2, -1, 0, 0, 1009, (817, 974, 1))
    const z = _find_zeros_mod_p_odd(1n, 2n, 2n, -1n, 0n, 0n, 1009n, [817n, 974n, 1n]);
    expect(z.length).toBe(1010);
    const sorted = z.slice().sort((a, b) => Number(a[0] - b[0] || a[1] - b[1] || a[2] - b[2]));
    expect(sorted[0]).toEqual([0n, 32n, 1n]);
    // Q((0, 32, 1)) == 2018 == 2 * 1009
    expect(evaluate(1n, 2n, 2n, -1n, 0n, 0n, [0n, 32n, 1n])).toBe(2018n);
    for (const w of z) {
      expect(evaluate(1n, 2n, 2n, -1n, 0n, 0n, w) % 1009n).toBe(0n);
    }
    // Every zero is distinct as a point of P^2(F_1009).
    const norm = z.map(normalizeModP(1009n)).map((w) => w.join(','));
    expect(new Set(norm).size).toBe(1010);
  });

  test('_find_zeros_mod_p_2 (ternary.pyx:508)', () => {
    expect(_find_zeros_mod_p_2(1n, 2n, 2n, -1n, 0n, 0n)).toEqual([
      [0n, 1n, 0n],
      [0n, 0n, 1n],
      [1n, 1n, 1n],
    ]);
  });
});

describe('zeros mod p and p-neighbours', () => {
  test('find_zeros_mod_p (ternary_qf.py:821)', () => {
    const Q = new TernaryQF([4n, 7n, 8n, -4n, -1n, -3n]);
    expect(Q.is_positive_definite()).toBe(true);
    expect(Q.disc()).toBe(741n); // 3 * 13 * 19
    expect(Q.find_zeros_mod_p(2n)).toEqual([
      [1n, 0n, 0n],
      [1n, 1n, 0n],
      [0n, 0n, 1n],
    ]);
    const zeros17 = Q.find_zeros_mod_p(17n);
    expect(zeros17.length).toBe(18);
    for (const v of zeros17) {
      expect((((Q.__call__(v) as bigint) % 17n) + 17n) % 17n).toBe(0n);
    }
    // The full set, up to the choice of representatives mod 17 (the order and
    // the exact representatives depend on the starting zero; see the deviation
    // note on primitive_zero_mod_p).  From sage:
    //   sorted(TernaryQF([4,7,8,-4,-1,-3]).find_zeros_mod_p(17))
    expect(sortTuples(zeros17.map(normalizeModP(17n)))).toEqual(
      sortTuples([
        [0n, 7n, 1n],
        [0n, 13n, 1n],
        [1n, 6n, 1n],
        [1n, 12n, 1n],
        [2n, 1n, 0n],
        [3n, 1n, 0n],
        [4n, 0n, 1n],
        [4n, 12n, 1n],
        [6n, 4n, 1n],
        [8n, 6n, 1n],
        [8n, 15n, 1n],
        [9n, 0n, 1n],
        [9n, 2n, 1n],
        [10n, 4n, 1n],
        [10n, 13n, 1n],
        [12n, 15n, 1n],
        [14n, 2n, 1n],
        [14n, 7n, 1n],
      ])
    );
  });

  test('find_p_neighbor_from_vec (ternary_qf.py:848)', () => {
    const Q = new TernaryQF([1n, 3n, 3n, -2n, 0n, -1n]);
    expect(Q.disc()).toBe(29n);
    const v: [bigint, bigint, bigint] = [9n, 7n, 1n];
    expect(sortTuples(Q.find_zeros_mod_p(11n))).toContain('9,7,1');
    const [Q11, M] = Q.find_p_neighbor_from_vec(11n, v, true);
    expect(Q11.toString()).toBe(
      'Ternary quadratic form with integer coefficients:\n[1 2 4]\n[-1 -1 0]'
    );
    expect(M.rows().map((r) => r.map((x) => x.toString()))).toEqual([
      ['-1', '-5/11', '7/11'],
      ['0', '-10/11', '3/11'],
      ['0', '-3/11', '13/11'],
    ]);
    expect((Q.__call__(M) as TernaryQF).equals(Q11)).toBe(true);
  });

  test('find_p_neighbor_from_vec with (0, 0, 1) (ternary_qf.py:884)', () => {
    const Q = new TernaryQF([1n, 3n, 3n, -2n, 0n, -1n]);
    expect(Q.find_p_neighbor_from_vec(3n, [0n, 0n, 1n]).toString()).toBe(
      'Ternary quadratic form with integer coefficients:\n[1 3 3]\n[-2 0 -1]'
    );
  });

  test('find_p_neighbors (ternary_qf.py:899)', () => {
    const Q0 = new TernaryQF([1n, 3n, 3n, -2n, 0n, -1n]);
    const neig = Q0.find_p_neighbors(5n);
    expect(neig.length).toBe(6);
    const Q1 = new TernaryQF([1n, 1n, 10n, 1n, 1n, 1n]);
    const Q2 = new TernaryQF([1n, 2n, 4n, -1n, -1n, 0n]);
    expect(neig.filter((q) => q.equals(Q0)).length).toBe(2);
    expect(neig.filter((q) => q.equals(Q1)).length).toBe(1);
    expect(neig.filter((q) => q.equals(Q2)).length).toBe(3);
    // Neighbours share the discriminant of the original form.
    for (const q of neig) {
      expect(q.disc()).toBe(Q0.disc());
    }
  });
});

describe('genus characters', () => {
  test('xi (ternary_qf.py:938)', () => {
    const Q1 = new TernaryQF([26n, 42n, 53n, -36n, -17n, -3n]);
    expect(Q1.omega()).toBe(3n);
    expect(Q1.xi(3n)).toBe(-1n);
    // find_zeros_mod_p(2) is deterministic in SageMath too, so index 1 matches.
    const Q2 = Q1.find_p_neighbors(2n)[1]!;
    expect(Q2.xi(3n)).toBe(-1n);
  });

  test('xi_rec (ternary_qf.py:975)', () => {
    const Q1 = new TernaryQF([1n, 1n, 7n, 0n, 0n, 0n]);
    expect(Q1.delta()).toBe(28n);
    expect(Q1.xi_rec(7n)).toBe(1n);
    // Order-independent form of the doctest: every 3-neighbour has xi_rec 1.
    for (const Q2 of Q1.find_p_neighbors(3n)) {
      expect(Q2.xi_rec(7n)).toBe(1n);
    }
  });

  test('xi rejects invalid characters', () => {
    expect(() => new TernaryQF([1n, 1n, 1n, 0n, 0n, 0n]).xi(5n)).toThrow('not a valid character');
  });
});

describe('symmetries and spin norms', () => {
  test('symmetry (ternary_qf.py:990)', () => {
    const Q = new TernaryQF([4n, 5n, 8n, 5n, 2n, 2n]);
    const M = Q.symmetry([1n, 1n, 1n]);
    expect(M.rows().map((r) => r.map((x) => x.toString()))).toEqual([
      ['7/13', '-17/26', '-23/26'],
      ['-6/13', '9/26', '-23/26'],
      ['-6/13', '-17/26', '3/26'],
    ]);
    // M * v == -v
    const v = [1n, 1n, 1n];
    for (let i = 0; i < 3; i++) {
      let s = Rational.zero();
      for (let j = 0; j < 3; j++) s = s.add(M.get(i, j).mul(v[j]!));
      expect(s.eq(-1n)).toBe(true);
    }
  });

  test('automorphism_symmetries and spin norm (ternary_qf.py:1026,1066)', () => {
    const Q = new TernaryQF([9n, 12n, 30n, -26n, -28n, 20n]);
    const A = [
      [9n, 10n, -10n],
      [-6n, -7n, 6n],
      [2n, 2n, -3n],
    ];
    expect((Q.__call__(A) as TernaryQF).equals(Q)).toBe(true);
    const [v1, v2] = Q.automorphism_symmetries(A);
    expect(v1!.map((x) => x.toString())).toEqual(['8', '-6', '2']);
    expect(v2!.map((x) => x.toString())).toEqual(['1', '-5/4', '-1/4']);
    expect(Q.automorphism_spin_norm(A)).toBe(7n);
    expect(
      Q.automorphism_symmetries([
        [1n, 0n, 0n],
        [0n, 1n, 0n],
        [0n, 0n, 1n],
      ])
    ).toEqual([]);
    expect(
      Q.automorphism_spin_norm([
        [1n, 0n, 0n],
        [0n, 1n, 0n],
        [0n, 0n, 1n],
      ])
    ).toBe(1n);
  });
});

describe('honest stubs', () => {
  test('automorphisms / number_of_automorphisms', () => {
    const Q = new TernaryQF([1n, 1n, 1n, 0n, 0n, 0n]);
    expect(() => Q.automorphisms()).toThrow(NotImplementedError);
    expect(() => Q.number_of_automorphisms()).toThrow(NotImplementedError);
    expect(() => Q.polynomial()).toThrow(/multivariate polynomial/);
  });

  test('find_(a|all)_ternary_qf_by_level_disc', () => {
    expect(() => find_all_ternary_qf_by_level_disc(44n, 11n)).toThrow(NotImplementedError);
    expect(() => find_a_ternary_qf_by_level_disc(44n, 11n)).toThrow(NotImplementedError);
  });
});
