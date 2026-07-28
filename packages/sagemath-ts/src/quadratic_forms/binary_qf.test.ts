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
      // binary_qf.py:905 doctest: the matrix is [[1, -1], [0, 1]]
      expect(M).toEqual([
        [1n, -1n],
        [0n, 1n],
      ]);
    });

    // binary_qf.py:867-870
    test('15x^2 + 15y^2 is already reduced', () => {
      const Q = new BinaryQF(15n, 0n, 15n);
      expect(Q.is_reduced()).toBe(true);
      expect(Q.reduced_form().toTuple()).toEqual([15n, 0n, 15n]);
    });

    // binary_qf.py:860 doctest: 33x^2 + 11xy + 5y^2 -> 5x^2 - xy + 27y^2
    test('33x^2 + 11xy + 5y^2 reduces to 5x^2 - xy + 27y^2', () => {
      const R = new BinaryQF(33n, 11n, 5n).reduced_form();
      expect(R.toTuple()).toEqual([5n, -1n, 27n]);
      expect(R.is_reduced()).toBe(true);
    });

    // C15: b == -a must be normalized into (-a, a]
    test('b = -a is normalized (2x^2 - 2xy + 5y^2 -> 2x^2 + 2xy + 5y^2)', () => {
      const Q = new BinaryQF(2n, -2n, 5n);
      expect(Q.is_reduced()).toBe(false);
      const R = Q.reduced_form();
      expect(R.toTuple()).toEqual([2n, 2n, 5n]);
      expect(R.is_reduced()).toBe(true);
    });

    // binary_qf.py:914-916 (issue 34229)
    test('negative definite BinaryQF(-225, -743, -743) reduces to a reduced form', () => {
      const Q = new BinaryQF(-225n, -743n, -743n);
      const R = Q.reduced_form();
      expect(R.is_reduced()).toBe(true);
      expect(R.discriminant()).toBe(Q.discriminant());
      const [R2, M] = Q.reduced_form({ transformation: true });
      expect(R2.equals(R)).toBe(true);
      expect(Q.matrix_action_right(M).equals(R)).toBe(true);
    });

    // H74: the negative definite branch must conjugate the base change
    test('negative definite transformation satisfies f * M == g', () => {
      const Q = new BinaryQF(-20n, -22n, -16n);
      const [R, M] = Q.reduced_form({ transformation: true });
      expect(Q.matrix_action_right(M).equals(R)).toBe(true);
      expect(M[0][0] * M[1][1] - M[0][1] * M[1][0]).toBe(1n);
    });

    describe('indefinite forms', () => {
      // binary_qf.py:876-884
      test('x^2 - 3y^2 -> x^2 + 2xy - 2y^2', () => {
        const Q = new BinaryQF(1n, 0n, -3n);
        expect(Q.is_reduced()).toBe(false);
        const R = Q.reduced_form();
        expect(R.toTuple()).toEqual([1n, 2n, -2n]);
        expect(R.is_reduced()).toBe(true);
      });

      // binary_qf.py:886-888 (square discriminant)
      test('x^2 - y^2 -> x^2 + 2xy', () => {
        expect(new BinaryQF(1n, 0n, -1n).reduced_form().toTuple()).toEqual([1n, 2n, 0n]);
      });

      // binary_qf.py:890-894
      test('BinaryQF(1, 9, 4) -> 4x^2 + 7xy - y^2 with M = [[0, -1], [1, 2]]', () => {
        const Q = new BinaryQF(1n, 9n, 4n);
        const [R, M] = Q.reduced_form({ transformation: true });
        expect(R.toTuple()).toEqual([4n, 7n, -1n]);
        expect(M).toEqual([
          [0n, -1n],
          [1n, 2n],
        ]);
        expect(Q.matrix_action_right(M).equals(R)).toBe(true);
      });

      // binary_qf.py:895-899
      test('BinaryQF(3, 7, -2) is already reduced', () => {
        const Q = new BinaryQF(3n, 7n, -2n);
        const [R, M] = Q.reduced_form({ transformation: true });
        expect(R.toTuple()).toEqual([3n, 7n, -2n]);
        expect(M).toEqual([
          [1n, 0n],
          [0n, 1n],
        ]);
      });

      // binary_qf.py:900-904
      test('BinaryQF(-6, 6, -1) -> -x^2 + 2xy + 2y^2 with M = [[0, -1], [1, -4]]', () => {
        const Q = new BinaryQF(-6n, 6n, -1n);
        const [R, M] = Q.reduced_form({ transformation: true });
        expect(R.toTuple()).toEqual([-1n, 2n, 2n]);
        expect(M).toEqual([
          [0n, -1n],
          [1n, -4n],
        ]);
        expect(Q.matrix_action_right(M).equals(R)).toBe(true);
      });

      // binary_qf.py:780-788 (_reduce_indef doctest)
      test('BinaryQF(-1, 0, 3) -> -x^2 + 2xy + 2y^2 with M = [[-1, 1], [0, -1]]', () => {
        const Q = new BinaryQF(-1n, 0n, 3n);
        const [R, M] = Q.reduced_form({ transformation: true });
        expect(R.toTuple()).toEqual([-1n, 2n, 2n]);
        expect(M).toEqual([
          [-1n, 1n],
          [0n, -1n],
        ]);
        expect(Q.matrix_action_right(M).equals(R)).toBe(true);
      });

      // binary_qf.py:790-792 (_reduce_indef doctest, c == 0 branch)
      test('BinaryQF(0, 5, 24) reduces with a valid transformation', () => {
        const Q = new BinaryQF(0n, 5n, 24n);
        const [R, M] = Q.reduced_form({ transformation: true });
        expect(R.is_reduced()).toBe(true);
        expect(Q.matrix_action_right(M).equals(R)).toBe(true);
      });

      // H72: square discriminants used to loop forever
      test('square discriminant reductions terminate', () => {
        expect(new BinaryQF(-3n, 4n, 0n).reduced_form().toTuple()).toEqual([1n, 4n, 0n]);
        expect(new BinaryQF(6n, 5n, -11n).reduced_form().toTuple()).toEqual([6n, 17n, 0n]);
      });
    });
  });

  describe('cycle', () => {
    const tuples = (forms: BinaryQF[]) => forms.map((f) => f.toTuple());

    // binary_qf.py:1114-1123
    test('cycle of 14x^2 + 17xy - 2y^2', () => {
      const Q = new BinaryQF(14n, 17n, -2n);
      expect(tuples(Q.cycle())).toEqual([
        [14n, 17n, -2n],
        [2n, 19n, -5n],
        [5n, 11n, -14n],
      ]);
      expect(tuples(Q.cycle({ proper: true }))).toEqual([
        [14n, 17n, -2n],
        [-2n, 19n, 5n],
        [5n, 11n, -14n],
        [-14n, 17n, 2n],
        [2n, 19n, -5n],
        [-5n, 11n, 14n],
      ]);
    });

    // binary_qf.py:1155-1180 (Example 6.10.7 of [BUVO2007])
    test('cycle of x^2 + 8xy - 3y^2', () => {
      const Q = new BinaryQF(1n, 8n, -3n);
      expect(tuples(Q.cycle())).toEqual([
        [1n, 8n, -3n],
        [3n, 4n, -5n],
        [5n, 6n, -2n],
        [2n, 6n, -5n],
        [5n, 4n, -3n],
        [3n, 8n, -1n],
      ]);
      expect(tuples(Q.cycle({ proper: true }))).toEqual([
        [1n, 8n, -3n],
        [-3n, 4n, 5n],
        [5n, 6n, -2n],
        [-2n, 6n, 5n],
        [5n, 4n, -3n],
        [-3n, 8n, 1n],
      ]);
    });

    // binary_qf.py:1185-1191 (a negative)
    test('proper cycle of -x^2 + 8xy + 3y^2', () => {
      expect(tuples(new BinaryQF(-1n, 8n, 3n).cycle({ proper: true }))).toEqual([
        [-1n, 8n, 3n],
        [3n, 4n, -5n],
        [-5n, 6n, 2n],
        [2n, 6n, -5n],
        [-5n, 4n, 3n],
        [3n, 8n, -1n],
      ]);
    });

    // binary_qf.py:1119-1121 (issue 28989): odd cycles get doubled
    test('proper cycle of x^2 + xy - y^2 (issue 28989)', () => {
      expect(tuples(new BinaryQF(1n, 1n, -1n).cycle({ proper: true }))).toEqual([
        [1n, 1n, -1n],
        [-1n, 1n, 1n],
      ]);
    });

    // binary_qf.py:1122-1134 (Example 6.10.6)
    test('cycle of x^2 + 7xy - 6y^2', () => {
      expect(tuples(new BinaryQF(1n, 7n, -6n).cycle())).toEqual([
        [1n, 7n, -6n],
        [6n, 5n, -2n],
        [2n, 7n, -3n],
        [3n, 5n, -4n],
        [4n, 3n, -4n],
        [4n, 5n, -3n],
        [3n, 7n, -2n],
        [2n, 5n, -6n],
        [6n, 7n, -1n],
      ]);
      expect(new BinaryQF(1n, 7n, -6n).cycle({ proper: true }).length).toBe(18);
    });

    test('rejects definite and non-reduced forms', () => {
      expect(() => new BinaryQF(1n, 0n, 1n).cycle()).toThrow(/indefinite and reduced/);
      expect(() => new BinaryQF(1n, 9n, 4n).cycle()).toThrow(/indefinite and reduced/);
    });

    test('rejects square discriminants', () => {
      // (1, 2, 0) is reduced with discriminant 4
      expect(() => new BinaryQF(1n, 2n, 0n).cycle()).toThrow(/non-square/);
    });

    test('every form in a cycle is reduced and of the same discriminant', () => {
      for (const f of [
        new BinaryQF(14n, 17n, -2n),
        new BinaryQF(1n, 8n, -3n),
        new BinaryQF(1n, 7n, -6n),
      ]) {
        for (const g of f.cycle({ proper: true })) {
          expect(g.is_reduced()).toBe(true);
          expect(g.discriminant()).toBe(f.discriminant());
        }
      }
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

    // C14: the order-2 class of discriminant -20 must square to the principal form
    test('(2, 2, 3)^2 is the principal form of discriminant -20', () => {
      const Q = new BinaryQF(2n, 2n, 3n);
      expect(Q.discriminant()).toBe(-20n);
      expect(Q.compose(Q).reduced_form().toTuple()).toEqual([1n, 0n, 5n]);
    });

    // binary_qf.py:235-244 doctests for the class group of discriminant -23
    test('class group of discriminant -23', () => {
      const R = BinaryQF_reduced_representatives(-23n, { primitive_only: false });
      expect(R.map((f) => f.toTuple())).toEqual([
        [1n, 1n, 6n],
        [2n, -1n, 3n],
        [2n, 1n, 3n],
      ]);
      expect(R[0]!.compose(R[0]!).toTuple()).toEqual([1n, 1n, 6n]);
      expect(R[1]!.compose(R[1]!).toTuple()).toEqual([4n, 3n, 2n]);
      expect(R[1]!.compose(R[1]!).reduced_form().toTuple()).toEqual([2n, 1n, 3n]);
      expect(R[1]!.compose(R[1]!).compose(R[1]!).reduced_form().toTuple()).toEqual([1n, 1n, 6n]);
    });

    test('composition makes the reduced forms an abelian group', () => {
      for (const D of [-20n, -23n, -47n, -56n, -95n, -104n]) {
        const reps = BinaryQF_reduced_representatives(D, { primitive_only: false });
        const principal = BinaryQF.principal(D).reduced_form();
        const idx = (f: BinaryQF) => reps.findIndex((r) => r.equals(f));
        const compose = (i: number, j: number) => idx(reps[i]!.compose(reps[j]!).reduced_form());
        const e = idx(principal);
        expect(e).toBeGreaterThanOrEqual(0);
        for (let i = 0; i < reps.length; i++) {
          expect(compose(i, e)).toBe(i);
          // inverse
          expect(reps[i]!.compose(reps[i]!.inverse()).reduced_form().equals(principal)).toBe(true);
          for (let j = 0; j < reps.length; j++) {
            expect(compose(i, j)).toBe(compose(j, i));
            expect(compose(i, j)).toBeGreaterThanOrEqual(0);
            for (let k = 0; k < reps.length; k++) {
              expect(compose(compose(i, j), k)).toBe(compose(i, compose(j, k)));
            }
          }
        }
      }
    });

    // binary_qf.py:213-216: the principal form is idempotent
    test('principal form is idempotent for many discriminants', () => {
      for (let D = -200n; D <= 200n; D++) {
        const m = ((D % 4n) + 4n) % 4n;
        if (D === 0n || (m !== 0n && m !== 1n)) continue;
        const Q = BinaryQF.principal(D);
        if (Q.discriminant() > 0n && Q.is_reducible()) continue;
        expect(Q.compose(Q).is_equivalent(Q)).toBe(true);
      }
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

    // binary_qf.py:1308-1319
    test('definite doctests', () => {
      expect(new BinaryQF(4n, -4n, 15n).is_equivalent(new BinaryQF(4n, 4n, 15n))).toBe(true);
      const a = new BinaryQF(33n, 11n, 5n);
      expect(a.is_equivalent(a.reduced_form())).toBe(true);
      expect(a.is_equivalent(new BinaryQF(3n, 4n, 5n))).toBe(false);
    });

    // binary_qf.py:1321-1328
    test('indefinite proper vs improper equivalence', () => {
      const Q1 = new BinaryQF(9n, 8n, -7n);
      const Q2 = new BinaryQF(9n, -8n, -7n);
      expect(Q1.is_equivalent(Q2, { proper: true })).toBe(false);
      expect(Q1.is_equivalent(Q2, { proper: false })).toBe(true);
    });

    // binary_qf.py:1332-1341 (issue 25888)
    test('issue 25888: (3, 4, -2) and (-2, 4, 3) are properly equivalent', () => {
      const Q1 = new BinaryQF(3n, 4n, -2n);
      const Q2 = new BinaryQF(-2n, 4n, 3n);
      expect(Q1.is_equivalent(Q2, { proper: true })).toBe(true);
      expect(Q1.is_equivalent(Q2)).toBe(Q2.is_equivalent(Q1));
      expect(Q1.is_equivalent(Q2, { proper: false })).toBe(Q2.is_equivalent(Q1, { proper: false }));
    });

    // binary_qf.py:1343-1352 (issue 29028) - square discriminant
    test('issue 29028: BinaryQF(0, 2, 0) is equivalent to itself', () => {
      const Q = new BinaryQF(0n, 2n, 0n);
      expect(Q.discriminant()).toBe(4n);
      expect(Q.is_equivalent(Q, { proper: true })).toBe(true);
      expect(Q.is_equivalent(Q, { proper: false })).toBe(true);
    });

    // binary_qf.py:1354-1360 (square discriminant, Conway-Sloane branch)
    test('(0, 4, 2) and (2, 4, 0) are improperly equivalent', () => {
      expect(
        new BinaryQF(0n, 4n, 2n).is_equivalent(new BinaryQF(2n, 4n, 0n), { proper: false })
      ).toBe(true);
    });

    // binary_qf.py:1362-1366 (issue 28989)
    test('issue 28989: (1, 1, -1) and (-1, 1, 1) are properly equivalent', () => {
      expect(
        new BinaryQF(1n, 1n, -1n).is_equivalent(new BinaryQF(-1n, 1n, 1n), { proper: true })
      ).toBe(true);
    });

    test('a form is equivalent to its reduction (random sweep)', () => {
      let seed = 20260727n;
      const rnd = (lo: bigint, hi: bigint) => {
        seed = (seed * 6364136223846793005n + 1442695040888963407n) & ((1n << 64n) - 1n);
        return lo + ((seed >> 17n) % (hi - lo + 1n));
      };
      let checked = 0;
      for (let i = 0; i < 300; i++) {
        const f = new BinaryQF(rnd(-40n, 40n), rnd(-40n, 40n), rnd(-40n, 40n));
        if (f.discriminant() === 0n) continue;
        checked++;
        const g = f.reduced_form();
        expect(g.is_reduced()).toBe(true);
        expect(g.discriminant()).toBe(f.discriminant());
        expect(f.is_equivalent(g)).toBe(true);
        const [g2, M] = f.reduced_form({ transformation: true });
        expect(g2.equals(g)).toBe(true);
        expect(f.matrix_action_right(M).equals(g)).toBe(true);
      }
      expect(checked).toBeGreaterThan(250);
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

  const tuples = (forms: BinaryQF[]) => forms.map((f) => f.toTuple());

  // binary_qf.py:1855-1866
  test('definite doctests', () => {
    expect(tuples(BinaryQF_reduced_representatives(-12n))).toEqual([
      [1n, 0n, 3n],
      [2n, 2n, 2n],
    ]);
    expect(tuples(BinaryQF_reduced_representatives(-16n))).toEqual([
      [1n, 0n, 4n],
      [2n, 0n, 2n],
    ]);
    expect(tuples(BinaryQF_reduced_representatives(-63n))).toEqual([
      [1n, 1n, 16n],
      [2n, -1n, 8n],
      [2n, 1n, 8n],
      [3n, 3n, 6n],
      [4n, 1n, 4n],
    ]);
  });

  // binary_qf.py:1883-1898
  test('D = -207 with and without primitive_only', () => {
    expect(tuples(BinaryQF_reduced_representatives(-23n * 9n))).toEqual([
      [1n, 1n, 52n],
      [2n, -1n, 26n],
      [2n, 1n, 26n],
      [3n, 3n, 18n],
      [4n, -1n, 13n],
      [4n, 1n, 13n],
      [6n, -3n, 9n],
      [6n, 3n, 9n],
      [8n, 7n, 8n],
    ]);
    expect(tuples(BinaryQF_reduced_representatives(-23n * 9n, { primitive_only: true }))).toEqual([
      [1n, 1n, 52n],
      [2n, -1n, 26n],
      [2n, 1n, 26n],
      [4n, -1n, 13n],
      [4n, 1n, 13n],
      [8n, 7n, 8n],
    ]);
  });

  // binary_qf.py:1902-1915 (indefinite, non-square discriminant)
  test('D = 73, 76 and 136', () => {
    expect(tuples(BinaryQF_reduced_representatives(73n))).toEqual([[4n, 3n, -4n]]);
    expect(tuples(BinaryQF_reduced_representatives(76n, { primitive_only: true }))).toEqual([
      [-3n, 4n, 5n],
      [3n, 4n, -5n],
    ]);
    expect(tuples(BinaryQF_reduced_representatives(136n))).toEqual([
      [-5n, 4n, 6n],
      [-2n, 8n, 9n],
      [2n, 8n, -9n],
      [5n, 4n, -6n],
    ]);
    expect(tuples(BinaryQF_reduced_representatives(136n, { proper: false }))).toEqual([
      [-2n, 8n, 9n],
      [2n, 8n, -9n],
      [5n, 4n, -6n],
    ]);
  });

  // binary_qf.py:1917-1928
  test('D = 148 with all four flag combinations', () => {
    expect(
      tuples(BinaryQF_reduced_representatives(148n, { proper: false, primitive_only: false }))
    ).toEqual([
      [1n, 12n, -1n],
      [4n, 6n, -7n],
      [6n, 2n, -6n],
    ]);
    expect(
      tuples(BinaryQF_reduced_representatives(148n, { proper: false, primitive_only: true }))
    ).toEqual([
      [1n, 12n, -1n],
      [4n, 6n, -7n],
    ]);
    expect(
      tuples(BinaryQF_reduced_representatives(148n, { proper: true, primitive_only: true }))
    ).toEqual([
      [-7n, 6n, 4n],
      [1n, 12n, -1n],
      [4n, 6n, -7n],
    ]);
    expect(
      tuples(BinaryQF_reduced_representatives(148n, { proper: true, primitive_only: false }))
    ).toEqual([
      [-7n, 6n, 4n],
      [1n, 12n, -1n],
      [4n, 6n, -7n],
      [6n, 2n, -6n],
    ]);
  });

  // binary_qf.py:1930-1954 (issue 29028): square discriminant D = 100
  test('D = 100 (square discriminant)', () => {
    expect(
      tuples(BinaryQF_reduced_representatives(100n, { proper: false, primitive_only: false }))
    ).toEqual([
      [-4n, 10n, 0n],
      [-3n, 10n, 0n],
      [-2n, 10n, 0n],
      [-1n, 10n, 0n],
      [0n, 10n, 0n],
      [1n, 10n, 0n],
      [2n, 10n, 0n],
      [5n, 10n, 0n],
    ]);
    expect(
      tuples(BinaryQF_reduced_representatives(100n, { proper: false, primitive_only: true }))
    ).toEqual([
      [-3n, 10n, 0n],
      [-1n, 10n, 0n],
      [1n, 10n, 0n],
    ]);
    expect(
      tuples(BinaryQF_reduced_representatives(100n, { proper: true, primitive_only: true }))
    ).toEqual([
      [-3n, 10n, 0n],
      [-1n, 10n, 0n],
      [1n, 10n, 0n],
      [3n, 10n, 0n],
    ]);
    expect(
      BinaryQF_reduced_representatives(100n, { proper: true, primitive_only: false }).length
    ).toBe(10);
  });

  // M85: the lower bound must floor, not truncate
  test('small square discriminants', () => {
    expect(tuples(BinaryQF_reduced_representatives(1n))).toEqual([[0n, 1n, 0n]]);
    expect(BinaryQF_reduced_representatives(9n).length).toBe(3);
    expect(BinaryQF_reduced_representatives(25n).length).toBe(5);
  });

  test('every representative is reduced and has the right discriminant', () => {
    for (const D of [-4n, -23n, -207n, 73n, 76n, 100n, 136n, 148n, 229n, 9n, 25n]) {
      for (const primitive_only of [false, true]) {
        for (const proper of [false, true]) {
          const reps = BinaryQF_reduced_representatives(D, { primitive_only, proper });
          for (const f of reps) {
            expect(f.discriminant()).toBe(D);
            expect(f.is_reduced()).toBe(true);
            if (primitive_only) expect(f.is_primitive()).toBe(true);
          }
          // representatives are pairwise inequivalent
          for (let i = 0; i < reps.length; i++) {
            for (let j = i + 1; j < reps.length; j++) {
              expect(reps[i]!.is_equivalent(reps[j]!, { proper })).toBe(false);
            }
          }
        }
      }
    }
  });

  test('rejects non-discriminants', () => {
    expect(() => BinaryQF_reduced_representatives(2n)).toThrow(/is not a discriminant/);
    expect(() => BinaryQF_reduced_representatives(-1n)).toThrow(/is not a discriminant/);
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
