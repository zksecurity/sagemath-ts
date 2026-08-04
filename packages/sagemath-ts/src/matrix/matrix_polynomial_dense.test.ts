/**
 * Tests for `sage/matrix/matrix_polynomial_dense.pyx`.
 *
 * Every expected value in the "doctest" blocks below is copied verbatim from
 * the docstrings of `reference/sage/src/sage/matrix/matrix_polynomial_dense.pyx`
 * (line numbers cited per block).  The remaining blocks are independent
 * oracles: exhaustive brute force over small matrices for the form predicates,
 * and algebraic identities (`U*A == P`, `U` unimodular, idempotence of the
 * canonical forms, Hermite = shifted Popov) for the algorithms.
 */

import { describe, expect, it } from 'bun:test';
import { GF } from '../rings/finite_rings/finite_field_constructor.js';
import { GFExtended } from '../rings/finite_rings/finite_field_extension.js';
import { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';
import { QQ } from '../rings/rational_field.js';
import { Matrix } from './matrix_generic.js';
import { determinant } from './matrix_operations.js';
import * as MP from './matrix_polynomial_dense.js';

// ============================================================================
// Helpers
// ============================================================================

/** Parse a Sage-printed univariate polynomial such as `"x^3 + 5*x^2 + 5*x + 1"`. */
function poly(pR: any, s: string): any {
  const src = s.replace(/\s+/g, '');
  const varName: string = pR.variable_name;
  const coeffs: bigint[] = [];
  const add = (c: bigint, d: number) => {
    while (coeffs.length <= d) coeffs.push(0n);
    coeffs[d] = coeffs[d]! + c;
  };
  let i = 0;
  let sign = 1n;
  if (src[0] === '+' || src[0] === '-') {
    sign = src[0] === '-' ? -1n : 1n;
    i = 1;
  }
  while (i < src.length) {
    let j = i;
    while (j < src.length && src[j] !== '+' && src[j] !== '-') j++;
    const term = src.slice(i, j);
    let c = 1n;
    let d = 0;
    const star = term.indexOf('*');
    let rest = term;
    if (star >= 0 && term.slice(star + 1).startsWith(varName)) {
      c = BigInt(term.slice(0, star));
      rest = term.slice(star + 1);
    } else if (term.startsWith(varName)) {
      rest = term;
    } else {
      c = BigInt(term);
      rest = '';
    }
    if (rest.startsWith(varName)) {
      const pw = rest.slice(varName.length);
      d = pw === '' ? 1 : Number(pw.replace('^', ''));
    }
    add(sign * c, d);
    if (j < src.length) {
      sign = src[j] === '-' ? -1n : 1n;
      i = j + 1;
    } else {
      i = j;
    }
  }
  if (coeffs.length === 0) coeffs.push(0n);
  return pR.__call__(coeffs.map((c) => pR.base_ring.__call__(c)));
}

/** Build a polynomial matrix from Sage-printed entries. */
function pmat(pR: any, rows: string[][]): any {
  return new Matrix(
    pR as any,
    rows.length,
    rows[0]!.length,
    rows.map((r) => r.map((s) => poly(pR, s)))
  );
}

/** Build the `m x n` zero polynomial matrix. */
function pzero(pR: any, m: number, n: number): any {
  return new Matrix(pR as any, m, n);
}

/** Render a matrix the way SageMath prints its rows, for readable assertions. */
function render(M: any): string {
  const rows: string[] = [];
  for (let i = 0; i < M.nrows; i++) {
    const r: string[] = [];
    for (let j = 0; j < M.ncols; j++) r.push(String(M.get(i, j)));
    rows.push('[' + r.join(' ') + ']');
  }
  return rows.join(' / ');
}

function matEq(A: any, B: any): boolean {
  if (A.nrows !== B.nrows || A.ncols !== B.ncols) return false;
  for (let i = 0; i < A.nrows; i++) {
    for (let j = 0; j < A.ncols; j++) {
      if (!A.get(i, j).eq(B.get(i, j))) return false;
    }
  }
  return true;
}

/** Deterministic small PRNG so the property tests are reproducible. */
function makeRng(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s >>> 1;
  };
}

const F7 = GF(7n) as any;
const pR = new PolynomialRing(F7, 'x');

// The running example of the module docstrings (pyx:240-246 and elsewhere).
const Mbig = pmat(pR, [
  ['x^3+5*x^2+5*x+1', '5', '6*x+4', '0'],
  ['6*x^2+3*x+1', '1', '2', '0'],
  ['2*x^3+4*x^2+6*x+4', '5*x + 1', '2*x^2+5*x+5', 'x^2+5*x+6'],
]);

// M = matrix(pR, [[3*x+1, 0, 1], [x^3+3, 0, 0]])  (pyx:123)
const Msmall = pmat(pR, [
  ['3*x+1', '0', '1'],
  ['x^3+3', '0', '0'],
]);

// ============================================================================
// _check_shift_dimension, degree, degree_matrix  (pyx:104-235)
// ============================================================================

describe('_check_shift_dimension (pyx:104-135)', () => {
  it('accepts a shift of the column dimension when row-wise', () => {
    expect(() => MP._check_shift_dimension(Msmall, [1, 3, 2])).not.toThrow();
  });

  it('rejects a shift of the wrong length when column-wise', () => {
    expect(() => MP._check_shift_dimension(Msmall, [1, 3, 2], false)).toThrowError(
      'shifts length should be the row dimension'
    );
  });
});

describe('degree (pyx:137-166)', () => {
  it('is the maximum entry degree', () => {
    expect(MP.degree(Msmall)).toBe(3);
  });

  it('is -1 for zero and empty matrices', () => {
    expect(MP.degree(pzero(pR, 2, 3))).toBe(-1);
    expect(MP.degree(pzero(pR, 3, 0))).toBe(-1);
  });
});

describe('degree_matrix (pyx:168-235)', () => {
  it('matches the unshifted doctest', () => {
    expect(MP.degree_matrix(Msmall)).toEqual([
      [1, -1, 0],
      [3, -1, -1],
    ]);
  });

  it('matches the shifted doctests', () => {
    expect(MP.degree_matrix(Msmall, { shifts: [0, 1, 2] })).toEqual([
      [1, -1, 2],
      [3, -1, -1],
    ]);
    expect(MP.degree_matrix(Msmall, { shifts: [-2, 1, 2] })).toEqual([
      [-1, -3, 2],
      [1, -3, -3],
    ]);
    expect(MP.degree_matrix(Msmall, { shifts: [-1, 2], row_wise: false })).toEqual([
      [0, -2, -1],
      [5, -2, -2],
    ]);
  });
});

// ============================================================================
// constant_matrix, is_constant, coefficient_matrix  (pyx:237-382)
// ============================================================================

describe('constant_matrix / is_constant (pyx:237-289)', () => {
  it('matches the constant_matrix doctest', () => {
    expect(render(MP.constant_matrix(Mbig))).toBe('[1 5 4 0] / [1 1 2 0] / [4 1 5 6]');
  });

  it('matches the is_constant doctests', () => {
    expect(MP.is_constant(Mbig)).toBe(false);
    expect(
      MP.is_constant(
        pmat(pR, [
          ['1', '5', '2'],
          ['3', '1', '5'],
        ])
      )
    ).toBe(true);
    expect(MP.is_constant(pzero(pR, 3, 5))).toBe(true);
  });
});

describe('coefficient_matrix (pyx:291-382)', () => {
  it('matches the integer-degree doctests', () => {
    expect(render(MP.coefficient_matrix(Mbig, 2))).toBe('[5 0 0 0] / [6 0 0 0] / [4 0 2 1]');
    expect(render(MP.coefficient_matrix(Mbig, 0))).toBe(render(MP.constant_matrix(Mbig)));
  });

  it('matches the row-wise and column-wise doctests', () => {
    expect(render(MP.coefficient_matrix(Mbig, [3, 2, 1]))).toBe(
      '[1 0 0 0] / [6 0 0 0] / [6 5 5 5]'
    );
    expect(render(MP.coefficient_matrix(Mbig, [2, 0, 1, 3], { row_wise: false }))).toBe(
      '[5 5 6 0] / [6 1 0 0] / [4 1 5 0]'
    );
  });

  it('gives zero coefficients for negative degrees', () => {
    expect(render(MP.coefficient_matrix(Mbig, [-1, 0, 1, 3], { row_wise: false }))).toBe(
      '[0 5 6 0] / [0 1 0 0] / [0 1 5 0]'
    );
  });

  it('checks the length of the degree list', () => {
    expect(() => MP.coefficient_matrix(Mbig, [2, 1, 1, 2])).toThrowError(
      'length of input degree list should be the row dimension of the input matrix'
    );
    expect(() => MP.coefficient_matrix(Mbig, [3, 2, 1], { row_wise: false })).toThrowError(
      'length of input degree list should be the column dimension of the input matrix'
    );
  });
});

// ============================================================================
// truncate, shift, reverse  (pyx:384-712)
// ============================================================================

describe('truncate (pyx:384-475)', () => {
  it('matches the doctests', () => {
    expect(render(MP.truncate(Mbig, 2))).toBe(
      '[5*x + 1 5 6*x + 4 0] / [3*x + 1 1 2 0] / [6*x + 4 5*x + 1 5*x + 5 5*x + 6]'
    );
    expect(render(MP.truncate(Mbig, 1))).toBe(render(MP.constant_matrix(Mbig)));
    expect(render(MP.truncate(Mbig, [3, 2, 1]))).toBe(
      '[5*x^2 + 5*x + 1 5 6*x + 4 0] / [3*x + 1 1 2 0] / [4 1 5 6]'
    );
    expect(render(MP.truncate(Mbig, [2, 1, 1, 2], { row_wise: false }))).toBe(
      '[5*x + 1 5 4 0] / [3*x + 1 1 2 0] / [6*x + 4 1 5 5*x + 6]'
    );
  });

  it('checks the length of the precision list', () => {
    expect(() => MP.truncate(Mbig, [2, 1, 1, 2])).toThrowError(
      'length of input precision list should be the row dimension of the input matrix'
    );
    expect(() => MP.truncate(Mbig, [3, 2, 1], { row_wise: false })).toThrowError(
      'length of input precision list should be the column dimension of the input matrix'
    );
  });
});

describe('shift (pyx:477-567)', () => {
  it('matches the doctests', () => {
    expect(render(MP.shift(Mbig, -2))).toBe('[x + 5 0 0 0] / [6 0 0 0] / [2*x + 4 0 2 1]');
    expect(render(MP.shift(Mbig, [-1, 2, -2]))).toBe(
      '[x^2 + 5*x + 5 0 6 0] / [6*x^4 + 3*x^3 + x^2 x^2 2*x^2 0] / [2*x + 4 0 2 1]'
    );
    expect(render(MP.shift(Mbig, [-1, 1, 0, 0], { row_wise: false }))).toBe(
      '[x^2 + 5*x + 5 5*x 6*x + 4 0] / [6*x + 3 x 2 0] / ' +
        '[2*x^2 + 4*x + 6 5*x^2 + x 2*x^2 + 5*x + 5 x^2 + 5*x + 6]'
    );
  });

  it('satisfies M.shift(-M.row_degrees()) == M.leading_matrix()', () => {
    expect(
      render(
        MP.shift(
          Mbig,
          MP.row_degrees(Mbig).map((d) => -d)
        )
      )
    ).toBe(render(MP.leading_matrix(Mbig)));
  });

  it('checks the length of the shift list', () => {
    expect(() => MP.shift(Mbig, [1, 3, 1, 4])).toThrowError(
      'length of input shift list should be the row dimension of the input matrix'
    );
    expect(() => MP.shift(Mbig, [5, 2, -1], { row_wise: false })).toThrowError(
      'length of input shift list should be the column dimension of the input matrix'
    );
  });
});

describe('reverse (pyx:569-712)', () => {
  it('matches the whole-matrix-degree doctests', () => {
    expect(render(MP.reverse(Mbig))).toBe(
      '[x^3 + 5*x^2 + 5*x + 1 5*x^3 4*x^3 + 6*x^2 0] / ' +
        '[x^3 + 3*x^2 + 6*x x^3 2*x^3 0] / ' +
        '[4*x^3 + 6*x^2 + 4*x + 2 x^3 + 5*x^2 5*x^3 + 5*x^2 + 2*x 6*x^3 + 5*x^2 + x]'
    );
    expect(render(MP.reverse(Mbig, { degree: 1 }))).toBe(
      '[x + 5 5*x 4*x + 6 0] / [x + 3 x 2*x 0] / [4*x + 6 x + 5 5*x + 5 6*x + 5]'
    );
    expect(render(MP.reverse(Mbig, { degree: 0 }))).toBe(render(MP.constant_matrix(Mbig)));
  });

  it('matches the entry-wise doctest', () => {
    expect(render(MP.reverse(Mbig, { entry_wise: true }))).toBe(
      '[x^3 + 5*x^2 + 5*x + 1 5 4*x + 6 0] / ' +
        '[x^2 + 3*x + 6 1 2 0] / ' +
        '[4*x^3 + 6*x^2 + 4*x + 2 x + 5 5*x^2 + 5*x + 2 6*x^2 + 5*x + 1]'
    );
  });

  it('matches the row-wise and column-wise doctests', () => {
    expect(render(MP.reverse(Mbig, { degree: [2, 3, 1] }))).toBe(
      '[x^2 + 5*x + 5 5*x^2 4*x^2 + 6*x 0] / ' +
        '[x^3 + 3*x^2 + 6*x x^3 2*x^3 0] / ' +
        '[4*x + 6 x + 5 5*x + 5 6*x + 5]'
    );
    expect(render(MP.reverse(Mbig, { degree: MP.column_degrees(Mbig), row_wise: false }))).toBe(
      '[x^3 + 5*x^2 + 5*x + 1 5*x 4*x^2 + 6*x 0] / ' +
        '[x^3 + 3*x^2 + 6*x x 2*x^2 0] / ' +
        '[4*x^3 + 6*x^2 + 4*x + 2 x + 5 5*x^2 + 5*x + 2 6*x^2 + 5*x + 1]'
    );
  });

  it('checks the length and the sign of the degree list', () => {
    expect(() => MP.reverse(Mbig, { degree: [1, 3, 1, 4] })).toThrowError(
      'length of input degree list should be the row dimension of the input matrix'
    );
    expect(() => MP.reverse(Mbig, { degree: [5, 2, 1], row_wise: false })).toThrowError(
      'length of input degree list should be the column dimension of the input matrix'
    );
    expect(() => MP.reverse(Mbig, { degree: [2, 3, -1] })).toThrowError(
      'degree argument must be a nonnegative integer, got -1'
    );
  });
});

// ============================================================================
// row_degrees / column_degrees  (pyx:1120-1266)
// ============================================================================

describe('row_degrees (pyx:1120-1197)', () => {
  it('matches the doctests', () => {
    expect(MP.row_degrees(Msmall)).toEqual([1, 3]);
    expect(MP.row_degrees(Msmall, [0, 1, 2])).toEqual([2, 3]);
  });

  it('identifies zero rows as min(shifts) - 1', () => {
    const M = pmat(pR, [
      ['3*x+1', '0', '1'],
      ['x^3+3', '0', '0'],
      ['0', '0', '0'],
    ]);
    expect(MP.row_degrees(M)).toEqual([1, 3, -1]);
    expect(MP.row_degrees(M, [-2, 1, 2])).toEqual([2, 1, -3]);
  });

  it('handles empty matrices', () => {
    expect(MP.row_degrees(pzero(pR, 0, 3))).toEqual([]);
    expect(MP.row_degrees(pzero(pR, 0, 3), [1, 2, 3])).toEqual([]);
    expect(MP.row_degrees(pzero(pR, 3, 0))).toEqual([-1, -1, -1]);
    expect(MP.row_degrees(pzero(pR, 3, 0), [])).toEqual([-1, -1, -1]);
  });
});

describe('column_degrees (pyx:1198-1266)', () => {
  it('matches the doctests', () => {
    expect(MP.column_degrees(Msmall)).toEqual([3, -1, 0]);
    expect(MP.column_degrees(Msmall, [0, 2])).toEqual([5, -1, 0]);
    expect(MP.column_degrees(Msmall, [-2, 1])).toEqual([4, -3, -2]);
  });

  it('handles empty matrices', () => {
    expect(MP.column_degrees(pzero(pR, 3, 0))).toEqual([]);
    expect(MP.column_degrees(pzero(pR, 3, 0), [1, 2, 3])).toEqual([]);
    expect(MP.column_degrees(pzero(pR, 0, 3))).toEqual([-1, -1, -1]);
    expect(MP.column_degrees(pzero(pR, 0, 3), [])).toEqual([-1, -1, -1]);
  });
});

// ============================================================================
// leading_matrix, leading_positions, is_reduced, _is_empty_popov
// ============================================================================

describe('leading_matrix (pyx:1268-1362)', () => {
  it('matches the doctests', () => {
    expect(render(MP.leading_matrix(Msmall))).toBe('[3 0 0] / [1 0 0]');
    expect(render(MP.leading_matrix(Msmall, { shifts: [0, 1, 2] }))).toBe('[0 0 1] / [1 0 0]');
    expect(render(MP.leading_matrix(Msmall, { row_wise: false }))).toBe('[0 0 1] / [1 0 0]');
    expect(render(MP.leading_matrix(Msmall, { shifts: [-2, 1], row_wise: false }))).toBe(
      '[0 0 1] / [1 0 0]'
    );
    expect(render(MP.leading_matrix(Msmall, { shifts: [2, 0], row_wise: false }))).toBe(
      '[3 0 1] / [1 0 0]'
    );
  });

  it('matches the empty-matrix doctests', () => {
    expect(render(MP.leading_matrix(pzero(pR, 0, 3), { shifts: [1, 2, 3] }))).toBe('');
    expect(render(MP.leading_matrix(pzero(pR, 0, 3), { shifts: [], row_wise: false }))).toBe('');
  });
});

describe('_is_empty_popov (pyx:1364-1415)', () => {
  it('matches the doctests', () => {
    expect(MP._is_empty_popov(pzero(pR, 0, 0))).toBe(true);
    expect(MP._is_empty_popov(pzero(pR, 0, 0), true, false)).toBe(true);
    expect(MP._is_empty_popov(pzero(pR, 0, 3), true, false)).toBe(true);
    expect(MP._is_empty_popov(pzero(pR, 0, 3), false)).toBe(true);
    expect(MP._is_empty_popov(pzero(pR, 0, 3), false, false)).toBe(false);
  });
});

describe('is_reduced (pyx:1417-1495)', () => {
  it('matches the doctests', () => {
    expect(MP.is_reduced(Msmall)).toBe(false);
    expect(MP.is_reduced(Msmall, { shifts: [0, 1, 2] })).toBe(true);
    expect(MP.is_reduced(Msmall, { shifts: [2, 0], row_wise: false })).toBe(true);
    expect(
      MP.is_reduced(Msmall, { shifts: [2, 0], row_wise: false, include_zero_vectors: false })
    ).toBe(false);
    const M = pmat(pR, [
      ['3*x+1', '0', '1'],
      ['x^3+3', '0', '0'],
      ['0', '1', '0'],
    ]);
    expect(MP.is_reduced(M, { shifts: [2, 0, 0], row_wise: false })).toBe(true);
  });
});

describe('leading_positions (pyx:1497-1621)', () => {
  it('matches the doctests', () => {
    expect(MP.leading_positions(Msmall)).toEqual([0, 0]);
    expect(MP.leading_positions(Msmall, { return_degree: true })).toEqual([
      [0, 0],
      [1, 3],
    ]);
    expect(MP.leading_positions(Msmall, { shifts: [0, 5, 2], return_degree: true })).toEqual([
      [2, 0],
      [0, 3],
    ]);
    expect(MP.leading_positions(Msmall, { row_wise: false, return_degree: true })).toEqual([
      [1, -1, 0],
      [3, -1, 0],
    ]);
    expect(
      MP.leading_positions(Msmall, { shifts: [1, 2], row_wise: false, return_degree: true })
    ).toEqual([
      [1, -1, 0],
      [3, -1, 0],
    ]);
  });

  it('chooses the rightmost (resp. bottommost) entry on ties', () => {
    expect(MP.leading_positions(Msmall, { shifts: [0, 5, 1], return_degree: true })).toEqual([
      [2, 0],
      [0, 3],
    ]);
    expect(
      MP.leading_positions(Msmall, { shifts: [2, 0], row_wise: false, return_degree: true })
    ).toEqual([
      [1, -1, 0],
      [3, -1, 0],
    ]);
  });

  it('matches the empty-matrix conventions', () => {
    expect(MP.leading_positions(pzero(pR, 0, 3), { return_degree: true })).toEqual([[], []]);
    expect(
      MP.leading_positions(pzero(pR, 0, 3), { shifts: [], row_wise: false, return_degree: true })
    ).toEqual([
      [-1, -1, -1],
      [-1, -1, -1],
    ]);
  });
});

// ============================================================================
// is_weak_popov / is_popov / is_hermite
// ============================================================================

describe('is_weak_popov (pyx:1623-1752)', () => {
  const M = pmat(pR, [
    ['x^3+3*x^2+6*x+6', '3*x^2+3*x+6', '4*x^2+x+3'],
    ['5', '1', '0'],
    ['2*x^2+2', '2*x+5', 'x^2+4*x+6'],
  ]);

  it('matches the base doctests', () => {
    expect(MP.is_weak_popov(M)).toBe(true);
    expect(MP.is_weak_popov(M, { ordered: true })).toBe(true);
  });

  it('detects a non-ordered form after swapping rows 1 and 2', () => {
    const N = pmat(pR, [
      ['x^3+3*x^2+6*x+6', '3*x^2+3*x+6', '4*x^2+x+3'],
      ['2*x^2+2', '2*x+5', 'x^2+4*x+6'],
      ['5', '1', '0'],
    ]);
    expect(MP.is_weak_popov(N)).toBe(true);
    expect(MP.is_weak_popov(N, { ordered: true })).toBe(false);
  });

  it('supports shifts and column-wise orientation', () => {
    expect(MP.is_weak_popov(M, { shifts: [2, 3, 1] })).toBe(false);
    expect(MP.is_weak_popov(M, { shifts: [0, 2, 0], row_wise: false, ordered: true })).toBe(true);
  });

  it('supports rectangular matrices', () => {
    expect(MP.is_weak_popov(Mbig, { shifts: [0, 2, 1, 3] })).toBe(true);
    expect(MP.is_weak_popov(Mbig, { shifts: [0, 2, 1, 3], ordered: true })).toBe(true);
  });

  it('can forbid zero rows/columns', () => {
    const Z = pmat(pR, [
      ['6*x+4', '0', '5*x+1', '0'],
      ['2', '5*x + 1', '6*x^2+3*x+1', '0'],
      ['2*x^2+5*x+5', '1', '2*x^3+4*x^2+6*x+4', '0'],
    ]);
    expect(MP.is_weak_popov(Z, { shifts: [2, 1, 0], row_wise: false, ordered: true })).toBe(true);
    expect(
      MP.is_weak_popov(Z, { shifts: [2, 1, 0], row_wise: false, include_zero_vectors: false })
    ).toBe(false);
  });
});

describe('is_popov (pyx:1754-1891)', () => {
  const M = pmat(pR, [
    ['x^4+6*x^3+4*x+4', '3*x+6', '3'],
    ['x^2+6*x+6', 'x^2+5*x+5', '2'],
    ['3*x', '6*x+5', 'x+5'],
  ]);
  const N = pmat(pR, [
    ['x^4+3*x^3+x^2+2*x+6', 'x^3+5*x^2+5*x+1'],
    ['6*x+1', 'x^2+4*x+1'],
    ['6', '6'],
  ]);

  it('matches the base doctests', () => {
    expect(MP.is_popov(M)).toBe(true);
    expect(MP.is_popov(M, { shifts: [0, 1, 2] })).toBe(true);
    expect(
      MP.is_popov(
        pmat(pR, [
          ['x^4+6*x^3+4*x+4', '3*x+6'],
          ['x^2+6*x+6', 'x^2+5*x+5'],
          ['3*x', '6*x+5'],
        ])
      )
    ).toBe(false);
    expect(
      MP.is_popov(
        pmat(pR, [
          ['x^4+6*x^3+4*x+4', '3*x+6', '3'],
          ['x^2+6*x+6', 'x^2+5*x+5', '2'],
        ]),
        { shifts: [0, 1, 2] }
      )
    ).toBe(true);
  });

  it('matches the column-wise doctests', () => {
    expect(MP.is_popov(N, { row_wise: false })).toBe(false);
    expect(MP.is_popov(N, { shifts: [0, 2, 3], row_wise: false })).toBe(true);
  });

  it('can forbid zero rows', () => {
    const Z = pmat(pR, [
      ['x^4+3*x^3+x^2+2*x+6', '6*x+1'],
      ['5*x^2+5*x+1', 'x^2+4*x+1'],
      ['0', '0'],
    ]);
    expect(MP.is_popov(Z)).toBe(true);
    expect(MP.is_popov(Z, { include_zero_vectors: false })).toBe(false);
  });

  it('supports the up_to_permutation option', () => {
    // N.swap_columns(0, 1)
    const Ns = pmat(pR, [
      ['x^3+5*x^2+5*x+1', 'x^4+3*x^3+x^2+2*x+6'],
      ['x^2+4*x+1', '6*x+1'],
      ['6', '6'],
    ]);
    expect(MP.is_popov(Ns, { shifts: [0, 2, 3], row_wise: false })).toBe(false);
    expect(MP.is_popov(Ns, { shifts: [0, 2, 3], row_wise: false, up_to_permutation: true })).toBe(
      true
    );

    // Z.swap_rows(0, 2)
    const Zs = pmat(pR, [
      ['0', '0'],
      ['5*x^2+5*x+1', 'x^2+4*x+1'],
      ['x^4+3*x^3+x^2+2*x+6', '6*x+1'],
    ]);
    expect(MP.is_popov(Zs)).toBe(false);
    expect(MP.is_popov(Zs, { up_to_permutation: true })).toBe(true);
  });
});

describe('is_hermite (pyx:1893-2008)', () => {
  const M = pmat(pR, [
    ['x^4+6*x^3+4*x+4', '3*x+6', '3'],
    ['0', 'x^2+5*x+5', '2'],
    ['0', '0', 'x+5'],
  ]);
  const N = pmat(pR, [
    ['x+5', '0', '0'],
    ['2', 'x^4+6*x^3+4*x+4', '0'],
    ['3', '3*x^3+6', 'x^2+5*x+5'],
  ]);

  it('matches the upper-echelon doctests', () => {
    expect(MP.is_hermite(M)).toBe(true);
    expect(MP.is_hermite(M, { row_wise: false })).toBe(true);
    expect(MP.is_hermite(M, { row_wise: false, lower_echelon: true })).toBe(false);
  });

  it('matches the lower-echelon doctests', () => {
    expect(MP.is_hermite(N)).toBe(false);
    expect(MP.is_hermite(N, { lower_echelon: true })).toBe(true);
    expect(MP.is_hermite(N, { row_wise: false })).toBe(false);
    expect(MP.is_hermite(N, { row_wise: false, lower_echelon: true })).toBe(false);
  });

  it('supports rectangular matrices with zero rows', () => {
    // N[:,1:]
    expect(
      MP.is_hermite(
        pmat(pR, [
          ['0', '0'],
          ['x^4+6*x^3+4*x+4', '0'],
          ['3*x^3+6', 'x^2+5*x+5'],
        ]),
        { lower_echelon: true }
      )
    ).toBe(false);
    // N[[1,2,0],1:]
    expect(
      MP.is_hermite(
        pmat(pR, [
          ['x^4+6*x^3+4*x+4', '0'],
          ['3*x^3+6', 'x^2+5*x+5'],
          ['0', '0'],
        ]),
        { lower_echelon: true }
      )
    ).toBe(true);
    // N[:2,:]
    const N2 = pmat(pR, [
      ['x+5', '0', '0'],
      ['2', 'x^4+6*x^3+4*x+4', '0'],
    ]);
    expect(MP.is_hermite(N2, { row_wise: false, lower_echelon: true })).toBe(true);
    expect(
      MP.is_hermite(N2, {
        row_wise: false,
        lower_echelon: true,
        include_zero_vectors: false,
      })
    ).toBe(false);
  });

  it('matches the empty-matrix doctests', () => {
    expect(MP.is_hermite(pzero(pR, 0, 3))).toBe(true);
    expect(MP.is_hermite(pzero(pR, 0, 3), { row_wise: false })).toBe(true);
    expect(MP.is_hermite(pzero(pR, 0, 3), { row_wise: false, include_zero_vectors: false })).toBe(
      false
    );
  });
});

// ============================================================================
// weak_popov_form  (pyx:2010-2343)
// ============================================================================

describe('weak_popov_form (pyx:2010-2227)', () => {
  const M = pmat(pR, [
    ['6*x+4', '5*x^3+5*x', '6*x^2+2*x+2'],
    ['4*x^2+5*x+2', 'x^4+5*x^2+2*x+4', '4*x^3+6*x^2+6*x+5'],
  ]);

  it('matches the transformation doctest', () => {
    const [P, U] = MP.weak_popov_form(M, { transformation: true });
    expect(render(P)).toBe('[4 x^2 6*x^2 + x + 2] / [2 4*x^2 + 2*x + 4 5]');
    expect(render(U)).toBe('[2*x^2 + 1 4*x] / [4*x 1]');
    expect(MP.is_weak_popov(P)).toBe(true);
    expect(matEq(U.mul(M), P)).toBe(true);
  });

  it('matches the ordered doctest', () => {
    const P = MP.weak_popov_form(M);
    expect(MP.leading_positions(P)).toEqual([2, 1]);
    const PP = MP.weak_popov_form(M, { ordered: true });
    expect(render(PP)).toBe('[2 4*x^2 + 2*x + 4 5] / [4 x^2 6*x^2 + x + 2]');
    expect(MP.leading_positions(PP)).toEqual([1, 2]);
  });

  it('matches the shifts doctest', () => {
    const P = MP.weak_popov_form(M, { shifts: [0, 2, 4] });
    expect(render(P)).toBe(
      '[6*x^2 + 6*x + 4 5*x^4 + 4*x^3 + 5*x^2 + 5*x 2*x + 2] / [2 4*x^2 + 2*x + 4 5]'
    );
    expect(matEq(P, MP.weak_popov_form(M, { shifts: [-10, -8, -6] }))).toBe(true);
  });

  it('is the transpose of the column-wise form of the transpose', () => {
    expect(
      matEq(
        MP.weak_popov_form(M),
        MP.weak_popov_form(M.transpose(), { row_wise: false }).transpose()
      )
    ).toBe(true);
  });

  it('matches the discard-zero-vectors doctests', () => {
    expect(render(MP.weak_popov_form(M, { row_wise: false }))).toBe('[x + 4 6 0] / [5 1 0]');
    const [P, U] = MP.weak_popov_form(M, {
      transformation: true,
      row_wise: false,
      include_zero_vectors: false,
    });
    expect(render(P)).toBe('[x + 4 6] / [5 1]');
    expect(render(U)).toBe(
      '[5*x + 2 5*x^2 + 4*x + 4 3*x^3 + 3*x^2 + 2*x + 4] / [1 1 2*x + 1] / [5*x + 5 2 6]'
    );
    const MU = M.mul(U);
    expect(MU.get(0, 0).eq(P.get(0, 0)) && MU.get(1, 1).eq(P.get(1, 1))).toBe(true);
    expect(MU.get(0, 2).isZero() && MU.get(1, 2).isZero()).toBe(true);
  });

  it('supports empty matrices', () => {
    expect(render(MP.weak_popov_form(pzero(pR, 0, 3)))).toBe('');
    const [E, EU] = MP.weak_popov_form(pzero(pR, 0, 3), { transformation: true });
    expect([render(E), render(EU)]).toEqual(['', '']);
    const [E2, EU2] = MP.weak_popov_form(pzero(pR, 0, 3), {
      transformation: true,
      row_wise: false,
    });
    expect([render(E2), render(EU2)]).toEqual(['', '[1 0 0] / [0 1 0] / [0 0 1]']);
  });

  it('matches the issue #41278 regression test (pyx:2150-2157)', () => {
    const F3 = GF(3n) as any;
    const pR3 = new PolynomialRing(F3, 'x');
    const A = pmat(pR3, [
      ['x^3 + x', '0', '0'],
      ['2*x^2', 'x', '0'],
      ['x', '0', 'x'],
      ['x^2 + 1', 'x^2 + 1', '0'],
      ['2*x + 2', '2*x + 2', 'x'],
      ['x^2 + x + 1', 'x^2 + 2*x + 1', '2*x^3 + 2*x^2'],
      ['0', '0', 'x^2 + 1'],
      ['x^2 + x', 'x^2 + 2*x', '2*x^3 + 2*x^2 + 2*x + 2'],
      ['2*x^4 + x^3 + 2*x^2 + 2', '2*x^4 + x^2 + 2', 'x^5 + 2*x^4 + x^3 + x^2 + 2*x + 1'],
    ]);
    expect(render(MP.weak_popov_form(A, { ordered: true, include_zero_vectors: false }))).toBe(
      '[x + 2 2 2] / [0 2*x 1] / [x 0 x]'
    );
  });
});

describe('_weak_popov_form (pyx:2229-2343)', () => {
  const pQ = new PolynomialRing(QQ as any, 'x');
  const A = pmat(pQ, [
    ['x', 'x^2', 'x^3'],
    ['x^2', 'x', '0'],
    ['x^3', 'x^3', 'x^3'],
  ]);

  it('matches the QQ doctest for weak_popov_form', () => {
    // [        x       x^2       x^3]
    // [      x^2         x         0]
    // [  x^3 - x x^3 - x^2         0]
    const expected = pmat(pQ, [
      ['x', 'x^2', 'x^3'],
      ['x^2', 'x', '0'],
      ['x^3-x', 'x^3-x^2', '0'],
    ]);
    expect(matEq(MP.weak_popov_form(A), expected)).toBe(true);
  });

  it('matches the QQ doctest with shifts=[16,8,0]', () => {
    const M = A.copy();
    const U = MP._weak_popov_form(M, { transformation: true, shifts: [16, 8, 0] })!;
    const expected = pmat(pQ, [
      ['x', 'x^2', 'x^3'],
      ['0', '-x^2+x', '-x^4+x^3'],
      ['0', '0', '-x^5+x^4+x^3'],
    ]);
    expect(matEq(M, expected)).toBe(true);
    expect(matEq(U.mul(A), M)).toBe(true);
  });

  it('matches the GF(2^4) doctest', () => {
    const F16 = GFExtended(16n, { name: 'a' }) as any;
    const pF = new PolynomialRing(F16, 'x');
    const a = F16.gen();
    const zero = F16.zero();
    const one = F16.one();
    const c2: any[] = new Array(18).fill(zero);
    c2[0] = one;
    c2[17] = a;
    const c3: any[] = new Array(12).fill(zero);
    c3[0] = one;
    c3[7] = a.mul(a);
    c3[11] = a;
    const A16 = new Matrix(pF as any, 2, 2, [
      [pF.one(), pF.__call__(c2)],
      [pF.zero(), pF.__call__(c3)],
    ] as any);
    const M = A16.copy();
    const U = MP._weak_popov_form(M, { transformation: true })!;
    expect(matEq(U.mul(A16), M)).toBe(true);
    expect(MP.is_weak_popov(M)).toBe(true);
    // U is unimodular: its determinant is a nonzero constant
    const du = determinant(U) as any;
    expect(du.isZero()).toBe(false);
    expect(du.degree()).toBe(0);
  });
});

describe('reduced_form (pyx:2542-2669)', () => {
  it('matches the GF(2^3) doctest', () => {
    const F8 = GFExtended(8n, { name: 'a' }) as any;
    const pF = new PolynomialRing(F8, 'x');
    const a = F8.gen();
    const z = F8.zero();
    const o = F8.one();
    // A = [[x^2 + a, x^4 + a], [x^3, a*x^4]]
    const A = new Matrix(pF as any, 2, 2, [
      [pF.__call__([a, z, o]), pF.__call__([a, z, z, z, o])],
      [pF.__call__([z, z, z, o]), pF.__call__([z, z, z, z, a])],
    ] as any);
    const [W, U] = MP.reduced_form(A, { transformation: true });
    // (
    // [          x^2 + a           x^4 + a]  [1 0]
    // [x^3 + a*x^2 + a^2               a^2], [a 1]
    // )
    const a2 = a.mul(a);
    const Wexp = new Matrix(pF as any, 2, 2, [
      [pF.__call__([a, z, o]), pF.__call__([a, z, z, z, o])],
      [pF.__call__([a2, z, a, o]), pF.__call__([a2])],
    ] as any);
    const Uexp = new Matrix(pF as any, 2, 2, [
      [pF.one(), pF.zero()],
      [pF.__call__([a]), pF.one()],
    ] as any);
    expect(matEq(W, Wexp)).toBe(true);
    expect(matEq(U, Uexp)).toBe(true);
    expect(MP.is_reduced(W)).toBe(true);
    expect(matEq(U.mul(W), A)).toBe(true);
  });
});

// ============================================================================
// popov_form  (pyx:2345-2540)
// ============================================================================

describe('popov_form (pyx:2345-2540)', () => {
  const M = pmat(pR, [
    ['6*x+4', '5*x^3+5*x', '6*x^2+2*x+2'],
    ['4*x^2+5*x+2', 'x^4+5*x^2+2*x+4', '4*x^3+6*x^2+6*x+5'],
  ]);

  it('matches the transformation doctest', () => {
    const [P, U] = MP.popov_form(M, { transformation: true });
    expect(render(P)).toBe('[4 x^2 + 4*x + 1 3] / [0 4*x + 1 x^2 + 6*x + 1]');
    expect(render(U)).toBe('[x 2] / [5*x^2 + x + 6 3*x + 2]');
    expect(MP.is_popov(P)).toBe(true);
    expect(matEq(U.mul(M), P)).toBe(true);
  });

  it('matches the shifts and Hermite doctests', () => {
    const P = MP.popov_form(M, { shifts: [0, 2, 4] });
    expect(render(P)).toBe(
      '[4*x^2 + 3*x + 4 x^4 + 3*x^3 + 5*x^2 + 5*x + 5 0] / [6 5*x^2 + 6*x + 5 1]'
    );
    expect(MP.is_popov(P, { shifts: [0, 2, 4] })).toBe(true);
    expect(matEq(P, MP.popov_form(M, { shifts: [-6, -4, -2] }))).toBe(true);

    const dd = MP.row_degrees(M).reduce((a, b) => a + b, 0) + 1;
    expect(matEq(MP.popov_form(M, { shifts: [2 * dd, dd, 0] }), MP.hermite_form(M))).toBe(true);
  });

  it('is the transpose of the column-wise form of the transpose', () => {
    expect(
      matEq(MP.popov_form(M), MP.popov_form(M.transpose(), { row_wise: false }).transpose())
    ).toBe(true);
  });

  it('matches the discard-zero-vectors doctests', () => {
    expect(render(MP.popov_form(M, { row_wise: false }))).toBe('[x + 2 6 0] / [0 1 0]');
    const [P, U] = MP.popov_form(M, {
      transformation: true,
      row_wise: false,
      include_zero_vectors: false,
    });
    expect(render(P)).toBe('[x + 2 6] / [0 1]');
    expect(render(U)).toBe(
      '[3*x^2 + 6*x + 3 5*x^2 + 4*x + 4 3*x^3 + 3*x^2 + 2*x + 4] / [3 1 2*x + 1] / [5*x + 2 2 6]'
    );
    const MU = M.mul(U);
    expect(
      MU.get(0, 0).eq(P.get(0, 0)) &&
        MU.get(0, 1).eq(P.get(0, 1)) &&
        MU.get(1, 0).eq(P.get(1, 0)) &&
        MU.get(1, 1).eq(P.get(1, 1))
    ).toBe(true);
    expect(MU.get(0, 2).isZero() && MU.get(1, 2).isZero()).toBe(true);
  });
});

// ============================================================================
// hermite_form  (pyx:2671-2743)
// ============================================================================

describe('hermite_form (pyx:2671-2743)', () => {
  const A = pmat(pR, [
    ['x', '1', '2*x'],
    ['x', '1+x', '2'],
  ]);
  const B = pmat(pR, [
    ['x', '1', '2*x'],
    ['2*x', '2', '4*x'],
  ]);

  it('matches the base doctests', () => {
    expect(render(MP.hermite_form(A))).toBe('[x 1 2*x] / [0 x 5*x + 2]');
    const [H, U] = MP.hermite_form(A, { transformation: true });
    expect(render(H)).toBe('[x 1 2*x] / [0 x 5*x + 2]');
    expect(render(U)).toBe('[1 0] / [6 1]');
  });

  it('matches the rank-deficient doctests', () => {
    const [H1, U1] = MP.hermite_form(B, { transformation: true, include_zero_rows: false });
    expect([render(H1), render(U1)]).toEqual(['[x 1 2*x]', '[0 4]']);
    expect(matEq(U1.mul(B), H1)).toBe(true);

    const [H2, U2] = MP.hermite_form(B, { transformation: true, include_zero_rows: true });
    expect([render(H2), render(U2)]).toEqual(['[x 1 2*x] / [0 0 0]', '[0 4] / [5 1]']);
    expect(matEq(U2.mul(B), H2)).toBe(true);
  });
});

// ============================================================================
// approximant bases  (pyx:3392-3919)
// ============================================================================

const F97 = GF(97n) as any;
const R97 = new PolynomialRing(F97, 'x');

function fromCoeffs(pRing: any, cs: number[]) {
  return pRing.__call__(cs.map((c) => pRing.base_ring.__call__(c)));
}

describe('is_minimal_approximant_basis (pyx:3392-3590)', () => {
  // Storjohann's "Notes on computing minimal approximant bases" example.
  const order = 8;
  const shifts = [1, 1, 0, 0, 0];
  const pmat97 = new Matrix(R97 as any, 5, 1, [
    [fromCoeffs(R97, [35, 0, 41, 87, 3, 42, 22, 90])],
    [fromCoeffs(R97, [80, 15, 62, 87, 14, 93, 24, 0])],
    [fromCoeffs(R97, [42, 57, 90, 87, 22, 80, 71, 53])],
    [fromCoeffs(R97, [37, 72, 74, 6, 5, 75, 23, 47])],
    [fromCoeffs(R97, [36, 10, 74, 1, 29, 44, 87, 74])],
  ] as any);
  const appbas = pmat(R97, [
    ['x+47', '57', '58*x+44', '9*x+23', '93*x+76'],
    ['15', 'x+18', '52*x+23', '15*x+58', '93*x+88'],
    ['17', '86', 'x^2+77*x+16', '76*x+29', '90*x+78'],
    ['44', '36', '3*x+42', 'x^2+50*x+26', '85*x+44'],
    ['2', '22', '54*x+94', '73*x+24', 'x^2+2*x+25'],
  ]);

  it('accepts the Storjohann example basis', () => {
    expect(
      MP.is_minimal_approximant_basis(appbas, pmat97, order, {
        shifts,
        row_wise: true,
        normal_form: true,
      })
    ).toBe(true);
  });

  it('rejects x^8 * Id_5, which generates a strictly smaller module', () => {
    const M = new Matrix(R97 as any, 5, 5);
    const x8 = R97.gen().pow(8);
    for (let i = 0; i < 5; i++) M.set(i, i, x8);
    expect(MP.is_minimal_approximant_basis(M, pmat97, 8)).toBe(false);
  });

  it('accepts [x^8] column-wise', () => {
    const M = new Matrix(R97 as any, 1, 1, [[R97.gen().pow(8)]] as any);
    expect(
      MP.is_minimal_approximant_basis(M, pmat97, 8, { row_wise: false, normal_form: true })
    ).toBe(true);
  });

  it('raises on unsound input dimensions', () => {
    expect(() => MP.is_minimal_approximant_basis(appbas, pmat97, [8, 8], { shifts })).toThrowError(
      'order length should be the column dimension of the input matrix'
    );
    expect(() =>
      MP.is_minimal_approximant_basis(appbas, pmat97, order, { shifts, row_wise: false })
    ).toThrowError('shifts length should be the column dimension of the input matrix');
    const M = new Matrix(R97 as any, 1, 1, [[R97.gen().pow(8)]] as any);
    expect(() => MP.is_minimal_approximant_basis(M, pmat97, 8)).toThrowError(
      'column dimension should be the row dimension of the input matrix'
    );
  });
});

describe('minimal_approximant_basis (pyx:3593-3770)', () => {
  const order = [4, 3];
  const shifts = [-1, 2, 0];
  const F = pmat(pR, [
    ['5*x^3 + 4*x^2 + 4*x + 6', '5*x^2 + 4*x + 1'],
    ['2*x^2 + 2*x + 3', '6*x^2 + 6*x + 3'],
    ['4*x^3 + x + 1', '4*x^2 + 2*x + 3'],
  ]);

  it('matches the base doctests', () => {
    const P = MP.minimal_approximant_basis(F, order, { shifts });
    expect(MP.is_minimal_approximant_basis(P, F, order, { shifts })).toBe(true);
    expect(MP.is_minimal_approximant_basis(P, F, order, { shifts, normal_form: true })).toBe(false);
    const Pn = MP.minimal_approximant_basis(F, order, { shifts, normal_form: true });
    expect(MP.is_minimal_approximant_basis(Pn, F, order, { shifts, normal_form: true })).toBe(true);
  });

  it('defaults shifts to zero and accepts a single integer order', () => {
    expect(
      matEq(
        MP.minimal_approximant_basis(F, 3),
        MP.minimal_approximant_basis(F, [3, 3], { shifts: [0, 0, 0] })
      )
    ).toBe(true);
  });

  it('matches the column-wise doctests', () => {
    const P = MP.minimal_approximant_basis(F, [5, 2, 2], { shifts: [0, 1], row_wise: false });
    expect(
      MP.is_minimal_approximant_basis(P, F, [5, 2, 2], { shifts: [0, 1], row_wise: false })
    ).toBe(true);
    expect(
      matEq(
        MP.minimal_approximant_basis(F, 3, { row_wise: true }),
        MP.minimal_approximant_basis(F.transpose(), 3, { row_wise: false }).transpose()
      )
    ).toBe(true);
  });

  it('ignores columns/rows with nonpositive order', () => {
    const P = MP.minimal_approximant_basis(F, [4, 0, 3], { row_wise: false });
    expect(matEq(P, MP.minimal_approximant_basis(F, [4, -2, 3], { row_wise: false }))).toBe(true);
    const Fsub = pmat(pR, [
      ['5*x^3 + 4*x^2 + 4*x + 6', '5*x^2 + 4*x + 1'],
      ['4*x^3 + x + 1', '4*x^2 + 2*x + 3'],
    ]);
    expect(matEq(P, MP.minimal_approximant_basis(Fsub, [4, 3], { row_wise: false }))).toBe(true);
  });

  it('raises on unsound input dimensions', () => {
    expect(() => MP.minimal_approximant_basis(F, [4], { shifts })).toThrowError(
      'order length should be the column dimension'
    );
    expect(() => MP.minimal_approximant_basis(F, order, { shifts: [0, 0, 0, 0] })).toThrowError(
      'shifts length should be the row dimension'
    );
  });
});

describe('_approximant_basis_iterative (pyx:3772-3919)', () => {
  const pm = pmat(pR, [
    ['5*x^3 + 4*x^2 + 4*x + 6', '5*x^2', '3*x^2 + 4'],
    ['2*x^3 + 2*x^2 + 2*x + 3', 'x^3 + 6', '6*x + 3'],
  ]);

  it('matches the arbitrary-shift doctest', () => {
    const [P, rdeg] = MP._approximant_basis_iterative(pm, [4, 1, 2], [-3, 4]);
    expect(MP.is_minimal_approximant_basis(P, pm, [4, 1, 2], { shifts: [-3, 4] })).toBe(true);
    expect(rdeg).toEqual(MP.row_degrees(P, [-3, 4]));
  });

  it('supports zero and negative orders', () => {
    const [P, rdeg] = MP._approximant_basis_iterative(pm, [4, 0, 2], [3, -1]);
    expect(MP.is_minimal_approximant_basis(P, pm, [4, 0, 2], { shifts: [3, -1] })).toBe(true);
    expect(rdeg).toEqual(MP.row_degrees(P, [3, -1]));
    const [P2] = MP._approximant_basis_iterative(pm, [4, -3, 2], [3, -1]);
    expect(matEq(P, P2)).toBe(true);
  });

  it('returns the identity for the zero matrix', () => {
    const [P, rdeg] = MP._approximant_basis_iterative(pzero(pR, 3, 2), [2, 5], [5, 0, -4]);
    expect(rdeg).toEqual([5, 0, -4]);
    expect(render(P)).toBe('[1 0 0] / [0 1 0] / [0 0 1]');
  });
});

// ============================================================================
// Independent oracles
// ============================================================================

describe('exhaustive brute force for the form predicates', () => {
  // All 2x2 matrices over GF(2)[x] with entries of degree <= 1, tested against
  // the definitions of weak Popov / Popov re-derived from the class docstring.
  const F2 = GF(2n) as any;
  const pR2 = new PolynomialRing(F2, 'x');
  const entries: any[] = [[], [1], [0, 1], [1, 1]].map((c) =>
    pR2.__call__(c.map((v) => F2.__call__(v)))
  );

  function defLeadingPositions(A: any, shifts: number[] | null): [number[], number[]] {
    const lpos: number[] = [];
    const pdeg: number[] = [];
    for (let i = 0; i < A.nrows; i++) {
      let best = Number.NEGATIVE_INFINITY;
      let bestj = -1;
      for (let j = 0; j < A.ncols; j++) {
        const e = A.get(i, j);
        if (e.isZero()) continue;
        const d = e.degree() + (shifts ? shifts[j]! : 0);
        if (d >= best) {
          best = d;
          bestj = j;
        }
      }
      lpos.push(bestj);
      pdeg.push(bestj < 0 ? -1 : A.get(i, bestj).degree());
    }
    return [lpos, pdeg];
  }

  function defIsWeakPopov(A: any, shifts: number[] | null, ordered: boolean): boolean {
    const [lpos] = defLeadingPositions(A, shifts);
    const nz = lpos.filter((l) => l >= 0);
    if (new Set(nz).size !== nz.length) return false;
    if (ordered) {
      let seenZero = false;
      for (const l of lpos) {
        if (l < 0) seenZero = true;
        else if (seenZero) return false;
      }
      for (let i = 0; i + 1 < nz.length; i++) if (nz[i]! >= nz[i + 1]!) return false;
    }
    return true;
  }

  function defIsPopov(A: any, shifts: number[] | null): boolean {
    if (!defIsWeakPopov(A, shifts, true)) return false;
    const [lpos, pdeg] = defLeadingPositions(A, shifts);
    for (let i = 0; i < A.nrows; i++) {
      if (lpos[i]! < 0) continue;
      if (!A.get(i, lpos[i]!).is_monic()) return false;
      for (let k = 0; k < A.nrows; k++) {
        if (k === i) continue;
        if (A.get(k, lpos[i]!).degree() >= pdeg[i]!) return false;
      }
    }
    return true;
  }

  it('agrees with the definitions on all 256 matrices and 4 shifts', () => {
    let cases = 0;
    for (const a of entries)
      for (const b of entries)
        for (const c of entries)
          for (const d of entries) {
            const A = new Matrix(pR2 as any, 2, 2, [
              [a, b],
              [c, d],
            ] as any);
            for (const shifts of [null, [0, 1], [2, 0], [-1, 3]] as Array<number[] | null>) {
              cases++;
              expect(MP.is_weak_popov(A, { shifts })).toBe(defIsWeakPopov(A, shifts, false));
              expect(MP.is_weak_popov(A, { shifts, ordered: true })).toBe(
                defIsWeakPopov(A, shifts, true)
              );
              expect(MP.is_popov(A, { shifts })).toBe(defIsPopov(A, shifts));
            }
          }
    expect(cases).toBe(1024);
  });
});

describe('algebraic identities on random matrices', () => {
  const rng = makeRng(20260728);

  function randomPolyMatrix(F: any, ring: any, m: number, n: number, dmax: number) {
    const rows: any[][] = [];
    for (let i = 0; i < m; i++) {
      const row: any[] = [];
      for (let j = 0; j < n; j++) {
        const d = rng() % (dmax + 2);
        const cs: any[] = [];
        for (let k = 0; k < d; k++) cs.push(F.__call__(rng() % Number(F.order ?? 7n)));
        row.push(ring.__call__(cs));
      }
      rows.push(row);
    }
    return new Matrix(ring, m, n, rows as any);
  }

  it('weak_popov_form / popov_form / hermite_form satisfy U*A == form with U unimodular', () => {
    for (const p of [2n, 3n, 5n, 7n, 11n]) {
      const F = GF(p) as any;
      const ring = new PolynomialRing(F, 'x');
      for (let t = 0; t < 12; t++) {
        const m = 1 + (rng() % 4);
        const n = 1 + (rng() % 4);
        const A = randomPolyMatrix(F, ring, m, n, rng() % 4);
        const shifts = t % 3 === 0 ? null : Array.from({ length: n }, () => (rng() % 7) - 3);

        const [W, UW] = MP.weak_popov_form(A, { transformation: true, shifts });
        expect(matEq(UW.mul(A), W)).toBe(true);
        expect(MP.is_weak_popov(W, { shifts })).toBe(true);
        const dW = determinant(UW) as any;
        expect(dW.isZero()).toBe(false);
        expect(dW.degree()).toBe(0);
        expect(
          MP.is_weak_popov(MP.weak_popov_form(A, { shifts, ordered: true }), {
            shifts,
            ordered: true,
          })
        ).toBe(true);

        const [P, UP] = MP.popov_form(A, { transformation: true, shifts });
        expect(matEq(UP.mul(A), P)).toBe(true);
        expect(MP.is_popov(P, { shifts })).toBe(true);
        const dP = determinant(UP) as any;
        expect(dP.isZero()).toBe(false);
        expect(dP.degree()).toBe(0);
        // the Popov form is canonical: applying it again changes nothing
        expect(matEq(MP.popov_form(P, { shifts }), P)).toBe(true);
        // column-wise form of the transpose
        expect(
          matEq(MP.popov_form(A.transpose(), { shifts, row_wise: false }).transpose(), P)
        ).toBe(true);

        const [H, UH] = MP.hermite_form(A, { transformation: true });
        expect(matEq(UH.mul(A), H)).toBe(true);
        expect(MP.is_hermite(H)).toBe(true);
        const dH = determinant(UH) as any;
        expect(dH.isZero()).toBe(false);
        expect(dH.degree()).toBe(0);
        // pyx:2434: the Hermite form is the shifted Popov form for
        // shifts = ((n-1)d, ..., d, 0) with d larger than every degree
        const dd = MP.row_degrees(A).reduce((a, b) => a + Math.max(b, 0), 0) + 1;
        const hshift = Array.from({ length: n }, (_, j) => (n - 1 - j) * dd);
        expect(matEq(MP.popov_form(A, { shifts: hshift }), H)).toBe(true);
      }
    }
  });

  it('minimal_approximant_basis outputs are verified by is_minimal_approximant_basis', () => {
    for (const p of [2n, 5n, 7n, 11n]) {
      const F = GF(p) as any;
      const ring = new PolynomialRing(F, 'x');
      for (let t = 0; t < 10; t++) {
        const m = 1 + (rng() % 3);
        const n = 1 + (rng() % 3);
        const A = randomPolyMatrix(F, ring, m, n, rng() % 4);

        const order = Array.from({ length: n }, () => rng() % 5);
        const shifts = t % 2 === 0 ? null : Array.from({ length: m }, () => (rng() % 7) - 3);
        for (const normal_form of [false, true]) {
          const P = MP.minimal_approximant_basis(A, order, { shifts, normal_form });
          expect(MP.is_minimal_approximant_basis(P, A, order, { shifts, normal_form })).toBe(true);
        }

        const corder = Array.from({ length: m }, () => rng() % 5);
        const cshifts = Array.from({ length: n }, () => (rng() % 5) - 2);
        for (const normal_form of [false, true]) {
          const P = MP.minimal_approximant_basis(A, corder, {
            shifts: cshifts,
            row_wise: false,
            normal_form,
          });
          expect(
            MP.is_minimal_approximant_basis(P, A, corder, {
              shifts: cshifts,
              row_wise: false,
              normal_form,
            })
          ).toBe(true);
        }
      }
    }
  });
});

// ============================================================================
// Independent verification of krylov_kernel_basis's shifted-Popov claim
// ============================================================================

describe('krylov_kernel_basis returns shifts-Popov matrices (matrix2.pyx:20343-20478)', () => {
  it('confirms the Popov property of the polynomial forms in the matrix2 doctests', async () => {
    const { krylov_kernel_basis } = await import('./matrix_decompositions.js');
    const F = GF(97n) as any;
    const ffmat = (entries: number[][]) =>
      new Matrix(
        F,
        entries.length,
        entries[0]!.length,
        entries.map((r) => r.map((x) => F.__call__(x)))
      );
    const E = ffmat([
      [27, 49, 29],
      [50, 58, 0],
      [77, 10, 29],
    ]);
    const Mn = ffmat([
      [0, 1, 0],
      [0, 0, 1],
      [0, 0, 0],
    ]);

    // Sage documents each of these outputs as being in shifts-Popov form.
    const P = krylov_kernel_basis(E, Mn, undefined, undefined, false, 'x') as any;
    expect(MP.is_popov(P, { shifts: [0, 0, 0] })).toBe(true);

    const H = krylov_kernel_basis(E, Mn, [0, 3, 6], undefined, false, 'x') as any;
    expect(MP.is_popov(H, { shifts: [0, 3, 6] })).toBe(true);

    const Q = krylov_kernel_basis(E, Mn, [3, 0, 2], undefined, false, 'Y') as any;
    expect(MP.is_popov(Q, { shifts: [3, 0, 2] })).toBe(true);
  });
});

// ============================================================================
// Degenerate shapes
// ============================================================================

describe('zero and empty matrices', () => {
  const F5 = GF(5n) as any;
  const R5 = new PolynomialRing(F5, 'x');

  it('returns canonical forms and identity transformations', () => {
    for (const [m, n] of [
      [2, 3],
      [3, 2],
      [0, 0],
      [1, 1],
      [0, 3],
      [3, 0],
    ] as Array<[number, number]>) {
      const Z = pzero(R5, m, n);
      const [P, U] = MP.popov_form(Z, { transformation: true });
      expect(P.nrows).toBe(m);
      expect(P.ncols).toBe(n);
      expect(P.is_zero()).toBe(true);
      expect(matEq(U, MP.weak_popov_form(Z, { transformation: true })[1])).toBe(true);
      expect(matEq(U.mul(Z), P)).toBe(true);
      expect(MP.is_popov(P)).toBe(true);

      const [H, UH] = MP.hermite_form(Z, { transformation: true });
      expect(H.is_zero()).toBe(true);
      expect(matEq(UH.mul(Z), H)).toBe(true);
      expect(MP.is_hermite(H)).toBe(true);

      // an approximant basis for the zero matrix is the identity
      const A = MP.minimal_approximant_basis(Z, 3);
      expect(MP.is_minimal_approximant_basis(A, Z, 3)).toBe(true);
      expect(A.nrows).toBe(m);
      for (let i = 0; i < m; i++) {
        for (let j = 0; j < m; j++) {
          expect(A.get(i, j).eq(i === j ? R5.one() : R5.zero())).toBe(true);
        }
      }
    }
  });

  it('accepts an empty shift on a matrix with zero columns', () => {
    // upstream would raise from ``min([])`` here; the shift is irrelevant since
    // there are no columns, so the empty shift is accepted
    const Z = pzero(R5, 3, 0);
    expect(MP.weak_popov_form(Z, { shifts: [] }).nrows).toBe(3);
  });
});
