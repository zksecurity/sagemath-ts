/**
 * sagemath-ts side of the `matrix` property-test area.
 *
 * Matrix results are transported as nested bigint arrays so the transcript
 * tests their entries without depending on either implementation's printer.
 *
 * Cases: tests/property/cases/matrix.cases.json
 * SageMath counterpart: tests/property/python/areas/matrix.py
 */

import {
  type IntegerMatrix,
  IntegerMatrixFromEntries,
  LLL,
  hermite_normal_form,
  smith_form_integer,
} from '../../../../packages/sagemath-ts/src/matrix/index.js';

function matrixRows(matrix: IntegerMatrix): bigint[][] {
  const rows: bigint[][] = [];
  for (let i = 0; i < matrix.nrows; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < matrix.ncols; j++) {
      row.push(matrix.get(i, j).value);
    }
    rows.push(row);
  }
  return rows;
}

function matrix2(a: bigint, b: bigint, c: bigint, d: bigint): IntegerMatrix {
  return IntegerMatrixFromEntries([
    [a, b],
    [c, d],
  ]);
}

function matrix3(
  a: bigint,
  b: bigint,
  c: bigint,
  d: bigint,
  e: bigint,
  f: bigint,
  g: bigint,
  h: bigint,
  i: bigint
): IntegerMatrix {
  return IntegerMatrixFromEntries([
    [a, b, c],
    [d, e, f],
    [g, h, i],
  ]);
}

export const functions = {
  determinant_2x2: (a: bigint, b: bigint, c: bigint, d: bigint): bigint =>
    matrix2(a, b, c, d).determinant().value,

  determinant_3x3: (
    a: bigint,
    b: bigint,
    c: bigint,
    d: bigint,
    e: bigint,
    f: bigint,
    g: bigint,
    h: bigint,
    i: bigint
  ): bigint => matrix3(a, b, c, d, e, f, g, h, i).determinant().value,

  rank_2x3: (a: bigint, b: bigint, c: bigint, d: bigint, e: bigint, f: bigint): number =>
    IntegerMatrixFromEntries([
      [a, b, c],
      [d, e, f],
    ]).rank(),

  rank_3x3: (
    a: bigint,
    b: bigint,
    c: bigint,
    d: bigint,
    e: bigint,
    f: bigint,
    g: bigint,
    h: bigint,
    i: bigint
  ): number => matrix3(a, b, c, d, e, f, g, h, i).rank(),

  hnf_2x2: (a: bigint, b: bigint, c: bigint, d: bigint): bigint[][] =>
    matrixRows(hermite_normal_form(matrix2(a, b, c, d)) as IntegerMatrix),

  hnf_3x3: (
    a: bigint,
    b: bigint,
    c: bigint,
    d: bigint,
    e: bigint,
    f: bigint,
    g: bigint,
    h: bigint,
    i: bigint
  ): bigint[][] =>
    matrixRows(hermite_normal_form(matrix3(a, b, c, d, e, f, g, h, i)) as IntegerMatrix),

  snf_2x2: (a: bigint, b: bigint, c: bigint, d: bigint): bigint[][] =>
    matrixRows(smith_form_integer(matrix2(a, b, c, d)) as IntegerMatrix),

  snf_3x3: (
    a: bigint,
    b: bigint,
    c: bigint,
    d: bigint,
    e: bigint,
    f: bigint,
    g: bigint,
    h: bigint,
    i: bigint
  ): bigint[][] =>
    matrixRows(smith_form_integer(matrix3(a, b, c, d, e, f, g, h, i)) as IntegerMatrix),

  lll_2x2: (a: bigint, b: bigint, c: bigint, d: bigint): bigint[][] =>
    matrixRows(LLL(matrix2(a, b, c, d)) as IntegerMatrix),

  lll_3x3: (
    a: bigint,
    b: bigint,
    c: bigint,
    d: bigint,
    e: bigint,
    f: bigint,
    g: bigint,
    h: bigint,
    i: bigint
  ): bigint[][] => matrixRows(LLL(matrix3(a, b, c, d, e, f, g, h, i)) as IntegerMatrix),
};
