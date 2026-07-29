/**
 * sagemath-ts side of the `matrix_ops` property-test area.
 *
 * Cases: tests/property/cases/matrix_ops.cases.json
 * SageMath counterpart: tests/property/python/areas/matrix_ops.py
 */

import { IntegerMatrix } from '../../../../packages/sagemath-ts/src/index.js';
import {
  LLL,
  elementary_divisors_integer,
  hermite_normal_form,
  rank_integer,
  smith_form_integer,
} from '../../../../packages/sagemath-ts/src/matrix/matrix_integer.js';

export const functions = {
  determinant_2x2: (a: bigint, b: bigint, c: bigint, d: bigint) => {
    const M = new IntegerMatrix(2, 2, [
      [a, b],
      [c, d],
    ]);
    return M.determinant().value;
  },
  determinant_3x3: (
    a11: bigint,
    a12: bigint,
    a13: bigint,
    a21: bigint,
    a22: bigint,
    a23: bigint,
    a31: bigint,
    a32: bigint,
    a33: bigint
  ) => {
    const M = new IntegerMatrix(3, 3, [
      [a11, a12, a13],
      [a21, a22, a23],
      [a31, a32, a33],
    ]);
    return M.determinant().value;
  },
  determinant_4x4: (...args: bigint[]) => {
    const entries: bigint[][] = [];
    for (let i = 0; i < 4; i++) {
      entries.push(args.slice(i * 4, (i + 1) * 4));
    }
    const M = new IntegerMatrix(4, 4, entries);
    return M.determinant().value;
  },
  rank_2x3: (a11: bigint, a12: bigint, a13: bigint, a21: bigint, a22: bigint, a23: bigint) => {
    const M = new IntegerMatrix(2, 3, [
      [a11, a12, a13],
      [a21, a22, a23],
    ]);
    return BigInt(rank_integer(M));
  },
  rank_3x3: (
    a11: bigint,
    a12: bigint,
    a13: bigint,
    a21: bigint,
    a22: bigint,
    a23: bigint,
    a31: bigint,
    a32: bigint,
    a33: bigint
  ) => {
    const M = new IntegerMatrix(3, 3, [
      [a11, a12, a13],
      [a21, a22, a23],
      [a31, a32, a33],
    ]);
    return BigInt(rank_integer(M));
  },
  hnf_2x2: (a: bigint, b: bigint, c: bigint, d: bigint) => {
    const M = new IntegerMatrix(2, 2, [
      [a, b],
      [c, d],
    ]);
    const H = hermite_normal_form(M);
    return [
      [H.get(0, 0).value, H.get(0, 1).value],
      [H.get(1, 0).value, H.get(1, 1).value],
    ];
  },
  hnf_3x3: (...args: bigint[]) => {
    const entries: bigint[][] = [];
    for (let i = 0; i < 3; i++) {
      entries.push(args.slice(i * 3, (i + 1) * 3));
    }
    const M = new IntegerMatrix(3, 3, entries);
    const H = hermite_normal_form(M);
    const result: bigint[][] = [];
    for (let i = 0; i < 3; i++) {
      result.push([H.get(i, 0).value, H.get(i, 1).value, H.get(i, 2).value]);
    }
    return result;
  },
  snf_2x2: (a: bigint, b: bigint, c: bigint, d: bigint) => {
    const M = new IntegerMatrix(2, 2, [
      [a, b],
      [c, d],
    ]);
    const result = smith_form_integer(M, true);
    const D = Array.isArray(result) ? result[0] : result;
    return [D.get(0, 0).value, D.get(1, 1).value];
  },
  snf_3x3: (...args: bigint[]) => {
    const entries: bigint[][] = [];
    for (let i = 0; i < 3; i++) {
      entries.push(args.slice(i * 3, (i + 1) * 3));
    }
    const M = new IntegerMatrix(3, 3, entries);
    const result = smith_form_integer(M, true);
    const D = Array.isArray(result) ? result[0] : result;
    return [D.get(0, 0).value, D.get(1, 1).value, D.get(2, 2).value];
  },
  lll_2x2: (a: bigint, b: bigint, c: bigint, d: bigint) => {
    const M = new IntegerMatrix(2, 2, [
      [a, b],
      [c, d],
    ]);
    const L = LLL(M);
    return [
      [L.get(0, 0).value, L.get(0, 1).value],
      [L.get(1, 0).value, L.get(1, 1).value],
    ];
  },
  lll_3x3: (...args: bigint[]) => {
    const entries: bigint[][] = [];
    for (let i = 0; i < 3; i++) {
      entries.push(args.slice(i * 3, (i + 1) * 3));
    }
    const M = new IntegerMatrix(3, 3, entries);
    const L = LLL(M);
    const result: bigint[][] = [];
    for (let i = 0; i < 3; i++) {
      result.push([L.get(i, 0).value, L.get(i, 1).value, L.get(i, 2).value]);
    }
    return result;
  },
  elementary_divisors_2x2: (a: bigint, b: bigint, c: bigint, d: bigint) => {
    const M = new IntegerMatrix(2, 2, [
      [a, b],
      [c, d],
    ]);
    const divs = elementary_divisors_integer(M);
    return divs.map((d) => d.value);
  },
  elementary_divisors_3x3: (...args: bigint[]) => {
    const entries: bigint[][] = [];
    for (let i = 0; i < 3; i++) {
      entries.push(args.slice(i * 3, (i + 1) * 3));
    }
    const M = new IntegerMatrix(3, 3, entries);
    const divs = elementary_divisors_integer(M);
    return divs.map((d) => d.value);
  },
};
