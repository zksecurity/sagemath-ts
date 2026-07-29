/**
 * @module sage/algebras/quatalg/quaternion_algebra_cython
 * @description Optimized code needed by quaternion algebras
 *
 * Port of: sage/algebras/quatalg/quaternion_algebra_cython.pyx
 * Reference: reference/sage/src/sage/algebras/quatalg/quaternion_algebra_cython.pyx
 *
 * These are the routines that convert between lists of rational quaternions
 * and matrices over `ZZ` / `QQ`.
 */

import { IntegerMatrix } from '../../matrix/matrix_integer.js';
import { Rational } from '../../rings/rational.js';
import { QQ } from '../../rings/rational_field.js';
import type { QuaternionAlgebra_ab } from './quaternion_algebra.js';
import {
  MatrixQQ,
  QuaternionAlgebraElement_rational_field,
  type RationalMatrix,
} from './quaternion_algebra_element.js';

function lcmBigInt(a: bigint, b: bigint): bigint {
  if (a === 0n || b === 0n) return 0n;
  let x = a < 0n ? -a : a;
  let y = b < 0n ? -b : b;
  const A = x;
  const B = y;
  while (y) {
    [x, y] = [y, x % y];
  }
  return (A / x) * B;
}

/**
 * Given a list of rational quaternions, return a matrix `A` over `ZZ` and a
 * denominator `d` such that the rows of `(1/d)A` are the coefficients of the
 * quaternions.
 *
 * @param v - list of quaternions in a rational quaternion algebra
 * @param reverse - whether the order of the coordinates as well as the order
 *   of the list `v` should be reversed
 *
 * @example
 * ```typescript
 * // A.<i,j,k> = QuaternionAlgebra(-4,-5)
 * // integral_matrix_and_denom_from_rational_quaternions([i/2, 1/3+j+k])
 * // ([[0, 3, 0, 0], [2, 0, 6, 6]], 6)
 * ```
 *
 * @see Reference: quaternion_algebra_cython.pyx:46
 *   (integral_matrix_and_denom_from_rational_quaternions)
 */
export function integral_matrix_and_denom_from_rational_quaternions(
  v: readonly QuaternionAlgebraElement_rational_field[],
  reverse: boolean = false
): [IntegerMatrix, bigint] {
  const n = v.length;
  const A = new IntegerMatrix(n, 4);
  if (n === 0) {
    return [A, 1n];
  }

  // Least common multiple of the denominators
  let d = 1n;
  for (const x of v) {
    d = lcmBigInt(d, x.denominator());
  }

  for (let i = 0; i < n; i++) {
    const x = v[i] as QuaternionAlgebraElement_rational_field;
    const q = d / x.denominator();
    const [, cx, cy, cz, cw] = x.denominator_and_integer_coefficient_tuple();
    if (reverse) {
      A.set(n - i - 1, 3, q * cx);
      A.set(n - i - 1, 2, q * cy);
      A.set(n - i - 1, 1, q * cz);
      A.set(n - i - 1, 0, q * cw);
    } else {
      A.set(i, 0, q * cx);
      A.set(i, 1, q * cy);
      A.set(i, 2, q * cz);
      A.set(i, 3, q * cw);
    }
  }
  return [A, d];
}

/**
 * Return the matrix over `QQ` whose rows have entries the coefficients of the
 * rational quaternions in `v`.
 *
 * @param v - list of quaternions in a rational quaternion algebra
 * @param reverse - whether the order of the coordinates as well as the order
 *   of the list `v` should be reversed
 *
 * @see Reference: quaternion_algebra_cython.pyx:126
 *   (rational_matrix_from_rational_quaternions)
 */
export function rational_matrix_from_rational_quaternions(
  v: readonly QuaternionAlgebraElement_rational_field[],
  reverse: boolean = false
): RationalMatrix {
  const n = v.length;
  const A = new MatrixQQ(QQ, n, 4);
  for (let i = 0; i < n; i++) {
    const t = (v[i] as QuaternionAlgebraElement_rational_field).coefficient_tuple();
    for (let j = 0; j < 4; j++) {
      if (reverse) {
        A.set(n - i - 1, 3 - j, t[j] as Rational);
      } else {
        A.set(i, j, t[j] as Rational);
      }
    }
  }
  return A;
}

/**
 * Given an integral matrix and a denominator, return the list of rational
 * quaternions given by its rows divided by `d`.
 *
 * @param A - rational quaternion algebra
 * @param H - matrix over the integers
 * @param d - the denominator
 * @param reverse - whether the order of the coordinates as well as the order
 *   of the list should be reversed
 *
 * @see Reference: quaternion_algebra_cython.pyx:184
 *   (rational_quaternions_from_integral_matrix_and_denom)
 */
export function rational_quaternions_from_integral_matrix_and_denom(
  A: QuaternionAlgebra_ab,
  H: IntegerMatrix,
  d: bigint,
  reverse: boolean = false
): QuaternionAlgebraElement_rational_field[] {
  const v: QuaternionAlgebraElement_rational_field[] = [];
  const idx: number[] = [];
  if (reverse) {
    for (let i = H.nrows - 1; i >= 0; i--) idx.push(i);
  } else {
    for (let i = 0; i < H.nrows; i++) idx.push(i);
  }
  for (const i of idx) {
    const coeffs: Rational[] = [];
    for (let j = 0; j < 4; j++) {
      const entry = reverse ? H.get(i, 3 - j).value : H.get(i, j).value;
      coeffs.push(new Rational(entry, d));
    }
    v.push(new QuaternionAlgebraElement_rational_field(A, coeffs));
  }
  return v;
}
