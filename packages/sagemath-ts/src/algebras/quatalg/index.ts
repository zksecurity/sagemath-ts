/**
 * @module sage/algebras/quatalg
 * @description Quaternion algebras (`sage.algebras.quatalg`)
 *
 * Port of: sage/algebras/quatalg/all.py
 */

export {
  INFINITE_PLACE_QQ,
  QuaternionAlgebra,
  QuaternionAlgebra_ab,
  QuaternionFractionalIdeal_rational,
  QuaternionOrder,
  ZZLattice,
  basis_for_quaternion_lattice,
  hilbert_symbol_QQ,
  intersection_of_row_modules_over_ZZ,
  maxord_solve_aux_eq,
  normalize_basis_at_p,
  type IntegerPowerSeries,
} from './quaternion_algebra.js';

export {
  integral_matrix_and_denom_from_rational_quaternions,
  rational_matrix_from_rational_quaternions,
  rational_quaternions_from_integral_matrix_and_denom,
} from './quaternion_algebra_cython.js';

export {
  MatrixQQ,
  QuaternionAlgebraElement_generic,
  QuaternionAlgebraElement_number_field,
  QuaternionAlgebraElement_rational_field,
  type QuaternionAlgebraElement,
  type QuaternionCoefficient,
  type QuaternionLike,
  type RationalMatrix,
  type RationalPolynomial,
} from './quaternion_algebra_element.js';
