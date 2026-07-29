/**
 * @module sage/quadratic_forms
 * @description Quadratic forms
 */

export {
  BinaryQF,
  BinaryQuadraticForm,
  BinaryQF_reduced_representatives,
  class_number,
} from './binary_qf.js';

export {
  QuadraticForm,
  DiagonalQuadraticForm,
  quadratic_form_from_invariants,
  QFEvaluateVector,
  QFEvaluateMatrix,
  matrixQQ,
  identityQQ,
  toMatrixQQ,
  isRationalMatrix,
  determinantQQ,
  inverseQQ,
  adjugateQQ,
  type QFBaseRing,
  type MatrixLike,
  type VectorLike,
  type RationalMatrix,
} from './quadratic_form.js';

export {
  rational_diagonal_form,
  _rational_diagonal_form_and_transformation,
  signature,
  signature_vector,
  hasse_invariant,
  hasse_invariant__OMeara,
  is_hyperbolic,
  is_anisotropic,
  is_isotropic,
  anisotropic_primes,
  compute_definiteness,
  compute_definiteness_string,
  compute_definiteness_string_by_determinants,
  is_positive_definite,
  is_negative_definite,
  is_indefinite,
  is_definite,
  hilbert_symbol_QQ,
  is_padic_square,
  qfgaussred,
} from './quadratic_form__local_field_invariants.js';

export {
  TernaryQF,
  red_mfact,
  _reduced_ternary_form_eisenstein_with_matrix,
  _reduced_ternary_form_eisenstein_without_matrix,
  primitivize,
  evaluate,
  _find_zeros_mod_p_2,
  _find_zeros_mod_p_odd,
  primitive_zero_mod_p,
  extend,
  _basic_lemma,
  _basic_lemma_vec,
  _find_p_neighbor_from_vec,
  find_all_ternary_qf_by_level_disc,
  find_a_ternary_qf_by_level_disc,
  type TernaryCoefficients,
} from './ternary_qf.js';
