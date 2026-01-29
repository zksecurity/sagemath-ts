/**
 * @module sage/matrix
 * @description Matrix operations and matrix spaces
 *
 * Port of: sage/matrix/
 */

// Generic matrix class
export {
  Matrix,
  zero_matrix,
  identity_matrix,
  diagonal_matrix,
  scalar_matrix,
} from './matrix_generic.js';

// Matrix space factory
export {
  MatrixSpace,
  MatrixSpaceClass,
  MatrixFromEntries,
  matrix,
} from './matrix_space.js';

// Integer matrices
export {
  IntegerMatrix,
  IntegerMatrixFromEntries,
  zero_integer_matrix,
  identity_integer_matrix,
  // Normal forms
  hermite_normal_form,
  smith_form_integer,
  elementary_divisors_integer,
  invariant_factors_integer,
  rank_integer,
  kernel_matrix,
  left_kernel_matrix as left_kernel_matrix_integer,
  cyclic_decomposition_integer,
  // Integer-specific operations
  frobenius_form_integer,
  saturation,
  index_in_saturation,
  pivots_integer,
  height,
  gcd_integer_matrix,
  is_primitive,
  LLL,
  BKZ,
  is_LLL_reduced,
  rational_reconstruction,
  symplectic_form_integer,
  p_minimal_polynomials,
  null_ideal,
  integer_valued_polynomials_generators,
  antitranspose,
  insert_row,
  augment_integer,
  delete_zero_columns,
  factor_out_common_factors_from_each_row,
} from './matrix_integer.js';

// Matrix operations (determinant, rank, kernel, etc.)
export {
  // Determinant and related
  determinant,
  det,
  quantum_determinant,
  pfaffian,
  permanent,
  permanental_minor,
  minors,
  // Inverse
  inverse,
  pseudoinverse,
  adjugate,
  adjoint_classical,
  // Rank and nullity
  rank,
  left_nullity,
  nullity,
  right_nullity,
  // Kernel and image
  right_kernel_matrix,
  left_kernel_matrix,
  right_kernel,
  left_kernel,
  kernel,
  image,
  row_space,
  column_space,
  // Linear system solving
  solve_left,
  solve_right,
  // Characteristic and minimal polynomials
  charpoly,
  characteristic_polynomial,
  minpoly,
  minimal_polynomial,
  // Eigenvalues and eigenvectors
  eigenvalues,
  singular_values,
  eigenvectors_left,
  eigenvectors_right,
  eigenspaces_left,
  eigenspaces_right,
  eigenmatrix_left,
  eigenmatrix_right,
  // Type
  type Eigenspace,
  // Matrix properties
  is_one,
  is_scalar,
  is_symmetric,
  is_diagonal,
  is_triangular,
  is_unitary,
  is_normal,
  is_nilpotent,
  is_diagonalizable,
  is_semisimple,
  is_positive_definite,
  is_positive_semidefinite,
  is_similar,
  // Norms
  norm,
  density,
  // Conjugate
  conjugate,
  conjugate_transpose,
  H,
  T,
  C,
  // Matrix construction operations
  direct_sum,
  stack,
  augment,
} from './matrix_operations.js';

// Matrix decompositions
export {
  // Echelon forms
  echelonize,
  echelon_form,
  rref,
  extended_echelon_form,
  pivot_rows,
  // LU decomposition
  LU,
  // QR decomposition
  QR,
  gram_schmidt_noscale,
  gram_schmidt,
  // Cholesky
  cholesky,
  inverse_positive_definite,
  // LDL
  indefinite_factorization,
  block_ldlt,
  // Smith normal form
  smith_form,
  elementary_divisors,
  // Hermite normal form
  hermite_form,
  // Hessenberg form
  hessenberg_form,
  hessenbergize,
  // Jordan and rational forms
  jordan_form,
  jordan_decomposition,
  diagonalization,
  rational_form,
  zigzag_form,
  // Symplectic form
  symplectic_form,
  // LLL
  LLL_gram,
  // Square root
  principal_square_root,
  // Exponential
  exp,
  // Decomposition
  decomposition,
  cyclic_subspace,
  maxspin,
  wiedemann,
  // Krylov methods
  krylov_matrix,
  krylov_basis,
  krylov_kernel_basis,
} from './matrix_decompositions.js';

// Special matrix constructors
export {
  // Random matrices
  random_matrix,
  random_rref_matrix,
  random_echelonizable_matrix,
  random_subspaces_matrix,
  random_unimodular_matrix,
  random_unitary_matrix,
  random_diagonalizable_matrix,
  // Special constructors
  column_matrix,
  ones_matrix,
  lehmer,
  elementary_matrix,
  circulant,
  block_matrix,
  block_diagonal_matrix,
  jordan_block,
  companion_matrix,
  hilbert,
  vandermonde,
  toeplitz,
  hankel,
  vector_on_axis_rotation_matrix,
  ith_to_zero_rotation_matrix,
  // Berlekamp-Massey
  berlekamp_massey,
  // Matrix operations
  tensor_product,
  elementwise_product,
  rook_vector,
  prod_of_row_sums,
  hadamard_bound,
  subdivide,
  subdivision,
  subdivisions,
  set_block,
  find,
  apply_map,
  apply_morphism,
  subs,
  derivative,
  numerical_approx,
  denominator,
  randomize,
  get_bandwidth,
  as_bipartite_graph,
  automorphisms_of_rows_and_columns,
  permutation_normal_form,
  is_permutation_of,
} from './matrix_special.js';

// GF(2) matrices
export {
  Matrix_mod2_dense,
  from_png,
  from_png_data,
  to_png,
  to_png_data,
  pluq,
  ple,
  zero_matrix_gf2,
  identity_matrix_gf2,
  random_matrix_gf2,
  matrix_gf2_from_entries,
} from './matrix_mod2.js';

// Z/nZ matrices
export {
  Matrix_modn_dense,
  zero_matrix_modn,
  identity_matrix_modn,
  random_matrix_modn,
  matrix_modn_from_entries,
} from './matrix_modn.js';

// Real matrix decompositions (IEEE 754 doubles)
export {
  // SVD
  SVD_double,
  SVD_reconstruct,
  pseudoinverse_SVD,
  rank_SVD,
  condition_number_SVD,
  frobenius_norm_SVD,
  spectral_norm_SVD,
  low_rank_approx_SVD,
  // QR (Householder)
  QR_double,
  solve_QR,
  // LU (with partial pivoting)
  LU_double,
  solve_LU,
  det_LU,
  inverse_LU,
  // Types
  type SVDResult,
  type QRResult,
  type LUResult,
} from './matrix_decompositions_additions.js';
