/**
 * @module sagemath-ts
 * @description TypeScript port of SageMath
 *
 * This is the main entry point for sagemath-ts.
 * It mirrors SageMath's top-level imports.
 */

// Type coercion utilities
export { type IntegerLike, type RationalLike, toBigInt, toRational } from './types/index.js';

// Ring structures
export { ZZ, Integer, IntegerRing } from './rings/index.js';
export { QQ, Q, RationalField } from './rings/index.js';
export { Rational } from './rings/index.js';
export { Zmod, Integers, IntegerModRing, IntegerMod, Mod } from './rings/index.js';

// Finite fields
export { GF, FiniteFieldPrime, FiniteFieldElement } from './rings/index.js';

// Arithmetic functions
export * as arith from './arith/index.js';

// Matrix operations
export {
  Matrix,
  MatrixSpace,
  MatrixFromEntries,
  matrix,
  zero_matrix,
  identity_matrix,
  diagonal_matrix,
  scalar_matrix,
  IntegerMatrix,
  IntegerMatrixFromEntries,
  zero_integer_matrix,
  identity_integer_matrix,
} from './matrix/index.js';

// Re-export commonly used functions at top level (like SageMath does)
export {
  gcd,
  GCD,
  lcm,
  LCM,
  xgcd,
  factor,
  is_prime,
  is_prime_power,
  next_prime,
  previous_prime,
  power_mod,
  inverse_mod,
  euler_phi,
  moebius,
  kronecker_symbol,
  legendre_symbol,
  jacobi_symbol,
  crt,
  CRT_list,
  CRT_basis,
  CRT_vectors,
  rational_reconstruction,
  half_gcd,
  continued_fraction,
  continued_fraction_value,
  convergents,
  isqrt,
  is_square,
  is_squarefree,
  squarefree_part,
  divisors,
  number_of_divisors,
  sigma,
  radical,
  prime_range,
  prime_factors,
  prime_divisors,
  valuation,
  sqrt_mod,
  formatFactorization,
  type Factorization,
} from './arith/index.js';

// Coding theory
export * as coding from './coding/index.js';

// Re-export commonly used coding functions at top level
export {
  ReedSolomonCode,
  DecodingError,
  createClassicalReedSolomonCode,
  createFRIReedSolomonCode,
} from './coding/index.js';

// Groups
export * as groups from './groups/index.js';

// Re-export commonly used group functions at top level
export {
  discrete_log as generic_discrete_log,
  bsgs,
  pohlig_hellman,
  order_from_multiple,
  multiple_of_order,
  has_order,
  multiples,
  multiple,
} from './groups/index.js';

// Cryptography
export * as crypto from './crypto/index.js';

// Free modules and lattices
export * as modules from './modules/index.js';

// Re-export commonly used module functions at top level
export {
  FreeModule,
  VectorSpace,
  span,
  vector,
  IntegerLattice,
} from './modules/index.js';

// Elliptic curves
export * as schemes from './schemes/index.js';
export {
  EllipticCurve,
  EllipticCurvePoint,
  EllipticCurveFiniteField,
  EllipticCurveGeneric,
  weil_pairing,
  tate_pairing,
  ate_pairing,
  discrete_log,
  abelian_group,
  embedding_degree,
  is_supersingular,
  is_ordinary,
} from './schemes/index.js';

// Re-export more module functions
export {
  gramSchmidt,
  lllReduce,
  isLLLReduced,
  genLattice,
  randomLattice,
} from './modules/index.js';

// Re-export crypto functions at top level
export {
  LWE,
  Regev,
  LindnerPeikert,
  UniformSampler,
  UniformNoiseLWE,
  samples,
  balance_sample,
} from './crypto/index.js';

// Statistics and distributions
export * as stats from './stats/index.js';

// Re-export commonly used distribution functions at top level
export {
  DiscreteGaussianDistributionIntegerSampler,
  DiscreteGaussianInteger,
  DiscreteGaussianDistributionLatticeSampler,
  DiscreteGaussianLattice,
  DiscreteGaussianDistributionPolynomialSampler,
  DiscreteGaussianPolynomial,
} from './stats/index.js';

// Errors
export * from './errors.js';
