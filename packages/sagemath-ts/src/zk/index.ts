/**
 * @module sage/zk
 * @description Zero-knowledge proof primitives and protocols
 *
 * This module provides implementations of fundamental ZK building blocks
 * including the sumcheck protocol and multilinear extensions.
 */

// Multilinear extension utilities
export {
  closestPowerOfTwo,
  intToBinary,
  binaryToInt,
  booleanHypercube,
  MAX_HYPERCUBE_DIM,
  eqPolynomial,
  multilinearExtension,
  sparseMultilinearExtension,
  isMultilinear,
  verifyMLEEvaluation,
} from './multilinear.js';

// Sumcheck protocol
export {
  type SumcheckProof,
  type SumcheckResult,
  sumcheckProve,
  sumcheckVerify,
  sumcheckRun,
  sumcheckRoundProver,
  sumcheckRoundVerifier,
  createPolyEvaluator,
} from './sumcheck.js';

// Polynomial commitment scheme helpers (KZG, FRI).
//
// Moved here from src/rings/polynomial/polynomial_commitment.ts: SageMath has no
// sage/rings/polynomial/polynomial_commitment.py, so the module does not belong
// in the mirrored sage tree. src/rings/polynomial/index.ts still re-exports
// every symbol below for backwards compatibility.
export {
  type FieldElement,
  type EvaluationPoint,
  // Core quotient operations (KZG)
  compute_quotient,
  computeQuotient,
  batch_quotient,
  batchQuotient,
  // Vanishing and interpolation
  compute_vanishing_polynomial,
  computeVanishingPolynomial,
  lagrange_interpolation,
  lagrangeInterpolation,
  // Linearization
  linearization,
  // Lagrange basis evaluation
  evaluate_basis,
  evaluateBasis,
  evaluate_basis_with_weights,
  barycentric_weights,
  barycentricWeights,
  // Degree operations
  degree_bound,
  degreeBound,
  check_degree_bound,
  checkDegreeBound,
  // Polynomial splitting
  split_poly,
  splitPoly,
  recombine_chunks,
  recombineChunks,
  // FRI operations
  compose_with_linear,
  composeWithLinear,
  fri_fold,
  friFold,
  // Utilities
  polynomial_derivative,
  polynomialDerivative,
  multi_evaluate,
  multiEvaluate,
  roots_of_unity_vanishing,
  rootsOfUnityVanishing,
  coset_vanishing,
  cosetVanishing,
  verify_quotient_proof,
  verifyQuotientProof,
  generate_powers,
  generatePowers,
} from './polynomial_commitment.js';
