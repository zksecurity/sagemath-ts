/**
 * @module sage/rings/polynomial
 * @description Polynomial rings and related structures
 */

// From polynomial_element.ts
export {
  Polynomial,
  type RingElement,
  type CoefficientRing,
  type PolynomialRingBase,
} from './polynomial_element.js';

// From polynomial_ring.ts
export { PolynomialRing, PolynomialRingConstructor } from './polynomial_ring.js';

// From quotient_ring.ts
export { QuotientRing, QuotientRingElement } from './quotient_ring.js';

// From multi_polynomial_element.ts
export {
  MPolynomial,
  type MPolynomialRingBase,
  type TermOrder,
  type Exponent,
  exponentToKey,
  keyToExponent,
  totalDegree,
  addExponents,
  subtractExponents,
  exponentDivides,
  getTermOrderComparator,
} from './multi_polynomial_element.js';

// From multi_polynomial_ring.ts
export { MPolynomialRing, MPolynomialRingConstructor } from './multi_polynomial_ring.js';

// From zk/polynomial_commitment.ts - ZK polynomial commitment operations.
//
// This module has NO SageMath counterpart (there is no
// sage/rings/polynomial/polynomial_commitment.py), so it now lives under
// src/zk/ next to sumcheck and multilinear. These re-exports exist purely for
// backwards compatibility: the symbols used to be defined in this directory and
// are part of the published `./rings` / `./rings/polynomial` surface.
// New code should import them from `src/zk/index.ts` instead.
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
} from '../../zk/polynomial_commitment.js';

// From convolution.ts - FFT/NTT operations
export {
  // NTT operations
  ntt,
  intt,
  ntt_multiply,
  // Generic FFT
  fft,
  ifft,
  fft_multiply,
  type FFTElement,
  // Multi-point evaluation/interpolation
  evaluate_on_domain,
  interpolate_from_domain,
  // Root of unity utilities
  find_primitive_root,
  max_fft_size,
  // Naive implementations for testing
  dft_naive,
  idft_naive,
  convolve_naive,
  // Prime utilities
  supports_ntt_size,
  find_ntt_prime,
  NTT_FRIENDLY_PRIMES,
} from './convolution.js';

// From polynomial_gf2x.ts - GF(2)[x] with bit-packed storage
export {
  GF2X,
  GF2XRing,
  GF2X_Ring,
  buildIrred,
  buildSparseIrred,
  buildRandomIrred,
  GF2X_BuildIrred_list,
  GF2X_BuildSparseIrred_list,
  GF2X_BuildRandomIrred_list,
  squareFreeDecomp,
  distinctDegreeFactorization,
  equalDegreeFactorization,
  factor as factorGF2X,
} from './polynomial_gf2x.js';
