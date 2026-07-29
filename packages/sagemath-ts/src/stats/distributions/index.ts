/**
 * @module sage/stats/distributions
 * @description Statistical distributions for sampling
 *
 * This module provides samplers for various discrete distributions commonly
 * used in lattice-based cryptography and related applications.
 *
 * @see Reference: sage/stats/distributions/
 */

// Discrete Gaussian over integers
export {
  DiscreteGaussianDistributionIntegerSampler,
  DiscreteGaussianInteger,
  statisticalDistance,
  klDivergence,
  sampleConvolution,
} from './discrete_gaussian_integer.js';

export type {
  DiscreteGaussianOptions,
  DiscreteGaussianAlgorithm,
  // The only way to pass a non-integer centre `c`, which is the public
  // `uniform+table` / `uniform+online` behaviour, so it belongs in the barrel.
  DiscreteGaussianOptionsInternal,
} from './discrete_gaussian_integer.js';

// Discrete Gaussian over lattices
export {
  DiscreteGaussianDistributionLatticeSampler,
  DiscreteGaussianLattice,
  DiscreteGaussianDistributionPolynomialSampler,
  DiscreteGaussianPolynomial,
  sampleShortVector,
  samplePreimage,
  _iter_vectors,
  // `_normalisation_factor_zz` returns a multiprecision real (SageMath's
  // RealField(prec) element), so its type is part of the public surface.
  RealNumberMP,
} from './discrete_gaussian_lattice.js';
// NOTE: `_mp` (test-only surface) and the module-local `qfrep` adapter over
// parigp-ts are deliberately NOT re-exported.

export type { DiscreteGaussianLatticeOptions } from './discrete_gaussian_lattice.js';
