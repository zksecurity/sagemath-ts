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
} from './discrete_gaussian_integer.js';

// Discrete Gaussian over lattices
export {
  DiscreteGaussianDistributionLatticeSampler,
  DiscreteGaussianLattice,
  DiscreteGaussianDistributionPolynomialSampler,
  DiscreteGaussianPolynomial,
  sampleShortVector,
  samplePreimage,
} from './discrete_gaussian_lattice.js';

export type { DiscreteGaussianLatticeOptions } from './discrete_gaussian_lattice.js';
