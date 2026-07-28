/**
 * @module sage/crypto
 * @description Cryptographic functions and classes
 *
 * This module provides cryptographic primitives including:
 * - Hard lattice generators
 * - Learning with Errors (LWE) oracles
 * - Ring-LWE oracles
 * - S-Box cryptographic analysis
 * - Boolean functions for cryptographic analysis
 *
 * @see Reference: sage/crypto/
 */

// Lattice generation
export { gen_lattice } from './lattice.js';
export type { GenLatticeOptions, LatticeType } from './lattice.js';

// LWE (Learning with Errors)
export {
  // Samplers
  UniformSampler,
  UniformPolynomialSampler,
  // LWE classes
  LWE,
  LWEVector,
  Regev,
  LindnerPeikert,
  UniformNoiseLWE,
  // Ring-LWE classes
  RingLWE,
  RingLindnerPeikert,
  RingLWEConverter,
  // Functions
  samples,
  balance_sample,
} from './lwe.js';

export type {
  SecretDistribution,
  UniformNoiseLWEInstance,
  OracleKeywords,
  SamplesOptions,
  LWESample,
  LWEOracle,
  DistributionSampler,
} from './lwe.js';

// Re-export actual discrete Gaussian implementations from stats module
export {
  DiscreteGaussianDistributionIntegerSampler,
  DiscreteGaussianInteger,
  DiscreteGaussianDistributionLatticeSampler,
  DiscreteGaussianLattice,
  DiscreteGaussianDistributionPolynomialSampler,
  DiscreteGaussianPolynomial,
} from '../stats/index.js';

// S-Box cryptographic analysis
export {
  SBox,
  // Standard S-boxes
  AES_SBOX,
  PRESENT_SBOX,
  DES_SBOX1,
  GIFT_SBOX,
  SKINNY_SBOX,
  // Construction functions
  feistel_construction,
  misty_construction,
} from './sbox.js';

export type { LATScale } from './sbox.js';

// Boolean functions
export {
  BooleanFunction,
  hammingWeight,
  walshHadamardInPlace,
  randomBooleanFunction,
  fromANF,
  createBentFunction,
  createAffineFunction,
  verifyParseval,
} from './boolean_function.js';

export type { TruthTableInput } from './boolean_function.js';
