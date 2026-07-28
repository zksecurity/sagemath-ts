/**
 * @module sage/modules
 * @description Free modules and lattices
 *
 * This module provides:
 * - Free modules over commutative rings
 * - Vector spaces over fields
 * - Integer lattices with LLL/BKZ reduction
 * - Tools for lattice-based cryptography
 */

// Free module elements (vectors)
export {
  FreeModuleElement,
  FreeModuleElementDense,
  FreeModuleElementSparse,
  vector,
  zeroVector,
  randomVector,
} from './free_module_element.js';

// Free modules and vector spaces
export {
  // Factory functions
  FreeModule,
  VectorSpace,
  span,
  // Options
  type FreeModuleOptions,
  // Base classes
  ModuleFreeAmbient,
  FreeModuleGeneric,
  FreeModuleDomain,
  FreeModulePID,
  FreeModuleField,
  // Ambient modules
  FreeModuleAmbient,
  FreeModuleAmbientPID,
  FreeModuleAmbientField,
  // Submodules with basis
  FreeModuleWithBasis,
  FreeModuleSubspaceWithBasis,
  // Submodules given by generators (basis = echelon form of the generators)
  FreeModuleSubmodule,
  FreeModuleSubmodulePID,
  FreeModuleSubspace,
  // Quotients
  FreeModuleQuotient,
  // Fraction field of a Euclidean base ring (e.g. QQ(x) for QQ[x])
  FractionFieldElement,
} from './free_module.js';

// Integer lattices
export {
  // Factory function
  IntegerLattice,
  // Main lattice class
  FreeModuleIntegerLattice,
  // Gram-Schmidt
  gramSchmidt,
  type GramSchmidtResult,
  // LLL reduction
  lllReduce,
  isLLLReduced,
  type LLLReduceOptions,
  // Lattice generation
  genLattice,
  standardLattice,
  randomLattice,
  qaryLattice,
  qaryDualLattice,
  // Cryptographic utilities
  discreteGaussianSample,
  smoothingParameter,
  hermiteFactor,
  estimateBKZBlockSize,
  // Options types
  type LLLOptions,
  type BKZOptions as BKZOptionsLattice,
  type ShortestVectorOptions,
  type ClosestVectorOptions,
  type GenLatticeOptions,
} from './free_module_integer.js';

// BKZ reduction
export {
  // Main BKZ function
  BKZ,
  BKZWithInfo,
  // HKZ reduction (BKZ with blockSize = rank)
  HKZ,
  // Verification
  isBKZReduced,
  // Utilities
  computeHermiteFactor,
  estimateBlockSize,
  // Types
  type BKZOptions,
  type BKZResult,
} from './bkz.js';
