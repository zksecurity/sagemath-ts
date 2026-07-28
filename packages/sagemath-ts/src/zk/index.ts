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
