/**
 * @module parigp-ts
 * @description TypeScript port of PARI/GP number theory functions
 *
 * This package provides TypeScript implementations of PARI/GP algorithms.
 * Reference: reference/pari/src/ and reference/cypari2/
 */

export const VERSION = '0.0.1';

// Core types
export {
  PariType,
  type GEN,
  type PariInt,
  type PariReal,
  type PariFfelt,
  type PariVec,
  type PariCol,
  type PariMat,
  // Type constructors
  mkInt,
  stoi,
  itos,
  mkFfeltFp,
  mkvec,
  mkcol,
  mkmat,
  // Type guards
  isInt,
  isReal,
  isFfelt,
  isVec,
  isCol,
  isMat,
  // Integer utilities
  signe,
  absi,
  abscmpii,
  equali1,
  is_pm1,
  isZero,
  mod2,
  mod4,
  mod8,
  vali,
  vals,
  // Constants
  gen_0,
  gen_1,
  gen_m1,
  gen_2,
} from './types.js';

// Finite field operations (Fp)
export {
  // Core Fp operations
  Fp_red,
  Fp_add,
  Fp_sub,
  Fp_neg,
  Fp_mul,
  Fp_sqr,
  Fp_inv,
  Fp_div,
  Fp_pow,
  // Square root and quadratic residues
  Fp_sqrt,
  Fp_issquare,
  kronecker,
  // Utility functions
  Fp_center,
  Fp_mulu,
  Fp_addmul,
  Fp_double,
  Fp_halve,
  Fp_eq,
  // GCD operations
  gcd,
  xgcd,
} from './ff.js';

// Elliptic curve point operations
export {
  // Types
  type EllipticPoint,
  type JacobianPoint,
  type ShortWeierstrassCurve,
  // Point at infinity
  ellinf,
  ell_is_inf,
  ellinf_FpJ,
  FpJ_is_inf,
  // Point creation
  mkpoint,
  // Point finding
  ellordinate,
  random_FpE,
  // Coordinate conversions
  FpE_to_FpJ,
  FpJ_to_FpE,
  // Curve membership
  FpE_isoncurve,
} from './elliptic/points.js';

// Elliptic curve initialization
export {
  ellinit,
  ellfromj,
  ellfromjFp,
  ellj,
  elldisc,
  ellcoeffs,
  ellisnonsingular,
  ellToShortWeierstrass,
  EllCurveType,
  EllipticCurveError,
  type EllipticCurve,
  type EllInitInput,
} from './elliptic/init.js';

// Elliptic curve group operations (ellcard, ellgroup, etc.)
export {
  // Types
  type EllipticCurveFp,
  type EllipticPointFp,
  // Point operations
  ellinf as ellinf_Fp,
  ellpoint,
  ell_is_inf as ell_is_inf_Fp,
  ellequal,
  FpE_neg,
  FpE_add,
  FpE_dbl,
  FpE_mul,
  FpE_random,
  // Curve operations
  ellinit_Fp,
  ellisoncurve,
  elllift_x,
  // Cardinality and trace
  ellcard,
  trace_of_frobenius,
  Fp_elltrace_naive,
  // Group structure
  ellgroup,
  ellgenerators,
  // Point order
  ellorder,
} from './elliptic/group.js';

// Elliptic curve point operations with Jacobian coordinates
// Source: FpE.c (Jacobian coordinate formulas)
export {
  // Jacobian coordinate operations
  FpJ_neg,
  FpJ_dbl,
  FpJ_add,
  // High-level point operations using points.ts types
  ellisoncurve as ellisoncurve_sw,
  ellneg,
  elladd,
  ellsub,
  ellmul,
} from './elliptic/point.js';

// Integer factorization
// Source: ifactor1.c, arith2.c
export {
  // Factorization
  Z_factor,
  factoru,
  formatFactorization,
  // Primality testing
  isPrime,
  // Types
  type Factorization,
} from './ifactor.js';

// Advanced elliptic curve functions (stubs)
// Source: elliptic.c, ellsea.c, ellisog.c, FpE.c
export {
  // Types
  type Isogeny,
  type IsogenyMap,
  // Discrete logarithm
  elllog,
  // Pairings
  elltatepairing,
  ellweilpairing,
  // SEA algorithm
  ellcard_sea,
  // Isogenies
  ellisogeny,
  ellisogenyapply,
  ellisogenycompose,
  // Frobenius
  ellfrobenius,
  // Internal/helper functions
  _FpE_Miller,
  ellembeddingdegree,
  // Division polynomials
  elldivpol,
  ellxn,
} from './elliptic/advanced.js';
