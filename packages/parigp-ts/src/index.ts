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

// Linear algebra over Z/dZ for arbitrary (composite) d: Howell normal form.
// Source: bb_hnf.c
export {
  // Matrix layout helpers (PARI column-major)
  type ZM,
  zm_from_rows,
  zm_to_rows,
  // Kernel / image / determinant / inverse mod d
  matkermod,
  matkermod_basis,
  matimagemod,
  matdetmod,
  matinvmod,
  // Shared PARI error kinds
  PariTypeError,
  PariDomainError,
  PariDimError,
  PariInvError,
  PariPrimeError,
  PariSqrtnError,
  PariFlagError,
} from './matkermod.js';

// Irreducible polynomials over F_p (Adleman-Lenstra).
// Source: polarit3.c, subcyclo.c
export {
  // FpX arithmetic
  type FpX,
  type FpXY,
  FpX_renormalize,
  FpX_red,
  FpX_degree,
  FpX_add,
  FpX_sub,
  FpX_neg,
  FpX_mul,
  FpX_Fp_mul,
  FpX_divrem,
  FpX_rem,
  FpX_normalize,
  FpX_gcd,
  FpXQ_mul,
  FpXQ_pow,
  FpX_is_irred,
  pol_xn,
  // Resultants / composed sums
  FpX_FpXY_resultant,
  FpX_composedsum,
  FpXV_composedsum,
  // Cyclotomic subfields and Artin-Schreier towers
  polsubcyclo_prime,
  fpinit_check,
  fpinit,
  ffinit_Artin_Schreier,
  ffinit_Artin_Schreier_2,
  // Main entry points
  ffinit,
  init_Fq,
  ffinit_rand,
  ffnbirred,
} from './ffinit.js';

// Binary quadratic forms: composition, reduction, representation.
// Source: Qfb.c, quad.c
export {
  // Type and constructors
  type Qfb as QfbForm,
  Qfb,
  mkqfb,
  qfb_disc,
  qfb_disc3,
  qfb_is_qfi,
  qfb_equal,
  qfb_1,
  qfbinv,
  qfb_apply,
  // Reduction
  qfbred,
  qfbredsl2,
  qfi_rho,
  qfi_red_fast,
  qfbred_withLimit,
  qfbredsl2_withLimit,
  // Composition and powering
  qfbcompraw,
  qfbcomp,
  qfbsqr,
  qfbsqrraw,
  qfbpow,
  qfbpowraw,
  // Prime forms and representation
  primeform,
  qfbsolve,
  cornacchia,
  cornacchia2,
  qfbcornacchia,
  // Supporting number theory
  sqrti,
  Z_issquareall,
  Z2_sqrt,
  Zp_sqrt,
  Zn_quad_roots,
} from './qfb.js';
