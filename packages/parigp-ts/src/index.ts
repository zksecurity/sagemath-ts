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
  // Extended t_QFB: a form carrying Shanks' logarithmic distance (Qfb.c:111-123)
  type QfbExt,
  type QfbLike,
  is_qfbext,
  // qfr3 / qfr5 containers and the distance (Qfb.c:396-430, 552-558)
  type Qfr3,
  type Qfr5,
  type QfrData,
  qfr5_dist,
  qfr_data_init,
  // PARI's t_REAL kernel (kernel/none/{mp.c,mp_indep.c,add.c}, basemath/trans1.c).
  // This is the package-level t_REAL kernel; see the note on './buch.js' below.
  type MpReal,
  nbits2prec,
  realprec,
  precision,
  expo,
  realsigne,
  real_0_bit,
  real_0,
  real_1,
  negr,
  absr,
  shiftr,
  setexpo,
  gequal1,
  itor,
  rtor,
  truncr,
  gcvtoi,
  addrr,
  subrr,
  addir,
  subir,
  addrs,
  subrs,
  mulrr,
  sqrr,
  mulir,
  mulri,
  mulsr,
  mulrs,
  divrr,
  divir,
  divri,
  divru,
  sqrtr_abs,
  sqrtr,
  mplog2,
  logr_abs,
  mpreal_to_frac,
} from './qfb.js';

// Multiple polynomial quadratic sieve
// Source: mpqs.c
export { mpqs, type MpqsOptions } from './mpqs.js';

// Modular polynomials Phi_L(X, Y) and Hilbert class polynomials.
// Source: polmodular.c, polclass.c, volcano.c
export {
  // Errors
  PariBugError,
  PariImplError,
  PariPriorityError,
  PariArchError,
  // Class invariants
  INV_J,
  INV_F,
  INV_F2,
  INV_F3,
  INV_F4,
  INV_G2,
  INV_W2W3,
  INV_F8,
  INV_W3W3,
  INV_W2W5,
  INV_W2W7,
  INV_W3W5,
  INV_W3W7,
  INV_W2W3E2,
  INV_W2W5E2,
  INV_W2W13,
  INV_W2W7E2,
  INV_W3W3E2,
  INV_W5W7,
  INV_W3W13,
  INV_ATKIN3,
  INV_ATKIN5,
  INV_ATKIN7,
  INV_ATKIN11,
  INV_ATKIN13,
  INV_ATKIN17,
  INV_ATKIN19,
  INV_ATKIN23,
  INV_ATKIN29,
  INV_ATKIN31,
  INV_LAST,
  check_modinv,
  modinv_level,
  modinv_degree,
  modinv_odd_conductor,
  modinv_height_factor,
  modinv_sparse_factor,
  modinv_pfilter,
  modinv_good_prime,
  modinv_good_disc,
  modinv_ramified,
  modinv_is_Weber,
  modinv_is_double_eta,
  modinv_max_internal_level,
  disc_best_modinv,
  qfb_nform,
  // Modular polynomial database
  type PolmodularDB,
  polmodular_db_init,
  polmodular_db_for_inv,
  polmodular_db_add_level,
  polmodular_db_add_levels,
  polmodular_db_getp,
  // Modular polynomials
  type ZXX,
  sympol_to_ZM,
  Flm_Fl_polmodular_evalx,
  polmodular0_ZM,
  polmodular_ZM,
  polmodular_ZXX,
  Fp_polmodular_evalx,
  polmodular,
  // Hilbert class polynomials
  polclass0,
} from './polmodular.js';
// NOTE: `_internal` is deliberately NOT re-exported (test-only surface).

// Schoof-Elkies-Atkin point counting over F_p.
// Source: ellsea.c (+ FpE.c for the supersingularity test)
export {
  Fp_ellcard_SEA,
  Fp_elljissupersingular,
  Fq_elldivpolmod,
  setSeaDebugLevel,
} from './elliptic/ellsea.js';
// NOTE: `_internal` is deliberately NOT re-exported (test-only surface).
// NOTE: ellsea.ts re-uses PariBugError (polmodular.js) and PariPrimeError
// (matkermod.js) rather than defining its own, so there is no name clash.

// Representation numbers (theta series) of a positive definite integral
// quadratic form.  Source: bibli1.c (qfrep0/minim0_dolll), alglin2.c
// (qfgaussred_positive), lll.c (lllgramint)
export {
  qfrep0,
  qfrep,
  lllgramint,
  qfgaussred_positive,
  qf_ZM_apply,
  ZM_det,
  PariPrecError,
  type Frac,
  type BoundLike,
} from './qfrep.js';
// NOTE: qfrep.ts re-exports PariDomainError / PariTypeError / type ZM for the
// convenience of direct module importers; they must NOT be re-exported here
// (index.ts already has them from './matkermod.js').  `isqrt` is a module-local
// duplicate of `sqrti` and is likewise not re-exported.

// Class group and unit group of quadratic fields (index calculus).
// Source: buch1.c, plus hnf_snf.c / Qfb.c (qfr3-qfr5) / alglin1.c (ZM_pivots)
// and the GRH check of buch2.c.
//
// NOTE: buch.ts carries its own transcription of PARI's t_REAL kernel and of
// the qfr3/qfr5 containers, structurally identical to (but distinct from) the
// ones exported from './qfb.js' above.  Only the buch names that do not clash
// with that kernel are re-exported here; the clashing ones (real_0, real_1,
// itor, shiftr, addrr, subrr, mulrr, sqrr, divrr, mulir, divri, divru, divir,
// truncr, gcvtoi, sqrtr, mplog2, logr_abs, Qfr3, Qfr5, QfrData, qfr5_dist,
// qfr_data_init) are reachable from './qfb.js'.  Merging the two kernels into a
// single module is follow-up work.  `sqrti` and the shared PARI error classes
// that buch.ts re-exports are likewise omitted (already exported above).
export {
  // Multiprecision reals (PARI t_REAL) used by the Shanks distances
  type Real,
  type CReal,
  DEFAULTPREC,
  real_neg,
  real_abs,
  real_sign,
  real_expo,
  setprec,
  mulur,
  cmprr,
  rtodbl,
  dbltor,
  expr as mpexp,
  // Integer matrices: HNF/SNF of a relation matrix
  type ZC,
  type ZMat,
  type zv,
  ZM_mul,
  ZM_det_triangular,
  ZM_pivots,
  ZM_hnflll,
  hnfspec_i,
  hnfadd_i,
  ZM_snf_group,
  // Indefinite forms with Shanks distance (buch.ts's own qfr3/qfr5 layer)
  type Qfr3 as BuchQfr3,
  type Qfr5 as BuchQfr5,
  type QfrData as BuchQfrData,
  qfr3_rho,
  qfr5_rho,
  qfr3_red,
  qfr5_red,
  qfr3_comp,
  qfr5_comp,
  qfr3_pow,
  primeform_u3,
  // Class group / unit group
  type QuadClassUnit,
  type Bnf,
  bnf_increase_LIMC,
  setBuchRandomSeed,
  Buchquad,
  quadclassunit0,
  quadclassno,
  bnfinit,
} from './buch.js';

// Galois groups of number fields: galoisinit/galoisgen, galoispermtopol,
// galoisfixedfield, galoissubgroups.  Source: galconj.c (+ perm.c, Zp.c, FpX.c)
//
// NOTE: galconj.ts re-exports NotImplementedError, PariBugError,
// PariDomainError, PariFlagError, PariImplError, PariInvError and
// PariTypeError for convenience; they are NOT re-exported here because
// index.ts already has them from ifactor/matkermod/polmodular.
export {
  // errors
  PariIrredpolError,
  // integer helpers
  ugcd,
  ulcm,
  logint,
  factoru_small,
  Forprime,
  // ZX
  type ZX,
  ZX_renormalize,
  ZX_degree,
  ZX_add,
  ZX_sub,
  ZX_neg,
  ZX_mul,
  ZX_Z_mul,
  ZX_deriv,
  ZX_equal,
  ZX_is_monic,
  ZX_xn,
  ZX_divrem_monic,
  ZX_resultant,
  ZX_disc,
  ZX_is_squarefree,
  indexpartial,
  // FpX additions
  FpX_deriv,
  FpX_eval,
  FpX_center,
  FpX_div_by_X_x,
  FpX_extgcd,
  FpXQ_powers,
  FpXQ_powBig,
  FpX_FpXQ_eval,
  FpX_Frobenius,
  FpXQ_autpow,
  FpXQ_autpowers,
  FpV_roots_to_pol,
  FpXQ_minpoly,
  FpX_is_squarefree,
  FpX_split_part,
  FpX_nbroots,
  FpX_is_totally_split,
  FpX_ddf,
  FpX_factor_squarefree,
  cmp_FpX,
  FpX_roots,
  FpX_nbfact_by_degree,
  FpV_invVandermonde,
  // Zp
  ZpX_liftroot,
  ZpX_roots,
  ZpX_liftroots,
  ZpX_liftfact,
  bezout_lift_fact,
  ZpX_ZpXQ_liftroot,
  FpXQ_inv,
  // permutations and groups
  type Perm,
  type Group,
  type Quotient,
  identity_perm,
  perm_mul,
  perm_sqr,
  perm_inv,
  perm_conj,
  perm_commute,
  perm_powu,
  perm_cycles,
  perm_orderu,
  perm_relorder,
  vecperm_orbits,
  cyc_pow,
  vecpermute,
  vecsmall_lexcmp,
  vecsmall_uniq,
  zv_equal,
  trivialgroup,
  cyclicgroup,
  dicyclicgroup,
  group_order,
  group_domain,
  group_elts,
  group_set,
  groupelts_set,
  group_leftcoset,
  group_rightcoset,
  perm_generate,
  group_perm_normalize,
  groupelts_quotient,
  group_quotient,
  quotient_perm,
  quotient_subgroup_lift,
  quotient_group,
  group_isA4S4,
  group_subgroups,
  // galconj proper
  type QPoly,
  type SymPol,
  type GaloisBorne,
  type GaloisGens,
  type GaloisInit,
  type FixedField,
  QPoly_normalize,
  QPoly_to_fractions,
  QPoly_to_FpX,
  initgaloisborne,
  galoisborne,
  galoisanalysis,
  listznstarelts,
  sympol_eval,
  fixedfieldsympol,
  fixedfieldorbits,
  permtopol,
  galoisgen,
  galoisinit,
  galoisvecpermtopol,
  galoispermtopol,
  galoisfixedfield,
  galois_group,
  galoissubgroups,
  galoisconj4,
} from './galconj.js';
