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
  Fp_order,
  znorder,
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
  ell_to_a4a6_bc,
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
  FpE_changepoint,
  FpE_changepointinv,
} from './elliptic/point.js';

// Integer factorization
// Source: ifactor1.c, arith2.c
export {
  // Factorization
  Z_factor,
  Z_isanypower,
  factoru,
  formatFactorization,
  // Primality testing
  isPrime,
  isprimepower,
  nextprime,
  precprime,
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
  type MpComplex,
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
// Buchmann real helpers delegate to the shared kernels. Its qfr3/qfr5 containers
// remain structurally identical but distinct. Only nonclashing Buchmann names are
// exported here; qfb supplies the common arithmetic, square-root and logarithm names.
// Buchmann's expr is the public mpexp entry and retains native BigInt exponents.
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
  numberofconjugates,
} from './galconj.js';

export { resultant } from './polarit2.js';

export { prime } from './prime.js';

export { sumdedekind } from './elltrans.js';

export { algdep } from './bibli1.js';

export { hilbert } from './arith1.js';

export { eulerphi, numdiv } from './arith2.js';

// Real factorial backend, basemath/trans2.c.
export { mpfactr } from './trans2.js';

export { PariError } from './errors.js';
export { pari_init_rand, pari_rand, setrand, getrand, random_bits, random_Fl, randomi, random_F2x, random_zv } from './random.js';
export { random_Flx } from './Flx.js';
export { random_FpX } from './FpX.js';
export { F2x_degree, F2x_add, F2x_mul, F2x_sqr, F2x_sqrt, F2x_rem, F2x_divrem, F2x_gcd, F2x_deriv, F2x_valrem, F2xq_mul, F2xq_sqr, F2xq_powers, F2x_Frobenius, F2x_matFrobenius } from './F2x.js';
export { F2m_ker_sp, F2m_ker } from './F2v.js';
export { F2x_factor_squarefree, F2x_ddf, F2x_factor } from './FpX_factor.js';

export { Flm_mul, FpM_mul } from './FpV.js';
export { F2m_mul } from './F2v.js';

export { brent_kung_optpow } from './RgX.js';
export { Flxq_powers, Flx_Flxq_eval, Flx_FlxqV_eval } from './Flx.js';
export { FpX_FpXQV_eval, FpXQ_auttrace } from './FpX.js';

export { ZX_sqr } from './ZX.js';
export { Flx_mul, Flx_sqr } from './Flx.js';
export { FpX_sqr } from './FpX.js';

export { FpX_invBarrett } from './FpX.js';
export { Flx_invBarrett, Flx_divrem, Flx_rem } from './Flx.js';

export { FpX_halfgcd, FpX_halfgcd_all } from './FpX.js';
export { Flx_gcd, Flx_extgcd, Flx_halfgcd, Flx_halfgcd_all, Flx_deriv, Flx_is_squarefree } from './Flx.js';
export type { PolynomialMatrix, HalfGcdResult } from './_polynomial_gcd.js';
export { gen_pow_i, gen_pow_fold, gen_order, type GroupOrder, type GroupOrderFactors } from './bb_group.js';

export { Flxq_minpoly } from './Flx.js';
export { Flx_ddf, Flx_nbfact_by_degree, FpX_nbfact, type PolynomialDegreeFactor } from './FpX_factor.js';

export { Fl_sqrt, Fp_sqrt_i } from './ff.js';

export { Flx_normalize } from './Flx.js';

export { FpX_factor, Flx_factor, type PolynomialFactor } from './FpX_factor.js';

export { Flx_roots, Flx_nbroots, Flx_is_totally_split } from './FpX_factor.js';

export { FpXQX_mul, FpXQX_sqr, FpXQX_red } from './FpXX.js';
export { FpXQX_normalize } from './polarit3.js';
export { FlxqX_mul, FlxqX_sqr, FlxqX_red, FlxqX_normalize } from './FlxX.js';
export { F2xqX_mul, F2xqX_sqr, F2xqX_red, F2xqX_normalize } from './F2x.js';
export type { ExtensionCoefficient, ExtensionPolynomial } from './_extension_polynomial.js';

export { FpXQX_divrem, FpXQX_rem, FpXQX_div, FpXQX_invBarrett, FpXQX_get_red } from './FpXX.js';
export { FlxqX_divrem, FlxqX_rem, FlxqX_div, FlxqX_invBarrett, FlxqX_get_red } from './FlxX.js';
export { F2xqX_divrem, F2xqX_rem, F2xqX_div, F2xqX_invBarrett, F2xqX_get_red } from './F2x.js';
export type {ExtensionReduction,ExtensionModulus} from './_extension_division.js';

export { FpXQX_gcd, FpXQX_extgcd, FpXQX_halfgcd } from './FpXX.js';
export { FlxqX_gcd, FlxqX_extgcd, FlxqX_halfgcd } from './FlxX.js';
export { F2xqX_gcd, F2xqX_extgcd, F2xqX_halfgcd } from './F2x.js';
export type { ExtensionMatrix } from './_extension_gcd.js';

export { FpXQXQ_mul, FpXQXQ_sqr, FpXQXQ_invsafe, FpXQXQ_inv, FpXQXQ_div, FpXQXQ_pow, FpXQXQ_powers } from './FpXX.js';
export { FlxqXQ_mul, FlxqXQ_sqr, FlxqXQ_invsafe, FlxqXQ_inv, FlxqXQ_div, FlxqXQ_pow, FlxqXQ_powers, FlxqXQ_powu } from './FlxX.js';
export { F2xqXQ_mul, F2xqXQ_sqr, F2xqXQ_invsafe, F2xqXQ_inv, F2xqXQ_pow, F2xqXQ_powers } from './F2x.js';

export { FlxqM_mul } from './alglin1.js';
export { FpXQX_FpXQXQ_eval, FpXQX_FpXQXQV_eval } from './FpXX.js';
export { FlxqX_FlxqXQ_eval, FlxqX_FlxqXQV_eval } from './FlxX.js';
export { F2xqX_F2xqXQ_eval, F2xqX_F2xqXQV_eval } from './F2x.js';

export { FpXY_FpXQ_evalx, FpXY_FpXQV_evalx } from './FpXX.js';
export { FlxY_Flxq_evalx, FlxY_FlxqV_evalx } from './FlxX.js';
export { F2xY_F2xq_evalx, F2xY_F2xqV_evalx, F2x_F2xq_eval, F2x_F2xqV_eval } from './F2x.js';

export { FpXQXQ_autpow, FpXQXQ_auttrace, FpXQXQ_autsum } from './FpXX.js';
export { FlxqXQ_autpow, FlxqXQ_auttrace, FlxqXQ_autsum } from './FlxX.js';
export { F2xqXQ_autpow, F2xqXQ_auttrace } from './F2x.js';

export { random_FpXQX, FpXQX_dotproduct, FpXQXn_mul, FpXQXn_sqr } from './FpXX.js';

export { random_FlxqX, FlxqXn_mul, FlxqXn_sqr } from './FlxX.js';

export { FlxqX_dotproduct } from './Flx.js';

export { FpXQXQ_minpoly } from './FpXX.js';
export { FlxqXQ_minpoly } from './FlxX.js';

export { FpXQX_Frobenius, FlxqX_Frobenius, F2xqX_Frobenius, FpXQXQ_halfFrobenius, FlxqXQ_halfFrobenius } from './FpXQX_factor.js';
export { Flx_Frobenius } from './Flx.js';

export { FpXQX_split_part, FpXQX_nbroots, FlxqX_nbroots, F2xqX_nbroots, FqX_nbroots, FpXQX_is_squarefree, FlxqX_is_squarefree } from './FpXQX_factor.js';
export { FpXX_deriv } from './FpXX.js';
export { FlxX_deriv } from './FlxX.js';

export { FF_issquareall } from './ff.js';

export { FF_issquare, FF_norm } from './ff.js';
export { Flx_resultant, Flxq_norm, Flxq_issquare } from './Flx.js';
export { FpX_resultant, FpXQ_norm, FpXQ_issquare } from './FpX.js';

export { FF_trace } from './ff.js';
export { FpXQ_trace } from './FpX.js';
export { Flxq_trace } from './Flx.js';
export { F2xq_trace } from './F2x.js';

export { FF_charpoly } from './ff.js';
export { FpXQ_charpoly, FpV_polint } from './FpX.js';
export { Flxq_charpoly, Flv_polint } from './Flx.js';
export { Flx_FlxY_resultant } from './polarit3.js';

export { F2xq_invsafe, F2xq_inv, F2xq_div } from './F2x.js';
export { type F2xqECoefficient, type F2xqEChange, F2xqE_changepoint,
  F2xqE_changepointinv, F2xqE_add, F2xqE_dbl, F2xqE_neg, F2xqE_sub,
  F2xqE_mul, F2xqE_order } from './F2xqE.js';

export {
  type FFEllipticScalar, type FFEllipticInvariants, type FFEllipticCurve,
  type BinaryFFEllipticCurve, type OddFFEllipticCurve,
  type FFEllipticPoint, type FFEllipticInputPoint,
  FF_ellinit, FF_ellmul, FF_ellorder, ellinit_Fq,
} from './_elliptic_finite_field.js';

export { type FqEllipticPoint, type FqEllipticChange, type FlxqECoefficient } from './_odd_elliptic.js';
export * from './FlxqE.js';
export * from './FpE.js';

export * from './FlE.js';

export { F2xq_autpow, F2xq_sqrt_fast, F2xq_sqrt } from './F2x.js';
export { Fl2_sqrt_pre, Flxq_sqrt_pre, Flxq_sqrt } from './Flx.js';

export { type ZpPolynomialTree, gen_ZpX_Dixon, gen_ZpX_Newton, ZpXQ_invlift, ZpXQ_inv, ZpXQ_div } from './Zp.js';
export { Flx_Teichmuller } from './Zp.js';

export { Zp_sqrtnlift, Zp_sqrtlift } from './Zp.js';
