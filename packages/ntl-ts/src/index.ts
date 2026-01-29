/**
 * @module ntl-ts
 * @description TypeScript port of NTL (Number Theory Library) functions
 *
 * This package provides TypeScript implementations of NTL algorithms.
 * Reference: reference/ntl/src/
 *
 * NTL is a high-performance, portable C++ library providing data structures
 * and algorithms for arbitrary length integers, vectors, matrices, polynomials
 * over integers, finite fields, and arbitrary precision floating point numbers.
 *
 * @see https://libntl.org/
 */

export const VERSION = '0.0.1';

// ============================================
// Core Integer Types
// ============================================

/**
 * Arbitrary precision integers.
 * @see ZZ.ts
 */
export {
  ZZ,
  GCD,
  LCM,
  XGCD,
  PowerMod,
  InvMod,
  ProbPrime,
  RandomPrime,
  NextPrime,
} from './ZZ.js';

// ============================================
// Integers Modulo p
// ============================================

/**
 * Integers mod p (prime modulus).
 * @see ZZ_p.ts
 */
export {
  ZZ_p,
  ZZ_pContext,
  ZZ_pBak,
  add as ZZ_p_add,
  sub as ZZ_p_sub,
  mul as ZZ_p_mul,
  div as ZZ_p_div,
  inv as ZZ_p_inv,
  negate as ZZ_p_negate,
  sqr as ZZ_p_sqr,
  power as ZZ_p_power,
  IsZero as ZZ_p_IsZero,
  IsOne as ZZ_p_IsOne,
  conv as ZZ_p_conv,
} from './ZZ_p.js';

// ============================================
// Polynomials over Z_p
// ============================================

/**
 * Polynomials over Z_p.
 * @see ZZ_pX.ts
 */
export {
  ZZ_pX,
  deg as ZZ_pX_deg,
  coeff as ZZ_pX_coeff,
  LeadCoeff as ZZ_pX_LeadCoeff,
  ConstTerm as ZZ_pX_ConstTerm,
  GCD as ZZ_pX_GCD,
  XGCD as ZZ_pX_XGCD,
  IsZero as ZZ_pX_IsZero,
  IsOne as ZZ_pX_IsOne,
  eval_ as ZZ_pX_eval,
  diff as ZZ_pX_diff,
  MulMod as ZZ_pX_MulMod,
  PowerMod as ZZ_pX_PowerMod,
  InvMod as ZZ_pX_InvMod,
} from './ZZ_pX.js';

// ============================================
// Binary Field GF(2)
// ============================================

/**
 * Binary field GF(2).
 * @see GF2.ts
 */
export {
  GF2,
  to_GF2,
  rep as GF2_rep,
  add as GF2_add,
  sub as GF2_sub,
  mul as GF2_mul,
  div as GF2_div,
  power as GF2_power,
  IsZero as GF2_IsZero,
  IsOne as GF2_IsOne,
  inv as GF2_inv,
} from './GF2.js';

// ============================================
// Polynomials over GF(2)
// ============================================

/**
 * Polynomials over GF(2).
 * @see GF2X.ts
 */
export {
  GF2X,
  deg as GF2X_deg,
  coeff as GF2X_coeff,
  LeadCoeff as GF2X_LeadCoeff,
  ConstTerm as GF2X_ConstTerm,
  GCD as GF2X_GCD,
  XGCD as GF2X_XGCD,
  IsZero as GF2X_IsZero,
  IsOne as GF2X_IsOne,
  IsX as GF2X_IsX,
  MulMod as GF2X_MulMod,
  SqrMod as GF2X_SqrMod,
  InvMod as GF2X_InvMod,
  PowerMod as GF2X_PowerMod,
  BuildIrred as GF2X_BuildIrred,
} from './GF2X.js';

// ============================================
// Extension Field GF(2^n)
// ============================================

/**
 * Extension field GF(2^n).
 * @see GF2E.ts
 */
export {
  GF2E,
  GF2EContext,
  GF2EBak,
  rep as GF2E_rep,
  add as GF2E_add,
  sub as GF2E_sub,
  mul as GF2E_mul,
  div as GF2E_div,
  sqr as GF2E_sqr,
  inv as GF2E_inv,
  power as GF2E_power,
  IsZero as GF2E_IsZero,
  IsOne as GF2E_IsOne,
  conv as GF2E_conv,
  trace as GF2E_trace,
} from './GF2E.js';
