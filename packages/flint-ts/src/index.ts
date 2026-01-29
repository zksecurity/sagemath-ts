/**
 * @module flint-ts
 * @description TypeScript port of FLINT library functions
 *
 * FLINT (Fast Library for Number Theory) is a C library for number theory
 * developed at the University of Warwick, UK. This package provides TypeScript
 * implementations of core FLINT algorithms.
 *
 * Reference: https://flintlib.org/
 * Source: reference/flint/src/
 *
 * @license GPL-2.0-or-later (following FLINT's LGPL license)
 */

export const VERSION = '0.0.1';

// ============================================================
// Core Types
// ============================================================

/**
 * Arbitrary precision integers (fmpz)
 *
 * The fmpz type is FLINT's main arbitrary precision integer type.
 * Uses a small/large representation for efficiency.
 *
 * @see Reference: flint/src/fmpz.h
 */
export {
  fmpz,
  fmpz_factor,
  // Functional API
  fmpz_init,
  fmpz_clear,
  fmpz_add,
  fmpz_sub,
  fmpz_mul,
  fmpz_tdiv_q,
  fmpz_gcd,
  fmpz_lcm,
  fmpz_is_prime,
  fmpz_factor_func,
} from './fmpz.js';

/**
 * Polynomials over the integers (fmpz_poly)
 *
 * Dense polynomials with fmpz coefficients.
 *
 * @see Reference: flint/src/fmpz_poly.h
 */
export {
  fmpz_poly,
  fmpz_poly_factor,
  // Functional API
  fmpz_poly_init,
  fmpz_poly_clear,
  fmpz_poly_set_coeff_fmpz,
  fmpz_poly_add,
  fmpz_poly_mul,
  fmpz_poly_gcd,
} from './fmpz_poly.js';

/**
 * Integers modulo n (fmpz_mod)
 *
 * Modular arithmetic in Z/nZ for arbitrary precision modulus.
 *
 * @see Reference: flint/src/fmpz_mod.h
 */
export {
  fmpz_mod_ctx,
  fmpz_mod_discrete_log_pohlig_hellman,
  // Functional API
  fmpz_mod_ctx_init,
  fmpz_mod_ctx_clear,
  fmpz_mod_ctx_modulus,
  fmpz_mod_add,
  fmpz_mod_sub,
  fmpz_mod_mul,
  fmpz_mod_neg,
  fmpz_mod_inv,
  fmpz_mod_pow_ui,
} from './fmpz_mod.js';

/**
 * Polynomials modulo n (nmod_poly)
 *
 * Polynomials with word-size modulus for fast arithmetic.
 *
 * @see Reference: flint/src/nmod_poly.h
 */
export {
  nmod_poly,
  nmod_t,
  nmod_poly_factor,
  // Functional API
  nmod_poly_init,
  nmod_poly_clear,
  nmod_poly_add,
  nmod_poly_sub,
  nmod_poly_mul,
  nmod_poly_gcd,
  nmod_poly_pow,
  nmod_poly_evaluate_nmod,
} from './nmod_poly.js';
