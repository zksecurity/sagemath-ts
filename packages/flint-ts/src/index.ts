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
export { fmpz_get_d } from './fmpz/get.js';

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
export type { fmpz_factor } from './fmpz.js';
export {
  fmpz,
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
export type { fmpz_poly_factor } from './fmpz_poly.js';
export {
  fmpz_poly,
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
export type { nmod_t } from './nmod_poly.js';
export {
  nmod_poly,
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

// Dense-array kernels; lengths are inferred and inputs are not mutated.
export {
  _fmpz_poly_gcd,
  _fmpz_poly_gcd_subresultant,
  _fmpz_poly_gcd_heuristic,
  _fmpz_poly_gcd_modular,
} from './fmpz_poly/gcd.js';
export { _nmod_poly_gcd } from './nmod_poly/gcd.js';
export { _nmod_poly_make_monic } from './nmod_poly/make_monic.js';
export { _fmpz_poly_divrem } from './fmpz_poly/divrem.js';
export { _nmod_poly_divrem } from './nmod_poly/divrem.js';
export { _fmpq_poly_gcd } from './fmpq_poly/gcd.js';
export { _nmod_poly_xgcd } from './nmod_poly/xgcd.js';
export { _fmpz_poly_resultant } from './fmpz_poly/resultant.js';
export { _nmod_poly_resultant } from './nmod_poly/resultant.js';
export { _fmpz_poly_xgcd } from './fmpz_poly/xgcd.js';
export { _fmpq_poly_xgcd, type DenseRationalPolynomial } from './fmpq_poly/xgcd.js';

export { _fmpq_poly_resultant } from './fmpq_poly/resultant.js';
export { _fmpz_poly_derivative } from './fmpz_poly/derivative.js';
export { _fmpq_poly_derivative } from './fmpq_poly/derivative.js';
export {
  _nmod_poly_mul,
  _nmod_poly_mul_classical,
  _nmod_poly_mul_KS,
  _nmod_poly_mul_KS2,
  _nmod_poly_mul_KS4,
} from './nmod_poly/mul.js';
export { _nmod_poly_pow } from './nmod_poly/pow.js';
export { _nmod_poly_add } from './nmod_poly/add.js';
export { _nmod_poly_sub } from './nmod_poly/sub.js';
export { _fmpz_poly_pow } from './fmpz_poly/pow.js';
export { _fmpz_poly_pow_small } from './fmpz_poly/pow_small.js';
export { _fmpz_poly_pow_binomial } from './fmpz_poly/pow_binomial.js';
export { _fmpz_poly_pow_multinomial } from './fmpz_poly/pow_multinomial.js';
export { _fmpz_poly_pow_binexp } from './fmpz_poly/pow_binexp.js';
export { _fmpq_poly_pow } from './fmpq_poly/pow.js';
export { _fmpz_poly_mullow } from './fmpz_poly/mullow.js';
export { _fmpq_poly_mullow } from './fmpq_poly/mullow.js';
export { _nmod_poly_mullow } from './nmod_poly/mullow.js';
export { _fmpz_poly_pow_trunc } from './fmpz_poly/pow_trunc.js';
export { _nmod_poly_pow_trunc } from './nmod_poly/pow_trunc.js';
export { _fmpz_poly_inv_series } from './fmpz_poly/inv_series.js';
export { _fmpq_poly_inv_series_newton } from './fmpq_poly/inv_series_newton.js';
export { _fmpz_poly_mul } from './fmpz_poly/mul.js';
export { _fmpq_poly_mul } from './fmpq_poly/mul.js';
export { _nmod_poly_inv_series_newton } from './nmod_poly/inv_series_newton.js';
export { _nmod_poly_powmod_ui_binexp } from './nmod_poly/powmod_binexp.js';
export { _nmod_poly_powmod_fmpz_binexp_preinv } from './nmod_poly/powmod_binexp_preinv.js';
export { _nmod_poly_powmod_x_fmpz_preinv } from './nmod_poly/powmod_x_preinv.js';
export { _fmpz_poly_evaluate_fmpz } from './fmpz_poly/evaluate_fmpz.js';
export { _fmpz_poly_evaluate_fmpq } from './fmpz_poly/evaluate_fmpq.js';
export { _fmpq_poly_evaluate_fmpz, _fmpq_poly_evaluate_fmpq } from './fmpq_poly/evaluate.js';
export { _nmod_poly_evaluate_nmod } from './nmod_poly/evaluate_nmod.js';
export { _fmpz_poly_taylor_shift } from './fmpz_poly/taylor_shift.js';
export { _fmpz_poly_compose } from './fmpz_poly/compose.js';
export { _fmpq_poly_compose } from './fmpq_poly/compose.js';
export { _nmod_poly_compose } from './nmod_poly/compose.js';

export { fmpq_dedekind_sum } from './fmpq/dedekind_sum.js';

export { _n_jacobi_unsigned, n_jacobi_unsigned, n_jacobi } from './ulong_extras/jacobi.js';
export { n_is_square } from './ulong_extras/is_square.js';
export { n_preinvert_limb } from './ulong_extras/preinvert_limb.js';
export { n_powmod2_ui_preinv } from './ulong_extras/powmod2_ui_preinv.js';
export { n_sqrtmod } from './ulong_extras/sqrtmod.js';
export { _gr_poly_inv_series_basecase } from './gr_poly/inv_series_basecase.js';
export { _gr_poly_sqrt_series_basecase } from './gr_poly/sqrt_series_basecase.js';
export { _gr_poly_rsqrt_series_basecase } from './gr_poly/rsqrt_series_basecase.js';
export { _gr_poly_sqrt_series_newton } from './gr_poly/sqrt_series_newton.js';
export { _nmod_poly_sqrt } from './nmod_poly/sqrt.js';
export { _nmod_poly_sqrt_series } from './nmod_poly/sqrt_series.js';
export { _nmod_poly_mulhigh } from './nmod_poly/mulhigh.js';
export { _nmod_poly_mulhigh_classical } from './nmod_poly/mulhigh_classical.js';
export { nmod_poly_factor_squarefree } from './nmod_poly_factor/factor_squarefree.js';
export { fmpq_poly_get_numerator, fmpq_poly_get_denominator } from './fmpq_poly/get_numerator_denominator.js';
export { _fmpz_poly_lcm } from './fmpz_poly/lcm.js';
export { _fmpq_poly_lcm } from './fmpq_poly/lcm.js';
export { flint_rand_init, flint_rand_set_seed, flint_rand_get_seed, flint_rand_clear, type flint_rand_t } from './flint.js';
export { n_randlimb, n_randint, n_urandint, n_randbits, n_randtest_bits, n_randtest, n_randtest_not_zero } from './ulong_extras/randomisation.js';
export { _nmod_vec_rand, _nmod_vec_randtest } from './nmod_vec/rand.js';
export { nmod_poly_deflation } from './nmod_poly/deflation.js';
export { nmod_poly_deflate } from './nmod_poly/deflate.js';
export { nmod_poly_inflate } from './nmod_poly/inflate.js';
export { nmod_poly_remove } from './nmod_poly/remove.js';
export { nmod_poly_randtest } from './nmod_poly/randtest.js';
export { nmod_mat_mul } from './nmod_mat/mul.js';
export { _nmod_poly_mod_matrix_rows_evaluate } from './nmod_poly/mod_matrix_rows_evaluate.js';
export { nmod_poly_precompute_matrix, _nmod_poly_reduce_matrix_mod_poly, nmod_poly_compose_mod_brent_kung_precomp_preinv } from './nmod_poly/compose_mod_brent_kung_precomp_preinv.js';
export { nmod_poly_compose_mod_brent_kung_vec_preinv } from './nmod_poly/compose_mod_brent_kung_vec_preinv.js';
export { nmod_poly_factor_equal_deg_prob } from './nmod_poly_factor/factor_equal_deg_prob.js';
export { nmod_poly_factor_equal_deg } from './nmod_poly_factor/factor_equal_deg.js';
export { nmod_poly_factor_cantor_zassenhaus } from './nmod_poly_factor/factor_cantor_zassenhaus.js';
export { nmod_poly_factor_distinct_deg } from './nmod_poly_factor/factor_distinct_deg.js';
export { nmod_poly_factor_kaltofen_shoup } from './nmod_poly_factor/factor_kaltofen_shoup.js';
export { nmod_poly_factor, nmod_poly_factor_with_cantor_zassenhaus, nmod_poly_factor_with_kaltofen_shoup } from './nmod_poly_factor/factor.js';
