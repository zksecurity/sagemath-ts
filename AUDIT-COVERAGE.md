# Port audit coverage inventory

Generated with `bun tests/audit/coverage.ts`. This is an inventory, not a claim of full behavioral equivalence.

- Modular-integer/ring member dispatch coverage: 39/39; factory aliases/zero/negative orders are compared; broader coercion/category domains remain open.
- Parser/tokenizer implemented callable dispatch coverage: 23/23; excludes 7 throw-only stubs. Function-call branches remain incomplete.
- Production TypeScript files: 432.
- Named functions, constructors and public methods with bodies: 8365 (includes internal functions).
- Explicit throw-only stubs: 555.
- Other callables containing an explicit unimplemented branch: 253.
- Comparative areas: 29; configured seeded cases: 2163617.
- Integer method dispatch coverage: 107/107. Adapters and input limits are explicit in the oracle and DEVIATIONS.md.
- Rational method, accessor and constructor dispatch coverage: 84/84. Symbolic/complex return branches remain explicit implementation gaps.
- RationalField member dispatch coverage: 39/39; random streams are compared in rand_stats.
- IntegerRing member dispatch coverage: 10/10. Integer construction is also compared through the zz_* dispatchers with both wrapper modes.
- Finite-field parent/element member dispatch coverage: 105/105 across both prime-field implementations and the extension backend. Backend adapters and unimplemented input branches remain explicit in DEVIATIONS.md.
- Binary matrix supported API dispatch coverage: 51/51; two file-I/O stubs are excluded. Numeric/container and pixel-data adapters are explicit in DEVIATIONS.md.
- Instrumented source-line execution coverage: 111541/133728 (83.41%). LCOV does not establish branch or behavioral coverage.

Execution totals retain the previously recorded measurements for unchanged source files.
The edited ell_generic record is replaced by the 24.54.1 caller LCOV
measurements; their older records are not merged. These per-file measurements describe the same current source.

The September source review and repairs are recorded in [AUDIT-2026-09.md](AUDIT-2026-09.md).
Every file below remains in the audit inventory, including files absent from the supplied execution reports.
Existing implementation gaps are tracked separately in DEVIATIONS.md; stubs are not treated as implemented behavior.

## Comparative areas

| Area | Dispatch functions | Seeded cases |
|---|---:|---:|
| arith | 21 | 77 |
| arith_extended | 4 | 26 |
| arith_special | 80 | 428310 |
| coding_crypto | 80 | 920 |
| core_fidelity | 21 | 89 |
| ec_advanced | 70 | 10202 |
| elliptic_curves | 22 | 262 |
| finite_fields | 98 | 37802 |
| function_fields | 121 | 750140 |
| groups_modn | 80 | 15268 |
| hyperelliptic | 3 | 436 |
| integers | 129 | 8230 |
| lattices | 39 | 159 |
| lwe | 5 | 15 |
| matrix | 10 | 16 |
| matrix_extended | 135 | 236507 |
| matrix_ops | 59 | 137510 |
| modular_integers | 32 | 9199 |
| mpfr | 138 | 815 |
| number_fields | 51 | 57596 |
| padics_series | 83 | 490 |
| parser | 5 | 566 |
| polynomial_ops | 53 | 42209 |
| polynomials | 82 | 285639 |
| quadratic_forms | 26 | 82 |
| quaternion_algebras | 5 | 855 |
| rand_stats | 58 | 20427 |
| rationals | 136 | 7436 |
| real_literals | 17 | 112334 |

## Production files

| File | Callables | Throw-only stubs | Partial stubs | Executed lines |
|---|---:|---:|---:|---:|
| [packages/flint-ts/src/flint.ts](packages/flint-ts/src/flint.ts) | 4 | 0 | 0 | 8/8 (100.00%) |
| [packages/flint-ts/src/fmpq/dedekind_sum.ts](packages/flint-ts/src/fmpq/dedekind_sum.ts) | 1 | 0 | 0 | 26/28 (92.86%) |
| [packages/flint-ts/src/fmpq_poly/compose.ts](packages/flint-ts/src/fmpq_poly/compose.ts) | 1 | 0 | 0 | 27/28 (96.43%) |
| [packages/flint-ts/src/fmpq_poly/derivative.ts](packages/flint-ts/src/fmpq_poly/derivative.ts) | 1 | 0 | 0 | 7/7 (100.00%) |
| [packages/flint-ts/src/fmpq_poly/evaluate.ts](packages/flint-ts/src/fmpq_poly/evaluate.ts) | 2 | 0 | 0 | 18/18 (100.00%) |
| [packages/flint-ts/src/fmpq_poly/gcd.ts](packages/flint-ts/src/fmpq_poly/gcd.ts) | 1 | 0 | 0 | 17/18 (94.44%) |
| [packages/flint-ts/src/fmpq_poly/get_numerator_denominator.ts](packages/flint-ts/src/fmpq_poly/get_numerator_denominator.ts) | 2 | 0 | 0 | 4/4 (100.00%) |
| [packages/flint-ts/src/fmpq_poly/inv_series_newton.ts](packages/flint-ts/src/fmpq_poly/inv_series_newton.ts) | 1 | 0 | 0 | 54/57 (94.74%) |
| [packages/flint-ts/src/fmpq_poly/lcm.ts](packages/flint-ts/src/fmpq_poly/lcm.ts) | 1 | 0 | 0 | 8/8 (100.00%) |
| [packages/flint-ts/src/fmpq_poly/mul.ts](packages/flint-ts/src/fmpq_poly/mul.ts) | 1 | 0 | 0 | 18/18 (100.00%) |
| [packages/flint-ts/src/fmpq_poly/mullow.ts](packages/flint-ts/src/fmpq_poly/mullow.ts) | 1 | 0 | 0 | 14/14 (100.00%) |
| [packages/flint-ts/src/fmpq_poly/pow.ts](packages/flint-ts/src/fmpq_poly/pow.ts) | 1 | 0 | 0 | 7/7 (100.00%) |
| [packages/flint-ts/src/fmpq_poly/rem.ts](packages/flint-ts/src/fmpq_poly/rem.ts) | 2 | 0 | 0 | 22/22 (100.00%) |
| [packages/flint-ts/src/fmpq_poly/resultant.ts](packages/flint-ts/src/fmpq_poly/resultant.ts) | 2 | 0 | 0 | 55/57 (96.49%) |
| [packages/flint-ts/src/fmpq_poly/xgcd.ts](packages/flint-ts/src/fmpq_poly/xgcd.ts) | 3 | 0 | 0 | 61/62 (98.39%) |
| [packages/flint-ts/src/fmpz.ts](packages/flint-ts/src/fmpz.ts) | 88 | 78 | 0 | 89/185 (48.11%) |
| [packages/flint-ts/src/fmpz/get.ts](packages/flint-ts/src/fmpz/get.ts) | 1 | 0 | 0 | 8/8 (100.00%) |
| [packages/flint-ts/src/fmpz_mod.ts](packages/flint-ts/src/fmpz_mod.ts) | 54 | 45 | 0 | 58/112 (51.79%) |
| [packages/flint-ts/src/fmpz_poly.ts](packages/flint-ts/src/fmpz_poly.ts) | 91 | 85 | 0 | 93/184 (50.54%) |
| [packages/flint-ts/src/fmpz_poly/compose.ts](packages/flint-ts/src/fmpz_poly/compose.ts) | 1 | 0 | 0 | 47/47 (100.00%) |
| [packages/flint-ts/src/fmpz_poly/derivative.ts](packages/flint-ts/src/fmpz_poly/derivative.ts) | 1 | 0 | 0 | 4/4 (100.00%) |
| [packages/flint-ts/src/fmpz_poly/divrem.ts](packages/flint-ts/src/fmpz_poly/divrem.ts) | 9 | 0 | 0 | 112/116 (96.55%) |
| [packages/flint-ts/src/fmpz_poly/evaluate_fmpq.ts](packages/flint-ts/src/fmpz_poly/evaluate_fmpq.ts) | 1 | 0 | 0 | 57/59 (96.61%) |
| [packages/flint-ts/src/fmpz_poly/evaluate_fmpz.ts](packages/flint-ts/src/fmpz_poly/evaluate_fmpz.ts) | 1 | 0 | 0 | 39/40 (97.50%) |
| [packages/flint-ts/src/fmpz_poly/gcd.ts](packages/flint-ts/src/fmpz_poly/gcd.ts) | 20 | 0 | 0 | 344/357 (96.36%) |
| [packages/flint-ts/src/fmpz_poly/inv_series.ts](packages/flint-ts/src/fmpz_poly/inv_series.ts) | 1 | 0 | 0 | 32/33 (96.97%) |
| [packages/flint-ts/src/fmpz_poly/lcm.ts](packages/flint-ts/src/fmpz_poly/lcm.ts) | 1 | 0 | 0 | 12/12 (100.00%) |
| [packages/flint-ts/src/fmpz_poly/mul.ts](packages/flint-ts/src/fmpz_poly/mul.ts) | 1 | 0 | 0 | 28/28 (100.00%) |
| [packages/flint-ts/src/fmpz_poly/mullow.ts](packages/flint-ts/src/fmpz_poly/mullow.ts) | 1 | 0 | 0 | 22/23 (95.65%) |
| [packages/flint-ts/src/fmpz_poly/pow.ts](packages/flint-ts/src/fmpz_poly/pow.ts) | 1 | 0 | 0 | 16/16 (100.00%) |
| [packages/flint-ts/src/fmpz_poly/pow_binexp.ts](packages/flint-ts/src/fmpz_poly/pow_binexp.ts) | 1 | 0 | 0 | 15/15 (100.00%) |
| [packages/flint-ts/src/fmpz_poly/pow_binomial.ts](packages/flint-ts/src/fmpz_poly/pow_binomial.ts) | 1 | 0 | 0 | 35/35 (100.00%) |
| [packages/flint-ts/src/fmpz_poly/pow_multinomial.ts](packages/flint-ts/src/fmpz_poly/pow_multinomial.ts) | 1 | 0 | 0 | 26/26 (100.00%) |
| [packages/flint-ts/src/fmpz_poly/pow_small.ts](packages/flint-ts/src/fmpz_poly/pow_small.ts) | 1 | 0 | 0 | 8/8 (100.00%) |
| [packages/flint-ts/src/fmpz_poly/pow_trunc.ts](packages/flint-ts/src/fmpz_poly/pow_trunc.ts) | 1 | 0 | 0 | 16/16 (100.00%) |
| [packages/flint-ts/src/fmpz_poly/pseudo_divrem_basecase.ts](packages/flint-ts/src/fmpz_poly/pseudo_divrem_basecase.ts) | 1 | 0 | 0 | 22/23 (95.65%) |
| [packages/flint-ts/src/fmpz_poly/pseudo_divrem_divconquer.ts](packages/flint-ts/src/fmpz_poly/pseudo_divrem_divconquer.ts) | 2 | 0 | 0 | 60/61 (98.36%) |
| [packages/flint-ts/src/fmpz_poly/pseudo_rem.ts](packages/flint-ts/src/fmpz_poly/pseudo_rem.ts) | 1 | 0 | 0 | 7/7 (100.00%) |
| [packages/flint-ts/src/fmpz_poly/resultant.ts](packages/flint-ts/src/fmpz_poly/resultant.ts) | 3 | 0 | 0 | 79/81 (97.53%) |
| [packages/flint-ts/src/fmpz_poly/taylor_shift.ts](packages/flint-ts/src/fmpz_poly/taylor_shift.ts) | 1 | 0 | 0 | 50/52 (96.15%) |
| [packages/flint-ts/src/fmpz_poly/xgcd.ts](packages/flint-ts/src/fmpz_poly/xgcd.ts) | 1 | 0 | 0 | 87/92 (94.57%) |
| [packages/flint-ts/src/gr_poly/inv_series_basecase.ts](packages/flint-ts/src/gr_poly/inv_series_basecase.ts) | 1 | 0 | 0 | 18/19 (94.74%) |
| [packages/flint-ts/src/gr_poly/rsqrt_series_basecase.ts](packages/flint-ts/src/gr_poly/rsqrt_series_basecase.ts) | 1 | 0 | 0 | 14/14 (100.00%) |
| [packages/flint-ts/src/gr_poly/sqrt_series_basecase.ts](packages/flint-ts/src/gr_poly/sqrt_series_basecase.ts) | 1 | 0 | 0 | 20/20 (100.00%) |
| [packages/flint-ts/src/gr_poly/sqrt_series_newton.ts](packages/flint-ts/src/gr_poly/sqrt_series_newton.ts) | 1 | 0 | 0 | 30/31 (96.77%) |
| [packages/flint-ts/src/index.ts](packages/flint-ts/src/index.ts) | 0 | 0 | 0 | 87/87 (100.00%) |
| [packages/flint-ts/src/nmod_mat/mul.ts](packages/flint-ts/src/nmod_mat/mul.ts) | 1 | 0 | 0 | 90/92 (97.83%) |
| [packages/flint-ts/src/nmod_poly.ts](packages/flint-ts/src/nmod_poly.ts) | 93 | 85 | 0 | 95/188 (50.53%) |
| [packages/flint-ts/src/nmod_poly/add.ts](packages/flint-ts/src/nmod_poly/add.ts) | 1 | 0 | 0 | 8/8 (100.00%) |
| [packages/flint-ts/src/nmod_poly/compose.ts](packages/flint-ts/src/nmod_poly/compose.ts) | 1 | 0 | 0 | 101/104 (97.12%) |
| [packages/flint-ts/src/nmod_poly/compose_mod_brent_kung_precomp_preinv.ts](packages/flint-ts/src/nmod_poly/compose_mod_brent_kung_precomp_preinv.ts) | 3 | 0 | 0 | 71/71 (100.00%) |
| [packages/flint-ts/src/nmod_poly/compose_mod_brent_kung_vec_preinv.ts](packages/flint-ts/src/nmod_poly/compose_mod_brent_kung_vec_preinv.ts) | 1 | 0 | 0 | 46/48 (95.83%) |
| [packages/flint-ts/src/nmod_poly/deflate.ts](packages/flint-ts/src/nmod_poly/deflate.ts) | 1 | 0 | 0 | 10/10 (100.00%) |
| [packages/flint-ts/src/nmod_poly/deflation.ts](packages/flint-ts/src/nmod_poly/deflation.ts) | 1 | 0 | 0 | 19/20 (95.00%) |
| [packages/flint-ts/src/nmod_poly/divrem.ts](packages/flint-ts/src/nmod_poly/divrem.ts) | 1 | 0 | 0 | 12/13 (92.31%) |
| [packages/flint-ts/src/nmod_poly/evaluate_nmod.ts](packages/flint-ts/src/nmod_poly/evaluate_nmod.ts) | 1 | 0 | 0 | 10/10 (100.00%) |
| [packages/flint-ts/src/nmod_poly/gcd.ts](packages/flint-ts/src/nmod_poly/gcd.ts) | 15 | 0 | 0 | 239/242 (98.76%) |
| [packages/flint-ts/src/nmod_poly/inflate.ts](packages/flint-ts/src/nmod_poly/inflate.ts) | 1 | 0 | 0 | 11/11 (100.00%) |
| [packages/flint-ts/src/nmod_poly/inv_series_newton.ts](packages/flint-ts/src/nmod_poly/inv_series_newton.ts) | 1 | 0 | 0 | 28/28 (100.00%) |
| [packages/flint-ts/src/nmod_poly/make_monic.ts](packages/flint-ts/src/nmod_poly/make_monic.ts) | 1 | 0 | 0 | 17/17 (100.00%) |
| [packages/flint-ts/src/nmod_poly/mod_matrix_rows_evaluate.ts](packages/flint-ts/src/nmod_poly/mod_matrix_rows_evaluate.ts) | 1 | 0 | 0 | 36/37 (97.30%) |
| [packages/flint-ts/src/nmod_poly/mul.ts](packages/flint-ts/src/nmod_poly/mul.ts) | 13 | 0 | 0 | 163/167 (97.60%) |
| [packages/flint-ts/src/nmod_poly/mulhigh.ts](packages/flint-ts/src/nmod_poly/mulhigh.ts) | 1 | 0 | 0 | 13/13 (100.00%) |
| [packages/flint-ts/src/nmod_poly/mulhigh_classical.ts](packages/flint-ts/src/nmod_poly/mulhigh_classical.ts) | 1 | 0 | 0 | 10/10 (100.00%) |
| [packages/flint-ts/src/nmod_poly/mullow.ts](packages/flint-ts/src/nmod_poly/mullow.ts) | 1 | 0 | 0 | 31/31 (100.00%) |
| [packages/flint-ts/src/nmod_poly/pow.ts](packages/flint-ts/src/nmod_poly/pow.ts) | 1 | 0 | 0 | 25/27 (92.59%) |
| [packages/flint-ts/src/nmod_poly/pow_trunc.ts](packages/flint-ts/src/nmod_poly/pow_trunc.ts) | 1 | 0 | 0 | 23/23 (100.00%) |
| [packages/flint-ts/src/nmod_poly/powmod_binexp.ts](packages/flint-ts/src/nmod_poly/powmod_binexp.ts) | 1 | 0 | 0 | 25/25 (100.00%) |
| [packages/flint-ts/src/nmod_poly/powmod_binexp_preinv.ts](packages/flint-ts/src/nmod_poly/powmod_binexp_preinv.ts) | 2 | 0 | 0 | 45/45 (100.00%) |
| [packages/flint-ts/src/nmod_poly/powmod_x_preinv.ts](packages/flint-ts/src/nmod_poly/powmod_x_preinv.ts) | 1 | 0 | 0 | 46/46 (100.00%) |
| [packages/flint-ts/src/nmod_poly/randtest.ts](packages/flint-ts/src/nmod_poly/randtest.ts) | 1 | 0 | 0 | 5/5 (100.00%) |
| [packages/flint-ts/src/nmod_poly/remove.ts](packages/flint-ts/src/nmod_poly/remove.ts) | 1 | 0 | 0 | 18/18 (100.00%) |
| [packages/flint-ts/src/nmod_poly/resultant.ts](packages/flint-ts/src/nmod_poly/resultant.ts) | 2 | 0 | 0 | 59/62 (95.16%) |
| [packages/flint-ts/src/nmod_poly/sqrt.ts](packages/flint-ts/src/nmod_poly/sqrt.ts) | 1 | 0 | 0 | 26/26 (100.00%) |
| [packages/flint-ts/src/nmod_poly/sqrt_series.ts](packages/flint-ts/src/nmod_poly/sqrt_series.ts) | 1 | 0 | 0 | 11/11 (100.00%) |
| [packages/flint-ts/src/nmod_poly/sub.ts](packages/flint-ts/src/nmod_poly/sub.ts) | 1 | 0 | 0 | 8/8 (100.00%) |
| [packages/flint-ts/src/nmod_poly/xgcd.ts](packages/flint-ts/src/nmod_poly/xgcd.ts) | 2 | 0 | 0 | 71/73 (97.26%) |
| [packages/flint-ts/src/nmod_poly_factor/factor.ts](packages/flint-ts/src/nmod_poly_factor/factor.ts) | 5 | 0 | 0 | 43/43 (100.00%) |
| [packages/flint-ts/src/nmod_poly_factor/factor_cantor_zassenhaus.ts](packages/flint-ts/src/nmod_poly_factor/factor_cantor_zassenhaus.ts) | 1 | 0 | 0 | 29/29 (100.00%) |
| [packages/flint-ts/src/nmod_poly_factor/factor_distinct_deg.ts](packages/flint-ts/src/nmod_poly_factor/factor_distinct_deg.ts) | 1 | 0 | 0 | 99/99 (100.00%) |
| [packages/flint-ts/src/nmod_poly_factor/factor_equal_deg.ts](packages/flint-ts/src/nmod_poly_factor/factor_equal_deg.ts) | 1 | 0 | 0 | 23/23 (100.00%) |
| [packages/flint-ts/src/nmod_poly_factor/factor_equal_deg_prob.ts](packages/flint-ts/src/nmod_poly_factor/factor_equal_deg_prob.ts) | 1 | 0 | 0 | 39/40 (97.50%) |
| [packages/flint-ts/src/nmod_poly_factor/factor_kaltofen_shoup.ts](packages/flint-ts/src/nmod_poly_factor/factor_kaltofen_shoup.ts) | 1 | 0 | 0 | 24/24 (100.00%) |
| [packages/flint-ts/src/nmod_poly_factor/factor_squarefree.ts](packages/flint-ts/src/nmod_poly_factor/factor_squarefree.ts) | 1 | 0 | 0 | 24/24 (100.00%) |
| [packages/flint-ts/src/nmod_vec/rand.ts](packages/flint-ts/src/nmod_vec/rand.ts) | 2 | 0 | 0 | 11/11 (100.00%) |
| [packages/flint-ts/src/ulong_extras/is_square.ts](packages/flint-ts/src/ulong_extras/is_square.ts) | 1 | 0 | 0 | 13/14 (92.86%) |
| [packages/flint-ts/src/ulong_extras/jacobi.ts](packages/flint-ts/src/ulong_extras/jacobi.ts) | 3 | 0 | 0 | 21/23 (91.30%) |
| [packages/flint-ts/src/ulong_extras/powmod2_ui_preinv.ts](packages/flint-ts/src/ulong_extras/powmod2_ui_preinv.ts) | 1 | 0 | 0 | 13/14 (92.86%) |
| [packages/flint-ts/src/ulong_extras/preinvert_limb.ts](packages/flint-ts/src/ulong_extras/preinvert_limb.ts) | 1 | 0 | 0 | 4/4 (100.00%) |
| [packages/flint-ts/src/ulong_extras/randomisation.ts](packages/flint-ts/src/ulong_extras/randomisation.ts) | 7 | 0 | 0 | 55/57 (96.49%) |
| [packages/flint-ts/src/ulong_extras/sqrtmod.ts](packages/flint-ts/src/ulong_extras/sqrtmod.ts) | 1 | 0 | 0 | 45/45 (100.00%) |
| [packages/gsl-ts/src/index.ts](packages/gsl-ts/src/index.ts) | 0 | 0 | 0 | 3/3 (100.00%) |
| [packages/gsl-ts/src/specfunc/exp.ts](packages/gsl-ts/src/specfunc/exp.ts) | 1 | 0 | 0 | 4/4 (100.00%) |
| [packages/gsl-ts/src/specfunc/log.ts](packages/gsl-ts/src/specfunc/log.ts) | 1 | 0 | 0 | 2/2 (100.00%) |
| [packages/gsl-ts/src/sys/pow_int.ts](packages/gsl-ts/src/sys/pow_int.ts) | 1 | 0 | 0 | 14/14 (100.00%) |
| [packages/m4ri-ts/src/brilliantrussian.ts](packages/m4ri-ts/src/brilliantrussian.ts) | 11 | 0 | 0 | 244/244 (100.00%) |
| [packages/m4ri-ts/src/echelonform.ts](packages/m4ri-ts/src/echelonform.ts) | 3 | 0 | 0 | 27/27 (100.00%) |
| [packages/m4ri-ts/src/index.ts](packages/m4ri-ts/src/index.ts) | 0 | 0 | 0 | 10/10 (100.00%) |
| [packages/m4ri-ts/src/mzd.ts](packages/m4ri-ts/src/mzd.ts) | 10 | 0 | 0 | 160/167 (95.81%) |
| [packages/m4ri-ts/src/ple.ts](packages/m4ri-ts/src/ple.ts) | 5 | 0 | 0 | 100/104 (96.15%) |
| [packages/m4ri-ts/src/ple_russian.ts](packages/m4ri-ts/src/ple_russian.ts) | 4 | 0 | 0 | 182/187 (97.33%) |
| [packages/m4ri-ts/src/solve.ts](packages/m4ri-ts/src/solve.ts) | 2 | 0 | 0 | 66/74 (89.19%) |
| [packages/m4ri-ts/src/strassen.ts](packages/m4ri-ts/src/strassen.ts) | 1 | 0 | 0 | 91/96 (94.79%) |
| [packages/m4ri-ts/src/triangular.ts](packages/m4ri-ts/src/triangular.ts) | 3 | 0 | 0 | 35/36 (97.22%) |
| [packages/m4ri-ts/src/triangular_russian.ts](packages/m4ri-ts/src/triangular_russian.ts) | 1 | 0 | 0 | 49/51 (96.08%) |
| [packages/mpfr-ts/src/add.ts](packages/mpfr-ts/src/add.ts) | 1 | 0 | 0 | 42/48 (87.50%) |
| [packages/mpfr-ts/src/arithmetic.ts](packages/mpfr-ts/src/arithmetic.ts) | 2 | 0 | 0 | 24/27 (88.89%) |
| [packages/mpfr-ts/src/cmp.ts](packages/mpfr-ts/src/cmp.ts) | 1 | 0 | 0 | 12/12 (100.00%) |
| [packages/mpfr-ts/src/cmp_si.ts](packages/mpfr-ts/src/cmp_si.ts) | 1 | 0 | 0 | 13/13 (100.00%) |
| [packages/mpfr-ts/src/frac.ts](packages/mpfr-ts/src/frac.ts) | 1 | 0 | 0 | 44/46 (95.65%) |
| [packages/mpfr-ts/src/get_d.ts](packages/mpfr-ts/src/get_d.ts) | 1 | 0 | 0 | 12/12 (100.00%) |
| [packages/mpfr-ts/src/get_str.ts](packages/mpfr-ts/src/get_str.ts) | 1 | 0 | 0 | 31/31 (100.00%) |
| [packages/mpfr-ts/src/get_z.ts](packages/mpfr-ts/src/get_z.ts) | 1 | 0 | 0 | 10/12 (83.33%) |
| [packages/mpfr-ts/src/index.ts](packages/mpfr-ts/src/index.ts) | 0 | 0 | 0 | 19/19 (100.00%) |
| [packages/mpfr-ts/src/init2.ts](packages/mpfr-ts/src/init2.ts) | 1 | 0 | 0 | 6/6 (100.00%) |
| [packages/mpfr-ts/src/isinf.ts](packages/mpfr-ts/src/isinf.ts) | 1 | 0 | 0 | 2/2 (100.00%) |
| [packages/mpfr-ts/src/isinteger.ts](packages/mpfr-ts/src/isinteger.ts) | 1 | 0 | 0 | 5/5 (100.00%) |
| [packages/mpfr-ts/src/isnan.ts](packages/mpfr-ts/src/isnan.ts) | 1 | 0 | 0 | 2/2 (100.00%) |
| [packages/mpfr-ts/src/isnum.ts](packages/mpfr-ts/src/isnum.ts) | 1 | 0 | 0 | 2/2 (100.00%) |
| [packages/mpfr-ts/src/mul.ts](packages/mpfr-ts/src/mul.ts) | 1 | 0 | 0 | 16/16 (100.00%) |
| [packages/mpfr-ts/src/rint.ts](packages/mpfr-ts/src/rint.ts) | 6 | 0 | 0 | 70/78 (89.74%) |
| [packages/mpfr-ts/src/round_raw.ts](packages/mpfr-ts/src/round_raw.ts) | 3 | 1 | 0 | 10/10 (100.00%) |
| [packages/mpfr-ts/src/set.ts](packages/mpfr-ts/src/set.ts) | 1 | 0 | 0 | 28/28 (100.00%) |
| [packages/mpfr-ts/src/set_d.ts](packages/mpfr-ts/src/set_d.ts) | 1 | 0 | 0 | 45/45 (100.00%) |
| [packages/mpfr-ts/src/set_str.ts](packages/mpfr-ts/src/set_str.ts) | 1 | 0 | 0 | 45/46 (97.83%) |
| [packages/mpfr-ts/src/set_z.ts](packages/mpfr-ts/src/set_z.ts) | 1 | 0 | 0 | 30/30 (100.00%) |
| [packages/mpfr-ts/src/sgn.ts](packages/mpfr-ts/src/sgn.ts) | 1 | 0 | 0 | 2/2 (100.00%) |
| [packages/mpfr-ts/src/strtofr.ts](packages/mpfr-ts/src/strtofr.ts) | 1 | 0 | 0 | 16/20 (80.00%) |
| [packages/mpfr-ts/src/types.ts](packages/mpfr-ts/src/types.ts) | 0 | 0 | 0 | unmeasured |
| [packages/ntl-ts/src/FFT.ts](packages/ntl-ts/src/FFT.ts) | 20 | 0 | 0 | 279/288 (96.88%) |
| [packages/ntl-ts/src/FFT_impl.ts](packages/ntl-ts/src/FFT_impl.ts) | 1 | 0 | 0 | 10/10 (100.00%) |
| [packages/ntl-ts/src/GF2.ts](packages/ntl-ts/src/GF2.ts) | 30 | 1 | 0 | 53/79 (67.09%) |
| [packages/ntl-ts/src/GF2E.ts](packages/ntl-ts/src/GF2E.ts) | 42 | 30 | 0 | 52/94 (55.32%) |
| [packages/ntl-ts/src/GF2X.ts](packages/ntl-ts/src/GF2X.ts) | 72 | 6 | 0 | 375/456 (82.24%) |
| [packages/ntl-ts/src/GF2X_irred_tab.ts](packages/ntl-ts/src/GF2X_irred_tab.ts) | 0 | 0 | 0 | 345/345 (100.00%) |
| [packages/ntl-ts/src/LLL.ts](packages/ntl-ts/src/LLL.ts) | 5 | 0 | 0 | 195/200 (97.50%) |
| [packages/ntl-ts/src/ZZ.ts](packages/ntl-ts/src/ZZ.ts) | 95 | 44 | 0 | 704/765 (92.03%) |
| [packages/ntl-ts/src/ZZX1.ts](packages/ntl-ts/src/ZZX1.ts) | 9 | 0 | 0 | 278/286 (97.20%) |
| [packages/ntl-ts/src/ZZXFactoring.ts](packages/ntl-ts/src/ZZXFactoring.ts) | 47 | 0 | 0 | 1418/1418 (100.00%) |
| [packages/ntl-ts/src/ZZ_p.ts](packages/ntl-ts/src/ZZ_p.ts) | 37 | 26 | 0 | 47/84 (55.95%) |
| [packages/ntl-ts/src/ZZ_pEX.ts](packages/ntl-ts/src/ZZ_pEX.ts) | 9 | 0 | 0 | 294/369 (79.67%) |
| [packages/ntl-ts/src/ZZ_pX.ts](packages/ntl-ts/src/ZZ_pX.ts) | 86 | 46 | 0 | 572/630 (90.79%) |
| [packages/ntl-ts/src/ZZ_pX1.ts](packages/ntl-ts/src/ZZ_pX1.ts) | 2 | 0 | 0 | 134/141 (95.04%) |
| [packages/ntl-ts/src/index.ts](packages/ntl-ts/src/index.ts) | 0 | 0 | 0 | 81/81 (100.00%) |
| [packages/ntl-ts/src/lip.ts](packages/ntl-ts/src/lip.ts) | 1 | 0 | 0 | 24/24 (100.00%) |
| [packages/ntl-ts/src/lzz_pX.ts](packages/ntl-ts/src/lzz_pX.ts) | 31 | 0 | 0 | 531/555 (95.68%) |
| [packages/ntl-ts/src/lzz_pX1.ts](packages/ntl-ts/src/lzz_pX1.ts) | 21 | 0 | 0 | 455/473 (96.19%) |
| [packages/ntl-ts/src/lzz_pXFactoring.ts](packages/ntl-ts/src/lzz_pXFactoring.ts) | 14 | 0 | 0 | 392/397 (98.74%) |
| [packages/ntl-ts/src/mat_ZZ.ts](packages/ntl-ts/src/mat_ZZ.ts) | 7 | 0 | 0 | 183/188 (97.34%) |
| [packages/ntl-ts/src/mat_ZZ_p.ts](packages/ntl-ts/src/mat_ZZ_p.ts) | 1 | 0 | 0 | 45/45 (100.00%) |
| [packages/ntl-ts/src/mat_lzz_p.ts](packages/ntl-ts/src/mat_lzz_p.ts) | 15 | 0 | 0 | 346/352 (98.30%) |
| [packages/parigp-ts/src/F2v.ts](packages/parigp-ts/src/F2v.ts) | 18 | 0 | 0 | 286/287 (99.65%) |
| [packages/parigp-ts/src/F2x.ts](packages/parigp-ts/src/F2x.ts) | 42 | 0 | 0 | 218/226 (96.46%) |
| [packages/parigp-ts/src/F3v.ts](packages/parigp-ts/src/F3v.ts) | 1 | 0 | 0 | 36/37 (97.30%) |
| [packages/parigp-ts/src/Flv.ts](packages/parigp-ts/src/Flv.ts) | 4 | 0 | 0 | 29/29 (100.00%) |
| [packages/parigp-ts/src/Flx.ts](packages/parigp-ts/src/Flx.ts) | 28 | 0 | 0 | 210/210 (100.00%) |
| [packages/parigp-ts/src/FlxX.ts](packages/parigp-ts/src/FlxX.ts) | 32 | 0 | 0 | 97/105 (92.38%) |
| [packages/parigp-ts/src/FpV.ts](packages/parigp-ts/src/FpV.ts) | 2 | 0 | 0 | 36/36 (100.00%) |
| [packages/parigp-ts/src/FpX.ts](packages/parigp-ts/src/FpX.ts) | 32 | 0 | 0 | 240/263 (91.25%) |
| [packages/parigp-ts/src/FpXQX_factor.ts](packages/parigp-ts/src/FpXQX_factor.ts) | 12 | 0 | 0 | 40/40 (100.00%) |
| [packages/parigp-ts/src/FpXX.ts](packages/parigp-ts/src/FpXX.ts) | 31 | 0 | 0 | 124/128 (96.88%) |
| [packages/parigp-ts/src/FpX_factor.ts](packages/parigp-ts/src/FpX_factor.ts) | 29 | 0 | 0 | 276/292 (94.52%) |
| [packages/parigp-ts/src/QX_factor.ts](packages/parigp-ts/src/QX_factor.ts) | 26 | 0 | 0 | 611/623 (98.07%) |
| [packages/parigp-ts/src/RgV.ts](packages/parigp-ts/src/RgV.ts) | 4 | 0 | 0 | 60/61 (98.36%) |
| [packages/parigp-ts/src/RgX.ts](packages/parigp-ts/src/RgX.ts) | 5 | 0 | 0 | 71/71 (100.00%) |
| [packages/parigp-ts/src/ZV.ts](packages/parigp-ts/src/ZV.ts) | 7 | 0 | 0 | 119/125 (95.20%) |
| [packages/parigp-ts/src/ZX.ts](packages/parigp-ts/src/ZX.ts) | 7 | 0 | 0 | 121/125 (96.80%) |
| [packages/parigp-ts/src/Zp.ts](packages/parigp-ts/src/Zp.ts) | 14 | 0 | 0 | 232/236 (98.31%) |
| [packages/parigp-ts/src/_binary64.ts](packages/parigp-ts/src/_binary64.ts) | 4 | 0 | 0 | 59/62 (95.16%) |
| [packages/parigp-ts/src/_binary_composition.ts](packages/parigp-ts/src/_binary_composition.ts) | 1 | 0 | 0 | 28/29 (96.55%) |
| [packages/parigp-ts/src/_coefficient_composition.ts](packages/parigp-ts/src/_coefficient_composition.ts) | 1 | 0 | 0 | 54/54 (100.00%) |
| [packages/parigp-ts/src/_ecm.ts](packages/parigp-ts/src/_ecm.ts) | 17 | 0 | 0 | 409/416 (98.32%) |
| [packages/parigp-ts/src/_error_display.ts](packages/parigp-ts/src/_error_display.ts) | 1 | 0 | 0 | 4/4 (100.00%) |
| [packages/parigp-ts/src/_extension_automorphism.ts](packages/parigp-ts/src/_extension_automorphism.ts) | 1 | 0 | 0 | 97/101 (96.04%) |
| [packages/parigp-ts/src/_extension_composition.ts](packages/parigp-ts/src/_extension_composition.ts) | 1 | 0 | 0 | 112/116 (96.55%) |
| [packages/parigp-ts/src/_extension_derivative.ts](packages/parigp-ts/src/_extension_derivative.ts) | 1 | 0 | 0 | 24/24 (100.00%) |
| [packages/parigp-ts/src/_extension_display.ts](packages/parigp-ts/src/_extension_display.ts) | 4 | 0 | 0 | 53/56 (94.64%) |
| [packages/parigp-ts/src/_extension_division.ts](packages/parigp-ts/src/_extension_division.ts) | 7 | 0 | 0 | 191/196 (97.45%) |
| [packages/parigp-ts/src/_extension_field.ts](packages/parigp-ts/src/_extension_field.ts) | 1 | 0 | 0 | 130/140 (92.86%) |
| [packages/parigp-ts/src/_extension_frobenius.ts](packages/parigp-ts/src/_extension_frobenius.ts) | 2 | 0 | 0 | 84/87 (96.55%) |
| [packages/parigp-ts/src/_extension_gcd.ts](packages/parigp-ts/src/_extension_gcd.ts) | 10 | 0 | 0 | 178/179 (99.44%) |
| [packages/parigp-ts/src/_extension_matrix.ts](packages/parigp-ts/src/_extension_matrix.ts) | 1 | 0 | 0 | 49/52 (94.23%) |
| [packages/parigp-ts/src/_extension_minpoly.ts](packages/parigp-ts/src/_extension_minpoly.ts) | 1 | 0 | 0 | 94/102 (92.16%) |
| [packages/parigp-ts/src/_extension_polynomial.ts](packages/parigp-ts/src/_extension_polynomial.ts) | 8 | 0 | 0 | 197/199 (98.99%) |
| [packages/parigp-ts/src/_extension_projection.ts](packages/parigp-ts/src/_extension_projection.ts) | 1 | 0 | 0 | 61/62 (98.39%) |
| [packages/parigp-ts/src/_extension_quotient.ts](packages/parigp-ts/src/_extension_quotient.ts) | 1 | 0 | 0 | 135/139 (97.12%) |
| [packages/parigp-ts/src/_extension_roots.ts](packages/parigp-ts/src/_extension_roots.ts) | 1 | 0 | 0 | 28/28 (100.00%) |
| [packages/parigp-ts/src/_finite_field_square_root.ts](packages/parigp-ts/src/_finite_field_square_root.ts) | 1 | 0 | 0 | 117/117 (100.00%) |
| [packages/parigp-ts/src/_lll_gso.ts](packages/parigp-ts/src/_lll_gso.ts) | 9 | 0 | 0 | 107/108 (99.07%) |
| [packages/parigp-ts/src/_lll_wrapper.ts](packages/parigp-ts/src/_lll_wrapper.ts) | 16 | 0 | 0 | 337/350 (96.29%) |
| [packages/parigp-ts/src/_matrix_inverse.ts](packages/parigp-ts/src/_matrix_inverse.ts) | 18 | 0 | 0 | 602/616 (97.73%) |
| [packages/parigp-ts/src/_matrix_kernel.ts](packages/parigp-ts/src/_matrix_kernel.ts) | 10 | 0 | 0 | 375/375 (100.00%) |
| [packages/parigp-ts/src/_matrix_mul.ts](packages/parigp-ts/src/_matrix_mul.ts) | 4 | 0 | 0 | 75/77 (97.40%) |
| [packages/parigp-ts/src/_modular_sqrt.ts](packages/parigp-ts/src/_modular_sqrt.ts) | 11 | 0 | 0 | 178/189 (94.18%) |
| [packages/parigp-ts/src/_mpqs_hash.ts](packages/parigp-ts/src/_mpqs_hash.ts) | 6 | 0 | 0 | 65/68 (95.59%) |
| [packages/parigp-ts/src/_polynomial_composition.ts](packages/parigp-ts/src/_polynomial_composition.ts) | 2 | 0 | 0 | 73/75 (97.33%) |
| [packages/parigp-ts/src/_polynomial_ddf.ts](packages/parigp-ts/src/_polynomial_ddf.ts) | 1 | 0 | 0 | 53/54 (98.15%) |
| [packages/parigp-ts/src/_polynomial_division.ts](packages/parigp-ts/src/_polynomial_division.ts) | 5 | 0 | 0 | 148/149 (99.33%) |
| [packages/parigp-ts/src/_polynomial_factor.ts](packages/parigp-ts/src/_polynomial_factor.ts) | 1 | 0 | 0 | 109/118 (92.37%) |
| [packages/parigp-ts/src/_polynomial_gcd.ts](packages/parigp-ts/src/_polynomial_gcd.ts) | 12 | 0 | 0 | 238/246 (96.75%) |
| [packages/parigp-ts/src/_polynomial_minpoly.ts](packages/parigp-ts/src/_polynomial_minpoly.ts) | 1 | 0 | 0 | 75/81 (92.59%) |
| [packages/parigp-ts/src/_polynomial_normalize.ts](packages/parigp-ts/src/_polynomial_normalize.ts) | 1 | 0 | 0 | 24/25 (96.00%) |
| [packages/parigp-ts/src/_polynomial_packing.ts](packages/parigp-ts/src/_polynomial_packing.ts) | 8 | 0 | 0 | 97/100 (97.00%) |
| [packages/parigp-ts/src/_polynomial_quotient.ts](packages/parigp-ts/src/_polynomial_quotient.ts) | 1 | 0 | 0 | 47/48 (97.92%) |
| [packages/parigp-ts/src/_polynomial_quotient_power.ts](packages/parigp-ts/src/_polynomial_quotient_power.ts) | 2 | 0 | 0 | 46/51 (90.20%) |
| [packages/parigp-ts/src/_polynomial_roots.ts](packages/parigp-ts/src/_polynomial_roots.ts) | 4 | 0 | 0 | 170/184 (92.39%) |
| [packages/parigp-ts/src/_rational_polynomial.ts](packages/parigp-ts/src/_rational_polynomial.ts) | 8 | 0 | 0 | 30/31 (96.77%) |
| [packages/parigp-ts/src/_real_matrix.ts](packages/parigp-ts/src/_real_matrix.ts) | 3 | 0 | 0 | 24/24 (100.00%) |
| [packages/parigp-ts/src/_scalar_power.ts](packages/parigp-ts/src/_scalar_power.ts) | 4 | 0 | 0 | 133/144 (92.36%) |
| [packages/parigp-ts/src/alglin1.ts](packages/parigp-ts/src/alglin1.ts) | 12 | 0 | 0 | 181/184 (98.37%) |
| [packages/parigp-ts/src/alglin2.ts](packages/parigp-ts/src/alglin2.ts) | 4 | 0 | 0 | 46/46 (100.00%) |
| [packages/parigp-ts/src/arith1.ts](packages/parigp-ts/src/arith1.ts) | 5 | 0 | 0 | 76/111 (68.47%) |
| [packages/parigp-ts/src/arith2.ts](packages/parigp-ts/src/arith2.ts) | 3 | 0 | 0 | 40/43 (93.02%) |
| [packages/parigp-ts/src/base1.ts](packages/parigp-ts/src/base1.ts) | 1 | 0 | 0 | 56/56 (100.00%) |
| [packages/parigp-ts/src/base2.ts](packages/parigp-ts/src/base2.ts) | 16 | 0 | 0 | 306/312 (98.08%) |
| [packages/parigp-ts/src/base3.ts](packages/parigp-ts/src/base3.ts) | 3 | 0 | 0 | 34/34 (100.00%) |
| [packages/parigp-ts/src/base4.ts](packages/parigp-ts/src/base4.ts) | 20 | 0 | 0 | 338/346 (97.69%) |
| [packages/parigp-ts/src/bb_group.ts](packages/parigp-ts/src/bb_group.ts) | 5 | 0 | 0 | 100/104 (96.15%) |
| [packages/parigp-ts/src/bern.ts](packages/parigp-ts/src/bern.ts) | 5 | 0 | 0 | 100/100 (100.00%) |
| [packages/parigp-ts/src/bibli1.ts](packages/parigp-ts/src/bibli1.ts) | 5 | 0 | 0 | 121/128 (94.53%) |
| [packages/parigp-ts/src/bibli2.ts](packages/parigp-ts/src/bibli2.ts) | 1 | 0 | 0 | 12/12 (100.00%) |
| [packages/parigp-ts/src/buch.ts](packages/parigp-ts/src/buch.ts) | 168 | 0 | 3 | 2166/2373 (91.28%) |
| [packages/parigp-ts/src/char.ts](packages/parigp-ts/src/char.ts) | 1 | 0 | 0 | 63/64 (98.44%) |
| [packages/parigp-ts/src/elliptic/advanced.ts](packages/parigp-ts/src/elliptic/advanced.ts) | 67 | 0 | 4 | 898/1189 (75.53%) |
| [packages/parigp-ts/src/elliptic/ellsea.ts](packages/parigp-ts/src/elliptic/ellsea.ts) | 171 | 0 | 1 | 2008/2157 (93.09%) |
| [packages/parigp-ts/src/elliptic/group.ts](packages/parigp-ts/src/elliptic/group.ts) | 60 | 0 | 0 | 704/868 (81.11%) |
| [packages/parigp-ts/src/elliptic/init.ts](packages/parigp-ts/src/elliptic/init.ts) | 20 | 0 | 0 | 285/412 (69.17%) |
| [packages/parigp-ts/src/elliptic/point.ts](packages/parigp-ts/src/elliptic/point.ts) | 12 | 0 | 0 | 162/253 (64.03%) |
| [packages/parigp-ts/src/elliptic/points.ts](packages/parigp-ts/src/elliptic/points.ts) | 15 | 0 | 0 | 121/187 (64.71%) |
| [packages/parigp-ts/src/elltrans.ts](packages/parigp-ts/src/elltrans.ts) | 3 | 0 | 0 | 44/45 (97.78%) |
| [packages/parigp-ts/src/errors.ts](packages/parigp-ts/src/errors.ts) | 1 | 0 | 0 | 4/4 (100.00%) |
| [packages/parigp-ts/src/ff.ts](packages/parigp-ts/src/ff.ts) | 32 | 0 | 0 | 263/270 (97.41%) |
| [packages/parigp-ts/src/ffinit.ts](packages/parigp-ts/src/ffinit.ts) | 37 | 0 | 0 | 347/396 (87.63%) |
| [packages/parigp-ts/src/galconj.ts](packages/parigp-ts/src/galconj.ts) | 181 | 0 | 3 | 2363/2395 (98.66%) |
| [packages/parigp-ts/src/gen2.ts](packages/parigp-ts/src/gen2.ts) | 4 | 0 | 0 | 76/76 (100.00%) |
| [packages/parigp-ts/src/gen3.ts](packages/parigp-ts/src/gen3.ts) | 3 | 0 | 0 | 33/34 (97.06%) |
| [packages/parigp-ts/src/hnf_snf.ts](packages/parigp-ts/src/hnf_snf.ts) | 16 | 0 | 0 | 463/473 (97.89%) |
| [packages/parigp-ts/src/ifactor.ts](packages/parigp-ts/src/ifactor.ts) | 41 | 0 | 1 | 1022/1139 (89.73%) |
| [packages/parigp-ts/src/ifactor1.ts](packages/parigp-ts/src/ifactor1.ts) | 3 | 0 | 0 | 155/164 (94.51%) |
| [packages/parigp-ts/src/index.ts](packages/parigp-ts/src/index.ts) | 0 | 0 | 0 | 92/92 (100.00%) |
| [packages/parigp-ts/src/ispower.ts](packages/parigp-ts/src/ispower.ts) | 2 | 0 | 0 | 50/53 (94.34%) |
| [packages/parigp-ts/src/kernel/gmp/mp.ts](packages/parigp-ts/src/kernel/gmp/mp.ts) | 8 | 0 | 0 | 190/208 (91.35%) |
| [packages/parigp-ts/src/kernel/gmp/sqrtrem.ts](packages/parigp-ts/src/kernel/gmp/sqrtrem.ts) | 4 | 0 | 0 | 89/89 (100.00%) |
| [packages/parigp-ts/src/kernel/none/add.ts](packages/parigp-ts/src/kernel/none/add.ts) | 1 | 0 | 1 | 21/21 (100.00%) |
| [packages/parigp-ts/src/kernel/none/cmp.ts](packages/parigp-ts/src/kernel/none/cmp.ts) | 2 | 0 | 0 | 12/14 (85.71%) |
| [packages/parigp-ts/src/kernel/none/halfgcd.ts](packages/parigp-ts/src/kernel/none/halfgcd.ts) | 10 | 0 | 0 | 152/152 (100.00%) |
| [packages/parigp-ts/src/kernel/none/level1.ts](packages/parigp-ts/src/kernel/none/level1.ts) | 3 | 0 | 0 | 25/25 (100.00%) |
| [packages/parigp-ts/src/kernel/none/mp_indep.ts](packages/parigp-ts/src/kernel/none/mp_indep.ts) | 13 | 0 | 0 | 212/222 (95.50%) |
| [packages/parigp-ts/src/language/forprime.ts](packages/parigp-ts/src/language/forprime.ts) | 3 | 0 | 0 | 57/60 (95.00%) |
| [packages/parigp-ts/src/lll.ts](packages/parigp-ts/src/lll.ts) | 23 | 0 | 0 | 817/851 (96.00%) |
| [packages/parigp-ts/src/matkermod.ts](packages/parigp-ts/src/matkermod.ts) | 70 | 0 | 0 | 778/847 (91.85%) |
| [packages/parigp-ts/src/mpqs.ts](packages/parigp-ts/src/mpqs.ts) | 45 | 0 | 0 | 1119/1218 (91.87%) |
| [packages/parigp-ts/src/nffactor.ts](packages/parigp-ts/src/nffactor.ts) | 1 | 0 | 0 | 2/2 (100.00%) |
| [packages/parigp-ts/src/polarit2.ts](packages/parigp-ts/src/polarit2.ts) | 5 | 0 | 0 | 115/126 (91.27%) |
| [packages/parigp-ts/src/polarit3.ts](packages/parigp-ts/src/polarit3.ts) | 3 | 0 | 0 | 26/26 (100.00%) |
| [packages/parigp-ts/src/polmodular.ts](packages/parigp-ts/src/polmodular.ts) | 184 | 0 | 9 | 2761/3640 (75.85%) |
| [packages/parigp-ts/src/prime.ts](packages/parigp-ts/src/prime.ts) | 3 | 0 | 0 | 101/108 (93.52%) |
| [packages/parigp-ts/src/qfb.ts](packages/parigp-ts/src/qfb.ts) | 167 | 0 | 4 | 1584/1693 (93.56%) |
| [packages/parigp-ts/src/qfrep.ts](packages/parigp-ts/src/qfrep.ts) | 24 | 0 | 0 | 283/338 (83.73%) |
| [packages/parigp-ts/src/random.ts](packages/parigp-ts/src/random.ts) | 11 | 0 | 0 | 91/91 (100.00%) |
| [packages/parigp-ts/src/subcyclo.ts](packages/parigp-ts/src/subcyclo.ts) | 2 | 0 | 0 | 25/25 (100.00%) |
| [packages/parigp-ts/src/subgroup.ts](packages/parigp-ts/src/subgroup.ts) | 3 | 0 | 0 | 166/171 (97.08%) |
| [packages/parigp-ts/src/trans1.ts](packages/parigp-ts/src/trans1.ts) | 19 | 0 | 0 | 182/238 (76.47%) |
| [packages/parigp-ts/src/trans2.ts](packages/parigp-ts/src/trans2.ts) | 8 | 0 | 0 | 104/107 (97.20%) |
| [packages/parigp-ts/src/types.ts](packages/parigp-ts/src/types.ts) | 24 | 0 | 0 | 61/101 (60.40%) |
| [packages/sagemath-ts/src/algebras/index.ts](packages/sagemath-ts/src/algebras/index.ts) | 0 | 0 | 0 | 1/1 (100.00%) |
| [packages/sagemath-ts/src/algebras/quatalg/index.ts](packages/sagemath-ts/src/algebras/quatalg/index.ts) | 0 | 0 | 0 | 3/3 (100.00%) |
| [packages/sagemath-ts/src/algebras/quatalg/quaternion_algebra.ts](packages/sagemath-ts/src/algebras/quatalg/quaternion_algebra.ts) | 164 | 0 | 12 | 1824/2020 (90.30%) |
| [packages/sagemath-ts/src/algebras/quatalg/quaternion_algebra_cython.ts](packages/sagemath-ts/src/algebras/quatalg/quaternion_algebra_cython.ts) | 4 | 0 | 0 | 79/86 (91.86%) |
| [packages/sagemath-ts/src/algebras/quatalg/quaternion_algebra_element.ts](packages/sagemath-ts/src/algebras/quatalg/quaternion_algebra_element.ts) | 35 | 2 | 0 | 217/248 (87.50%) |
| [packages/sagemath-ts/src/arith/index.ts](packages/sagemath-ts/src/arith/index.ts) | 0 | 0 | 0 | 2/2 (100.00%) |
| [packages/sagemath-ts/src/arith/misc.ts](packages/sagemath-ts/src/arith/misc.ts) | 103 | 0 | 3 | 2034/2887 (70.45%) |
| [packages/sagemath-ts/src/arith/power.ts](packages/sagemath-ts/src/arith/power.ts) | 2 | 0 | 0 | 24/24 (100.00%) |
| [packages/sagemath-ts/src/categories/fields.ts](packages/sagemath-ts/src/categories/fields.ts) | 1 | 0 | 1 | 21/21 (100.00%) |
| [packages/sagemath-ts/src/coding/bch_code.ts](packages/sagemath-ts/src/coding/bch_code.ts) | 25 | 0 | 1 | 592/808 (73.27%) |
| [packages/sagemath-ts/src/coding/goppa_code.ts](packages/sagemath-ts/src/coding/goppa_code.ts) | 22 | 0 | 0 | 541/752 (71.94%) |
| [packages/sagemath-ts/src/coding/index.ts](packages/sagemath-ts/src/coding/index.ts) | 0 | 0 | 0 | 6/6 (100.00%) |
| [packages/sagemath-ts/src/coding/reed_muller_code.ts](packages/sagemath-ts/src/coding/reed_muller_code.ts) | 43 | 0 | 0 | 407/526 (77.38%) |
| [packages/sagemath-ts/src/coding/reed_solomon.ts](packages/sagemath-ts/src/coding/reed_solomon.ts) | 23 | 0 | 0 | 428/582 (73.54%) |
| [packages/sagemath-ts/src/crypto/boolean_function.ts](packages/sagemath-ts/src/crypto/boolean_function.ts) | 41 | 0 | 0 | 526/670 (78.51%) |
| [packages/sagemath-ts/src/crypto/index.ts](packages/sagemath-ts/src/crypto/index.ts) | 0 | 0 | 0 | 5/5 (100.00%) |
| [packages/sagemath-ts/src/crypto/lattice.ts](packages/sagemath-ts/src/crypto/lattice.ts) | 4 | 0 | 0 | 165/251 (65.74%) |
| [packages/sagemath-ts/src/crypto/lwe.ts](packages/sagemath-ts/src/crypto/lwe.ts) | 43 | 0 | 0 | 516/672 (76.79%) |
| [packages/sagemath-ts/src/crypto/sbox.ts](packages/sagemath-ts/src/crypto/sbox.ts) | 52 | 0 | 0 | 440/538 (81.78%) |
| [packages/sagemath-ts/src/errors.ts](packages/sagemath-ts/src/errors.ts) | 9 | 0 | 0 | 54/55 (98.18%) |
| [packages/sagemath-ts/src/groups/generic.ts](packages/sagemath-ts/src/groups/generic.ts) | 30 | 0 | 0 | 883/1211 (72.91%) |
| [packages/sagemath-ts/src/groups/index.ts](packages/sagemath-ts/src/groups/index.ts) | 0 | 0 | 0 | 1/1 (100.00%) |
| [packages/sagemath-ts/src/index.ts](packages/sagemath-ts/src/index.ts) | 0 | 0 | 0 | 30/30 (100.00%) |
| [packages/sagemath-ts/src/matrix/args.ts](packages/sagemath-ts/src/matrix/args.ts) | 2 | 0 | 0 | 71/75 (94.67%) |
| [packages/sagemath-ts/src/matrix/change_ring.ts](packages/sagemath-ts/src/matrix/change_ring.ts) | 1 | 0 | 0 | 12/12 (100.00%) |
| [packages/sagemath-ts/src/matrix/index.ts](packages/sagemath-ts/src/matrix/index.ts) | 0 | 0 | 0 | 17/17 (100.00%) |
| [packages/sagemath-ts/src/matrix/matrix0.ts](packages/sagemath-ts/src/matrix/matrix0.ts) | 1 | 0 | 0 | 111/118 (94.07%) |
| [packages/sagemath-ts/src/matrix/matrix_decompositions.ts](packages/sagemath-ts/src/matrix/matrix_decompositions.ts) | 88 | 1 | 8 | 2215/3492 (63.43%) |
| [packages/sagemath-ts/src/matrix/matrix_decompositions_additions.ts](packages/sagemath-ts/src/matrix/matrix_decompositions_additions.ts) | 15 | 0 | 0 | 449/586 (76.62%) |
| [packages/sagemath-ts/src/matrix/matrix_generic.ts](packages/sagemath-ts/src/matrix/matrix_generic.ts) | 26 | 0 | 0 | 294/323 (91.02%) |
| [packages/sagemath-ts/src/matrix/matrix_integer.ts](packages/sagemath-ts/src/matrix/matrix_integer.ts) | 101 | 0 | 0 | 2333/2960 (78.82%) |
| [packages/sagemath-ts/src/matrix/matrix_mod2.ts](packages/sagemath-ts/src/matrix/matrix_mod2.ts) | 58 | 2 | 1 | 1054/1187 (88.80%) |
| [packages/sagemath-ts/src/matrix/matrix_modn.ts](packages/sagemath-ts/src/matrix/matrix_modn.ts) | 30 | 0 | 2 | 604/818 (73.84%) |
| [packages/sagemath-ts/src/matrix/matrix_operations.ts](packages/sagemath-ts/src/matrix/matrix_operations.ts) | 86 | 1 | 9 | 1591/2374 (67.02%) |
| [packages/sagemath-ts/src/matrix/matrix_polynomial_dense.ts](packages/sagemath-ts/src/matrix/matrix_polynomial_dense.ts) | 47 | 0 | 0 | 1058/1250 (84.64%) |
| [packages/sagemath-ts/src/matrix/matrix_space.ts](packages/sagemath-ts/src/matrix/matrix_space.ts) | 12 | 0 | 0 | 114/130 (87.69%) |
| [packages/sagemath-ts/src/matrix/matrix_special.ts](packages/sagemath-ts/src/matrix/matrix_special.ts) | 78 | 4 | 8 | 1738/2615 (66.46%) |
| [packages/sagemath-ts/src/misc/derivative.ts](packages/sagemath-ts/src/misc/derivative.ts) | 2 | 0 | 0 | 38/39 (97.44%) |
| [packages/sagemath-ts/src/misc/mrange.ts](packages/sagemath-ts/src/misc/mrange.ts) | 1 | 0 | 0 | 68/82 (82.93%) |
| [packages/sagemath-ts/src/misc/parser.ts](packages/sagemath-ts/src/misc/parser.ts) | 30 | 7 | 2 | 249/257 (96.89%) |
| [packages/sagemath-ts/src/misc/randstate.ts](packages/sagemath-ts/src/misc/randstate.ts) | 29 | 0 | 0 | 433/465 (93.12%) |
| [packages/sagemath-ts/src/modules/bkz.ts](packages/sagemath-ts/src/modules/bkz.ts) | 13 | 0 | 0 | 291/600 (48.50%) |
| [packages/sagemath-ts/src/modules/free_module.ts](packages/sagemath-ts/src/modules/free_module.ts) | 241 | 0 | 12 | 2345/3027 (77.47%) |
| [packages/sagemath-ts/src/modules/free_module_element.ts](packages/sagemath-ts/src/modules/free_module_element.ts) | 46 | 0 | 2 | 482/612 (78.76%) |
| [packages/sagemath-ts/src/modules/free_module_integer.ts](packages/sagemath-ts/src/modules/free_module_integer.ts) | 79 | 0 | 3 | 1462/2128 (68.70%) |
| [packages/sagemath-ts/src/modules/index.ts](packages/sagemath-ts/src/modules/index.ts) | 0 | 0 | 0 | 4/4 (100.00%) |
| [packages/sagemath-ts/src/quadratic_forms/binary_qf.ts](packages/sagemath-ts/src/quadratic_forms/binary_qf.ts) | 40 | 0 | 2 | 489/570 (85.79%) |
| [packages/sagemath-ts/src/quadratic_forms/index.ts](packages/sagemath-ts/src/quadratic_forms/index.ts) | 0 | 0 | 0 | 4/4 (100.00%) |
| [packages/sagemath-ts/src/quadratic_forms/quadratic_form.ts](packages/sagemath-ts/src/quadratic_forms/quadratic_form.ts) | 77 | 16 | 2 | 608/655 (92.82%) |
| [packages/sagemath-ts/src/quadratic_forms/quadratic_form__local_field_invariants.ts](packages/sagemath-ts/src/quadratic_forms/quadratic_form__local_field_invariants.ts) | 22 | 0 | 3 | 301/354 (85.03%) |
| [packages/sagemath-ts/src/quadratic_forms/ternary_qf.ts](packages/sagemath-ts/src/quadratic_forms/ternary_qf.ts) | 61 | 5 | 0 | 958/1049 (91.33%) |
| [packages/sagemath-ts/src/rings/complex_mpfr.ts](packages/sagemath-ts/src/rings/complex_mpfr.ts) | 73 | 0 | 1 | 720/953 (75.55%) |
| [packages/sagemath-ts/src/rings/fast_arith.ts](packages/sagemath-ts/src/rings/fast_arith.ts) | 10 | 0 | 0 | 68/68 (100.00%) |
| [packages/sagemath-ts/src/rings/finite_rings/conway_polynomials.ts](packages/sagemath-ts/src/rings/finite_rings/conway_polynomials.ts) | 4 | 0 | 0 | 265/267 (99.25%) |
| [packages/sagemath-ts/src/rings/finite_rings/element_base.ts](packages/sagemath-ts/src/rings/finite_rings/element_base.ts) | 1 | 0 | 0 | 149/175 (85.14%) |
| [packages/sagemath-ts/src/rings/finite_rings/finite_field_base.ts](packages/sagemath-ts/src/rings/finite_rings/finite_field_base.ts) | 1 | 0 | 0 | 28/28 (100.00%) |
| [packages/sagemath-ts/src/rings/finite_rings/finite_field_constructor.ts](packages/sagemath-ts/src/rings/finite_rings/finite_field_constructor.ts) | 4 | 0 | 0 | 52/69 (75.36%) |
| [packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.ts](packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.ts) | 75 | 0 | 2 | 843/865 (97.46%) |
| [packages/sagemath-ts/src/rings/finite_rings/finite_field_prime.ts](packages/sagemath-ts/src/rings/finite_rings/finite_field_prime.ts) | 37 | 0 | 0 | 200/235 (85.11%) |
| [packages/sagemath-ts/src/rings/finite_rings/gf2.ts](packages/sagemath-ts/src/rings/finite_rings/gf2.ts) | 28 | 0 | 0 | 101/108 (93.52%) |
| [packages/sagemath-ts/src/rings/finite_rings/index.ts](packages/sagemath-ts/src/rings/finite_rings/index.ts) | 0 | 0 | 0 | 11/11 (100.00%) |
| [packages/sagemath-ts/src/rings/finite_rings/integer_mod.ts](packages/sagemath-ts/src/rings/finite_rings/integer_mod.ts) | 30 | 0 | 0 | 374/437 (85.58%) |
| [packages/sagemath-ts/src/rings/finite_rings/integer_mod_ring.ts](packages/sagemath-ts/src/rings/finite_rings/integer_mod_ring.ts) | 20 | 0 | 1 | 211/233 (90.56%) |
| [packages/sagemath-ts/src/rings/finite_rings/roots_of_unity.ts](packages/sagemath-ts/src/rings/finite_rings/roots_of_unity.ts) | 31 | 0 | 0 | 363/577 (62.91%) |
| [packages/sagemath-ts/src/rings/finite_rings/tower_field.ts](packages/sagemath-ts/src/rings/finite_rings/tower_field.ts) | 6 | 0 | 0 | 44/78 (56.41%) |
| [packages/sagemath-ts/src/rings/fraction_field.ts](packages/sagemath-ts/src/rings/fraction_field.ts) | 22 | 0 | 0 | 208/228 (91.23%) |
| [packages/sagemath-ts/src/rings/fraction_field_FpT.ts](packages/sagemath-ts/src/rings/fraction_field_FpT.ts) | 20 | 0 | 1 | 179/182 (98.35%) |
| [packages/sagemath-ts/src/rings/fraction_field_element.ts](packages/sagemath-ts/src/rings/fraction_field_element.ts) | 30 | 0 | 2 | 270/289 (93.43%) |
| [packages/sagemath-ts/src/rings/function_field/constant_field.ts](packages/sagemath-ts/src/rings/function_field/constant_field.ts) | 8 | 0 | 0 | 53/92 (57.61%) |
| [packages/sagemath-ts/src/rings/function_field/constructor.ts](packages/sagemath-ts/src/rings/function_field/constructor.ts) | 2 | 1 | 0 | 9/14 (64.29%) |
| [packages/sagemath-ts/src/rings/function_field/divisor.ts](packages/sagemath-ts/src/rings/function_field/divisor.ts) | 37 | 1 | 1 | 344/395 (87.09%) |
| [packages/sagemath-ts/src/rings/function_field/element.ts](packages/sagemath-ts/src/rings/function_field/element.ts) | 14 | 4 | 0 | 95/115 (82.61%) |
| [packages/sagemath-ts/src/rings/function_field/element_rational.ts](packages/sagemath-ts/src/rings/function_field/element_rational.ts) | 39 | 0 | 2 | 284/316 (89.87%) |
| [packages/sagemath-ts/src/rings/function_field/function_field.ts](packages/sagemath-ts/src/rings/function_field/function_field.ts) | 20 | 8 | 0 | 98/111 (88.29%) |
| [packages/sagemath-ts/src/rings/function_field/function_field_polymod.ts](packages/sagemath-ts/src/rings/function_field/function_field_polymod.ts) | 2 | 2 | 0 | 2/10 (20.00%) |
| [packages/sagemath-ts/src/rings/function_field/function_field_rational.ts](packages/sagemath-ts/src/rings/function_field/function_field_rational.ts) | 31 | 1 | 0 | 191/220 (86.82%) |
| [packages/sagemath-ts/src/rings/function_field/ideal.ts](packages/sagemath-ts/src/rings/function_field/ideal.ts) | 23 | 0 | 1 | 121/130 (93.08%) |
| [packages/sagemath-ts/src/rings/function_field/ideal_rational.ts](packages/sagemath-ts/src/rings/function_field/ideal_rational.ts) | 28 | 0 | 0 | 96/96 (100.00%) |
| [packages/sagemath-ts/src/rings/function_field/index.ts](packages/sagemath-ts/src/rings/function_field/index.ts) | 0 | 0 | 0 | 15/15 (100.00%) |
| [packages/sagemath-ts/src/rings/function_field/order.ts](packages/sagemath-ts/src/rings/function_field/order.ts) | 12 | 0 | 1 | 40/44 (90.91%) |
| [packages/sagemath-ts/src/rings/function_field/order_rational.ts](packages/sagemath-ts/src/rings/function_field/order_rational.ts) | 22 | 1 | 1 | 151/176 (85.80%) |
| [packages/sagemath-ts/src/rings/function_field/place.ts](packages/sagemath-ts/src/rings/function_field/place.ts) | 19 | 0 | 0 | 60/64 (93.75%) |
| [packages/sagemath-ts/src/rings/function_field/place_rational.ts](packages/sagemath-ts/src/rings/function_field/place_rational.ts) | 6 | 0 | 1 | 76/91 (83.52%) |
| [packages/sagemath-ts/src/rings/function_field/valuation_ring.ts](packages/sagemath-ts/src/rings/function_field/valuation_ring.ts) | 6 | 0 | 0 | 42/45 (93.33%) |
| [packages/sagemath-ts/src/rings/generic.ts](packages/sagemath-ts/src/rings/generic.ts) | 10 | 0 | 0 | 87/90 (96.67%) |
| [packages/sagemath-ts/src/rings/index.ts](packages/sagemath-ts/src/rings/index.ts) | 0 | 0 | 0 | 15/15 (100.00%) |
| [packages/sagemath-ts/src/rings/integer_ring.ts](packages/sagemath-ts/src/rings/integer_ring.ts) | 128 | 0 | 3 | 1068/1290 (82.79%) |
| [packages/sagemath-ts/src/rings/laurent_series_ring.ts](packages/sagemath-ts/src/rings/laurent_series_ring.ts) | 68 | 0 | 0 | 526/594 (88.55%) |
| [packages/sagemath-ts/src/rings/number_field/class_group.ts](packages/sagemath-ts/src/rings/number_field/class_group.ts) | 50 | 0 | 3 | 469/569 (82.43%) |
| [packages/sagemath-ts/src/rings/number_field/galois_group.ts](packages/sagemath-ts/src/rings/number_field/galois_group.ts) | 70 | 0 | 5 | 650/1087 (59.80%) |
| [packages/sagemath-ts/src/rings/number_field/index.ts](packages/sagemath-ts/src/rings/number_field/index.ts) | 0 | 0 | 0 | 8/8 (100.00%) |
| [packages/sagemath-ts/src/rings/number_field/number_field.ts](packages/sagemath-ts/src/rings/number_field/number_field.ts) | 146 | 2 | 11 | 1544/1661 (92.96%) |
| [packages/sagemath-ts/src/rings/number_field/number_field_element.ts](packages/sagemath-ts/src/rings/number_field/number_field_element.ts) | 8 | 0 | 1 | 13/88 (14.77%) |
| [packages/sagemath-ts/src/rings/number_field/number_field_embeddings.ts](packages/sagemath-ts/src/rings/number_field/number_field_embeddings.ts) | 64 | 0 | 1 | 567/618 (91.75%) |
| [packages/sagemath-ts/src/rings/number_field/number_field_ideal.ts](packages/sagemath-ts/src/rings/number_field/number_field_ideal.ts) | 43 | 0 | 4 | 478/573 (83.42%) |
| [packages/sagemath-ts/src/rings/number_field/order.ts](packages/sagemath-ts/src/rings/number_field/order.ts) | 39 | 0 | 4 | 227/353 (64.31%) |
| [packages/sagemath-ts/src/rings/number_field/pari_nf.ts](packages/sagemath-ts/src/rings/number_field/pari_nf.ts) | 80 | 0 | 1 | 1376/1440 (95.56%) |
| [packages/sagemath-ts/src/rings/number_field/unit_group.ts](packages/sagemath-ts/src/rings/number_field/unit_group.ts) | 36 | 0 | 7 | 367/488 (75.20%) |
| [packages/sagemath-ts/src/rings/padics/index.ts](packages/sagemath-ts/src/rings/padics/index.ts) | 0 | 0 | 0 | 2/2 (100.00%) |
| [packages/sagemath-ts/src/rings/padics/padic_generic.ts](packages/sagemath-ts/src/rings/padics/padic_generic.ts) | 48 | 0 | 0 | 179/251 (71.31%) |
| [packages/sagemath-ts/src/rings/padics/padic_generic_element.ts](packages/sagemath-ts/src/rings/padics/padic_generic_element.ts) | 64 | 0 | 0 | 995/1389 (71.63%) |
| [packages/sagemath-ts/src/rings/polynomial/convolution.ts](packages/sagemath-ts/src/rings/polynomial/convolution.ts) | 34 | 0 | 1 | 560/767 (73.01%) |
| [packages/sagemath-ts/src/rings/polynomial/index.ts](packages/sagemath-ts/src/rings/polynomial/index.ts) | 0 | 0 | 0 | 10/10 (100.00%) |
| [packages/sagemath-ts/src/rings/polynomial/multi_polynomial_element.ts](packages/sagemath-ts/src/rings/polynomial/multi_polynomial_element.ts) | 71 | 19 | 0 | 624/744 (83.87%) |
| [packages/sagemath-ts/src/rings/polynomial/multi_polynomial_ideal.ts](packages/sagemath-ts/src/rings/polynomial/multi_polynomial_ideal.ts) | 19 | 0 | 1 | 345/346 (99.71%) |
| [packages/sagemath-ts/src/rings/polynomial/multi_polynomial_ring.ts](packages/sagemath-ts/src/rings/polynomial/multi_polynomial_ring.ts) | 14 | 0 | 1 | 224/281 (79.72%) |
| [packages/sagemath-ts/src/rings/polynomial/polynomial_compiled.ts](packages/sagemath-ts/src/rings/polynomial/polynomial_compiled.ts) | 4 | 0 | 1 | 158/167 (94.61%) |
| [packages/sagemath-ts/src/rings/polynomial/polynomial_element.ts](packages/sagemath-ts/src/rings/polynomial/polynomial_element.ts) | 185 | 0 | 12 | 5078/5532 (91.79%) |
| [packages/sagemath-ts/src/rings/polynomial/polynomial_gf2x.ts](packages/sagemath-ts/src/rings/polynomial/polynomial_gf2x.ts) | 71 | 0 | 0 | 404/448 (90.18%) |
| [packages/sagemath-ts/src/rings/polynomial/polynomial_ring.ts](packages/sagemath-ts/src/rings/polynomial/polynomial_ring.ts) | 33 | 0 | 0 | 798/970 (82.27%) |
| [packages/sagemath-ts/src/rings/polynomial/quotient_ring.ts](packages/sagemath-ts/src/rings/polynomial/quotient_ring.ts) | 27 | 0 | 0 | 196/268 (73.13%) |
| [packages/sagemath-ts/src/rings/power_series_ring.ts](packages/sagemath-ts/src/rings/power_series_ring.ts) | 109 | 0 | 3 | 1300/1655 (78.55%) |
| [packages/sagemath-ts/src/rings/rational.ts](packages/sagemath-ts/src/rings/rational.ts) | 87 | 0 | 2 | 699/878 (79.61%) |
| [packages/sagemath-ts/src/rings/rational_field.ts](packages/sagemath-ts/src/rings/rational_field.ts) | 44 | 0 | 0 | 335/499 (67.13%) |
| [packages/sagemath-ts/src/rings/real_double.ts](packages/sagemath-ts/src/rings/real_double.ts) | 18 | 0 | 1 | 87/93 (93.55%) |
| [packages/sagemath-ts/src/rings/real_mpfr.ts](packages/sagemath-ts/src/rings/real_mpfr.ts) | 140 | 0 | 1 | 1244/1567 (79.39%) |
| [packages/sagemath-ts/src/rings/real_mpfr_dd.ts](packages/sagemath-ts/src/rings/real_mpfr_dd.ts) | 26 | 0 | 0 | 285/380 (75.00%) |
| [packages/sagemath-ts/src/rings/sum_of_squares.ts](packages/sagemath-ts/src/rings/sum_of_squares.ts) | 7 | 0 | 0 | 107/127 (84.25%) |
| [packages/sagemath-ts/src/schemes/elliptic_curves/cm.ts](packages/sagemath-ts/src/schemes/elliptic_curves/cm.ts) | 14 | 1 | 8 | 530/687 (77.15%) |
| [packages/sagemath-ts/src/schemes/elliptic_curves/constructor.ts](packages/sagemath-ts/src/schemes/elliptic_curves/constructor.ts) | 3 | 0 | 0 | 71/102 (69.61%) |
| [packages/sagemath-ts/src/schemes/elliptic_curves/ell_curve_isogeny.ts](packages/sagemath-ts/src/schemes/elliptic_curves/ell_curve_isogeny.ts) | 135 | 0 | 10 | 1980/2910 (68.04%) |
| [packages/sagemath-ts/src/schemes/elliptic_curves/ell_finite_field.ts](packages/sagemath-ts/src/schemes/elliptic_curves/ell_finite_field.ts) | 102 | 1 | 6 | 1197/1709 (70.04%) |
| [packages/sagemath-ts/src/schemes/elliptic_curves/ell_generic.ts](packages/sagemath-ts/src/schemes/elliptic_curves/ell_generic.ts) | 62 | 2 | 7 | 922/960 (96.04%) |
| [packages/sagemath-ts/src/schemes/elliptic_curves/ell_point.ts](packages/sagemath-ts/src/schemes/elliptic_curves/ell_point.ts) | 43 | 0 | 10 | 614/797 (77.04%) |
| [packages/sagemath-ts/src/schemes/elliptic_curves/ell_torsion.ts](packages/sagemath-ts/src/schemes/elliptic_curves/ell_torsion.ts) | 29 | 0 | 2 | 575/854 (67.33%) |
| [packages/sagemath-ts/src/schemes/elliptic_curves/formal_group.ts](packages/sagemath-ts/src/schemes/elliptic_curves/formal_group.ts) | 64 | 0 | 1 | 445/625 (71.20%) |
| [packages/sagemath-ts/src/schemes/elliptic_curves/index.ts](packages/sagemath-ts/src/schemes/elliptic_curves/index.ts) | 0 | 0 | 0 | 19/19 (100.00%) |
| [packages/sagemath-ts/src/schemes/elliptic_curves/isogeny_class.ts](packages/sagemath-ts/src/schemes/elliptic_curves/isogeny_class.ts) | 62 | 0 | 5 | 1334/1657 (80.51%) |
| [packages/sagemath-ts/src/schemes/elliptic_curves/padic_lseries.ts](packages/sagemath-ts/src/schemes/elliptic_curves/padic_lseries.ts) | 77 | 4 | 10 | 806/1198 (67.28%) |
| [packages/sagemath-ts/src/schemes/elliptic_curves/types.ts](packages/sagemath-ts/src/schemes/elliptic_curves/types.ts) | 1 | 0 | 0 | 20/20 (100.00%) |
| [packages/sagemath-ts/src/schemes/elliptic_curves/weierstrass_morphism.ts](packages/sagemath-ts/src/schemes/elliptic_curves/weierstrass_morphism.ts) | 32 | 0 | 3 | 490/517 (94.78%) |
| [packages/sagemath-ts/src/schemes/hyperelliptic_curves/constructor.ts](packages/sagemath-ts/src/schemes/hyperelliptic_curves/constructor.ts) | 4 | 1 | 1 | 97/134 (72.39%) |
| [packages/sagemath-ts/src/schemes/hyperelliptic_curves/field_ops.ts](packages/sagemath-ts/src/schemes/hyperelliptic_curves/field_ops.ts) | 17 | 0 | 8 | 213/347 (61.38%) |
| [packages/sagemath-ts/src/schemes/hyperelliptic_curves/hyperelliptic_finite_field.ts](packages/sagemath-ts/src/schemes/hyperelliptic_curves/hyperelliptic_finite_field.ts) | 35 | 7 | 0 | 598/658 (90.88%) |
| [packages/sagemath-ts/src/schemes/hyperelliptic_curves/hyperelliptic_g2.ts](packages/sagemath-ts/src/schemes/hyperelliptic_curves/hyperelliptic_g2.ts) | 27 | 3 | 0 | 77/104 (74.04%) |
| [packages/sagemath-ts/src/schemes/hyperelliptic_curves/hyperelliptic_generic.ts](packages/sagemath-ts/src/schemes/hyperelliptic_curves/hyperelliptic_generic.ts) | 42 | 7 | 3 | 340/395 (86.08%) |
| [packages/sagemath-ts/src/schemes/hyperelliptic_curves/hyperelliptic_rational_field.ts](packages/sagemath-ts/src/schemes/hyperelliptic_curves/hyperelliptic_rational_field.ts) | 2 | 2 | 0 | 12/13 (92.31%) |
| [packages/sagemath-ts/src/schemes/hyperelliptic_curves/index.ts](packages/sagemath-ts/src/schemes/hyperelliptic_curves/index.ts) | 0 | 0 | 0 | 11/11 (100.00%) |
| [packages/sagemath-ts/src/schemes/hyperelliptic_curves/invariants.ts](packages/sagemath-ts/src/schemes/hyperelliptic_curves/invariants.ts) | 16 | 0 | 1 | 231/252 (91.67%) |
| [packages/sagemath-ts/src/schemes/hyperelliptic_curves/jacobian_g2.ts](packages/sagemath-ts/src/schemes/hyperelliptic_curves/jacobian_g2.ts) | 1 | 1 | 0 | 8/8 (100.00%) |
| [packages/sagemath-ts/src/schemes/hyperelliptic_curves/jacobian_generic.ts](packages/sagemath-ts/src/schemes/hyperelliptic_curves/jacobian_generic.ts) | 13 | 2 | 1 | 69/71 (97.18%) |
| [packages/sagemath-ts/src/schemes/hyperelliptic_curves/jacobian_homset.ts](packages/sagemath-ts/src/schemes/hyperelliptic_curves/jacobian_homset.ts) | 9 | 0 | 0 | 88/95 (92.63%) |
| [packages/sagemath-ts/src/schemes/hyperelliptic_curves/jacobian_morphism.ts](packages/sagemath-ts/src/schemes/hyperelliptic_curves/jacobian_morphism.ts) | 18 | 0 | 0 | 241/275 (87.64%) |
| [packages/sagemath-ts/src/schemes/index.ts](packages/sagemath-ts/src/schemes/index.ts) | 0 | 0 | 0 | 2/2 (100.00%) |
| [packages/sagemath-ts/src/stats/distributions/discrete_gaussian_integer.ts](packages/sagemath-ts/src/stats/distributions/discrete_gaussian_integer.ts) | 31 | 0 | 1 | 452/526 (85.93%) |
| [packages/sagemath-ts/src/stats/distributions/discrete_gaussian_lattice.ts](packages/sagemath-ts/src/stats/distributions/discrete_gaussian_lattice.ts) | 89 | 0 | 2 | 1190/1352 (88.02%) |
| [packages/sagemath-ts/src/stats/distributions/index.ts](packages/sagemath-ts/src/stats/distributions/index.ts) | 0 | 0 | 0 | 2/2 (100.00%) |
| [packages/sagemath-ts/src/stats/index.ts](packages/sagemath-ts/src/stats/index.ts) | 0 | 0 | 0 | 1/1 (100.00%) |
| [packages/sagemath-ts/src/types/coercion.ts](packages/sagemath-ts/src/types/coercion.ts) | 3 | 0 | 0 | 22/23 (95.65%) |
| [packages/sagemath-ts/src/types/gmp.ts](packages/sagemath-ts/src/types/gmp.ts) | 1 | 0 | 0 | 29/31 (93.55%) |
| [packages/sagemath-ts/src/types/gmp_factorial.ts](packages/sagemath-ts/src/types/gmp_factorial.ts) | 6 | 0 | 0 | 166/171 (97.08%) |
| [packages/sagemath-ts/src/types/index.ts](packages/sagemath-ts/src/types/index.ts) | 0 | 0 | 0 | 1/1 (100.00%) |
| [packages/sagemath-ts/src/types/python_float.ts](packages/sagemath-ts/src/types/python_float.ts) | 3 | 0 | 0 | 97/101 (96.04%) |
| [packages/sagemath-ts/src/types/python_float_data.ts](packages/sagemath-ts/src/types/python_float_data.ts) | 0 | 0 | 0 | 720/720 (100.00%) |
| [packages/sagemath-ts/src/types/python_sort.ts](packages/sagemath-ts/src/types/python_sort.ts) | 7 | 0 | 0 | 258/272 (94.85%) |
| [packages/sagemath-ts/src/types/real_format.ts](packages/sagemath-ts/src/types/real_format.ts) | 3 | 0 | 0 | 35/36 (97.22%) |
| [packages/sagemath-ts/src/zk/index.ts](packages/sagemath-ts/src/zk/index.ts) | 0 | 0 | 0 | 3/3 (100.00%) |
| [packages/sagemath-ts/src/zk/multilinear.ts](packages/sagemath-ts/src/zk/multilinear.ts) | 10 | 0 | 0 | 143/180 (79.44%) |
| [packages/sagemath-ts/src/zk/polynomial_commitment.ts](packages/sagemath-ts/src/zk/polynomial_commitment.ts) | 22 | 0 | 0 | 325/448 (72.54%) |
| [packages/sagemath-ts/src/zk/sumcheck.ts](packages/sagemath-ts/src/zk/sumcheck.ts) | 8 | 0 | 0 | 210/308 (68.18%) |
| [packages/zksecurity-cheatsheets/src/curves.ts](packages/zksecurity-cheatsheets/src/curves.ts) | 2 | 0 | 0 | unmeasured |
| [packages/zksecurity-cheatsheets/src/ethereum.ts](packages/zksecurity-cheatsheets/src/ethereum.ts) | 0 | 0 | 0 | unmeasured |
| [packages/zksecurity-cheatsheets/src/fields.ts](packages/zksecurity-cheatsheets/src/fields.ts) | 0 | 0 | 0 | unmeasured |
| [packages/zksecurity-cheatsheets/src/gcm-key-commitment.ts](packages/zksecurity-cheatsheets/src/gcm-key-commitment.ts) | 29 | 0 | 0 | 317/366 (86.61%) |
| [packages/zksecurity-cheatsheets/src/index.ts](packages/zksecurity-cheatsheets/src/index.ts) | 0 | 0 | 0 | unmeasured |
| [packages/zksecurity-cheatsheets/src/poseidon.ts](packages/zksecurity-cheatsheets/src/poseidon.ts) | 0 | 0 | 0 | unmeasured |
| [packages/zksecurity-cheatsheets/src/security.ts](packages/zksecurity-cheatsheets/src/security.ts) | 0 | 0 | 0 | unmeasured |

## Integer methods missing a comparative dispatcher

None. Method dispatch coverage is distinct from complete input-domain or algorithm coverage.

## Binary matrix API missing a comparative dispatcher

None. Dispatch coverage is distinct from full input-domain coverage.

## Modular-integer/ring members missing a comparative dispatcher

None. Broader input-domain and category gaps remain open.

## Rational members missing a comparative dispatcher

None. Dispatch coverage does not establish full input-domain coverage.

## RationalField members missing a comparative dispatcher

None.

## IntegerRing members missing a comparative dispatcher

None.

## Finite-field members missing a comparative dispatcher

None. Dispatch coverage is not full input-domain coverage.
