# Project Scope

This document tracks implementation progress. Update this file when completing modules.

## Status Legend

- ⬜ Not started
- 🟡 In progress (note: include agent/person working on it)
- ✅ Complete (with test coverage %)
- 🔴 Blocked (note: reason)

**Maintenance notes:**
- 2026-01-30 consolidated deviations log into root `DEVIATIONS.md`.
- 2026-07-28 reconciled every status marker against `AUDIT-2026-07.md` (370 confirmed findings)
  and the fix pass that followed. Percentages in this file previously tracked *API surface*,
  not fidelity: several modules marked ✅ 100% had reducible Conway entries, wrong Vélu
  formulas, or no test file at all. Coverage figures below now mean "ported **and** verified
  against upstream doctests or an executed oracle". Modules whose ✅ was contradicted by the
  audit are downgraded here even where the audited defect is fixed, when a documented gap
  remains.
- 2026-07-28 (0.0.12) updated for the **deferred-work pass**: the items the audit fix pass had
  left as honest `NotImplementedError` stubs. Several dependencies that were stubbed are now
  implemented (PARI `ffinit`, `matkermod`, the `Qfb` family, the real `Z_factor` chain, GMP's
  MT19937 seeding, `dgs_bern.c`), so the modules that were blocked on them moved up. Anything
  still throwing is recorded as throwing — no percentage was raised for work that only *changed
  the error message*.

---

## Phase 1: Core Number Theory (Cryptography Focus)

### `sage.rings.integer` - Arbitrary Precision Integers
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | Reviewed vs reference (2026-02-04); re-audited 2026-07-27 (`AUDIT-2026-07.md`) |
| `integer.py` | ✅ 98% | ✅ | Core Integer class - 47+ methods |
| `integer_ring.py` | ✅ 98% | ✅ | ZZ ring with full number-theoretic operations. 2026-07: floor division / `ZeroDivisionError` messages, `nth_root_mod`, `is_discriminant`, `is_fundamental_discriminant`, `real_log` above 2^53, `ndigits(0)` fixed; `multiplicative_order` no longer takes a modulus (Sage's signature) |
| `rational.py` | ✅ 98% | ✅ | Rational numbers - 60+ methods. 2026-07: `integerNthRoot` non-termination, `period()` O(order) loop, `norm`/`trace`/`list` fixed |
| `rational_field.py` | ✅ 98% | ✅ 322 tests | QQ field with iteration, Selmer groups, quadratic defect. 2026-07: `quadratic_defect` rewritten; `rational_field.test.ts` created (the module had **no** test file) |

**Integer Methods Implemented:**
- **Roots:** nth_root, exact_log, sqrtrem, is_perfect_power, is_prime_power
- **Combinatorial:** binomial, factorial, bell_number, catalan_number, fibonacci, lucas_number
- **Divisibility:** divisors, prime_divisors, number_of_divisors, is_squarefree, squarefree_part, radical, core
- **Arithmetic Functions:** euler_phi, sigma, moebius, carmichael_lambda
- **Modular:** inverse_mod, powermod, sqrt_mod, multiplicative_order, is_primitive_root, primitive_root
- **Symbols:** jacobi, kronecker, legendre_symbol, valuation, is_discriminant
- **Primality:** is_prime, is_pseudoprime, nth_prime, prime_pi
- **Partitions:** number_of_partitions

**Rational Methods Implemented:**
- **Construction:** from bigint, number, string ("n/d" format), decimal strings
- **Arithmetic:** add, sub, mul, div, neg, inv, pow, abs
- **Comparison:** eq, lt, le, gt, ge, cmp
- **Conversion:** toString, toNumber, floor, ceil, round (6 modes), trunc
- **Predicates:** isZero, isOne, isInteger, isPositive, isNegative, is_unit, is_integral
- **Continued Fractions:** continued_fraction_list (std, hj), continued_fraction
- **Valuations:** valuation, ord, local_height, global_height, padic_valuation
- **Roots:** is_square, sqrt, is_nth_power, nth_root
- **S-units:** support, prime_to_S_part, val_unit, is_S_unit, is_S_integral
- **Algebraic:** minpoly, charpoly, norm, trace, real, imag, conjugate
- **Combinatorial:** factorial, gamma
- **Misc:** period, content, rational_gcd, rational_lcm, height, ndigits, nbits

**RationalField (QQ) Methods Implemented:**
- **Field Properties:** is_field, is_ring, is_integral_domain, is_prime_field, is_absolute, is_finite
- **Algebraic:** characteristic, degree, absolute_degree, ngens, gens, gen, order
- **Discriminants:** discriminant, absolute_discriminant, relative_discriminant
- **Number Field:** maximal_order, ring_of_integers, number_field, power_basis, class_number, signature
- **Iteration:** Symbol.iterator (by height), range_by_height, primes_of_bounded_norm_iter
- **Selmer Groups:** selmer_generators, selmer_group_iterator
- **Quadratic Forms:** quadratic_defect
- **Elements:** zero, one, __call__ (coercion), random_element, an_element, some_elements, zeta

### `sage.rings.real_mpfr` - Real Numbers
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | 2026-07 (`AUDIT-2026-07.md`): 17 findings, all fixed |
| `real_mpfr.ts` | ✅ 95% | ✅ 165 tests | RealField, RealNumber with IEEE 754 double precision. 2026-07: `round` half-away-from-zero, `ulp`, `epsilon`, `fp_rank`, `simplest_rational`, `nearby_rational`, `zeta` on (0,1), `erfc` tail, Bessel `jn` accuracy fixed |
| `complex_mpfr.ts` | ✅ 90% | ✅ | ComplexField. 2026-07: `sqrt` was returning √(conj z); `dilog`, `gamma_inc`, `algebraic_dependency` (complex branch) implemented — see DEVIATIONS.md |

**Real Number Features (using JavaScript Math):**
- Trigonometric: sin, cos, tan, arcsin, arccos, arctan, sincos
- Hyperbolic: sinh, cosh, tanh, arcsinh, arccosh, arctanh, coth, sech, csch
- Exponential/logarithmic: exp, exp2, exp10, expm1, log, log2, log10, log1p
- Power functions: sqrt, cube_root, nth_root, pow
- Constants: pi, euler_constant, catalan_constant, log2
- Special: agm, erf, erfc, gamma, log_gamma, zeta, eint
- Bessel: j0, j1, jn, y0, y1, yn
- Float representation: sign_mantissa_exponent, fp_rank, ulp, epsilon
- Rational conversion: exact_rational, nearby_rational
- Navigation: nextabove, nextbelow, nexttoward
- 85 tests, all passing

**Deviation:** Uses IEEE 754 double precision (53-bit mantissa) instead of MPFR arbitrary precision. See DEVIATIONS.md for details.

### `sage.arith` - Basic Arithmetic
| Function | Status | Tests | Notes |
|----------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | Reviewed 2026-02-04; re-audited 2026-07-27 — 18 findings, all fixed (246 tests in `arith/misc.test.ts`, up from 187) |
| `gcd` | ✅ 100% | ✅ | Binary GCD (Stein's algorithm) |
| `lcm` | ✅ 100% | ✅ | |
| `xgcd` | ✅ 100% | ✅ | Extended Euclidean algorithm |
| `factor` | ✅ 95% | ✅ | Delegates to PARI `Z_factor`, which is now the real `ifac_crack` chain (trial division, pure powers, SQUFOF, Pollard-Brent, ECM). **MPQS is still unported**, so a hard semiprime with a >= 25-digit smallest factor raises `NotImplementedError` naming `mpqs.c` instead of the previous behaviour of returning the composite as if prime |
| `is_prime` | ✅ 100% | ✅ | Delegates to parigp-ts BPSW — **probabilistic**, not the APRCL/ECPP proof Sage's `proof=True` gives |
| `is_pseudoprime` | ✅ 100% | ✅ | Same BPSW entry point as `is_prime`, so the two coincide here |
| `is_prime_power` | ✅ 100% | ✅ | PARI `isprimepower`. A real port now lives in `parigp-ts/src/ifactor.ts` (exact integer k-th roots + BPSW, never factors `n`), but it is **not re-exported from that package's barrel**, so `arith/misc.ts:836` still keeps a duplicate local copy — see DEVIATIONS.md |
| `next_prime` | ✅ 100% | ✅ | |
| `previous_prime` | ✅ 100% | ✅ | |
| `prime_range` | ✅ 100% | ✅ | List of primes in range |
| `euler_phi` | ✅ 100% | ✅ | Euler's totient |
| `moebius` | ✅ 100% | ✅ | Mobius function |
| `carmichael_lambda` | ✅ 100% | ✅ | Carmichael function lambda(n) |
| `primitive_root` | ✅ 100% | ✅ | Primitive root modulo n |
| `quadratic_residues` | ✅ 100% | ✅ | List of QRs mod n |
| `hilbert_symbol` | ✅ 100% | ✅ | Hilbert symbol |
| `two_squares` | ✅ 100% | ✅ | Sum of two squares |
| `four_squares` | ✅ 100% | ✅ | Lagrange's theorem |
| `power_mod` | ✅ 100% | ✅ | Binary exponentiation |
| `inverse_mod` | ✅ 100% | ✅ | Via xgcd |
| `crt` | ✅ 100% | ✅ | Chinese Remainder Theorem |
| `CRT_list` | ✅ 100% | ✅ | CRT for list of residues/moduli |
| `kronecker` | ✅ 100% | ✅ | Kronecker symbol |
| `jacobi_symbol` | ✅ 100% | ✅ | |
| `legendre_symbol` | ✅ 100% | ✅ | |
| `sqrt_mod` | ✅ 100% | ✅ | Tonelli-Shanks algorithm with all_roots option |
| `isqrt` | ✅ 100% | ✅ | Integer square root (Newton) |
| `is_square` | ✅ 100% | ✅ | With optional root extraction |
| `is_squarefree` | ✅ 100% | ✅ | |
| `divisors` | ✅ 100% | ✅ | |
| `number_of_divisors` | ✅ 100% | ✅ | |
| `sigma` | ✅ 100% | ✅ | Sum of k-th powers of divisors |
| `radical` | ✅ 100% | ✅ | |
| `trial_division` | ✅ 100% | ✅ | With optional bound |
| `is_strong_probable_prime` | ✅ 100% | ✅ | Miller-Rabin witness test |
| `rational_reconstruction` | ✅ 100% | ✅ | Reconstruct p/q from a mod m |
| `CRT_basis` | ✅ 100% | ✅ | CRT basis elements e_i |
| `CRT_vectors` | ✅ 100% | ✅ | CRT for vectors element-wise |
| `half_gcd` | ✅ 100% | ✅ | Fast transformation matrix |
| `continued_fraction` | ✅ 100% | ✅ | Continued fraction expansion |
| `continued_fraction_value` | ✅ 100% | ✅ | Evaluate continued fraction |
| `convergents` | ✅ 100% | ✅ | Convergents of continued fraction |

### `sage.rings.finite_rings` - Finite Fields
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | Reviewed 2026-02-04; re-audited 2026-07-27 — 13 findings, all fixed (286 tests pass in this directory) |
| `finite_field_constructor.ts` | ✅ 100% | ✅ | `GF()` for prime fields; `FiniteField` now aliased to `GFExtended` so both names accept prime powers |
| `integer_mod.ts` | ✅ 100% | ✅ new | Z/nZ elements with Mod(). `log(b, order)` rewritten to Sage's CRT algorithm; test file created (module previously had **none**) |
| `integer_mod_ring.ts` | ✅ 100% | ✅ new | Zmod() ring constructor with iteration; test file created (module previously had **none**) |
| `finite_field_prime.ts` | ✅ 100% | ✅ | GF(p) with sqrt, multiplicative_generator |
| `finite_field_extension.ts` | ✅ 95% | ✅ 12 new | GF(p^n) via polynomial quotient rings. `irreducible_element` is now a faithful port of `polynomial_ring.py:3560-3626` **delegating to ntl-ts `GF2X_BuildSparseIrred` and parigp-ts `ffinit`**; 45/45 explicit-algorithm and 0/130 genuinely-mismatching default moduli against a SageMath 10.3 oracle. Remaining gaps: `algorithm='ffprimroot'` throws (needs PARI `ffgen`/`ffprimroot`/`charpoly`), `algorithm='random'` uses rejection sampling rather than NTL `BuildRandomIrred`, and `minimalPolynomial()` is still the simplified constant-coefficient version |
| `gf2.ts` | ✅ 100% | ✅ new | Optimized GF(2) singleton; test file created (module previously had **none**) |
| `conway_polynomials.ts` | ✅ 100% | ✅ new | Regenerated by porting FLINT's `conway.c` decoder against the vendored bit-packed table. 7 entries were **reducible** and one GF(2^128) entry was fabricated; every entry is now verified irreducible, primitive, normalised and subfield-compatible. p = 2 to n=64, 3 to 24, 5 to 18, 7 to 14, 11/13 to 12, 17/19/23/29/31 to 10 |
| `tower_field.ts` | ✅ 100% | ✅ | Binary tower fields (Binius): Ti(i) = GF(2^(2^i)) |
| `roots_of_unity.ts` | ✅ 100% | ✅ | FFTDomain, CosetDomain, primitive roots |

**Extension Field Features:**
- Frobenius automorphism
- Trace and norm functions
- Minimal polynomial computation
- Integer representation round-trip
- Primitive element finding

### `sage.rings.polynomial` - Polynomials
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | Reviewed 2026-02-04; re-audited 2026-07-27 — 37 findings across univariate + multivariate, all fixed (524 tests in `rings/polynomial`, up from 423) |
| `polynomial_ring.ts` | ✅ 100% | ✅ new | PolynomialRing with lagrange, vanishing, cyclotomic; Neville/divided-difference corrected to Sage's row+table semantics |
| `polynomial_element.ts` | ✅ 95% | ✅ | Full arithmetic + `factor()`, `roots()`, `is_irreducible()`. **Fixed in 0.0.12:** factoring over QQ works (`QQ.__call__` rewritten), and the underlying integer factorization is now a **real Zassenhaus** (DDF + Cantor-Zassenhaus, multifactor Hensel lifting, exact Landau-Mignotte bound, subset recombination) — the old code only peeled off integer roots for degree <= 10 and returned the rest as one "irreducible" factor. Verified against an independent Kronecker oracle on 750 polynomials. Remaining: no van Hoeij/LLL recombination (raises after a 200 000-subset budget) |
| `polynomial_ring_constructor.ts` | ✅ 100% | ✅ | PolynomialRingConstructor() returning [R, x] |
| `quotient_ring.ts` | ✅ 100% | ✅ new | R[x]/<f(x)> for field extensions |
| `convolution.ts` | ✅ 100% | ✅ | Full port of `convolution.py`'s ring-agnostic Schönhage algorithm, plus the FFT/NTT fast path |
| `multi_polynomial_ring.ts` | ✅ 95% | ✅ | Multivariate polynomials R[x,y,z,...]. `__call__` from a univariate polynomial throws (needs `_mpoly_dict_recursive`) |
| `multi_polynomial_element.ts` | 🟡 60% | ✅ | Sparse multivariate + sumcheck/GKR methods. ~19 of Sage's `MPolynomial` methods are now honest `SAGE_NOT_IMPLEMENTED` stubs rather than absent; only 3 of Sage's 12 term orders |
| `multi_polynomial_ideal.ts` | 🟡 70% | ✅ new | Buchberger only. Raises rather than truncating or hanging; `dimension()` now the real Cox-Little-O'Shea algorithm. `multi_polynomial_ideal.test.ts` created, pinning Sage's Katsura-3 lex basis |
| `polynomial_gf2x.ts` | ✅ 100% | ✅ 52 tests | GF(2)[x] bit-packed. **Now delegates to ntl-ts** for the whole arithmetic and irreducibility layer (`buildSparseIrred` uses NTL's real minimal-weight table); old-vs-new values identical for every n ∈ [2,160]. The four factoring routines stay local because ntl-ts's still throw |

> `polynomial_commitment.ts` **moved to `src/zk/` in 0.0.12** — it has no SageMath counterpart, so
> it no longer sits inside the mirrored Sage tree. See the ZK row in Cross-Cutting Infrastructure.
> `rings/polynomial/index.ts` keeps a backwards-compatible re-export block, so the `./rings` and
> `./rings/polynomial` package subpaths are unchanged.

**Polynomial Features:**
- Generic polynomials over any CoefficientRing
- Polynomial division (quo_rem) over fields
- Polynomial GCD via extended Euclidean algorithm
- **Factorization:** factor(), roots(), is_irreducible(), squarefree_decomposition()
- **Interpolation:** lagrange_polynomial(), newton_interpolation(), barycentric
- **FFT/NTT:** O(n log n) multiplication, domain evaluation/interpolation
- **Multivariate:** MPolynomialRing, term orders (lex, deglex, degrevlex)
- **Sumcheck/GKR:** degrees(), variables(), subs(), args(), partialEvaluate()
- **ZK Helpers:** vanishing_polynomial(), cyclotomic_polynomial(), FRI fold

### `sage.matrix` - Matrices
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | 2026-07 (`AUDIT-2026-07.md`): 57 findings across the matrix modules, all fixed. **749 tests pass in `matrix/`** (was 691) |
| `matrix_space.ts` | ✅ 90% | ✅ new 17 | MatrixSpace, Matrix class, scalar-matrix `__call__`. `matrix_space.test.ts` created |
| `matrix_generic.ts` | ✅ 90% | ✅ new 14 | Constructors, arithmetic, `pow` with negative exponents. `matrix_generic.test.ts` created |
| `matrix_operations.ts` | ✅ 93% | ✅ 196 | `minpoly` (works over QQ now that `factor()` does; Sage's `x^3-30x^2-80x` doctest reproduced), `is_semisimple`, exact `is_positive_(semi)definite`, `right_kernel_matrix`, `density`, `eigenvalues`. **New in 0.0.12:** `norm(A,2)` (exact isolation of the largest eigenvalue of `A^H A` — matches three Sage doctests bit-for-bit; raises for inexact base rings), `is_similar(transformation=true)` (verified against exhaustive brute force over 134 785 matrix pairs), a generic `change_ring`, and `is_diagonalizable(base_field)` |
| `matrix_integer.ts` | ✅ 97% | ✅ | HNF, SNF, elementary divisors, kernel, exact integral LLL (`delta` 0.99), symplectic form. **New in 0.0.12:** `frobenius_form` flags 0/1/**2** are a verbatim port of PARI's `RgM_Frobenius` (`alglin2.c:428-720`); `flag=2` returns `[F, B]` over QQ with `B^-1 F B == A`, verified against SageMath on 300 random matrices. ⚠️ flags 0/1 **changed block order** to PARI's (minimal polynomial first) |
| `matrix_modn.ts` | ✅ 95% | ✅ | charpoly, determinant, echelonize, `right_kernel_matrix` (all three basis formats). **Composite modulus now delegates to parigp-ts `matkermod`** — 300 random matrices match SageMath entry for entry, 2190 brute-force cases confirm the returned rows generate the *full* kernel |
| `matrix_mod2.ts` | ✅ 95% | ✅ | GF(2) matrices; `pluq`/`ple` now use M4RI's transposition-list convention for both P and Q |
| `matrix_decompositions.ts` | ✅ 97% | ✅ 121 | RREF echelon form, LU, QR, Cholesky, Bunch-Kaufman `block_ldlt`, Smith, Hermite, LLL_gram. **New in 0.0.12:** `jordan_form(transformation=true)` (reproduces Sage's *exact* `P`, not merely a valid one) and `krylov_kernel_basis(variable=…)` (every matrix in Sage's docstring reproduced verbatim). `pivots` is now re-exported from `matrix/index.ts`. Still ignored: `jordan_form`'s `subdivide` (the port's `Matrix` has no subdivision support) |
| `matrix_decompositions_additions.ts` | ✅ 100% | ✅ | SVD_double, QR_double, LU_double for IEEE 754 real matrices |
| `matrix_special.ts` | ✅ 95% | ✅ 98 | `companion_matrix`, `toeplitz`, `hankel`, `elementary_matrix`, `block_matrix`, `rook_vector`, `berlekamp_massey`, `is_permutation_of`, `permutation_normal_form`, random matrix constructors. **Its divergences are now consolidated in DEVIATIONS.md** (audit item L44); five constructors that need `sqrt`/trigonometry over an inexact ring throw |

**Matrix Integer Features:**
- `hermite_normal_form()` - Hermite Normal Form for integer matrices
- `smith_form_integer()` - Smith Normal Form (D, U, V) where D = U·A·V
- `elementary_divisors_integer()` - Diagonal of SNF
- `rank_integer()` - Rank via HNF
- `kernel_matrix()` / `left_kernel_matrix()` - Null space computation

---

## Phase 2: Elliptic Curves

### `sage.schemes.elliptic_curves`
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | Reviewed 2026-02-04; re-audited 2026-07-27 — 48 findings, 46 fixed. 585 pass / 13 skip / 0 fail in this directory (was 520) |
| `constructor.ts` | ✅ 100% | ✅ new 8 | EllipticCurve() - delegates to parigp-ts ellinit. `constructor.test.ts` created (module previously had **none**) |
| `ell_generic.ts` | ✅ 95% | ✅ | Invariants, torsion_points, is_on_curve, `multiplication_by_m`, `_isomorphisms` (all char 2/3/p branches), `montgomery_model`, `lift_x`/`is_x_coord`, bivariate `division_polynomial` |
| `ell_finite_field.ts` | ✅ 96% | ✅ | cardinality, trace, generators, twists, torsion_basis, `abelian_group`, `set_order`/`has_order`, `frobenius_order`, Vélu. `is_j_supersingular` skips Sage's precomputed j-polynomial table (exact anyway under the default `proof=True`) |
| `ell_point.ts` | ✅ 96% | ✅ | Point arithmetic, weil/tate/ate pairings, `division_points` (Sage's `_multiple_x_numerator` algorithm), `is_divisible_by`, `point_log`. No p-adic shortcut for anomalous curves |
| `ell_curve_isogeny.ts` | 🟡 88% | ✅ ~35 new | Vélu (y-coordinate sign error fixed), real Kohel implementation, real fastElkies' BMSS, `compute_intermediate_curves`, `dual()`. Stark's algorithm throws (needs `ell_wp.py`); odd-degree kernel-polynomial validation skipped (needs `isogeny_small_degree.py`) |
| `formal_group.ts` | ✅ 95% | ✅ 41 | `differential()`, `log()`, `inverse()`, `group_law()`, `mult_by_n()`, `sigma()`. **0.0.12 audit finding: these were already correct at 0.0.11 — the deferral note calling `differential()` hardcoded was stale.** This pass proved it by running Sage's own doctests verbatim (incl. the 35-term `w(35)`, the 16-term `mult_by_n(100,20)` for 37a, and the `# long time` GF(17) `mult_by_n(10,50)`) and by checking defining identities: the Weierstrass equation for `x,y` with a negative control, `log`/`exp` mutual inversion, `log(F(t1,t2)) = log t1 + log t2`, `log([n]t) = n·log t`, `[m+n] = F([m],[n])`, `F(t, i(t)) = 0`. Added `x_list`/`y_list` because `x()`/`y()` returned objects whose coefficients no caller could read. Only gap: Sage's characteristic-zero shortcut inside `mult_by_n` (needs Laurent arithmetic + curve base-change) — **speed only**, outputs verified identical |
| `ell_torsion.ts` | ✅ 90% | ✅ | `_p_primary_torsion_basis` replaced with Sage's division-polynomial algorithm. Torsion over number fields still throws |
| `weierstrass_morphism.ts` | ✅ 100% | ✅ | `order()`, automorphism enumeration, isomorphism composition |
| `isogeny_class.ts` | 🟡 80% | ✅ | IsogenyClass, IsogenyClassNumberField, IsogenyClassRational, `fill_isogeny_matrix`. **0.0.12:** `Frobenius_filter` (`gal_reps_number_field.py:492-586`) ported and applied over QQ — both of Sage's doctests reproduce. Two real defects fixed in `isogeny_degrees_cm`: the horizontal (class-group) primes were **missing entirely**, which made the list potentially *too small* (a soundness gap), and the downward-ramified test used `n` instead of `n/(2h)`. Sage's `d = -23` verbose trace now reproduces line for line. **Still deferred:** the filter is inapplicable over number fields (needs `K.primes_above`/`E.reduction(P)`), so `isogeny_degrees_cm` returns the unfiltered superset there; `possible_isogeny_degrees` still uses Mazur's list |
| `cm.ts` | ✅ 95% | ✅ new 18 | `cm_j_invariants`, `cm_orders`, `is_cm_j_invariant`, `discriminants_with_bounded_class_number` (6.26 s -> 47 ms), `largest_(fundamental_)disc_with_class_number`. All four return shapes corrected to Sage's. `cm.test.ts` created (module previously had **none** despite being listed as delivered) |
| `padic_lseries.ts` | 🔴 40% | ✅ 82 tests | **Still blocked on modular symbols** (`sage.modular.modsym`), which gates `series`, `measure`, `modular_symbol`, `order_of_vanishing`, `_c_bound` and the three `Dp_valued_*` methods. **New in 0.0.12:** `bernardi_sigma_function` (Sage's 14a doctest reproduced verbatim, plus an independent Weierstrass-℘ recursion check) and `alpha` at a **supersingular** prime, via new `pAdicEisensteinQuadraticExtension`/`Element` classes implementing the ramified quadratic extension `K[x]/(x^2 - a_p x + p)` |

**Elliptic Curve Features:**
- Curve initialization (ellinit, ellfromj) via parigp-ts
- Point operations (elladd, ellsub, ellmul, ellneg)
- Point constructor accepts both `E.point(x, y)` and `E.point([x, y])` for Sage parity
- Group operations (ellcard, ellgroup, ellorder, ellgenerators)
- **Pairings:** weil_pairing(), tate_pairing(), ate_pairing() via Miller's algorithm
- **elllog** (Pohlig-Hellman ECDLP) - discrete log on elliptic curves
- **Torsion:** torsion_basis(), torsion_subgroup(), division_points()
- **Twists:** quadratic_twist(), twists(), curves_with_j_0(), curves_with_j_1728()
- **Isogeny graph:** j_invariant_neighbors(), isogenies_prime_degree(), isogeny_class()
- **Isogeny class:** IsogenyClass with matrix(), get(), index(), reorder(), copy(), contains()
- **p-adic L-series:** Basic structure (elliptic_curve, prime), subclasses (Ordinary, Supersingular)
- **Formal group:** EllipticCurveFormalGroup with w(), x(), y(), differential(), log(), inverse(), group_law(), mult_by_n(), sigma()
- **CM functions:** hilbert_class_polynomial, cm_j_invariants, cm_orders, is_cm_j_invariant, discriminants_with_bounded_class_number
- **Point discrete log:** point_log() using baby-step giant-step
- **Point divisibility:** is_divisible_by() for checking if m|P
- **Hyperelliptic polynomials:** hyperelliptic_polynomials() returns (g(x), h(x)) for y^2 + h(x)y = g(x)
- **Weierstrass morphisms:** baseWI (u,r,s,t) transformations, WeierstrassIsomorphism class, _isomorphisms generator, identity_morphism, negation_morphism

### `sage.schemes.elliptic_curves.weierstrass_morphism`
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| `weierstrass_morphism.ts` | ✅ 100% | ✅ 18 tests | baseWI, WeierstrassIsomorphism, _isomorphisms, identity_morphism, negation_morphism |

---

## Phase 3: Additional Crypto Primitives

### `sage.crypto`
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | 2026-07 (`AUDIT-2026-07.md`): 15 findings, all fixed. 335 crypto tests pass |
| `lattice.ts` | ✅ 96% | ✅ | `gen_lattice` with all four types and the `ntl`/`lattice` output flags. **0.0.12:** the `modular`/`random`/`dual` doctests now reproduce **exactly** — Sage draws that block with `rstate.c_random() % q` row-major (`matrix_modn_dense_template.pxi:2843`), not `mpz_urandomm`. `ideal`/`cyclotomic` still differ: they go through `PolynomialRing.random_element`, whose draw order is not ported, and the expected value is Sage-version dependent |
| `lwe.ts` | ✅ 98% | ✅ | LWE, Regev, LindnerPeikert, RingLWE, RingLindnerPeikert, RingLWEConverter. Sage's doctest parameters now reproduce exactly |
| `boolean_function.ts` | ✅ 100% | ✅ | Walsh transform, nonlinearity, ANF, correlation immunity, algebraic immunity, hex truth tables |
| `sbox.ts` | ✅ 100% | ✅ | DDT, LAT (per-mask Walsh-Hadamard; AES 175 ms -> 6.5 ms), APN detection, `min_degree`, MISTY/Feistel constructions |

### `sage.modules` - Lattices and Modules
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | 2026-07 (`AUDIT-2026-07.md`): 27 findings, all fixed. 255 tests pass in `modules/` |
| `free_module.ts` | ✅ 95% | ✅ 282 | FreeModule hierarchy matching Sage's, exact rank, echelonized bases, coordinates, kernels, saturation, discriminant, cardinality. **New in 0.0.12:** free modules over **non-ZZ PIDs** (`QQ[x]`, `GF(p)[x]`) via an exact `QQ(x)` fraction-field layer plus ports of `_echelon_form_PID`, `_generic_clear_column` and Sage's generic Smith normal form — 700 random spans match the vendored SageMath exactly; real `quotient` lift/project (Sage's `__quotient_matrices` over a field, `FGP_Module` over ZZ, both agreeing with Sage on 391 random cases); and a real embedded `tensor_product` (the previous `@see Reference` citation was **fabricated** — Sage's `free_module.py` has no such method). Remaining: `intersection()` over a non-ZZ PID returns the correct *module* but its echelon basis can differ from Sage's by a unit constant |
| `free_module_element.ts` | ✅ 95% | ✅ new 23 | Inner-product-matrix pairing, 7-D cross product, true p-norms, exact normalization, Python-style indexing. `free_module_element.test.ts` created |
| `free_module_integer.ts` | ✅ 90% | ✅ new 41 | Exact integral LLL, exact CVP/SVP enumeration, Voronoi relevant vectors and cell, q-ary lattices. `free_module_integer.test.ts` created. Note: this and `matrix_integer.LLL` are two independent LLL implementations where SageMath has one |
| `bkz.ts` | ✅ 100% | ✅ 34 | BKZ reduction, HKZ, Schnorr-Euchner enumeration (no audit findings) |

### `sage.stats.distributions` - Sampling Distributions
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| `discrete_gaussian_integer.ts` | ✅ 95% | ✅ | **All four** algorithms now, including `uniform+logtable` and `sigma2+logtable`: `dgs_bern.c` and `dgs_disc_gauss_sigma2p` are ported, and the acceptance test is dgs's real `mpfr_urandomb` comparison (0.0.11 used a different test, so even the two implemented algorithms consumed randomness differently). All four reproduce SageMath's **seeded sample streams** bit-for-bit, and `repr()` matches Sage's format. **Not implemented:** the `precision='dp'` mode (Sage documents its results as not reproducible) |
| `discrete_gaussian_lattice.ts` | ✅ 90% | ✅ | GPV algorithm, coset sampling, exact (`Rational` basis/centre/GSO) with `sampleExact()`. **New in 0.0.12:** non-spherical Σ (matrix sigma, `sigma_basis`, Peikert's `r` by power iteration, Cholesky, offline samples, `_call_non_spherical`), `set_c`/`c()`/`sigma()`/`f()`, and `_normalisation_factor_zz` with a local `qfrep`. Every doctest in the vendored source reproduces. **Remaining:** `_normalisation_factor_zz` runs in double precision, so Sage's `prec=100` doctest matches only to 15 significant digits; `qfrep` belongs in parigp-ts; the local LLL is not fpLLL |

### `sage.coding` - Error-Correcting Codes
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | 2026-07 (`AUDIT-2026-07.md`): 18 findings, all fixed. 211 tests (was 173) |
| `reed_solomon.ts` | ✅ 100% | ✅ | RS codes with encode/decode; syndrome/Forney reproduce Sage's `GRSKeyEquationSyndromeDecoder` doctests. FRI fold/query are port-only additions |
| `bch_code.ts` | ✅ 95% | ✅ | BCH codes, PGZ decoding, Chien search, real field embedding with section. `minimum_distance()` enumerates exactly and raises above `q^k > 2^17` (it previously returned the *designed* distance) |
| `goppa_code.ts` | ✅ 95% | ✅ | Goppa codes, Patterson algorithm (McEliece-ready) plus a key-equation decoder for the non-binary case. `distance_bound()` is now Sage's `1 + deg(g)` |
| `reed_muller_code.ts` | ✅ 100% | ✅ | RM(r,m) codes, Plotkin construction, majority decoding (now decodes `u` from both halves); monomial order matches Sage's `Subsets` enumeration |

**Coding Features (for ZK and Post-Quantum):**
- Reed-Solomon encoding/decoding with FRI operations
- BCH codes for classical error correction
- Goppa codes for McEliece post-quantum cryptosystem
- Reed-Muller codes for recursive ZK constructions
- Error correction via Gao, PGZ, Patterson algorithms

---

## Cross-Cutting Infrastructure

| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| `misc/randstate.ts` | ✅ 100% | ✅ | Centralized RNG + `set_random_seed` parity. **0.0.12: seeded streams now match SageMath exactly.** GMP 6.3.0 is not vendored under `reference/`, so the sources were obtained and `rand/randmts.c` (`mangle_seed`, `randseed_mt`), `rand/randmt.c` (incl. the 624-word `default_state` that seed 0 lands on) and `mpz/urandomm.c` were ported verbatim. CPython's `random.Random` is ported too, as `PythonRandom`, since `randstate.python_random()` is Sage's second generator. Verified against a C oracle linked to libgmp 6.3.0 and against SageMath 10.3's own doctest values |
| `quadratic_forms/binary_qf.ts` | ✅ 97% | ✅ 99 tests | BinaryQF. **0.0.12: composition and reduction now delegate to parigp-ts `qfb.ts`** (~170 lines of transcribed `Qfb.c` deleted), and a `solve_integer` the port did not have was added. Equivalence with the pre-delegation code proven against a side-by-side HEAD import: 29 944 forms, 59 280 compositions and 500 class-group Cayley tables with 0 differences. A fidelity bug was fixed en route — every `D > 0` form went through Sage's `_reduce_indef`, where Sage uses it only for **square** discriminants |
| `groups/generic.ts` | ✅ 95% | ✅ 106 tests | Sage's `discrete_log` loop verbatim (incl. repair of a non-minimal `ord` and the `<30` linear branch of `bsgs`); `order_from_multiple` honours `check=True` |
| `zk/sumcheck.ts`, `zk/multilinear.ts`, `zk/polynomial_commitment.ts` | ✅ 90% | ✅ 172 tests | Ports of `reference/sage_blueprints/`, **not** of SageMath. `polynomial_commitment.ts` (KZG/FRI helpers) **moved here from `rings/polynomial/` in 0.0.12** because it has no SageMath counterpart; `rings/polynomial/index.ts` keeps a compatibility re-export. ⚠️ `package.json` has no `./zk` subpath export yet, so the symbols are reachable from the package root and via `./rings` but not as `@sagemath-ts/sagemath-ts/zk` |

---

## Dependency Libraries Progress

### parigp-ts (port of cypari2 / PARI/GP)

**Total: 655 tests passing, 0 failing** (was 414; +241 in the 0.0.12 deferred-work pass)

| Feature | Status | Tests | Notes |
|---------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | Reviewed 2026-02-04; re-audited 2026-07-27 — 11 findings, all fixed. The 2026-02 audit had marked `ellcard`/`ellgroup` fidelity-checked; they were in fact returning **wrong values** at primes as small as p ≈ 100 |
| Core types (GEN, t_INT, etc.) | ✅ | - | types.ts - PariType enum, PariInt, PariFfelt, PariVec, etc. |
| Fp arithmetic | ✅ | 55 | ff.ts - Fp_add, Fp_sub, Fp_mul, Fp_sqr, Fp_neg, Fp_inv, Fp_div, Fp_pow |
| Fp_sqrt (Tonelli-Shanks) | ✅ | - | ff.ts - includes Fp_issquare, kronecker symbol; returns the canonical smallest root as PARI does |
| ellinit | ✅ | 35 | elliptic/init.ts - Short/general Weierstrass, from j-invariant. `j`/`ellj` return an exact `Ratio` when non-integral |
| Point operations (elladd, ellmul) | ✅ | 65 | elliptic/point.ts - Jacobian coordinates for efficiency |
| ellcard, ellgroup, ellorder | ✅ | 36+ | elliptic/group.ts - faithful `Fp_ellcard_Shanks` and `gen_ellgroup`/`gen_ellgens` with the real Weil pairing. Verified against exhaustive point-count oracles (4542 curves), a brute-force group-structure oracle (476 runs, 0 wrong — was 85 wrong), and real PARI 2.15.4 via Sage 10.3 for primes up to 2^32 |
| ellordinate, random_FpE | ✅ | 32 | elliptic/points.ts - Find y from x, random point generation |
| ellgenerators, trace_of_frobenius | ✅ | - | elliptic/group.ts - 10 068-curve sweep |
| `Fp_ellcard_CM` (full CM table) | ✅ 100% | ✅ new | **New in 0.0.12.** All **thirteen** class-number-one discriminants (`FpE.c:624-666`, `:1282-1421`), delegating to `qfb.ts`'s `cornacchia2`. Verified against brute-force point counting (15 392 curves), Shanks (3744), a counting-independent `[#E]P = O` oracle (936 at 64/80/96 bits) and the **published SECG group orders** of secp160k1/192k1/224k1/256k1 |
| `Fp_ellcard_Schoof` / `ellcard_sea` | 🟡 60% | ✅ new | **New in 0.0.12.** Schoof's base algorithm — the "S" of SEA — so `ellcard_sea` no longer throws and returns exact cardinalities at every size. Verified exhaustively on all 121 104 curves over every prime `5 <= p <= 120` and against **PARI's own `ellsea` regression vectors** at 65, 70 and 101 bits. **Elkies and Atkin are NOT ported** (see Not Yet Implemented) |
| `Z_factor` (real factoring chain) | ✅ 90% | ✅ 49 | **Rewritten in 0.0.12** (746 -> 1609 lines): `tridiv_bound` + gcd-with-primorial trial division, then PARI's `ifac_crack` order — pure powers, SQUFOF, Pollard-Brent rho, ECM — driven by an `ifac_decomp` worklist. Plus `ispower.c`'s `Z_issquareall`, `is_357_power`, `is_kth_power`, `is_pth_power`, `Z_isanypower` and a real `isprimepower` that never factors `n`. **MPQS is the one missing stage** |
| `ffinit` (Adleman-Lenstra) | ✅ 100% | ✅ 24 | **New in 0.0.12.** `polarit3.c`'s full chain plus the supporting `FpX` layer, `FpX_composedsum`, the bivariate resultant and `polsubcyclo` for prime conductor. Reproduces PARI **coefficient for coefficient** for all of the first 60 primes × n ∈ [2,12] (660/660), each independently re-verified irreducible |
| `matkermod` / `matimagemod` / `matdetmod` / `matinvmod` | ✅ 100% | ✅ 41 | **New in 0.0.12.** `bb_hnf.c` specialised to the `Z/dZ` Hermite ring (Howell form, `gen_kernel`, `gen_matimage`, `gen_inv`, `gen_detops`). All 24 golden values decoded from PARI's own regression suite reproduced verbatim; kernels confirmed complete by exhaustive enumeration |
| `Qfb` family | 🟡 85% | ✅ 33 | **New in 0.0.12.** `qfbred`, `qfbredsl2`, `qfbcomp(raw)`, `qfbsqr(raw)`, `qfbpow(raw)`, `qfbsolve` (all 4 flags), `primeform`, `cornacchia`/`cornacchia2`, `Zp_sqrt`/`Z2_sqrt`/`Zn_quad_roots`, and the Schoenhage fast reduction. Verified against real PARI 2.15.4 on ~2500 golden values. **Not ported:** the Shanks-distance (`qfr5_*`, `t_REAL`) variants |

**Implemented PARI Functions:**
- **Types:** mkInt, stoi, itos, mkFfeltFp, mkvec, mkcol, mkmat, gen_0/gen_1/gen_2/gen_m1
- **Fp operations:** Fp_red, Fp_add, Fp_sub, Fp_neg, Fp_mul, Fp_sqr, Fp_inv, Fp_div, Fp_pow, Fp_sqrt, Fp_issquare, Fp_center, Fp_halve, Fp_double, Fp_addmul, kronecker
- **Elliptic curves (init):** ellinit, ellfromj, ellfromjFp, ellj, elldisc, ellisnonsingular, ellcoeffs, ellToShortWeierstrass
- **Elliptic curves (points):** ellinf, ell_is_inf, mkpoint, ellordinate, random_FpE, FpE_to_FpJ, FpJ_to_FpE, FpE_isoncurve
- **Elliptic curves (point ops):** elladd, ellsub, ellneg, ellmul, FpJ_add, FpJ_dbl, FpJ_neg (Jacobian arithmetic)
- **Elliptic curves (group):** ellcard, ellgroup, ellorder, ellgenerators, trace_of_frobenius, FpE_random, ellinit_Fp, ellisoncurve, elllift_x

**Advanced Features (elliptic/advanced.ts):**
- elllog (Pohlig-Hellman discrete logarithm)
- Weil pairing (Miller's algorithm)
- Tate pairing (Miller's algorithm)
- Division polynomials (psi_n)

**Not Yet Implemented:**
- **Elkies and Atkin** (`ellsea.c`'s `find_trace`, `find_trace_Elkies_powerell`,
  `find_trace_Atkin`, `match_and_sort`, `champion`). Both need the modular polynomials `Phi_l`,
  which PARI reads from the separately distributed **`seadata` package** — `reference/pari` ships
  the reader (`ellsea.c:47-101`) but `reference/pari/data` is **empty**. Base Schoof is ported in
  their place and is exact; only the complexity differs (`O(log^5 p)` vs `O(log^4 p)`)
- **MPQS** (`mpqs.c`, ~2600 lines) — the one missing stage of `Z_factor`'s chain
- `ffgen` / `ffprimroot` / `charpoly` over `F_q` — blocks `irreducible_element(algorithm='ffprimroot')`
- The Shanks-distance (`qfr5_*`) variants of the `Qfb` family — need an arbitrary-precision
  float kernel that CLAUDE.md forbids
- `nf` module (nfbasis, nfdisc, idealprimedec, nfgaloisconj, quadunit) — currently ported inside sagemath-ts
- `matfrobenius` — currently ported inside sagemath-ts's `matrix_integer.ts`
- `qfbclassno` / `quadclassunit`, `qfrep` — currently ported inside sagemath-ts
- Transcendental functions (`dilog`, `incgam`) — currently ported inside sagemath-ts
- Extension field elliptic curves (t_ELL_Fq)
- Polynomial operations (t_POL)
- ⚠️ **Barrel gap:** `src/index.ts` does not re-export `isprimepower`, `Z_isanypower`, `Z_iroot`,
  `squfof`, `pollardbrent`, `ellfacteur`, `forprime`, `FactorOptions`, `Fp_ellcard_CM`,
  `Fp_ellj_get_CM`, `ec_ap_cm`, `Fp_ellcard_Schoof`, `Fp_elldivpol`, `Ratio` or `isRatio`.
  Nothing is broken (tests import by path), but `arith/misc.ts` still keeps a duplicate local
  `isprimepower` because of it. Adding names to that barrel is exactly how a runtime-only
  `export`-of-a-type bug shipped in `flint-ts`, so this was deliberately not done speculatively

### flint-ts (port of FLINT)
| Feature | Status | Tests | Notes |
|---------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | 2026-07: **all APIs are still stubs** (FLINT_NOT_IMPLEMENTED). Consequence: sagemath-ts reimplements FLINT primitives rather than delegating, which the audit identifies as a recurring structural cause of defects |
| Package barrel (`src/index.ts`) | ✅ | ✅ new | **0.0.12 bug fix:** the barrel re-exported five *interfaces* (`nmod_t`, `fmpz_factor`, …) through a **value** `export` clause. `tsc --noEmit` passed (it elides them) but at runtime `import … from '@sagemath-ts/flint-ts'` threw `export 'nmod_t' not found`, i.e. **the package could not be imported at all** — pre-existing since the first commit. Fixed with `export type`. `index.test.ts` created; flint-ts had **zero** tests before, which is why it survived |
| fmpz (integers) | ⬜ | - | |
| fmpz_poly | ⬜ | - | |
| fmpz_mod_poly | ⬜ | - | |
| nmod_poly | ⬜ | - | |

### ntl-ts (port of NTL)
| Feature | Status | Tests | Notes |
|---------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | 2026-07: `GF2`/`GF2X` implemented; the rest are still stubs (NTL_NOT_IMPLEMENTED). 15 tests / 2645 assertions — the package had **zero** tests before |
| ZZ | ⬜ | - | |
| ZZ_p | ⬜ | - | |
| ZZ_pX | ⬜ | - | |
| GF2 | ✅ 100% | ✅ | |
| GF2X | ✅ 85% | ✅ 15 | Bit-packed bigint representation. `IterIrredTest`, `BuildIrred`, `BuildSparseIrred` ported line-for-line from `GF2XFactoring.cpp` over a vendored copy of NTL's 2049-row `GF2X_irred_tab`; verified against Sage's `polynomial_gf2x.pyx` doctests and exhaustive brute force for degree <= 10. **Consumers wired up in 0.0.12:** both `finite_field_extension.irreducible_element` and `polynomial_gf2x.ts` now delegate here. **Not implemented:** `random`, `BuildRandomIrred`, `factor`, `SquareFreeDecomp`, `DistinctDegFactor`, `EqualDegFactor`, `BerlekampFactor` (need NTL's ChaCha `RandomStream`, `IrredPolyMod`/`GF2XModulus`, or randomized factoring) — which is why `polynomial_gf2x.ts` keeps four local factoring routines |

---

## Future Cryptography Modules (Not Yet Implemented)

This section tracks cryptography-relevant SageMath functionality that could be added.

### `sage.groups.generic` - Generic Group Operations
Generic algorithms for any group with compatible operations.

| Function | Status | Priority | Notes |
|----------|--------|----------|-------|
| Audit (algorithm fidelity) | ✅ | - | Reviewed 2026-02-04; re-audited 2026-07-27 — 7 findings, all fixed |
| `discrete_log` | ✅ 90% | HIGH | Sage's Pohlig-Hellman + BSGS loop verbatim; no bounds/algorithm options |
| `discrete_log_rho` | ✅ 90% | HIGH | Pollard's rho; requires explicit prime order |
| `discrete_log_lambda` | ✅ 100% | MEDIUM | Pollard's kangaroo (bounded DLP) |
| `bsgs` | ✅ 100% | HIGH | Baby-step giant-step algorithm, incl. Sage's `<30` linear branch |
| `pohlig_hellman` | ✅ 100% | HIGH | Reduce DLP to prime power subgroups |
| `order_from_multiple` | ✅ 95% | HIGH | Compute element order given a multiple; `check`/`plist` now honoured (in a trailing options object, not Sage's argument positions) |
| `multiple_of_order` | ✅ 100% | HIGH | Find multiple of element order |
| `has_order` | ✅ 95% | HIGH | Check if element has given order (integer input only) |
| `multiple` | ✅ 100% | HIGH | Generic scalar multiplication |
| `multiples` | ✅ 100% | HIGH | Compute [0·g, 1·g, ..., n·g] |
| `linear_relation` | ⬜ | MEDIUM | Not implemented |
| `merge_points` | ⬜ | LOW | Not implemented |
| `structure_description` | ⬜ | LOW | GAP-dependent |

**Use Cases:** Finite field DLP, elliptic curve DLP, generic cyclic group attacks

### `sage.rings.finite_rings` - Additional Finite Field Features
| Feature | Status | Priority | Notes |
|---------|--------|----------|-------|
| `element.log(base, order?)` | ✅ 100% | HIGH | Discrete log in Z/nZ via Pohlig-Hellman |
| `multiplicative_order()` on GF(p^n) | ⬜ | HIGH | Currently only on GF(p) |
| `GaloisGroup_GF` | ⬜ | MEDIUM | Galois group (cyclic, Frobenius) |
| `FiniteFieldHomset` | ⬜ | LOW | Field homomorphism enumeration |
| `ResidueField` | ⬜ | LOW | Residue fields of DVRs |

### `sage.crypto` - Classical and Symmetric Cryptography
| Module | Status | Priority | Notes |
|--------|--------|----------|-------|
| `SBox` | ✅ 100% | MEDIUM | S-box with differential_uniformity, linearity, DDT, LAT |
| `BooleanFunction` | ✅ 100% | MEDIUM | Walsh transform, nonlinearity, algebraic immunity |
| `RingLWE` | ✅ 100% | HIGH | Ring-LWE, RingLindnerPeikert, RingLWEConverter |
| `DiscreteGaussianSampler` | ✅ 100% | HIGH | Discrete Gaussian over integers/lattices (GPV algorithm) |
| `DiffieHellman` | ⬜ | LOW | DH key exchange class |
| `BlumGoldwasser` | ⬜ | LOW | Probabilistic public-key encryption |

**Classical ciphers (low priority for ZK):**
- Affine, Hill, Shift, Substitution, Transposition, Vigenere
- LFSR, Shrinking Generator (stream ciphers)
- DES, Mini-AES, PRESENT (block ciphers)

### `sage.quadratic_forms` - Quadratic Forms (MEDIUM PRIORITY)
Relevant for lattice-based crypto and class group crypto.

| Module | Status | Priority | Notes |
|--------|--------|----------|-------|
| `BinaryQF` | ✅ 97% | MEDIUM | Binary quadratic forms ax^2 + bxy + cy^2 — now **delegating to parigp-ts `qfb.ts`**, plus a new `solve_integer`. See the Cross-Cutting table above |
| `TernaryQF` | ⬜ | LOW | Ternary quadratic forms |
| `BQFClassGroup` | 🟡 60% | MEDIUM | Class group of binary QFs — reduced representatives form a genuine class group (28 discriminants verified against the literature); no dedicated `BQFClassGroup` class |
| `qfbsolve` | ✅ 95% | MEDIUM | **New in 0.0.12** as `BinaryQF.solve_integer` and parigp-ts `qfbsolve` (all four PARI flags). Inherits `Z_factor`'s MPQS gap for hard `n`, with an optional precomputed-factorisation escape hatch |
| `qfsolve` (Simon's algorithm, rational solutions) | ⬜ | LOW | Solve quadratic equations |
| `least_quadratic_nonresidue` | ⬜ | LOW | Find smallest QNR mod p |

**Use Cases:** Class group cryptography, lattice reduction, quadratic sieve

### `sage.rings.number_field` - Algebraic Number Theory (PARTIAL)
Needed for advanced algebraic crypto constructions.

| Module | Status | Priority | Notes |
|--------|--------|----------|-------|
| Audit (algorithm fidelity) | ✅ | - | 2026-07 (`AUDIT-2026-07.md`): 17 findings, 16 fixed. **343 tests** (was 228) |
| `pari_nf.ts` | ✅ 92% | ✅ | The number-field kernel: `nfbasis`/`nfdisc` (Round 2), `idealprimedec`, `polisirreducible`. **0.0.12:** `nfgaloisconj`'s degree-8 cap is **gone** (replaced by LLL reconstruction with a two-sided proof — exact `g(beta)=0` verification and a Gram-Schmidt non-existence certificate); `idealprimedec` gained Buchmann-Lenstra **round 4** (`base2.c:2248`), so inessential discriminant divisors work; and `quadunit`/`quadunitnorm` (`quad.c:281`) were added. Belongs in parigp-ts — see DEVIATIONS.md |
| `NumberField` | ✅ 85% | LOW | Basic operations, element arithmetic, norm, trace, real field discriminants, integral bases |
| `NumberFieldElement` | ✅ 85% | LOW | Arithmetic, charpoly, minpoly, is_integral, `is_unit` (field semantics) + `is_integral_unit` |
| `CyclotomicField` | ✅ 85% | MEDIUM | Q(zeta_n) - cyclotomic polynomials, degree, exact `automorphisms()` at any degree |
| `QuadraticField` | ✅ 85% | LOW | Q(sqrt(D)) using `x^2 - D` verbatim as Sage does (it previously rewrote the polynomial to the squarefree part), discriminant |
| `RationalPolynomial` | ✅ 90% | - | Supporting polynomial arithmetic over Q |
| `class_group` | 🟡 60% | LOW | **Quadratic** class groups computed exactly (imaginary and real, incl. Sage's C38 x C2 doctest), guarded at \|D\| <= 2e6. **0.0.12:** degree > 2 answers in the one rigorous sub-case — when the Minkowski bound provably admits **no** prime ideal, `h = 1` (gives Sage's answer for `x^3-x-1`, `x^4-x-1`, `Q(zeta_5)`, `Q(zeta_7)`, `Q(zeta_12)`). Everything else, **including fields whose true class number is 1** (`Q(2^(1/3))`, `Q(zeta_8)`, `Q(zeta_23)`), still throws — we will not return an unproved answer. Needs PARI `bnfinit` |
| `unit_group` | 🟡 75% | LOW | **0.0.12:** real quadratic fundamental units implemented (PARI `quadunit_uv_basecase`), so `fundamental_units()`, `units()` and `regulator()` answer; `UnitGroup.log()` for rank 1 is exact (sign test + `O(log k)` doubling/binary search). `regulator()` is overflow- and cancellation-safe (the ~250-digit unit of `Q(sqrt(1000003))` gives a finite 576.646). Degree > 2 still needs `bnfinit` |
| `galois_group` | 🟡 70% | LOW | Built from genuine automorphisms — and **now at any degree**, since `nfgaloisconj`'s cap is gone: `is_galois()` for degree > 2 used to return `false` unconditionally. Still throws for non-Galois fields (Sage returns the Galois closure's group via `galoisinit`), as do `fixed_field`/`frobenius`/`inertia_group` above degree 2 |
| `number_field_ideal` | ✅ 90% | LOW | Ideal arithmetic in HNF: rational norms, `smallest_integer`, `is_prime`, canonical equality, `prime_above`/`primes_above`/`decomposition`. **0.0.12:** `decomposition(p)` branches exactly as PARI does — Dedekind-Kummer when `p` is prime to the index, Buchmann-Lenstra round 4 otherwise. Validated on 2000 decompositions with the exact lattice identity `prod P^e = pO_K` |
| `order` | ✅ 80% | LOW | Maximal order vs equation order, `different`/`codifferent`, `is_maximal`, `index_in_maximal_order`, `conductor` |

**Implemented (Pure TypeScript):**
- Element arithmetic (add, sub, mul, div, pow, inv)
- Norm, trace via characteristic polynomial
- Characteristic polynomial (Faddeev-LeVerrier)
- Minimal polynomial
- Signature via Sturm's theorem
- is_integral, is_unit / is_integral_unit
- QuadraticField with discriminant; CyclotomicField with cyclotomic polynomials
- Maximal orders, integral bases, field discriminants, prime decomposition, ideal HNF arithmetic
- Quadratic class groups; automorphisms; Galois groups of Galois fields

**Not Implemented (Requires PARI bnfinit / galoisinit):**
- Class group for degree > 2, except the provably-trivial Minkowski case
- Archimedean embeddings (`embeddings()`, `real_embeddings()`, `complex_embeddings()`) — the same
  gap as `bnfinit`, which needs LLL reduction of ideal lattices under the `T2` form
- Galois group of a non-Galois field; `fixed_field`/`frobenius`/`inertia_group` above degree 2

*Resolved in 0.0.12:* fundamental units and regulators of real quadratic fields (`quadunit`),
`nfgaloisconj` at any degree (LLL), and `decomposition(p)` at inessential discriminant divisors
(Buchmann-Lenstra round 4).

**Use Cases:** Ring-LWE over cyclotomic fields, class group crypto, NTRU

### `sage.rings.padics` - p-adic Numbers
| Module | Status | Priority | Notes |
|--------|--------|----------|-------|
| Audit (algorithm fidelity) | ✅ | - | 2026-07 (`AUDIT-2026-07.md`): 16 findings, all fixed. 148 tests (was 86 assertions) |
| `Zp(p, prec)` | ✅ 90% | LOW | p-adic integers with capped-relative precision. **`add()` multiplied by p^v twice for operands of unequal valuation** until 0.0.11 |
| `Qp(p, prec)` | ✅ 90% | LOW | p-adic field |
| `pAdicGenericElement` | ✅ 90% | LOW | `nth_root`, `square_root`, `is_square`, `multiplicative_order`, `additive_order`, `expansion`, `slice`, repr all ported from upstream |
| `power_series_ring.ts` | ✅ 90% | LOW | `log` (integral of the logarithmic derivative), `nth_root` (Newton via `_nth_root_series`), `pade` (returns a `PadeApproximant` — no `Frac(R[z])` type here) |
| `Zq(q, prec)` | ⬜ | LOW | Unramified extension of Zp |
| `extension()` | ⬜ | LOW | p-adic field extensions |
| Hensel lifting | ✅ | MEDIUM | Used in sqrt, teichmuller |

**Implemented p-adic Features:**
- Element arithmetic (add, sub, mul, div, neg, pow, inv)
- Valuation, unit_part, normalized_valuation
- Precision handling (absolute, relative, add_bigoh, lift_to_precision)
- Expansion and coefficients (list, expansion, residue, __getitem__, slice)
- Predicates (is_zero, is_one, is_unit, is_integral, is_square)
- Square root (Tonelli-Shanks + Hensel lifting)
- Teichmuller lift (Newton iteration)
- Log (for units, series expansion)
- Exp (for convergent inputs, Horner's method)
- Norm, trace (trivial for base field)
- Multiplicative order, abs, equality

**Not Yet Implemented:**
- Extension fields (Zq, ramified/unramified extensions) — `nth_root`'s p-th-root extraction is
  written only for absolute degree 1, which is exactly SageMath's result there
- minimal_polynomial, charpoly (meaningful only with extension fields)

**Use Cases:** p-adic methods in point counting, Hensel lifting for root finding

### `sage.rings.polynomial` - Advanced Polynomial Features
| Feature | Status | Priority | Notes |
|---------|--------|----------|-------|
| `groebner_basis()` | ✅ 70% | MEDIUM | Buchberger only; no Singular/FGb backends |
| `BooleanPolynomialRing` | ⬜ | MEDIUM | GF(2) polynomials with Grobner |
| `SkewPolynomialRing` | ⬜ | LOW | Ore/skew polynomials F[x;sigma] |
| `WeilPolynomials` | ⬜ | LOW | Iterator over q-Weil polynomials |

**Use Cases:** Algebraic cryptanalysis, multivariate crypto (MQ), code-based crypto

### `sage.modules` - Additional Lattice Features
| Feature | Status | Priority | Notes |
|---------|--------|----------|-------|
| `BKZ` | ✅ 100% | HIGH | Block Korkin-Zolotarev reduction with HKZ |
| `shortest_vector` | ✅ 95% | HIGH | SVP solver: **exact** Fincke-Pohst enumeration for rank ≤ 30, LLL first row above (was a float Schnorr-Euchner walk with a hard coefficient cap of 15) |
| `closest_vector` | ✅ 95% | HIGH | CVP solver: **exact** Fincke-Pohst enumeration seeded with Babai (was centred on the origin with coefficients in [−3,3], so it degraded to Babai for distant targets while claiming exactness); `approximateClosestVector` (nearest_plane, rounding_off, embedding) |
| `voronoi_cell` / `voronoi_relevant_vectors` | ✅ 90% | MEDIUM | Voronoi's L/2L theorem; H-representation instead of a `Polyhedron` (rank ≤ 24) |
| `IntegralLattice` | ⬜ | MEDIUM | Lattices with inner product |
| `TorsionQuadraticModule` | ⬜ | LOW | Finite quadratic modules |

**Use Cases:** Lattice attacks, LWE parameter selection, NTRU cryptanalysis

---

## Reference Repositories

| Repository | Cloned | Path |
|------------|--------|------|
| sagemath/sage | ✅ | `reference/sage/` |
| sagemath/cypari2 | ✅ | `reference/cypari2/` |
| pari/pari | ✅ | `reference/pari/` |
| flintlib/flint | ✅ | `reference/flint/` |
| libntl/ntl | ✅ | `reference/ntl/` |

---

## Notes

### Behavioral Differences Log
`DEVIATIONS.md` at the project root is the single source of truth (**64 sections as of 0.0.12**).
Headlines:
- `GF()` returns `FiniteFieldPrime` for primes; `FiniteField` is aliased to `GFExtended` and accepts prime powers
- Element coercion is explicit via `__call__()` rather than automatic
- Uses TypeScript's `bigint` for arbitrary precision integers
- **The RNG now matches SageMath's seeded streams exactly** (GMP `randseed_mt` + CPython's `Random` both ported); the 0.0.11 caveat about a different seeding step is gone
- Several PARI routines are still ported *in place* rather than delegated (`nf`, `matfrobenius`, `qfbclassno`, `qfrep`, transcendentals), because the target package lacks the module — but `ffinit`, `matkermod`, the `Qfb` family and `GF2X` are now genuinely delegated

### Known Limitations
- Conway database covers p = 2 to n=64, 3 to 24, 5 to 18, 7 to 14, 11/13 to 12, 17/19/23/29/31 to 10. Outside it the default modulus is now **SageMath's own** `ffinit` / `BuildSparseIrred` choice, so it agrees with Sage wherever our Conway coverage matches Sage's — but Sage's table is much larger (it has entries for `37^2`, `97^2`, `2^100`, `19^21`, …), so element representations of those fields are not interoperable and the generator need not be primitive
- **`Z_factor` has no MPQS**, so a hard semiprime with a >= 25-digit smallest factor raises. This is a *source-breaking* change from 0.0.11, where it returned the composite as if prime
- **`ellcard_sea` is Schoof, not SEA** — Elkies and Atkin need the `seadata` modular polynomials, which are not vendored. Exact but `O(log^5 p)`, so `ellcard` keeps Shanks below `expi(p) = 96`
- **Class number / class group of a degree > 2 field** needs `bnfinit`; only the provably-trivial Minkowski case answers, and fields whose true `h` is 1 still throw when we cannot prove it
- **Polynomial factorization has no van Hoeij/LLL recombination** — raises after a 200 000-subset budget rather than returning a partial answer
- **No discrete_log on GF(p^n) extension field elements** - extension field elements lack log() method (Z/nZ has log())
- **BKZ uses Schnorr-Euchner enumeration** - pure TypeScript implementation (no fpylll/NTL bindings)
- **`flint-ts` is 100% stubs and `ntl-ts` is stubs apart from GF2/GF2X**, so sagemath-ts reimplements FLINT/NTL primitives instead of delegating — the audit identifies this as a recurring structural cause of defects
- **`packages/sagemath-ts` has no `tsconfig.json`** (`git ls-files` confirms it never has), so its `typecheck` script has never actually typechecked the package. Against a synthesized strict config it reports ~1670 errors, **all** pre-dating the 0.0.12 pass (mostly the `PolynomialRing<T>` / `PolynomialRingBase<RingElement>` variance issue in `polynomial_element.ts`, plus `vitest` imports). Needs a dedicated task
- **`packages/parigp-ts` reports 29 pre-existing type errors** (27 `EllipticPoint` union narrowing, 2 missing `vitest` types), byte-identical to the 0.0.11 baseline; all pass at runtime under Bun
- **`tests/property/typescript/elliptic.test.ts` is an inert oracle.** `tests/property/transcripts/` is gitignored and empty, so every PARI comparison sits behind `if (pariOutputAvailable && pari)` and never runs. Generating the transcript from a real `gp` 2.15.4 produces 40/74 failures — all *formatting* mismatches (PARI prints `a1=Mod(0, 101)` where `formatCurveInfo` prints `a1=0`) plus test-number drift, not numeric ones. The harness's formatters were never reconciled with `tests/property/pari/elliptic.gp`. Wiring this up is a real project
- **Three tests are timeout-flaky under load** at the default 5000 ms (`ell_point` `point_log`, `discrete_gaussian_lattice` `_normalisation_factor_zz` consistency, and one other). In isolation they take 4.97 s and 3.89 s. Timeouts were **not** raised
- **`types.ts:FieldRing.characteristic` is a property while `rational_field.ts` declares `characteristic()` as a method**, and `constructor.ts:125` reads it as a property. For QQ this falls through to the char-0 branch, which happens to be correct, and every finite field in the repo declares it as a property — so there is no live bug, but it is a type-level trap
- **`packages/parigp-ts/src/matkermod.ts`** has a latent `matDims()`/`transposeMat()` round-trip issue for `n × 0` matrices (it faithfully mirrors PARI's own limitation); worked around from the caller side with `wantIm = true`

### Implementation Summary

**Phase 1 Progress: ~96% complete**
**Phase 2 Progress: ~96% complete**
**Phase 3 Progress: ~95% complete**

> Percentages were revised downward on 2026-07-28 for 0.0.11 (the July 2026 audit found 370
> confirmed defects behind the previous figures) and nudged back up for 0.0.12, where the
> deferred items were actually implemented rather than merely documented. They still mean
> "ported **and** verified against upstream doctests or an executed oracle".

| Category | Implemented | Remaining |
|----------|-------------|-----------|
| sage.arith | 47+ functions (incl. rational_reconstruction, CRT_basis, continued_fraction, trial_division, prime_powers, smooth_part, sum_of_squares) | MPQS inside `Z_factor`; proven primality (APRCL/ECPP) |
| sage.rings.finite_rings | 9 modules (incl. tower_field, roots_of_unity, log in Z/nZ) | element.log() for GF(p^n); real `minimalPolynomial`; a wider Conway table; `ffprimroot` |
| sage.rings.polynomial | 10 modules (incl. FFT, Schönhage convolution, multivariate, ideals, GF2X) | 9 term orders; van Hoeij/LLL recombination; Singular-grade Gröbner |
| sage.matrix | 9 modules (incl. HNF, SNF, mod2, modn, special, decompositions, Jordan + transformation, Frobenius flag 2, `matkermod`) | `rational_form`; `norm(2)` over inexact rings; `jordan_form` subdivisions; symbolic base rings |
| sage.schemes.elliptic_curves | 11 modules (incl. pairings, torsion, twists, isogeny graph, formal group, CM, Bernardi sigma, supersingular alpha) | Stark's algorithm, `ell_wp`, `isogeny_small_degree`, **modular symbols** (gates most of `padic_lseries`), `Frobenius_filter` over number fields |
| sage.crypto | 4 modules (lattice, lwe, boolean_function, sbox) | Seeded parity for `gen_lattice`'s `ideal`/`cyclotomic` branches |
| sage.modules | 4 modules (free_module incl. non-ZZ PIDs / quotients / tensor products, free_module_element, LLL, BKZ, CVP, SVP, Voronoi) | `intersection()` basis normalisation over non-ZZ PIDs (module is correct; representative differs by a unit) |
| sage.stats.distributions | 2 modules, all four integer algorithms + non-spherical Σ | Arbitrary-precision `_normalisation_factor_zz`; `precision='dp'`; `qfrep` in parigp-ts |
| sage.coding | 4 modules (RS, BCH, Goppa, Reed-Muller) | Brouwer-Zimmermann minimum distance |
| sage.groups.generic | ✅ 95% | `bounds`/`algorithm`/`verify`; `linear_relation`, `merge_points`, `structure_description` |
| sage.quadratic_forms | ✅ 97% BinaryQF (incl. `solve_integer`) | TernaryQF, `BQFClassGroup` class, Shanks distance forms |
| sage.rings.number_field | ✅ 90% | `bnfinit` (degree > 2 class groups), archimedean embeddings, `galoisinit` |
| sage.rings.real_mpfr / complex_mpfr | ✅ 95% / 90% | IEEE 754 double precision (see DEVIATIONS.md) |
| sage.rings.padics | ✅ 90% | Extension fields (Zq) as a first-class module; minimal_polynomial/charpoly |

**Test Coverage (2026-07-28, 0.0.12, verbatim from the runners):**
- `bun test`: **6216 pass, 32 skip, 0 fail**, 2 738 804 expect() calls across 106 files
  (0.0.11 baseline: 6208 pass / 33 skip / 0 fail across 105 files)
- `bun run test:property` (SageMath transcript comparison): **433/433 passed, 0 failed, 0 errors**
- Typecheck: `flint-ts` **0 errors**, `ntl-ts` **0 errors**, `parigp-ts` **29 errors byte-identical
  to the HEAD baseline** (all in test files: 2 `vitest` module resolution, 27 `EllipticPoint`
  union narrowing). `sagemath-ts` has no `tsconfig.json` so it is out of scope; under an ad-hoc
  strict config its error set is **identical before and after** this pass
- Test files created in 0.0.12 for packages/modules that had **none**: `flint-ts/index`
- Test files created during the 0.0.11 audit fix pass:
  `randstate`, `conway_polynomials`, `gf2`, `integer_mod`, `integer_mod_ring`, `rational_field`,
  `matrix_generic`, `matrix_space`, `matrix_integer`, `matrix_modn`, `free_module_element`,
  `free_module_integer`, `unit_group`, `pari_nf`, `polynomial_element`, `polynomial_ring`,
  `quotient_ring`, `multi_polynomial_ideal`, `formal_group`, `constructor`, `cm`,
  `ntl-ts/GF2X`

**ZK-Specific Features:**
- FFT/NTT with O(n log n) polynomial multiplication
- Lagrange interpolation, vanishing polynomials
- Multivariate polynomials with term orders
- KZG helpers: quotient polynomials, batch openings
- FRI helpers: folding, coset domains, proximity testing
- Reed-Solomon codes with error correction

## Tooling

| Tooling | Status | Tests | Notes |
|---------|--------|-------|-------|
| Benchmarks | ✅ | - | SageMath vs sagemath-ts harness in tests/bench |
| zksecurity-cheatsheets | ✅ | - | Curve parameter catalog for benchmarks/tests |
| Tutorials | ✅ | ✅ | Bun-only runner/tests; imports normalized to sagemath-ts |
| Playground | ✅ | ✅ (manual) | Browser bundle + CodeMirror highlighting; lessons generated from TS with snippet validation; docs sidebar with search/filters |

### Dependencies Between Modules
```
sage.rings.integer
    └── sage.arith (uses Integer)
        └── sage.rings.finite_rings (uses arith functions)
            └── sage.rings.polynomial (uses finite fields)
                └── sage.schemes.elliptic_curves (uses everything above)
```
