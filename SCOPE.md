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
| `factor` | 🟡 85% | ✅ | Delegates to PARI `Z_factor`, which still gives up on composites past its trial-division bound (see DEVIATIONS.md) |
| `is_prime` | ✅ 100% | ✅ | Delegates to parigp-ts BPSW — **probabilistic**, not the APRCL/ECPP proof Sage's `proof=True` gives |
| `is_pseudoprime` | ✅ 100% | ✅ | Same BPSW entry point as `is_prime`, so the two coincide here |
| `is_prime_power` | ✅ 100% | ✅ | PARI `isprimepower` algorithm (ported in place — see DEVIATIONS.md) |
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
| `finite_field_extension.ts` | 🟡 92% | ✅ | GF(p^n) via polynomial quotient rings. Default modulus outside the Conway table is `first_lexicographic`, not NTL `BuildSparseIrred` / PARI `ffinit` (both delegation targets pending — see DEVIATIONS.md). `minimalPolynomial()` is still the simplified constant-coefficient version |
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
| `polynomial_element.ts` | 🟡 90% | ✅ new | Full arithmetic + `factor()`, `roots()`, `is_irreducible()` (now Rabin's test — the old `x^(p^n)=x` check passed every fully-split polynomial). **Known bug:** factoring over QQ throws because `QQ.__call__` has no object form; this blocks `minpoly` over QQ for non-squarefree charpolys |
| `polynomial_ring_constructor.ts` | ✅ 100% | ✅ | PolynomialRingConstructor() returning [R, x] |
| `quotient_ring.ts` | ✅ 100% | ✅ new | R[x]/<f(x)> for field extensions |
| `convolution.ts` | ✅ 100% | ✅ | Full port of `convolution.py`'s ring-agnostic Schönhage algorithm, plus the FFT/NTT fast path |
| `multi_polynomial_ring.ts` | ✅ 95% | ✅ | Multivariate polynomials R[x,y,z,...]. `__call__` from a univariate polynomial throws (needs `_mpoly_dict_recursive`) |
| `multi_polynomial_element.ts` | 🟡 60% | ✅ | Sparse multivariate + sumcheck/GKR methods. ~19 of Sage's `MPolynomial` methods are now honest `SAGE_NOT_IMPLEMENTED` stubs rather than absent; only 3 of Sage's 12 term orders |
| `multi_polynomial_ideal.ts` | 🟡 70% | ✅ new | Buchberger only. Raises rather than truncating or hanging; `dimension()` now the real Cox-Little-O'Shea algorithm. `multi_polynomial_ideal.test.ts` created, pinning Sage's Katsura-3 lex basis |
| `polynomial_commitment.ts` | ✅ 100% | ✅ | KZG/FRI helpers — **port-only, no SageMath counterpart**; belongs under `src/zk/` (see DEVIATIONS.md) |
| `polynomial_gf2x.ts` | ✅ 100% | ✅ 46 tests | GF(2)[x] bit-packed: mul, div, gcd, factorization. Should delegate to ntl-ts rather than reimplement |

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
| Audit (algorithm fidelity) | ✅ | - | 2026-07 (`AUDIT-2026-07.md`): 57 findings across the matrix modules, all fixed. 691 tests pass in `matrix/` |
| `matrix_space.ts` | ✅ 90% | ✅ new 17 | MatrixSpace, Matrix class, scalar-matrix `__call__`. `matrix_space.test.ts` created |
| `matrix_generic.ts` | ✅ 90% | ✅ new 14 | Constructors, arithmetic, `pow` with negative exponents. `matrix_generic.test.ts` created |
| `matrix_operations.ts` | 🟡 85% | ✅ 202 | `minpoly` (was returning the minpoly of e₀), division-free `determinant`/`adjugate` over Z/8, exact `is_positive_(semi)definite`, `is_similar` via elementary divisors, `right_kernel_matrix`, `density`, `eigenvalues`. `norm(2)` (SVD), `is_similar` transformation and foreign `base_field` raise |
| `matrix_integer.ts` | ✅ 95% | ✅ new 23 | HNF, SNF (no longer loops forever), elementary divisors, kernel, exact integral LLL (`delta` 0.99), symplectic form, Frobenius form flags 0/1. `matrix_integer.test.ts` created |
| `matrix_modn.ts` | 🟡 85% | ✅ new 20 | charpoly (Faddeev-LeVerrier sign fixed), determinant, echelonize, `right_kernel_matrix` (all three basis formats). Composite modulus raises — needs PARI `matkermod`. `matrix_modn.test.ts` created |
| `matrix_mod2.ts` | ✅ 95% | ✅ | GF(2) matrices; `pluq`/`ple` now use M4RI's transposition-list convention for both P and Q |
| `matrix_decompositions.ts` | ✅ 95% | ✅ 104 | RREF echelon form, LU, QR, Cholesky, Bunch-Kaufman `block_ldlt`, Smith, Hermite, Jordan, LLL_gram, Krylov. `pivot_rows` now returns rows (and a new `pivots` returns columns) |
| `matrix_decompositions_additions.ts` | ✅ 100% | ✅ | SVD_double, QR_double, LU_double for IEEE 754 real matrices |
| `matrix_special.ts` | ✅ 95% | ✅ 98 | `companion_matrix`, `toeplitz`, `hankel`, `elementary_matrix`, `block_matrix`, `rook_vector`, `berlekamp_massey`, `is_permutation_of`, `permutation_normal_form`, random matrix constructors. (This module was absent from this table entirely.) |

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
| `formal_group.ts` | ✅ 90% | ✅ new 20 | `differential()`, `log()`, `inverse()`, `group_law()`, `mult_by_n()`, `sigma()` implemented for real (were placeholders emitting four hardcoded coefficients). `formal_group.test.ts` created (module previously had **none**) |
| `ell_torsion.ts` | ✅ 90% | ✅ | `_p_primary_torsion_basis` replaced with Sage's division-polynomial algorithm. Torsion over number fields still throws |
| `weierstrass_morphism.ts` | ✅ 100% | ✅ | `order()`, automorphism enumeration, isomorphism composition |
| `isogeny_class.ts` | 🟡 70% | ✅ new | IsogenyClass, IsogenyClassNumberField, IsogenyClassRational, `fill_isogeny_matrix`. `isogeny_degrees_cm` omits Sage's `Frobenius_filter` and `possible_isogeny_degrees` uses Mazur's list — both **supersets** of Sage's answer (see DEVIATIONS.md) |
| `cm.ts` | ✅ 95% | ✅ new 18 | `cm_j_invariants`, `cm_orders`, `is_cm_j_invariant`, `discriminants_with_bounded_class_number` (6.26 s -> 47 ms), `largest_(fundamental_)disc_with_class_number`. All four return shapes corrected to Sage's. `cm.test.ts` created (module previously had **none** despite being listed as delivered) |
| `padic_lseries.ts` | 🔴 30% | ✅ 68 tests | Blocked: requires modular symbols (sage.modular.modsym). Has: class structure, `teichmuller` (correct), `alpha` (good-ordinary + multiplicative), `_e_bounds` (matches Sage's doctest vectors), validation. `bernardi_sigma_function` and `_c_bound` now **throw** rather than return wrong values |

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
| `lattice.ts` | ✅ 95% | ✅ | `gen_lattice` with all four types (`modular`, `ideal`, `cyclotomic`, `random`) and the `ntl`/`lattice` output flags. Seeded output still differs from Sage's (two different upstream generators — see DEVIATIONS.md) |
| `lwe.ts` | ✅ 98% | ✅ | LWE, Regev, LindnerPeikert, RingLWE, RingLindnerPeikert, RingLWEConverter. Sage's doctest parameters now reproduce exactly |
| `boolean_function.ts` | ✅ 100% | ✅ | Walsh transform, nonlinearity, ANF, correlation immunity, algebraic immunity, hex truth tables |
| `sbox.ts` | ✅ 100% | ✅ | DDT, LAT (per-mask Walsh-Hadamard; AES 175 ms -> 6.5 ms), APN detection, `min_degree`, MISTY/Feistel constructions |

### `sage.modules` - Lattices and Modules
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | 2026-07 (`AUDIT-2026-07.md`): 27 findings, all fixed. 255 tests pass in `modules/` |
| `free_module.ts` | ✅ 90% | ✅ | FreeModule hierarchy matching Sage's (`ambient_pid`/`submodule_pid`/`submodule_field`), exact rank, echelonized bases, coordinates, kernels, saturation, discriminant, cardinality. The float RREF/kernel/determinant helpers are gone. Non-ZZ PIDs (e.g. `QQ[x]`) raise |
| `free_module_element.ts` | ✅ 95% | ✅ new 23 | Inner-product-matrix pairing, 7-D cross product, true p-norms, exact normalization, Python-style indexing. `free_module_element.test.ts` created |
| `free_module_integer.ts` | ✅ 90% | ✅ new 41 | Exact integral LLL, exact CVP/SVP enumeration, Voronoi relevant vectors and cell, q-ary lattices. `free_module_integer.test.ts` created. Note: this and `matrix_integer.LLL` are two independent LLL implementations where SageMath has one |
| `bkz.ts` | ✅ 100% | ✅ 34 | BKZ reduction, HKZ, Schnorr-Euchner enumeration (no audit findings) |

### `sage.stats.distributions` - Sampling Distributions
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| `discrete_gaussian_integer.ts` | 🟡 70% | ✅ | Rejection sampling, `uniform+table` and `uniform+online`. **Not implemented:** `uniform+logtable` / `sigma2+logtable` (need `dgs_bern.c`), the `dp` precision mode, `_flush_cache`. Rounding of a non-integer centre corrected to dgs's round-half-to-even |
| `discrete_gaussian_lattice.ts` | 🟡 60% | ✅ | GPV algorithm, coset sampling, now exact (`Rational` basis/centre/GSO) with `sampleExact()`. **Not implemented:** non-spherical Σ (matrix sigma, Peikert `r`, Cholesky, offline samples, `_call_non_spherical`), `set_c`/`c()`/`sigma()`/`f()`, `_normalisation_factor_zz` (PARI `qfrep`) — so Sage's `DGL(ZZ^3, Matrix(...), [7,2,5])` doctest is inexpressible |

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
| `misc/randstate.ts` | ✅ 95% | ✅ 24 tests | Centralized RNG + `set_random_seed` parity. **2026-07:** the 64-bit MMIX LCG (whose low bits had period 2^(k+1), so every consumer emitted short deterministic cycles) was replaced with MT19937 — the generator family GMP's `gmp_randinit_default` uses — with GMP-compatible per-call bit consumption. The seeding step still differs from GMP's `randseed_mt`, so streams differ from Sage's for identical seeds. `randstate.test.ts` created; the module had **zero** tests before, which is why the defect survived |
| `quadratic_forms/binary_qf.ts` | ✅ 95% | ✅ 77 tests | BinaryQF: composition (PARI `qfb_comp`/`qfb_sqr`), reduction (`qfi_redsl2_basecase`, `_reduce_indef`), cycles, `is_equivalent`, `BinaryQF_reduced_representatives`. Gauss squaring and `reduced_form` were both wrong before 2026-07 |
| `groups/generic.ts` | ✅ 95% | ✅ 106 tests | Sage's `discrete_log` loop verbatim (incl. repair of a non-minimal `ord` and the `<30` linear branch of `bsgs`); `order_from_multiple` honours `check=True` |
| `zk/sumcheck.ts`, `zk/multilinear.ts` | ✅ 90% | ✅ 106 tests | Ports of `reference/sage_blueprints/`, **not** of SageMath. Round polynomials are now built symbolically at arbitrary degree; `sumcheckVerify` requires `numVars` (a short proof previously verified nothing) |

---

## Dependency Libraries Progress

### parigp-ts (port of cypari2 / PARI/GP)

**Total: 414 tests passing, 0 failing**

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
- SEA algorithm for large prime cardinality (`ellcard` uses a faithful Shanks BSGS for all p above 2048; PARI switches to SEA at `expi(p) >= 56`, so 56..126-bit primes have an `O(p^(1/4))` table here). PARI's own "install the seadata package" overflow guard is reproduced
- `Fp_ellcard_CM` beyond the `a6 = 0` (j = 1728) branch — performance only; Shanks is verified correct on those curves
- `nf` module (nfbasis, nfdisc, idealprimedec, nfgaloisconj) — currently ported inside sagemath-ts
- `Qfb` module (qfbcompraw, qfbred, qfbredsl2, qfbsolve, qfbclassno) — currently ported inside sagemath-ts
- Matrix routines (`matkermod`, `matfrobenius`) — blocks two sagemath-ts matrix paths
- `ffinit` — blocks the GF(p^n) default-modulus delegation
- Transcendental functions (`dilog`, `incgam`) — currently ported inside sagemath-ts
- Extension field elliptic curves (t_ELL_Fq)
- Integer arithmetic (factor, is_prime, etc.) - use sagemath-ts/arith instead. `Z_factor` still gives up on composites past its trial-division bound, which constrains `is_prime_power`, `IntegerMod.log()` and `arith.factor`
- Polynomial operations (t_POL)

### flint-ts (port of FLINT)
| Feature | Status | Tests | Notes |
|---------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | 2026-07: **all APIs are still stubs** (FLINT_NOT_IMPLEMENTED). Consequence: sagemath-ts reimplements FLINT primitives rather than delegating, which the audit identifies as a recurring structural cause of defects |
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
| GF2X | ✅ 85% | ✅ 15 | Bit-packed bigint representation. `IterIrredTest`, `BuildIrred`, `BuildSparseIrred` ported line-for-line from `GF2XFactoring.cpp` over a vendored copy of NTL's 2049-row `GF2X_irred_tab`; verified against Sage's `polynomial_gf2x.pyx` doctests and exhaustive brute force for degree <= 10. **Not implemented:** `random`, `factor`, `SquareFreeDecomp`, `DistinctDegFactor`, `EqualDegFactor`, `BerlekampFactor` (need NTL's ChaCha `RandomStream` / randomized factoring). **Consumers have not yet been wired up:** `finite_field_extension.getDefaultModulus` and `polynomial_gf2x.ts` should delegate here |

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
| `BinaryQF` | ✅ 95% | MEDIUM | Binary quadratic forms ax^2 + bxy + cy^2 — implemented; see the Cross-Cutting table above. (This row said "not started" while the module existed with two critical defects.) |
| `TernaryQF` | ⬜ | LOW | Ternary quadratic forms |
| `BQFClassGroup` | 🟡 60% | MEDIUM | Class group of binary QFs — reduced representatives now form a genuine class group (28 discriminants verified against the literature); no dedicated `BQFClassGroup` class |
| `qfsolve` | ⬜ | LOW | Solve quadratic equations |
| `least_quadratic_nonresidue` | ⬜ | LOW | Find smallest QNR mod p |

**Use Cases:** Class group cryptography, lattice reduction, quadratic sieve

### `sage.rings.number_field` - Algebraic Number Theory (PARTIAL)
Needed for advanced algebraic crypto constructions.

| Module | Status | Priority | Notes |
|--------|--------|----------|-------|
| Audit (algorithm fidelity) | ✅ | - | 2026-07 (`AUDIT-2026-07.md`): 17 findings, 16 fixed. 348 tests (was 228) |
| `pari_nf.ts` | ✅ 85% | - | **New in 0.0.11.** The number-field kernel: `nfbasis`/`nfdisc` (Round 2), `idealprimedec` (Dedekind-Kummer), `nfgaloisconj` (p-adic reconstruction, degree <= 8), `polisirreducible` (Cantor-Zassenhaus / Zassenhaus). Belongs in parigp-ts — see DEVIATIONS.md |
| `NumberField` | ✅ 85% | LOW | Basic operations, element arithmetic, norm, trace, real field discriminants, integral bases |
| `NumberFieldElement` | ✅ 85% | LOW | Arithmetic, charpoly, minpoly, is_integral, `is_unit` (field semantics) + `is_integral_unit` |
| `CyclotomicField` | ✅ 85% | MEDIUM | Q(zeta_n) - cyclotomic polynomials, degree, exact `automorphisms()` at any degree |
| `QuadraticField` | ✅ 85% | LOW | Q(sqrt(D)) using `x^2 - D` verbatim as Sage does (it previously rewrote the polynomial to the squarefree part), discriminant |
| `RationalPolynomial` | ✅ 90% | - | Supporting polynomial arithmetic over Q |
| `class_group` | 🟡 55% | LOW | **Quadratic** class groups computed exactly (imaginary and real, incl. Sage's C38 x C2 doctest) by enumerating binary quadratic forms, guarded at \|D\| <= 2e6. Degree > 2 still needs PARI `bnfinit` |
| `unit_group` | 🟡 40% | LOW | Torsion subgroup and log embedding correct (the embedding was wrong before); real quadratic fundamental units still need PARI `quadunit`, so `regulator()` throws unless a unit is supplied. `unit_group.test.ts` created — the module had **no** test file, which is how two defects shipped |
| `galois_group` | 🟡 60% | LOW | Built from genuine automorphisms. Throws for non-Galois fields (Sage returns the Galois closure's group via `galoisinit`) |
| `number_field_ideal` | ✅ 85% | LOW | Ideal arithmetic in HNF: rational norms, `smallest_integer`, `is_prime`, canonical equality, `prime_above`/`primes_above`/`decomposition` |
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

**Not Implemented (Requires PARI bnfinit / galoisinit / quadunit):**
- Class group for degree > 2
- Fundamental unit of a real quadratic field, and hence its regulator
- Galois group of a non-Galois field
- `nfgaloisconj` above degree 8; `decomposition(p)` at inessential discriminant divisors

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
`DEVIATIONS.md` at the project root is the single source of truth (62 sections as of 0.0.11).
Headlines:
- `GF()` returns `FiniteFieldPrime` for primes; `FiniteField` is aliased to `GFExtended` and accepts prime powers
- Element coercion is explicit via `__call__()` rather than automatic
- Uses TypeScript's `bigint` for arbitrary precision integers
- The RNG is MT19937 with GMP-compatible bit consumption but a different seeding step, so seeded streams differ from Sage's
- Several PARI/NTL routines are ported *in place* rather than delegated, because the target package lacks the module

### Known Limitations
- Conway database covers p = 2 to n=64, 3 to 24, 5 to 18, 7 to 14, 11/13 to 12, 17/19/23/29/31 to 10; outside it the default modulus is `first_lexicographic`, not NTL/PARI's choice
- **No discrete_log on GF(p^n) extension field elements** - extension field elements lack log() method (Z/nZ has log())
- **Pollard's rho for DLP** - implemented for prime order groups (sage.groups.generic.discrete_log_rho)
- **BKZ uses Schnorr-Euchner enumeration** - pure TypeScript implementation (no fpylll/NTL bindings)
- **`flint-ts` is 100% stubs and `ntl-ts` is stubs apart from GF2/GF2X**, so sagemath-ts reimplements FLINT/NTL primitives instead of delegating — the audit identifies this as a recurring structural cause of defects
- **`parigp-ts` covers only elliptic curves, finite fields and integer factorisation**, so PARI's `nf`, `Qfb`, matrix and transcendental routines are ported inside sagemath-ts (see DEVIATIONS.md)
- **Factoring over QQ is broken** (`QQ.__call__` has no object form), which blocks `minpoly` over QQ for non-squarefree characteristic polynomials
- **`packages/sagemath-ts` has no `tsconfig.json`**, so its `typecheck` script has never actually typechecked the package; against a synthesized strict config it reports ~1284 errors, the large majority pre-dating this pass (mostly `FiniteFieldElement`/`IntegerMod` not satisfying the `this`-typed `RingElement` constraint). Needs a dedicated task
- **`packages/parigp-ts` reports 29 pre-existing type errors** (27 `EllipticPoint` union narrowing, 2 missing `vitest` types); all pass at runtime under Bun

### Implementation Summary

**Phase 1 Progress: ~95% complete**
**Phase 2 Progress: ~95% complete**
**Phase 3 Progress: ~93% complete**

> Percentages were revised downward on 2026-07-28. Nothing regressed — the July 2026 audit found
> 370 confirmed defects behind the previous figures, and the fix pass raised *fidelity* while
> making the remaining gaps explicit (honest `SAGE_NOT_IMPLEMENTED` stubs where there used to be
> silently wrong answers).

| Category | Implemented | Remaining |
|----------|-------------|-----------|
| sage.arith | 47+ functions (incl. rational_reconstruction, CRT_basis, continued_fraction, trial_division, prime_powers, smooth_part, sum_of_squares) | Full `Z_factor` (ECM/MPQS); proven primality |
| sage.rings.finite_rings | 9 modules (incl. tower_field, roots_of_unity, log in Z/nZ) | element.log() for GF(p^n); real `minimalPolynomial`; NTL/PARI modulus delegation |
| sage.rings.polynomial | 10 modules (incl. FFT, Schönhage convolution, multivariate, ideals, commitment, GF2X) | 9 term orders; factoring over QQ; Singular-grade Gröbner |
| sage.matrix | 9 modules (incl. HNF, SNF, mod2, modn, special, decompositions) | SVD/`norm(2)`, `rational_form`, `jordan_form` transformation, `matkermod`, `matfrobenius` flag 2 |
| sage.schemes.elliptic_curves | 11 modules (incl. pairings, torsion, twists, isogeny graph, formal group, CM) | Stark's algorithm, `ell_wp`, `isogeny_small_degree`, modular symbols |
| sage.crypto | 4 modules (lattice, lwe, boolean_function, sbox) | Seeded parity with Sage's RNG streams |
| sage.modules | 4 modules (free_module, free_module_element, LLL, BKZ, CVP, SVP, Voronoi) | Free modules over non-ZZ PIDs |
| sage.stats.distributions | 2 modules (discrete Gaussian integer/lattice) | Logtable algorithms, non-spherical Σ, `_normalisation_factor_zz` |
| sage.coding | 4 modules (RS, BCH, Goppa, Reed-Muller) | Brouwer-Zimmermann minimum distance |
| sage.groups.generic | ✅ 95% | `bounds`/`algorithm`/`verify`; `linear_relation`, `merge_points`, `structure_description` |
| sage.quadratic_forms | ✅ 95% BinaryQF | TernaryQF, `BQFClassGroup` class, `qfsolve` |
| sage.rings.number_field | ✅ 85% | `bnfinit` (degree > 2 class groups), `quadunit`, `galoisinit` |
| sage.rings.real_mpfr / complex_mpfr | ✅ 95% / 90% | IEEE 754 double precision (see DEVIATIONS.md) |
| sage.rings.padics | ✅ 90% | Extension fields (Zq), minimal_polynomial/charpoly |

**Test Coverage (2026-07-28, verbatim from the runners):**
- `bun test`: **5782 pass, 33 skip, 0 fail**, 1 253 493 expect() calls across 102 files
- `bun run test:property` (SageMath transcript comparison): **433/433 passed, 0 failed, 0 errors**
- Test files created during the audit fix pass for modules that had **none**:
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
