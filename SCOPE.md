# Project Scope

This document tracks implementation progress. Update this file when completing modules.

## Status Legend

- ⬜ Not started
- 🟡 In progress (note: include agent/person working on it)
- ✅ Complete (with test coverage %)
- 🔴 Blocked (note: reason)

**Maintenance notes:** 2026-01-30 consolidated deviations log into root `DEVIATIONS.md`.

---

## Phase 1: Core Number Theory (Cryptography Focus)

### `sage.rings.integer` - Arbitrary Precision Integers
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | Reviewed vs reference (2026-02-04) |
| `integer.py` | ✅ 98% | ✅ | Core Integer class - 47+ methods |
| `integer_ring.py` | ✅ 98% | ✅ | ZZ ring with full number-theoretic operations |
| `rational.py` | ✅ 98% | ✅ | Rational numbers - 60+ methods |
| `rational_field.py` | ✅ 98% | ✅ | QQ field with iteration, Selmer groups, quadratic defect |

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
| `real_mpfr.ts` | ✅ 95% | ✅ | RealField, RealNumber with IEEE 754 double precision |

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
| Audit (algorithm fidelity) | ✅ | - | Reviewed vs reference (2026-02-04) |
| `gcd` | ✅ 100% | ✅ | Binary GCD (Stein's algorithm) |
| `lcm` | ✅ 100% | ✅ | |
| `xgcd` | ✅ 100% | ✅ | Extended Euclidean algorithm |
| `factor` | ✅ 95% | ✅ | Trial division (TODO: add Pollard rho, ECM) |
| `is_prime` | ✅ 100% | ✅ | Deterministic Miller-Rabin |
| `is_prime_power` | ✅ 100% | ✅ | With optional data extraction |
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
| Audit (algorithm fidelity) | ✅ | - | Reviewed vs reference (2026-02-04) |
| `finite_field_constructor.ts` | ✅ 100% | ✅ | GF() constructor for prime fields |
| `integer_mod.ts` | ✅ 100% | ✅ | Z/nZ elements with Mod() function |
| `integer_mod_ring.ts` | ✅ 100% | ✅ | Zmod() ring constructor with iteration |
| `finite_field_prime.ts` | ✅ 100% | ✅ | GF(p) with sqrt, multiplicative_generator |
| `finite_field_extension.ts` | ✅ 100% | ✅ | GF(p^n) via polynomial quotient rings |
| `gf2.ts` | ✅ 100% | ✅ | Optimized GF(2) singleton |
| `conway_polynomials.ts` | ✅ 100% | ✅ | Database for p=2,3,5,7,11,13,17,19,23,29,31 |
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
| Audit (algorithm fidelity) | ✅ | - | Reviewed vs reference (2026-02-04) |
| `polynomial_ring.ts` | ✅ 100% | ✅ | PolynomialRing with lagrange, vanishing, cyclotomic |
| `polynomial_element.ts` | ✅ 100% | ✅ | Full arithmetic + factor(), roots(), is_irreducible() |
| `polynomial_ring_constructor.ts` | ✅ 100% | ✅ | PolynomialRingConstructor() returning [R, x] |
| `quotient_ring.ts` | ✅ 100% | ✅ | R[x]/<f(x)> for field extensions |
| `convolution.ts` | ✅ 100% | ✅ | FFT/NTT, fast multiplication O(n log n) |
| `multi_polynomial_ring.ts` | ✅ 100% | ✅ | Multivariate polynomials R[x,y,z,...] |
| `multi_polynomial_element.ts` | ✅ 100% | ✅ | Sparse multivariate + sumcheck/GKR methods |
| `polynomial_commitment.ts` | ✅ 100% | ✅ | KZG/FRI helpers: quotients, linearization |
| `polynomial_gf2x.ts` | ✅ 100% | ✅ 46 tests | GF(2)[x] bit-packed: mul, div, gcd, factorization |

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
| `matrix_space.ts` | ✅ 80% | ✅ | MatrixSpace, Matrix class |
| `matrix_operations.ts` | ✅ 70% | ✅ | Basic operations, transpose |
| `matrix_integer.ts` | ✅ 95% | ✅ | HNF, SNF, elementary divisors, kernel |
| `matrix_decompositions.ts` | ✅ 95% | ✅ | LU, QR, Cholesky, Smith, Hermite, Jordan, LLL_gram, Krylov methods |
| `matrix_decompositions_additions.ts` | ✅ 100% | ✅ | SVD_double, QR_double, LU_double for IEEE 754 real matrices |
| `matrix_mod2_dense.py` | ⬜ | - | GF(2) matrices |

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
| Audit (algorithm fidelity) | ✅ | - | Reviewed vs reference (2026-02-04) |
| `constructor.ts` | ✅ 100% | ✅ | EllipticCurve() - delegates to parigp-ts ellinit |
| `ell_generic.ts` | ✅ 95% | ✅ | Invariants, torsion_points, is_on_curve |
| `ell_finite_field.ts` | ✅ 98% | ✅ | cardinality, trace, generators, twists, torsion_basis; FpE_mul uses Jacobian |
| `ell_point.ts` | ✅ 98% | ✅ | Point arithmetic, weil_pairing, tate_pairing |
| `ell_curve_isogeny.ts` | ✅ 95% | ✅ | Velu's formulas, isogeny_class, is_isogenous |
| `isogeny_class.ts` | ✅ 70% | ✅ | IsogenyClass, IsogenyClassNumberField, IsogenyClassRational |
| `padic_lseries.ts` | 🔴 35% | ✅ 68 tests | Blocked: requires modular symbols (sage.modular.modsym). Has: class structure, teichmuller, _e_bounds, validation |

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
| `lattice.ts` | ✅ 80% | ✅ | IntegerLattice with LLL reduction |
| `lwe.ts` | ✅ 98% | ✅ | LWE, Regev, LindnerPeikert, RingLWE, RingLindnerPeikert, RingLWEConverter |
| `boolean_function.ts` | ✅ 100% | ✅ | Walsh transform, nonlinearity, ANF, correlation immunity |
| `sbox.ts` | ✅ 100% | ✅ | DDT, LAT, APN detection, standard S-boxes (AES, PRESENT) |

### `sage.modules` - Lattices and Modules
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| `free_module.ts` | ✅ 70% | ✅ | FreeModule base class |
| `free_module_integer.ts` | ✅ 90% | ✅ | LLL algorithm, Gram-Schmidt |
| `bkz.ts` | ✅ 100% | ✅ | BKZ reduction, HKZ, Schnorr-Euchner enumeration |

### `sage.stats.distributions` - Sampling Distributions
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| `discrete_gaussian_integer.ts` | ✅ 100% | ✅ | Rejection sampling, table/online algorithms |
| `discrete_gaussian_lattice.ts` | ✅ 100% | ✅ | GPV algorithm, coset sampling |

### `sage.coding` - Error-Correcting Codes
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| `reed_solomon.ts` | ✅ 100% | ✅ | RS codes with encode/decode, FRI fold/query |
| `bch_code.ts` | ✅ 100% | ✅ | BCH codes, PGZ decoding, Chien search |
| `goppa_code.ts` | ✅ 100% | ✅ | Goppa codes, Patterson algorithm (McEliece-ready) |
| `reed_muller_code.ts` | ✅ 100% | ✅ | RM(r,m) codes, Plotkin construction, majority decoding |

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
| `misc/randstate.ts` | ✅ 100% | ⚠️ not run | Centralized RNG + set_random_seed parity |

---

## Dependency Libraries Progress

### parigp-ts (port of cypari2 / PARI/GP)

**Total: 376 tests passing (2001 expect() calls)**

| Feature | Status | Tests | Notes |
|---------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | Reviewed vs reference (2026-02-04) |
| Core types (GEN, t_INT, etc.) | ✅ | - | types.ts - PariType enum, PariInt, PariFfelt, PariVec, etc. |
| Fp arithmetic | ✅ | 55 | ff.ts - Fp_add, Fp_sub, Fp_mul, Fp_sqr, Fp_neg, Fp_inv, Fp_div, Fp_pow |
| Fp_sqrt (Tonelli-Shanks) | ✅ | - | ff.ts - includes Fp_issquare, kronecker symbol |
| ellinit | ✅ | 35 | elliptic/init.ts - Short/general Weierstrass, from j-invariant |
| Point operations (elladd, ellmul) | ✅ | 65 | elliptic/point.ts - Jacobian coordinates for efficiency |
| ellcard, ellgroup, ellorder | ✅ | 36 | elliptic/group.ts - Baby-step giant-step for medium primes |
| ellordinate, random_FpE | ✅ | 32 | elliptic/points.ts - Find y from x, random point generation |
| ellgenerators, trace_of_frobenius | ✅ | - | elliptic/group.ts |

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
- SEA algorithm for large prime cardinality (ellcard uses BSGS up to ~10^6, then slower BSGS for larger)
- Extension field elliptic curves (t_ELL_Fq)
- Integer arithmetic (factor, is_prime, etc.) - use sagemath-ts/arith instead
- Polynomial operations (t_POL)

### flint-ts (port of FLINT)
| Feature | Status | Tests | Notes |
|---------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | All APIs are stubs (FLINT_NOT_IMPLEMENTED) |
| fmpz (integers) | ⬜ | - | |
| fmpz_poly | ⬜ | - | |
| fmpz_mod_poly | ⬜ | - | |
| nmod_poly | ⬜ | - | |

### ntl-ts (port of NTL)
| Feature | Status | Tests | Notes |
|---------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | All APIs are stubs (NTL_NOT_IMPLEMENTED) |
| ZZ | ⬜ | - | |
| ZZ_p | ⬜ | - | |
| ZZ_pX | ⬜ | - | |
| GF2 | ⬜ | - | |
| GF2X | ⬜ | - | |

---

## Future Cryptography Modules (Not Yet Implemented)

This section tracks cryptography-relevant SageMath functionality that could be added.

### `sage.groups.generic` - Generic Group Operations
Generic algorithms for any group with compatible operations.

| Function | Status | Priority | Notes |
|----------|--------|----------|-------|
| Audit (algorithm fidelity) | ✅ | - | Reviewed vs reference (2026-02-04) |
| `discrete_log` | ✅ 90% | HIGH | Pohlig-Hellman + BSGS; no bounds/algorithm options |
| `discrete_log_rho` | ✅ 90% | HIGH | Pollard's rho; requires explicit prime order |
| `discrete_log_lambda` | ✅ 100% | MEDIUM | Pollard's kangaroo (bounded DLP) |
| `bsgs` | ✅ 100% | HIGH | Baby-step giant-step algorithm |
| `pohlig_hellman` | ✅ 100% | HIGH | Reduce DLP to prime power subgroups |
| `order_from_multiple` | ✅ 90% | HIGH | Compute element order given a multiple; no check/plist |
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
| `BinaryQF` | ⬜ | MEDIUM | Binary quadratic forms ax^2 + bxy + cy^2 |
| `TernaryQF` | ⬜ | LOW | Ternary quadratic forms |
| `BQFClassGroup` | ⬜ | MEDIUM | Class group of binary QFs |
| `qfsolve` | ⬜ | LOW | Solve quadratic equations |
| `least_quadratic_nonresidue` | ⬜ | LOW | Find smallest QNR mod p |

**Use Cases:** Class group cryptography, lattice reduction, quadratic sieve

### `sage.rings.number_field` - Algebraic Number Theory (PARTIAL)
Needed for advanced algebraic crypto constructions.

| Module | Status | Priority | Notes |
|--------|--------|----------|-------|
| `NumberField` | ✅ 70% | LOW | Basic operations, element arithmetic, norm, trace |
| `NumberFieldElement` | ✅ 70% | LOW | Arithmetic, charpoly, minpoly, is_integral, is_unit |
| `CyclotomicField` | ✅ 80% | MEDIUM | Q(zeta_n) - cyclotomic polynomials, degree |
| `QuadraticField` | ✅ 80% | LOW | Q(sqrt(D)), discriminant |
| `RationalPolynomial` | ✅ 90% | - | Supporting polynomial arithmetic over Q |
| `class_group` | 🟡 20% | LOW | Structure only, computation requires PARI |
| `unit_group` | 🟡 20% | LOW | Structure only, fundamental units require PARI |
| `galois_group` | ⬜ | LOW | Galois group of number field |
| `number_field_ideal` | ⬜ | LOW | Ideal arithmetic, factorization |
| `order` | ⬜ | LOW | Orders in number fields |

**Implemented (Pure TypeScript):**
- Element arithmetic (add, sub, mul, div, pow, inv)
- Norm, trace via characteristic polynomial
- Characteristic polynomial (Faddeev-LeVerrier)
- Minimal polynomial
- Signature via Sturm's theorem
- is_integral, is_unit
- QuadraticField with discriminant
- CyclotomicField with cyclotomic polynomials
- 50 tests passing

**Not Implemented (Requires PARI bnfinit):**
- Class group computation
- Unit group / fundamental units
- Ideal factorization
- Ring of integers / maximal order

**Use Cases:** Ring-LWE over cyclotomic fields, class group crypto, NTRU

### `sage.rings.padics` - p-adic Numbers
| Module | Status | Priority | Notes |
|--------|--------|----------|-------|
| `Zp(p, prec)` | ✅ 90% | LOW | p-adic integers with capped-relative precision |
| `Qp(p, prec)` | ✅ 90% | LOW | p-adic field |
| `pAdicGenericElement` | ✅ 85% | LOW | 35 of 39 methods implemented |
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
- nth_root for n != 2
- log with branch parameter for non-units
- artin_hasse_exp
- minimal_polynomial, charpoly
- Extension fields (Zq, ramified/unramified extensions)

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
| `shortest_vector` | ✅ 100% | HIGH | SVP solver: Schnorr-Euchner enumeration for dim≤50, simple enumeration for dim≤4 |
| `closest_vector` | ✅ 100% | HIGH | CVP solver: closestVector (exact enum for rank≤4, Babai for larger), approximateClosestVector (nearest_plane, rounding_off, embedding) |
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
Document any intentional differences from SageMath here:
- GF() constructor currently returns FiniteFieldPrime for primes, requires GFExtended() for extension fields
- Element coercion is explicit via __call__() rather than automatic
- Uses TypeScript's bigint for arbitrary precision integers

### Known Limitations
- Large Conway polynomials (degree > 20 for p=2) not included in database
- Irreducibility testing uses probabilistic methods for large degrees
- **No discrete_log on GF(p^n) extension field elements** - extension field elements lack log() method (Z/nZ has log())
- **Pollard's rho for DLP** - implemented for prime order groups (sage.groups.generic.discrete_log_rho)
- **BKZ uses Schnorr-Euchner enumeration** - pure TypeScript implementation (no fpylll/NTL bindings)

### Implementation Summary

**Phase 1 Progress: ~98% complete**
**Phase 2 Progress: ~98% complete**
**Phase 3 Progress: ~97% complete**

| Category | Implemented | Remaining |
|----------|-------------|-----------|
| sage.arith | 47+ functions (incl. rational_reconstruction, CRT_basis, continued_fraction) | - |
| sage.rings.finite_rings | 9 modules (incl. tower_field, roots_of_unity, log in Z/nZ) | element.log() for GF(p^n) extensions |
| sage.rings.polynomial | 9 modules (incl. FFT, multivariate, commitment, GF2X) | Grobner |
| sage.matrix | 4 modules (incl. HNF, SNF) | matrix_mod2_dense |
| sage.schemes.elliptic_curves | 5 modules (incl. pairings, torsion, twists, isogeny graph) | weierstrass_morphism (partial) |
| sage.crypto | 4 modules (lattice, lwe, boolean_function, sbox) | gen_lattice (stub) |
| sage.modules | 3 modules (free_module, LLL, BKZ, CVP, SVP) | - |
| sage.stats.distributions | 2 modules (discrete Gaussian integer/lattice) | - |
| sage.coding | 4 modules (RS, BCH, Goppa, Reed-Muller) | - |
| sage.groups.generic | ✅ complete | discrete_log, bsgs, pohlig_hellman, order_from_multiple |
| sage.quadratic_forms | ⬜ not started | BinaryQF, class groups |
| sage.rings.number_field | ✅ 70% | Basic ops (50 tests), advanced features require PARI |
| sage.rings.real_mpfr | ✅ 95% | IEEE 754 double precision (85 tests, see DEVIATIONS.md) |
| sage.rings.padics | ⬜ not started | p-adic numbers |

**Test Coverage:**
- parigp-ts: 376 tests, 2001 expect() calls
- sagemath-ts: 2,503 tests, 16,146 expect() calls
- property tests: 504 comparison tests
- **Total: 3,383+ tests passing**

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
