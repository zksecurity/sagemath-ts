# Changelog

All notable changes to this project will be documented in this file.

## 0.0.16 - 2026-07-29

The **silent-answer cleanup**: the three remaining known wrong-value paths from the 0.0.15 audit
now compute faithfully, the new crypto module families have live Sage coverage, the last 57
vacuous property cases are real comparisons, and routine verification is split from the
long-running research-grade vectors.

### Fidelity fixes

- `IsogenyClass._compute` now performs the full breadth-first isogeny traversal, records actual
  maps and prime/filled degree matrices, discovers characteristic-zero kernel polynomials, and
  normalizes rational codomains to Sage's global minimal models. The complete `11a1` class matches
  Sage curve-for-curve and matrix-for-matrix. Degrees 43/67/163 now fail explicitly on the one
  genuinely missing input: Sage's precomputed exceptional kernel table.
- `p_minimal_polynomials`, `null_ideal` and `integer_valued_polynomials_generators` implement the
  exact Smith-congruence/J-ideal computation, including exponent lifts, composite and negative
  moduli, and rational generators. `is_LLL_reduced` uses exact Rational Gram-Schmidt with Sage's
  `delta=0.99` default and validation instead of IEEE-754 plus an epsilon.
- Extension-field `minpoly`/`minimal_polynomial` delegates to parigp-ts `FpXQ_minpoly`, matching
  Sage's PARI architecture for generators, constants and proper-subfield elements;
  `minimalPolynomial` remains as a compatibility alias.
- The general `GF(p)` element used by rational function fields now exposes `is_square` and `sqrt`;
  this was found by the new function-field oracle rather than hidden by narrowing its cases.

### Differential oracle and test tiers

- Added live areas for hyperelliptic curves (7), quaternion algebras (8), and rational function
  fields (9); extended `matrix_extended` to 913 cases and `ec_advanced` to 250.
- Added real Python/TypeScript area modules for `arith_extended` (26), `lwe` (15), and `matrix`
  (16). The harness now requires both area modules and treats unknown module/function dispatch as
  an error, so a missing implementation cannot pass because both runners failed.
- `bun run test:fast` / `test:slow` form an exhaustive 102/23-file unit partition; the measured fast
  tier completes in about 40 seconds. `test:property:fast` / `test:property:slow` partition the 23
  live areas, while `test:property` still runs everything.
- Live SageMath 10.3 differential result: **4673/4673 passed, 0 failed, 0 errors**.
- Exhaustive fast + slow unit result: **7122 pass, 32 skip, 0 fail**, 3,025,438 `expect()` calls
  across 125 files.

## 0.0.15 - 2026-07-29

The **differential-oracle pass**: nine new property-test areas comparing against a real SageMath
10.3 process, four new ported module families wired into the package surface, and every
port-vs-SageMath disagreement those areas found, fixed. `bun run test:property` goes from **433**
cases to **4643**, all passing. The complete unit suite is **7101 pass / 32 skip / 0 fail** across
125 files.

### New property-test areas (tests/property)

`mpfr` (815), `matrix_extended` (908), `coding_crypto` (926), `padics_series` (462),
`ec_advanced` (249), `groups_modn` (335), `rand_stats` (274), `lattices` (159),
`quadratic_forms` (82).

### New modules wired into the package surface

- `sage.schemes.hyperelliptic_curves` — curves, Jacobians, Mumford divisors, Igusa/Clebsch
  invariants. Exported from `schemes/index.ts` and as the `./schemes/hyperelliptic_curves`
  subpath.
- `sage.algebras.quatalg` — quaternion algebras, orders and fractional ideals over QQ. Exported
  from the package root and as `./algebras` / `./algebras/quatalg`.
- `sage.rings.function_field` — rational function fields, orders, ideals, places, divisors and
  Riemann-Roch. Exported from `rings/index.ts` and as `./rings/function_field`.
- `sage.quadratic_forms` — `QuadraticForm`, the local-field invariants and `TernaryQF`, exported
  as the namespaced `quadratic_forms` from the package root (`RationalMatrix`, `evaluate`,
  `extend` and `primitivize` are too generic to flatten).
- `rings/laurent_series_ring.ts` is now re-exported from `rings/index.ts`.

### sagemath-ts — correctness fixes found by the oracle

**Real and complex numbers (`rings/real_mpfr.ts`, `rings/complex_mpfr.ts`, new
`rings/real_mpfr_dd.ts`)**

- `RealNumber.str()` is now a port of `real_mpfr.pyx:1897` for base 10, and `toString()` is
  `str(truncate=True)`; `ComplexNumber.toString()` composes the parts as
  `complex_mpfr.pyx:1311-1326` does. `RR(1)` prints `1.00000000000000`, not `1`.
- `exact_rational()` reduces to lowest terms, which also repairs `nearby_rational`.
- `sqrt`/`log`/`log2`/`log10`/`log1p` of a negative real widen to the complex field instead of
  returning NaN; `pow` retries over CC on NaN and honours IEEE `pow(1, y) == 1`.
- The field's rounding mode is applied to `div` and `sqrt` (exact error-sign oracle plus a
  one-ulp nudge).
- `sign_mantissa_exponent` branches on the sign BIT and no longer raises on NaN/infinity;
  `is_square(NaN)` is True; `frac(-0.0)` keeps its sign.
- **Bessel and error functions rewritten in double-double arithmetic** (`real_mpfr_dd.ts`): the
  Numerical Recipes rational approximations they used were single-precision (`RR(1).j0()` was
  wrong from the 9th significant digit). Series below `|x| = 17.5`, Hankel asymptotic above.
- `gamma` reproduces factorials exactly; `log_gamma` is exact at the integers, `+infinity` at 0
  and raises at the poles; `zeta` at the negative integers uses `-B_{n+1}/(n+1)`.
- Complex `abs`/`sqrt`/`log` use `Math.hypot` (no more overflow at `1e300`); `arccos`, `arcsin`,
  `arctan`, `arccosh`, `arcsinh`, `arctanh` are exact on their real/imaginary axes and raise at
  PARI's branch-cut endpoints; `arccosh` takes the principal branch for `z < -1`; `gamma` returns
  the unsigned infinity at its poles and uses the `g = 607/128` Lanczos coefficients; `zeta` uses
  Borwein's Algorithm 2; `is_imaginary` is "real part is zero", as upstream.

**p-adics (`rings/padics/`)**

- Division, inversion and negative powers move to the fraction field, as
  `padic_generic_element.pyx:449` does — which also fixes the repr of negative-valuation elements.
- `is_unit` is True for every nonzero element of a field; `lift()` returns a Rational for negative
  valuation; `lift_to_precision` enforces the precision cap; `__getitem__` accepts negative
  indices.
- `nth_root` seeds the Newton iteration with SageMath's residue-field root (a port of
  `element_base.pyx:_nth_root_common`), so it returns SageMath's root, not just *a* root.
- `log()` of an exact 1-unit returns `O(p^aprec)`, not an exact zero.
- `artin_hasse_exp` handles `p = 2` with `x = 2 mod 4`, where `AH(x) = -exp(...)`.

**Power and Laurent series**

- `_repr_` does upstream's string surgery instead of deciding on `coefficient == -1` (over GF(5)
  the coefficient 4 satisfies that), and emits the short `O(1)` / `O(x)` forms.
- Division/inversion by a non-unit lands in the Laurent series ring; division by a unit stays
  exact; `truncate()` returns a polynomial.
- `__getitem__` past the precision raises `IndexError`; `LaurentSeriesRing.characteristic()`
  works over `GF(p)`; `LaurentSeriesElement.__call__` accepts an argument of negative valuation
  when the unit part is exact.

**Matrices**

- `Matrix_mod2_dense.right_kernel_matrix` echelonizes (the default basis over a field).
- `hermite_form` is `_echelon_form_PID`, not the RREF.
- `is_primitive` is Perron-Frobenius primitivity, not "every row has gcd 1".
- `minpoly` refuses characteristic 2 and composite moduli, as `matrix_modn_dense_template.pxi`
  does; `inverse` raises `ZeroDivisionError`; `jordan_form` raises `RuntimeError`; the permanent
  and the backend-specific inverse messages match.
- `Polynomial.roots()` returns SageMath's order (`Factorization.sort`: multiplicity ascending,
  then the root descending), which is what Jordan block ordering keys off.
- `LLL` of a linearly dependent generating set reproduces fpLLL's basis when the independent
  prefix already spans the row lattice.

**Lattices**

- `shortestVector` tracks the argmin row, so it can no longer return a non-shortest vector.
- BKZ/HKZ now run a genuine Schnorr-Euchner tour with **exact** block SVP enumeration and a
  unimodular insertion, so `HKZ()` really does realise `lambda_1`.
- `volume()` of a non-full-rank lattice returns an exact `sqrt(N)` instead of a floored integer,
  which also repairs `isUnimodular()`.
- The constructor rejects a linearly dependent basis; `BKZ({blockSize: 1})` is accepted.

**Elliptic curves**

- `compute_isogeny_stark` and the `weierstrass_p` it needs are ported, and
  `compute_isogeny_kernel_polynomial` follows SageMath's `ell < 10` dispatch. BMSS and Stark do
  NOT agree for even degrees, so `dual()` of every even-degree isogeny was wrong.
- `is_kernel_polynomial` (odd-degree validation) is ported.
- `EllipticCurveHom.formal()` is implemented on top of the formal group; it used to return the
  series `t` for every isogeny.
- `hilbert_class_polynomial` delegates to parigp-ts's `polclass0` instead of a nine-entry lookup
  table that returned an unprintable object.
- `_equation_string` implements the `±1` special cases, so every repr and every error message
  embedding the equation matches.
- `EllipticCurveTorsionSubgroup.invariants()` returns increasing invariant factors.
- Point construction raises SageMath's `TypeError` with the projective coordinates and the curve.

**Random state and samplers**

- `ZZ.random_element` burns the unconditional `den` draw (`integer_ring.pyx:801`), so every
  seeded stream lines up with SageMath again; `distribution='mpz_rrandomb'` is now GMP's runs
  generator (`randstate.ts:random_bits_rrandomb`), not `mpz_urandomb`; the gaussian branch pins
  `algorithm='uniform+logtable'`.
- `RandState.ZZ_seed()` / `long_seed()` added; `python_random(seed)` no longer reseeds a cached
  generator.
- `discrete_gaussian_integer`'s `upper_bound` rounds like MPFR; `_maximal_r` no longer runs an
  extra power-iteration step; `_iter_vectors` added.

**Coding and crypto**

- `SBox`'s output size is `ZZ(max(S)).nbits()` (exact, and 0 when every output is 0);
  `is_involution` raises for a non-permutation.
- `ReedMullerCode.length()`/`minimum_distance()` return `bigint`.
- `GoppaCode.generator_matrix()` returns the echelon basis.
- `gen_lattice(type='random')` draws from `IntegerModRing.random_element` above the dense-template
  modulus bound, instead of bounding every entry by 2^31.

**Finite fields, groups and arithmetic**

- `GF(p)` reports `inverse of Mod(0, p) does not exist`; `IntegerMod` inversion, `isOne`,
  `multiplicative_order`, `multiplicative_generator` and `units()` all handle `Z/1Z` and match
  SageMath's exception classes; `cyclotomic_polynomial(n <= 0)` raises `ArithmeticError`;
  `factor(0)` raises `ArithmeticError`; `bsgs` and `discrete_log_rho` messages match; the
  polynomial quotient ring enumerates with the constant coefficient varying fastest; exhausted
  Cantor-Zassenhaus splitting raises `AssertionError`, as SageMath does.

### parigp-ts

- `Z_factor(0)` returns the matrix `0^1` as PARI's `ifactor` does
  (`basemath/ifactor1.c:4459-4463`), instead of throwing — which is what made
  `BinaryQF.solve_integer(0)` raise where SageMath returns None.

### Quadratic forms

- `BinaryQF.compose` reproduces PARI's `Qfb` domain validation (negative definite, square
  discriminant); `solve_integer` on the zero form returns None for `n != 0`.

### Errors

- `IndexError`, `AssertionError` and `PariError` added to `errors.ts`.

### Property-test framework

- `tests/property/typescript/mersenne-twister.ts` `getrandbits(k)` packs the 32-bit words
  LITTLE-endian and masks the LAST one, as CPython's `_random_getrandbits` does. It used to pack
  them big-endian and mask the first, which agrees only for `k <= 32`; above that the two runners
  would have generated DIFFERENT arguments, which `compare.ts` reports as a missing test rather
  than as a failure. Verified against CPython at k = 1, 8, 31, 32, 33, 64, 65, 128 and through
  `_randbelow` at seven widths up to 2^100.

### Docs

- `DEVIATIONS.md`: the "Isogeny Kernel-Polynomial Algorithms", "Real and Complex Printing" and
  "Power Series — truncate()" gaps are closed and rewritten; the false "BMSS and Stark return the
  same kernel polynomial" claim is retracted; accepted return-type and numeric-backend differences
  are separated from the remaining open fidelity gaps.
- `SCOPE.md`: new sections for `sage.schemes.hyperelliptic_curves`, `sage.algebras.quatalg`,
  `sage.rings.function_field` and the general quadratic forms modules.

## 0.0.14 - 2026-07-28

The **upstream-porting pass**: the last large PARI/SageMath modules the earlier passes had
listed as out of reach, plus their consumers.

### parigp-ts

- **MPQS** (`src/mpqs.ts`, port of `mpqs.c`): self-initialising multiple polynomial quadratic
  sieve with PARI's size-indexed parameter tables, full/large-prime relation stores and GF(2)
  elimination. Wired into `ifac_crack` in PARI's own position, so `Z_factor`'s chain is complete
  and the five call sites that used to throw on hard semiprimes now answer. The relation store is
  in memory rather than disk-backed; MPQS still declines above 107 decimal digits, as PARI does.
- **Modular and Hilbert class polynomials** (`src/polmodular.ts`, port of `polmodular.c`,
  `polclass.c`, `volcano.c`): `polmodular_ZM`/`polmodular_ZXX`/`Fp_polmodular_evalx`, the
  `polmodular_db_*` cache, the class-invariant predicates and `polclass0`. This replaces PARI's
  separately distributed `seadata` package, which is not vendored.
- **SEA** (`src/elliptic/ellsea.ts`, port of `ellsea.c`): `Fp_ellcard_SEA` with Elkies, Atkin,
  `match_and_sort`, the CM branch, `Fp_elljissupersingular` and `Fq_elldivpolmod`.
- **Class and unit groups of quadratic fields** (`src/buch.ts`, port of `buch1.c` plus
  `hnf_snf.c`, `Qfb.c`'s `qfr3`/`qfr5` and `alglin1.c`'s `ZM_pivots`): `Buchquad`,
  `quadclassunit0`, `quadclassno`, `bnfinit`, with PARI's `t_REAL` kernel for Shanks distances.
- **Galois groups** (`src/galconj.ts`, port of `galconj.c` with the `perm.c`, `Zp.c` and `FpX.c`
  support it needs): `galoisinit`, `galoisgen`, `galoispermtopol`, `galoisfixedfield`,
  `galoissubgroups`, `galoisconj4`.
- **Theta series** (`src/qfrep.ts`, port of `bibli1.c`'s `qfrep0`/`minim0_dolll` with
  `lllgramint` and `qfgaussred_positive`).
- **`qfb.ts` gained PARI's `t_REAL` kernel and the Shanks-distance forms**: `QfbExt`, the
  `qfr3`/`qfr5` containers and a second overload on `qfbred`/`qfbcomp`/`qfbcompraw`/`qfbsqr`/
  `qfbsqrraw`/`qfbpow`/`qfbpowraw`, so `flag | qf_NOD` works.
- **Fixed `common_nbr` (`volcano.c:407-427`)**: a double root of the degree-2 gcd was reported as
  *two* candidates rather than one, so `surface_parallel_path` took the ambiguous branch and
  failed unconditionally when `n[0] == 2`. Every j-invariant `polclass0` drew was rejected and
  the routine never terminated for non-fundamental discriminants such as `D = -288`. Regression
  tests added for eight non-fundamental discriminants against PARI's `polclass(D)`.
- **Restored `polclass_roots_modp`'s `endo_cert` handling** (`polclass.c:1748-1777`): the port had
  dropped PARI's `if (!res && endo_cert) pari_err_BUG(...)` and its `vecsmall_isin_skip` repeat
  test, turning a diagnosable bug into an infinite loop.
- **`src/index.ts` now exports** `mpqs`, the `polmodular`/`polclass` surface, `Fp_ellcard_SEA`,
  `qfrep0`/`qfrep`, the `buch` class-group surface, the `galconj` surface and the extended `qfb`
  surface. `buch.ts` carries a second, independent copy of PARI's `t_REAL` kernel; only `qfb.ts`'s
  is re-exported from the package root, and the clashing `buch` names are omitted with a note.

### sagemath-ts

- **Laurent series** (`src/rings/laurent_series_ring.ts`): `LaurentSeriesRing` /
  `LaurentSeriesElement` with full arithmetic, `__getitem__`, `list`, `coefficients`,
  `exponents`. `PowerSeriesElement.__call__`, `.inv()` and the constructor now follow SageMath's
  precision rules; `MPowerSeriesRing` / `MPowerSeries` added.
- **Polynomial matrices** (`src/matrix/matrix_polynomial_dense.ts`): shifted reduced / weak Popov
  / Popov / Hermite forms and minimal approximant bases. Exported from `matrix/index.ts` with the
  generic names (`degree`, `truncate`, `shift`, `reverse`, `hermite_form`) aliased.
- **van Hoeij / LLL recombination** in `rings/polynomial/polynomial_element.ts`, replacing the
  200 000-subset Zassenhaus budget.
- **Number field embeddings** (`src/rings/number_field/number_field_embeddings.ts`);
  `embeddings` / `real_embeddings` / `complex_embeddings` / `places` are implemented rather than
  throwing, and `galois_group.ts` delegates to parigp-ts's `galoisinit`/`galoisfixedfield`.
- **`discrete_gaussian_lattice.ts` delegates its theta series to parigp-ts `qfrep0`** instead of a
  local floating-point Fincke-Pohst enumeration, and `_normalisation_factor_zz` returns a
  multiprecision `RealNumberMP` (exported from `stats/distributions/index.ts`).
- **`matrix_operations.norm(A, 2)` now follows `matrix2.pyx:16460-16471`** (`change_ring(CDF)` +
  SVD) and accepts RR/CC entries; `jordan_form` honours `subdivide`.
- **`elliptic_curves/formal_group.ts`** works over a genuine Laurent series ring:
  `mult_by_n`'s characteristic-zero branch (`formal_group.py:644-665`) is ported line for line
  and `group_law` computes in `MPowerSeriesRing`, so Sage's three-variable associativity TESTS
  block over `GF(7)[[x,y,z]]` is verified. `isogeny_class.ts`'s `Frobenius_filter` works over an
  arbitrary number field and is corrected against `isogeny_class.py:1202-1203`.
- **`number_field.ts`** proves `class_number() == 1` by exhibiting a generator for every
  Minkowski-bound prime (Sage's `[1, 1, 1]` Hecke doctest, `ZZ[2^(1/3)]`, `Q(sqrt2, sqrt3)`);
  `nfrootsof1` proves the number of roots of unity, so `unit_group()` no longer claimed torsion
  order 2 for every degree > 2 field. `NumberFieldIdeal.valuation` works at any residue degree.
- **`modules/free_module.ts`**: `intersection()` over `K[x]` returned a basis off by a unit of
  `K[x]`. The cause was `Matrix.integer_kernel`'s missing `self.denominator()` scaling — over
  `QQ[x]` that is the lcm of the rational *coefficient* denominators, and the port only cleared
  fraction-function denominators, so it never scaled. 914 cases now match SageMath exactly, where
  22 of 250 intersections were wrong.

### Still throws (unchanged, or newly narrowed — with the upstream routine named)

- **`bnfinit` for degree > 2**, and therefore the class group structure for `h > 1`,
  `regulator()` and `fundamental_units()` over a number field. `Buchall_param`
  (`buch2.c:3946`) needs `nfinit_basic`/`nfmaxord`/`idealprimedec`/ideal HNF arithmetic/the `T2`
  form/`nfrootsof1` — ~20k lines of `base1.c`–`base5.c` that `parigp-ts` does not have — and the
  algorithm is circular (a rigorous principality test needs the units, which come out of the same
  relation search). `buch.ts` implements degrees 0, 1 and 2 and throws above.
- **Galois group of a non-Galois field**: both routes (`splitting_field.py:371`, PARI
  `nfsplitting0` at `base1.c:1413`) need `nffactor` + `rnfequation`.
- **`find_isogenous_from_Atkin` / `find_isogenous_from_canonical`** (`ellsea.c:900`, `:964`) —
  the only routines here that genuinely need `seadata`, and PARI itself never reaches them
  without it. **`Fq_ellcard_SEA` with `T != NULL`** (extension fields) is absent.
- **Weber / double-eta / Atkin class invariants** in `polmodular`'s CM path (need ~1500 lines of
  double-eta tables plus `polclass.c`'s orientation machinery). SEA never asks for them.
  `polmodular0_powerup_ZM` is transcribed but unreachable, hence untested.
- **`s4galoisgen` / `f36galoisgen`** (`galconj.c:1519`, `:1698`) need `FpX_ffisom`/
  `FpXV_chinese`. `findpsi` and `galoisgenlift_nilp` are transcribed but unreachable from any
  tested input, and are flagged as untested rather than claimed.
- **`mpqs_class_init` / `mpqs_class_rels`** (`mpqs.c:1775`, `:1815`) are absent; their only
  caller is `buch2.c`. `buch.ts` uses upstream's own `imag_relations` fallback, so results are
  identical and only slower. The shared `MPQS_MODE_CLASSGROUP` branches are untested.
- **`MPowerSeries` division by a non-unit** (needs `quo_rem`, `# needs sage.libs.singular`);
  Laurent ordering comparisons; `polredbest` in `fixed_field`; `precision='dp'` in the integer
  Gaussian sampler.

### Known gaps left open deliberately

- `ellcard` (`parigp-ts/src/elliptic/group.ts:1318`, `:1357`) still routes >= 96 bits to base
  Schoof rather than the new `Fp_ellcard_SEA`. Correct, and orders of magnitude slower than it
  needs to be; switching it is a behaviour/performance change that needs its own test pass.
- `MPowerSeries.inv()` in `power_series_ring.ts` does not match upstream's precision (upstream
  inverts the background univariate series, `multi_power_series_ring_element.py:725`).
  `formal_group.ts` works around it locally: 12.1 s -> 0.2 s on `group_law(50)`.
- `Matrix.toString` (`matrix_generic.ts:429`) is not subdivision-aware and pads per column rather
  than to Sage's single global width; `jordan_form` attaches the new `matrix_str` per instance.
- `Polynomial.roots()` does not use Sage's `Factorization` order (multiplicity ascending, then
  root descending), so `jordan_form`'s block order differs for multi-eigenvalue matrices.
- `class_group.ts` and `cm.ts` still enumerate reduced forms behind
  `CLASS_GROUP_DISC_BOUND = 2e6` instead of delegating to the new `Buchquad`.
- `buch.ts` and `qfb.ts` each carry an independent transcription of PARI's `t_REAL` kernel (24
  clashing names). Only `qfb.ts`'s is exported at package level; they should be merged.

### Documentation

- `DEVIATIONS.md`: new section 64, *Newly Ported Upstream Modules (0.0.14)*, with the residual
  deviations of each ported module. Ten sections were rewritten and several rows **deleted**
  rather than softened, because the deviation existed only while a dependency was missing —
  MPQS, SEA/`seadata`, `polmodular`, `qfrep`, `galoisinit`, quadratic `bnfinit`, Laurent series,
  Popov/approximant bases and van Hoeij are all ported now. Table of contents synced.
- `SCOPE.md`: rows added for the six new parigp-ts modules and the three new sagemath-ts modules;
  the stale "remaining" notes on `polynomial_element`, `formal_group`, `discrete_gaussian_lattice`,
  `isogeny_class`, `free_module`, `class_group`, `unit_group`, `galois_group`,
  `matrix_operations` and `matrix_decompositions` corrected.

### Tests

- The `blockedByPolclass288` skip in `parigp-ts/src/elliptic/ellsea.test.ts` was removed now that
  `polclass0(-288)` works; the exhaustive `Fp_ellcard_SEA` sweeps no longer exclude any curve.
- `isogeny_class.test.ts`'s "candidate set for d = -23" pinned `[2, 3, 5]`; SageMath's doctest
  (`isogeny_class.py:1202-1203`) says `[2, 3]`. Corrected to build the doctest's curve and assert
  SageMath's value.
- `number_field.test.ts`'s "throws rather than guessing when the criterion is inconclusive"
  asserted that `Q(2^(1/3)).class_number()` throws. Sage's `order.py:1181` gives 1, and the new
  principality certificate proves it, so the assertion pinned a limitation rather than a value.
  Rewritten to expect `1n`, with the "must still throw" half moved onto `x^3-19` (`h = 3`).
- `ifactor.test.ts`'s "reports failure instead of declaring a composite prime" asserted a throw
  because MPQS was missing. MPQS now factors that number, so the test was pinning wrong
  behaviour. Split into the exact factorization through the MPQS stage plus the same throw-path
  test driven with `mpqsMaxPolys: 1`, so both paths stay covered.
- `galois_group.test.ts` pinned the port's own wording `"Prime 2 is ramified"`; Sage's
  `galois_group.py:767` names the ideal. Corrected to Sage's message.
- Eight non-fundamental discriminants added to `polmodular.test.ts` against PARI's
  `polclass(D)`. Every pre-existing `polclass0` case used a *fundamental* discriminant, which is
  why the `common_nbr` bug survived.
- No test was deleted, no assertion weakened, no tolerance loosened, and no skip added. Final
  state: `bun test` 6781 pass / 32 skip / 0 fail across 115 files (the 32 skips are the
  pre-existing ones); `bun run test:property` 433/433.

## 0.0.13 - 2026-07-28
- `DiscreteGaussianDistributionIntegerSampler` accepts Sage's `precision` keyword: `'mp'` (default) works, `'dp'` throws naming the unported `dgs_gauss_dp.c`, and any other value raises Sage's exact `ValueError("Parameter precision '...' not supported")`. Previously the keyword was absent, so an unsupported precision was silently ignored.

## 0.0.12 - 2026-07-28

The **deferred-work pass**: the items 0.0.11 left as honest `SAGE_NOT_IMPLEMENTED` stubs, taken
across 12 work units. The theme is that most of those stubs existed because a *dependency* was
missing, not because the algorithm was hard — so the dependency was ported and the stub deleted.

### Dependencies that stopped being stubs

- **`parigp-ts` `ifactor.ts` is now PARI's real factoring chain** (746 -> 1609 lines): `tridiv_bound`
  + gcd-with-primorial trial division, then `ifac_crack`'s order — pure powers, SQUFOF,
  Pollard-Brent rho, Lenstra-Montgomery ECM — driven by an `ifac_decomp` worklist, plus
  `ispower.c`'s perfect-power machinery and a real `isprimepower` that never factors its
  argument. Verified against brute force exhaustively on 1..5000 and 2..20000, 2000 random
  semiprimes, and the published factorizations of `F6`, `M67`, `M71`, `M101` and `F7 = 2^128+1`
  (the last found by the ported ECM in 16 s).
- **`parigp-ts` gained `ffinit`, `matkermod` and the `Qfb` family.** `ffinit` reproduces PARI
  coefficient-for-coefficient for all of the first 60 primes × n ∈ [2,12] (660/660);
  `matkermod` reproduces all 24 golden values decoded from PARI's own regression suite; `qfb`
  reproduces real PARI 2.15.4 on ~2500 golden values.
- **`misc/randstate.ts` now matches SageMath's seeded streams exactly.** GMP is not vendored, so
  GMP 6.3.0's sources were obtained and `randmts.c` (`mangle_seed`, `randseed_mt`), `randmt.c`
  (including the `default_state` buffer that seed 0 lands on) and `urandomm.c` were ported
  verbatim; CPython's `random.Random` was ported as `PythonRandom`. Checked against a C oracle
  linked to libgmp and against SageMath 10.3's own doctest values.
- **`dgs_bern.c` ported**, so all four discrete-Gaussian integer algorithms work — and all four
  reproduce SageMath's *seeded sample streams* bit-for-bit, not merely its distributions.

### Now implemented (previously threw)

- `factor()` over ZZ/QQ: a **real Zassenhaus** (DDF + Cantor-Zassenhaus, multifactor Hensel
  lifting, exact Landau-Mignotte bound, subset recombination) replacing code that peeled off
  integer roots for degree <= 10 and returned the rest as one "irreducible" factor. `QQ.__call__`
  was rewritten as a port of `Rational.__set_value`. Together these unblock `minpoly` over QQ.
- `matrix_integer.frobenius_form(2)`, via a verbatim port of PARI's `RgM_Frobenius`.
- `matrix_modn.right_kernel_matrix` for composite moduli, via `matkermod`.
- `matrix_decompositions.jordan_form(transformation=true)` — reproducing Sage's *exact* `P` —
  and `krylov_kernel_basis(variable=…)`.
- `matrix_operations.norm(A, 2)`, `is_similar(transformation=true)`, `is_diagonalizable(base_field)`
  and a generic `change_ring`.
- `finite_field_extension.irreducible_element` delegating to NTL `BuildSparseIrred` and PARI
  `ffinit`; `polynomial_gf2x.ts` delegating its whole arithmetic layer to `ntl-ts`.
- Number fields: real quadratic **fundamental units** and regulators (`quadunit`),
  `nfgaloisconj` at **any degree** (the degree-8 cap is gone), and `decomposition(p)` at
  **inessential discriminant divisors** (Buchmann-Lenstra round 4).
- Free modules over **non-ZZ PIDs** (`QQ[x]`, `GF(p)[x]`), real quotient `lift`/`project`, and a
  real embedded `tensor_product`.
- Elliptic curves: the **full 13-discriminant `Fp_ellcard_CM` table**, `bernardi_sigma_function`,
  supersingular `alpha` (with a new ramified quadratic p-adic extension), `Frobenius_filter` over
  QQ, and `BinaryQF.solve_integer`.
- Non-spherical Σ for the lattice Gaussian sampler (Peikert's `r`, Cholesky, offline samples,
  `_normalisation_factor_zz` with a local `qfrep`).

### Landed only partially

- **`ellcard_sea` is Schoof, not SEA.** Elkies and Atkin both need the modular polynomials
  `Phi_l`, which PARI reads from the separately distributed `seadata` package —
  `reference/pari` ships the reader but `reference/pari/data` is empty. Schoof's base algorithm
  is ported in full and is exact (verified against PARI's own `ellsea` regression vectors at 65,
  70 and 101 bits), but it is `O(log^5 p)`, so `ellcard` keeps Shanks below `expi(p) = 96` — a
  **measured** crossover, not PARI's 56.
- **Class groups of degree > 2 fields** answer only in the provably-trivial Minkowski case.
  Fields whose true class number is 1 (`Q(2^(1/3))`, `Q(zeta_8)`, `Q(zeta_23)`) still throw,
  because we cannot *prove* it without `bnfinit`.
- **`Frobenius_filter` works only over QQ.** Sage's headline `d = -23` example lives over a
  degree-6 field, so `isogeny_degrees_cm` returns the unfiltered (still sufficient) superset
  there rather than a wrong answer.
- **`_normalisation_factor_zz` runs in double precision**, matching Sage's `prec=100` doctest to
  15 significant digits rather than 28.
- **The formal group needed no work** — the deferral note calling `differential()` hardcoded was
  stale as of 0.0.11. This pass proved it correct against Sage's doctests and defining
  identities, and added `x_list`/`y_list` because `x()`/`y()` returned objects whose coefficients
  no caller could read.

### Still genuinely deferred, with reasons

- **MPQS** (`mpqs.c`, ~2600 lines of sieving and GF(2) linear algebra) is the one missing stage
  of `Z_factor`. **This is a behavioural change, not just a gap:** 0.0.11 `console.warn`-ed and
  returned an unfactored composite *as if it were prime*; it now throws `NotImplementedError`
  naming `mpqs.c`. Five call sites that previously got a wrong answer now get an exception.
- **`seadata` modular polynomials** — see above.
- **`bnfinit`** and, equivalently, archimedean embeddings of number fields.
- **van Hoeij/LLL recombination** for polynomial factorization: the classical subset search
  raises after a 200 000-subset budget rather than returning a partial factorization. Nothing
  constructible reached it (Swinnerton-Dyer degree 32 finishes in 324 ms).
- **Modular symbols** (`sage.modular.modsym`), which still gate `padic_lseries`'s `series`,
  `measure`, `modular_symbol`, `order_of_vanishing`, `_c_bound` and the `Dp_valued_*` methods.
- **Shanks distance forms** (`qfr5_*`) — a `t_REAL` quantity needing an arbitrary-precision float
  kernel that CLAUDE.md forbids.
- **NTL's randomized routines** (`BuildRandomIrred`, `SquareFreeDecomp`, `DistinctDegFactor`,
  `EqualDegFactor`, `BerlekampFactor`), so `polynomial_gf2x.ts` keeps four local factoring
  routines and `algorithm='random'` uses rejection sampling — which is the same fallback
  SageMath takes when its own NTL import fails.
- **PARI `ffgen`/`ffprimroot`/`charpoly` over `F_q`**, so `algorithm='ffprimroot'` throws rather
  than returning a polynomial that is irreducible but silently not primitive.
- **`flint-ts` remains 100% stubs.**

### Bugs found and fixed en route

- **`packages/flint-ts/src/index.ts` could not be imported at all.** The barrel re-exported five
  *interfaces* through a value `export` clause; `tsc` elides them so it typechecked, but at
  runtime `import … from '@sagemath-ts/flint-ts'` threw `export 'nmod_t' not found`.
  Pre-existing since the first commit, found by a sweep that imports all 128 source modules.
  Fixed with `export type`; `flint-ts` gained its first test file.
- **`parigp-ts/src/qfb.ts` hung forever** in `qfbsolve` for indefinite forms with negative `n` —
  two independent fidelity bugs: `normforms` used `a/N` where PARI's `itou` takes `|a|/N`
  (`Qfb.c:1766`), and `Zn_quad_roots` kept the `-1` that PARI's `clean_Z_factor` drops, so
  `Z_pvalrem(D, -1n)` looped. Validated against a live PARI oracle on 220 random forms:
  0 disagreements.
- **`crypto/lattice.ts` `gen_lattice` used the wrong generator** — Sage's modular block goes
  through `rstate.c_random() % q` row-major, not `mpz_urandomm`. All three doctests now
  reproduce exactly.
- **`isogeny_degrees_cm` was unsound**, not merely non-minimal: the horizontal (class-group)
  primes step was missing entirely, so the returned list could be *too small*.
- `errors.ts` gained the `RuntimeError` class SageMath uses.
- Several barrel gaps closed: `matrix/index.ts` (`change_ring`, `pivots`),
  `schemes/elliptic_curves/index.ts` (`Frobenius_filter` and the two p-adic extension classes),
  `modules/index.ts` (five module classes).

### Housekeeping

- `polynomial_commitment.ts` **moved from `rings/polynomial/` to `src/zk/`** — it has no
  SageMath counterpart, so it no longer sits inside the mirrored Sage tree.
  `rings/polynomial/index.ts` keeps a byte-identical compatibility re-export block, verified by
  object identity across all five import surfaces. `playground/docs-data.json` regenerated.
- `DEVIATIONS.md` is now 64 sections: 2 new (**Matrix Special Constructors**, closing audit item
  L44, and **Bounded Search Budgets and Measured Thresholds**) and 18 rewritten — most of them
  because the deviation they described no longer exists.
- Tests: **6216 pass, 32 skip, 0 fail**, 2 738 804 expect() calls across 106 files (up from
  6208/33/0 across 105). Property transcripts 433/433. One `test.skip` with an empty body
  (`ellcard_sea`) was un-skipped and given three PARI golden values. Across the whole 12-unit
  diff, 316 tests were added and 9 removed; all 9 removals are renames where a test pinning the
  port's own wrong or stubbed value was replaced by golden SageMath values.
- Typecheck: `flint-ts` and `ntl-ts` 0 errors; `parigp-ts` 29 errors **byte-identical to the
  HEAD baseline** (all pre-existing, all in test files). `packages/sagemath-ts` still has no
  `tsconfig.json`, so it cannot be typechecked; under an ad-hoc strict config its error set is
  unchanged by this pass.

## 0.0.11 - 2026-07-28

Adversarial audit of the whole port against the vendored upstream (`AUDIT-2026-07.md`, 370
confirmed findings) and the fix pass that followed.

**355 of the 370 confirmed findings fixed outright** across 23 work units (audit severities:
17 Critical, 120 High, 153 Medium, 80 Low). All 17 Criticals were addressed; C8 (`minpoly` over
QQ) and C11 (`right_kernel_matrix` over composite Z/nZ) retain a narrowed, documented residual
gap. The doc-only findings (M38, M108, M121, M142, M146, H119, L31, L44, L60) are closed by this
commit. Genuinely still open: H100 (`bernardi_sigma_function`, needs an exact formal-group log)
and the delegation half of H120 (`getDefaultModulus` -> NTL `BuildSparseIrred` / PARI `ffinit`).

Headline defects:
- **RNG period** — `misc/randstate.ts` took the *low* bits of a 64-bit MMIX LCG, so bit *k* had period 2^(k+1) and every consumer emitted short deterministic cycles (`random_below(2n)` alternated `1,0,1,0,…`). Replaced with MT19937, the generator family GMP's `gmp_randinit_default` uses, with GMP-compatible per-call bit consumption.
- **Reducible Conway entries** — 7 entries in the Conway polynomial table were reducible (so `GF(29^2)` was not a field) and a `GF(2^128)` pentanomial was fabricated. The table is now regenerated by porting FLINT's `conway.c` decoder; every entry is verified irreducible, primitive, normalised and subfield-compatible.
- **p-adic addition** — `padic_generic_element.add()` multiplied by `p^v` twice for operands of unequal valuation, so most p-adic sums were wrong.
- **Vélu y-coordinate** — a sign error in `ell_curve_isogeny.ts` put isogeny images off the codomain.
- **minpoly** — `matrix_operations.minpoly` returned the minimal polynomial of `e_0` rather than of the matrix.
- **PARI `ellcard`** — returned wrong cardinalities at primes as small as p ≈ 100, and `ellgroup` was wrong in 85 of 476 brute-force-checked cases. Both now port PARI's real `Fp_ellcard_Shanks` / `gen_ellgroup` and are verified against exhaustive oracles and real PARI 2.15.4.

Also in this release:
- Exact arithmetic replaces floating point in LLL, free-module `coordinates`/`echelonize`/`discriminant`, lattice Gram-Schmidt, CVP/SVP enumeration, binary quadratic form reduction and `real_log`.
- `ntl-ts` gains a real `GF2`/`GF2X` (`IterIrredTest`, `BuildIrred`, `BuildSparseIrred` over a vendored `GF2X_irred_tab`), replacing 100% stubs.
- A new number-field kernel (`rings/number_field/pari_nf.ts`) ports `nfbasis`/`nfdisc`, `idealprimedec`, `nfgaloisconj` and `polisirreducible`, giving real maximal orders, prime decomposition and quadratic class groups.
- Paths that cannot match SageMath now raise `SAGE_NOT_IMPLEMENTED` instead of returning placeholders or silently wrong answers.
- Test suite: 5782 pass / 33 skip / 0 fail (1 253 493 expect() calls, 102 files); property transcripts 433/433. 22 modules that had **no test file** — including `randstate`, `conway_polynomials`, `cm` and `formal_group` — now have one; ~20 tests that pinned the port's own wrong values were corrected to SageMath's.
- `DEVIATIONS.md` grew to 62 sections: 18 new, and 16 existing entries corrected (several claimed functions throw when they did not, or claimed parity that did not hold).
- `SCOPE.md` status markers reconciled with the audit; percentages now mean "ported **and** verified", so several were revised downward without any regression.

## 0.0.10 - 2026-07-27
- Fixed `FiniteFieldExtension.isIrreducible` accepting reducible polynomials that split completely (`x^{p^k} = x mod f`), which produced non-field `GF(p^n)` for p outside the Conway table.
- Added regression tests asserting default moduli past the Conway table are irreducible.

## 0.0.9 - 2026-02-04
- Audited Gröbner basis implementation and documented simplified algorithm deviation.
- Linked multivariate ideal Groebner docstrings to deviations.

## 0.0.8 - 2026-02-04
- Completed `sage.groups.generic` audit and documented API/algorithm limitations.
- Aligned Pohlig-Hellman digit solving with Sage behavior.

## 0.0.7 - 2026-02-04
- Completed polynomial audit and documented roots/factorization and ideal-dimension deviations.
- Added deviation links in polynomial and multivariate ideal docstrings.

## 0.0.6 - 2026-02-04
- Completed finite_rings audit and documented Conway/minimal polynomial deviations.
- Refreshed DEVIATIONS.md table of contents.

## 0.0.5 - 2026-02-04
- Completed elliptic curve audit and documented isogeny/torsion deviations and PARI linkage notes.

## 0.0.4 - 2026-02-04
- Completed `integer_ring` audit and documented deviation links and cleanup.

## 0.0.3 - 2026-02-04
- Completed `sage.arith` audit and documented Bernoulli and Dedekind sum deviations.
- Added deviation references in `arith/misc.ts`.

## 0.0.2 - 2026-02-04
- Audited parigp-ts algorithms and documented deviations for factorization and elliptic curve advanced algorithms.
- Fixed PARI Fp_ellpoint mapping to match the reference implementation.
- Documented arith/misc deviations (algebraic_dependency approximation, gauss_sum numeric-only, hilbert_symbol direct-only) and linked docstrings.
- Recorded audit status for flint-ts and ntl-ts.

## 0.0.1 - 2026-01-30
- Initial public snapshot.
