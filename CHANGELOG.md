# Changelog

All notable changes to this project will be documented in this file.

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
