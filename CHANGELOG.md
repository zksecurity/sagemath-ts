# Changelog

All notable changes to this project will be documented in this file.

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
