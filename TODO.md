# Audit handoff — updated 2026-10-03

Latest elliptic checkpoint (28.12.0): lift_isogeny, bilinear/power-derivative
helpers, getc2 and all four Kohel coefficient tables are ported. Comparisons
check complete Newton evaluation/correction states for synthetic and Kohel
polynomials, in addition to exact results and native errors. 25,079 checks,
405 docs/storage tests and eight builds pass; all 568 complete baseline
diagnostics are unchanged.
Next: inverse-Frobenius precomputation (Flxq_lroot_pre and its automorphism
powering dependencies), get_norm's small extension-polynomial power, and
zx_is_pcyc; then connect Flxq_ellcard_Kohel. The ternary branch uses a4=[a2]
(native vector tag), not a short-model a4. Harley also needs modular-polynomial
construction and get_trace_Robert. Satoh/Kedlaya, extension SEA, subfield
descent, general FF dispatch and group integration remain open.
Case definitions total 999,756 bytes; compact repeated constant recipes while
preserving their inputs before adding more controls. No generated data belongs
in the repository.

Previous checkpoint (28.11.0): integer cvtop, Qp_exp precision records,
ZpXQ_sqrtnorm, Flxq_lroot_fast_pre and Teichmuller_lift are ported. Norms retain
relative precision e-1+valuation(s); element lifts use native Newton/Dixon
linearized Frobenius solves. Full GEN p-adic arithmetic and noninteger cvtop
remain unimplemented outside the standalone adapter contract.

Previous checkpoint (28.10.0): scalar inverse/division, Zp_exp and ZpXQ_log
were ported with native digit/binary splitting and atanh evaluation. Native
precision-one quotient behavior is preserved. 15,032 comparisons passed.

Previous checkpoint (28.9.0): scalar root lifts, Frobenius and cyclotomic norms
are ported. Zp_sqrt preserves native inverse exceptions and -2 binary dispatch.
The binary n=1 native zero-exponent bug is corrected with a guarded original-body
regression and independent odd-residue comparisons.

Previous checkpoint (28.8.0): native polynomial Newton/Dixon drivers,
quotient inverse/division and Flx_Teichmuller are ported. The four-argument
FpXQ_inv now delegates to the original precision schedule; comparative traces
retain the p=3,e=5 regression. 6,011 comparisons pass.
Next p-adic dependencies: ZpXQ Frobenius, cyclotomic/general square-root norms,
logarithm and scalar lift helpers; then isogeny/Teichmuller element lifting and
Kohel/Harley. Satoh special cases, Kedlaya, extension SEA, subfield descent,
general FF dispatch and group/default-order integration remain open.

Previous checkpoint (28.7.0): native extension-field Shanks–Mestre
counting and order search are ported. 17,157 native comparisons pass, including
1,002 Shanks cases and order-search traces; the generated GF(49) j=1728 input
that stalls direct native Shanks is retained under its correct constant-j backend.
Next: small-characteristic/p-adic backends, subfield descent and general field
dispatch; then FF/Sage cardinality, group and default-order integration.

Previous checkpoint (28.6.0): native Flxq/FpXQ constant-j twist-count
helpers are ported and compared, including j=0, j=1728 and other base-field j.
Next: extension Shanks, remaining small-characteristic/p-adic backends and
subfield descent/dispatcher integration.

Previous checkpoint (28.5.0): the ternary supersingular counting
formula is ported with a documented correction for a bundled PARI bug: an
odd-degree random root can cause the native code to return the twist count.
The GF(27) 19/37 seed regression and independent enumeration are retained.
Ordinary ternary and general field cardinality backends remain open.

Previous checkpoint (28.4.0): native direct Flxq/Fl2 square-root
branches and binary automorphism/fast-root dependencies are ported and compared
for exact roots and RNG fingerprints. These support ternary special counting;
field cardinality dispatch and the other counting backends remain open.

Previous checkpoint (28.3.0): native word-field point kernels and
Fl_ellcard_Shanks are ported; default prime counting now selects word Shanks and
SEA in the native ranges. Raw point, count and RNG-fingerprint comparisons pass.
Exact search bounds intentionally replace native floating approximations.

Previous checkpoint (28.2.1): default prime counting selects SEA at the native
threshold after CM; comparative branch traces and actual large-prime counts pass.
The machine-word Shanks kernel remains to port.

Previous checkpoint (28.2.0): native elltrace_extension and Fp_ffellcard
helpers are ported and compared, including logarithmic degree powering. They
support the base-field-model branch of the future FF cardinality dispatcher.

Cardinality dependency order for the next pass:
1. Prime counter dispatch and word Shanks are implemented (28.3.0); continue
   field-counting dependencies below. Composite inputs remain outside the kernel
   contract; exact-bound adaptation is documented in DEVIATIONS.md.
2. Wire base-field-model extension counts through Fp_ffellcard (FpE.c:2187,
   FlxqE.c:1457). Native j/minimal-polynomial descent has separate branches.
3. Port the binary and ternary special-count branches, preserving the native
   square-root/root-selection dependencies (Flxq_sqrt differs from FF_issquareall;
   the native direct-root dependency is now available in 28.4.0).
4. Complete Harley/Kohel/Satoh/Kedlaya/Shanks/extension-SEA dependencies where the
   native dispatcher selects them, then FF_ellcard and general Sage callers.
5. Complete native group/exponent/generator dependencies and default FF point
   order. Supplied-bound order kernels are already available.
Remaining coercion, torsion and full module inventory still follow that work.

Previous checkpoint (28.1.0): gen_order and F2xqE/FlxqE/FpXQE_order are
ported and compared, including native callback schedules. FF_ellorder now handles
a supplied annihilating multiple/factorization. Prime order uses the same native
recursion and errors. Default field order still needs cardinality/group exponent;
next is native field cardinality dispatch, then group/default caller integration.

Previous checkpoint (28.0.0): integer-domain ellinit for p=2/3 now uses
native finite-field records (null for singular models); optimized ternary scalar
multiplication uses the same FF backend. Field accessors preserve their types,
and nonsingularity tests field-valued zero correctly. See AUDIT-2026-10.md.
Next: finite-field cardinality/group backends; remaining coercion/API/torsion/module inventory. The older singular
ellinit exception outside p=2/3 also needs migration to native empty-result behavior.

Previous checkpoint (27.0.1): ellinit_Fq accepts j-invariant and short
coefficient tuples, with live binary/ternary/large-field comparisons.

Previous checkpoint (27.0.0): native odd-extension/characteristic-three
kernels, FF model initialization and general Sage scalar callers are integrated
and compared, including the ordinary ternary x² term and 127-bit prime fields.
FFEllipticCurve now has binary and odd variants. Next: legacy integer-only
field-record API boundaries,
finite-extension cardinality/order/group backends, and the remaining full audit.

Previous checkpoint (26.0.0): binary model conversion, field initialization and
Sage scalar callers use native PARI kernels and preserve original coordinates.

Previous checkpoint (25.5.0): optimized scalar multiplication now coerces
through ZZ, delegates to PARI ellmul and propagates known orders. The short-record
PARI adapter uses native prime-model conversion for p > 3.

Previous checkpoint (25.4.0): optimized prime-field points now support
Sage's order algorithms/additive alias and share the curve cardinality cache.

Previous checkpoint (25.3.0): hybrid point order now follows the native
bounded-search/partial-factorization schedule. Default general prime-curve
cardinality/order and bounded Integer.factor are integrated and compared.

Previous checkpoint (25.2.0): general prime-field point order now uses
PARI and caches curve cardinality; native general-model ellcard/ellorder are
compared. Additional cardinality algorithms/options, groups and other parent backends remain open.

Previous checkpoint (25.1.0): general prime-field models now delegate scalar
multiplication to PARI, preserving coordinates, integer coercion and known orders.
Point order/cardinality/group backends and small-characteristic/extension scalar
backends remain open. Constructor model/coefficient fixes are in 25.0.0.
Typecheck now has 568 existing diagnostics (two constructor errors removed).


The user resumed the behavioral audit on 2026-09-18 after repository cleanup.
Completed: is_x_coord input conversion, GF(2) constructor parity, and lift_x
canonical coercion/promotion over the tested scalar fields, and lift_x extend=True
over QQ, prime fields and explicit finite extensions.
Optimized prime-curve coordinate callers now share that implementation (24.59.0).
Hyperelliptic root selection and scalar dispatch are also repaired (24.59.1).
Extension trace now uses the native dependency path (24.60.0).
Extension norm/charpoly now follows Sage/PARI, including native resultants (24.61.0).
Next: remaining coordinate parents, torsion and irreducibility. The September 9–14 commits
were squashed to reduce Git storage. Their messages
and former IDs are preserved in `tests/audit/squashed-commits-2026-09.txt`; old
IDs in audit notes are historical labels and may no longer resolve in Git.
Completed implementation code is preserved; the now-resolved partial patch and
its bulk coordinate research records have been retired. The user authorized
retiring unclassified bulk input sweeps on September 14; historical case counts below
describe the earlier audit checkpoints, not the current randomized suite.

The last verified audit commit before the squash was `e34ae87` (24.54.1); the earlier
square-predicate batch is `5e4f8d0` (24.54.0).

## Randomized input migration

The user explicitly chose fresh generators and accepted losing unclassified old
inputs. Area tests and native tests now generate from small constrained domains,
compare live, and print replay seeds. Only named regression inputs and focused
hand-written controls remain. Do not rebuild the old Cartesian corpora.

Constructor choices are currently deliberately bounded: valid field moduli,
packed native command layouts and matrix dimensions remain small finite domains.
Payloads, scalar values and native random seeds vary where their contracts permit.
Expanding these domains is future testing work; the migration does not preserve
all historical input coverage or prove full behavioral equivalence.

## Findings from the generator migration

Retained `replaySeeds` entries run these comparisons alongside fresh trials; no
bulk inputs or expected outputs are saved. These findings remain open.
Use generator `mt19937-recipes-v2` and seed 456 with the current case definitions:

- `--case lattices --runs 2`: two `lattice_volume` radical-format mismatches
  (`4*sqrt(14)` versus `sqrt(224)`, `9*sqrt(3)` versus `sqrt(243)`).
- `--case rationals --runs 2`: symbolic negative rational gamma is unimplemented;
  negative nonintegral factorial raises ValueError instead of TypeError.
- `--case quaternion_algebras --runs 2`: an indefinite lattice reduction exhausts
  the original PARI stack and the port's BigInt limit, producing different errors.
- `--case polynomials --runs 2`: the same constant-unit discrepancy in `poly_factor`.
- `--case polynomial_ops --runs 8`: `poly_factor_ff` includes a constant unit that
  Sage omits from the factor list (three mismatches). The zero-divisor generator
  issue was corrected by requiring a nonzero leading divisor coefficient.

Reproduce with `bun tests/property/compare.ts` plus the options above and
`--seed 456`. Check existing deviations/source contracts before repairing these
behaviors. Seed retention means affected area suites currently fail visibly;
this migration does not claim a green full behavioral suite.

The `function_fields` area also exceeded the normal 120-second Sage subprocess
deadline at seed 456 / runs 2 while compiling/executing native helpers. Its full
comparison remains unverified; investigate the slow helper before increasing the
limit (`SAGEMATH_TEST_TIMEOUT_MS`) or broadening that generator domain.

The broader ec_advanced run at seed 583475302 / runs 25 hit the 120-second
Sage deadline during this batch. A bounded native diagnostic progressed through
1,561 rows in 50 seconds, reaching polynomial-root isomorphism checks; this does
not identify a single hanging input. Large runs remain a performance/validation
gap, not evidence that their mathematical comparisons pass. Replay with the
same seed/generator and investigate cumulative native-helper overhead.

## What is verified

- September 29 norm/charpoly batch (24.61.0): 2,420 live new/dependency trace
  comparisons and 570 caller/docs/storage checks pass, as do eight builds. Twelve
  compact regression/control inputs accompany fresh generators. Typechecking has
  the same 570 diagnostics; 69 existing member-list messages gain the added
  charpoly member, with no other diagnostic changes. The compressed baseline is
  refreshed only after verifying that exact difference.

- September 29 trace batch (24.60.0): 908 live dependency/caller comparisons, 538
  existing caller/docs/storage tests and eight builds pass. Eight small inputs
  cover routing and kernel edge cases; no bulk test data. All 570 complete baseline
  type diagnostics remain unchanged.

- September 28 optimized caller/root batches (24.59.0–24.59.1): 1,854 coordinate,
  103 polynomial and 207 hyperelliptic comparisons pass. Nineteen small regression
  inputs accompany fresh generators. All 95 optimized curve tests, 121 hyperelliptic
  tests, 375 executed-docs/storage checks and final 41 case/storage checks pass.
  Eight builds pass; all 570 complete baseline diagnostic messages remain unchanged.

- September 28 extension batch (24.58.0): 1,242 fresh coordinate comparisons and
  103 named-polynomial dependency comparisons pass, including exact parent and
  root displays, point equations, doubling and negation. Sixteen small regressions
  accompany fresh generators; no bulk data was added. Existing caller/dependency
  tests: 245 pass; executed docs/case-format/storage: 373 pass. Eight builds pass;
  all 570 baseline type diagnostics remain unchanged after normalizing positions.
- September 18 lifting batch (24.57.2): 3,500 historical lift comparisons,
  929 fresh focused and 307 smaller advanced-area comparisons pass, along with
  298 caller tests and 374 focused/docs/storage checks. Eight builds pass; 14
  existing type diagnostics are removed with no additions (570 remain).
- September 18 coercion batch (24.57.1): 1,750 historical predicate comparisons,
  309 focused fresh comparisons, 221 fresh advanced-area comparisons, and 740
  focused/caller/docs/storage tests pass. Eight builds pass; full diagnostic
  messages match the 584-entry baseline after normalizing shifted positions.
- Generic curve coordinate roots, Montgomery selection and coefficient formatting:
  5,822 new shared native comparisons; all 10,202 advanced-area rows match.
  Broader callers/docs/case-format: 21,571 pass, 13 existing skips.
- Extension square roots: 4,324 comparisons, including exact PARI root choice/state.
- Extension square predicates and native resultant/norm dependencies: 4,832
  comparisons; all 37,802 field-area rows matched at that checkpoint.
- Earlier completed polynomial, modular-ring, isomorphism, PARI and NTL batches
  are recorded in AUDIT-2026-09.md and SCOPE.md.
- Eight package builds pass. There are 570 remaining TypeScript diagnostics (14 removed by base-change fixes);
  compare full messages/continuations, not just the count.
- Recorded source-line execution: 111541/133728 (83.41%). This is NOT full
  behavioral or input-domain coverage. AUDIT-COVERAGE.md inventories the gaps.

## Coordinate coercion resumed on September 18

The saved `is_x_coord` patch is now applied and completed; do not apply it again.
The predicate always calls K(x) and replaces only conversion TypeError, matching
ell_generic.py:710. The specialized GF2 constructor now accepts integer strings,
null/default input and matches nonintegral-float errors through the shared prime
field constructor. Unrecognized plain-object integer errors use the native type
name. Nine small comparative regressions cover these repairs; fresh inputs include
QQ, prime/extension fields, wrapped integers, strings, null, booleans and invalid
scalar values. No new bulk input or output records were committed.

All 1,750 historical is_x_coord and all 3,500 historical lift_x cases now match
live original execution. The expanded fresh profile passes 929 comparisons,
including whole-coefficient constructor traces, field promotion and parent
identity. Twenty new small regressions cover lift/base-change behavior; the nine
predicate/GF2 regressions remain. The completed historical patch and compressed
bulk records have been removed. No new generated corpus was committed.

The canonical-map subset covers ZZ/QQ, prime/explicit extension fields and
quotient-ring inputs such as Z/6Z and Z/9Z. The port retains integral JavaScript
numbers as integer inputs. Nonintegral floating-point/RDF promotion and roots
need their own profile; successful element conversion is never evidence of a
canonical parent map.

## Active goal: finish the elliptic-curve audit

Completion requires reviewing the implemented elliptic-curve API and reachable
branches against bundled Sage/PARI, fixing discrepancies with comparative
regressions, and verifying remaining items below. Passing one profile or closing
one batch does not complete this goal. Keep unrelated module audits separate.

- [ ] Close coordinate-parent, construction/model and coercion gaps (items 1–3).
- [ ] Repair and compare torsion callers and native dependency routing (item 4).
- [ ] Review other implemented elliptic modules and documented open discrepancies
      (isogenies, point/group algorithms, formal maps, rational/number-field paths).
- [ ] Resolve elliptic comparison timeouts and run the relevant complete suites.
- [ ] Refresh the elliptic source/API/branch coverage inventory and account for
      every remaining stub or unsupported branch before claiming completion.

## Next source and behavior work

1. **Remaining lift_x parents.** extend=True now adjoins y over QQ, prime fields
   and explicit finite extensions, using number fields, finite fields and quotient
   rings respectively. The returned points' equations, doubling and negation match
   original execution. Starting from general number-field/quotient parents, further
   lifting on the newly returned algebraic parents, cross-call field/curve identity
   and extension-construction RNG/backend parity remain unverified. The default
   Givaro/NTL root ordering is matched without claiming full backend equivalence.
2. **Other coordinate parents.** Audit real/complex inputs, other field/ring
   families and global coercion/factory metadata. Explicit finite extensions do
   not automatically coerce by degree divisibility: finite_field_base.pyx:1243
   requires compatible `_prefix` metadata for non-prime canonical embeddings.
   Add those dependencies before claiming general parent promotion. The current
   base_extend/change_ring methods preserve whole coefficients and unchanged-field
   curve identity. General ring morphisms and global curve caching remain open.
3. **Separate optimized finite-curve class.** The public default EllipticCurve
   export comes from ell_finite_field.ts, not constructor.ts/EllipticCurveGeneric.
   Its prime-field coordinate predicates/lifting now share the audited generic
   caller, with twelve regressions and fresh dependency-trace tests. Its broader
   types/backends remain prime-field oriented. Starting extension-field inputs
   bypass declared types and can fail in ellinit_Fp for two-coefficient calls.
   Five-coefficient construction now preserves the original equation through the
   generic implementation, including QQ/finite extensions and characteristic 2/3
   (25.0.0); coefficient coercion, tuple ownership, validation and tested point
   arithmetic match. This changes that runtime path to the generic point API.
   General scalar multiplication now delegates over prime fields p > 3 (25.1.0)
   and supported binary finite fields (26.0.0), plus characteristic-three and
   explicit odd extensions (27.0.0).
   Remaining finite-field backend integration is still required, along with two-term
   extension dispatch and native factory identity/caching. Promoted/extended lift
   results also use generic points.
4. **Torsion callers.** Generic torsion_points currently enumerates p constants
   (not all p^n elements), with only y=0/1 in characteristic two. Sage has no
   finite-field method by that name: the native rational implementation is in
   ell_rational_field.py:4126. Existing finite property checks compare cardinality,
   so document the alias and repair its domain/algorithm. It is not the same as
   _p_primary_torsion_basis: its relation dependency, negative bounds, selected
   basis ordering and division-point reduction/cache behavior are now repaired
   and compared (24.62.0). General point-order PARI routing is repaired for prime fields p > 3 (25.2.0),
   including cardinality/point caches and general-model backend conversion. Default general prime-curve cardinality/order and hybrid scheduling are repaired
   (25.3.0). Additional cardinality algorithms/options, group APIs, scalar/order routing for
   small characteristics/extensions, and finite torsion subgroup basis correction
   still need work. Optimized prime points now expose the same order algorithms
   and additive alias, including shared curve cardinality caching (25.4.0). Their
   prime scalar/coercion/cache path is repaired and compared in 25.5.0; native
   general characteristic-three/odd-extension scalar backends are integrated in
   27.0.0; legacy optimized characteristic-three scalars are integrated in 28.0.0. The optimized prime curve's
   points() now uses native group-basis enumeration, sorting and immutable caching;
   abelian_group caches its result and updates gens (24.63.0). This does not repair
   the separate generic torsion_points path. Group-record/tuple container semantics
   and exact native generator RNG behavior remain to be audited.
5. **Remaining hyperelliptic parents.** Root callers now preserve distinct/default
   root order, use PARI element comparison, and call scalar sqrt directly. Seven
   regressions and fresh comparisons check scalar state, binary lift root options,
   odd-degree models and Cantor reduction. Audit lift_x common-parent promotion
   (currently K(x)), other parent families and custom-field sqrt fallback; the
   same-parent root profile does not establish their correctness.
6. **Extension norm/charpoly — completed for the tested contracts (24.61.0).**
   Public norm calls charpoly('x'), which delegates to FF_charpoly. Native word
   interpolation/subresultants and arbitrary-prime interpolation are ported;
   nonmonic resultant moduli work. Fresh comparisons cover values, parent/variable,
   dependency calls, random state and PARI's degree-one binary charpoly error.
   The separate FpX_composedsum Newton/Laplace route remains open; its callers use
   the repaired resultant but still differ algorithmically from native composed sum.
7. **Finite polynomial irreducibility.** Word FLINT Shoup/DDF versus Rabin routing,
   extension NTL Iter/Det/Prob algorithms and constructor options remain open.
   See preserved audit-finite-irreducibility notes.
8. **Other coverage gaps.** Continue the source/API inventory in SCOPE.md and
   AUDIT-COVERAGE.md, including field embeddings, global coercion/factory metadata,
   unimplemented backend paths and remaining stubs. Eight lines in the shared
   half-GCD helper were unexecuted at 24.54.0: deferred constant-divisor resultant
   transition and asymmetric large-GCD remainder. Use real native cases, not
   artificial thresholds merely to increase line coverage.

## Rules when resuming

Read AGENTS.md, DESIGN.md, DEVIATIONS.md and SCOPE.md. Read bundled Sage and
native dependency source before changing algorithms. For every bug, retain a
comparative regression using the original implementation, including errors and
state where relevant. Preserve full integer precision when parsing native args.
Do not weaken assertions to accept another mathematically equivalent ordering.
Run appropriate callers and shared comparisons; update scope, deviations/API docs,
changelog and synchronized versions; commit completed batches. Do not claim full
coverage from dispatcher presence or raw line execution alone.

Small research handoffs formerly inside ignored transcripts are preserved in
`tests/audit/notes/`. Completed generated logs and duplicate transcripts are
storage artifacts, not the authoritative regression suite.


## Storage state

All 29 area suites use compact generator definitions and focused regressions.
97 bulk expected-output files have been replaced with live native calls.
Successful runs save no transcript caches. Failed runs save replay inputs in the
OS temporary directory; promote minimized counterexamples when fixing bugs.
Colocated native tests require SageMath and the bundled native toolchains.

Scratch output, unused audit JSON observations/baselines and duplicate transcripts
were removed. Historical audit notes may name those retired data files; their
behavioral conclusions remain in the notes. The 568-diagnostic baseline
remains at `tests/audit/baseline-typecheck.log.gz`. September audit commits were
squashed above `72b1e25`; historical commit IDs in audit notes are labels only.

## Local native build note (September 18)

An incomplete temporary PARI source cache initially lacked Configure. Rebuilding
from intact reference/pari succeeded after using `LD=cc`,
`LDFLAGS="-Wl,-search_paths_first -L/opt/homebrew/opt/readline/lib -L/opt/homebrew/lib"`
and `LIBS=-lreadline` for make. The repaired build remains outside the repository.
A future cold build may need the same local toolchain adjustment; the source and
comparative assertions were not changed to accommodate build failures.
