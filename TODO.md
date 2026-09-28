# Audit handoff — updated 2026-09-28

The user resumed the behavioral audit on 2026-09-18 after repository cleanup.
Completed: is_x_coord input conversion, GF(2) constructor parity, and lift_x
canonical coercion/promotion over the tested scalar fields, and lift_x extend=True
over QQ, prime fields and explicit finite extensions.
Next: remaining coordinate parents and optimized finite-curve callers. The September 9–14 commits
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
   Its types/backends are prime-field oriented; its lift_x calls elllift_x while
   Sage inherits the generic scalar-root caller. Extension inputs bypass the
   declared types and currently fail in ellinit_Fp. Audit this API separately;
   the executed extension example correctly uses exported EllipticCurveGeneric.
4. **Torsion callers.** Generic torsion_points currently enumerates p constants
   (not all p^n elements), with only y=0/1 in characteristic two. Sage has no
   finite-field method by that name: the native rational implementation is in
   ell_rational_field.py:4126. Existing finite property checks compare cardinality,
   so document the alias and repair its domain/algorithm. It is not the same as
   _p_primary_torsion_basis, whose native division-point algorithm needs review.
5. **Hyperelliptic root callers.** field_ops.compare_elements uses integer
   representation instead of PARI field order. odd_degree_model selects the
   first DISTINCT root; cantor_reduction selects the first DEFAULT root. Both
   currently apply a wrong extra sort. sqrt_all_unsorted still adds predicate/
   binary-power work beyond the scalar method. Compare order, state and errors.
6. **Extension norm and trace.** Sage norm is charpoly-based
   (element_base.pyx:632, element_pari_ffelt.pyx:982 -> FF_charpoly); do not merely
   substitute FF_norm. Port/check native FpXQ/Flxq charpoly/resultant dependencies.
   Trace delegates to FF_trace. Current public methods use Frobenius sums/products.
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
behavioral conclusions remain in the notes. The 584-diagnostic baseline
remains at `tests/audit/baseline-typecheck.log.gz`. September audit commits were
squashed above `72b1e25`; historical commit IDs in audit notes are labels only.

## Local native build note (September 18)

An incomplete temporary PARI source cache initially lacked Configure. Rebuilding
from intact reference/pari succeeded after using `LD=cc`,
`LDFLAGS="-Wl,-search_paths_first -L/opt/homebrew/opt/readline/lib -L/opt/homebrew/lib"`
and `LIBS=-lreadline` for make. The repaired build remains outside the repository.
A future cold build may need the same local toolchain adjustment; the source and
comparative assertions were not changed to accommodate build failures.
