# Audit handoff — paused 2026-09-14

The user paused the behavioral audit to reduce repository size. Resume only when
requested. The September 9–14 commits were squashed to reduce Git storage. Their messages
and former IDs are preserved in `tests/audit/squashed-commits-2026-09.txt`; old
IDs in audit notes are historical labels and may no longer resolve in Git.
Current implementation code and the pending patch are preserved. The user authorized
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
bulk inputs or expected outputs are saved. The mathematical audit remains paused.
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

## What is verified

- Generic curve coordinate roots, Montgomery selection and coefficient formatting:
  5,822 new shared native comparisons; all 10,202 advanced-area rows match.
  Broader callers/docs/case-format: 21,571 pass, 13 existing skips.
- Extension square roots: 4,324 comparisons, including exact PARI root choice/state.
- Extension square predicates and native resultant/norm dependencies: 4,832
  comparisons; all 37,802 field-area rows matched at that checkpoint.
- Earlier completed polynomial, modular-ring, isomorphism, PARI and NTL batches
  are recorded in AUDIT-2026-09.md and SCOPE.md.
- Eight package builds pass. There are 584 pre-existing TypeScript diagnostics;
  compare full messages/continuations, not just the count.
- Recorded source-line execution: 111541/133728 (83.41%). This is NOT full
  behavioral or input-domain coverage. AUDIT-COVERAGE.md inventories the gaps.

## Unfinished work preserved before cleanup

`tests/audit/pending/curve-coordinate-coercion/` contains a patch and compressed
shared/native/old/partial-result records. The active source tree was restored to
the verified implementation before storage work so it does not contain a partial
coercion repair. Apply the saved patch only when resuming the audit.

- 5,250 coordinate-parent comparisons cover QQ, prime fields, explicit extensions,
  wrapped integers, bigint, strings, null, plain objects and booleans.
- Old code differs on 4,261 rows: 1,375 is_x_coord, 1,441 all-point lifts and
  1,445 single-point lifts.
- Partial patch makes is_x_coord always call K(x), preserving Sage's TypeError
  wrapper. Eight is_x_coord mismatches remain: the specialized GF2 constructor
  rejects integer strings and null that native GF(2) accepts. Fix at that field
  conversion boundary, with direct comparative regressions.
- lift_x parent promotion is not yet fixed. After the partial predicate change,
  2,894 total rows still differ. No full tests/builds were run for this patch.
- Nonintegral floating-point/RDF coordinates were deliberately not included in
  this initial matrix; their promotion and scalar hooks need a separate profile.

## Next source and behavior work

1. **Coordinate coercion and promotion.** Read ell_generic.py:710 and 913.
   is_x_coord always converts via K(x), replacing only TypeError with
   `x must be coercible into the base ring of the curve`. lift_x instead calls
   py_scalar_to_element (from sage.structure.coerce), then K.coerce_map_from(L),
   or promotes E to L if L canonically receives K. Constructor conversion is not
   canonical coercion: a QQ value can numerically convert to GF(p) yet lift_x
   rejects that parent. A constant extension element still promotes the curve.
   Preserve exact errors and returned parents. Audit base_extend/change_ring
   alongside this; it currently strips `.value` before conversion.
2. **Remaining lift_x options and inputs.** Port extend=True's fraction-field
   extension construction; support scalar wrappers and appropriate return types.
   Explicit finite extensions do not automatically coerce by degree divisibility:
   finite_field_base.pyx:1243 requires compatible `_prefix` metadata for the
   non-prime canonical embedding. Unknown fields/global factory metadata remain
   documented gaps. Update LLM.md for exported signature changes.
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
