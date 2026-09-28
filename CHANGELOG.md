# Changelog

All notable changes to this project will be documented in this file.

## 24.59.1 - 2026-09-28

- Match hyperelliptic scalar square-root dispatch and PARI extension ordering;
  preserve first-root selection in odd-degree models and Cantor reduction, and
  request distinct roots for characteristic-two lifting.
- Add fresh live original/port comparisons, scalar state and dependency traces,
  and seven compact regressions. Six retained inputs fail against the prior
  implementation; the additional binary-root regression detects the final fix.
- Validate 207 live comparisons, all 121 existing hyperelliptic tests, 41 storage/
  case checks and eight builds; preserve all 570 baseline type diagnostics. Keep
  coordinate-parent promotion and custom-field fallback gaps explicit.

## 24.59.0 - 2026-09-28

- Route optimized prime-curve coordinate predicates/lifting through Sage's generic
  caller logic; repair conversion, errors, wrapped integers, promotion and ignored
  extension flags. Preserve optimized points over the original curve.
- Add fresh public-API and dependency-trace comparisons with twelve compact
  regressions, plus executable examples for both coordinate APIs.
- Validate 1,854 coordinate and 103 polynomial comparisons, 95 existing optimized
  curve tests and 375 docs/storage checks. Remaining constructor/model and broader
  audit gaps stay explicit; no bulk inputs or output snapshots added.

## 24.58.0 - 2026-09-28

- Implement generic curve lift_x's third positional extend flag over QQ and
  finite fields, adjoining y through number-field, finite-field or quotient-ring
  constructors while preserving native conjugate order and coordinate parents.
- Preserve rational-polynomial variable names through arithmetic and field
  definitions; bridge NumberFieldElement.isZero to its existing is_zero method.
- Add fresh original/port profiles and sixteen small regressions for extension
  creation, existing roots, disabled extension, promotion and dependency behavior.
  Compare returned point equations, parents, doubling and negation. Update the
  public API examples and remaining audit scope without adding bulk test data.
- Validate 1,242 live coordinate and 103 dependency comparisons, 245 existing
  caller/dependency tests, 373 docs/storage checks and eight package builds.
  Preserve all 570 existing type diagnostics; broader audit gaps remain open.

## 24.57.2 - 2026-09-18

- Match generic lift_x canonical scalar maps and curve promotion, including exact
  missing-parent/errors, extension coordinates and wrapped/modular integers.
- Preserve whole coefficients in base_extend/change_ring and original curve
  identity for an unchanged field; expose promotion-aware point overloads.
- Add fresh live comparisons and twenty small regressions, including constructor
  call traces. Remove the completed bulk coordinate research records and patch.
- Fix fourteen existing type diagnostics; the remaining 570 messages match the
  previous baseline. Update the baseline, executed API docs and audit handoff.
- Validate 3,500 historical lift comparisons, 929 fresh focused comparisons,
  307 smaller advanced-area comparisons, 298 caller tests and 374 focused/docs/
  storage checks; eight builds pass. Record the larger area's native timeout
  and remaining extend/real-coordinate/global-parent gaps without weakening tests.

## 24.57.1 - 2026-09-18

- Resume the behavior audit: make generic curve `is_x_coord` convert every input
  through its base field and preserve Sage's conversion error handling.
- Fix specialized GF2 string/default conversions and nonintegral-float errors;
  match the integer constructor's bare-object error message.
- Add fresh live original/port comparisons with replay seeds and nine small
  regression inputs. Update executed API examples and the remaining audit work.
  No new stored transcripts or bulk input corpus.
- Validate 1,750 historical predicate comparisons, 309 focused fresh comparisons,
  221 advanced-area comparisons and 740 focused/caller/docs/storage tests. Eight
  builds pass and all 584 existing TypeScript diagnostic messages are unchanged.
  General lift_x promotion and the broader audit gaps remain open in TODO.md.

## 24.57.0 - 2026-09-14

- Add `bun run quiet -- COMMAND` to stream full command output into temporary
  logs and return bounded summaries with exit status, replay settings, reported
  totals and diagnostic excerpts. Clip giant single lines as well as total output.
- Add bounded streaming log search/paging, literal argument forwarding, exit-code
  preservation and optional process-group deadlines. Document the wrapper as the
  default for noisy agent commands; do not rerun tests just to inspect their output.
- Validate large-output retention, early-failure visibility, UTF-8/ANSI handling,
  argument boundaries, inspection limits, missing commands and child cancellation.
  Six CLI tests and strict TypeScript checks pass; the full typecheck log retains
  all 584 existing diagnostics while its command summary is capped at 2,200 chars.

## 24.56.0 - 2026-09-14

- Retire 2,162,729 legacy explicit area inputs with the user's authorization.
  Generate fresh scalar, polynomial, byte, parser and constrained matrix inputs;
  retain compact constructor domains and named handwritten regression checks.
- Generate inputs for 88 native suites and execute the originals live. Replace
  ordinal regression lookups with stable IDs, print native replay settings, and
  retain failure-only artifacts for area mismatches and process failures.
- Remove the remaining alternate-named fixtures and unused audit baselines;
  reconstruct the 50,100-bit LLL regression from shifts and offsets. Fix native
  modular-root test dispatch to inspect the modulus instead of its random seed.
- Bound comparative subprocess lifetime, validate generator domains, and enforce
  a storage budget that rejects bulk input corpora and new output snapshots.
- Retain replay seeds for newly exposed lattice, rational, quaternion and finite
  polynomial mismatches; keep assertions strict and document the open findings.
- Validate 59 generator/replay/storage checks, two independent 340-comparison
  arithmetic runs, 164 final targeted native tests and the unchanged 584-diagnostic
  type baseline;
  all eight package builds pass. Full area comparisons still expose the open
  behavior differences and function-field timeout recorded in TODO.md.
- Keep the behavioral audit paused and document reduced historical input coverage.
  No production mathematical implementation or public API changes.

## 24.55.0 - 2026-09-14

- Generate property inputs once from a fresh replayable seed and run both original
  and port live. Retire cached-output CLI modes and passing transcript writes.
  Save only failing inputs with seed, generator/source version and runtime metadata.
- Replace 86 bulk native output snapshots (225,669,133 bytes) with live reference
  calls in 87 test modules. Preserve all 249,220 original inputs and their order,
  verified against independent fingerprints; keep positional regression checks.
- Reject missing, extra, ambiguously keyed and empty comparison result streams.
  Preserve structured fast-check counterexamples with seed/path metadata, use fresh
  default seeds, and prevent custom failure reporters from swallowing failures.
- Validate 55,890 live results against 16 former snapshots, 5,729 focused native
  tests, 3,261 final live-loader tests, 41 randomized arithmetic properties and two
  seeded live arithmetic runs (217 comparisons each). All 54 framework, storage
  and input-fingerprint checks pass; failure artifacts fail again on exact replay. All 584
  baseline type diagnostics and their continuation text are unchanged; eight builds pass.
- Retain legacy explicit input sweeps pending the user's choice of regression
  classification policy; this release does not claim their conversion to generators.

## 24.54.3 - 2026-09-14

- Squash the 265 September 9–14 audit and storage commits above `72b1e25`
  and reclaim superseded Git objects. Preserve current source and regression
  files, the paused audit handoff, and a compact inventory of former commits.
- Verify that changes from the previous tip are limited to history documentation
  and synchronized version metadata; production code and tests are unchanged.

## 24.54.2 - 2026-09-14

- Pause the audit at user request; write TODO.md and preserve the unfinished
  coercion patch, comparison records and research handoffs outside scratch output.
- Losslessly compress large shared case files and native fixtures. Support gzip
  in property discovery, normalization, comparison and coverage inventory, and
  write retained transcripts as atomic compressed streams by default.
- Verify all 52 compressed files against their original hashes, retaining
  2,304,091 case/fixture records while reducing 1,058,853,220 bytes to 214,538,563.
  Replay checks pass 15,853 tests; selected cached property checks pass. All
  canonical/storage checks, eight builds and unchanged type-diagnostic checks
  pass (584 existing errors).
- Remove completed scratch output and regenerable transcript caches; repack Git
  without rewriting history. The audit remains paused; TODO.md explains cache
  regeneration and how to resume the saved work.

## 24.54.1 - 2026-09-14

- Repair curve scalar-root delegation over QQ and extension fields, native
  y-coordinate order, distinct polynomial-root modes and Montgomery selection.
- Fix parentheses and spacing for extension coefficients in curve equations
  and their error messages. Add 5,822 shared native comparisons, including
  scalar random state, dependency traces and large-field subprocess watches.
- All 10,202 advanced-area comparisons match. Final broader callers/docs/case
  format pass 21,571 tests (13 existing skips); eight builds pass and all 584
  baseline type diagnostics are unchanged. Raw execution is 111541/133728
  (83.41%); generic curve execution records 922/960 lines.
- Correct the Curve25519 unit to native -486662 selection and use the exported
  generic curve class in the executed extension-field documentation example.

## 24.54.0 - 2026-09-14

- Add the native extension square predicate and its PARI resultant/norm
  dependencies. Preserve the existing GCD basecase dispatch while sharing
  native half-resultant accumulation and threshold selection.
- Add 4,832 distinct native comparisons: 4,414 dependency calls and 418 scalar
  predicates, including retained random state, nonmonic norms, common factors,
  recursive resultant thresholds and 276 scalar-record/default-modulus cases.
  Recognize the new scalar adapters in the audit inventory.
- All 37,802 field-area comparisons match. Broader callers/docs/case format
  pass 26,555 tests (13 existing skips); eight builds pass. The 584 existing
  type errors remain, with member details updated for is_square. Raw execution
  is 111577/133764 (83.41%); full-port behavioral coverage remains incomplete.

## 24.53.0 - 2026-09-14

- Add extension-element scalar sqrt and delegate to the bundled PARI
  FF_issquareall n=2 path, preserving exact root choice and random state.
- Add 4,324 distinct native comparisons (1,398 dependency and 2,926 API),
  including word/large/binary branches, option errors and boundary fields.
  Permanent watches cover large primes and degree-130 binary extensions.
- All 32,970 field-area comparisons match; final callers/docs/case format
  pass 15,341 tests (13 existing skips). Eight builds pass; 584 existing type
  errors remain, with member-detail changes from adding sqrt. Raw execution
  is 111392/133602 (83.38%); the new dependency records 117/117 lines.
- Correct the executed API example to import GFpn from the finite-rings subpath;
  document the output-pointer adapter and remaining predicate/caller gaps.

## 24.52.2 - 2026-09-14

- Share the Weierstrass generator with generic-curve isomorphism APIs; preserve
  first-result selection, native sorted lists and early exit for boolean queries.
- Implement all six morphism comparisons, native parameter ordering and the
  NotImplemented/null adapter. Check parent definitions and coerce equivalent
  prime-field implementations before curve arithmetic and comparison.
- Add 2,371 shared native comparisons, including full comparison matrices,
  invalid arguments, custom prime generators, parent definitions and root-call
  traces. All 4,380 advanced-area comparisons match native records.
- Broader callers/docs/case format pass 5,428 tests (13 existing skips). Eight
  builds pass; all 584 baseline type diagnostics are unchanged. Raw execution
  is 111330/133673 (83.29%). Update API examples, deviations and audit tracking.

## 24.52.1 - 2026-09-14

- Delegate all Weierstrass isomorphism root equations to distinct polynomial
  roots. Repair the characteristic-two quartic, extension roots/order, QQ
  failures and invalid-argument errors; remove the obsolete local root helpers.
- Add 1,759 shared native comparisons across finite fields and QQ, including
  large-prime/binary-extension boundaries, tuple validity and first-isomorphism
  selection. Fix 577 discrepancies in the initial/auxiliary comparison sets.
- Complete advanced area: 2,009 comparisons match. Broader elliptic callers,
  documentation and case-format checks pass 3,055 tests (13 existing skips).
  Eight builds pass; all 584 baseline type diagnostics are unchanged. Current
  raw execution is 111328/134058 (83.04%).

## 24.52.0 - 2026-09-14

- Match modular polynomial root parents, zero/constant errors and the bundled
  CRT/Hensel algorithms. Add cached field/factored-order methods and the residue
  root-lifting helper; preserve the bundled nonunit-linear bug explicitly.
- Add 5,126 shared native comparisons, including direct hook/caching/lifting
  cases and watched roots over composite moduli exceeding 10^12. Fix 3,132
  differences in the initial root comparison set; all new comparisons match.
- Preserve separate integer-lift and residue-root-lift test dispatchers, with a
  native-record runner regression. All 9,199 modular area comparisons match.
- Fast 239,012 tests pass (32 existing skips), affected slow callers 434 pass,
  and final focused checks 5,496 pass. Eight builds pass; 584 existing type
  diagnostics remain with display-order differences only. Raw execution is
  111416/134573 (82.79%). Update API examples, deviations and audit tracking.

## 24.51.0 - 2026-09-14

- Extract finite-field roots through factorization, preserving zero errors,
  factor order and binary factorization's PARI random-state effects.
- Add PARI's finite-field universal comparator and use it for extension factors,
  including low-first coefficients, packed binary words and signed word slots.
- Support distinct-root lists and use Sage's finite-field GCD route before
  factoring. Update elliptic division_points to request distinct roots; retain
  a watched degree-5,100 caller regression.
- Add 6,007 shared comparisons, including 236 PARI-state traces and 211 elliptic
  point lists. All new live comparisons match their native records.
- All 10,211 focused/caller tests and eight builds pass. No new type errors;
  584 existing diagnostics remain, with 39 private-member-count details changed.
  Updated raw execution coverage is 111281/134035 (83.02%).

## 24.50.2 - 2026-09-14

- Match Sage's ZZ/QQ factor and irreducibility backend dispatch, integer
  content ordering, zero/constant results and cached boolean predicates.
- Add 1,118 shared comparative regressions, including 304 PARI-state traces
  and watched degree-boundary cases. All comparisons and 4,917 polynomial,
  number-field caller, documentation and case-format tests pass.
- Eight builds pass; all 584 baseline type diagnostics are unchanged.
  Raw execution is 111247/135036 (82.38%); finite-field roots and
  irreducibility remain under audit.

## 24.50.1 - 2026-09-14

- Replace integer/rational root divisor enumeration with Sage's sparse-gap,
  FLINT GCD and NTL/PARI factorization routes. Preserve sparse root order and
  the bundled source's distinct ZZ/QQ zero errors.
- Add 699 shared comparisons, including 18 watched large-coefficient cases.
  Live 699/699 and all 1,625 broader polynomial/doc/case-format tests pass.
- Full fast: 226,757 pass, 32 skips, one existing NTL timeout. The exact case
  passes alone in 1.85 seconds at its unchanged 10-second deadline. Eight
  builds pass; all 584 baseline type diagnostics are unchanged.
- Current root helpers execute 87/87 recorded lines; polynomial_element
  5075/5306, aggregate recorded execution 111250/134402 (82.77%). Public
  factor/irreducibility and finite-field roots remain under audit.

## 24.50.0 - 2026-09-13

- Add PARI conjugate-count bounds, word derivatives and word squarefree predicates;
  delegate the legacy bound to native distinct-degree counting and prime iteration.
- Fix 40 negative-start bounds and a signed -2^63 scan hang. Add 1,500 native
  comparisons, watched regressions and five executed API expressions.
- Live 1,500/1,500, fast 226,059 tests (32 skips), 62 separate Galois tests and
  eight builds pass. All 584 baseline type diagnostics remain unchanged.
  Flx executes 193/193 recorded lines, galconj 2363/2395; aggregate recorded
  execution is 111293/135116 (82.37%). The whole-port audit continues.

## 24.49.2 - 2026-09-13

- Match Sage's unembedded number-field automorphism order, including signed
  integral-model scaling and the natural embedding's swap into first place.
- Add 190 comparisons using the bundled Sage scaling method and complete native
  PARI, fixing 94 ordered-result discrepancies. Exclude installed Sage 10.3's
  stale polredbest conversion behavior from the oracle.
- Live 190/190, broader 11,372 tests and eight builds pass; all 584 baseline
  type diagnostics are unchanged. Current number_field execution is 1544/1661
  lines; aggregate recorded execution is 111298/135150 (82.35%).

## 24.49.1 - 2026-09-13

- Delegate legacy integer polynomial factor helpers to PARI, preserving native
  factor order and random-state advancement and removing the local 200-bit cap.
- Add 224 native comparative regressions through degree 512 and 512-bit
  coefficients, including repeated inputs for the irreducibility predicate.
- Live 224/224 and focused 21,665 tests pass. Eight builds pass and the 584
  baseline type diagnostics are unchanged. Current pari_nf execution is
  1402/1465 lines; aggregate recorded execution is 111282/135135 (82.35%).

## 24.49.0 - 2026-09-13

- Complete NTL integer factorization with native deflation, adaptive van Hoeij
  recombination, additional lifting, exact lattice reduction and shared state.
- Route number-field irreducibility checks through NTL outside Sage's PARI degree
  window, fixing 16 comparative constructor failures above degree 300.
- Add 274 full native factorization traces, 100 Sage constructor comparisons and
  three full-word-prime precision comparisons. Correct degree-two overlap and
  bigint prime handling in comparative adapters; document the exported driver APIs.
- All 377 live comparisons and eight builds pass; the 584 baseline type diagnostics
  are unchanged. Full fast tier: 224,144 pass / 32 skips / zero failures, with
  2,086,767 assertions. The native factor driver executes 1418/1418 measured lines;
  aggregate recorded execution is 111396/135262 (82.36%).

## 24.48.0 - 2026-09-13

- Add explicit FFT-prime word coefficient contexts with native transform profiles,
  direct residue reconstruction and algorithm crossovers.
- Preserve FFT initialization and shared state through integer GCD, both exact
  certifications and squarefree decomposition. Add 519 native comparisons, fixing
  491 old trace differences; native diagnostics cover reconstruction retry paths.
- Live comparisons, 156,758 NTL tests, 143 Sage caller tests, five fallback
  comparisons and eight builds pass. All 584 baseline type diagnostics remain
  unchanged; two executed API examples cover the new context/state arguments.
  Aggregate raw execution is 111137/135197 (82.20%); the whole-port audit stays open.

## 24.47.0 - 2026-09-13

- Preserve Hensel factor-tree delegation, quotient/multiplier construction order,
  shared integer/FFT state and degree-dependent additional-lifting contexts.
- Fix missing FFT-limit and constant-quotient errors, incorrect word-context
  errors and nonunit-inverse messages. Add 362 shared native traces; the old
  implementation differs on 267 traces / 882 operations.
- Shared comparisons and eight builds pass; all 584 type diagnostics remain
  unchanged. Full fast tier: 223,231 pass / 32 skips / two timeouts, both passing
  unchanged in isolation. Add two executed API examples and an ownership check.
  Recorded execution is 111101/135118 (82.23%).

## 24.46.0 - 2026-09-13

- Add native word polynomial ordinary products and extended GCD, preserving word
  crossovers, recursive matrix transforms and pre-product FFT-size checks.
- Route polynomial inverses through word extended GCD and retain caller context
  in negative powers. Restoring the old call fails 23 traces / 85 operations.
- Add 780 shared native state comparisons, an ownership check and two executed
  API examples. All 72,026 focused tests and eight builds pass; the 584 baseline
  type diagnostics remain unchanged. Recorded execution is 111086/135104
  (82.22%); the whole-port audit remains open.

## 24.45.1 - 2026-09-13

- Preserve native word-multiplier errors when a rebuild changes the allocated FFT
  prime count, including retained buffers after plain builds and failures.
- Preserve original cached FFT residues during partial inversion and normalize
  complete cached products before packed arithmetic in word and arbitrary-modulus
  multiplication. Add 740 word and 88 arbitrary-modulus native traces, with
  regressions proven to fail before each repair. All 828 live comparisons and
  37,754 focused tests pass; eight builds pass and all 584 type diagnostics are
  unchanged. Recorded execution is 111010/135045 (82.20%).

## 24.45.0 - 2026-09-13

- Add native arbitrary-modulus polynomial quotient/multiplier array APIs, retained
  failed-build state, Newton/block reduction and the 800-prime CRT algorithm switch.
- Fix cached multiplier rebuilds accepting incompatible FFT prime counts; all 28
  comparative regressions fail before the guard and pass afterward.
- Add 1,499 native fixtures (1,458 distinct shared rows), including 360 watched
  constant-modulus boundary traces and two executed API examples. Focused 2,134
  tests pass; broad 221,256 pass/32 skips with one timeout passing unchanged in
  isolation. Eight builds pass and all 584 baseline type diagnostics are unchanged.
  Recorded current-source execution is 111008/135043 (82.20%).

## 24.44.2 - 2026-09-13

- Preserve ordinary FFT-size validation in cached-word multiplication after a
  failed oversized multiplier build, before coefficient cancellation or reduction.
- Add 309 shared native comparisons; 36 fail before the fix and all pass afterward.
  Focused 18,357 tests and 19,547 assertions pass; eight builds pass and all 584
  baseline type diagnostics are unchanged. Raw execution is 110644/134682 (82.15%).

## 24.44.1 - 2026-09-13

- Fix an infinite remainder loop after a failed constant-modulus rebuild by
  preserving NTL's negative-buffer-length error and its validation order.
- Add 240 shared native comparisons and a subprocess regression with a fixed
  watchdog, proven to fail before the fix and pass afterward. Focused 18,048
  tests and 19,238 assertions pass; eight builds pass and all 584
  baseline type diagnostics are unchanged. Raw execution is 110638/134676 (82.15%).

## 24.44.0 - 2026-09-13

- Preserve retained NTL transforms and native errors in projection after modulus or
  multiplier rebuilds; fix the defined plain-multiplier right-shift boundary.
- Port both transposed FFT primitives and cache seeded root powers and retained
  forward representations. Add 1,603 native fixtures covering 1,579 distinct shared
  rows, an ownership check and two executed API examples.
- Existing/new projection tests pass 10,461; final focused tests pass 7,165; fast
  tests pass 219,802 (32 existing skips), with two timeouts passing unchanged
  isolated checks. Eight builds pass and all 584 baseline
  type diagnostics are unchanged. Raw execution is 110637/134683 (82.15%); the audit remains open.

## 24.43.0 - 2026-09-13

- Preserve independently updated NTL word-modulus and multiplier caches after
  failed rebuilds, including cyclic remainders and coefficient-cancellation errors.
- Port forward/inverse truncated transforms and prefix rounding; reconstruct stale
  cache interpolation with native roots, exact one-prime centering and fused CRT.
- Add 5,210 native comparison fixtures spanning 5,186 distinct shared rows, three
  unit checks and four executed examples. Focused 5,523 tests pass; the broad run
  passes 218,159 (32 skips), and three timeouts pass unchanged isolated checks.
  Eight builds pass and all 584 type diagnostics are unchanged.
  Raw execution is 110541/134581 (82.14%). The whole-port audit remains open.

## 24.42.0 - 2026-09-13

- Add a compact `rows` form to property-test case files (`[seed, arg1, ...]` per case) and
  `tests/property/normalize-cases.ts`, which converts all-`fixedValue` cases to rows, merges
  adjacent same-function row cases, drops exact duplicates and writes one row per line. Both
  runners accept rows alongside generator cases; `case-format.test.ts` rejects non-canonical files.
- Normalize all 29 case files: 1083 MB to 687 MB on disk, 2,170,093 to 2,122,406 executions
  (47,687 exact duplicates removed, no distinct case added or lost; transcripts verified
  identical for TypeScript and SageMath runners). Per-case `description`/`trials` metadata,
  which no runner read, is dropped from row cases.

## 24.41.0 - 2026-09-13

- Preserve NTL integer polynomial product/square FFT initialization through optional
  shared state. Export and delegate to the original ChooseSS dispatch predicate.
- Add 1,497 shared comparisons (1,490 defined native and seven arithmetic guards),
  two unit checks and three executed API examples. Fix 69 prior cache/random-byte
  discrepancies; native diagnostics confirm every multiplication/squaring branch
  and every SS selection criterion with unchanged arithmetic results.
- Focused 8,736 and fast 212,945 tests pass (32 existing skips); eight builds pass
  and all 584 baseline type diagnostics remain unchanged. Raw execution is 110250/134267 (82.11%).

## 24.40.0 - 2026-09-13

- Propagate shared NTL cache/stream state through local factor-information updates,
  selected/complementary products, exact division and cardinality search.
- Add 4,344 shared comparisons (4,339 native and five explicit indexing guards),
  one ownership check and three executed API examples. Fix 112 earlier trace
  discrepancies, including 80 with later output or retained-information differences.
- Focused 12,084 and fast 211,443 tests pass (32 existing skips); eight builds pass
  and all 584 baseline type diagnostics are unchanged. Combined execution is
  110217/134233 (82.11%); Hensel state and full van Hoeij iteration remain under audit.

## 24.39.0 - 2026-09-13

- Preserve cold NTL word-context initialization through optional shared FFT cache
  and random-stream state, including factor recovery order and failed context
  construction. SmallPrimeFactorization accepts the shared cache as an option.
- Add 3,331 native comparisons, two unit checks and two executed API examples.
  Fix 2,786 trace cache/stream discrepancies, including 864 traces with different
  factor outputs or retained information. Focused 12,318 and fast 207,095 tests
  pass (32 existing skips); eight builds pass and all 584 baseline type diagnostics
  remain unchanged. Combined raw execution is 110206/134223 (82.11%).

## 24.38.0 - 2026-09-13

- Port NTL GotThem factor acceptance with native row order, retained prefixes,
  constant-term and coefficient-bound filters, exact trial division, diagnostics
  and shared FFT cache/random state.
- Add 1,669 shared comparisons (1,667 defined native cases and two adapter guards),
  three unit checks and two executed API examples. All 8,771 focused and 203,760
  fast tests pass (32 existing skips); all eight package builds pass and all 584
  baseline type diagnostics remain unchanged.
- Combined raw execution coverage is 110188/134208 (82.10%). Cold word-context initialization
  and the complete van Hoeij iteration remain under audit.

## 24.37.0 - 2026-09-13

- Correct full-library oracle AES profile guards; unsupported-profile checks and
  fresh native replays of 3,547 division/product traces pass.

- Preserve NTL exact polynomial division's degree-first zero-numerator rejection,
  retained outputs, and distinct PlainDivide/HomDivide behavior. Add explicit shared
  prime/cache state for reconstruction and the original retry/certification order.
- Add 2,870 native comparisons, one ownership check and three executed API examples.
  Fix 181 old value/status discrepancies and 172 first-call state discrepancies.
  Native diagnostics confirm skipped primes and all targeted retry/rejection paths.
- Focused 10,793 and fast 202,086 tests pass (32 existing skips); eight builds pass
  and all 584 baseline type diagnostics remain unchanged. Combined raw execution
  coverage is 110123/134143 (82.09%).

## 24.36.0 - 2026-09-13

- Preserve NTL modular-product FFT cache initialization and random-byte consumption
  through explicit shared state, including native alias-square and SS dispatch.
- Add 677 shared comparisons (670 defined native cases, seven labeled guards),
  two unit checks and three executed API examples. Fix 99 old state mismatches.
  Build and validate the complete bundled NTL for future comparisons.
- Focused 4,223 tests and fast 199,212 tests pass (32 existing skips); all eight
  builds pass and all 584 baseline type diagnostics remain unchanged.
  Combined raw execution coverage is 110019/134067 (82.06%).

## 24.35.0 - 2026-09-13

- Port NTL's certified integer elimination helper with native modular-rank,
  inverse-zero and exact-form retries and full cache/stream fidelity.
- Add 3,250 native comparisons, one ownership check and two executed API examples.
  A diagnostic native replay confirms all three seeded retry paths without changing
  their outputs. Focused 7,138 and fast 198,530 tests pass (32 existing skips).
- Builds pass; all 584 baseline type diagnostics remain unchanged. All 46 new
  emitted lines execute; combined raw execution coverage is 109978/134028 (82.06%).

## 24.34.0 - 2026-09-13

- Port NTL integer determinant/inverse reconstruction with native FFT-cache,
  random-stream, CRT, probabilistic and deterministic stopping behavior.
- Add 3,595 native comparisons, two ownership checks and four executed API
  examples. Focused 12,790 and fast 195,277 tests pass (32 existing skips);
  builds pass and all 584 baseline type diagnostics are unchanged.
- New source execution is 109/114 lines; post-certificate rescaling remains
  unexercised. Combined raw execution coverage is 109932/133962 (82.06%).

## 24.33.0 - 2026-09-13

- Port NTL scalar/matrix CRT, signed reconstruction intervals, integer determinant
  bounds and delayed-reduction determinants over arbitrary-size moduli.
- Add 8,881 native comparisons, 22 labeled input guards, two ownership checks and
  four executed API examples. Focused 13,814 and fast 191,676 tests pass
  (32 existing skips); builds pass and all 584 baseline type diagnostics are unchanged.
- All 158 new emitted lines execute; combined raw execution coverage is 109823/133833 (82.06%).

## 24.32.0 - 2026-09-13

- Port NTL FFT-prime recognition, descending candidate search and lazy prime/root
  tables with explicit cache and random-stream state.
- Add 3,230 native comparisons, 29 labeled input guards, two ownership checks
  and four executed API examples, including a seeded rare special-factor rejection.
- Focused 8,166 and fast 182,767 tests pass (32 existing skips); builds pass and
  all 584 baseline type diagnostics are unchanged. All 147 new emitted lines
  execute; combined raw execution coverage is 109665/133675 (82.04%).

## 24.31.0 - 2026-09-13

- Port NTL word-matrix elimination, left-kernel bases, inverse/determinant and
  left/right solves, preserving strict/relaxed pivots and retained singular outputs.
- Add 5,086 native comparisons, eight labeled guards, two ownership checks and
  four executed API examples. Focused 6,620 and fast 179,502 tests pass
  (32 existing skips); builds pass and all 584 baseline type diagnostics are unchanged.
- All 235 new emitted lines execute; the matrix module is 346/346 and combined
  raw execution coverage is 109518/133528 (82.02%).

## 24.30.0 - 2026-09-13

- Port NTL probable-prime tests, probability bounds and sequential/seeded-block
  prime generators with explicit random streams and retained legacy overloads.
- Add 4,621 native comparisons, three labeled word guards, one compatibility
  check and three executed API examples. Focused 12,631 and fast 174,402 tests
  pass (32 existing skips); package builds pass and 584 type diagnostics are unchanged.
- All 206 newly emitted lines execute; combined raw coverage is
  109283/133298 (81.98%). Rare candidate-retry/exhaustion paths remain source-reviewed.

## 24.29.0 - 2026-09-13

- Port NTL cardinality-based integer-factor recombination, including unsigned
  ratio filters, iterative pruning tables, cached tests and local-information refresh.
- Add 2,879 matching native comparisons, five labeled vector guards and four
  ownership/API tests. Focused 13,844 and fast 169,774 tests pass (32 existing skips).
- Builds pass; all 584 baseline type diagnostics remain unchanged. All 426 new
  emitted lines execute; combined raw coverage is 109077/133107 (81.95%).

## 24.28.0 - 2026-09-13

- Port NTL retained factor-information updates, including partial pattern writes,
  saved contexts, suffix-degree caches and optional additional-prime selection.
- Add 1,216 matching native stateful comparisons and four ownership/API tests.
  Focused 10,956 and fast 166,886 tests pass (32 existing skips); builds pass.
- All 584 baseline type diagnostics remain unchanged. The edited module executes
  570/570 emitted lines; overall raw coverage is 108651/132651 (81.91%).

## 24.27.0 - 2026-09-13

- Port NTL's stateful small-prime selection for integer polynomials, preserving
  ordered modular factors, retained vector storage, contexts and random-stream bytes.
- Fix 102 prototype signed-bit discrepancies with permanent comparative regressions.
- Add 2,175 native stateful comparisons and four ownership/API tests. Focused 9,736
  tests pass; the broad suite records 165,664 passes, 32 skips and two existing
  test timeouts. Both tests pass unchanged in isolation; builds pass.
- All 584 baseline type diagnostics remain unchanged. The edited module executes
  494/494 emitted lines; overall raw coverage is 108575/132556 (81.91%).

## 24.26.0 - 2026-09-13

- Add explicit degree-dependent NTL word-context settings and propagate them through
  polynomial GCD, minimum polynomials, root products and modular factor recovery.
- Preserve native context thresholds, FFT-size errors and inverse-error ordering;
  fix 66 observed prototype discrepancies with permanent comparative regressions.
- Add 6,500 defined native comparisons, six labeled numeric guards and four ownership/API
  tests. Focused 6,772 and fast 163,487 tests pass (32 existing skips); builds pass.
- All 584 baseline type diagnostics remain unchanged. Updated modules execute
  1007/1007 emitted lines; overall raw coverage is 108479/132460 (81.90%).

## 24.25.0 - 2026-09-13

- Fix word-polynomial division errors leaking the capital-modulus error message;
  retain 70 observed failing cases and 20 passing controls as comparative regressions.
- Port NTL root splitting, ordered equal-degree recovery and squarefree factorization.
  Add 4,433 defined native comparisons, 120 labeled guard traces and four API/ownership tests.
- Focused 7,298 and fast 156,977 tests pass (32 existing skips); builds pass and
  584 baseline type diagnostics are unchanged. Updated word modules execute 543/543 lines.
- Raw execution 108380/132361 (81.88%). Integer-factorization dependencies remain open.

## 24.24.0 - 2026-09-13

- Port NTL cached bounded samplers, random word polynomials and monic root products,
  preserving sampler assignment, native iteration/tree cutovers and zero roots.
- Add 2,933 defined native comparisons and 406 explicitly labeled guard traces,
  plus five ownership/API tests. Compare outputs, cache state, errors and stream bytes.
- Focused 22,018 and fast 152,420 tests pass (32 existing skips); builds pass and all
  584 baseline type diagnostics are unchanged. All 82 new routine lines execute.
- Raw execution 108238/132219 (81.86%). Root splitting and modular factor recovery continue.

## 24.23.0 - 2026-09-13

- Port NTL quotient-element minimum polynomials, including native projection,
  probabilistic candidates, certification and exact random-stream retry behavior.
- Add 11,713 native comparisons covering polynomial results, exception types/messages
  and subsequent stream bytes, plus four ownership/API example checks.
- Focused 15,218 and fast 149,076 tests pass (32 existing skips); builds pass and
  all 584 baseline type diagnostics are unchanged. Updated module executes 334/334 lines.
- Raw execution 108155/132136 (81.85%). Root splitting and modular factor recovery continue.

## 24.22.0 - 2026-09-13

- Port NTL word and arbitrary-integer random samplers with explicit stream contexts,
  preserving native overloads, rejection order, early returns and error timing.
- Add 2,626 native comparisons and six explicit signed-word guard controls,
  checking sampled values and subsequent stream bytes. Add four ownership/API tests.
- Focused 7,988 and fast 137,359 tests pass (32 existing skips); builds pass and
  all 584 baseline type diagnostics are unchanged. New source executes 81/81 lines.
- Raw execution 108104/132085 (81.84%). Element minimum polynomials and factor recovery continue.

## 24.21.0 - 2026-09-13

- Port NTL's generic deterministic byte stream, SHA/HMAC helpers and key derivation,
  preserving copy/assignment state, nonce resets and counter carry/wrap.
- Preserve the native HMAC 32-bit input-length conversion found during source review;
  add 216 boundary comparisons without multi-gigabyte allocations.
- Add 4,038 native comparisons and 28 explicit short-buffer guard controls, plus
  five ownership/API example tests. Focused 5,352 and fast 134,723 tests pass
  (32 existing skips); builds pass and all 584 baseline type diagnostics are
  unchanged. New source executes 198/198 emitted lines.
- Raw execution 108023/132004 (81.83%). Integer sampling and factor recovery continue.

## 24.20.0 - 2026-09-13

- Port cached NTL polynomial multipliers, transposed multiplication and projected
  powers, including native correction caching and failed-build state transitions.
- Add 8,432 native comparisons plus 423 explicit guard controls for a recorded
  upstream prepared-zero projection crash. Recheck that native boundary whenever
  source/compiler/header profiles change; preserve earlier validation errors.
- Add 8,860 unit/example checks. Focused 36,282 and final fast 130,652 tests pass
  (32 existing skips); builds pass and all 584 baseline type diagnostics are
  unchanged. Both updated modules execute 467/467 emitted lines. An existing
  factorization timeout passes isolated/full reruns without code or deadline changes.
- Raw execution 107825/131806 (81.81%). Modular factor recovery and the audit continue.

## 24.19.0 - 2026-09-13

- Port NTL PrimeSeq with its native segmented sieve, independent iterator state,
  block resets and finite exhaustion/restart behavior.
- Add 1,031 shared operation traces and two executed examples, including 40,000
  consecutive primes. Focused 1,279 and fast 121,792 tests pass (32 existing
  skips); builds pass and all 584 baseline type diagnostics are unchanged.
  New code executes 80/80 emitted lines. No mismatch was observed in this batch.
- Raw execution 107713/131694 (81.79%). Modular factor recovery and the audit continue.

## 24.18.0 - 2026-09-13

- Port NTL sequence minimum-polynomial reconstruction with native Berlekamp–Massey
  and half-GCD cutovers, sharing the existing word-polynomial GCD transformation.
- Add 3,245 shared comparisons and 3,249 unit/example checks. Focused 16,538
  and fast 120,759 tests pass (32 existing skips); builds pass and all 584 baseline
  type diagnostics are unchanged. The updated module executes 212/212 lines.
  No mismatch was observed in this batch.
- Raw execution 107633/131614 (81.78%). Modular factor recovery and the audit continue.

## 24.17.0 - 2026-09-13

- Port NTL word-polynomial distinct-degree factorization and the first squarefree
  Cantor–Zassenhaus stage, retaining baby/giant steps, batched GCDs and factor order.
- Add 2,481 shared comparisons and 2,485 unit/example checks through degree 1024.
  Focused 13,289 and fast 117,510 tests pass (32 existing skips); builds pass
  and all 584 baseline type diagnostics are unchanged. The updated factorization
  module executes 211/211 emitted lines. No mismatch was observed in this batch.
- Raw execution 107590/131571 (81.77%). Equal-degree splitting and the audit continue.

## 24.16.0 - 2026-09-13

- Port NTL word polynomial block composition, cache reduction, shared compositions,
  trace maps, iterated composition and native half-GCD scheduling.
- Fix native allocation-error precedence in the initial PowerCompose adapter;
  preserve 120 failing prototype comparisons and add 45 boundary controls.
- Add 10,559 shared comparisons and 10,565 unit/example checks. Focused 10,804
  and fast 115,025 tests pass (32 existing skips); builds pass and all 584
  baseline type diagnostics are unchanged. New modules execute 217/217 lines.
- Raw execution 107427/131408 (81.75%). Full factorization and the audit continue.

## 24.15.0 - 2026-09-13

- Port NTL word matrix multiplication with native Strassen–Winograd scheduling,
  exact packed dot products, odd-dimension corrections and empty-shape adapters.
- Add 1,246 shared comparisons and 1,250 unit/example checks, including recursive
  products through 899 by 899. Focused 1,487 and final fast 104,460 tests pass
  (32 existing skips); builds pass and all 584 baseline type diagnostics remain
  unchanged. New code executes 111/111 emitted lines. Two existing constructor
  timeouts pass isolated/full reruns without changing their limits or code.
- Raw execution 107210/131191 (81.72%). Full factorization and the audit continue.

## 24.14.0 - 2026-09-13

- Port NTL word-polynomial quotient construction, cached reduction, products,
  inverses and binary powering. Preserve native reciprocal construction timing,
  fixing 135 observed mismatches in the initial adapter on composite moduli.
- Add 10,880 shared comparisons and 10,884 unit/example checks. Focused 11,119
  and final fast 103,210 tests pass (32 existing skips); builds pass and all 584
  baseline type diagnostics remain unchanged. The new module executes 143/143
  emitted lines. An existing automorphism timeout passes both isolated and full
  reruns without changing its timeout or implementation.
- Raw execution 107099/131080 (81.71%). Full factorization and the audit continue.

## 24.13.0 - 2026-09-12

- Port eleven NTL factor selection helpers, including degree bitsets, cached
  constant-term divisibility, balanced copying and minimum-degree products.
- Add 7,291 shared comparisons and 7,295 unit/example checks. Focused 7,528
  and fast 92,326 tests pass (32 existing skips); builds pass and all 584 baseline
  type diagnostics remain unchanged. ZZXFactoring executes 398/398 emitted lines.
- Raw execution 106956/130937 (81.69%). Full factorization and the audit continue.

## 24.12.0 - 2026-09-12

- Port NTL factor recovery evaluation/root bounds, lattice row removal and
  additional Hensel lifting, retaining native signed division and error behavior.
- Add 3,848 shared comparisons and 3,852 unit/example checks. Focused 4,083
  and fast 85,031 tests pass (32 existing skips); builds pass and all 584 baseline
  type diagnostics remain unchanged. ZZXFactoring executes 278/278 emitted lines.
- Raw execution 106836/130817 (81.67%). Full factorization and the audit continue.

## 24.11.0 - 2026-09-12

- Add native exact NTL LLL/LLL_plus, integer image and lattice solving, with
  transformations, dependent rows and ordered integer Gram determinants.
- Add the signed extended-GCD primitive and preserve NTL's zero convention,
  rounding ties and unsuccessful solver output behavior.
- Add 5,576 shared comparisons and 5,581 unit/example checks. Focused 5,810
  and fast 81,179 tests pass (32 existing skips); builds pass and all 584 baseline
  type diagnostics remain unchanged. New kernels execute 219/219 lines.
- Raw execution 106759/130742 (81.66%). Factor recovery and the full audit continue.

## 24.10.0 - 2026-09-12

- Add native NTL multifactor Hensel lifting with its minimum-degree factor tree,
  staged inverse corrections and restoration of the original factor order.
- Supply exact integer products through the existing portable NTL multiplication
  boundary. Preserve native overflow validation at 2^60 and error precedence.
- Add 4,058 shared comparisons and 4,062 unit/example checks. Focused 8,582,
  fast 75,598 and 729 existing GCD/squarefree callers pass (32 existing skips).
  Builds pass; all 584 baseline type diagnostics remain unchanged.
- Raw execution 106540/130523 (81.63%). NTL lattice factor recovery, constructor
  routing outside the PARI degree window and the broader audit remain open.

## 24.9.0 - 2026-09-12

- Add eight native NTL integer-factor trace dependencies: cached Newton sums,
  sparse/dense chopping, precision updates and factor-lattice construction.
- Preserve native validation order, balanced rounding, output dimensions and
  retained vector tails; return independent vectors and row matrices.
- Add 4,289 shared comparisons against bundled original bodies and 4,293
  unit/example checks. Focused 4,520 and fast 71,536 tests pass (32 existing skips).
  Builds pass; all 584 baseline type diagnostics are unchanged.
- Updated NTL module executes 96/96 lines; combined raw coverage is
  106409/130346 (81.64%). Full NTL factorization and the broader audit continue.

## 24.8.0 - 2026-09-12

- Complete PARI integer polynomial factorization stages: Newton traces, lattice
  recombination, prime selection, deflation, squarefree and full factor output.
- Route degree-30–300 number-field irreducibility checks through Sage's PARI
  factorization path, repairing reproduced timeouts and unsupported boundaries.
  Cache repeated checks with Sage's boolean/random-state behavior.
- Add 5,599 shared comparisons and 5,605 unit/example checks. Focused: 84,850 pass;
  fast: 67,219 pass, 32 existing skips; 24 final knapsack controls also pass.
  Builds/native suites pass; 584 baseline type diagnostics unchanged.
- Raw execution 106329/130122 (81.71%). Native NTL constructor routing outside
  the PARI degree window and the broader audit remain in progress.

## 24.7.0 - 2026-09-12

- Add native ZM_hnf, ZM_hnfall_i, ZM_hnfall and hnfall, retaining exact basis
  transformations, all zero-column removal modes and the column-count dispatch.
- Share the native Bézout/sign/remainder column operations with permuted HNF;
  keep caller inputs independent of returned matrices and transformations.
- Add 11,200 shared comparisons and 11,202 unit/example checks. Focused: 32,459
  pass; fast: 61,638 pass with 32 existing skips. Builds/native suites pass;
  584 baseline type diagnostics remain unchanged. Execution 106009/129864 (81.63%).
- Preserve separate private-kernel stress observations; public dispatch handles
  those same wide matrices. Constructor factorization and the audit continue.

## 24.6.0 - 2026-09-12

- Port bounded PARI modular-factor recombination with its native trace filters,
  candidate order, precision selection and unresolved-group stopping decisions.
- Add exact bounded quotients, centered scalar remainders, the native subset cap
  and balanced modular polynomial products. Preserve native leading-bound and
  singleton/empty-product behavior.
- Add 7,437 shared comparisons and 7,441 unit/example checks. Native and caller
  comparisons pass; two new documentation import paths were corrected and all
  220 examples pass. Fast: 50,436 pass, 32 existing skips. Builds/native suites
  pass; 584 baseline type diagnostics remain unchanged. Execution 105947/129805 (81.62%).

## 24.5.0 - 2026-09-12

- Add native permuted HNF, knapsack recognition, LLL factor progress and
  coefficient/root bounds, with their binomial, sparse evaluation, real ceiling,
  half-power and mixed-comparison dependencies.
- Match Qdivii zero-denominator errors, fixing five reproduced existing mismatches.
- Add 17,092 shared original comparisons and 17,099 unit/example checks. All
  46,543 focused and 42,995 fast tests pass (32 existing skips); builds and native
  suites pass, with 584 baseline type diagnostics unchanged.
- Execution 105787/129606 (81.62%). Constructor factorization and the full audit remain open.

## 24.4.0 - 2026-09-12

- Add native ZM_lll_norms output to the adaptive lattice wrapper, preserving
  unset norm pointers, precision escalation and caught upstream errors. Existing
  ZM_lll keeps its signature and result shape.
- Add 2,318 shared native probes and permanent regressions, including separate
  resource and bounded nontermination observations; keep 66 two-second timeout
  controls in the explicitly listed slow unit file.
- Validate 123,811 focused tests plus the single isolated deadline replay, and
  25,896 fast tests (32 existing skips). Builds/native suites pass and 584
  baseline type diagnostics remain unchanged. Execution 105574/129259 (81.68%); audit continues.

## 24.3.0 - 2026-09-12

- Honor supported trace/norm base-field arguments: QQ, the parent and degree-one
  fields. Preserve Sage's field-valued results and degree-one norm exception.
- Match quadratic trace arity, generic argument counts, primitive-input errors
  and impossible embedding degrees. Keep other-field/morphism dependencies open.
- Add 672 shared Sage comparisons and one executed API example, covering 276
  reproduced old mismatches. All 66,757 focused and 23,667 fast tests pass (32
  existing skips); builds pass and 584 baseline type errors remain unchanged.
  Execution 105557/129017 (81.82%); the full audit continues.

## 24.2.0 - 2026-09-12

- Replace the full characteristic-polynomial computation used by absolute
  number-field trace/norm with Sage's quadratic formulas and dedicated PARI
  rational quotient trace/norm routes, including scalar shortcuts.
- Add native QQ multiplication, monic/generic remainders and quotient arithmetic
  dependencies. Cover nonmonic equations, signed storage and trailing zero input.
- Add 3,561 shared original comparisons, 12 bounded operation regressions and an
  executed API example. All 95,430 focused and 22,994 fast tests pass
  (32 existing skips). Builds/native pol suite pass; 584 baseline type errors remain.
- Record the separate constructor irreducibility slowdown and remaining optional
  bases/characteristic-polynomial gaps. Execution 105490/128950 (81.81%); audit remains open.

## 24.1.0 - 2026-09-12

- Implement general integral/fractional ideal factorization and NumberField.factor
  inputs, preserving negative exponents, cached results and Sage factor order.
  Remove false-unit results for fractional ideals, including norm-one nonunits.
- Delegate full Kummer/Buchmann--Lenstra prime decomposition and ideal primality to
  PARI. Add finite-field image, supplementary basis, inverse and prime-HNF kernels.
- Replace unbounded ideal HNF construction with PARI principal-HNF and idealadd
  routines; fix BigInt exhaustion on degree-16 ideals with small generators.
- Add 7,137 permanent native/Sage comparisons plus an executed API example.
  All 122,642 combined and 19,420 fast tests pass (32 existing skips); native
  suites/builds pass and 584 baseline type errors remain unchanged.
  Execution 105315/128779 (81.78%); broader audit coverage remains open.

## 24.0.0 - 2026-09-12

- **Breaking:** correct ideal valuation to `I.valuation(P): bigint | 'Infinity'`;
  migrate Galois callers and document the signature with an executed API example.
- Delegate fractional HNF valuation and integral-vector valuation to native PARI
  algorithms; finalize prime data with the original uniformizer/kernel route.
  Remove the fractional zero loop, arbitrary exponent cap and prime_below trial search.
- Add 8,753 shared comparisons and 2,493 unit/example checks. Preserve raw display/
  string discrepancies and the native low-level zero-vector hang as separate controls.
- All 46,645 number-field comparisons, 49,961 focused tests and 12,282 fast tests
  pass (32 existing skips). Builds pass; 584 baseline type errors remain unchanged.
  Execution 104754/128227 (81.69%). Factorization and broader coverage remain in progress.

## 23.0.13 - 2026-09-12

- Match native class-specific ideal method errors and nonprime error messages;
  construct prime-decomposition results through the fractional-ideal factory.
- Add 1,948 shared comparative cases, 260 native unit fixtures and an executed
  API example. Repair nine legacy test setups that constructed the wrong class.
- Full number-field comparisons 37,892, focused tests 38,805 and fast tests 9,789
  pass (32 existing skips). Builds pass and 584 baseline type errors are unchanged.
  Combined execution 104604/128122 (81.64%); valuation and factorization remain under audit.

## 23.0.12 - 2026-09-12

- Delegate Gram LLL and algebraic relations through their distinct native PARI
  routes. Preserve HNF kernel columns and native Smith representatives, including
  half-tie rounding, diagonal residues and machine-word reduction thresholds.
- Match quaternion reduced-basis orientation and integer quadratic-form errors;
  pin the bundled source's exact conjugators and generator images.
- Add 2,938 shared cases and 2,298 unit/example checks. Native PARI 2.18.1 replaces
  installed-version relation and quaternion dependency oracles; 27 resource-limit
  failures remain separate controls rather than exact matches.
- Matrix comparisons 58,456, algebraic dependencies 10,386 and quaternions 855 pass;
  focused 51,833 and fast 9,528 tests pass (32 skips). Builds pass, 584 baseline
  type diagnostics are unchanged, and execution coverage is 104577/128070 (81.66%).

## 23.0.11 - 2026-09-12

- Match Sage's distinct membership, sum/product, quotient, divisibility and
  coprimality coercion rules; support exact-length field coefficient vectors.
- Restore PARI's centered maximal-order basis normalization and exact half-tie
  integer rounding; delegate rational ideal products with primitive-content and
  dimension handling. Keep existing class-number certificates basis-independent.
- Add 9,176 shared comparisons, 128 unit/example checks and native kernel fixtures.
  The first input/operator/basis probes contain 5,049 pre-fix discrepancy records.
- All 35,944 live comparisons, 36,720 focused and 7,230 fast tests pass (32 skips);
  builds pass and 584 type diagnostics remain unchanged. Coverage 104657/128068 (81.72%).
  Gram/HNF/SNF dependencies and the broader audit continue.

## 23.0.10 - 2026-09-12

- Route ideal intersection through the PARI LLL-kernel implementation, preserving
  fractional HNF normalization and supporting Sage's scalar/list input coercion.
- Fix ideal factory and direct constructor class, identity, generator and zero-error
  semantics, plus rational constants from foreign abstract number fields.
- Add 5,126 shared Sage/native comparisons covering 3,478 observed old mismatches,
  native denominator normalization, input preservation and random-state controls.
- All 26,768 number-field comparisons, 27,416 focused tests and 7,102 fast
  tests pass (32 existing skips); builds and API examples pass; 584 pre-existing type
  diagnostics are unchanged. Coverage 104535/128028 (81.65%); the broader audit continues.

## 23.0.9 - 2026-09-12

- Port native adaptive integer LLL, typed lllfp and recursive FLATTER, preserving
  tuning tables, rank augmentation, stage precision schedules and transformations.
- Fix all four LLL stages to multiply by the word/truncated mantissa before shifting.
  A 50,100-bit regression drops from 15.2s to 2.1s while matching native results.
- Add native-compatible definite-form reduction at the large-input sign boundary,
  3,935 exact comparative cases, 78 unit tests and an executed API example. Record
  three native timeouts separately; bound the corresponding Gram augmentation loop.
- All 56,382 matrix records, 56,727 focused and 7,096 fast tests pass (32 existing
  skips); builds pass; 584 type diagnostics unchanged. Coverage 104464/128162 (81.51%).
  Ideal routing/coercion and the wider audit remain open.

## 23.0.8 - 2026-09-12

- Port native word-prime pivots and matrix solves, modular integer rank with exact
  certificates, and batched CRT/rational lifting with exact solve certification.
- Add 5,418 exact bundled-PARI comparisons, 42 unit tests and an executed API example,
  covering singular reductions, later primes, zero-column RHS and rank retries.
- All 52,447 matrix records, 52,694 focused and 7,017 fast tests pass (32 existing
  skips); builds pass; 584 type diagnostics unchanged. Coverage 104068/127850 (81.40%).
  FLATTER/full LLL and ideal routing remain open.

## 23.0.7 - 2026-09-12

- Port PARI adaptive QR/Cholesky precision selection, exponent statistics and the
  upper-triangular shortcut, preserving retry schedules and native real records.
- Add 1,230 exact native comparisons, 32 unit tests and an executed API example.
  Clarify the previous product probe's measured maximum coefficient size (4,102 bits).
- All 47,029 matrix records, 47,265 focused and 6,974 fast tests pass (32 existing
  skips); builds pass; 584 type diagnostics unchanged. Coverage 103752/127539 (81.35%); closing
  delimiters remain counted. Rank/FLATTER/full LLL and ideal routing remain open.

## 23.0.6 - 2026-09-12

- Port native integer/real dot, Gram and matrix products, preserving summation order,
  integer backend delegation and real-to-integer rescaling/rounding boundaries.
- Add 5,903 exact native comparisons, 39 unit tests and an executed API example,
  including large multiplication cutoffs and all-zero real exponent overflow.
- All 45,799 matrix records, 46,041 focused and 6,941 fast tests pass (32 existing
  skips); builds pass; 584 type diagnostics unchanged. Changed files execute 117/117 lines.
  Coverage 103644/127430 (81.33%); adaptive Gram-Schmidt/FLATTER and ideal routing remain open.

## 23.0.5 - 2026-09-12

- Port real PARI positive Gaussian reduction, Cholesky and upper-triangular inversion,
  including exact-zero numeric branches and native error names/accuracy payloads.
- Add 10,445 exact native comparisons, 39 unit tests and an executed API example;
  keep 19 undefined large allocated-zero inverse probes outside passing counts.
- All 39,896 matrix records, 40,137 focused and 6,901 fast tests pass (32 existing
  skips); builds pass; 584 type diagnostics unchanged. Changed files execute 154/154 lines.
  Overall coverage 103529/127335 (81.30%); generic coefficients/FLATTER and ideal routing remain open.

## 23.0.4 - 2026-09-12

- Port native Householder QR, upper R and Gaussian reduction with mixed integer/real
  cells, exact precision state, public reflectors and precision-failure outputs.
- Add 4,752 exact native comparisons, 23 unit tests and an executed API example.
- All 29,451 matrix records, 29,675 focused and 6,861 fast tests pass (32 existing
  skips); builds pass; 584 type diagnostics unchanged. bibli1 executes 123/123 lines.
  Overall coverage 103460/127269 (81.29%); Cholesky/FLATTER and ideal-intersection routing remain open.

## 23.0.3 - 2026-09-12

- Port PARI's heuristic and arbitrary-precision proved LLL stages, safe real rounding
  and absolute real comparison, preserving precision, native coefficients and partial state.
- Add 13,338 exact native comparisons, 42 unit tests and an executed API example;
  retain 12 native/TypeScript timeout probes separately from passing comparisons.
- All 24,699 matrix records, 24,984 focused and 6,837 fast tests pass (32 existing
  skips); builds pass; 584 type diagnostics unchanged. Coverage 103370/127185 (81.28%) overall.
  FLATTER/full-wrapper and ideal-intersection routing remain under audit.

## 23.0.2 - 2026-09-12

- Port PARI's DPE LLL stage with exact Gram updates, supplied/incremental Gram modes,
  optional B/U/norms and partial failure state. Guard unrepresentable shifts explicitly.
- Add 3,588 exact native comparisons, nine separate resource-boundary checks, 26
  unit tests and an executed API example; retain 75 native timeout probes separately.
- All 11,361 matrix records, 11,603 focused and 6,794 fast tests pass (32 existing
  skips); build passes; 584 type diagnostics unchanged. Coverage 102999/126816 (81.22%) overall.
  Heuristic/arbitrary-precision/FLATTER wrapper and ideal-intersection routing stay open.

## 23.0.1 - 2026-09-12

- Port PARI's internal fast LLL/Babai stage, including partial failure transforms,
  optional U, keep-first and the audited native binary64 contraction profile.
- Add exact binary64 FMA/scaling/decomposition dependencies and 5,345 native
  comparisons, 17 unit tests and an executed API example. Retain 94 native timeout
  probes separately; full certified LLL and ideal-intersection routing remain open.
- All 7,764 matrix comparisons, 7,979 focused and 6,767 fast tests pass (32 existing
  skips). Build passes; 584 type diagnostics unchanged. New code executes 272/272
  lines; expanded instrumentation reports 102779/126594 (81.19%) overall. Audit remains open.

## 23.0.0 - 2026-09-12

- Breaking type correction: ideal denominator() now returns a cached integral ideal,
  matching Sage's coprime numerator/denominator decomposition. Fix fractional coprimality.
- Preserve rational ideal intersections, principal/PARI/monoid division dispatch and
  zero-product errors. Normalize empty/all-zero field ideal generators to a fresh zero ideal.
- Add 9,090 comparative cases, eight unit regressions and an executed API example;
  repair 1,595 prepatch wrong-behavior records. Keep native intersection LLL routing open.
- All 21,642 field comparisons pass; 85,064 focused, 6,749 fast (32 existing skips)
  and 417 caller tests pass. Build passes; 584 existing type diagnostics unchanged.
  Audit remains open at 102507/126010 (81.35%) instrumented line coverage.

## 22.0.0 - 2026-09-12

- Breaking type correction: NumberFieldIdeal.gens_two() returns two NumberFieldElement
  objects, preserving fractional intersections with QQ. Cache the pair and use PARI's
  certified generator reduction; rational ideals have a zero second generator.
- Build integral bases and free-module records from the maximal-order HNF, including
  fractional/nonprincipal and zero ideals. Cache free-module records like Sage.
- Add 1,878 native comparative cases, six unit regressions and an executed API example.
  Repair 910 wrong-behavior records and 742 missing-implementation records.
- All 12,552 field comparisons pass; 75,965 focused, 6,740 fast (32 existing skips)
  and 417 caller tests pass. Build passes; 584 existing type diagnostics unchanged.
  Audit remains open at 102483/126020 (81.32%) instrumented line coverage.

## 21.0.32 - 2026-09-12

- Delegate general fractional-ideal inversion, bounded HNF products and maximal-order
  differents to PARI trace/ideal dependencies, fixing shifted quadratic conjugation
  and nonmonogenic different ideals. Match ideal IntegerLike powers and zero errors.
- Port modular matrix adjoints/inversion, recursive half-GCD, rational lifting and
  integral triangular solves. Document the remaining higher-degree generator-search specialization.
- Add 3,574 native comparative cases, 12 unit regressions and an executed API example.
  Repair 278 wrong-behavior records and 588 missing-implementation records.
- All 10,674 field and 2,419 matrix comparisons pass; 76,499 focused, 6,733 fast
  (32 existing skips) and 417 caller tests pass. Build passes; 584 existing type
  diagnostics unchanged. Audit remains open at 81.30% instrumented line coverage.

## 21.0.31 - 2026-09-12

- Reduce number-field coefficient vectors through the new bundled-FLINT pseudo-division
  and rational-remainder ports. Preserve defining-polynomial scaling and degree-one roots.
- Match Sage integer-power scheduling and exponent-one object reuse; accept IntegerLike.
  Normalize arithmetic callers explicitly for quadratic units, different ideals and Frobenius filtering.
- Add 6,002 comparative cases, 12 regressions and an executed API-reference example.
  Repair 1,218 original discrepancies; 288 scaled-caller controls catch integration regressions.
- All 8,594 number-field and 600 new FLINT records match native results. 71,987 focused,
  6,720 fast (32 existing skips) and 417 caller tests pass; build passes and 584 type
  diagnostics are unchanged. Audit remains open at 81.12% instrumented line coverage.

## 21.0.30 - 2026-09-12

- Preserve Pollard parent resolution, native set_immutable calls and rho memory
  equality/eviction for colliding keys. Restore Sage's singleton lambda error and RNG state.
- Add 692 comparisons exposing 504 former discrepancies, including actual binary
  matrices, nine regressions and an executed API-reference assertion.
- All 15,454 group records, 65,972 focused, 6,707 fast (32 existing skips) and 355
  caller tests pass. Build passes; all 584 type diagnostics remain unchanged.

## 21.0.29 - 2026-09-12

- Resolve generic logarithm/order parent identities once and preserve them through
  nested BSGS, exponential bounds and order reduction. Use source equality tests.
- Add 3,870 parent/equality/error comparisons exposing 795 former discrepancies
  and seven regressions. All 14,762 group records match.
- All 65,270 focused, 6,697 fast tests (32 existing skips) and 355 caller tests pass.
  Build passes; all 584 existing type diagnostics remain unchanged.

## 21.0.28 - 2026-09-12

- Preserve Sage's BSGS, has_order, order_from_bounds and discrete-log operation
  schedules, including eager recursive powers and native negative-power calls.
- Resolve BSGS string-key collisions with element equality and preserve RDF's
  distinct ordinary and f-string error representations.
- Add 3,604 comparisons exposing 1,575 former discrepancies and ten regressions.
  All 10,892 group records match; 61,393 focused, 6,690 fast (32 existing skips)
  and 355 caller tests pass. Build passes; 584 type diagnostics are unchanged.

## 21.0.27 - 2026-09-12

- Coerce RDF arithmetic/equality scalar operands and preserve native integer/bool
  object-reuse shortcuts. Apply Sage's binary algorithm to standard group multiple.
- Add 2,922 comparisons exposing 1,789 former discrepancies, eight regressions and
  an executed example. All 278,897 polynomial and 7,288 group records match.
- All 57,779 focused, 6,680 fast tests (32 existing skips) and 355 callers pass.
  Build passes; four old type errors are removed and 584 remain unchanged.

## 21.0.26 - 2026-09-12

- Compare Rational with host numbers through RDF coercion, preserving fractional
  values, native real rounding, underflow, overflow, signed zero and NaN behavior.
- Add 1,170 comparisons (812 former discrepancies), five regressions and an
  executed example. All 7,472 rational records match.
- All 54,848 focused, 6,671 fast tests (32 existing skips) and 355 caller tests pass.
  Build passes; all 588 existing type diagnostics remain unchanged.

## 21.0.25 - 2026-09-12

- Route Integer-left rational arithmetic/comparisons through QQ kernels; retain
  Integer-only return types and expose rational result overloads.
- Recognize native __invert__ in generic groups and derive missing parent identities
  from neutral powers/actions. Type group equality against this and expose the
  rational completion of Integer multiplicative operations.
- Add 1,680 comparisons (333 former discrepancies), eight regressions and an executed
  example. All 8,230 Integer and 6,892 group records match; 47,370 focused, 6,665 fast
  tests (32 existing skips) and 355 caller tests pass. Build passes; all 588 existing
  type diagnostics remain unchanged.

## 21.0.24 - 2026-09-12

- Preserve number-field coefficient-index bounds, quadratic negative indices,
  float-index diagnostics and trimmed-coefficient behavior. Retain the original
  polynomial scaling and bounded discriminant normalization used to select the
  native quadratic representation.
- Add 1,638 comparisons (1,332 former discrepancies), eight regressions, and an
  executed example. All 3,192 number-field and 6,508 group records match.
- All 38,747 focused, 6,656 fast tests (32 existing skips), and 355 caller tests
  pass. Build passes; all 588 existing type diagnostics are unchanged.

## 21.0.23 - 2026-09-12

- Coerce exact scalars in number-field addition, subtraction, division and equality;
  accept Integer wrappers in field construction and convert real inputs through
  the existing simplest-rational kernel instead of truncating them.
- Preserve native inversion and scalar/field zero-division errors. Add 1,530
  comparisons (850 former discrepancies), eight regressions and an executed example.
- All 1,554 number-field and 6,508 group records, 37,100 focused tests, 6,647 fast
  tests (32 existing skips), and 355 caller tests pass. Build passes; all 588
  existing type diagnostics remain unchanged.

## 21.0.22 - 2026-09-12

- Restore exact integer/rational scalar multiplication of number-field elements;
  recognize Sage identity predicates in generic group algorithms.
- Add 821 native comparisons, six regressions, and an executed public example.
  These repair 746 formerly discrepant scalar/group records.
- All 6,508 group comparisons, 35,537 focused tests, 6,638 fast tests (32 existing
  skips), and 355 number-field/elliptic caller tests pass. Build passes; all 588
  existing type diagnostics remain unchanged.

## 21.0.21 - 2026-09-12

- Restore eager multiple-iterator validation/copying, advance-before-return, and
  recovery after callback errors; resolve both stored and callable group parents.
- Add 608 comparative cases (576 former discrepancies), eight unit regressions,
  and an executed iterator example. Document host copy and generator adapters.
- All 5,687 group comparisons, 34,709 focused tests, 6,631 fast tests (32 existing
  skips), and 212 elliptic caller tests pass. Build passes; all 588 existing type
  diagnostics are unchanged.

## 21.0.20 - 2026-09-12

- Treat empty group factor/prime lists as absent and preserve factor precedence.
  Delegate valuations to the existing Integer/GMP path; preserve terminating
  zero-multiple cases and native infinite-cost diagnostics.
- Add 480 native comparisons, including 63 former port hangs captured in isolated
  processes, five focused regressions and an executed empty-list example.
- All 5,079 group comparisons, 34,092 focused, 6,622 fast tests (32 existing skips)
  and 212 elliptic caller tests pass. Build passes; 588 existing type diagnostics
  are unchanged. Original nonterminating inputs remain recorded audit limits.

## 21.0.19 - 2026-09-12

- Preserve native custom-group validation text and precedence across group APIs.
- Restore the custom multiple operation schedule and idempotence shortcut;
  divide order bounds exactly and disable the redundant post-BSGS order check.
- Add 1,620 comparisons of errors, results and callback traces plus five focused
  regressions. All 4,599 group comparisons, 33,606 focused, 6,616 fast tests
  (32 existing skips) and 212 elliptic caller tests pass. Build passes and all
  588 existing type diagnostics are unchanged.

## 21.0.18 - 2026-09-12

- Restore Pollard rho's modular-ring random draws in native order, identity and
  retry behavior, inclusive small-order BSGS bound and exact error messages.
- Add 1,304 comparisons against the bundled original, checking results/errors
  and subsequent random state with the documented shared string hash. Add five
  focused regressions and an executed prime-order subgroup example; correct two
  unit expectations that asserted the former port-only errors.
- All 2,979 group-area comparisons, 31,981 focused and 6,611 fast tests pass
  (32 existing skips). Build passes; 588 existing type diagnostics are unchanged.

## 21.0.17 - 2026-09-12

- Restore Sage's CPython step-size stream, nonnegative signed-hash indices and
  exact invalid-bound errors in Pollard lambda.
- Add 1,340 live Sage comparisons of results, errors and complete walk traces,
  three focused regressions and an executed public example. All 1,675 group-area
  comparisons, 30,671 focused and 6,605 fast tests pass (32 existing skips).
- Build passes and all 588 existing type diagnostics are unchanged. Default
  object hashing and the documented width-zero workaround remain adapters.

## 21.0.16 - 2026-09-12

- Preserve PARI's nonfatal MPQS warning after Gauss and continue factor extraction.
- Add 240 warning comparisons, 232 class-candidate controls covering 7,642
  relations, and 11 driver factor-base early returns. Three focused regressions
  and an executed example protect warning text, count and continued results.
- Replay, fast, MPQS and caller suites pass; build passes and all 588 existing
  type diagnostics are unchanged. The broader audit remains open.

## 21.0.15 - 2026-09-12

- Preserve native MPQS relation hashing, duplicate checks, collision chains and
  incremental resizing so matrix columns, ordered factors and final RNG state
  match PARI.
- Add 94 seeded driver comparisons and 136 native hash/table controls, four
  focused regressions and an executed example. Replay, fast, MPQS and caller
  suites pass; all 65 instrumented lines of the new backend execute. Build
  passes and the 588 existing type diagnostics are unchanged.

## 21.0.14 - 2026-09-12

- Restore native MPQS word-inverse errors through the audited scalar backend.
- Add 1,093 native comparisons covering reduced inverses, class initialization,
  complete sieve bytes/candidates, factor-base branches, prime-table expansion,
  and the candidate capacity boundary. Four focused tests and an executed example
  pass. Replay, fast, MPQS and caller suites pass; build passes and the 588 existing
  type diagnostics are unchanged.

## 21.0.13 - 2026-09-12

- Delegate MPQS Kronecker symbols to the exact backend, preserving signed
  numerators and high word bits and terminating at power-of-two denominators.
- Restore native signed modular relation products, inverse powers and debug
  relation error types/messages. Add 4,802 native comparisons, seven regressions
  and an executed example. Replay, fast, MPQS and caller suites pass; build
  passes and the 588 existing type diagnostics are unchanged.

## 21.0.12 - 2026-09-12

- Restore native sparse binary kernel dispatch: singleton elimination and seeded
  block Lanczos above 640 rows, preserving retries, returned vectors and RNG state.
  MPQS delegates to the binary-matrix backend through its existing bitset adapter.
- Add 766 native comparisons, five focused regressions and an executed example.
  Replay, fast, MPQS and caller suites pass; build passes and the 588 existing
  type diagnostics are unchanged.

## 21.0.11 - 2026-09-12

- Preserve native negative k-th-power errors and avoid giant intermediate powers
  when the exact bit bound already determines the integer root.
- Use the magnitude of signed MPQS discriminants for relation inverses and the
  native canonical word square root for factor-base construction.
- Add 2,732 native comparisons, five focused regressions and an executed source
  example. Replay, fast, MPQS and caller suites pass; build passes and the 588
  existing type diagnostics are unchanged.

## 21.0.10 - 2026-09-12

- Replace serial ECM with native batched inversions, PRAC multiplication and
  helix/baby-step continuation. Preserve native failure handling, persistent
  state, prime/gcd schedules and exact large-seed increments.
- Add 800 native comparisons, six focused regressions and an executed source
  example. Replay, fast and affected slow suites pass; build passes and type
  diagnostics remain unchanged. All instrumented lines of the new module execute.

## 21.0.9 - 2026-09-12

- Restore signed perfect-power results, native word/multiword search order,
  residue-mask updates and the exponent-search failure sentinel.
- Add 8,914 native comparisons, four focused tests and executed public examples.
  Replay, fast and all slow unit files pass; build passes and type diagnostics
  remain unchanged. Keep the older magnitude-only factory adapter explicit.

## 21.0.8 - 2026-09-12

- Restore native SQUFOF/Pollard-Brent stage gates and exact Pollard-Brent budget
  counters, seed arithmetic and retry-exhaustion errors.
- Add 1,888 native factor-stage comparisons, 210 fixed-field prime-replacement
  comparisons, six focused tests and an executed source-subpath example.
  Replay, fast/caller suites and build pass; type diagnostics remain unchanged.

## 21.0.7 - 2026-09-12

- Delegate the Galois matrix kernel to native binary, packed ternary, word and
  generic-field backends. Preserve recursive echelon/triangular-solve schedules,
  native thresholds and signed generic-field kernel coefficients.
- Add 3,920 native comparisons through 144-by-144 matrices, three focused tests
  and an executed source-subpath example. Replay, fast/caller suites and build
  pass; type diagnostics remain unchanged.

## 21.0.6 - 2026-09-12

- Route the Galois word totient through the exact native word backend, preserving
  magnitude conversion and PARI's value 2 at zero.
- Add 4,838 comparisons covering Frobenius selection/random state, permutation
  filter decisions/caches, subgroup filtering, polynomial helpers and word totients.
  Replay, fast/caller suites, docs and build pass; type diagnostics are unchanged.

## 21.0.5 - 2026-09-12

- Preserve the native symmetric-polynomial search range with an exact BigInt
  limit/counter; JavaScript's 32-bit shift had prematurely exhausted the search.
- Add 557 native comparisons, three focused tests and an executed source-subpath
  helper example. Replay, fast/caller suites and build pass; type errors are unchanged.

## 21.0.4 - 2026-09-12

- Preserve PARI's positive safe-rounding boundary for fixed-field precision,
  repairing valid flag-two factorizations. Match fixed-field flag/permutation
  errors and malformed generator/order-pair validation in subgroup operations.
- Add 2,084 native comparisons using 24 identical Galois fixtures and four focused
  tests. Replay, fast/caller suites and build pass; type diagnostics are unchanged.

## 21.0.3 - 2026-09-12

- Match native Galois validation order, error classes/messages, generic linear
  conjugates and the flag-four identity fallback. Return a non-null conjugate list.
- Add 322 native comparisons, four focused tests and an executed API example.
  Replay, fast/caller suites and build pass; type diagnostics remain unchanged.

## 21.0.2 - 2026-09-12

- Complete scalar-bound and unbounded subgroup enumeration through the shared
  Birkhoff engine. Port strict partial factorization, including native word and
  multiword power rules, cache boundaries and unresolved cofactors; delegate the
  existing default helper to it.
- Add 2,781 native comparisons and an executed source-subpath example. Validate
  406,536 replay checks, fast/caller suites, build and unchanged type diagnostics.

## 21.0.1 - 2026-09-12

- Fix 555 unit-subgroup failures: 543 native ordering differences and 12 incorrect
  subgroup selections for signed order arguments. Restore the cyclic-factor,
  exact-index Birkhoff, modular-HNF and coset dependency route.
- Add 3,505 native comparisons covering decomposition, subgroup matrices, word
  boundaries and modular HNF. Validate 403,755 replay checks, fast/caller suites,
  build and unchanged type diagnostics. Remaining bound modes are explicit.

## 21.0.0 - 2026-09-12

- Fix 586 symmetric-polynomial evaluation failures: combine reduced Newton sums
  over the integers, preserve scalar zero for zero weights, and delegate unsigned
  powers to PARI's logarithmic schedules and native small-exponent shortcuts.
- Breaking API correction: sympol_eval returns bigint[] | bigint. Callers with
  all-zero weights now receive 0n and must narrow before array operations.
- Add 2,571 native comparisons and executed API examples; validate 400,249 replay
  checks, fast/caller suites, build and unchanged type diagnostics.

## 20.27.13 - 2026-09-12

- Preserve PARI's exception class and exact message in unsupported quotient
  construction and permutation paths, fixing six comparative failures.
- Add 3,515 native comparisons covering subgroup presentations, ordering, cosets,
  quotient maps/lifts and the A4/S4/F36 exceptional branches. Validate 397,674
  replay checks, fast/caller suites, build and unchanged type diagnostics.

## 20.27.12 - 2026-09-12

- Fix 1,926 permutation/vector failures: checked word exponents, exact order
  intermediates, three-way prefix comparison and native trivial-group errors.
  Restore counting-sort dispatch and enforce zv_prod's no-overflow precondition.
- Add 10,696 direct native comparisons and eight explicit guard checks; validate
  394,163 replay checks, fast and arithmetic/Galois suites, native suites, build
  and unchanged type diagnostics.

## 20.27.11 - 2026-09-12

- Fix 1,792 Vandermonde error discrepancies by porting product/remainder trees,
  native batch inversion and denominator ordering; preserve public matrix storage.
- Add 4,227 native comparisons and focused/API regressions. Validate 383,456 replay
  checks, fast and arithmetic/Galois suites, native suites, build and unchanged
  type diagnostics; share the existing balanced-product tree scheme.

## 20.27.10 - 2026-09-12

- Fix 2,234 reproduced Hensel-lifting failures: use native monic normalization,
  factor trees, precision doubling, Bezout propagation and root/inverse updates.
- Add 13,016 native comparisons and 42 explicit precondition-guard checks; preserve
  existing one-indexed public vectors. Validate 379,231 replay checks, fast and
  arithmetic/Galois suites, native suites, build and unchanged type diagnostics.

## 20.27.9 - 2026-09-12

- Fix 855 reproduced integer-polynomial failures: normalize trailing zeros,
  recognize squarefree constants, and replace indexpartial's approximation with
  strict partial factorization and native p-adic reduced-resultant refinement.
- Port modular integer GCD and column echelon dependencies; add 8,161 native
  comparisons and focused/API regressions. Validate 366,171 replay checks,
  arithmetic/Galois callers, fast and native suites, and build.

## 20.27.8 - 2026-09-12

- Fix 242 reproduced Galois integer-helper failures, including 27 time-limit
  cases: preserve word overflow, delegate factorization, port binary integer logs
  and use exact prime-iterator state with native cache/sieve/successor phases.
- Add 1,291 native comparisons plus focused/API regressions; validate 358,007
  replay checks and arithmetic/Galois caller suites.

## 20.27.7 - 2026-09-12

- Fix 3,771 reproduced rational Galois polynomial discrepancies in denominator
  normalization, per-coefficient finite-field conversion, centering and errors.
- Add 7,003 original-PARI comparisons plus focused/API regressions; validate
  356,715 replay checks and complete arithmetic/Galois caller suites.

## 20.27.6 - 2026-09-12

- Fix 5,872 reproduced predicate/order discrepancies: native Kronecker square
  testing, word-bound defaults, generic validation/raw bases, coprimality errors
  and zero Bezout coefficients. Route word-order factorization through factoru.
- Add 28,782 native comparisons and focused/API regressions; validate 349,769
  replay checks plus arithmetic and extended elliptic-curve callers.

## 20.27.5 - 2026-09-11

- Fix 479 reproduced scalar exponentiation discrepancies in zero powers, inverse
  dispatch and signed reductions; port native Montgomery/Barrett selection and
  window/fused schedules. Check a negative-modulus stack-exhaustion boundary.
- Add 12,137 native comparisons and focused/API regressions; validate 320,985
  replay checks plus arithmetic and extended elliptic-curve callers.

## 20.27.4 - 2026-09-11

- Fix 6,460 reproduced scalar arithmetic discrepancies, including incorrect
  products and inverses of negative operands and non-native error payloads.
- Restore native word-divisor/exact-division branches, zero shortcuts and the
  distinct squaring remainder rule; preserve single-comparison centering/doubling.
- Add 21,712 native comparisons and focused/API regressions; validate 308,848
  replay checks plus arithmetic and extended elliptic-curve callers.

## 20.27.3 - 2026-09-11

- Fix 3,379 reproduced derivative, evaluation, centering and linear-quotient
  discrepancies at signed/zero-modulus and raw coefficient-storage boundaries.
- Restore native sparse evaluation, cancellation and constant shortcuts; omit
  unrequested remainders and preserve native centering output storage.
- Add 10,703 shared comparisons and focused/API regressions; all 287,077 final
  replay checks pass, along with fast tests, Galois callers and native suites.

## 20.27.2 - 2026-09-11

- Fix signed power-table and composition routing, including cached trace contexts
  and native matrix empty/zero-product boundaries previously rejected by a guard.
- Preserve binary matrix dispatch for negative two, signed word magnitudes and
  conversion of supplied reciprocals; correct zero-polynomial oracle guard order.
- Add 44,838 native comparisons, refresh 18 guard checks against C, and execute
  six regressions/API examples; all 276,374 final replay checks pass.

## 20.27.1 - 2026-09-11

- Fix signed-word polynomial division, GCD and half-GCD backend selection,
  resolving 1,902 observed discrepancies while preserving earlier shortcuts.
- Route multiplication and squaring by native unsigned magnitude, and frame
  native oracle errors as JSON without changing their text or exception classes.
- Add 28,490 native comparisons and five focused regression tests; validate
  231,516 replay checks, fast tests, arithmetic callers and native suites.

## 20.27.0 - 2026-09-11

- Add nine extension split-part, root-count, squarefree and derivative APIs using
  the existing native Frobenius/GCD backends and coefficient-tag conventions.
- Fix scalar modular inverse/residue behavior on signed and zero moduli,
  resolving 274 observed discrepancies and the affected prime root-count caller.
- Add 15,230 permanent native comparisons, focused regressions and executed API
  examples; validate 203,027 replay checks and shared arithmetic callers.

## 20.26.1 - 2026-09-11

- Fix prime-polynomial linear operations at signed/zero modulus boundaries,
  resolving 5,280 observed discrepancies in results, shortcuts and exceptions.
- Preserve native cancellation and zero-scalar returns before division; remove
  trailing zero storage before arithmetic and retain exact dvmdii diagnostics.
- Add 18,414 shared native comparisons, regressions and executed examples;
  continue extension factor dependencies and full behavioral coverage.

## 20.26.0 - 2026-09-11

- Add six native prime and extension Frobenius APIs using the original
  composition/powering threshold and shared reciprocal caches.
- Preserve characteristic-two shortcuts, unsigned-word half-Frobenius conversion,
  coefficient tags, constant quotients and native error phases.
- Add 6,195 shared comparisons, cache/branch regressions and executed API examples;
  continue extension factorization and full source/behavioral coverage.

## 20.25.1 - 2026-09-11

- Fix prime quotient power reciprocal preparation, inverse error payloads and
  signed/unsigned word backend selection, resolving 1,753 observed discrepancies.
- Add 16,228 native comparisons and strengthen 2,820 older comparisons to check
  exact errors; include cache reuse and large-payload regressions.
- Validate 163,185 comparative/regression/API cases, fast tests and arithmetic
  callers; continue Frobenius, factor routing and full port coverage.

## 20.25.0 - 2026-09-11

- Add generic and word extension minimal polynomials using native Shoup
  projection, transposed products, half-GCD and shared power tables.
- Preserve random-state consumption, repeated-factor multiplicity, coefficient
  tags and outer-before-inner reciprocal preparation.
- Add 5,695 shared comparisons, cache/error regressions and executed examples;
  continue quotient repairs, Frobenius, factor routing and full coverage work.

## 20.24.0 - 2026-09-11

- Add native extension random projections, delayed-reduction dot products and
  generic/word truncated products and squares.
- Preserve random-state consumption, coefficient tags, integer-only square
  coefficients and backend-specific truncation/error ordering.
- Add 20,281 shared comparisons, regressions and executed examples for all eight
  APIs; continue extension minimal polynomials and full port coverage.

## 20.23.0 - 2026-09-11

- Add eight native extension automorphism power, additive trace and multiplicative
  aggregate APIs with the original tuple layouts and window powering schedules.
- Preserve signed exponent casting, unit-count reciprocal preparation, raw
  coefficient tags and shared inner/outer caches.
- Add 10,686 shared native comparisons, cache/count regressions and executed API
  examples; continue extension minimal-polynomial and full coverage work.

## 20.22.0 - 2026-09-11

- Add six native coefficient-substitution APIs and binary scalar Brent–Kung
  composition, preserving scalar tags, raw power tables and native cache phases.
- Preserve safe table guards after native power preparation and earlier
  coefficient errors; retain prepared inner reciprocals for internal callers.
- Add 11,437 shared native comparisons, cache/schedule regressions and executed
  examples for the eight new exports.
- Continue extension automorphism, minimal-polynomial and full port coverage work.

## 20.21.2 - 2026-09-11

- Fix native reciprocal preparation and canonical storage in automorphism powers,
  successive automorphism tables and additive traces.
- Preserve generic-to-word cache reuse, the direct word zero shortcut and native
  count-zero/count-one remainder ordering.
- Add 7,161 permanent comparisons, including all 1,532 observed discrepancies;
  retain exact table-reuse assertions and add cache-routing regressions.
- Update executed API examples and continue the full port audit.

## 20.21.1 - 2026-09-11

- Fix native reciprocal preparation in public prime-field power tables and
  composition, including inverse errors before zero/unit and constant results.
- Canonicalize trailing zero storage before degree checks and raw first-power
  copies; preserve unreduced coefficient values and zero-polynomial shortcuts.
- Route word power tables through Flx products and reuse prepared reducers.
- Add 8,130 shared comparisons, regression/routing tests and executed examples.
- Continue the full port audit and extension-factor dependencies.

## 20.21.0 - 2026-09-11

- Add native generic, word and binary extension composition, preserving
  Brent–Kung blocks, raw tables, coefficient tags and source-order errors.
- Add packed word-extension matrix multiplication through the existing native
  integer-matrix backend; retain the shared inner reciprocal during evaluation.
- Add 22,386 comparative cases, delegation/cache tests and executed API examples.
- Continue coefficient substitution, Frobenius, minimal-polynomial and full
  extension-factor dependencies; full audit coverage remains open.

## 20.20.0 - 2026-09-11

- Add 21 native PARI extension quotient operations for products, inverses,
  powers and power tables, with reusable inner and outer reciprocal caches.
- Fix long polynomial and modular-integer inverse errors to match native display;
  preserve multiline oracle errors without corrupting subsequent comparisons.
- Add 24,023 shared native comparisons and executed API examples; preserve the
  bundled generic word-power cache behavior and guard its unsafe input domain.
- Continue composition, Frobenius and extension-factor dependencies; full coverage
  remains open.

## 20.19.0 - 2026-09-11

- Add native PARI extension GCD, Bézout coefficients and recursive half-GCD,
  preserving unscaled outputs, coefficient tags and the native matrix schedules.
- Fix binary packing for raw coefficients within native word allocations,
  including overlapping blocks and XOR cancellation.
- Add 24,532 shared comparisons, regression/routing tests and executed API examples;
  isolate the bundled binary GCD cycle in bounded termination comparisons.
- Continue extension quotient and factor dependencies; full coverage remains open.

## 20.18.0 - 2026-09-11

- Add native PARI extension-polynomial division, cached reductions and Newton
  reciprocals, preserving coefficient tags, algorithm cutoffs and native errors.
- Fix the binary product guard at the modulus-degree boundary; preserve raw
  word/binary reciprocal terms and linear copying during Barrett traversal.
- Add 27,367 shared native comparisons and regression/routing tests; document
  reciprocal containers and the bundled binary remainder behavior.
- Continue the extension GCD/factor audit; full behavioral coverage remains open.

## 20.17.0 - 2026-09-11

- Add native PARI extension-polynomial products, squares, coefficient reduction
  and normalization as dependencies for the extension-factor backend.
- Preserve integer/polynomial coefficient tags, native packing and scalar paths,
  distinct word/binary normalization and exact inverse-error payloads.
- Add 10,575 shared comparisons and routing/API tests. Node and Bun comparisons
  catch and prevent argument-spread overflow for large coefficient arrays.
- Document representations and guards; full extension-factor routing remains open.

## 20.16.0 - 2026-09-11

- Restore PARI deterministic polynomial roots, word-field cutoffs and native
  root ordering; fix zero-polynomial errors and raw-degree splitting predicates.
- Add word root/count/predicate entry points, share the native nonsquare helper,
  and preserve quadratic and binary root-count shortcuts.
- Add 17,709 shared comparisons, routing checks and executed API examples;
  document native ordering quirks, guards and remaining extension-factor work.

## 20.15.0 - 2026-09-11

- Restore PARI full polynomial factorization, native equal-degree splitting and
  random-state consumption; delegate large-prime Sage and Galois factor callers.
- Fix normalization errors and modulus-one leading coefficients, and route word
  squarefree arithmetic through Flx. Add 19,566 shared native/Sage comparisons.
- Reduce packing conversion costs with hexadecimal strings for aligned blocks;
  retain native packed values, algorithm cutoffs and the existing memory budget.
- Document the new factor/normalization APIs with executed examples and track
  remaining root and extension-factor audit work.

## 20.14.0 - 2026-09-11

- Restore PARI modular square-root word dispatch, signed-small/Gauss-sum
  shortcuts, Cipolla and native composite/zero boundaries; eliminate eight
  observed nonresidue-search hangs.
- Add 21,551 shared native comparisons, fused powering and supplied-generator
  entry points, routing regressions and executed API documentation.

## 20.13.0 - 2026-09-11

- Restore PARI Shoup distinct-degree factorization, word-factor scales, native
  factor counts, binary dispatch and reusable quotient reduction contexts.
- Add 4,165 shared comparisons and caller/dispatch regressions; move packed
  products into recursive multiplication earlier to reduce large-case runtime.

## 20.12.0 - 2026-09-11

- Replace repeated linear-algebra minimal-polynomial searches with PARI's Shoup
  projections, recursive half-GCD reconstruction and shared Barrett reductions.
- Preserve exact random-state consumption and word dispatch; add 6,203 shared
  comparisons, native-contract guards, routing checks and executed API examples.

## 20.11.0 - 2026-09-11

- Restore PARI unscaled GCD/Bézout results, recursive half-GCD and native cutoffs;
  repair dependent inverse, Hensel and split-part normalization.
- Restore arbitrary-exponent quotient power schedules and signed-word dispatch;
  add 10,613 shared native comparisons, routing tests and executed API examples.

## 20.10.0 - 2026-09-11

- Restore PARI division/remainder early returns, error payloads, word dispatch,
  reciprocal Newton inversion and Barrett reduction with native thresholds.
- Add 8,541 shared comparisons, including 630 prior discrepancies, modulus-one
  inverse representation, degree collapse and routing/cutoff regressions.

## 20.9.0 - 2026-09-11

- Restore PARI word/integer polynomial multiplication and separate square kernels,
  native packing thresholds and signed coefficient recovery; route existing callers.
- Bound packed intermediates with recursive products, fixing three engine-limit
  failures; add 2,710 bundled-native comparisons and backend/cutoff regressions.

## 20.8.0 - 2026-09-11

- Restore PARI Brent–Kung modular composition, reusable power tables and native
  automorphism trace schedules; preserve word-table conversion boundaries.
- Canonically reduce signed polynomial/scalar products and guard undefined native
  vector-copy contracts; add 15,018 shared comparisons and routing regressions.

## 20.7.0 - 2026-09-11

- Restore PARI integer matrix multiplication thresholds, Winograd products and
  modular reconstruction; share binary/word/arbitrary-modulus matrix backends.
- Add 1,269 direct native comparisons, a modular-routing regression and executed
  column-array API examples, including empty shapes and valid modulus one.

## 20.6.1 - 2026-09-11

- Preserve PARI's unreduced initial quotient-power and automorphism-vector entries.
- Replace linear automorphism powering with native binary composition order and
  guard unrepresentable counts; add 3,775 direct bundled-PARI comparisons.

## 20.6.0 - 2026-09-11

- Restore binary finite-field polynomial factor delegation to PARI, including
  packed arithmetic, native kernel bases, Cantor/Berlekamp dispatch and RNG state.
- Add 11,139 shared native/Sage comparisons, routing regressions and executed API
  examples; correct the old documentation attributing binary factoring to NTL.

## 20.5.0 - 2026-09-11

- Port PARI's native 64-bit XORGEN state, save/restore, rejection sampling and
  polynomial/vector sampling needed by the remaining finite-field factor routes.
- Add 2,408 shared comparisons against an isolated bundled-PARI executable,
  automatic temporary source builds, native-contract guards and executed API docs.

## 20.4.0 - 2026-09-11

- Restore word-prime polynomial factor delegation to FLINT, including deflation,
  native CZ/KS dispatch, equal-degree splits and baby/giant distinct-degree splitting.
- Port modular Brent–Kung composition, matrix multiplication and row evaluation;
  compare factor ordering and exact random-state consumption against bundled C.
- Add 5,359 shared comparisons, backend-routing regressions, temporary native-build
  tooling, safe native-contract guards and executed public API documentation.

## 20.3.0 - 2026-09-11

- Port FLINT's 64-bit word random state and sampling helpers required by native
  polynomial factorization, preserving exact output streams and seed consumption.
- Add 3,036 shared comparisons against compiled bundled C routines, portable
  invalid-shift guards, unit regressions and an executed public API example.

## 20.2.4 - 2026-09-11

- Restore variable-map source conversion, identity-map coercion and function-field
  singleton-list delegation to the fraction-field constructor.
- Restore polynomial coefficient dictionaries with backend-specific coefficient,
  exponent, validation-order and overflow behavior; retain documented FLINT guards.
- Add 2,134 shared comparisons, focused regression tests and executed API examples.

## 20.2.3 - 2026-09-11

- Restore rational function-field matrix base arguments, native parent-cache warming
  and argument hashability errors without caching the returned matrix arrays.
- Add 2,205 shared comparisons for cold/warm fields, failed calls, raw fractions,
  scalar/container bases and cache isolation, plus executed API documentation.

## 20.2.2 - 2026-09-11

- Restore function-field predicate parent/category dispatch and native error boundaries,
  including scalar wrappers, generic matrices and dynamically supplied methods.
- Add 426 shared comparisons for predicate values, exceptions and attribute access counts.

## 20.2.1 - 2026-09-11

- Restore native function-field matrix, trace and norm operations on unreduced
  fractions, preserving normalization and element identity behavior.
- Return immutable matrix arrays with native assignment errors and a mutable-copy adapter.
- Add 6,892 shared comparisons for raw matrices, divisors, evaluations and copies.

## 20.2.0 - 2026-09-11

- Restore lazy polynomial and place enumeration, including native iterator order,
  empty products, degree bounds and zero-degree place assertion errors.
- Add inherited QQ cardinality and correct backend-specific negative monomials;
  document safe exceptions for undefined native FLINT inputs.
- Correct infinite-cardinality handling in hyperelliptic callers, reject invalid
  low-degree projective models and restore balanced prime-field multivariate display.
- Add 7,136 shared comparisons covering large-field prefixes, iterator restart/copy
  behavior, metadata failures, rational enumeration and curve/polynomial callers.

## 20.1.0 - 2026-09-11

- Restore ideal parents, cached monoid identities and the native binary power
  schedule, retaining exponent-one objects and raw generators.
- Coerce scalar membership inputs and preserve native zero-ideal division errors.
- Add 11,220 shared ideal comparisons and executed API examples, including
  raw/normalized ideals, factors, valuations and powers beyond machine-word size.

## 20.0.0 - 2026-09-11

- Restore divisor parent/zero identity, place conversion, zero-support cleanup and
  cached Riemann-Roch maps with native vector and coordinate errors.
- Correct `_format` to take Sage's formatter argument and preserve dependent/zero
  echelon rows, including original exception behavior.
- Add 5,640 shared comparisons, vector-modulus boundary cases and executed examples.

## 19.1.0 - 2026-09-11

- Restore place-set/valuation-ring identities, residue-map caching and domain checks,
  and existing-element identity during function-field coercion.
- Preserve native QQ residue reduction through polynomial coefficient numerator
  clearing; add polynomial numerator, denominator and LCM methods with FLINT routing.
- Add 6,278 comparisons for places, maps, nested coefficients and LCM kernels, plus
  executed API examples and source documentation.

## 19.0.0 - 2026-09-11

- Correct PARI `FpX_factor_squarefree` to return native multiplicity-indexed
  components, including repeated factors, placeholders, units and constant dispatch.
- Preserve Galois irreducible-factor behavior under a separate support helper.
- Add 844 native comparisons exposing 781 old discrepancies, public examples and
  a defensive guard for unsafe native large-prime zero input.

## 18.3.0 - 2026-09-11

- Restore squarefree units, factor order, zero/composite errors and backend dispatch;
  port FLINT squarefree factors and native field/Frobenius algorithms.
- Add dedicated GF2 coefficient square roots and native canonical-element caching,
  preserving directly constructed zero negation and existing-element conversion.
- Add 6,787 shared comparisons and executed API examples; document the bundled
  zero guard and explicit native-input adapters.

## 18.2.0 - 2026-09-11

- Route fraction and function-field roots through their native generic/FpT backends;
  preserve unreduced squares, zero identity, root choices, lists and errors.
- Port FLINT polynomial roots, basecase/Newton series and high-product dependencies.
- Preserve integer polynomial root parents through NTL squarefree decomposition,
  modular GCD, balanced CRT and exact-division certification.
- Add 50,082 shared native comparisons, root-parent regressions and executed API
  examples; document dense-kernel contracts and remaining named-extension limits.

## 18.1.0 - 2026-09-11

- Port FLINT word modular square roots and their Jacobi, square-check, modular-power
  and limb-preinversion dependencies, preserving root choices and native branches.
- Add 69,725 direct native comparisons, word-boundary regressions and executed
  public examples; document BigInt arithmetic and native-precondition adaptations.

## 18.0.0 - 2026-09-11

- Match function-field modular inversion assertions and preserve valid inverses
  modulo zero, including nonmonic native constant gcds.
- Preserve order construction errors, conversion between finite and infinite ideals,
  unique ideal monoids and identity when constructing an existing monoid element.
- Return the bundled infinite-order integer basis `[1n]`; document and test the
  installed Sage runtime's different implementation.
- Add 3,982 comparative cases plus focused regressions and executed public examples.

## 17.0.0 - 2026-09-11

- Restore actual fraction parents and elements from function-field `field()` and
  `element()`; retain identity, lazy normalization and native arithmetic dispatch.
- Match cross-parent conversion, unreduced fraction predicates and display, generic
  negation, and power identity without unwanted normalization.
- Restore fraction string fallback, denominator expressions, Boolean conversion,
  target-parent identity and syntax-before-name errors; use the bundled preparser
  in comparative tests where installed Sage's backslash behavior differs.
- Add 3,745 permanent comparisons, focused regressions and executed public examples.

## 16.0.1 - 2026-09-11

- Prevent rational-function valuation from looping on constants; convert complete
  fraction arguments and preserve divisor validation for zero inputs.
- Propagate bundled Sage polynomial zero-factorization errors through function
  fields and finite ideals; retain native infinite-ideal degree behavior.
- Add 1,023 comparisons, including 188 formerly looping input cases and 196 other
  discrepancies, with explicit compatibility for the older Sage runtime's missing guard.

## 16.0.0 - 2026-09-11

- Restore rational function-field factory identity and the field-plus-two-maps
  result of variable renaming. Preserve field and order generator `IndexError`.
- Match native fraction ordering, numerical rational coefficient comparison,
  characteristic-two square-root lists and backend-specific root/zero errors.
- Add 3,648 comparative cases exposing 415 old discrepancies, focused regressions
  and an executed public renaming example. Track remaining dependency/coercion work.

## 15.0.0 - 2026-09-11

- Preserve BigInt exponents and native overflow in direct `exp1r_abs`; keep native
  mantissa multiplication and range-reduced exponential behavior.
- Add 1,040 original-source comparisons exposing 270 old discrepancies, and focused
  regressions. Cross-check the oracle against an isolated build of bundled PARI.
- Document the native series allocation gap separately from valid-input equivalence.

## 14.0.0 - 2026-09-11

- Route Buchmann square roots and transcendentals through native kernels; restore
  complex roots, zero accuracy and exponential precision growth.
- Preserve BigInt exponents in public `mpexp` results and accept them in `rtodbl`.
- Restore native progressive-precision logarithm series, final precision clamping
  and binary splitting for log(2); add direct native helper comparisons.
- Add 8,061 original comparison cases, including 1,055 old result mismatches, and
  executed public examples; update API, architecture, deviation and audit records
  and synchronize workspace versions.

## 13.0.0 - 2026-09-11

- Replace independent Buchmann elementary real arithmetic with shared native kernels,
  including constructors, precision conversion, cancellation and truncated products.
- Restore the native one-argument `real_0_bit` signature; round requested working bits
  up to whole words and accept BigInt for unsigned-word multiplication/division.
- Restore the native `truncr` precision guard while preserving `gcvtoi` error estimates;
  route algebraic dependencies through PARI's distinct scaled truncation helper.
- Add 70,678 original comparison cases, focused regressions and public examples;
  update architecture/deviation records and synchronize workspace versions.

## 12.0.2 - 2026-09-11

- Route Buchmann binary64 conversions and real comparisons through native kernels;
  restore zero accuracy, whole-word working precision, tie rounding and overflow errors.
- Compare signs/exponents before mantissas, preserving approximate-zero semantics and
  avoiding allocations proportional to exponent gaps.
- Export the native `PariError` and preserve Sage's resultant fallback at that boundary.
- Add 16,728 direct comparative cases exposing 5,489 old mismatches, focused regressions
  and an executed error-class example; synchronize workspace versions.

## 12.0.1 - 2026-09-11

- Restore native `gcvtoi` error exponents for exact cancellation and allocated/canonical
  zeros in both real helpers; preserve the source validation and precision order.
- Add 4,012 original comparisons exposing 798 old mismatches, focused regressions and
  an executed public example; synchronize workspace versions.

## 12.0.0 - 2026-09-11

- **Breaking:** PARI `sqrtr` returns `MpReal | MpComplex`; negative roots preserve an
  exact zero real component. Low-level `sqrti` reads magnitude for negative integers.
- Restore native square-root guard rounding and full zero-exponent halving; delegate
  integer root/remainder to GMP's normalized Karatsuba algorithm.
- Add 6,276 shared native cases, including 66,048 hashed root/remainder values, exposing
  874 old mismatches. Update component types, executed examples and workspace versions.

## 11.0.3 - 2026-09-11

- Restore PARI's truncated real multiplication/squaring, native integer conversion
  precision and guard rounding, including pointer-identity kernel selection.
- Add 21,466 native comparisons exposing 3,404 old mismatches, focused regressions
  and an executed identity-dispatch example.
- Update stale real-arithmetic documentation and synchronize workspace versions.

## 11.0.2 - 2026-09-11

- Preserve native PARI zero allocation and accuracy through integer/real conversions;
  restore allocation-based zero division and inverse dispatch.
- Add 2,052 native comparisons exposing 1,089 old mismatches, plus focused regressions
  and an executed public example. Correct the oracle's zero serializer, revealing
  90 previously hidden rational-conversion mismatches.
- Document unspecified native zero payloads and synchronize workspace versions.

## 11.0.1 - 2026-09-11

- Restore native PARI real division, quotient correction, denominator truncation and
  guard rounding across single-word, short Knuth, GMP and reciprocal paths.
- Accept full 64-bit bigint divisors in `divru`; document the native unsigned domain.
- Add 11,374 native comparisons, exposing 943 old mismatches in the initial 9,784
  cases, plus colocated regressions and an executed public example.
- Correct stale dependency documentation and synchronize workspace versions.

## 11.0.0 - 2026-09-11

- **Breaking:** `factorial(n, 'pari')` now returns the original rounded real result as
  `{s,e,m,p}`, with bigint mantissa/exponent; default/GMP factorial still returns bigint.
- Port PARI's exact and precision-aware factorial products, positive-integer gamma,
  Bernoulli cache, pi, exponential and reciprocal Newton branches, and AGM logarithm.
- Preserve native signed-word/exponent errors, tuning thresholds and product/power
  operation order. Add 10,085 original comparisons through 8,192-bit precision;
  demonstrate 1,012 old free-PARI mismatches in 1,032 safe calls.
- Document the native real representation, internal C domains and precision-context
  limits; synchronize public examples and package versions.

## 10.0.2 - 2026-09-11

- Restore PARI real addition/subtraction word precision, carry handling, cancellation
  and zero accuracy; select the original precision for mixed integer operands.
- Add 57,825 native comparisons, including the final guard-rounding carry, with
  23,993 old mismatches demonstrated in the initial 57,820 cases.
- Document the deterministic zero guard used for an upstream out-of-range padding read.

## 10.0.1 - 2026-09-11

- Route Integer and free GMP/default factorials to GMP's split-product/prime-swing backend;
  restore exact negative-input, unsigned-word overflow and unknown-selector errors.
- Port native odd-factorial flags and the blocked prime sieve, retaining 64-bit thresholds
  and masks; add 3,098 native comparisons through 60,000! and five million sieve entries,
  exposing 197 mismatches in 984 safe old calls.
- Document runtime allocation boundaries and the remaining free-factorial/PARI backend gap.

## 10.0.0 - 2026-09-10

- **Breaking:** free integer valuation now returns `bigint | 'Infinity'`, restoring Sage's
  zero-before-base-validation behavior instead of throwing. Generic method dispatch now
  supports Rational values and preserves AttributeError-to-Integer fallback.
- Restore exact base error messages and native GMP factor extraction in Integer/Rational
  valuation and unit paths; accept IntegerLike bases in Rational convenience aliases.
- Cache immutable Integer magnitude metadata at construction and restore the nbits alias
  direction, replacing repeated shifts with constant-time bit-length queries.
- Add 32,307 original comparisons for scalar/protocol dispatch, signed native units, large
  factor powers and bit lengths through 262,144 bits; 3,364 old mismatches are demonstrated.

## 9.1.0 - 2026-09-10

- Restore complex display ordering for tiny imaginary parts, native binary rounding,
  nonfinite values, tuple records and empty-list identity. Keep MPFR keys to distinguish
  rounded finite values from infinity at the binary64 limit.
- Add a CPython-compatible stable powersort adapter, retaining original comparison order
  for NaN keys with natural runs, balanced merging and galloping; preserve its license.
- Add 3,120 original display comparisons exposing 2,969 previous mismatches, and 3,248
  direct sort result/trace comparisons through 16,384 entries.

## 9.0.1 - 2026-09-10

- Preserve exact rational division in maximal-quotient reconstruction, matching the
  executable original rather than its stale floor-division comment.
- Normalize IntegerLike inputs before zero/threshold handling; retain source selection
  and sign behavior. Replace permissive reconstruction unit tests with exact outputs.
- Add 44,332 original comparisons, exposing 18,629 previous mismatches, including 4,973
  within the documented input domain and boundaries through 1024 bits.

## 9.0.0 - 2026-09-10

- **Breaking:** two/three/four/k-square decompositions now raise Sage's ValueError for
  impossible representations instead of returning null; return types are nonnullable.
- Normalize IntegerLike inputs before branching, fixing wrapped-zero loops and missed
  perfect-square shortcuts. Support exact and wrapped square counts while retaining
  numeric counts with Python truncation and nonfinite conversion errors.
- Add the mirrored native sum_of_squares module with uint32 bounds and exact tuple order;
  restore low-bit extraction for large powers of four.
- Add 11,616 original comparisons. The safe old-input subset exposes 5,415 mismatches;
  former nonterminating inputs and all eight supported entry points have regressions.

## 8.18.0 - 2026-09-10

- Restore native PARI totient and divisor-count delegation, with distinct native/Sage
  zero conventions and original prime-power and balanced-product paths.
- Share PARI's divide-and-conquer valuation extraction in factorization and quadratic
  roots instead of repeated-division copies.
- Add 2,490 original comparisons for native arithmetic, signed/word boundaries,
  high prime-power exponents and direct C quadratic-root results.

## 8.17.0 - 2026-09-10

- Restore the original extended-LCM gcd reduction, signed multiple/factors and zero
  division behavior; accept IntegerLike inputs consistently.
- Preserve negative Dedekind psi values and the zero factorization error; normalize
  its and Carmichael lambda's IntegerLike boundaries.
- Add 19,624 original comparisons for extended LCM and six factor-based scalar APIs;
  13,172 old discrepancies are fixed. Reuse the documented source-version adapter for
  radical(0), where the current port already matches the bundled implementation.

## 8.16.0 - 2026-09-10

- Add integer ProductTree layers, remainders, cached CRT interpolation and balanced
  streaming product/derivative evaluation through the rings entry point.
- Route smooth/coprime factor-base extraction through cached or constructed product
  trees; preserve sorted/combined factors, overlap quirks and signed floor quotients.
- Restore zero-base errors and report source-proven nontermination explicitly instead
  of returning misleading results; document the timed-original regression adapter.
- Add 12,956 original comparisons, including tree boundaries through 129 leaves;
  a 2,332-case safe old-input subset exposes 300 mismatches.

## 8.15.0 - 2026-09-10

- Add CRT list and modular-element overloads, preserving singleton parent identity and
  original balanced combination order, errors and input arrays.
- Restore signed CRT basis residues, empty-option output and bundled partial-prefix
  behavior in the non-coprime fallback; preserve ragged-vector error order.
- Add 11,278 bundled-source comparisons covering both prime-field backends, native/wrapped
  inputs and malformed shapes; the previous implementation differs in 3,610 cases.

## 8.14.0 - 2026-09-10

- Implement the mirrored native-width gcd/inverse arithmetic classes and restore
  factory thresholds, bound-method selection and IntegerLike callback arguments.
- Preserve native conversion overflow, signed-modulus shortcuts and Cython error
  sentinels, documenting unstable diagnostics and native undefined arithmetic limits.
- Add 13,704 original comparisons; the initial 13,320 cases expose 2,902 old mismatches.

## 8.13.0 - 2026-09-10

- Preserve negative continuant orders, empty-input errors, truncated large orders and
  bigint returns for IntegerLike entries.
- Restore zero quadratic-residue errors and zero fundamental discriminants; normalize
  IntegerLike scalar inputs and remove wrapped-zero loops in odd/coprime-part helpers.
- Replace the 32-bit squarefree-divisor mask with Sage's progressive arbitrary-width
  subset iteration, preserving lazy prefixes and original zero-factorization errors.
- Add 8,468 original comparisons including large signed inputs, mapped return types and
  prefixes for integers with up to 64 prime factors; safe old cases expose 2,046 mismatches.

## 8.12.0 - 2026-09-10

- Delegate Hilbert symbols to PARI by default, preserve direct/all algorithm dispatch,
  reject unknown algorithms and normalize rational coefficients by their square classes.
- Fix wrapped-zero symbol loops and wrapped-unit conductor-inverse loops by normalizing
  IntegerLike inputs before control flow and arithmetic.
- Port native integer Hilbert valuation-parity formulas and divide-and-conquer unit
  extraction with original real-place,
  signed-unit, zero and modulus-error behavior.
- Add 22,724 original comparisons, including 5,274 native/large-valuation/wrapped-inverse
  cases; a safe old-code subset exposes 4,231 prior discrepancies.

## 8.11.0 - 2026-09-10

- Delegate prime-power extraction to the existing PARI port, removing duplicate root
  code and fixing false negatives for wrapped large pseudoprime powers.
- Port PARI's inclusive successor/predecessor residue wheel and delegate Integer prime
  traversal; preserve the distinct free-function and Integer previous-prime errors.
- Restore Sage's odd-candidate prime-power traversal and educational Eratosthenes sieve;
  accept IntegerLike consistently across the affected range/predicate functions.
- Add 9,870 original comparisons for wrapped values, powers, ranges, wheel residues,
  and 64-bit boundaries, exposing 32 prior discrepancies.

## 8.10.0 - 2026-09-10

- Delegate floating algebraic dependencies to PARI and restore exact integer/rational
  shortcuts, irreducible factors, height-bound nulls, proof errors and hint precedence.
- Preserve real-field input precision through relation powers and factor selection;
  real-field method aliases now forward the original value instead of a double.
- Add nearest-even MPFR addition/multiplication with native state, operand aliasing,
  singular signs and exact native midpoint status; expose and document both primitives.
- Add 33,642 permanent original comparisons, including direct PARI and MPFR calls,
  precision/proof boundaries, nonfinite values and 64-bit conversion limits.
- Record remaining native reduction/proof-interval and complex/p-adic audit gaps.

## 8.9.0 - 2026-09-10

- Dispatch integer floor/ceiling through exact numeric methods before Python-style float
  conversion, accepting the existing FloatInput union and preserving original errors.
- Make truncation choose floor/ceiling after the original sign comparison, including
  real-field NaN/infinity conventions and rejection of text/bytes/null comparisons.
- Preserve float fallback after method AttributeError and propagate other method errors.
- Add 1,569 original comparisons exposing 1,242 prior discrepancies, with exact large values,
  rational/real/modular inputs, numeric strings, subnormal bits and method-error protocols.

## 8.8.0 - 2026-09-10

- Delegate Dedekind sums to FLINT by default and PARI when selected, rejecting unknown
  algorithms before input coercion and accepting IntegerLike arguments.
- Port the FLINT continued-fraction matrix and PARI Knuth recurrences, preserving signed
  moduli, noncoprime reduction, word/generic boundary conventions and zero-division errors.
- Add 8,298 original comparisons and 59 fast original-vector/boundary tests; replace stale
  deviation claims with the remaining native quotient-batching performance gap.

## 8.7.0 - 2026-09-10

- Delegate nth-prime lookup to a PARI port with the original 64-bit ordinal table,
  segmented sieve/next-prime successor search and exact interval sizing.
- Accept IntegerLike nth-prime inputs and preserve native index/overflow errors;
  use the original nth-prime bound followed by prime_range for the first n primes.
- Truncate positive fractional first-prime counts after the original boundary checks,
  retaining the numeric-count API and matching nonfinite validation.
- Add 795 original comparisons, including every PARI table anchor and two successors
  through ordinal 100 billion; document cache/strategy adaptations and execute API examples.

## 8.6.0 - 2026-09-10

- Restore Sage's single-argument/default prime iterator bounds and lazy finite iteration;
  accept integer wrappers and the explicit Infinity stop, forwarding proof selection.
- Fix successive-difference validation and exhaustion, and preserve the original
  subfactorial negative-input error; both accept IntegerLike inputs.
- Replace the duplicated subfactorial oracle with Sage's actual function and add
  1,228 original comparisons for sequence boundaries, wrappers and enormous lazy intervals.
- Correct stale DESIGN guidance that had incorrectly claimed finite JavaScript numbers
  were accepted by IntegerLike coercion; execute the updated public API examples.

## 8.5.0 - 2026-09-10

- Forward Gaussian copy-helper overrides through constructor coercion/validation,
  preserving large/fractional centers, precision options and explicit null handling.
- Allow numeric centers and tail cutoffs in the public Gaussian options type, matching
  the constructor and copy helper; nonintegral numeric tails retain original validation.
- Add 432 original constructor/stream comparisons that exposed 225 copy discrepancies.

## 8.4.0 - 2026-09-10

- Preserve native fractional parts, signed-zero/infinity conventions and real comparisons
  after the original common-parent coercion, including literal reparsing and Python floats.
- Add MPFR value copying, comparison and fractional extraction with native inexact status,
  including the 64-bit limb paths for halfway cases.
- Add 31,500 original comparisons for fractional values, cross-precision/literal/float
  comparisons, native state transitions, aliasing and halfway rounding.

## 8.3.0 - 2026-09-10

- Evaluate real-number predicates, signs, multiplicative orders and integer rounding
  from native MPFR state, preserving large/small finite values and original singular errors.
- Add MPFR classification, signed-integer comparison, integer rounding and extraction,
  including directed modes, in-place operations and native inexact status.
- Add 23,484 original comparisons: 3,600 Sage observations and 19,884 direct MPFR cases,
  including extreme-exponent regressions that verify precision-bounded rounding.

## 8.2.0 - 2026-09-10

- Convert nearest-even real-field numeric inputs at the requested precision through
  MPFR double/integer setters and parse ordinary strings with the native conversion path.
- Add `mpfr_set_z` with integer normalization, nearest-even rounding and native inexact status.
- Add 6,267 original comparisons: 1,995 real constructor/negation/absolute-value cases
  exposing 1,661 old discrepancies and 4,272 direct native integer setter cases.

## 8.1.1 - 2026-09-10

- Restore integer Gaussian constructor validation order, missing/null argument errors,
  and rejection of fractional tail cutoffs.
- Round large integer centers through the 53-bit real field and compute exposed bounds
  with exact integer addition/subtraction.
- Add 806 original comparisons for combined invalid arguments, all algorithm selections,
  large-center rounding, bounds and seeded sample sequences.

## 8.1.0 - 2026-09-10

- Add the native MPFR double setter with precision rounding, signed zeros, subnormal
  normalization, inexact status and destination-state preservation across NaNs.
- Format Gaussian parameters through MPFR digits and Sage's Decimal rounding, preserving
  halfway cases, negative-zero centers, large fixed values and scientific notation.
- Add 8,265 original comparisons, including 6,720 native MPFR conversion/state cases,
  160 direct numeric formatting cases and 1,385 sampler representations that exposed 810 old discrepancies.

## 8.0.0 - 2026-09-10

- Default lattice Gaussian sigma to 1 and preserve original sigma/basis/center validation
  order, empty/dependent bases, row-space checks and failures deferred until sampling.
- Preserve `c=None` and the original center-change cache state. The center accessors now
  return null when deferred; the low-level `_call()` returns zero for an empty basis.
- Format nonpositive integer Gaussian sigma errors with Python's six-place, nearest-even
  rounding, including negative zero, negative infinity and large finite values.
- Compare exception messages as well as classes in the shared property harness. Adapt
  two older Sage lattice errors to the bundled source's spelling.
- Add 2,913 comparative regressions for constructor inputs, state changes and seeded streams.

## 7.2.0 - 2026-09-10

- Convert binary randomization density through Python float rules, including exact numeric
  wrappers, Unicode decimal text, byte strings, overflow, signed NaN and conversion errors.
- Preserve empty-matrix shortcuts and factory default handling before density conversion.
- Match the original empty PNG image error and compare binary image pixels and decoding.
- Add 10,961 original comparisons using CPython 3.12/Unicode 15 and bundled Gamma predicates;
  track comparative dispatch for all 51 implemented binary matrix API callables.

## 7.1.0 - 2026-09-10

- Use Sage's GMP-backed random state for binary matrices, preserving two-word draws,
  native bit order, padding consumption and sparse updates with replacement.
- Preserve randomization cache/error behavior for NaN, empty matrices and immutable inputs.
- Restore explicit-density random factories to select nonzero entries; preserve uniform
  random bits when density is omitted. Factory dimensions accept IntegerLike inputs.
- Construct identities after dimension conversion, and match inferred-row shape errors
  before coefficient conversion. Entry factories accept the constructor's coefficient types.
- Add 10,032 original comparisons, including random-stream advancement after repeated calls.

## 7.0.0 - 2026-09-10

- Restore binary cached row identity and immutability for `row(index, true)`. Assigning
  to these rows now raises Sage's immutable-vector error; default row copies remain mutable.
- Validate binary constructor dimensions, flat/nested sequence lengths and scalar shapes
  in the original order; accept IntegerLike dimensions and scalar/flat entry inputs.
- Coerce constructor and assigned entries through the prime-field backend, preserving
  large integers, rational reduction and original invalid-conversion errors.
- Add 11,554 comparative regressions covering constructor/assignment coercion, sequence
  validation, dimension limits and row cache invalidation.

## 6.1.0 - 2026-09-10

- Restore binary string mapping precedence and dictionary mutation, callable mappings,
  global alignment, Unicode code-point widths and original invalid-mapping errors.
- Add subdivision setters/getters, preserving sorted duplicate lines, coercion errors,
  mutability checks and cache invalidation in the original order.
- Preserve subdivision metadata through copies, negation, transpose and augmented
  boundaries; retain zero-size shortcuts and fresh metadata for arithmetic/slices.
- Add 23,178 original comparisons and a shared matrix0 string renderer. The installed
  Sage transpose oracle applies the metadata-preservation step present in bundled source.
  Boundary probes compile the bundled string method to avoid older C integer limits.
- Stream ASCII-safe property transcripts in chunks so large Unicode-bearing suites
  do not exceed the JavaScript engine's single-string allocation limit.

## 6.0.0 - 2026-09-10

- Restore binary column caching: column vectors are immutable, default copies share
  vectors, and `columns(false)` returns the mutable cached outer list. This corrects
  the previous independently mutable vector arrays and changes mutation behavior.
- Delegate binary addition, transpose and augmentation to M4RI; transpose uses recursive
  blocks and word-parallel bit swaps instead of a per-entry matrix loop.
- Preserve original addition/subtraction parent errors, determinant errors and
  augmentation dimension errors, with 10,773 new original comparisons.
- Expose native transpose/concatenation adapters, retain tested native empty-transpose
  faults, and document source buffer/caching representations with an executed example.

## 5.3.0 - 2026-09-10

- Restore rectangular and singular binary linear solves, default consistency checking,
  square-system dispatch, dimension errors and original cache effects.
- Delegate binary solves and computed kernel bases to M4RI PLUQ and triangular routines;
  preserve native mutated outputs, padding checks and tested empty-buffer faults.
- Restore kernel algorithm/basis/proof option validation, generic/pivot dispatch and
  immutable identity kernels for matrices with no rows.
- Add 12,222 original comparisons, compiling the bundled solve method absent in the
  installed Sage version; document the native adapters and expanded kernel options.

## 5.2.0 - 2026-09-10

- Restore binary swap index conversion before mutability/cache checks and preserve
  signed index boundaries, rational conversion and original errors.
- Apply binary permutations in original cycle order, including smaller degrees,
  fixed points outside the matrix, partial mutations and invalid image lists.
- Delegate binary inversion to M4RI's aligned augmented-window algorithm and repair
  nonsquare inverse errors; expose the native inverse kernel with its original k behavior.
- Add 18,153 permanent original comparisons and an executed public API example.

## 5.1.0 - 2026-09-10

- Delegate standalone binary PLE/PLUQ to their original standard, Russian and naive
  M4RI algorithms; preserve copied factors and native permutation lists.
- Restore C-int parameter conversion and validation order, including fractional,
  nonfinite, boolean and IntegerLike inputs, with 3,588 new original comparisons.
- Correct the PLE backend to the pinned seven-table schedule; preserve tested larger
  table parameters and translate native zero-column faults into exceptions.
- Remove direct private-storage access from the standalone wrappers, resolving 16
  pre-existing type diagnostics. Update API examples and documented backend bounds.

## 5.0.0 - 2026-09-10

- Restore binary elimination algorithm dispatch, cached forms and validation order;
  computed forms are immutable and empty in-place elimination returns self.
- Delegate rank and elimination to M4RI panel/table and recursive PLE/PLUQ algorithms,
  including native triangular solves and density-based crossover.
- Preserve classical reduction, keyword parameter errors, mutation cache invalidation,
  swap errors and original lexical-ordering cache and immutability behavior.
- Add 12,924 permanent original comparisons, with pinned native kernels and explicit
  adapters for differences in the installed oracle version.

## 4.0.0 - 2026-09-10

- Restore binary matrix density's exact Rational/integer and approximate RealLiteral return
  types, native M4RI sampling, word-boundary behavior and safe native fault translation.
- Add a bounded MPFR decimal conversion dependency and original RealLiteral factory,
  precision changes, formatting, exact observations and bundled float conversion behavior.
- Fix invalid real precision error classes and mutable native NaN parsing; add permanent
  original comparisons for each discrepancy, including binary rounding boundaries.
- Add 24,161 permanent original comparisons; all 47,751 live affected-area checks pass.
- Document conversion limits and the existing real arithmetic limitations; execute the
  public API examples and keep all workspace versions synchronized.

## 3.7.0 - 2026-09-10

- Restore binary matrix negative indexing, remaining-extent submatrix sizes and original
  bounds/conversion errors, including validation order for empty matrices.
- Accept IntegerLike, Rational and boolean indices alongside existing numeric indices;
  preserve Sage's truncating rational reads and exact rational row/slice/write indices.
- Delegate row/slice extraction to a packed M4RI submatrix port without scanning unrelated
  source rows or columns; document and execute the expanded public signatures.
- Add 19,522 original comparisons covering scalar index conversions, signed C boundaries,
  aligned/unaligned word slices, zero dimensions, padding bits and allocation identity.

## 3.6.0 - 2026-09-10

- Support specialized GF(2) matrices in polynomial evaluation and preserve matrix results
  in direct, mixed and argument-list TypeScript overloads.
- Replace binary matrix scalar cubic multiplication with a packed M4RI dependency port:
  native classical dispatch, Gray-code tables, M4RM and Strassen-Winograd/Bodrato squares.
- Add original binary matrix multiplication methods, parameter/dimension validation and
  safe translations of native cutoff faults; compare against the checksum-pinned source.
- Add 5,823 original comparisons, including large matrices, native algorithm cutoffs,
  allocation identity and compiler-checked evaluation return types.
- Verify all 374,423 configured comparisons; update the audit inventory with fresh
  execution coverage and preserve explicit open implementation gaps.

## 3.5.0 - 2026-09-10

- Restore polynomial evaluation at generic, integer and modular matrices, with canonical
  coefficient parents, scalar actions, sparse compiled powers and original dimension errors.
- Port FLINT's truncating integer-to-double conversion and Sage's integer matrix conversion;
  preserve rounding boundaries instead of applying scalar integer conversion to matrices.
- Add 21,017 original comparisons, including real result bits, nested coefficients, repeated
  calls, rectangular matrices, native conversions and degree cutoffs. Document oracle adapters.

## 3.4.0 - 2026-09-10

- Restore polynomial evaluation's canonical parents, argument unpacking, native dispatch
  and sparse compiled operation order; preserve IEEE outputs and validation errors.
- Port FLINT integer/rational evaluation, composition and Taylor shifts, modular evaluation
  and composition, and NTL prime/extension scalar evaluation with original comparisons.
- Cache polynomial generators and preserve native generator tests, zero-ring representation
  and large-prime first-power allocation; route composition through original evaluation.
- Add 34,922 comparative regressions for every discrepancy and document native raw-denominator
  oracle compilation, new public APIs and remaining evaluation-action gaps.

## 3.3.0 - 2026-09-10

- Restore polynomial modulus handling, canonical coercion, signed-boundary dispatch and
  native error types; add FLINT modular powers/inverse series and NTL extension powers/XGCD.
- Add the GSL numeric subset for RDF integer powers and inversion. Preserve real division
  operation order, bounded NaN remainders and unchecked quotient representation.
- Bound the modular half-GCD/Newton product helper; 16 large product, division and GCD
  regressions previously raised BigInt out-of-memory errors.
- Add 16,736 original comparisons for every discovered discrepancy and new kernel;
  document exported APIs, oracle version adapters and remaining portable backend limits.

## 3.2.0 - 2026-09-10

- Delegate full exact polynomial multiplication to FLINT/NTL and implement SageMath's
  generic Karatsuba recursion, replacing quadratic public-backend products.
- Expose native integer/rational and NTL modular/extension full-product array adapters;
  preserve rational cross-content cancellation and aliased square representation.
- Add 1,544 original comparisons for native products, all exact coefficient backends,
  large products, unequal degrees and noncommutative matrix coefficients.

## 3.1.2 - 2026-09-10

- Correct truncated multiplication's canonical scalar parents, nested coefficient
  promotion, nullable internal operands and original validation errors.
- Preserve real-polynomial product operation order, signed zeros and optimized squares;
  print scientific-notation and infinite coefficients as SageMath does.
- Add 14,012 comparative regressions for parent/type matrices and real coefficient bits.

## 3.1.1 - 2026-09-10

- Share bounded FLINT products with integer divide-and-conquer division and modular
  KS/KS2/KS4 multiplication to avoid oversized BigInt temporaries.
- Add 38 native comparative regressions for large quotients, products and squares;
  28 of those cases failed before the repair.

## 3.1.0 - 2026-09-10

- Restore rational polynomial exponents, exact roots and Newton truncated-series methods,
  including backend-specific typed argument validation and identity rules.
- Add FLINT truncated products, powers and inverses, plus NTL extension inverse series
  and its half-GCD coefficient kernel. Replace p-adic trial-division prime validation.
- Fix polynomial factor ordering by numeric coefficients and multiplicity, including
  deterministic finite-ring root selection and original factorization errors.
- Add 39,926 permanent original comparisons for every discovered discrepancy and new
  kernel; document the remaining backend, parent-category and audit limits.

## 3.0.0 - 2026-09-10

- Correct polynomial integer powers, negative fraction results, native overflow/zero
  behavior, identities and coefficient-ring preservation; widen the dynamic/negative
  return type to a polynomial/fraction union (TypeScript breaking change).
- Add FLINT integer/rational and NTL binary/prime/extension polynomial power adapters.
  Bound packed BigInt products with Karatsuba to prevent temporary-size failures.
- Add 7,856 original comparisons, including all-coefficient native checks, large-product
  regressions and strict rational coefficient adapters; document portable backend limits.

## 2.19.0 - 2026-09-10

- Restore generic fraction and polynomial scalar hooks, native FpT scalar sections,
  exact modular fallback conversions and GF2 scalar constructor support.
- Preserve native zero allocation during fraction reduction, constant-polynomial section
  precedence, ring unit-inverse errors and the bundled source's operand-prefixed errors.
- Add 7,338 original comparisons for scalar conversions, direct hooks, mutation and
  source identities; retain a narrow installed-Sage message adapter.

## 2.18.0 - 2026-09-10

- Restore polynomial/fraction partial conversions, native normalization side effects,
  mixed-parent fraction construction and standalone fraction-module imports.
- Expose polynomial parent coercion queries; repair nested coefficient inverses,
  constant-polynomial sections and polynomial/fraction equality.
- Add 16,937 comparative regressions checking values, errors, canonical parents,
  accessor identity and mutation; document the noninterned-parent oracle mapping.

## 2.17.0 - 2026-09-10

- Restore generic and specialized polynomial fraction classes, nonconstant denominators,
  reduction flags, inversion, powers, native equality, accessor identity and cached parents.
- Delegate small-prime fraction arithmetic to FLINT classical/KS/KS2/KS4 products,
  modular addition/subtraction, powers, GCD and division.
- Add 16,017 comparative fraction protocols and direct native-kernel regressions, including
  coefficient carries, aliasing, full-word moduli and safe zero-inversion oracle handling.

## 2.16.0 - 2026-09-10

- Restore pseudo-division scalar protocols, coefficient membership, negative scaling
  powers, mixed-parent arithmetic, error ordering and fraction-field zero quotients.
- Add the constant-denominator polynomial fraction-field subset and repair zero-ring
  scalar equality without introducing a synthetic polynomial coercion.
- Add 15,994 comparative regressions covering pseudo-division and direct fraction
  arithmetic/data protocols; document remaining nonconstant fraction and power backends.

## 2.15.0 - 2026-09-10

- Restore polynomial derivative variables, repetition counts, aliases, gradients,
  nested coefficient differentiation and native zero-result identity.
- Delegate QQ derivatives to FLINT ports and use native scalar coefficient multiplication,
  correcting RDF rounding and removing repeated-addition overhead.
- Add 21,025 comparative regressions for the public protocol, parser, native kernels
  and isolated real-coefficient rounding; expose the shared derivative utility module.

## 2.14.0 - 2026-09-10

- Port PARI's single-word real addition/division and general real Sylvester elimination,
  including native exponent-based pivots, approximate-zero checks and cancellation precision.
- Preserve inexact leading zeros at the PARI array boundary and Sage's earlier normalization.
- Add 7,062 comparative cases for native primitives, composed rounding, general resultants
  and explicit zero controls, with bit-level native and Sage outputs.

## 2.13.0 - 2026-09-10

- Port FLINT rational resultants and delegate Sage ZZ, QQ and word-prime resultants
  to native ports; restore composite matrix, scalar-parent and proof-keyword behavior.
- Add the RDF numeric coefficient subset and PARI constant-operand real resultants,
  retaining native binary64 conversion boundaries, rounding and exceptional values.
- Add 23,968 comparative regressions across native kernels, parent combinations,
  scalar inputs, aliasing and real coefficient bit patterns; document remaining backends.

## 2.12.0 - 2026-09-10

- Repair Sylvester-matrix canonical parents, scalar operands, zero-input errors,
  optional variable coercion and native small PARI-field matrix failures.
- Restore explicit finite-element polynomial representatives, specialized GF2 coefficient
  conversion and the distinct QQ singleton-list constructor path.
- Add 18,353 original-Sage comparisons, direct constructor regressions and executable
  API examples; retain existing resultant/discriminant behavior during separate backend work.

## 2.11.0 - 2026-09-10

- Port FLINT integer XGCD with resultant-scaled modular reconstruction and rational
  XGCD with native content/GCD removal and denominator normalization.
- Implement Sage integer XGCD, preserving Integer constant/zero returns, native
  input-order signs and the denominator-cleared QQ fallback for shared factors.
- Delegate QQ XGCD to its FLINT port. Add 1,791 comparative regressions covering
  coefficients, parents, aliasing, scalar types, reconstruction and noncoprime inputs.

## 2.10.0 - 2026-09-10

- Port FLINT integer resultants with signed content, subresultant/Euclidean dispatch,
  norm bounds, skipped word primes and balanced modular reconstruction.
- Port modular resultants with small formulas and recursive half-GCD accumulation,
  including deferred degree-drop factors and native pointer-alias behavior.
- Add 1,597 direct native FLINT comparisons across zeros, signs, composite errors,
  dispatch boundaries, large coefficients, skipped primes and common factors.

## 2.9.0 - 2026-09-10

- Port FLINT modular extended GCD with its Euclidean/recursive half-GCD cutoff tables,
  sharing fast arithmetic with GCD. Delegate Sage word-modular and binary XGCD to ports.
- Repair extended-GCD canonical parents, native zero conventions, composite errors and
  output aliasing. Add mixed coefficient overloads and nonconstant extension comparisons.
- Add 6,556 Sage comparisons and 1,328 direct native FLINT comparisons covering
  coefficients, zero/equal inputs, native errors, cutoffs and unbalanced shapes.

## 2.8.0 - 2026-09-10

- Add a FLINT rational GCD numerator kernel, delegating primitive integer GCD to the
  existing FLINT port. Delegate Sage word-modular and binary GCD to dependency ports.
- Repair mixed-parent GCD coercion, native finite zero/equal shortcuts, allocation
  identity and composite-modulus errors; stop modular Euclid at a nonzero constant.
- Add 5,654 Sage GCD comparisons and 301 rational-backend comparisons, with a narrow
  signal-guard adapter for composite GCD errors in the installed older Sage version.

## 2.7.1 - 2026-09-10

- Repair polynomial equality with integer, rational and finite scalars, including
  binary64 rounding, underflow and infinity in comparisons over rational bases.
- Embed coefficient polynomials into outer rings for mixed arithmetic and preserve
  generic zero-polynomial identity during scalar multiplication.
- Add 9,790 permanent comparative regressions covering coefficients, parents and identity.

## 2.7.0 - 2026-09-10

- Add dense FLINT integer and modular quotient/remainder kernels, sharing exact
  integer division with GCD and preserving recursive/Newton algorithm dispatch.
- Repair native ZZ partial-quotient behavior, polynomial division parent coercion,
  backend errors/identity and exact nested coefficient division; delegate binary
  polynomial division to NTL. Export NTLError and mixed-coefficient division overloads.
- Add 7,002 permanent comparative cases, including direct native FLINT checks;
  correct the prior unit assertion that required native ZZ division to be exact.

## 2.6.2 - 2026-09-09

- Select a common coefficient parent before polynomial addition, subtraction,
  multiplication and equality, including quotient-ring reductions and nested bases.
- Restore Sage polynomial signs, unit coefficients and parentheses in string output;
  correct unit expectations that had preserved the port's old formatting bugs.
- Add 7,554 shared Sage/TypeScript comparisons, mixed-type arithmetic overloads,
  colocated regressions and executable API examples.

## 2.6.1 - 2026-09-09

- Repair polynomial coefficient indexing, shifts, truncation and reversal across the
  implemented QQ, ZZ, finite and nested coefficient parents, including native integer
  bounds, scalar-parent arithmetic, fractional conversion, errors and result identity.
- Add 9,608 permanent Sage/TypeScript comparisons and colocated regression tests;
  document the widened index signatures and narrowly scoped Sage 10.3 oracle adapters.

## 2.6.0 - 2026-09-09

- Make nonmonic integer polynomials normalize over QQ without floor division, and retain
  identity for already-monic integer polynomials. Compare monicity with the coefficient-ring one.
- Preserve backend-specific zero/nonunit errors and modular allocation behavior; delegate
  word-modular normalization to the new FLINT `_nmod_poly_make_monic` array kernel.
- Add 573 comparative cases for coefficients, result parents, identity, zero and nonunit
  errors, extension-valued leading coefficients and specialized GF2 errors. All 2,244 area
  cases and 573 saved-Sage replays pass. Correct the old GF(2) zero-normalization unit assertion.
- Align specialized GF2 division/inversion errors with Sage and delegate powers to the
  existing IntegerMod implementation, including fractional exponents and the GMP cutoff.
- Document Integer/Rational result-type overloads, the new FLINT export and the remaining
  Sage-version and generic-base-extension boundaries. Final full saved-Sage comparison:
  54,674 pass; fast: 5,931 pass/32 skips; slow: 946 pass; type diagnostics unchanged at 609.
  Current measured source-line coverage is 76.58%.

## 2.5.0 - 2026-09-09

- Replace Sage-local integer-polynomial GCD with FLINT dependency delegation, including
  squarefree factorization's primitive GCD. Port subresultant, packed heuristic and modular
  kernels with coefficient/content normalization, CRT bounds and unlucky-prime handling.
- Add a modular polynomial backend with fast packed products, Newton division and recursive
  half-GCD using FLINT's cutoff tables. Integer divisibility uses recursive exact division.
- Add 896 comparisons against FLINT C/Sage, including heuristic failure, degree-3,500 and
  4,300-bit inputs, prime-table fallback, skipped primes and modular degree changes.
- Document the exported dense-array APIs and remaining native dispatch/performance boundaries.
  Existing mutable FLINT APIs remain stubbed. Full comparison: 54,060 pass; final area:
  1,142 pass; new replay: 896 pass; fast tests: 5,926 pass/32 skips; slow tests: 946 pass.
  Type diagnostics remain 609; current measured line coverage is 76.57%.

## 2.4.1 - 2026-09-09

- Replace test-local integer-polynomial GCD and evaluation algorithms with calls to the
  production Polynomial methods. Stop stripping content from Sage's expected GCD results.
- Preserve Sage's integer-polynomial zero/one GCD shortcuts, including negative signs and
  operand identity. Correct an old unit assertion that mislabeled positive normalization
  on the zero branch as Sage behavior.
- Add 216 comparative cases for zero, content, signs, trailing zeros and large exact values.
  All 53,205 full-suite cases pass, including 246 polynomial_ops cases.
  Fast tests: 5,921 pass with 32 existing skips;
  focused polynomial tests: 46 pass; saved-Sage replay: 216 pass. Type-check diagnostics
  remain unchanged at 609. Current measured lines: 80,781/106,358 (75.95%).
- Keep the separate integer-polynomial GCD backend dispatch gap open for the next repair.

## 2.4.0 - 2026-09-09

- Coerce polynomial coefficient arrays before normalization, preserve exact-parent identity
  and polynomial-base constant embedding, and validate generator indices. Follow QQ's
  singleton-list and GF2X's recursive-list constructor paths and broaden input signatures.
- Add rational-lift hooks to modular and prime elements. Accept compatible modular inputs
  in extension parents and preserve extension-to-prime conversion errors. Repair the
  Hilbert-class-polynomial Integer coefficient-ring adapter exposed by constructor validation.
- Add 1,626 comparative regressions. All 52,989 full-suite cases pass, including 1,671 in
  the polynomial area; all new cases also pass a saved-Sage transcript replay. Fast tests:
  5,920 pass with 32 existing skips before the final API example; focused source/replay:
  1,920 pass; final API/constructor tests: 41 pass. All 946 affected slow tests pass.
  Current measured source-line coverage: 80,543/106,155 (75.87%).
- Existing type-check diagnostics decrease from 763 to 609 after aligning polynomial input
  interfaces. The broader source audit and polynomial coercion/algorithm coverage remain open.

## 2.3.0 - 2026-09-09

- Delegate extension arithmetic to PARI-port FpX/FpXQ kernels, removing local multiplication,
  inversion and power implementations. Preserve exponent comparison/error order before ZZ
  coercion, including zero to negative fractions and signed nonfinite values.
- Share signed quotient powers across FpXQ_pow and FpXQ_powBig, using PARI's binary/sliding-window
  thresholds. Support three-argument FpXQ_inv while retaining its four-argument Hensel API.
- Add 3,940 comparative cases, including direct installed PARI C-kernel calls and exponent
  boundaries through 1,024 bits. All 51,363 full-suite cases pass; fast tests: 5,914 pass,
  32 existing skips before the final API example; affected Sage slow tests: 946 pass.
  Direct Galois tests: 62 pass; focused source tests: 219 pass; API examples: 35 pass;
  all 3,940 new cases also pass a final transcript replay.
- Keep native polynomial packing/multiplication optimizations and low-level PARI diagnostic
  payloads explicit as remaining boundaries. Type-check diagnostics unchanged at 763.
  Current instrumented lines: 81,254/106,262 (76.47%).

## 2.2.0 - 2026-09-09

- Validate finite-field generator indices with Python truthiness and quotient/QQ indices by
  equality to zero. Preserve Sage's distinct errors and broaden gen input signatures.
- Coerce both prime parent orders through ZZ before validation. Preserve the legacy
  check:false positive-order rule, including order one. Add modern prime/extension _integer_
  hooks and reject nonconstant extension-to-integer conversions with Sage's exact error.
- Add 1,209 comparative cases. Full suite: 47,380 pass before the final 43 QQ-index cases;
  final affected areas: 31,066 pass. Fast units: 5,909 pass, 32 existing skips before final
  QQ/API additions; focused units/docs: 269 pass; all 946 affected slow tests pass with a
  15-second per-test timeout. An additional 1,613-case transcript replay passes.
- Type-check diagnostics remain unchanged at 763. Updated API documentation and coverage:
  101/101 finite-field members dispatched; instrumented lines 81,193/106,081 (76.54%).

## 2.1.0 - 2026-09-09

- Share canonical finite-ring coercion across modular elements, both prime-field classes and
  extensions. Preserve compatible result parents and operand order; reject constructor-only
  scalar conversions in arithmetic and return false for incompatible equality.
- Support Python string/list repetition in native modular and prime multiplication, including
  Unicode, shallow list identity and signed index bounds. Correct an earlier comparative
  adapter that called .parent() on valid Python sequences and mislabeled them as errors.
- Add 7,979 comparative regressions and update arithmetic return overloads and API examples.
  All 46,214 full-suite comparisons pass; fast units: 5,907 pass, 32 existing skips; focused
  field/API tests: 251 pass. Existing type-check failures remain (763 diagnostic instances,
  down from 1,397, with some changed overload/inference diagnostics). Affected slow suite:
  945 pass and one timeout, which passes in isolation in 2.76 s. Final focused/transcript
  replay: 252 pass. Instrumented line coverage: 81,148/105,804 (76.70%).

## 2.0.1 - 2026-09-09

- Share IntegerMod power dispatch across both prime-field implementations, preserving native
  versus GMP error paths, original exponent types and exact scalar conversion diagnostics.
- Broaden pow input signatures to the tested Integer/Rational/number/boolean/string/null
  inputs; retain results in the original parent. Update API examples and architecture notes.
- Add 897 comparative regressions. All 38,235 full-suite comparisons and 946 affected slow
  tests pass. Fast suite: 5,902 pass, 32 existing skips before the final colocated regression;
  final focused units/docs: 119 pass. Type-check diagnostic instances are unchanged.

## 2.0.0 - 2026-09-09

- Breaking correction: Zmod/Integers return ZZ for zero/default orders; Mod with zero returns
  its original input unchanged. Conditional return types expose possibly-zero dynamic orders.
- Normalize negative factory orders, preserve signed cache keys with weak references, and
  honor explicit modular or prime-field parents in Mod. Keep direct ring constructors finite.
- Update positive-modulus consumers, API examples and stale factory assertions. Preserve LWE's
  existing nonzero-modulus requirement and document remaining category/coercion limits.
- Add 373 comparative regressions. Full suite: 37,338 pass; modular area: 3,890 pass. Fast units:
  5,902 pass, 32 existing skips; affected slow units: 946 pass; focused units: 454 pass.
  Type-check diagnostics are unchanged from 1.2.0. This version does not mark audit completion.

## 1.2.0 - 2026-09-09

- Repair modular-integer scalar conversion, arithmetic coercion, quotient/prime-field result
  parents, equality and native/GMP negative-power errors. Give Mod elements complete parents.
- Use Sage's Python random stream for modular rings and support the symmetric bound argument.
- Delegate modular multiplicative order to the new PARI znorder port, preserving its
  prime-power valuation and lcm algorithm. Update public signatures and examples.
- Add 3,517 comparative cases and focused regressions; dispatch all 34 public modular members.
  Full suite: 36,929 pass before 36 final string controls; final modular area passes.
  Fast units: 5,895 pass, 32 existing skips; affected slow units: 946 pass; focused: 112 pass.
- Document remaining factory gaps. Existing type-check diagnostic kinds persist: overloads
  repeat 18 prior RingElement diagnostics and broader inputs remove 149 diagnostic instances.

## 1.1.0 - 2026-09-09

- Parse finite-field polynomial strings through a port of Sage's tokenizer/arithmetic grammar;
  preserve implicit multiplication, precedence, cancellation, coercion and syntax errors.
- Reuse the existing rational-function normalization for intermediate fractions; require a
  polynomial before reducing an extension-field element. Preserve backend-specific diagnostics
  and polynomial-ring display labels, including GF2X and NTL.
- Add 867 field-string comparisons and a 566-case direct parser/tokenizer area. Cover Unicode
  Digit characters, Python whitespace, token state, C-int/character bounds, callback dictionaries
  and each implemented grammar entry point. Function/sequence/matrix parser stubs remain explicit.
- Validation: 33,443 full-suite comparisons pass; final parser area 566 pass, including five
  later boundary cases. Fast tests: 5,895 pass, 32 existing skips; 946 affected slow tests pass.
  Existing type-check diagnostics are unchanged after normalizing the widened string-input union.

## 1.0.0 - 2026-09-09

- Breaking correction: prime-field `sqrt()` now follows Sage's `extend=True` default,
  returning an extension element for nonsquares. Use `sqrt({extend:false})` to retain
  the original field and narrow the return type. This major increment follows the
  repository's breaking-change rule; it does not mark the ongoing audit complete.
- Implement all-root options, canonical root ordering, extension naming/moduli and exact
  keyword errors in both prime-element classes. Share the PARI-backed implementation.
- Preserve Sage's small-modulus empty-root shortcut and its general extension/all
  NotImplementedError. Document the newer reference's error-message oracle adapter.
- Keep elliptic-curve point enumeration in its original field explicitly. Add 1,572
  comparative regressions across small, machine-word and 127-bit fields and update API docs.
- Validation: 32,003 full-suite comparisons pass before twelve final keyword cases; final
  finite-field area: 14,037 pass. Fast units: 5,888 pass and 32 existing skips; type-check
  diagnostics unchanged. All 946 affected slow tests and 232 focused field/docs tests pass.

## 0.0.26 - 2026-09-09

- Make GF/FiniteField honor name and modulus options, full coefficient lists and polynomial
  inputs; normalize moduli, preserve custom linear generators, and validate names and degrees.
- Match Python's generator-list whitespace and exact name diagnostics; add the missing
  degree-one Conway polynomials for all existing characteristics.
- Replace the factory's exhaustive exponent scan with the existing PARI Z_isanypower port
  and export that backend. Compare factory boundaries, Conway entries and order helpers.
- Add 1,328 shared comparative cases. Full suite: 30,437 passed; final finite-field rerun:
  12,465 passed including six later whitespace cases. Fast units: 5,888 pass, 32 existing skips;
  946 affected slow tests pass; no new type-check diagnostics. Document Sage 10.3/10.9
  oracle adapters and remaining gaps.

## 0.0.25 - 2026-09-09

- Compare all 99 public parent/element members in both prime-field implementations and
  the polynomial-basis extension backend; replace two formula-only comparative dispatchers.
- Repair Integer/Rational/boolean/default field construction, Integer/boolean arithmetic
  and equality, incompatible extension conversions, foreign polynomial coefficient rings,
  integer-representation bounds, and exact error classes/messages.
- Match Sage's Python random stream for finite fields and the bundled reference's newer
  extension-field default draw. Record the Sage 10.3 oracle adaptation explicitly.
- Add PARI Fp_order's word and recursive algorithms, use it for prime-field element order,
  and share it with elliptic-curve group code. Add comparative tests for both paths.
- Validate 29,103 full-suite comparisons plus twelve new matrix cases, 734 affected slow
  tests and 109 focused PARI tests. Both fast-tier timeouts pass in isolation; no new
  type-check diagnostic kinds.
- Correct a stale matrix test: change_ring applies explicit entry conversion even when no
  ring morphism exists. Add live comparisons for quotient-ring matrices converted to GF(p).

## 0.0.24 - 2026-09-09

- Put Integer and IntegerRing on a shared Sage-compatible coercion path: base-aware strings,
  whitespace/prefix/underscore handling, numeric boundary errors, None, booleans, digit lists,
  Integer/Rational inputs, conversion hooks and modular lifts. Compare exact error messages.
- Make Rational's default and scalar factories accept Sage's zero, boolean and Integer-wrapper
  cases. Delegate QQ tuple components to Integer coercion and preserve the distinction between
  constructor-level None and a singleton list containing None.
- Add comparative regressions for every discrepancy, including arbitrary-size list bases,
  long string conversion, C int base overflow and the Q_to_Z extra-argument error.
- Update public examples, design/deviation records and IntegerRing coverage tracking.
- Validate all 18,044 comparative cases, 5,888 fast tests (32 existing skips), 734 affected
  slow tests and 32 public-example tests; existing type-check diagnostic kinds are unchanged.

## 0.0.23 - 2026-09-09

- Audit all 38 RationalField members; fix Sage's 100-element sample list, generator errors,
  integral tuple coercion, fractional prime bounds, random defaults and seeded draw order.
- Replace duplicated full-size binary root searches with GMP's precision-doubling Newton
  algorithm and extra-limb optimization over BigInt. Rational roots delegate to Integer roots.
  Preserve signed C exponent conversion, including OverflowError and minimum-negative behavior.
- Add comparative regressions for every discrepancy, large root powers and their neighbors;
  update public API examples, source/deviation records and member coverage tracking.
  Validation: 17,259 comparisons, 5,888 fast units and 734 affected slow units pass; 32
  existing unit skips and the repository type-check failures remain.

## 0.0.22 - 2026-09-09

- Compare all 84 public Rational methods, accessors and construction paths with Sage in
  5,757 shared cases, including binary64 bit patterns and exact valuation error messages.
- Fix rational and matrix conversion rounding, subnormal/overflow results, large logarithmic
  heights, zero valuations, gamma poles, zero denominator errors, simplest-rational float
  coercion, GMP string parsing and negative decimal rounding precision. Route Gaussian
  lattice numeric entries through the corrected coercion with six downstream regressions.
- Update public examples, source/deviation records and the coverage inventory. Symbolic
  rational-result branches and the broader repository audit remain open. Fast units: 5,886
  pass, 32 skips; affected slow tests: 734 pass after repairs and isolated timeout rerun.

## 0.0.21 - 2026-09-09

- Repair all nine September findings: custom field moduli, incompatible parents, polynomial
  square-free units, signed CRT, inverse Frobenius, additive field logarithms, negative integer
  powers and zero exact logarithms. Add permanent shared-case SageMath comparisons.
- Extend the audit to all 107 Integer methods (6,027 comparisons); fix ideal reduction,
  infinity results, reciprocal inversion, zero/base ordering, extended gcd and modular-inverse
  boundaries, combinatorial negative indices, large natural logarithms and characteristic-two
  residuosity. Replace the quadratic class-number table with the existing PARI backend.
- Make comparative exceptions check their classes and recognize empty messages. Repair the
  newly exposed invalid-curve-point and dependent-lattice-basis paths; fix curve equation strings.
- Update public API examples, deviation/source records and the reproducible coverage inventory.
  All 10,805 comparisons pass; fast units: 5,886 pass, 32 skipped; full initial slow run: 1,274 pass;
  affected slow files rerun after repairs: 246 pass. Measured production lines: 79.82%.
  Repository-wide API audit remains in progress; the existing type-check failures are not closed.

## 0.0.20 - 2026-09-09

- Audit core arithmetic, finite fields, polynomials and generic groups against the bundled
  SageMath source and a live SageMath oracle; record nine additional open discrepancies in
  `AUDIT-2026-09.md` and `DEVIATIONS.md`. No production implementation changes.
- Add an isolated differential audit with shared deterministic inputs and Python/TypeScript
  probes: 51 cases, 22 matches, 29 mismatches reproducing the findings. Its command exits 1
  while discrepancies remain and is separate from the ordinary test tiers.
- Verify the existing differential tiers: fast 539/539, slow 4134/4134, zero failures/errors;
  public API documentation tests: 32/32 pass. Synchronize package and lockfile workspace versions.

## 0.0.19 - 2026-08-04

**`bun run typecheck` now actually runs.** It was failing on a configuration fault that
had been masking the fact that `packages/sagemath-ts` had never been type-checked at all.
Three of the five packages are now clean; the fourth's remaining errors are real and
newly visible.

### The configuration fault

`packages/sagemath-ts` and `packages/zksecurity-cheatsheets` had no `tsconfig.json`, so
`tsc` fell back to the repo root's — `rootDir: "src"` with no `include`, which pulls in
every `.ts` in the tree (`tutorial/`, `tests/`, `playground/`) and emits TS6059 for each.
502 of the 531 workspace errors were this noise, and it exited before reaching anything
real. Both packages now have a scoped config matching their siblings; TS6059 count is 0.

### Fixed

- **44 test files imported `describe/expect/it` from `vitest`**, which is not a dependency
  of this repo anywhere. Bun aliases it to its own runner, so the tests did run — but
  `tsc` could not resolve the module, and the other 74 test files already used
  `bun:test`. All 44 migrated; the convention is now uniform at 118 files.
- **`parigp-ts` is at zero type errors** (was 29). Its elliptic tests mixed the two point
  representations — `EllipticPointFp` (`group.ts`, nullable coordinates plus a `boolean`
  flag) and `EllipticPoint` (`points.ts`, a discriminated union). `pari-tests.test.ts` now
  takes the ShortWeierstrass-family predicates from `point.js` under an explicit alias
  instead of borrowing `group.js`'s.
- **`mkpoint` now returns `EllipticPointFinite`**, a new exported alias for the finite
  branch of `EllipticPoint`. It can never produce the point at infinity, so declaring the
  full union forced callers to re-narrow before reaching `.x`/`.y`. Strictly narrower and
  assignable everywhere the union was.
- **`zksecurity-cheatsheets` is at zero errors in its own source.** `createPrimeFieldCurve`
  was calling `GF(p, { check: false })` — see below; the options object was silently being
  used as a generator name.

### `GF` is not what the docs said (corrected in `LLM.md`)

The publicly reachable `GF` — identical object at the package root and at
`sagemath-ts/rings/finite_rings` — is `GFExtended`, whose own comment describes it as an
alias kept "for backward compatibility in tests". Consequences, all verified:

- It returns `PrimeField` for prime `q` and `FiniteFieldExtension` for a prime power, not
  the single `FiniteFieldPrime` that 0.0.17's `LLM.md` claimed. Neither class is
  re-exported from any index, so the return type cannot be named without importing
  `finite_field_extension.ts` by path.
- Its second parameter is a generator **name**, not an options object. `GF(7n, { check:
  false })` does not disable the primality check.
- The prime-only `GF(order, { check })` in `finite_field_constructor.ts` is shadowed by
  the alias and unreachable through any entry point.
- `PrimeField` is not assignable to `FiniteFieldPrime`, the type `EllipticCurve` declares
  for its base field. The two prime-field classes have diverged. Invisible at runtime —
  `EllipticCurve(GF(101n), …)` works — but a type error anywhere the code is checked.

`LLM.md` also now names the two exports bound to different implementations depending on
entry point: `rational_reconstruction` (root vs `matrix`) and `order_from_multiple` (root
vs `schemes/elliptic_curves`).

### Known and now measurable

`packages/sagemath-ts` reports **1361 type errors — 216 in the shipped library, 1145 in
its tests**. These are pre-existing and were never surfaced, because the config fault meant
`tsc` never checked this package. They are not addressed here. The dominant cause is
variance in the ring-element hierarchy: `RingElement`'s methods take `RingElement`
parameters, so `FiniteFieldElement` and `PrimeFieldElement` do not satisfy the constraint,
which cascades through `PolynomialRing`, the matrix modules and the coding modules. Fixing
it is a design change, not a cleanup pass.

`bun run lint` also still fails (133 errors, 64 warnings), overwhelmingly
`noExplicitAny` and formatting. The 15 `noPrecisionLoss` hits are Lanczos gamma
coefficients in `complex_mpfr.ts` and are deliberate.

### Testing

- `tests/llm-doc.test.ts` grows to 34 tests / 122 assertions, pinning the `GF` alias
  behaviour (identity, prime-power dispatch, generator-name parameter) and asserting that
  exactly two names diverge between entry points.

## 0.0.18 - 2026-08-04

A **documentation audit** across the rest of the repo, prompted by the `LLM.md` rewrite in
0.0.17: if the one unreferenced, untested document had rotted, the question was what else
had. Every claim below was checked by executing it or reading the current source.

### The package root is a curated subset (newly documented)

The largest finding, and the one most likely to strand a consumer. `import { … } from
'sagemath-ts'` exposes 159 names, but the subpaths expose far more, and a missing name
fails at import time with `SyntaxError: Export named 'nth_prime' not found`:

| Subpath | Exports | Not re-exported at the root |
|---|---:|---:|
| `rings` | 257 | 241 |
| `matrix` | 239 | 226 |
| `arith` | 100 | 58 |
| `schemes/elliptic_curves` | 99 | 86 |
| `rings/polynomial` | 82 | **all 82** |
| `rings/finite_rings` | 43 | 35 |
| `crypto` | 35 | 22 |
| `stats` | 13 | 7 |

`nth_prime`, `binomial`, `factorial`, `fibonacci`, `kronecker`, `hilbert_symbol`,
`bernoulli`, `primitive_root`, `CRT`, `PolynomialRing` and `NumberField` are all
subpath-only. `LLM.md` now carries the table, the full `arith` gap list, and the rule that
the subpath import always works. No exports were added — that is an API decision, not a
documentation one.

`LLM.md` also now records that `IntegerLike` is not applied uniformly: ~37 exported `arith`
functions declare a bare `bigint` and never call `toBigInt`, so
`nth_prime(new Integer(5n))` throws `ValueError: nth prime not found` instead of returning
`11n`.

### Corrections

- `README.md` documented a `tests/unit/` directory that does not exist — unit tests are
  colocated as `packages/*/src/**/*.test.ts`. The structure block now shows the real tree,
  including `tutorial/` and `playground/`.
- `README.md` claimed 77 tutorial lessons in two places; the playground generates 90.
- `AGENTS.md` described the property-test layout as `python/rings/test_integer.py` +
  `typescript/rings/test_integer.ts`. Tests are organized by *area*: `cases/<area>.cases.json`
  plus `python/areas/<area>.py` and `typescript/areas/<area>.ts`.
- `AGENTS.md` and `DESIGN.md` both advertised `factor(n, options?: {algorithm})`. `factor`
  takes no options; both examples are renamed so they read as the illustrations they are.
- `DESIGN.md` told readers to wrap results with `ZZ(result)`, and its "TypeScript
  Limitations" section still used `F(3)`, `R(3n)`, `QQ(1n, 2n)` and
  `Matrix(F, [[F(1), …]])`. Rings are not callable — every one of those throws. Rewritten
  against the real `__call__` convention, with the `matrix()` factory and the reason the
  generic `Matrix<R>` works over `GF(p)` but not over `ZZ`.
- `DESIGN.md`'s overload example used a `sqrt(n, all)` that does not exist; it now shows the
  real `sqrt_mod(a, p, all_roots)` shape.
- `DESIGN.md` pointed at a DEVIATIONS.md section called "No JavaScript Number Coercion";
  the section is "Language and Type-System Adaptations".
- `DEVIATIONS.md` located five `DESIGN.md` cross-references by line number, all shifted by
  the edits above. Replaced with section names. Its "Callable parents" row also implied
  `QQ()` was a factory helper; `ZZ`/`QQ` are singleton ring objects.
- Repointed four dangling `@see Deviation:` tags at their real sections:
  `no-number-coercion` (`types/coercion.ts`, and its test),
  `Number Field Kernel Ported Locally Instead of parigp-ts` (`pari_nf.ts`,
  `number_field_embeddings.ts`), `Finite Fields — Conway Table and Minimal Polynomials`
  (`finite_field_extension.ts`), and `port-only API, not a SageMath method`
  (`goppa_code.ts`).

### Verified as accurate, no change

`SOURCES.md` (all five upstream commits match the checkouts), `SCOPE.md` (current through
0.0.16, with honest 🟡/🔴 rows), `DEVIATIONS.md` content, `docs/PORTING_GUIDE.md`,
`docs/style-guide.md`, and every command in `README.md`/`AGENTS.md` including
`bun test --filter`. The benchmark tables match `tests/bench/results/*.json` exactly and
carry their generation date, so they are old but not misleading.

### Testing

- `tests/deviation-refs.test.ts` guards `@see Deviation:` tags against the rot just fixed.
  It does not require every tag to name a section — many are deliberately inline prose —
  but flags kebab-case slugs that match no heading, and tags that are a near-miss for one
  (real renames score ~0.71 similarity; the closest prose tag scores 0.53; threshold 0.65).
  Verified to fail when either bug shape is reintroduced.
- `tests/llm-doc.test.ts` grows to 30 tests / 111 assertions, pinning the export-gap table,
  the `matrix()`-over-a-field-but-not-`ZZ` split, and the `nth_prime(Integer)` trap.

## 0.0.17 - 2026-08-03

The **LLM.md rewrite**: the API quick reference that agents and vendored consumers read
had drifted far enough to mislead. It is now generated from executed code and pinned by a
test, so it cannot silently rot again.

### Documentation corrections

`LLM.md` was an orphan — nothing linked it and no test touched it. Ten of its examples
were wrong. Every claim below was verified by running it:

- **Installation was fabricated.** `bun add sagemath-ts` and `npm install sagemath-ts`
  both 404; the package is unpublished and the root workspace is `private: true`. The doc
  now documents the only real path: `git clone` + `bun install`.
- **Import specifiers.** Added a table covering all three consumer positions. The bare
  `sagemath-ts` alias resolves from anywhere inside the checkout;
  `@sagemath-ts/sagemath-ts` (the real package name) resolves *only* from within
  `packages/sagemath-ts/`, not from the repo root or sibling packages. Consumers outside
  the tree — vendored copies, read-only Docker mounts — must import
  `<checkout>/packages/sagemath-ts/src/index.ts` by path, since no bare specifier can
  resolve from there. Subpath entry points and the `exports` allowlist are now listed.
- **`factor()` is not trial division.** It delegates to `@sagemath-ts/parigp-ts`, a port
  of PARI's `ifactor1.c` with the full cascade — trial division, SQUFOF, Pollard-Brent
  rho, ECM and MPQS. Replaced the omission with a measured size/time table: ~56 ms at 128
  bits, ~1.5 s at 160, ~32 s at 200. Reaching for an external CAS is only warranted past
  that range.
- **Ring construction.** `GF(p)`, `Zmod(n)` and `ZZ` are not callable — every ring uses
  `.__call__(x)`. The old `F(5n)` examples all threw `TypeError`.
- **Corrected signatures**: `xgcd` returns `[g, s, t]` not an object; `factor` exponents
  are `bigint`; `crt(a, b, m, n)` takes four scalars (`CRT_list` takes the arrays);
  `sigma(n, k)` needs a bigint `k`; `is_prime_power` returns a boolean unless passed
  `get_data`; `bsgs` is positional; `lllReduce`/`IntegerLattice` take an `IntegerMatrix`
  or `bigint[][]`, not `vector()`s; `LWE`/`Regev` are positional oracles sampled with
  `.call()`; `DiscreteGaussianInteger(sigma)` is a factory, and the class behind it takes
  an options object; Reed-Solomon length must divide `q - 1`.
- **Removed fabrications**: `Regev.keygen/encrypt/decrypt` do not exist (Regev is Sage's
  LWE parameter set, not an encryption scheme); `Matrix(...).det()` and `.inverse()` do
  not exist — integer linear algebra goes through `IntegerMatrixFromEntries` and
  `determinant()`; `toBigInt('123')` and `toRational('3/4')` throw, as strings are not
  coerced; `gcd(12, 8)` throws, since `IntegerLike` excludes `number` by design.
- Documented where `number` *is* accepted — non-integer parameters, array/matrix
  dimensions, ring constructors, and ring-element arithmetic — which the old blanket
  "functions accept `number`" claim got backwards.
- `DESIGN.md`'s "Ring Coercion" section described a callable-class pattern that was never
  implemented; it now describes the real `__call__` convention and the element-level
  `number` coercion.

### Testing

- New `tests/llm-doc.test.ts` executes every `LLM.md` example (27 tests, 90 assertions),
  including the negative cases the doc relies on: rings not being callable, `number` and
  string rejection, the absent `det`/`inverse`/`keygen` methods, and a factoring budget
  that a trial-division regression could never meet. Runs in the fast unit tier.
- `AGENTS.md` gains a same-commit requirement for `LLM.md`, and `README.md` now links it,
  so it is no longer an unreferenced file.

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
