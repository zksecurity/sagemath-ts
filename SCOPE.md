# Project Scope

This document tracks implementation progress. Update this file when completing modules.

## Status Legend

- ✅ 2026-09-29 (Codex): general prime-model point-order delegation, curve/point
  caches and general-model PARI ellcard/ellorder (25.2.0). Eleven compact rows;
  1,820 expanded native comparisons, 565 existing caller comparisons, 185 point/
  kernel tests and 381 docs/storage checks pass. Eight builds pass; all 568 full
  baseline diagnostics unchanged. Other finite-field backends, public generic
  cardinality/groups, hybrid order scheduling and full elliptic audit remain open.

- 🟡 2026-09-29 (Codex): full elliptic-curve audit remains active; investigation
  covers remaining finite point/group backends, hybrid order scheduling,
  torsion enumeration and the remaining coordinate-parent domains.

- ✅ 2026-09-29 (Codex): general prime-model PARI scalar delegation and cached
  model conversion (25.1.0). Nine compact regressions/controls; 809 expanded native
  comparisons, 1,225 existing caller comparisons, 132 point/kernel and 381 docs/
  storage checks pass. Eight builds pass; all 568 baseline diagnostics unchanged.
  Point order, cardinality/groups and small-characteristic/extension scalar paths
  remain open under the full elliptic audit.

- ✅ 2026-09-29 (Codex): constructor model preservation and coefficient coercion
  (25.0.0). Nine compact regressions/controls; 2,466 main and 505 follow-up live
  comparisons, 193 caller tests and 380 docs/storage checks pass. Eight builds pass;
  two old type errors removed, 568 unchanged. Five-coefficient default calls now
  use the generic point API; full finite-field backend integration remains open.

- ✅ 2026-09-29 (Codex): prime finite-curve group-based point enumeration,
  exact sorted list, immutable cache and corrected-generator caching (24.63.0).
  Four compact regressions/controls; 304 fresh comparisons, 95 finite-curve and
  379 docs/storage tests pass. Eight builds pass; 570 baseline type diagnostics
  unchanged. Generic/extension torsion and broader elliptic work remain open.

- ✅ 2026-09-29 (Codex): p-primary basis dependency, division-point ordering,
  reduced polynomials, signed indices and cached-order propagation (24.62.0).
  Sixteen compact regressions/controls; 3,172 initial plus 516 final live
  comparisons and 786 caller/docs/storage tests pass. Eight builds pass; all
  570 baseline type diagnostics unchanged. Full elliptic audit stays in progress.

- ✅ 2026-09-29 (Codex): native bivariate interpolation/subresultants and public
  norm/charpoly delegation (24.61.0). Twelve compact regressions/controls; 2,420
  live comparisons, 570 caller/docs/storage checks and eight builds pass. No new
  type diagnostics; 69 existing missing-member counts reflect the charpoly API.
  Composed-sum Newton/Laplace routing and the broader audit remain open.

- ✅ 2026-09-29 (Codex): extension trace now delegates to FF_trace and native
  binary/word/large-prime derivative-remainder kernels (24.60.0). Eight compact
  regression/control rows; 908 live comparisons and 538 caller/docs/storage tests
  pass, as do eight builds. All 570 full baseline type diagnostics are unchanged.
  Norm dependencies were completed in the subsequent 24.61.0 batch.

- ✅ 2026-09-28 (Codex): hyperelliptic scalar roots and first-root selection
  repaired in odd-degree models, Cantor reduction and binary lifting (24.59.1).
  Seven compact regressions; 207 live dependency/state comparisons, all 121
  existing hyperelliptic tests and 41 case/storage checks pass. Eight builds pass;
  all 570 full baseline type diagnostics are unchanged. Same-parent root coverage
  does not close parent promotion or custom-field fallback routing.

- ✅ 2026-09-28 (Codex): optimized prime-curve coordinate callers now reuse
  inherited Sage logic for scalar coercion, promotion and extensions (24.59.0).
  Twelve compact regressions; 1,854 live coordinate comparisons (including scalar
  dependency traces), 103 polynomial comparisons, 95 existing curve tests and
  375 docs/storage checks pass. Constructor/model conversion remains open.

- ✅ 2026-09-28 (Codex): lift_x extend=True construction, ordering and returned
  parents over QQ, prime fields and explicit finite extensions (24.58.0). Preserve
  rational-polynomial variable names and expose the number-field zero protocol.
  Sixteen small regressions plus fresh live generators; 1,242 coordinate and 103
  polynomial comparisons, 245 caller/dependency tests and 373 docs/storage checks
  pass. Eight builds pass; all 570 full baseline type diagnostics are unchanged.
  Remaining algebraic coordinate parents, further lifting, global factory identity
  and native constructor backend/RNG parity remain open; see TODO.md.

- ✅ 2026-09-18 (Codex): generic lift_x canonical coercion and finite-field
  promotion; whole-coefficient base change and unchanged-parent curve identity
  (24.57.2). Twenty new small regressions plus fresh live generators cover values,
  errors, parents and constructor arguments. All 3,500 historical lift comparisons,
  929 focused fresh comparisons, 307 smaller advanced-area comparisons, 298 caller
  tests and 374 focused/docs/storage checks pass. Eight builds pass; fourteen old
  type diagnostics are removed and the remaining 570 full messages are unchanged.
  Retire 151,903 bytes of completed coordinate research records. The larger area
  run hit its native 120-second deadline; extend=True, real coordinates, other
  parents and the broader audit remain open in TODO.md.

- ✅ 2026-09-18 (Codex): generic curve is_x_coord input conversion and
  specialized GF(2) scalar constructors (24.57.1). Match native conversion/error
  handling with nine small comparative regressions and fresh shared generators.
  All 1,750 historical predicate comparisons, 309 focused fresh comparisons,
  221 advanced-area comparisons and 740 focused/caller/docs/storage tests pass.
  Eight builds pass; all 584 baseline type diagnostics (including continuation
  messages, allowing shifted source positions) are unchanged. lift_x canonical
  coercion/promotion and the other TODO.md coverage gaps remain open.

- ✅ 2026-09-14 (Codex): bounded command summaries and streaming log inspection
  (24.57.0). Six CLI tests and strict TypeScript checks pass, including multi-MB
  lines, early failures, replay settings and process-group cancellation. A real
  318,225-byte typecheck log produces a 2,200-character summary; all 584 existing
  diagnostics and exit status are preserved. AGENTS.md defaults to this wrapper.

- ✅ 2026-09-14 (Codex): seeded live runner and native output removal (24.55.0).
  Replace 86 snapshots in 87 test modules; fingerprint all 249,220 original inputs.
  55,890 reference replays, 5,729 focused tests, 41 fresh arithmetic properties,
  two 217-row live seed runs and framework checks pass. All 584 baseline type
  diagnostics remain identical. The behavioral audit remains paused.
- ✅ 2026-09-14 (Codex): fresh constrained generators replace 2,162,729 legacy
  area inputs (24.56.0), with explicit user acceptance of lost unclassified inputs.
  97 native profiles run live; retain 24 named inputs/compact constructions and
  failure seeds, remove unused bulk audit records, and replace seed/ordinal-based
  test dispatch with input checks. 59 generator/replay/storage and 164 final targeted
  native checks pass; two fresh
  arithmetic runs each pass 340 comparisons, eight builds pass, and all 584 type
  diagnostics remain unchanged. Open randomized behavior mismatches and the slow
  function-field area are recorded in TODO.md; the implementation audit is paused.

- ✅ 2026-09-14 (Codex): consolidate 265 September 9–14 commits above
  `72b1e25` into one commit (24.54.3), preserving current implementation and
  every regression file. Retain the old commit-message inventory and paused
  TODO handoff; verify tree contents before pruning superseded Git history.

- ✅ 2026-09-13 (Codex): stateful NTL modular polynomial products. Add 677
  shared comparisons, two unit checks and three executed API examples; fix 99
  old state mismatches. Full bundled NTL builds and passes its own checks.
  Focused 4,223 and fast 199,212 tests pass (32 existing skips); eight builds
  pass and all 584 type diagnostics are unchanged. Version 24.36.0; combined
  raw execution 110019/134067 (82.06%).
- ✅ 2026-09-13 (Codex): native integer polynomial exact division. Add 2,870
  native comparisons, one ownership check and three executed API examples; fix
  181 value/status and 172 first-call state discrepancies. Native diagnostics
  confirm skipped primes, reconstruction retries and rejection paths. Focused
  10,793 and fast 202,086 tests pass (32 skips); builds pass and all 584 type
  diagnostics remain unchanged. Version 24.37.0; raw execution 110123/134143 (82.09%).
- ✅ 2026-09-13 (Codex): native GotThem factor-recombination acceptance.
  Add 1,669 shared comparisons, three unit checks and two executed API examples.
  Focused 8,771 and fast 203,760 tests pass (32 skips); eight builds pass and all
  584 baseline type diagnostics are unchanged. Version 24.38.0; execution 110188/134208 (82.10%).
- ✅ 2026-09-13 (Codex): cold word-context initialization and small-prime
  factor recovery. Add 3,331 native comparisons, two unit checks and two executed
  API examples; fix 2,786 cache/stream discrepancies including 864 factor-output
  or retained-information discrepancies. Focused 12,318 and fast 207,095 tests
  pass (32 skips), builds pass and all 584 type diagnostics are unchanged.
  Version 24.39.0; combined raw execution 110206/134223 (82.11%).
- ✅ 2026-09-13 (Codex): shared state through local factor information,
  cardinality search, complementary products and exact division. Add 4,344 shared
  comparisons, one ownership check and three executed API examples; repair 112
  old trace differences. Focused 12,084 and fast 211,443 tests pass (32 skips),
  builds pass and all 584 type diagnostics are unchanged. Version 24.40.0;
  combined raw execution 110217/134233 (82.11%).
- ✅ 2026-09-13 (Codex): native integer product/square cache state and
  ChooseSS predicate. Add 1,497 shared comparisons, two unit checks and three
  executed API examples; fix 69 earlier cache/random-byte discrepancies. Native
  diagnostics confirm all multiplication/squaring and SS-selection branches.
  Focused 8,736 and fast 212,945 tests pass (32 skips); builds pass and all 584
  type diagnostics are unchanged. Version 24.41.0; raw execution 110250/134267 (82.11%).
- ✅ 2026-09-13 (Codex): word-modulus/multiplier failed-rebuild state,
  native truncated transforms and exact/fused CRT. Add 5,210 native fixtures
  covering 5,186 distinct shared rows,
  three unit checks and four executed API examples. Focused 5,523 tests pass;
  broad 218,159 pass (32 skips), and three timeouts pass unchanged isolated checks.
  Builds pass and all 584 type diagnostics are unchanged.
  Version 24.43.0; raw execution 110541/134581 (82.14%).
- ✅ 2026-09-13 (Codex): transposed word projection with retained quotient
  transforms, native error order and the defined plain right-shift boundary.
  Add 1,603 native fixtures (1,579 distinct shared rows), one ownership check and
  two executed API examples. Existing/new projection 10,461, final focused 7,165
  and fast 219,802 tests pass (32 skips); two timeouts pass unchanged isolated
  checks. Eight builds pass and all 584 baseline
  type diagnostics are unchanged. Version 24.44.0; raw execution 110637/134683 (82.15%).
- ✅ 2026-09-13 (Codex): constant-modulus remainder nontermination. Add
  240 native boundary comparisons and a watched subprocess regression that fails
  before the fix and passes afterward. Focused 18,048 tests pass; eight
  builds pass and all 584 baseline type diagnostics are unchanged. Version 24.44.1;
  raw execution 110638/134676 (82.15%). The broader quotient audit continues below.
- ✅ 2026-09-13 (Codex): cached-word multiplier fallback FFT validation.
  Add 309 native comparisons, fixing 36 differing traces and 68 operation results.
  Focused 18,357 tests pass; eight builds pass and all 584 baseline type diagnostics
  are unchanged. Version 24.44.2; raw execution 110644/134682 (82.15%).
- ✅ 2026-09-13 (Codex): arbitrary-modulus quotient/multiplier caches and
  native CRT through the 800-prime algorithm switch. Add 1,499 native fixtures
  (1,458 distinct shared rows), fixing 28 incompatible-prime-count rebuild traces;
  include 360 watched constant-modulus traces, one ownership check and two API
  examples. Focused 2,134 pass; broad 221,256 pass/32 skips with one timeout that
  passes unchanged in isolation. Eight builds pass and all 584 baseline type
  diagnostics are unchanged. Version 24.45.0; current-file execution 572/630,
  aggregate recorded execution 111008/135043 (82.20%).
- ✅ 2026-09-13 (Codex): word multiplier allocation and cross-context cached
  coefficient fidelity. Add 828 distinct shared/native traces (740 word, 88 big),
  fixing 200 missing word allocation errors and 28 numeric traces per backend.
  Regressions fail before each repair. Live 828/828 and focused 37,754 tests pass;
  eight builds pass and all 584 baseline type diagnostics are unchanged. Version
  24.45.1; word execution 478/478, big 572/630, aggregate 111010/135045 (82.20%).
- ✅ 2026-09-13 (Codex): native word ordinary products and extended GCD, with
  word inverse delegation and context propagation through negative powers.
  Add 780 distinct shared/native traces, one ownership check and two executed
  API examples; restoring old inverse delegation fails 23 traces / 85 operations.
  Live 780/780 and focused 72,026 tests pass; eight builds pass and all 584 baseline
  type diagnostics are unchanged. Version 24.46.0; word modules execute 506/506
  and 454/454 lines, aggregate recorded execution 111086/135104 (82.22%).
- ✅ 2026-09-13 (Codex): shared state and native error order through Hensel
  lifting. Add 362 shared/native traces, two executed API examples and an ownership
  check; the old implementation differs on 267 traces / 882 operations. Live
  362/362 and eight builds pass, with 584 unchanged baseline type diagnostics.
  Full fast tier: 223,231 pass / 32 skips / two timeouts, both passing unchanged
  in isolation. Version 24.47.0; edited module 1144/1144 lines, aggregate recorded
  execution 111101/135118 (82.23%).
- ✅ 2026-09-13 (Codex): native FFT-prime word coefficient contexts and shared
  integer GCD/squarefree state. Add 519 shared traces, fixing 491 old mismatches;
  native diagnostics exercise all twelve targeted reconstruction outcomes. Live
  519/519, broad NTL 156,758, Sage callers 143, and five coefficient-only fallback
  comparisons pass. Eight builds pass and all 584 baseline type diagnostics remain
  unchanged. Version 24.48.0; aggregate raw execution 111137/135197 (82.20%).
- ✅ 2026-09-13 (Codex): complete NTL integer factorization, adaptive van Hoeij
  recombination, deflation and the remaining number-field constructor degree
  routing. Add 377 shared comparisons, three executed API examples and an
  ownership check. Live 377/377, fast 224,144 tests (32 skips) and eight builds
  pass; all 584 baseline type diagnostics are unchanged. Native diagnostics
  exercise every targeted lattice outcome. Version 24.49.0; factor-driver
  execution 1418/1418, aggregate raw execution 111396/135262 (82.36%).
- ✅ 2026-09-13 (Codex): legacy integer factor helpers now delegate to PARI,
  preserving factor order and random state and removing their 200-bit limit.
  Add 224 shared/native regressions; pre-fix tests reproduce 90 small-input
  mismatches and two size-limit errors. Live 224/224, focused 21,665 tests and
  eight builds pass; all 584 baseline type diagnostics are unchanged. Version
  24.49.1; edited module 1402/1465 lines, aggregate 111282/135135 (82.35%).
- ✅ 2026-09-13 (Codex): unembedded number-field automorphism ordering,
  signed integral scaling and the natural-map swap. Add 190 comparisons using
  the bundled Sage model and complete native PARI; 94 fail before repair and
  all pass afterward. Live 190/190, broader 11,372 tests and eight builds pass;
  584 baseline type diagnostics are unchanged. Version 24.49.2; current module
  1544/1661 lines, aggregate recorded execution 111298/135150 (82.35%).
- ✅ 2026-09-13 (Codex): native conjugate-count bound and word derivative/
  squarefree helpers. Add 1,500 shared/native comparisons and five executed
  API expressions; fix 40 negative-start values and a watched -2^63 scan hang.
  Live 1,500/1,500, fast 226,059 tests (32 skips), the separate 62-test Galois
  suite and eight builds pass; all 584 baseline type diagnostics are unchanged.
  Version 24.50.0; Flx 193/193, galconj 2363/2395, aggregate raw execution
  111293/135116 (82.37%). Full conjugate-kernel delegation remains open.
- ✅ 2026-09-14 (Codex): integer/rational polynomial roots. Source-based dense/
  sparse dispatch, FLINT block GCD, NTL/PARI factorization, zero errors and
  root ordering now match 699 live comparisons (18 watched large-coefficient
  cases). Broader 1,625 tests and eight builds pass; 584 baseline diagnostics
  unchanged. Full fast 226,757 pass/32 skips/one existing NTL timeout, whose
  unchanged isolated rerun passes in 1.85s. Version 24.50.1; root helpers 87/87
  recorded lines, polynomial_element 5075/5306, aggregate 82.77% raw execution.
- ✅ 2026-09-14 (Codex): public ZZ/QQ factor and irreducibility dispatch,
  content order, constant/zero results and cached predicates. Add 1,118 shared
  comparisons (304 PARI-state traces); all match. Broader polynomial/caller/doc
  checks pass 4,917 tests. Eight builds pass; all 584 baseline type diagnostics
  unchanged. Version 24.50.2; polynomial_element 5072/5940 measured lines,
  aggregate 111247/135036 (82.38%). Raw counters include unexecuted non-code lines.
- ✅ 2026-09-14 (Codex): finite-field root extraction, distinct-root options,
  elliptic division-point caller and extension factor comparison. All 6,007
  new native comparisons match, including 236 PARI-state traces and the watched
  degree-5,100 caller. Broader 10,211 tests and eight builds pass; 584 existing
  type errors remain, with 39 member-count details reflecting removal of an
  unused private helper. Version 24.51.0; execution 111281/134035 (83.02%).
- ✅ 2026-09-14 (Codex): IntegerModRing root contracts: 5,126 new live
  comparisons and the complete 9,199-case modular area match. Fast 239,012
  tests pass (32 existing skips), affected slow callers 434 pass, final focused
  tests 5,496 pass. Eight builds pass; 584 baseline type errors remain with
  type-display ordering changes only. Version 24.52.0; raw execution
  111416/134573 (82.79%). A runner regression preserves both lift dispatchers.
- ✅ 2026-09-14 (Codex): Weierstrass generator root equations and argument
  validation (24.52.1): 1,759 new native comparisons and 32 existing Weierstrass
  comparisons pass. Initial/auxiliary port discrepancies: 577. Final focused
  2,146 tests pass; eight builds pass, all 584 type diagnostics unchanged.
  The complete 2,009-case advanced area matches native records. Broader caller
  validation passes 3,055 tests (13 existing skips); raw execution is
  111328/134058 (83.04%), with weierstrass_morphism 356/422 recorded lines.
- ✅ 2026-09-14 (Codex): generic-curve isomorphism selection, ordering and
  parent guards (24.52.2): 2,371 new native comparisons, including comparison
  matrices and early-exit traces; all 4,380 advanced-area comparisons match.
  Broader callers/docs/case format: 5,428 pass, 13 existing skips. Eight builds
  pass and all 584 baseline type diagnostics are unchanged. Raw execution
  111330/133673 (83.29%); remaining scalar-root callers are still under audit.
- ✅ 2026-09-14 (Codex): scalar extension square roots and PARI
  FF_issquareall (24.53.0): 4,324 distinct native comparisons and 32,970
  complete-area comparisons match. Final callers/docs/case format: 15,341
  pass, 13 existing skips. Eight builds pass; 584 existing type errors remain
  with 32 member-detail changes from the new sqrt method. Raw execution
  111392/133602 (83.38%); new dependency 117/117 recorded lines.
- ✅ 2026-09-14 (Codex): finite-field square predicates and native PARI
  resultant/norm dependencies (24.54.0). Add 4,832 distinct comparisons; all
  37,802 field-area rows match. Broader callers/docs/case format pass 26,555
  tests (13 existing skips); eight builds pass. The 584 existing type errors
  remain, with member details updated for is_square. Raw execution records
  111577/133764 (83.41%); full-port behavioral coverage remains incomplete.
- ✅ 2026-09-14 (Codex): generic curve coordinate-root repairs (24.54.1).
  Add 5,822 shared comparisons; all 10,202 advanced-area rows match. Final
  broader callers/docs/case format pass 21,571 tests (13 existing skips).
  Eight builds pass and all 584 baseline type diagnostics are unchanged.
  Raw execution 111541/133728 (83.41%); generic curve 922/960 lines. General
  coordinate coercion/promotion and remaining caller audits continue.
- ⏸ 2026-09-14 (Codex): behavioral audit paused at user request. Remaining
  work and the preserved partial coercion patch are indexed in TODO.md.
- ⏸ 2026-09-14 (Codex): downstream elliptic/hyperelliptic roots and finite-field
  irreducibility remain open; resume from TODO.md.
- ✅ 2026-09-14 (Codex): test-storage cleanup (24.54.2). All 52 migrated files
  decompress to their original hashes: 2,304,091 case/fixture records retained.
  Native replay checks pass 15,853 tests; selected property replay passes;
  canonical/storage checks and eight builds pass. All 584 type diagnostics are
  unchanged. TODO.md and compressed pending-work records preserve the audit
  handoff; generated transcripts are removed and future output defaults to gzip.


- ✅ 2026-09-11 (Codex): word-prime polynomial factor delegation and FLINT dependencies:
  5,359 new shared comparisons against bundled C/Sage; live area 189,701 pass plus
  four final constant-boundary cases; caller/docs/units 236,034 pass; final kernel
  replay 1,015 pass; fast 6,283 pass/32 skips; all 588 type diagnostics unchanged.
  Native FLINT module tests and cached source build pass. Combined execution:
  94438/118850 (79.46%); versions 20.4.0. Other factor backends remain under audit.
- ✅ 2026-09-11 (Codex): PARI random-state and polynomial sampling dependency:
  2,408 new bundled-C comparisons; live 192,113 pass; callers/docs/units 238,447
  pass; fast 6,287 pass/32 skips; 588 type diagnostics unchanged; bundles pass.
  Native random suites and temporary-cache build pass. New kernels execute
  105/105 lines; combined 94546/118958 (79.48%); versions 20.5.0.
  Generic finite-field factor routing and unrelated local generators remain open.
- ✅ 2026-09-11 (Codex): binary PARI factor routing and arithmetic/kernel dependencies:
  11,139 new comparisons; live 203,252 pass; callers/docs/units 249,591 pass;
  fast 6,292 pass/32 skips; 588 type diagnostics unchanged; bundles pass.
  Native factor suites pass. New arithmetic/kernel files execute 196/196 lines;
  combined 94862/119066 (79.67%); versions 20.6.0. Other factor routes remain open.
- ✅ 2026-09-11 (Codex): quotient power entries and binary automorphism powering:
  3,775 new bundled-PARI comparisons fix 845 observed mismatches; live 207,027
  pass; caller/docs/units 253,114 pass; Galois callers 62 pass; fast 6,295 pass/
  32 skips; 588 type diagnostics unchanged; bundles pass. Combined execution:
  94881/119180 (79.61%); versions 20.6.1. Modular composition remains open.
- ✅ 2026-09-11 (Codex): PARI matrix multiplication dependencies and dispatch:
  1,269 new comparisons; live 208,296 pass; callers/docs/units 254,384 pass;
  Buchmann/Galois callers 97 pass; fast 6,298 pass/32 skips; native mat/linear
  suites pass; 588 type diagnostics unchanged; bundles pass. Combined execution:
  95090/119396 (79.64%); versions 20.7.0. Modular composition is next.
- ✅ 2026-09-11 (Codex): PARI Brent–Kung composition, trace and signed products:
  15,018 new comparisons; live 223,314 pass; callers/docs/units 269,404 pass;
  arithmetic/Galois/modular-polynomial callers 160 pass; fast 6,302 pass/32 skips;
  native ff/pol suites pass; 588 type diagnostics unchanged; bundles pass.
  Combined execution: 95167/119532 (79.62%); versions 20.8.0.
- ✅ 2026-09-11 (Codex): PARI word/integer polynomial multiplication and squaring:
  2,710 new bundled-native comparisons; live 226,024 pass; caller/docs/units
  272,116 pass; arithmetic/Galois callers 160 pass; fast 6,307 pass/32 skips.
  All 588 type diagnostics unchanged; bundles pass. Execution: 95398/119766
  (79.65%); versions 20.9.0. Initial broad timeouts pass without timeout changes.
- ✅ 2026-09-11 (Codex): PARI quotient division/reduction and Barrett inverse:
  8,541 new comparisons, including 630 old discrepancies; live 234,529 pass
  plus 36 final boundary controls; caller/docs/units 280,657 pass; arithmetic/
  Galois callers 160 pass; fast 6,311 pass/32 skips. All 588 type diagnostics
  unchanged; bundles pass. Division helper: 135/135 executed lines; combined
  95623/119991 (79.69%); versions 20.10.0. Full behavioral coverage remains open.
- ✅ 2026-09-11 (Codex): PARI unscaled GCD/Bézout, recursive half-GCD and
  dependent inverse/Hensel/split-part normalization; quotient power schedules:
  10,613 new comparisons; live 245,178 pass; caller/docs/units 291,272 pass;
  arithmetic/Galois callers 160 pass; schedule replay 931 pass; fast 6,316 pass/
  32 skips. All 588 type diagnostics unchanged; bundles pass. New GCD helper:
  138/138 lines; combined 95815/121120 (79.11%); versions 20.11.0.
- ✅ 2026-09-11 (Codex): native randomized Shoup minimal polynomials and
  shared reciprocal reduction: 6,203 new comparisons; live 251,381 pass;
  caller/docs/units 297,478 pass; arithmetic/Galois callers 160 pass; fast
  6,319 pass/32 skips. All 588 type diagnostics unchanged; bundles pass.
  New kernel executes 71/73 lines (native defensive restart remains unhit);
  combined 95854/121194 (79.09%); versions 20.12.0.
- ✅ 2026-09-11 (Codex): native distinct-degree factorization/counts and shared
  quotient contexts: 4,165 new comparisons; cached-native caller replay 301,648
  pass; fresh helper replay 2,941 pass; arithmetic/Galois callers 160 pass; fast
  6,324 pass/32 skips. All 588 type diagnostics unchanged; bundles pass.
  Combined execution: 95977/121323 (79.11%); versions 20.13.0.
- ✅ 2026-09-11 (Codex): PARI modular square-root algorithms and native
  boundaries: 21,551 new comparisons fix 402 observed discrepancies (eight
  hangs). Cached-native caller replay 323,203 pass; fast 6,328 pass/32 skips;
  extended arithmetic/elliptic callers 329 pass, including P-256/Curve25519.
  Native arith/sqrtn suites pass; all 588 type diagnostics unchanged; bundles
  pass. New helper 177/182 lines; combined 96108/121382 (79.18%); v20.14.0.
- ✅ 2026-09-11 (Codex): native full PARI factorization/equal-degree splitting,
  monic normalization and large-prime Sage delegation: 19,566 new shared cases.
  Cached-native caller replay 342,679 pass, plus 96 final error/state controls;
  fast 6,334 pass/32 skips; arithmetic/Galois callers 256 pass. Native factorff/
  ff/pol suites and bundles pass; all 588 type diagnostics unchanged. New helpers:
  133/133 lines; combined 96307/121595 (79.20%); versions 20.15.0.
- ✅ 2026-09-11 (Codex): native polynomial roots, total-splitting predicates and
  word root-count shortcuts: 17,709 new shared comparisons fix 91 baseline
  discrepancies. Cached-native caller replay: 360,490 pass; fast: 6,340 pass/
  32 skips; arithmetic/Galois callers: 256 pass. Native pol suites and bundles
  pass; all 588 type diagnostics unchanged. Root helper: 170/170 lines;
  combined: 96436/121735 (79.22%); versions 20.16.0.
- ✅ 2026-09-11 (Codex): PARI extension-polynomial products, squares, reduction
  and normalization: 10,575 new shared comparisons, including Node/Bun large-array
  regressions. New/native/docs replay: 10,721 pass; existing binary/quotient replay:
  15,064 pass; fast: 6,346 pass/32 skips; arithmetic/Galois callers: 256 pass.
  Native ff/factorff/pol suites and bundles pass; 588 type diagnostics unchanged.
  Shared helper: 195/195 lines; combined: 96672/122020 (79.23%); versions 20.17.0.
- ✅ 2026-09-11 (Codex): extension polynomial division, cached reduction and
  reciprocal inversion: 27,367 new shared comparisons; combined arithmetic/
  division/docs replay: 38,097 pass; binary/quotient: 15,064 pass; fast: 6,355 pass/
  32 skips; arithmetic/Galois callers: 256 pass. Native ff/factorff/pol suites and
  bundles pass; all 588 type diagnostics unchanged. Division/coefficient helpers:
  306/306 lines; combined: 96961/122318 (79.27%); versions 20.18.0.
- ✅ 2026-09-11 (Codex): extension polynomial GCD, unscaled Bézout coefficients
  and recursive half-GCD: 24,532 new comparisons, including 6,696 previously
  rejected native binary packing cases and isolated upstream cycle controls.
  Arithmetic/division/GCD/docs: 62,634 pass; binary/quotient: 15,064 pass; fast:
  6,360 pass/32 skips; arithmetic/Galois callers: 256 pass. Native ff/factorff/pol
  suites and bundles pass; 588 type diagnostics unchanged. GCD helper: 177/177
  lines; combined: 97188/122546 (79.31%); versions 20.19.0.
- ✅ 2026-09-11 (Codex): 21 extension quotient operations, native inverse display
  and reusable inner reciprocals: 24,023 new comparisons. Combined native replay
  and docs: 86,645 pass; binary/regressions: 15,081 pass; fast: 6,365 pass/32 skips;
  arithmetic/Galois callers: 256 pass. Native suites and bundles pass; 588 type
  diagnostics unchanged. Quotient helper: 135/135 lines; combined execution:
  97451/122820 (79.34%). Versions 20.20.0.
- ✅ 2026-09-11 (Codex): six extension composition APIs and packed word-extension
  matrix multiplication: 22,386 new comparisons. All values/errors match; the
  109,052-test aggregate had one transient timeout, passing with 39 neighboring
  cases under the unchanged limit. Fast: 6,369 pass/32 skips; binary: 15,064 pass;
  arithmetic/Galois: 256 pass; native suites/build pass; 588 type diagnostics
  unchanged. New helpers: 161/161 lines; combined: 97640/123016 (79.37%).
  Versions 20.21.0; see AUDIT-2026-09.md for the timing evidence.
- ✅ 2026-09-11 (Codex): public prime-field reciprocal preparation and canonical
  storage: 8,130 new comparisons cover all 29 saved discrepancies and controls.
  Focused: 8,286 pass; broad native replay: 469,210 pass without timeouts;
  fast: 6,373 pass/32 skips; arithmetic/Galois callers: 256 pass. Native suites
  and bundles pass; all 588 type diagnostics unchanged. Combined execution:
  97660/123041 (79.37%); versions 20.21.1.
- ✅ 2026-09-11 (Codex): automorphism reciprocal preparation and canonical
  storage: 7,161 permanent comparisons repair all 1,532 reproduced mismatches.
  Focused: 34,243 pass; fast: 6,378 pass/32 skips; arithmetic/Galois: 256 pass.
  Native suites and bundles pass; 588 type diagnostics unchanged. Combined
  execution: 97685/123069 (79.37%); versions 20.21.2.
- ✅ 2026-09-11 (Codex): eight coefficient-substitution/binary composition APIs:
  11,437 new comparisons; combined checks: 56,383 pass; binary: 15,064 pass;
  fast: 6,383 pass/32 skips; arithmetic/Galois: 256 pass. Native suites and
  bundles pass; 588 type diagnostics unchanged. New helpers: 54/54 and 28/29
  instrumented lines (one closing brace remains in the denominator). Combined:
  97778/123185 (79.37%); versions 20.22.0.
- ✅ 2026-09-11 (Codex): eight extension automorphism power/trace/aggregate APIs:
  10,686 new comparisons; combined checks: 67,073 pass; binary: 15,064 pass;
  fast: 6,388 pass/32 skips; arithmetic/Galois: 256 pass. Native suites/builds
  pass; 588 type diagnostics unchanged. New helper: 97/97 executed lines;
  combined: 97943/123354 (79.40%); versions 20.23.0.
- ✅ 2026-09-11 (Codex): eight extension random, dot-product and truncated-product
  APIs: 20,281 new comparisons; combined checks: 112,578 pass; binary: 15,064
  pass; fast: 6,393 pass/32 skips; arithmetic/Galois: 256 pass. Native suites/
  builds pass; 588 type diagnostics unchanged. Projection helper: 61/61 lines;
  extension arithmetic including truncation: 197/197. Combined: 98078/123483
  (79.43%); versions 20.24.0.
- ✅ 2026-09-11 (Codex): two native extension minimal-polynomial APIs with Shoup
  projection: 5,695 new comparisons; combined checks: 118,278 pass; binary:
  15,064 pass; fast: 6,398 pass/32 skips; arithmetic/Galois: 256 pass. Native
  suites/builds pass; 588 type diagnostics unchanged. New helper: 94/96 lines
  (two native defensive-restart lines unhit); combined: 98190/123603 (79.44%).
  Versions 20.25.0.
- ✅ 2026-09-11 (Codex): prime quotient cache/inverse repairs: 1,753 observed
  discrepancies fixed, 16,228 new comparisons and 2,820 strengthened comparisons.
  Combined checks: 163,185 pass; fast: 6,405 pass/32 skips; arithmetic/Galois:
  256 pass. Native suites/builds pass; 588 type diagnostics unchanged. New helper:
  46/46 lines; combined: 98225/123611 (79.46%). Versions 20.25.1.
- ✅ 2026-09-11 (Codex): six native Frobenius APIs: 6,195 new comparisons;
  combined checks: 124,449 pass; fast: 6,412 pass/32 skips; arithmetic/Galois:
  256 pass. Native suites/builds pass; 588 type diagnostics unchanged. Wrappers:
  16/16 lines; helper: 84/87 (three zero-hit instrumented comments retained).
  Combined: 98332/123736 (79.47%); versions 20.26.0.
- ✅ 2026-09-11 (Codex): prime-polynomial linear boundary repairs: 5,280
  observed discrepancies fixed with 18,414 new shared comparisons. Combined
  checks: 187,794 pass; fast: 6,417 pass/32 skips; arithmetic/Galois: 256 pass.
  Native suites/builds pass; 588 type diagnostics unchanged. Combined execution:
  98348/123737 (79.48%); versions 20.26.1.
- ✅ 2026-09-11 (Codex): nine extension root-count/split/squarefree/derivative
  APIs and scalar inverse/residue repair: 274 actual baseline discrepancies fixed,
  15,230 new shared comparisons; frozen replay 203,027 pass; fast 6,424 pass/
  32 skips; callers 256 pass. Native suites/build pass; 588 type diagnostics
  unchanged. Combined: 98446/123834 (79.50%); versions 20.27.0.
- ✅ 2026-09-11 (Codex): signed polynomial word dispatch: 1,902 observed
  division/GCD/half-GCD discrepancies fixed with 28,490 new native comparisons.
  Frozen replay 231,516 pass; fast 6,429 pass/32 skips; callers 256 pass. Native
  suites/build pass; 588 type diagnostics unchanged. Combined: 98457/123842
  (79.50%); versions 20.27.1.
- ✅ 2026-09-11 (Codex): signed composition, power tables and matrix boundaries:
  44,838 new native comparisons, 18 refreshed guard cases and cached trace repairs.
  Frozen replay 276,374 pass; fast 6,436 pass/32 skips; callers 256 pass. Native
  suites/build pass; 588 type diagnostics unchanged. Combined: 98461/123845
  (79.50%); versions 20.27.2.
- ✅ 2026-09-11 (Codex): polynomial observation helpers: 3,379 reproduced
  discrepancies fixed with 10,703 new native comparisons. Frozen replay 287,077
  pass; fast 6,442 pass/32 skips; Galois/arithmetic callers 256 pass. Native suites/
  build pass; 588 type diagnostics unchanged. Combined: 98482/123866 (79.51%);
  versions 20.27.3.
- ✅ 2026-09-11 (Codex): exported scalar arithmetic: 6,460 reproduced records
  repaired with 21,712 new native comparisons. Frozen replay 308,848 pass; fast
  6,448 pass/32 skips; arithmetic callers 256 pass; extended elliptic tests 73 pass.
  Native suites/build pass; 588 type diagnostics unchanged. Combined execution:
  98491/123860 (79.52%); versions 20.27.4.
- ✅ 2026-09-11 (Codex): scalar exponentiation: 479 reproduced records repaired
  with 12,137 new native comparisons, including nine stack-boundary controls.
  Frozen replay 320,985 pass; fast 6,453 pass/32 skips; arithmetic callers 256 pass;
  extended elliptic tests 73 pass. Native suites/build pass; 588 type diagnostics
  unchanged apart from equivalent union ordering. Combined: 98608/123979 (79.54%); versions 20.27.5.
- ✅ 2026-09-12 (Codex): scalar predicates/order: 5,872 reproduced records
  repaired with 28,782 new native comparisons. Frozen replay 349,769 pass; fast
  6,459 pass/32 skips; arithmetic callers 256 pass; extended elliptic tests 73 pass.
  Native suites/build pass; 588 type diagnostics unchanged. Combined: 98622/123987 (79.54%);
  versions 20.27.6.
- ✅ 2026-09-12 (Codex): rational Galois polynomials: 3,771 reproduced records
  repaired with 7,003 new native comparisons. Frozen replay 356,715 pass; fast
  6,465 pass/32 skips; arithmetic/Galois callers 256 pass. Native arithmetic and
  galoisinit suites/build pass; optional native galois suite skips without galdata.
  Type diagnostics remain 588. Combined: 98634/123999 (79.54%); versions 20.27.7.
- ✅ 2026-09-12 (Codex): Galois integer helpers: 242 baseline failures repaired,
  including 27 time-limit cases, with 1,291 new native comparisons. Frozen replay
  358,007 pass; fast 6,471 pass/32 skips; arithmetic/Galois callers 256 pass. Native
  arithmetic/galoisinit suites and build pass; 588 type diagnostics unchanged.
  Combined: 98705/124072 (79.55%); versions 20.27.8.
- ✅ 2026-09-12 (Codex): integer-polynomial observations and indexpartial: 855
  baseline failures repaired, 8,161 native comparisons, and source-based strict
  partial factorization/modular GCD/p-adic echelon dependencies. Frozen replay
  366,171 pass; fast 6,479 pass/32 skips; callers 256 pass. Native suites/build
  pass; 588 type diagnostics unchanged. Combined: 98975/124223 (79.68%); versions 20.27.9.
- ✅ 2026-09-12 (Codex): native Hensel root/factor lifting and Bezout propagation:
  2,234 baseline failures repaired; 13,016 direct native comparisons and 42
  precondition-guard checks. Frozen replay 379,231 pass; fast 6,488 pass/32 skips;
  callers 256 pass. Native suites/build pass; 588 type diagnostics unchanged.
  Combined: 99099/124353 (79.69%); versions 20.27.10.
- ✅ 2026-09-12 (Codex): Vandermonde interpolation and batch inverses: 1,792
  baseline error discrepancies repaired; 4,227 native comparisons. Frozen replay
  383,456 pass; fast 6,494 pass/32 skips; callers 256 pass. Native suites/build
  pass; 588 type diagnostics unchanged. Combined: 99155/124407 (79.70%); versions 20.27.11.
- ✅ 2026-09-12 (Codex): permutation/vector helpers and word orders: 1,926
  defined-execution failures repaired; 10,696 native records/eight guard checks.
  Frozen replay 394,163 pass; fast 6,502 pass/32 skips; callers 256 pass. Native
  suites/build pass; 588 type diagnostics unchanged. Combined: 99191/124442 (79.71%); 20.27.12.
- ✅ 2026-09-12 (Codex): subgroup enumeration, cosets and quotients: six error
  mismatches repaired; 3,515 native comparisons across 46 group presentations.
  Frozen replay 397,674 pass; fast 6,505 pass/32 skips; callers 256 pass. Build
  passes; 588 type diagnostics unchanged. Combined: 99191/124438 (79.71%); 20.27.13.
- ✅ 2026-09-12 (Codex): fixed-field symmetric evaluation: 586 baseline failures
  repaired; 2,571 native comparisons including unsigned-power dependency controls.
  Frozen replay 400,249 pass; fast 6,512 pass/32 skips; callers 256 pass. Build
  passes; 588 type diagnostics unchanged. Combined: 99163/124438 (79.69%); 21.0.0.
- ✅ 2026-09-12 (Codex): unit-subgroup enumeration and exact-index dependencies:
  555 baseline failures repaired; 3,505 native comparisons. Frozen replay 403,755
  pass; fast 6,519 pass/32 skips; callers 256 pass. Build passes; 588 type
  diagnostics unchanged. Combined: 99521/124798 (79.75%); 21.0.1.
- ✅ 2026-09-12 (Codex): scalar-bound/unbounded subgroup modes and strict partial
  factorization: 2,781 native comparisons; the explicit bound-mode gap is closed.
  Frozen replay 406,536 pass; fast 6,525 pass/32 skips; callers 256 pass. Build
  passes; 588 type diagnostics unchanged. Combined: 99678/124940 (79.78%); 21.0.2.
- ✅ 2026-09-12 (Codex): Galois validation/conjugate results: 130 discrepancies
  fixed and 322 native comparisons. Replay 406,858 pass; fast 6,530 pass/32 skips;
  callers 256 pass. Build passes; 588 type diagnostics unchanged. 99683/124940 (79.78%); 21.0.3.
- ✅ 2026-09-12 (Codex): Galois fixed-field/subgroup operations: 355 discrepancies
  fixed with 2,084 native comparisons. Replay 408,942 pass; fast 6,534 pass/32 skips;
  callers 256 pass. Build passes; 588 type diagnostics unchanged. 99687/124925 (79.80%); 21.0.4.
- ✅ 2026-09-12 (Codex): symmetric-polynomial search bounds: 17 failing records
  repaired; 557 native comparisons. Replay 409,499 pass; fast 6,538 pass/32 skips;
  callers 256 pass. Build passes; 588 type diagnostics unchanged. 99687/124925 (79.80%); 21.0.5.
- ✅ 2026-09-12 (Codex): Galois helper coverage: 33 totient discrepancies fixed;
  4,838 native comparisons, including stateful filters and Frobenius selection.
  Replay 414,339 pass; fast 6,543 pass/32 skips; callers 256 pass. Build passes;
  588 type errors unchanged. 99785/124938 (79.87%); 21.0.6.
- ✅ 2026-09-12 (Codex): finite-field kernel algorithm/output fidelity: 168
  discrepancies fixed; 3,920 native comparisons. Replay 418,259 pass; fast 6,547
  pass/32 skips; callers 256 pass. Build passes; 588 type diagnostics unchanged.
  Combined: 100044/125150 (79.94%); 21.0.7.
- ✅ 2026-09-12 (Codex): factor-stage gates, budgets, seed arithmetic and retry
  errors: 277 mismatching records repaired; 1,888 native stage comparisons and 210
  fixed-field prime-replacement comparisons. Replay 420,361 pass; fast 6,555
  pass/32 skips; callers 308 pass. Build passes; 588 type diagnostics unchanged.
  Combined: 100044/125146 (79.94%); 21.0.8.
- ✅ 2026-09-12 (Codex): perfect-power signs, masks, search order and failure
  sentinel: 1,016 discrepancies repaired; 8,914 native comparisons. Replay 429,274
  pass; fast 6,560 pass/32 skips; all 23 slow files: 1,274 pass. Build passes;
  588 type diagnostics unchanged. Combined: 100672/125203 (80.41%); 21.0.9.
- ✅ 2026-09-12 (Codex): native ECM batching, PRAC and continuation: 181
  discrepancies repaired; 800 native comparisons. Replay 430,077 pass; fast 6,567
  pass/32 skips; 22 affected slow files: 1,247 pass. New ECM module: 409/409
  instrumented lines. Build passes; 588 type diagnostics unchanged.
  Combined: 100941/125458 (80.46%); 21.0.10.
- ✅ 2026-09-12 (Codex): k-th-power errors/bounds and MPQS roots/relations:
  957 mismatching records repaired; 2,732 native comparisons. Replay 432,809 pass;
  fast 6,573 pass/32 skips; MPQS 27 and callers 308 pass. Build passes; 588 type
  diagnostics unchanged. Combined: 100927/125468 (80.44%); 21.0.11.
- ✅ 2026-09-12 (Codex): sparse binary kernels and MPQS block Lanczos:
  420 mismatching records repaired; 766 native comparisons. Replay 433,576 pass;
  fast 6,579 pass/32 skips; MPQS 27 and callers 308 pass. Build passes; 588 type
  diagnostics unchanged. Combined: 101078/125618 (80.46%); 21.0.12.
- ✅ 2026-09-12 (Codex): MPQS Kronecker symbols and debug relation checks:
  1,427 mismatching records repaired; 4,802 native comparisons. Replay 438,381
  pass; fast 6,587 pass/32 skips; MPQS 27 and callers 308 pass. Build passes;
  588 type diagnostics unchanged. Combined: 101054/125589 (80.46%); 21.0.13.
- ✅ 2026-09-12 (Codex): MPQS inverse errors and initialization/sieve controls:
  200 error discrepancies repaired; 1,093 native comparisons. Replay 439,472
  pass; fast 6,592 pass/32 skips; MPQS 27 and callers 308 pass. Build passes;
  588 type diagnostics unchanged. Combined: 101075/125588 (80.48%); 21.0.14.
- ✅ 2026-09-12 (Codex): MPQS relation hashing and final factor/RNG output:
  45 mismatching records repaired; 230 native comparisons. Replay 439,703 pass;
  fast 6,597 pass/32 skips; MPQS 27 and callers 308 pass. New module: 65/65
  instrumented lines. Build passes; 588 type diagnostics unchanged.
  Combined: 101137/125650 (80.49%); 21.0.15.
- ✅ 2026-09-12 (Codex): MPQS candidate reconstruction and warning semantics:
  216 mismatching records repaired; 483 native comparisons. Replay 440,186 pass;
  fast 6,601 pass/32 skips; MPQS 27 and callers 308 pass. Build passes;
  588 type diagnostics unchanged. Combined: 101144/125652 (80.50%); 21.0.16.
- ✅ 2026-09-12 (Codex): Pollard lambda signed hashing, CPython steps and bounds
  errors: 1,038 discrepant records repaired; 1,340 new and all 1,675 group-area
  comparisons pass. Focused 30,671; fast 6,605 pass/32 skips. Build passes;
  588 type diagnostics unchanged. Combined: 101197/125374 (80.72%); 21.0.17.
- ✅ 2026-09-12 (Codex): Pollard rho ring draws, native shortcuts and errors:
  887 discrepant records repaired; 1,304 new and all 2,979 group-area comparisons
  pass. Focused 31,981; fast 6,611 pass/32 skips. Build passes; 588 diagnostics
  unchanged. Combined: 101202/125356 (80.73%); 21.0.18.
- ✅ 2026-09-12 (Codex): Shared group parsing, custom multiples and exact bounds:
  1,053 discrepant records repaired; 1,620 new and all 4,599 group-area comparisons
  pass. Focused 33,606; fast 6,616 pass/32 skips; elliptic callers 212 pass.
  Build passes; 588 diagnostics unchanged. Combined: 101291/125629 (80.63%); 21.0.19.
- ✅ 2026-09-12 (Codex): Order-from-multiple factor lists and valuations:
  173 discrepant records repaired (63 former port hangs); 480 new and all 5,079
  group comparisons pass. Focused 34,092; fast 6,622 pass/32 skips; elliptic
  callers 212 pass. Build passes; 588 diagnostics unchanged. Combined: 101313/125651 (80.63%);
  21.0.20. Native nontermination is recorded and excluded from execution cases.
- ✅ 2026-09-12 (Codex): Generic multiple-iterator copying and state:
  576 discrepant records repaired; 608 new and all 5,687 group comparisons pass.
  Focused 34,709; fast 6,631 pass/32 skips; elliptic callers 212 pass.
  Build passes; 588 diagnostics unchanged. Combined: 101362/125690 (80.64%); 21.0.21.
- ✅ 2026-09-12 (Codex): Number-field scalar/group interoperability:
  746 discrepant records repaired; 821 new and all 6,508 group comparisons pass.
  Focused 35,537; fast 6,638 pass/32 skips; number-field/elliptic callers 355 pass.
  Build passes; 588 diagnostics unchanged. Combined: 101377/125656 (80.68%); 21.0.22.
- ✅ 2026-09-12 (Codex): Number-field scalar coercion and zero division:
  850 discrepant records repaired; 1,530 new and all 1,554 number-field comparisons
  pass, plus all 6,508 group records. Focused 37,100; fast 6,647 pass/32 skips;
  callers 355 pass. Build passes; 588 diagnostics unchanged. Combined: 101392/125656 (80.69%);
  21.0.23. Exact scalar equality retains generic-group type compatibility.
- ✅ 2026-09-12 (Codex): Number-field coefficient indexing:
  1,332 discrepant records repaired; 1,638 new and all 3,192 number-field comparisons
  pass, plus all 6,508 group records. Focused 38,747; fast 6,656 pass/32 skips;
  callers 355 pass. Build passes; 588 diagnostics unchanged. Combined: 101443/125712 (80.69%);
  21.0.24. Original quadratic scaling and bounded discriminant normalization retained.
- ✅ 2026-09-12 (Codex): Integer/rational arithmetic and group interoperability:
  333 discrepant records repaired; 1,680 new comparisons, all 8,230 Integer and
  6,892 group records pass. Focused 47,370; fast 6,665 pass/32 skips; callers 355
  pass. Build passes; 588 diagnostics unchanged. Combined: 101652/125699 (80.87%); 21.0.25.
  Typed arithmetic overloads and equality against this preserve concrete numeric types.
- ✅ 2026-09-12 (Codex): Rational equality with real inputs:
  812 discrepant records repaired; 1,170 new and all 7,472 rational comparisons pass.
  Focused 54,848; fast 6,671 pass/32 skips; callers 355 pass. Build passes;
  588 diagnostics unchanged. Combined: 101679/125657 (80.92%); 21.0.26. Comparisons delegate to RDF.
- ✅ 2026-09-12 (Codex): RDF scalar arithmetic and standard binary multiple:
  1,789 discrepant records repaired; 2,922 new comparisons pass, together with all
  278,897 polynomial and 7,288 group records. Focused 57,779; fast 6,680 pass/32 skips;
  callers 355 pass. Build passes; four old type errors removed, 584 remain unchanged.
  Combined: 101670/125647 (80.92%); 21.0.27.
- ✅ 2026-09-12 (Codex): Generic group schedules, BSGS collision equality and RDF errors:
  1,575 discrepant records repaired; 3,604 new and all 10,892 group comparisons pass.
  Focused 61,393; fast 6,690 pass/32 skips; callers 355 pass. Build passes; 584 type
  diagnostics unchanged. Combined: 101640/125592 (80.93%); 21.0.28.
- ✅ 2026-09-12 (Codex): Generic group parent identities and equality dispatch:
  795 discrepant records repaired; 3,870 new and all 14,762 group comparisons pass.
  Focused 65,270; fast 6,697 pass/32 skips; callers 355 pass. Build passes; 584 type
  diagnostics unchanged. Combined: 101621/125550 (80.94%); 21.0.29.
- ✅ 2026-09-12 (Codex): Pollard parent/mutability handling, collision equality and
  singleton bounds: 504 discrepant records repaired; 692 new and all 15,454 group
  comparisons pass. Focused 65,972; fast 6,707 pass/32 skips; callers 355 pass.
  Build passes; 584 type diagnostics unchanged. Combined: 101598/125518 (80.94%); 21.0.30.
- ✅ 2026-09-12 (Codex): Number-field coefficient reduction, defining equations and
  integer powers: 1,218 discrepant records repaired; 6,002 new comparisons, including
  288 scaled-caller controls. All 8,594 field and 600 new bundled-FLINT records match.
  Focused 71,987; fast 6,720/32 skips; callers 417 pass. Build passes; 584 existing
  type diagnostics unchanged. Combined: 101741/125419 (81.12%); version 21.0.31.
- ✅ 2026-09-12 (Codex): General ideal inversion, maximal-order differents and ideal
  powers delegate through PARI trace/HNF dependencies. 3,574 new comparative cases;
  866 prepatch mismatches (278 wrong behaviors and 588 missing-implementation records).
  All 10,674 field and 2,419 matrix-area records match fresh native runs. Focused
  76,499; fast 6,733/32 skips; callers 417 pass; build passes; 584 type diagnostics
  unchanged. Combined 102486/126053 (81.30%); version 21.0.32.
- ✅ 2026-09-12 (Codex): Ideal integral bases, zero modules and two-generator
  representations use HNF/PARI and cache the native results. 1,878 new comparisons;
  1,652 discrepant records (910 wrong behaviors and 742 missing implementations).
  All 12,552 field records match fresh Sage. Focused 75,965; fast 6,740/32 skips;
  callers 417; build pass; 584 type diagnostics unchanged. Combined 102483/126020 (81.32%).
  Version 22.0.0 corrects gens_two's first-entry type to NumberFieldElement.
- ✅ 2026-09-12 (Codex): Fractional ideal arithmetic, numerator/denominator ideals,
  zero dispatch and factory normalization: 1,595 wrong-behavior records repaired;
  9,090 new comparisons, including 1,215 direct PARI division controls. All 21,642
  field records match fresh Sage. Focused 85,064; fast 6,749/32 skips; callers 417;
  build passes; 584 type diagnostics unchanged. Combined 102507/126010 (81.35%).
  Version 23.0.0 corrects denominator's return type to NumberFieldIdeal.
- ✅ 2026-09-12 (Codex): Native PARI fast LLL stage and binary64 dependencies:
  5,345 new comparisons (1,970 stages; 3,375 scalar primitives), 17 unit tests and
  one executed API example. Full matrix area 7,764 matches; focused 7,979; fast
  6,767/32 skips; build passes; 584 type diagnostics unchanged. New files 272/272
  lines; combined 102779/126594 (81.19%). Native compiler profile and 94 timeout probes explicit.
- ✅ 2026-09-12 (Codex): PARI DPE stage: 3,588 exact comparisons plus nine separately
  normalized resource-boundary records; 26 unit tests and an executed API example.
  Full matrix area 11,361 matches; focused 11,603; fast 6,794/32 skips; build passes;
  584 type diagnostics unchanged. LLL file 433/435 emitted lines (two closing-call
  delimiters unhit; both expression branches execute). Combined 102999/126816 (81.22%).
  Keep 75 native timeout probes outside passing counts. Version 23.0.2.
- ✅ 2026-09-12 (Codex): PARI heuristic/proved real LLL stages and scalar dependencies:
  13,338 new exact native comparisons, 42 unit tests and an executed API example.
  Full matrix area 24,699 matches; focused 24,984; fast 6,837/32 skips; build passes;
  584 type diagnostics unchanged. LLL file 793/795 emitted lines (two closing-call
  delimiters retained); new scalar code fully executed. Combined 103370/127185 (81.28%).
  Twelve native/TypeScript bounded probes retained outside passing counts. Version 23.0.3.
- ✅ 2026-09-12 (Codex): PARI Householder QR, upper R and Gaussian reduction:
  4,752 exact native comparisons, 23 unit tests and an executed API example.
  Full matrix area 29,451 matches; focused 29,675; fast 6,861/32 skips; build passes;
  584 type diagnostics unchanged. bibli1 executes 123/123 emitted lines; combined
  103460/127269 (81.29%). Version 23.0.4. Full behavioral coverage remains open.
- ✅ 2026-09-12 (Codex): real PARI Cholesky, triangular inverse and numeric branches:
  10,445 exact native comparisons, 39 unit tests and an executed API example.
  Full matrix area 39,896 matches; focused 40,137; fast 6,901/32 skips; build passes;
  584 type diagnostics unchanged. All three changed files execute 154/154 emitted lines.
  Combined 103529/127335 (81.30%). Nineteen native undefined-payload probes retained separately.
  Version 23.0.5; generic coefficient domains and broader audit remain open.
- ✅ 2026-09-12 (Codex): PARI integer/real dot, Gram, matrix product and rescaling:
  5,903 exact native comparisons, 39 unit tests and an executed API example.
  Full matrix area 45,799 matches; focused 46,041; fast 6,941/32 skips; build passes;
  584 type diagnostics unchanged. Changed files execute 117/117 emitted lines.
  Combined 103644/127430 (81.33%). Version 23.0.6; broader audit remains open.
- ✅ 2026-09-12 (Codex): PARI adaptive QR/Cholesky and precision statistics:
  1,230 exact native comparisons, 32 unit tests and an executed API example.
  Full matrix area 47,029 matches; focused 47,265; fast 6,974/32 skips; build passes;
  584 type diagnostics unchanged. New helper 107/108 emitted lines (closing loop
  delimiter retained); LLL 794/796. Combined 103752/127539 (81.35%). Version 23.0.7.
- ✅ 2026-09-12 (Codex): PARI modular pivots, integer rank and rational matrix solve:
  5,418 exact native comparisons (2,802 word / 2,616 integer), 42 unit tests and
  an executed API example. Matrix area 52,447 matches; focused 52,694; fast 7,017
  with 32 skips; builds pass; 584 type diagnostics unchanged. Changed files execute
  750/753 lines; combined 104068/127850 (81.40%). Version 23.0.8.
- ✅ 2026-09-12 (Codex): typed native adaptive LLL and FLATTER wrapper, including
  a multiplication-order performance fix in all four stages. 3,935 exact native
  cases; 78 unit tests; three separate native-timeout controls; executed API example.
  Matrix area 56,382 matches; focused 56,727; fast 7,096/32 skips; builds pass;
  584 type diagnostics unchanged. Wrapper 320/322, stages 817/819; unreachable
  return/closing delimiters remain counted. Combined 104464/128162 (81.51%). Version 23.0.9.
- ✅ 2026-09-12 (Codex): native LLL-kernel ideal intersection, scalar/list ideal
  construction, fractional-class/zero/identity semantics and rational foreign-field
  coercion. 5,126 new shared cases, including 3,478 observed pre-fix mismatches;
  26,768 live number-field matches; 27,416 focused tests; fast 7,102/32 skips;
  builds and executed docs pass; 584 type diagnostics unchanged. Execution 104535/128028 (81.65%).
  Version 23.0.10. The 32 Sage-layer random-state controls also pass pre-fix.
- ✅ 2026-09-12 (Codex): ideal membership/operator coercion, exact-length field
  coefficient vectors, native centered integral bases, integer quotient rounding
  and primitive-content ideal multiplication. 9,176 new shared comparisons;
  35,944 live matches; 36,720 focused and 7,230 fast tests pass (32 existing skips);
  builds/examples pass; 584 type diagnostics unchanged. Execution 104657/128068 (81.72%).
  Version 23.0.11. Bounded class-number certificates retain their covered domain.
- ✅ 2026-09-12 (Codex): native Gram LLL, HNF kernel retention, Smith residues/
  limb thresholds/rounding, direct relation-basis reduction and quaternion callers.
  2,938 new shared comparisons; 58,456 matrix, 10,386 algebraic-dependency and
  855 quaternion comparisons pass. Focused replay 51,833; fast 9,528/32 skips;
  builds/examples pass; 584 type diagnostics unchanged. Execution 104577/128070 (81.66%).
  Version 23.0.12. Twenty-seven native/host resource failures are separate controls.
- ✅ 2026-09-12 (Codex): direct base-ideal method errors and fractional prime
  factory classes; 424 observed discrepancies fixed. 1,948 new shared comparisons,
  260 native unit fixtures and an executed API example. Full number fields 37,892,
  focused 38,805 and fast 9,789/32 skips pass; builds pass; 584 baseline type errors
  unchanged. Version 23.0.13. Combined execution 104604/128122 (81.64%).
- ✅ 2026-09-12 (Codex): native ideal valuation, prime finalization and rational
  prime lookup. Correct I.valuation(P), negative/infinite results, native validation
  precedence, and two captured hangs. Version 24.0.0 fixes the breaking signature.
  8,753 new shared comparisons; full number fields 46,645 match with old prefix
  unchanged; focused 49,961 and fast 12,282/32 skips pass. Builds/examples pass;
  584 type diagnostics unchanged. Combined execution 104754/128227 (81.69%).
- ✅ 2026-09-12 (Codex): full PARI prime decomposition, ideal factorization,
  maximality, modular HNF construction and finite-field matrix dependencies.
  7,137 new shared comparisons and 7,138 unit/example checks cover 273 old factor
  discrepancies and ten BigInt-growth failures. Native-backed areas: 51,830 number
  fields and 60,408 matrices; combined 122,642 pass; fast 19,420/32 skips pass.
  Native suites/builds/examples pass; 584 type diagnostics unchanged. Version 24.1.0;
  combined execution 105315/128779 (81.78%).
- ✅ 2026-09-12 (Codex): no-argument number-field trace/norm and PARI QQ
  quotient arithmetic. 3,561 new shared comparisons and 3,574 unit/example checks
  include 12 bounded child processes; five old three-second operation timeouts
  are repaired. All 55,391 number-field records match; combined 95,430
  and fast 22,994 pass (32 existing skips). Builds/examples/native pol
  suite pass; 584 type diagnostics unchanged. Version 24.2.0; execution 105490/128950 (81.81%).
- ✅ 2026-09-12 (Codex): trace/norm QQ, own-field and degree-one bases,
  result parents and argument errors. 276 old mismatches repaired; 672 new shared
  comparisons and 673 unit/example checks. Full number fields 56,063 match;
  focused 66,757 and fast 23,667 pass (32 existing skips). Builds/examples pass;
  584 type diagnostics unchanged. Version 24.3.0; execution 105557/129017 (81.82%).
- ✅ 2026-09-12 (Codex): adaptive PARI LLL norm output for polynomial
  factorization. 2,318 new shared probes: 2,161 exact result/error comparisons,
  66 documented resource boundaries and 91 bounded observations. Full matrix
  native prefix unchanged; 62,726 matrix and 56,063 number-field replay records.
  Focused 123,811 pass plus one unchanged isolated timeout replay (6.1 seconds);
  fast 25,896 pass/32 existing skips. Builds and native suites pass;
  584 type diagnostics unchanged. Version 24.4.0; execution 105574/129259 (81.68%).
- ✅ 2026-09-12 (Codex): 14 native PARI factorization dependencies: permuted
  HNF/knapsack recognition, LLL progress, coefficient/root bounds and scalar
  kernels. Fix five reproduced rational zero-denominator error mismatches.
  Add 17,092 shared comparisons and 17,099 unit/example checks; focused 46,543
  pass, fast 42,995 pass/32 existing skips. Builds/native suites pass; 584 type
  diagnostics unchanged. Version 24.5.0; execution 105787/129606 (81.62%).
- ✅ 2026-09-12 (Codex): native bounded integer-factor recombination, exact
  trial quotients, lifting precision and balanced modular products. Add 7,437
  shared comparisons and 7,441 unit/example checks. Focused 89,482 pass, then all
  220 docs pass after correcting two new import paths; fast 50,436 pass/32 skips.
  Builds/native suites pass; 584 type diagnostics unchanged. Version 24.6.0;
  execution 105947/129805 (81.62%). Larger unresolved factor groups are preserved.
- ✅ 2026-09-12 (Codex): native general Hermite forms, full transformations
  and zero-column removal. Preserve the original seven/eight-column dispatch;
  share exact column operations with permuted HNF. Add 11,200 shared native
  comparisons and 11,202 unit/example checks. Focused 32,459 and fast 61,638 pass
  (32 existing skips); builds/native suites pass; 584 type diagnostics unchanged.
  Version 24.7.0; execution 106009/129864 (81.63%); edited HNF module 463/465 lines.
- ✅ 2026-09-12 (Codex): full PARI integer factorization stages and Newton
  traces; number-field constructor routing for degree 30–300, with cached
  irreducibility results. Add 5,599 shared comparisons and 5,605 unit/example
  checks. Focused 84,850 and fast 67,219 pass (32 existing skips), plus 24 final
  knapsack controls. Builds/native suites pass; 584 type diagnostics unchanged.
  Version 24.8.0; raw execution 106329/130122 (81.71%).
- ✅ 2026-09-12 (Codex): eight native NTL trace/precision/lattice dependencies;
  4,289 shared original-body comparisons and 4,293 unit/example checks. Focused
  4,520 and fast 71,536 pass (32 existing skips); builds pass and 584 baseline
  type diagnostics are unchanged. Updated NTL module executes 96/96 lines.
  Version 24.9.0; combined raw execution 106409/130346 (81.64%).
- ✅ 2026-09-12 (Codex): native NTL multifactor Hensel lifting and portable
  exact integer products. Add 4,058 shared comparisons and 4,062 unit/example
  checks. Focused 8,582, fast 75,598 and 729 existing native callers pass (32
  existing skips). Builds pass; 584 baseline type diagnostics are unchanged.
  Version 24.10.0; raw execution 106540/130523 (81.63%).
- ✅ 2026-09-12 (Codex): exact NTL LLL/LLL_plus, integer image, lattice
  solving and extended-GCD dependency. Add 5,576 shared comparisons and 5,581
  unit/example checks. Focused 5,810 and fast 81,179 pass (32 existing skips);
  builds pass and 584 baseline type diagnostics are unchanged. New kernels
  execute 219/219 lines. Version 24.11.0; raw coverage 106759/130742 (81.66%).
- ✅ 2026-09-12 (Codex): NTL factor recovery bounds, lattice row removal and
  additional lifting. Add 3,848 shared comparisons and 3,852 unit/example checks.
  Focused 4,083 and fast 85,031 pass (32 existing skips); builds pass and 584
  baseline type diagnostics are unchanged. ZZXFactoring executes 278/278 lines.
  Version 24.12.0; raw coverage 106836/130817 (81.67%).
- ✅ 2026-09-12 (Codex): NTL degree patterns and factor selection helpers.
  Add 7,291 shared comparisons and 7,295 unit/example checks. Focused 7,528
  and fast 92,326 pass (32 existing skips); builds pass and 584 baseline type
  diagnostics are unchanged. ZZXFactoring executes 398/398 emitted lines.
  Version 24.13.0; raw coverage 106956/130937 (81.69%).
- ✅ 2026-09-13 (Codex): NTL word-polynomial quotient and power dependencies.
  Fix 135 initial-adapter mismatches and add 10,880 shared comparisons plus
  10,884 unit/example checks. Focused 11,119 and final fast 103,210 pass (32
  existing skips); builds pass and 584 baseline type diagnostics are unchanged.
  The new module executes 143/143 lines. Version 24.14.0; raw coverage
  107099/131080 (81.71%). An existing timeout passed isolated/full reruns.
- ✅ 2026-09-13 (Codex): native NTL word matrix multiplication schedule.
  Add 1,246 shared comparisons and 1,250 unit/example checks. Focused 1,487
  and final fast 104,460 pass (32 existing skips); builds pass and 584 baseline
  type diagnostics are unchanged. New code executes 111/111 emitted lines.
  Version 24.15.0; raw coverage 107210/131191 (81.72%). Existing constructor
  timeouts pass isolated/full reruns without timeout or implementation changes.
- ✅ 2026-09-13 (Codex): NTL word composition, trace maps and native GCD schedule.
  Fix 120 initial-adapter allocation-error mismatches, retain their comparative
  baseline and add 45 controls. Add 10,559 shared and 10,565 unit/example checks.
  Focused 10,804 and fast 115,025 pass (32 existing skips); builds pass and 584
  baseline type diagnostics are unchanged. New modules execute 217/217 lines.
  Version 24.16.0; raw coverage 107427/131408 (81.75%).
- ✅ 2026-09-13 (Codex): NTL word distinct-degree factorization and SFCanZass1.
  Add 2,481 shared comparisons and 2,485 unit/example checks. Focused 13,289
  and fast 117,510 pass (32 existing skips); builds pass and 584 baseline type
  diagnostics are unchanged. Updated factorization executes 211/211 lines.
  Version 24.17.0; raw coverage 107590/131571 (81.77%). No mismatch observed.
- ✅ 2026-09-13 (Codex): NTL sequence minimum-polynomial reconstruction.
  Add 3,245 shared comparisons and 3,249 unit/example checks. Focused 16,538
  and fast 120,759 pass (32 existing skips); builds pass and 584 baseline type
  diagnostics are unchanged. Updated composition/sequence module executes
  212/212 lines. Version 24.18.0; raw coverage 107633/131614 (81.78%).
- ✅ 2026-09-13 (Codex): NTL segmented PrimeSeq iterator.
  Add 1,031 shared operation traces and two executed examples. Focused 1,279
  and fast 121,792 pass (32 existing skips); builds pass and 584 baseline type
  diagnostics are unchanged. New code executes 80/80 lines. Version 24.19.0;
  raw coverage 107713/131694 (81.79%). No mismatch observed.
- ✅ 2026-09-13 (Codex): NTL cached multipliers and transposed projections.
  Add 8,432 native comparisons, 423 documented native-crash guard controls and
  8,860 unit/example checks. Focused 36,282 and final fast 130,652 pass (32 skips);
  builds pass and 584 baseline type diagnostics are unchanged. Updated modules
  execute 467/467 lines. Version 24.20.0; raw coverage 107825/131806 (81.81%).
  One existing timeout passes isolated/full reruns without code/deadline changes.
- ✅ 2026-09-13 (Codex): NTL deterministic byte streams and key derivation.
  Add 4,038 native records, 28 explicit short-buffer guard records and five
  ownership/API examples. Focused 5,352 and fast 134,723 pass (32 existing skips);
  builds pass and 584 baseline type diagnostics are unchanged. New source executes
  198/198 lines. Version 24.21.0; raw coverage 108023/132004 (81.83%).
- ✅ 2026-09-13 (Codex): NTL integer samplers with explicit stream contexts.
  Add 2,626 native traces, six signed-word guard controls and four ownership/API
  tests. Focused 7,988 and fast 137,359 pass (32 existing skips); builds pass and
  584 baseline type diagnostics are unchanged. New source executes 81/81 lines.
  Version 24.22.0; raw coverage 108104/132085 (81.84%).
- ✅ 2026-09-13 (Codex): NTL quotient-element minimum polynomials.
  Add 11,713 native result/error/stream traces and four ownership/API checks.
  Focused 15,218 and fast 149,076 pass (32 existing skips); builds pass and 584
  baseline type diagnostics are unchanged. Updated module executes 334/334 lines.
  Version 24.23.0; raw coverage 108155/132136 (81.85%).
- ✅ 2026-09-13 (Codex): NTL cached samplers, random polynomials and root products.
  Add 2,933 native records, 406 labeled adapter-guard traces and five ownership/API
  checks. Focused 22,018 and fast 152,420 pass (32 existing skips); builds pass and
  584 baseline type diagnostics are unchanged. All 82 new routine lines execute.
  Version 24.24.0; raw coverage 108238/132219 (81.86%).
- ✅ 2026-09-13 (Codex): NTL root splitting and equal-degree modular factor recovery.
  Fix 70 observed division-error discrepancies with 20 passing controls. Add 4,433
  native comparisons, 120 labeled guards and four ownership/API checks. Focused
  7,298 and fast 156,977 pass (32 existing skips); builds pass; 584 baseline type
  diagnostics unchanged. Updated word modules execute 543/543 lines.
  Version 24.25.0; raw coverage 108380/132361 (81.88%).
- ✅ 2026-09-13 (Codex): NTL degree-dependent word contexts and propagation.
  Add 6,500 native comparisons, six labeled numeric guards and four ownership/API
  checks; fix 66 prototype size-error discrepancies. Focused 6,772 and fast
  163,487 pass (32 existing skips); builds pass; 584 type diagnostics unchanged.
  Updated modules execute 1007/1007 lines. Version 24.26.0; overall raw execution
  108479/132460 (81.90%).
- ✅ 2026-09-13 (Codex): NTL stateful small-prime factor selection.
  Add 2,175 native stateful comparisons and four ownership/API checks; fix 102
  prototype signed-bit discrepancies. Focused 9,736 pass; fast 165,664 pass,
  32 skip and two timeouts. Both timeout cases pass unchanged in isolation.
  Builds pass; 584 type diagnostics unchanged. The edited
  module executes 494/494 lines. Version 24.27.0; overall raw execution
  108575/132556 (81.91%).
- ✅ 2026-09-13 (Codex): NTL retained local factor-information updates.
  Add 1,216 native stateful comparisons and four ownership/API checks. Focused
  10,956 and fast 166,886 pass (32 existing skips); builds pass; 584 type
  diagnostics unchanged. The edited module executes 570/570 lines. Version
  24.28.0; overall raw execution 108651/132651 (81.91%).
- ✅ 2026-09-13 (Codex): NTL cardinality-based integer-factor recombination.
  Add 2,879 native comparisons, five labeled guards and four ownership/API checks.
  Focused 13,844 and fast 169,774 pass (32 existing skips); builds pass; 584 type
  diagnostics unchanged. All 426 new emitted lines execute. Version 24.29.0;
  overall raw execution 109077/133107 (81.95%).
- ✅ 2026-09-13 (Codex): NTL probable-prime tests, probability bounds and
  sequential/seeded-block prime generators. Add 4,621 native comparisons, three
  labeled guards, one legacy-overload check and three executed API examples.
  Focused 12,631 and fast 174,402 pass (32 existing skips); builds pass; 584 type
  diagnostics unchanged. All 206 new emitted lines execute. Version 24.30.0;
  overall raw execution 109283/133298 (81.98%). Rare block retries and counter
  exhaustion are source-reviewed; full behavioral coverage remains open.
- ✅ 2026-09-13 (Codex): NTL word-matrix elimination, left kernel, inverse,
  determinant and left/right solving. Add 5,086 native comparisons, eight
  labeled guards, two ownership checks and four executed API examples.
  Focused 6,620 and fast 179,502 pass (32 existing skips); builds pass; all 584
  type diagnostics unchanged. All 235 new emitted lines execute; module 346/346.
  Version 24.31.0; overall raw execution 109518/133528 (82.02%).
- ✅ 2026-09-13 (Codex): NTL FFT-prime recognition, search and lazy root tables.
  Add 3,230 native comparisons, 29 labeled guards, two ownership checks and four
  executed API examples. Focused 8,166 and fast 182,767 pass (32 existing skips);
  builds pass; all 584 type diagnostics unchanged. All 147 emitted lines execute.
  Version 24.32.0; overall raw execution 109665/133675 (82.04%).
- ✅ 2026-09-13 (Codex): NTL scalar/matrix CRT, signed intervals, row-norm
  determinant bounds and arbitrary-modulus determinants. Add 8,881 native
  comparisons, 22 labeled guards, two ownership checks and four executed API
  examples. Focused 13,814 and fast 191,676 pass (32 existing skips); builds pass;
  all 584 type diagnostics unchanged. All 158 new emitted lines execute.
  Version 24.33.0; overall raw execution 109823/133833 (82.06%).
- ✅ 2026-09-13 (Codex): NTL integer determinant/inverse reconstruction.
  Add 3,595 native comparisons, two ownership checks and four executed API
  examples. Focused 12,790 and fast 195,277 pass (32 existing skips); builds
  pass; all 584 type diagnostics unchanged. New lines 109/114; post-certificate
  rescaling remains unexercised. Module 183/188; version 24.34.0; overall raw
  execution 109932/133962 (82.06%).
- ✅ 2026-09-13 (Codex): NTL certified integer elimination. Add 3,250 native
  comparisons, one ownership check and two executed API examples. Native
  diagnostics confirm each retry path and unchanged outputs. Focused 7,138
  and fast 198,530 pass (32 existing skips); builds pass; all 584 type diagnostics
  unchanged. All 46 new emitted lines execute. Version 24.35.0; overall raw
  execution 109978/134028 (82.06%).
- 🟡 2026-09-13 (Codex): constructor NTL integer-polynomial factorization
  outside degree 30–300: recombination acceptance checks, van Hoeij iteration
  and complete native factorization routing. Inverse rescaling branch remains
  under source/coverage review.
  Proper-subfield/morphism matrices, characteristic/minimal-polynomial routing,
  BNF generator display/string input, higher-degree ideal approximation, MPQS
  and full behavioral coverage remain open.
- ⬜ Not started
- 🟡 In progress (note: include agent/person working on it)
- ✅ 2026-09-11 (Codex): native FLINT word random-state dependency:
  3,036 new comparisons against compiled bundled C; live: 184,346 pass;
  caller/docs/units: 230,673 pass; fast: 6,278 pass/32 skips; all 588 type
  diagnostics unchanged; bundles pass. New kernels: 63/63 lines; combined:
  93839/117787 (79.67%). Versions 20.3.0; factor routing remains in progress.
- ✅ 2026-09-11 (Codex): variable-map domains and polynomial dictionaries:
  2,134 new comparisons; live: 181,310 pass; caller/docs/units: 227,634 pass;
  fast: 6,275 pass/32 skips; all 588 type diagnostics unchanged; bundles pass.
  Combined: 93774/117722 (79.66%); 1,240,889 configured cases. Versions 20.2.4.
  Full behavioral coverage remains open; finite-field factor routing is next.
- ✅ 2026-09-11 (Codex): matrix base arguments and parent-cache semantics:
  2,205 new comparisons; live: 179,176 pass; caller/docs/units: 225,496 pass;
  fast: 6,271 pass/32 skips; all 588 type diagnostics unchanged; bundles pass.
  Rational element: 284/284 lines; combined: 93670/117618 (79.64%), with
  1,238,755 configured cases. Versions 20.2.3; variable-map domains are next.
- ✅ 2026-09-11 (Codex): function-field parent/category predicates: 426 new
  comparisons; live: 176,971 pass; caller/docs/units: 223,289 pass; fast:
  6,269 pass/32 skips; all 588 type diagnostics unchanged; bundles pass.
  Combined execution: 93661/117609 (79.64%); 1,236,550 configured cases.
  Versions 20.2.2; matrix base arguments and broader input-domain coverage continue.
- ✅ 2026-09-11 (Codex): raw element matrices, trace/norm, divisors and evaluation:
  6,892 new comparisons expose 1,122 matrix/trace/norm discrepancies. Live:
  176,545 pass; caller/docs/units: 222,860 pass; fast: 6,266 pass/32 skips;
  all 588 type diagnostics unchanged; bundles pass. Rational element source:
  275/275 lines; combined: 93576/117530 (79.62%), with 1,236,124 configured
  cases. Versions 20.2.1; broader input-domain and category coverage continues.
- ✅ 2026-09-11 (Codex): lazy polynomial/place enumeration and caller dependencies:
  7,136 new comparisons; latest live areas: 169,653 function-field and 473
  hyperelliptic pass; caller/docs/units: 215,965 pass; fast: 6,263 pass/32 skips.
  All 588 type diagnostics unchanged; bundles pass. Cartesian iterator: 68/68
  lines; combined current-source execution: 93568/117523 (79.62%), with
  1,229,232 configured cases. Versions 20.2.0. Full behavioral coverage remains open.
- ✅ 2026-09-11 (Codex): ideal parents, powers and scalar/zero membership:
  11,220 new comparisons expose 1,816 old discrepancies. Live: 162,983 pass;
  caller/docs/units: 208,654 pass; fast: 6,259 pass/32 skips; all 588 type
  diagnostics unchanged; bundles pass. Ideal execution: 217/218 lines; divisor:
  344/344. Combined: 93429/117624 (79.43%); 1,222,096 configured cases.
  Versions 20.1.0; full behavioral coverage remains open.
- ✅ Complete (with test coverage %)
- 🔴 Blocked (note: reason)

**Maintenance notes:**
- 🟡 2026-09-11 (Codex): expanding rational-function-field parent, element, order,
  place and divisor comparisons; finite-field factor backend routing and
  dependency coverage are in progress. Higher-degree residues remain unimplemented.
- ✅ 2026-09-11 (Codex): divisor arithmetic, parent coercion, formatting and
  function-space maps: 5,640 new comparisons expose 2,842 old discrepancies. Live:
  151,763 pass; caller/docs/units: 197,431 pass; fast: 6,256 pass/32 skips; all
  588 type diagnostics unchanged; bundles pass. Divisor execution: 343/344 lines.
  Combined: 93415/117619 (79.42%); 1,210,876 configured cases. Versions 20.0.0
  for the native formatter signature; ideal and higher-degree residue coverage continue.
- ✅ 2026-09-11 (Codex): place/valuation-ring identities, residue maps and polynomial
  coefficient numerator/denominator/LCM dependencies: 6,278 new comparisons. Live:
  146,123 pass; caller/docs/units: 191,788 pass; fast: 6,251 pass/32 skips. All
  588 type diagnostics unchanged; bundles pass. New FLINT kernels: 24/24 lines;
  place_rational: 76/76; valuation_ring: 42/42. Combined: 93385/117639 (79.38%);
  1,205,236 configured cases. Versions 19.1.0; higher-degree residue fields remain open.
- ✅ 2026-09-11 (Codex): PARI squarefree component contract: 844 new comparisons
  expose 781 old discrepancies. Live area: 139,845; caller/Galois/docs replay:
  185,508 pass; fast: 6,244 pass/32 skips; all 588 type diagnostics unchanged;
  bundles pass. New kernel execution: 53/53 lines. Combined: 93238/117546 (79.32%);
  1,198,958 configured cases. Versions 19.0.0; full-factor routing remains open.
- ✅ 2026-09-11 (Codex): squarefree dispatch and GF2 roots/caching: 6,787 new
  comparisons. Live area: 139,001; caller/docs/units: 184,664; fast: 6,241 pass/32
  skips; all 588 type diagnostics unchanged; bundles pass. Native decomposition
  kernels execute 73/73 lines; GF2 101/101. Current source: 93184/117487 (79.31%),
  1,198,114 configured comparisons. Versions 18.3.0; broader factor routing remains open.
- ✅ 2026-09-11 (Codex): polynomial/fraction square roots and integer squarefree
  routing: 50,082 new native comparisons. Live area: 132,214; caller/docs/units:
  175,153; fast: 6,238 pass/32 skips; all 588 type diagnostics unchanged; bundles
  pass. Native roots retain raw fraction behavior and integer polynomial parents.
  Current-source coverage: 93143/117492 (79.28%); 1,191,327 configured comparisons.
  Versions 18.2.0; named nonsquare extensions and broader backend routing remain open.
- ✅ 2026-09-11 (Codex): FLINT word root dependencies: seven kernels now match
  69,725 direct native comparisons through unsigned/signed word boundaries.
  Live function-field/dependency area: 82,132; callers/docs/units: 125,064; fast:
  6,231 pass/32 skips; all 588 type diagnostics unchanged; bundles pass. New kernel
  execution: 96/96 lines. Combined: 92836/117829 (78.79%), 1,141,245 configured
  comparisons. Versions 18.1.0. Polynomial series and root routing are next.
- ✅ 2026-09-11 (Codex): function-field inversion and orders: 3,982 new shared
  cases cover assertion guards, zero moduli, ideal conversions, unique monoids,
  error propagation and the bundled infinite-order integer basis. Live: 12,407;
  caller/docs/units: 55,335; fast: 6,227 pass/32 skips; all 588 type diagnostics
  unchanged; bundles pass. Execution 92735/117728 (78.77%); 1,071,520 configured
  comparisons. Versions 18.0.0 for the corrected infinite-order basis entry type.
- ✅ 2026-09-11 (Codex): function-field fraction integration: 3,745 new shared
  cases cover actual parents/elements, cross-parent conversion, string fallback,
  unreduced arithmetic, predicates, display and power identity. Live: 8,425 pass;
  affected caller/docs/units: 51,348 pass; fast: 6,222 pass/32 skips; all 588 type
  diagnostics unchanged; bundles pass. Execution 92712/117730 (78.75%);
  1,067,538 configured cases. Versions 17.0.0; broader function-field audit active.
- ✅ 2026-09-11 (Codex): function-field valuation and zero-factorization batch:
  1,023 new cases cover 196 old terminating discrepancies and 188 formerly looping
  inputs. Complete fraction conversion, constant guards and bundled polynomial zero
  errors now match. Live comparisons: 4,680 pass; affected callers/docs/units: 47,599;
  fast: 6,218 pass/32 skips; type diagnostics unchanged at 588; bundles pass.
  Execution 92703/117800 (78.70%); 1,063,793 configured cases. Versions 16.0.1.
- ✅ 2026-09-11 (Codex): first rational-function-field parent/element batch: 3,648
  new comparisons expose and repair 415 discrepancies. Factory identity, renaming
  maps, backend ordering, roots and error behavior now match. Live: 3,657 pass;
  affected caller/docs/unit replay: 46,573 pass; fast: 6,215 pass/32 skips; type
  diagnostics unchanged at 588; bundles pass. Execution 92694/117797 (78.69%);
  1,062,770 configured cases. Versions 16.0.0; broader function-field audit active.
- 🟡 2026-09-11 (Codex): continuing native real audit: wide-exponent arithmetic,
  remaining source/branch coverage and native allocation edge cases.
- ✅ 2026-09-11 (Codex): direct `exp1r_abs` retains BigInt exponents and native overflow.
  Added 1,040 comparisons exposing 270 old discrepancies; the original-source oracle
  also matches the full bundled PARI build. Arithmetic/docs/units: 430,875 pass;
  fast: 6,210 pass/32 skips; live adjacent: 195 pass. Type diagnostics unchanged at
  588; bundles pass. Execution 92640/117808 (78.64%); 1,059,122 configured cases.
  Versions 15.0.0. Native invalid-allocation cases remain a documented limitation.
- ✅ 2026-09-11 (Codex): Buchmann transcendentals now share native kernels; progressive
  logarithm precision and binary-split constants restored. Added 8,061 comparisons,
  including 1,055 old result mismatches. Arithmetic/docs/units: 429,417 pass; direct
  split plus qfb/buch/units/docs: 628 pass; fast: 6,208 pass/32 skips; live adjacent:
  195 pass. All 588 type diagnostics unchanged; bundles pass. New native helper lines
  all execute; combined execution 92626/117796 (78.63%), 1,058,082 configured cases.
  Versions 14.0.0; general wide exponents and huge exponential behavior remain open.
- ✅ 2026-09-11 (Codex): Buchmann elementary real arithmetic delegates to shared
  native kernels; truncation precision guards and scaled-truncation dispatch repaired.
  Added 70,678 comparisons exposing 36,306 old mismatches. Arithmetic/docs/units:
  421,776 pass; fast: 6,203 pass/32 skips; qfb/buch plus units/docs: 208 pass; live
  adjacent: 195 pass. All 588 type diagnostics unchanged; bundles pass. Current-source
  execution: 92688/117878 (78.63%); 1,050,021 configured cases. Versions 13.0.0.
- ✅ 2026-09-11 (Codex): binary64 conversions and native real comparison repaired,
  with 16,728 direct cases exposing 5,489 old mismatches. Arithmetic/docs/conversion
  units: 342,837 pass; added comparison replay: 8,258 pass; final units/docs: 112 pass;
  fast: 6,197 pass/32 skips; qfb/buch plus units/docs: 208 pass; live adjacent: 195 pass;
  live resultant controls: 4,889 pass. All 588 type diagnostics unchanged; bundles pass.
  Current-source execution: 92716/118020 (78.56%); 979,343 configured cases. Versions
  12.0.2. Stale edited-file coverage removed; broader real arithmetic remains active.
- ✅ 2026-09-11 (Codex): native real-to-integer error metadata repaired in both
  helpers, with 4,012 comparisons exposing 798 old mismatches. Arithmetic/docs/units:
  334,366 pass; fast: 6,191 pass/32 skips; qfb/buch: 96 pass; live adjacent areas:
  195 pass. All 588 type diagnostics unchanged; bundles pass. Combined execution:
  93325/117921 (79.14%); 962,615 configured cases. Versions 12.0.1.
- ✅ 2026-09-11 (Codex): native real/complex square roots and zero-exponent halving
  restored; both low-level integer roots now use magnitude and a GMP Karatsuba
  dependency. Added 6,276 shared cases (including 66,048 hashed root/remainder
  values), exposing 874 old mismatches. Arithmetic/docs/units: 330,355 pass; fast:
  6,188 pass/32 skips; qfb/buch: 96 pass; live adjacent areas: 195 pass. All 588 type
  diagnostics unchanged; bundles pass. GMP root helper: 89/89 instrumented lines.
  Combined execution: 93333/117931 (79.14%); 958,603 configured cases. Versions 12.0.0.
- ✅ 2026-09-11 (Codex): native truncated multiplication and square kernels restored,
  with 21,466 direct comparisons exposing 3,404 old mismatches. All 110 instrumented
  multiplication-kernel lines execute. Arithmetic/docs/units: 324,078 pass; fast:
  6,183 pass/32 skips; qfb/buch: 96 pass; live adjacent areas: 195 pass. All 588 type
  diagnostics unchanged; bundles pass. Execution: 93249/117916 (79.08%); 952,327
  configured cases. Versions 11.0.3; broader source/branch audit remains open.
- ✅ 2026-09-11 (Codex): zero conversions now retain native allocation and accuracy;
  division/inverse dispatch follows allocation. Added 2,052 native comparisons
  (1,089 old mismatches) and corrected 90 hidden oracle-serialization mismatches.
  Arithmetic/docs/units: 302,611 pass; fast: 6,178 pass/32 skips; qfb/buch: 96 pass;
  adjacent replay: 195 pass. All 588 type diagnostics unchanged; bundles pass.
  Execution: 93153/117821 (79.06%); 930,861 configured cases. Versions 11.0.2.
- ✅ 2026-09-11 (Codex): native real division restored with 11,374 exact native
  comparisons (943 old mismatches in the first 9,784). Arithmetic/docs/units:
  300,558 pass; fresh real-kernel replay: 78,603 pass; fast: 6,173 pass/32 skips;
  qfb/buch: 96 pass; live quadratic forms: 82 pass. All 588 type diagnostics
  unchanged; bundles pass. Combined execution: 93151/117838 (79.05%); 928,809
  configured cases. Versions 11.0.1; remaining branch/source audit stays open.
- ✅ 2026-09-11 (Codex): free PARI factorial now returns the native rounded real with
  bigint exponent; exact/medium/gamma backends and dependency algorithms restored.
  Added 10,085 comparisons, exposing 1,012 old selector mismatches in 1,032 safe calls.
  Arithmetic/docs/units: 289,183 pass; final dependency replay: 10,089 pass; fast:
  6,168 pass/32 skips; qfb/buch: 96 pass; live quadratic: 82; adjacent replay: 113.
  All 588 type diagnostics unchanged; bundles pass. Combined execution: 92966/117649
  (79.02%); 917,435 configured cases. Versions 11.0.0; broader audit remains open.
- ✅ 2026-09-11 (Codex): PARI real/mixed-integer addition and subtraction repaired
  with 57,825 exact native comparisons (23,993 old mismatches in the first 57,820).
  All instrumented addition-kernel lines execute. Arithmetic/docs: 279,088 pass plus
  five boundary cases; fast: 6,163 pass/32 skips; qfb/buch: 96 pass; live quadratic
  forms: 82 pass. All 588 type diagnostics unchanged; bundles pass. Combined execution:
  92466/117153 (78.93%); 907,350 configured cases. Versions 10.0.2.
- ✅ 2026-09-11 (Codex): Integer and free GMP/default factorial routing, exact sign/word
  errors and selector validation restored. Added 3,098 comparisons, exposing 197 old
  mismatches in 984 safe calls. Native kernels: 166/168 lines (unhit closing braces only).
  Current arithmetic replay/docs/units: 221,609 pass; all 6159 fast cases covered across
  the initial run, unchanged timeout retries and the added public example; 32 skips.
  All 588 type diagnostics unchanged; bundles pass. All-area replay: 846,337 pass;
  the 90 added valuation and 3,098 factorial cases also pass, totaling 849,525.
  Combined execution: 92398/116208 (79.51%). Versions 10.0.1.
- ✅ 2026-09-10 (Codex): valuation dispatch/infinity, native GMP factor extraction and
  Integer metadata repaired with 32,307 comparisons exposing 3,364 old mismatches.
  Affected live: 231,409 pass; final replay/docs/quaternion: 218,211 pass; fast:
  6,157 pass/32 skips. Slow: all 1,274 pass across full run and two unchanged timeout
  retries. All 588 type diagnostics unchanged; bundles pass. Fresh execution:
  92230/116628 (79.08%); 846,427 configured cases. Versions 10.0.0.
- ✅ 2026-09-10 (Codex): complex display ordering repaired for exact-zero classification,
  binary rounding, tuple metadata, identity and nonfinite values. Added 6,368 comparisons:
  3,120 display cases (2,969 old mismatches) and 3,248 native Python sort result/trace cases.
  Live: 185,866 pass; replay/docs: 185,861 pass; fast: 6,155 pass/32 skips; all 588 type
  diagnostics unchanged; bundle passes. Sorting adapter: 258/266 lines (eight unhit brace
  lines); combined execution 92236/116124 (79.43%). Versions 9.1.0.
- ✅ 2026-09-10 (Codex): maximal-quotient reconstruction now preserves exact source
  division and IntegerLike zero/threshold behavior. Added 44,332 comparisons; 18,629
  old mismatches, including 4,973 in the stated domain. Live: 179,498 pass; replay/docs:
  179,492 pass; fast: 6,150 pass/32 skips; all 588 type diagnostics unchanged; bundle passes.
  Fresh execution: 91979/115882 (79.37%); workspace/lock versions 9.0.1.
- ✅ 2026-09-10 (Codex): exact square-decomposition errors, wrapped inputs, count dispatch
  and mirrored uint32 methods repaired. Added 11,616 original comparisons; 5,415 old
  mismatches in 9,308 safe calls. All eight supported entry points compared; native
  sum_of_squares execution 107/107 lines. Live: 142,100 pass; replay/docs: 135,159 pass;
  fast: 6,150 pass/32 skips; all 588 type diagnostics unchanged; bundles pass.
  Fresh execution: 91986/115902 (79.37%). Versions 9.0.0 for the null-to-error change.
- ✅ 2026-09-10 (Codex): native totient/divisor delegation and shared valuation extraction
  restored. Added 2,490 native comparisons; affected live: 130,566 pass; replay/docs:
  123,542 pass. Full saved integration: 751,804 pass; fast: 6,148 pass/32 skips;
  slow: 1,274 pass. All 588 type diagnostics unchanged; bundles pass. Fresh execution:
  91906/116300 (79.02%); expanded slow instrumentation increases the measured denominator.
  Workspace/lock versions 8.18.0.
- ✅ 2026-09-10 (Codex): signed/zero extended LCM and Dedekind psi behavior repaired;
  IntegerLike scalar signatures normalized. Added 19,624 comparisons; 13,172 old
  mismatches (the zero-radical 10.3 difference uses the existing bundled-source adapter).
  Affected live: 128,226 pass; replay/docs: 121,051 pass; fast: 6,145 pass/32 skips;
  all 588 type diagnostics unchanged; bundle passes. Fresh execution: 91900/115881
  (79.31%). Workspace/lock versions 8.17.0.
- ✅ 2026-09-10 (Codex): integer ProductTree and streaming product-rule support added;
  smooth/coprime extraction now follows original remainder trees, sorting, repeated
  factors and signed floor quotients. Added 12,956 original cases (300 old mismatches
  in 2,332 safe existing-input probes). Live: 101,333 pass; replay/docs: 101,426 pass;
  fast: 6,144 pass/32 skips; all 588 type diagnostics unchanged; bundles pass.
  Fresh execution: 91887/115925 (79.26%). Versions 8.16.0. General-ring domains remain open.
- ✅ 2026-09-10 (Codex): CRT list/modular dispatch, balanced combination, signed basis
  residues, source fallback prefixes and ragged-vector errors repaired. Added 11,278
  bundled-original comparisons; 3,610 old mismatches. Affected live: 88,815 pass;
  replay/docs: 88,469 pass; fast: 6,143 pass/32 skips. All 588 type diagnostics unchanged;
  bundles pass. Fresh execution: 91784/115833 (79.24%). Versions 8.15.0.
- ✅ 2026-09-10 (Codex): native-width gcd/inverse factory dispatch implemented with
  original conversion bounds, signed remainder, errors and bound-method identity.
  Added 13,704 original comparisons; live arithmetic: 77,099 pass; replay/docs: 77,190
  pass; fast: 6,142 pass/32 skips. All 588 type diagnostics unchanged; bundles pass.
  Fresh execution: 91750/115825 (79.21%). Workspace/lock versions 8.14.0.
- ✅ 2026-09-10 (Codex): continuant/zero scalar edges, wrapped returns and lazy squarefree
  divisor masks repaired; 8,468 new original comparisons. Affected live: 70,353 pass;
  replay/docs: 63,485 pass; fast: 6,141 pass/32 skips. All 588 existing type diagnostics
  unchanged; arithmetic bundle passes. Fresh execution: 91659/115738 (79.20%).
  Workspace/lock versions synchronized at 8.13.0.
- ✅ 2026-09-10 (Codex): all-area saved integration replay passed 660,560 comparisons.
  Hilbert rational/algorithm dispatch, wrapped-input loops and native valuation parity
  repaired with another 22,724 original cases. Affected live: 55,017 pass; replay/docs/native
  units: 55,019 pass; fast: 6,139 pass/32 skips. All 588 existing type diagnostics unchanged;
  arithmetic/PARI bundles pass. Fresh execution: 91613/115794 (79.12%). Versions 8.12.0.
- ✅ 2026-09-10 (Codex): prime-power dependency routing, wrapped pseudoprime powers,
  predecessor error dispatch and native successor/predecessor wheel repaired; 9,870 new
  original cases. Affected live: 39,240 pass; replay/docs: 32,291 pass; fast: 6,135 pass/32
  skips. All 588 existing type diagnostics unchanged; two entry points bundle. Fresh
  execution: 91557/115797 (79.07%). Workspace/lock versions synchronized at 8.11.0.
- ✅ 2026-09-10 (Codex): algebraic-dependency parent/options dispatch and native MPFR
  addition/multiplication repaired; 33,642 new original cases. Affected live: 135,488 pass;
  replay/docs: 134,760 pass; fast: 6,132 pass/32 skips. All 588 existing type diagnostics
  unchanged; three entry points bundle. Fresh execution: 91635/116114 (78.92%). Versions
  8.10.0. Native LLL/interval-boundary and complex/p-adic domains remain audit gaps.
- ✅ 2026-09-10 (Codex): integer floor/ceiling/truncation method and float fallback
  dispatch repaired; 1,569 new original cases. Affected live area: 11,947 pass;
  replay/docs: 12,032 pass; fast: 6,127 pass/32 skips. All 588 existing type diagnostics
  unchanged; arithmetic bundle passes. Fresh execution: 91441/116024 (78.81%).
  Workspace/lock versions synchronized at 8.9.0.
- ✅ 2026-09-10 (Codex): Dedekind backend selection, signed moduli, native word/generic
  conventions and zero errors repaired; 8,298 new original cases. Affected live area:
  10,378 pass; replay/docs: 10,462 pass; fast: 6,126 pass/32 skips. All 588 existing
  type diagnostics unchanged; arithmetic/PARI/FLINT bundles pass. Fresh execution:
  91409/116014 (78.79%). Versions 8.8.0; native FLINT quotient-batching optimization remains open.
- ✅ 2026-09-10 (Codex): PARI prime ordinal lookup/delegation, Integer coercion and
  fractional/nonfinite first-prime counts repaired; 795 new original cases. Live affected
  comparisons: 9,117 pass; replay/docs: 2,163 pass; fast: 6,067 pass/32 skips.
  All 588 existing type diagnostics unchanged; arithmetic/PARI bundles pass. Fresh execution:
  91376/116012 (78.76%). Workspace/lock versions synchronized at 8.7.0.
- ✅ 2026-09-10 (Codex): prime iterator bounds/laziness, differences exhaustion and
  subfactorial validation repaired; 1,228 new original cases. Live arithmetic/integer
  comparisons: 8,322 pass; replay/docs: 1,368 pass; fast: 6,067 pass/32 skips.
  All 588 existing type diagnostics unchanged; arithmetic bundle passes. Fresh execution:
  91273/115925 (78.73%). Workspace/lock versions synchronized at 8.6.0.
- ✅ 2026-09-10 (Codex): Gaussian copy overrides and public numeric option types repaired;
  432 new original comparisons. Affected live areas: 5,863 pass; replay/docs: 5,930 pass;
  fast: 6,066 pass/32 skips; Gaussian units: 157 pass. Type diagnostics fall from 589 to
  588 by correcting the mixed numeric/integer options type; no new diagnostics. Stats
  bundle passes. Fresh execution: 91,197/116,306 (78.41%). Versions synchronized at 8.5.0.
- ✅ 2026-09-10 (Codex): native fractional parts, value copying and common-parent real
  comparisons repaired; 31,500 new original comparisons. Affected live areas: 95,330 pass.
  Replay/docs: 89,165 pass; fast: 6,065 pass/32 skips; Gaussian tests: 157 pass. All 589
  existing type diagnostics unchanged; rings/MPFR bundles pass. Fresh execution:
  91,198/116,307 (78.41%). Workspace/lock versions synchronized at 8.4.0.
- ✅ 2026-09-10 (Codex): native real predicates, integer rounding and MPFR primitives
  repaired; 23,484 new original comparisons. Affected original records: 63,830 pass
  (62,246 live plus 1,584 native extreme cases). Replay/docs: 57,664 pass; fast: 6,064
  pass/32 skips; Gaussian tests: 157 pass. All 589 existing type diagnostics unchanged;
  rings/MPFR bundles pass. Fresh execution: 91,073/116,172 (78.39%). Version 8.3.0.
- ✅ 2026-09-10 (Codex): nearest-even real numeric/string constructors and native MPFR
  integer setter repaired; 6,267 new original comparisons. Affected live areas: 40,346 pass.
  Replay/docs: 34,179 pass; fast: 6,063 pass/32 skips; Gaussian tests: 157 pass. All 589
  existing type diagnostics unchanged; rings/MPFR bundles pass. Fresh execution:
  90,902/115,992 (78.37%). Versions synchronized at 8.2.0. Wider rerun confirms all
  543,475 preceding cases; 6,267 cases added during that run pass their separate live run.
- ✅ 2026-09-10 (Codex): integer Gaussian validation precedence, fractional-tail errors,
  large-center rounding and exact integer bounds repaired; 806 new original comparisons.
  Affected live areas: 5,431 pass. Replay/docs: 5,494 pass; fast: 6,062 pass/32 skips;
  affected samplers: 157 pass. All 589 existing type diagnostics unchanged; stats bundle
  passes. Fresh execution: 90,866/115,956 (78.36%). Versions synchronized at 8.1.1.
- ✅ 2026-09-10 (Codex): Gaussian parameter repr and native MPFR double conversion repaired;
  added 8,265 original comparisons. Affected four live areas: 33,273 pass; the other 25
  retain the preceding full-suite baseline. Replay/docs: 32,520 pass; fast: 6,061 pass/32
  skips; affected sampler tests: 157 pass. All 589 existing type diagnostics unchanged;
  stats/MPFR bundles pass. Fresh execution: 90,848/115,938 (78.36%). Direct public dependency
  import is executed; workspace and lock versions synchronized at 8.1.0.
- ✅ 2026-09-10 (Codex): lattice Gaussian defaults, deferred centers, cache transitions,
  dependent bases and delayed sampling errors repaired; fixed-six sigma errors match Python.
  Added 2,913 original comparisons. All 534,404 cases in 29 areas pass with strict exception
  messages; 159 lattice cases also pass live with the bundled-message adapter. Replay/docs:
  3,301 pass; fast: 6,060 pass/32 skips; affected slow: 946 pass. All 589 existing type
  diagnostics unchanged; stats bundle passes. Fresh execution: 90,767/115,859 (78.34%).
  Nullable center accessors and the empty low-level return are documented in version 8.0.0.
- ✅ 2026-09-10 (Codex): binary density coercion and PNG errors repaired; original Gamma
  predicates verified. Added 10,961 comparisons; all 51 implemented binary API callables
  have comparative dispatch. Matrix: 139,368 saved-original pass; polynomial/extended:
  277,296 pass. Replay/docs: 145,244 pass; fast: 6,058 pass/32 skips; affected slow: 946 pass.
  All 589 existing type diagnostics unchanged; matrix bundle passes. Fresh execution:
  90,738/115,930 (78.27%). Versions synchronized at 7.2.0.
- ✅ 2026-09-10 (Codex): binary seeded words, sparse/nonzero random schedules and factory
  density/dimension/shape dispatch repaired; 10,032 original comparisons. Matrix: 128,407
  saved-original pass; polynomial/extended matrix: 277,296 pass. Replay/docs: 134,282 pass.
  Fast: 6,057 pass/32 skips; affected slow: 946 pass. All 589 existing type diagnostics
  unchanged; matrix bundle passes. Fresh execution: 89,919/115,107 (78.12%). Released as 7.1.0.
- ✅ 2026-09-10 (Codex): binary constructor shapes/dimensions, field entry coercion and
  cached-row ownership repaired; 11,554 fresh original comparisons. Matrix: 118,375
  saved-original pass; polynomial/extended matrix: 277,296 pass. Replay/docs: 124,249 pass.
  Fast: 6,056 pass/32 skips; affected slow: 946 pass. All 589 existing type diagnostics
  unchanged; matrix bundle passes. Fresh execution: 89,905/115,100 (78.11%). Cached row
  assignment now raises the original immutable-vector error; released as 7.0.0.
- ✅ 2026-09-10 (Codex): binary formatting, subdivisions, coercion/order and metadata
  propagation repaired; 23,178 new original comparisons. Matrix: 106,821 saved-original
  pass, including all new observed cases; polynomial/extended matrix: 277,296 pass.
  Replay/docs: 112,694 pass; fast: 6,055 pass/32 skips. All 589 existing type diagnostics
  unchanged; matrix bundle passes. Fresh execution: 89,769/114,905 (78.12%). Large
  transcript serialization now uses bounded ASCII-safe chunks. Released as 6.1.0.
- ✅ 2026-09-10 (Codex): binary basic-operation errors, native transpose/add/concatenation
  delegation and column cache/ownership repaired; 10,773 new original comparisons. Matrix:
  83,643 live pass; polynomial/extended matrix: 277,296 saved-original pass; replay/docs:
  89,515 pass. Fast: 6,053 pass/32 skips; affected slow: 946 pass. All 589 existing type
  diagnostics unchanged; bundles pass. Fresh execution: 89,565/114,763 (78.04%). Column
  vectors are now immutable with original assignment errors; released as 6.0.0.
- ✅ 2026-09-10 (Codex): binary rectangular/singular solves, native solve/kernel delegation,
  basis/options/proof validation and cache/mutability behavior repaired; 12,222 new original
  comparisons. Matrix: 72,870 live pass; polynomial/extended matrix: 277,296 saved-original
  pass; replay/docs: 78,741 pass. Fast: 6,052 pass/32 skips; affected slow: 946 pass. All 589
  existing type diagnostics unchanged; bundles pass. Fresh execution: 89,485/114,712 (78.01%).
- ✅ 2026-09-10 (Codex): binary swap index conversion, permutation degrees/cycle order
  and native inverse delegation repaired; 18,153 new original comparisons. Matrix:
  60,648 live pass; polynomial/extended matrix: 277,296 saved-original pass; replay/docs:
  66,518 pass. Fast: 6,051 pass/32 skips; affected slow: 946 pass. All 589 existing type
  diagnostics unchanged; bundles pass. Fresh execution: 89,328/114,576 (77.96%).
- ✅ 2026-09-10 (Codex): standalone PLE/PLUQ delegation, parameter conversion and native
  seven-table PLE schedule repaired; 3,588 new original comparisons. Matrix: 42,495 live
  pass; polynomial/extended matrix: 277,296 saved-original pass; replay/docs: 48,364 pass.
  Fast: 6,050 pass/32 skips; affected slow: 946 pass. Type diagnostics reduced from 605
  to 589 with no additions; bundles pass. Fresh execution: 89,334/114,594 (77.96%).
- ✅ 2026-09-10 (Codex): binary elimination dispatch, native M4RI PLE/PLUQ and triangular
  algorithms, immutable cached forms, keyword errors and mutation cache behavior repaired.
  Added 12,924 original comparisons; matrix area: 38,907 live pass; full polynomials:
  276,371 saved-original pass. Replay/docs: 62,497 pass; fast: 6,049 pass/32 skips; affected
  slow: 946 pass. All 605 existing type diagnostics unchanged; bundles pass. Fresh line
  execution: 89,327/114,636 (77.92%). Standalone factorization delegation remains open.
- ✅ 2026-09-10 (Codex): exact/approximate binary density and bounded native real-literal
  conversion repaired, with 24,161 new original comparisons. Affected live areas: 47,751
  pass; full polynomials: 276,371 saved-original pass. Replay/docs: 49,572 pass; fast: 6,046
  pass/32 skips; affected slow: 946 pass. All 605 existing type diagnostics unchanged;
  affected bundles pass. Fresh line execution: 88,709/114,041 (77.79%). Documented native
  conversion limits and existing real arithmetic gaps remain explicit. Version 4.0.0.
- ✅ 2026-09-10 (Codex): binary matrix indices, rational conversion, bounds and slices
  repaired, with 19,522 new original comparisons and native M4RI extraction delegation.
  Matrix area: 22,775 live pass; full polynomials: 276,371 saved-original pass.
  New/previous binary replay and 65 docs: 25,410 pass; fast: 6,042 pass/32 skips;
  affected slow: 946 pass. All 605 existing type diagnostics remain unchanged; matrix/M4RI
  bundles pass. Fresh execution coverage: 88,334/113,755 instrumented lines (77.65%).
- ✅ 2026-09-10 (Codex): specialized GF(2) polynomial matrix evaluation, mixed/list return
  types and M4RI multiplication delegation repaired; 5,823 new original comparisons.
  All 374,423 configured comparisons pass (279,624 live; 94,799 saved originals).
  New replay/docs: 5,887 pass; fast: 6,041 pass
  with 32 skips; affected slow: 946 pass. Full polynomial replay: 276,326 pass including
  one targeted timeout retry; native GCD replay: 3,214 pass. Type diagnostics remain 605
  unchanged; matrix/polynomial/M4RI bundles pass. Fresh recorded source lines:
  88,227/113,670 (77.62%). Other binary methods remain under source review.
- ✅ 2026-09-10 (Codex): generic/integer/modular matrix evaluation actions and native
  integer-to-real conversion repaired, with 21,017 new original comparisons. Live polynomial:
  273,779 pass; live native: 43,505 pass; other 26 areas: 51,316 pass. New replay/docs: 21,080
  pass; full polynomial replay: 273,734 pass; native GCD replay: 3,214 pass. Fast: 6,036 pass/32
  skips; affected slow: 946 pass. All 605 existing type diagnostics remain unchanged. Matrix
  and polynomial entry points bundle successfully. Fresh coverage: 87,939/113,365 instrumented
  source lines (77.57%); both new production files have 100% recorded line execution.
- ✅ 2026-09-10 (Codex): polynomial scalar/composition dispatch, sparse compiled evaluation,
  generator flags/cache, native FLINT/NTL evaluation/composition/Taylor kernels repaired
  in the compared domain. Added 34,922 original comparisons. Live polynomial: 252,971 pass;
  native: 43,296 pass; other 26 areas: 51,316 pass. New replay/docs: 34,984 pass.
  Fast: 6,032 pass/32 skips; affected slow: 946 pass. Type diagnostics: 605, no additions
  and four previously invalid call sites resolved. Full polynomial replay: 252,926 pass.
  Coverage: 87,704/112,938 instrumented lines (77.66%); broader audit remains open.
- ✅ 2026-09-10 (Codex): modular polynomial powers, RDF/GSL integer powers, real division
  and bounded modular GCD/Newton products repaired. Added 16,736 original comparisons.
  Live polynomial: 226,648 pass; live native: 34,697 pass; other 26 areas: 51,316 pass.
  New replay/docs: 16,796 pass; full polynomial replay: 226,603 pass; GCD replay: 3,214 pass.
  Fast: 6,030 pass/32 skips; affected slow: 946 pass; type diagnostics unchanged at 609.
  Coverage: 87,005/112,335 instrumented lines (77.45%). Broader audit remains open.
- ✅ 2026-09-10 (Codex): full exact polynomial products now use FLINT/NTL dispatch and
  generic Karatsuba. Added 1,544 original comparisons, including noncommutative coefficients.
  Live polynomial: 211,714 pass; live native: 32,895 pass; other 26 areas: 51,316 pass.
  New replay/docs: 1,603 pass; full polynomial replay: 211,669 pass. Fast: 6,029 pass/32
  skips; affected slow: 946 pass; type diagnostics unchanged at 609. Coverage:
  86,395/111,524 lines (77.47%). Broader source/behavioral coverage remains open.
- ✅ 2026-09-10 (Codex): truncated scalar/nested-parent coercion, nullable operands,
  ordered real products and formatting repaired in the compared domain. Added 14,012
  original comparisons. Live polynomial: 210,690 pass; other 26 Sage areas: 51,316 pass;
  new replay/docs: 14,070 pass. Fast: 6,028 pass/32 skips; affected slow: 946 pass;
  type diagnostics unchanged at 609. Coverage: 86,224/111,363 lines (77.43%).
  Exact/native public multiplication and broader source/behavioral coverage remain open.
- ✅ 2026-09-10 (Codex): bounded native integer division and modular KS/KS2/KS4
  products repaired. Added 38 original regressions, 28 failing before the fix. Live native
  comparisons: 13,424 pass; other 27 areas: 247,994 pass. Fast: 6,027 pass/32 skips;
  affected slow: 946 pass; type diagnostics unchanged at 609. Both edited kernels have
  100% recorded line execution. Expanded merged instrumentation: 86,099/111,402 (77.29%).
  Mixed-parent truncation and broader source/behavioral coverage remain open.
- ✅ 2026-09-10 (Codex): rational polynomial powers, roots, truncated-series dispatch,
  finite coefficient root selection, numeric factor order and p-adic prime validation
  repaired in the compared domain. Added 39,926 original comparisons; all 280,331
  configured comparisons pass across live affected and saved remaining areas. Fast:
  6,027 pass/32 skips; affected slow: 946 pass; new replay/docs: 39,983 pass.
  Existing type diagnostics remain 609. Current-source execution coverage:
  86,100/110,707 lines (77.77%). Mixed-parent truncation, modular powers and broader
  behavioral/source coverage remain open.
- ✅ 2026-09-10 (Codex): integer polynomial powers, negative fraction results, backend
  dispatch, coefficient preservation and bounded native products repaired in the compared
  domain. Added 7,856 original comparisons; live areas: 163,780 Sage / 25,337 native pass;
  full saved-original suite: 240,405 pass. Fast: 6,023 pass/32 existing skips plus two
  final runtime/type tests; affected slow: 946 pass. Type diagnostics unchanged at 609.
  Current-source execution coverage: 85,218/110,190 lines (77.34%). Fractional/modular
  powers, root-series dependencies and broader source coverage remain open.
- ✅ 2026-09-10 (Codex): direct fraction/polynomial scalar conversions and hooks,
  GF2 object constructors, unit inverse errors and native zero allocation repaired in the
  compared domain. Added 7,338 comparisons; live polynomial area: 160,456 pass; full saved-original
  comparison: 232,549 pass. Fast:
  6,020 pass/32 existing skips; affected slow: 946 pass; replay/focused/docs: 7,399 pass.
  Type diagnostics remain 609 (one equivalent union reordered). Current-source execution
  coverage: 84,681/109,645 lines (77.23%). Public polynomial power remains open.
- ✅ 2026-09-10 (Codex): fraction-to-polynomial sections, mixed fraction construction,
  isolated imports, canonical parent queries and nested coefficient closure repaired in
  the compared domain. Added 16,937 comparisons; full live suite: 225,211 pass. Fast:
  6,013 pass/32 existing skips; affected slow: 946 pass; final replay/focused/docs: 16,994
  pass. Type diagnostics unchanged at 609. Current-source execution coverage:
  84,626/109,685 lines (77.15%). Direct scalar hooks and public polynomial power remain open.
- ✅ 2026-09-10 (Codex): generic/field/FpT polynomial fractions, nonconstant denominators,
  reduction flags, inversion/powers, parent/accessor identity and portable native modular
  products repaired in the compared domain. Added 16,017 comparative regressions;
  current 208,274 distinct cases pass. Fast: 6,013 pass/32 existing skips; affected slow:
  946 pass; new-case replays: 16,017 pass. Primary type diagnostics remain 609.
  Current-source execution coverage: 84,459/109,565 lines (77.09%). Complete fraction
  coercion, coefficient power closure and public polynomial multiplication/power remain open.
- ✅ 2026-09-10 (Codex): pseudo-division degree gaps, scalar protocols, mixed parents,
  fraction-field zero quotients and zero-ring scalar equality repaired in the compared
  domain. Added 15,994 comparisons; full live suite: 192,257 pass; all new-case replays
  pass. Fast units: 6,007 pass/32 existing skips; affected slow: 946 pass. Primary type
  diagnostics remain 609. Current-source execution coverage: 83,914/108,987 lines (76.99%).
  Nonconstant fraction denominators and full native polynomial power dispatch remain open.
- ✅ 2026-09-10 (Codex): derivative variables/counts, aliases, gradient, nested coefficients,
  zero identity and RDF scalar rounding repaired; QQ delegates to FLINT. Added 21,025
  comparisons; full live suite: 176,263 pass; all new-case replays pass. Fast units: 6,002
  pass/32 existing skips; affected slow: 946 pass. Primary type diagnostics remain 609.
  Current-source execution coverage: 83,622/108,761 lines (76.89%).
- ✅ 2026-09-10 (Codex): general RDF resultants use native PARI single-word elimination;
  corrected inexact leading-zero handling at the dependency boundary. Added 7,062 original
  comparisons; full live suite: 155,238 pass; all new-case replays pass. Fast units: 5,993
  pass/32 existing skips; primary type diagnostics unchanged at 609. Current-source
  execution coverage: 83,530/108,675 lines (76.86%). SciPy overflow, multiword real and
  native large-prime/binary/extension resultant backends remain open.
- ✅ 2026-09-10 (Codex): rational FLINT resultant and Sage ZZ/QQ/word-prime delegation,
  composite matrix behavior, mixed/scalar parents and proof keywords repaired. Added RDF
  numeric coefficients and PARI constant-operand real resultants. Full live suite: 148,176
  pass. All 23,968 new comparisons
  and replays pass; polynomial area 88,120 pass/native backend 8,768 pass. Fast units:
  5,988 pass/32 existing skips; affected slow: 946 pass. Primary type diagnostics remain
  609. Current-source execution coverage: 83,446/108,482 lines (76.92%). General real and
  native large-prime/binary/extension resultants remain open.
- ✅ 2026-09-10 (Codex): Sylvester mixed/scalar parents, zero errors, optional variable
  coercion and native finite matrix errors repaired, with direct scalar-constructor fixes.
  Added 18,353 comparisons; final polynomial area: 66,105 pass; full live suite: 124,144
  pass before 64 final degree-one cases. All 18,353 new-case replays pass. Fast units:
  5,978 pass/32 existing skips (final run also included 18,289 replays); affected slow:
  946 pass. Primary type diagnostics remain 609. Current-source execution coverage:
  83,171/108,265 lines (76.82%).
- ✅ 2026-09-10 (Codex): integer/rational XGCD delegated to FLINT ports; Integer constant
  returns and denominator-cleared noncoprime fallback implemented. All 1,791 new comparisons
  and replays pass, including native coefficient/sign and allocation checks. Fast: 5,974
  pass/32 existing skips; affected slow: 946 pass; primary type diagnostics unchanged at 609.
  Full saved-Sage/native suite: 105,855 pass. Current-source execution coverage:
  83,064/108,133 lines (76.82%).
- ✅ 2026-09-09 (Codex): modular-integer/ring repair batch: 34/34 public members dispatched,
  3,517 modular comparisons pass. Full suite: 36,929 pass before 36 final string controls;
  fast units: 5,895 pass, 32 existing skips; affected slow units: 946 pass; final focused
  units/docs: 112 pass. Measured line coverage: 81,166/105,683 (76.80%). No new primary
  type-check diagnostic kinds; widened overloads repeat 18 existing RingElement diagnostics,
  while 149 previous diagnostic instances disappear through the broader coercion signatures.
- ✅ 2026-09-09 (Codex): modular factories now handle default/zero/negative orders, signed
  weak caching and Mod's original-value/explicit-parent paths, including both prime classes.
  Added 373 comparisons; modular area: 3,890 pass; full suite: 37,338 pass. Final fast tests:
  5,902 pass, 32 existing skips; affected slow tests: 946 pass; focused: 454 pass; corrected
  legacy factory assertions/docs: 102 pass. Type-check diagnostics unchanged from 1.2.0.
  Instrumented lines: 81,217/105,513 (76.97%). Version 2.0.0 marks the necessary factory return
  type widening, not audit completion.
- ✅ 2026-09-09 (Codex): prime/modular powers now share native/GMP dispatch and preserve
  the original exponent type before coercion. Added 897 comparisons: field area 15,576 pass,
  modular area 4,115 pass, full suite 38,235 pass. Fast suite: 5,902 pass, 32 existing skips
  before the final colocated regression; final focused units/docs: 119 pass. All 946 affected
  slow tests pass. Type-check diagnostic instances unchanged. Instrumented lines:
  81,201/105,665 (76.85%).
- ✅ 2026-09-09 (Codex): mixed finite-ring coercion and Python sequence multiplication:
  7,979 new comparisons, 46,214 full-suite cases pass (23,555 field, 4,115 modular). Corrected
  the old sequence-result oracle and exposed promoted/sequence return overloads. Fast tests:
  5,907 pass, 32 existing skips. Slow tests: 945 pass, one timeout; both matching cases pass
  in isolation, including the timed-out test in 2.76 s. Final focused/replay tests: 252 pass.
  Type checking remains failing: 763 diagnostic instances versus 1,397 before this batch;
  overloads alter some existing inference/recursive-interface diagnostics. Instrumented lines:
  81,148/105,804 (76.70%).
- ✅ 2026-09-09 (Codex): finite/QQ generator indices, prime parent integer coercion and
  modern prime/extension _integer_ hooks: 1,209 new cases, 101/101 field members dispatched.
  Full suite: 47,380 pass before 43 final QQ-index cases; final field/rational areas: 31,066
  pass. Fast: 5,909 pass, 32 existing skips before final QQ/API additions; focused: 269 pass;
  transcript replay: 1,613 cases pass. All 946 affected slow tests pass (15-second per-test
  timeout). Type-check diagnostic instances unchanged at 763. Instrumented lines:
  81,193/106,081 (76.54%).
- ✅ 2026-09-09 (Codex): extension arithmetic now delegates to PARI-port FpX/FpXQ kernels;
  signed quotient powers, inverse signatures and exponent validation repaired. Added 3,940
  comparisons; all 51,363 full-suite cases pass, including 28,704 field cases. Fast tests:
  5,914 pass with 32 skips before the final API example. All 946 affected Sage slow tests and
  62 direct Galois tests pass; focused source: 219 pass; API examples: 35 pass; final new-case
  replay: 3,940 pass. Type-check diagnostics unchanged at 763. Instrumented lines:
  81,254/106,262 (76.47%). Native polynomial packing and fast multiplication remain open.
- ✅ 2026-09-09 (Codex): polynomial constructor/index and coefficient conversion batch:
  1,626 new comparisons; all 52,989 full-suite cases pass, including 1,671 polynomial cases.
  Exact-parent identity, polynomial-base constants, backend list rules, rational hooks and
  compatible extension inputs repaired. Fast: 5,920 pass, 32 skips before final API addition;
  all 946 affected slow tests pass; focused/replay: 1,920 pass; final docs/constructor: 41 pass.
  Type-check diagnostics decrease from 763 to 609. Current measured lines: 80,543/106,155
  (75.87%); finite-field member dispatch 103/103, modular dispatch 35/35.
- ✅ 2026-09-09 (Codex): integer-polynomial GCD/evaluation comparisons now call production
  methods and retain GCD content. Added 216 cases; all 53,205 full-suite cases, 246 area
  cases and all 216 replayed cases pass. Fixed 12 zero-GCD sign/identity mismatches and corrected the old unit assertion.
  Fast tests: 5,921 pass, 32 skips; focused: 46 pass; type diagnostics unchanged at 609.
  Current measured lines: 80,781/106,358 (75.95%).
- ✅ 2026-09-09 (Codex): FLINT integer/modular GCD kernels and Sage delegation repaired;
  added 896 comparisons. Full run: 54,060 pass; final area: 1,142 pass; new replay: 896 pass.
  Fast tests: 5,926 pass, 32 existing skips; affected slow tests: 946 pass; docs/backend: 43 pass.
  Type diagnostics unchanged at 609. Current measured lines: 81,857/106,906 (76.57%).
- ✅ 2026-09-09 (Codex): monic parent/coefficient/error/identity repair and specialized GF2
  division/inverse/power errors; 573 new comparisons. Final area: 2,244 pass; full saved-Sage
  comparison: 54,674 pass; new replay: 573 pass. Fast: 5,931 pass/32 skips; slow: 946 pass;
  docs/focused: 44 pass. Type diagnostics remain 609. Measured lines: 81,883/106,919 (76.58%).
- ✅ 2026-09-09 (Codex): polynomial indexing/shift/truncation/reversal repair; 9,608 new
  comparisons, all 64,282 live Sage comparisons pass. All 9,608 replay cases pass. Fast
  tests: 5,935 pass/32 skips plus one new-doc missing-import failure, repaired and verified
  in the final 43 focused tests. All 946 affected slow tests pass. Type diagnostics remain
  609. Measured lines: 82,104/107,132 (76.64%). Allocation limits remain documented.
- ✅ 2026-09-09 (Codex): common-parent polynomial add/sub/mul/equality and representation
  repaired; 7,554 new comparisons, all 71,836 live Sage comparisons pass. New replay: 7,554
  pass; final fast: 5,941 pass/32 skips; final affected slow: 946 pass; focused: 488 pass;
  API examples: 40 pass. Type diagnostics remain 609. Measured lines: 82,200/107,243 (76.65%).
- ✅ 2026-09-10 (Codex): native integer/modular division kernels and Sage quotient/remainder
  repair; 7,002 new comparisons, all 78,838 live Sage comparisons pass. New replay: 7,002
  pass; fast: 5,946 pass/32 skips before the last API example; affected slow: 946 pass;
  focused: 36 pass; final API examples: 41 pass. Type diagnostics remain 609. Measured
  lines: 82,283/107,393 (76.62%). Native QQ/large-prime/extension fast division remains open.
- ✅ 2026-09-10 (Codex): polynomial scalar equality and coefficient-polynomial embedding
  repaired; 9,790 new comparisons and replays pass. All 88,628 live Sage comparisons,
  5,950 fast tests (32 existing skips) and 946 affected slow tests pass. Type diagnostics
  remain 609. Current measured lines: 82,345/107,461 (76.63%).
- ✅ 2026-09-10 (Codex): polynomial GCD coercion/native shortcuts and FLINT rational,
  word-modular and NTL binary delegation repaired. Added 5,955 comparisons/replays;
  all 94,583 live comparisons pass. Fast: 5,956 pass/32 skips; affected slow: 946 pass;
  primary type diagnostics unchanged at 609. Measured lines: 82,457/107,550 (76.67%).
  Nested Singular GCD and fast large-prime/extension NTL GCD remain open.
- ✅ 2026-09-10 (Codex): finite extended GCD backend/coercion/zero/alias repair:
  7,884 new comparisons/replays pass; all 102,467 live comparisons pass. Fast: 5,962
  pass/32 skips; affected slow: 946 pass. Primary type diagnostics unchanged at 609.
  Measured lines: 82,635/107,682 (76.74%). Native modular Euclid/half-GCD and binary
  NTL dispatch are covered; integer/QQ and large-prime/extension fast backends remain open.
- ✅ 2026-09-10 (Codex): FLINT integer/modular resultant dependency kernels ported;
  all 1,597 new native comparisons/replays and 5,952 backend-area cases pass. Full saved
  comparison: 104,064 pass. Fast: 5,967 pass/32 skips; affected slow: 946 pass. Primary
  type diagnostics remain 609. Measured lines: 82,884/108,031 (76.72%).
- 🟡 2026-09-10 (Codex): continuing Sage resultant/Sylvester integration and remaining
  mixed-parent method coverage;
  specialized GF2 coercion and remaining
  finite-element conversion protocols, unchecked nonprime parents and dependency algorithm
  dispatch. General categories remain incomplete.
- ✅ 2026-09-09 (Codex): Rational repair batch: 84/84 public members dispatched in 5,757
  Sage comparisons, plus six Gaussian numeric-coercion regressions. Full comparative run:
  16,506 passed before the final boundaries; final affected areas: 5,757 rationals and 288
  rand_stats cases passed. Fast units: 5,886 pass, 32 existing skips. All 734 affected slow
  tests passed after fixing Gaussian coercion and rerunning a point-counting timeout in isolation.
  Current merged line coverage: 80,454/101,207 (79.49%). Existing type-check failures persist.
- ✅ 2026-09-09 (Codex): RationalField has 38/38 members in comparative coverage; repaired
  its enumeration, coercion, prime bounds and seeded randomness. Integer/Rational roots now
  share GMP’s precision-doubling Newton algorithm and preserve C exponent bounds. Validation:
  all 17,259 comparisons pass; final root/rational rerun 12,451 pass; final fast units 5,888 pass
  (32 existing skips), affected slow units 734 pass. Type-check diagnostic kinds are unchanged.
- ✅ 2026-09-09 (Codex): IntegerRing has 10/10 members in comparative coverage; both Integer
  constructors and scalar factories now follow Sage coercion and error paths. All 18,044
  comparisons pass (6,934 Integer, 6,302 Rational); fast units: 5,888 pass, 32 existing skips;
  affected slow units: 734 pass. Current line coverage: 80,525/101,343 (79.46%).
  Type-check diagnostic kinds are unchanged from the baseline.
- ✅ 2026-09-09 (Codex): finite-field repair batch has 99/99 public parent/element members
  dispatched across 11,137 comparative cases. Full suite: 29,103 comparisons pass; twelve
  additional matrix regressions also pass (925/925 matrix area). All 734 affected slow tests
  and 109 focused PARI tests pass. Fast run: 5,886 pass, 32 existing skips, two timeouts;
  both timed-out tests pass in isolation (0.17 s and 2.01 s). No new type-check diagnostic
  kinds. Instrumented line coverage: 80,562/105,221 (76.56%); more lines are now instrumented.
- ✅ 2026-09-09 (Codex): finite-field factory name/modulus repairs, degree-one Conway entries
  and PARI perfect-power delegation: 12,465 field comparisons pass; full suite 30,437 pass
  before six final whitespace regressions. Fast units: 5,888 pass, 32 existing skips; affected
  slow units: 946 pass; focused field/PARI/docs: 291 pass. No new type-check diagnostics.
  Current measured line coverage: 80,618/105,177 (76.65%). Sage-version adapters are documented.
- ✅ 2026-09-09 (Codex): prime-field square-root defaults/options repaired with 1,572 new
  comparisons. All 32,003 full-suite cases passed before twelve final keyword regressions;
  final finite-field area: 14,037 pass. Fast units: 5,888 pass, 32 existing skips; affected
  slow units: 946 pass; focused field/docs: 232 pass. No new type-check diagnostics.
  Current measured line coverage: 80,628/105,138 (76.69%). Version 1.0.0 reflects the wider
  default return type under the project's breaking-change rule, not completion of the audit.
- ✅ 2026-09-09 (Codex): finite-field arithmetic-string construction repaired: 14,904 field
  comparisons and 566 direct parser comparisons. Full suite 33,443 pass before five final
  parser-boundary cases; final parser area passes. Fast units: 5,895 pass, 32 existing skips;
  affected slow units: 946 pass; focused parser/field/polynomial/docs: 184 pass. Existing
  type-check diagnostic kinds are unchanged after normalizing the string-input union. Measured
  line coverage is 81,019/105,652 (76.68%); the parser has 23/23 implemented callables dispatched.
- 🟡 2026-09-09 (Codex): continuing finite-ring/polynomial source and input-domain coverage.
  Broader factory options and parser function/sequence/matrix branches remain open.
- 🟡 2026-09-09 (Codex): broader source/API audit in progress; `AUDIT-COVERAGE.md` inventories
  all 175 production files and distinguishes comparative coverage from executed lines.
- ✅ 2026-09-09 (Codex): September F1–F9 and the follow-up integer/coercion/curve/lattice
  defects are repaired with comparative regressions. The isolated audit is 51/51 matching;
  Integer has 107/107 implemented methods represented by 6,027 comparisons. This does not
  close the separately documented backend limits or establish full input-domain coverage.
  Validation: 10,805 comparisons pass; 5,886 fast units pass (32 skipped), initial slow tier
  1,274 pass, changed slow files rerun 246 pass. Executed production lines: 79.82%.
- ✅ 2026-09-09 (Codex): completed a focused core fidelity audit (`AUDIT-2026-09.md`):
  nine additional open findings, reproduced by 29 mismatches among 51 standalone audit cases
  (22 controls/matching cases). Existing fast/slow differential tiers: 4673/4673 pass.
  This describes the initial audit-only commit; the repairs are recorded in the newer entry above.
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
- 2026-07-28 (0.0.14) updated for the **upstream-porting pass**: the items earlier passes had
  deferred as "research-grade" or "needs an unavailable database" were ported from the vendored
  upstream instead — MPQS, `polmodular`/`polclass`/`volcano`, SEA, `buch1`, `galconj`, `qfrep`,
  the Shanks-distance `t_QFB`, Laurent and multivariate power series, polynomial matrices
  (Popov/Hermite/approximant bases) and van Hoeij. In every case the premise was wrong: nothing
  needed an external database, and PARI's own `seadata`-less fallback *is* `polmodular_ZXX`.
  New rows: `mpqs`, `polmodular`/`polclass`, `Fp_ellcard_SEA`, `buch1`, `galconj`, `qfrep`
  (parigp-ts); `laurent_series_ring.ts`, `matrix_polynomial_dense.ts`,
  `number_field_embeddings.ts` (sagemath-ts). Percentages here continue to mean "ported **and**
  verified against upstream doctests or an executed oracle", and where a module still throws,
  the row says so and names the upstream routine.
- 2026-07-29 (0.0.15) updated for the **differential-oracle pass**: nine new live-SageMath areas
  raised the oracle from 433 to 4643 cases; hyperelliptic curves, quaternion algebras, rational
  function fields, general quadratic forms and ternary quadratic forms were added; and the
  deviations register was split into accepted differences and open fidelity gaps.
- 2026-07-29 (0.0.16) closes the remaining silent-wrong-answer paths in `IsogenyClass`, integer
  matrices and extension-field minimal polynomials; adds live differential areas for the 0.0.15
  module families; converts the last 57 vacuous property cases into real Sage comparisons; and
  splits fast versus slow verification.

---

## Phase 1: Core Number Theory (Cryptography Focus)

### `sage.rings.integer` - Arbitrary Precision Integers
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | Reviewed vs reference (2026-02-04); re-audited 2026-07-27 (`AUDIT-2026-07.md`) |
| `integer.py` | 🟡 backend gaps | 107/107 method dispatch | Core Integer class; 6,027 live comparative cases with documented adapters (2026-09). Existing partition/prime-count backend limits remain |
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

### `sage.rings.generic` / `sage.rings.fast_arith`
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| `sum_of_squares.ts` | ✅ Supported IntegerLike domain audited | 11,616 comparisons with arithmetic wrappers; 107/107 native lines | All four uint32 entry points, exact errors/order, wrapped zero and upper/lower conversion bounds |
| `generic.ts` | 🟡 Integer domain audited | 12,956 product-tree/factor-base comparisons | ProductTree layers, iteration, remainders, cached interpolation and streaming product/derivative; general Euclidean-domain coefficients remain outside the current signatures |
| `fast_arith.ts` | 🟡 Factory methods audited | 13,704 factory comparisons | 32/64-bit gcd/inverse methods; separate xgcd/rational-reconstruction wrappers are not ported; native undefined behavior is documented |

### `sage.rings.real_mpfr` - Real Numbers
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| Audit (algorithm fidelity) | ✅ | - | 2026-07 (`AUDIT-2026-07.md`): 17 findings, all fixed |
| `real_mpfr.ts` | 🟡 Codex audit | 89,084 real-literal/native cases plus MPFR area | Native conversion, formatting, predicates, integer rounding, fractional parts and comparison coercion are audited. General arithmetic, directed conversion and remaining methods retain documented limitations. |
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
| `factor` | ✅ 98% | ✅ | Delegates to PARI `Z_factor`, the complete `ifac_crack` chain (trial division, pure powers, SQUFOF, Pollard-Brent, ECM, **MPQS**, insisting ECM). MPQS landed in 0.0.14, so hard semiprimes now factor; only inputs above MPQS's own 107-decimal-digit ceiling (`mpqs.h:400`) raise |
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
| `integer_mod_ring.ts` | ✅ Codex audit | ✅ | Zmod construction/iteration and cached field/factored-order/root hooks; 5,126 additional native root/lifting/cache comparisons (24.52.0). Complete modular area 9,199 comparisons and broad tests pass; category/proof framework remains outside the port profile. |
| `finite_field_prime.ts` | ✅ 100% | ✅ | GF(p) with sqrt, multiplicative_generator |
| `finite_field_extension.ts` | 🟡 Codex audit | ✅ comparative | Scalar sqrt/FF_issquareall: 4,324 native comparisons (24.53.0). Scalar is_square and native resultant/norm dependencies add 4,832 comparisons (24.54.0). Sage norm/trace and remaining caller audit continues.  GF(p^n) via polynomial quotient rings. `irreducible_element` faithfully ports `polynomial_ring.py:3560-3626`, delegating to ntl-ts `GF2X_BuildSparseIrred` and parigp-ts `ffinit`. `minpoly`/`minimal_polynomial` now delegate to parigp-ts `FpXQ_minpoly` (the same PARI branch Sage uses), with `minimalPolynomial` retained as an alias; generators, constants and proper-subfield elements are independently checked by Frobenius-orbit degree and annihilation. Remaining gaps: `algorithm='ffprimroot'` throws (needs PARI `ffgen`/`ffprimroot`/`charpoly`), `algorithm='random'` uses rejection sampling rather than NTL `BuildRandomIrred` |
| `gf2.ts` | 🟡 Codex audit | Unit tests and 85 direct comparisons | Specialized division/inverse/power errors repaired; constructor, generator and mixed-parent protocols remain open |
| `conway_polynomials.ts` | ✅ 100% | ✅ new | Regenerated by porting FLINT's `conway.c` decoder against the vendored bit-packed table. 7 entries were **reducible** and one GF(2^128) entry was fabricated; every entry is now verified irreducible, primitive, normalised and subfield-compatible. p = 2 to n=64, 3 to 24, 5 to 18, 7 to 14, 11/13 to 12, 17/19/23/29/31 to 10 |
| `tower_field.ts` | ✅ 100% | ✅ | Binary tower fields (Binius): Ti(i) = GF(2^(2^i)) |
| `roots_of_unity.ts` | ✅ 100% | ✅ | FFTDomain, CosetDomain, primitive roots |

**Extension Field Features:**
- Frobenius automorphism
- Trace and norm functions
- Minimal polynomial computation
- Integer representation round-trip
- Primitive element finding

### `sage.rings.fraction_field` - Polynomial fraction fields

| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| `fraction_field.ts` | 🟡 Codex | ✅ comparative + units | Cached generic/field/FpT parent selection, metadata and univariate arithmetic. General factory/coercion and nested coefficient domains remain open. |
| `fraction_field_element.ts` | 🟡 Codex | ✅ comparative + units | Nonconstant denominators, cross-cancellation, reduction flags, inversion/powers, accessor identity and subclasses; complete conversion and coefficient closure remain open. |
| `fraction_field_FpT.ts` | 🟡 Codex | ✅ comparative + units | Native modular arithmetic, reduction, powers and stored-pair equality/accessors; full coercion, iterators and remaining methods remain open. |

### `sage.rings.polynomial` - Polynomials
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| Audit (algorithm fidelity) | 🟡 Codex | ✅ comparative | September source audit continues; public native multiplication/power dispatch, general fraction coercion and other polynomial protocols remain open. Earlier July review fixed 37 findings but was not exhaustive. |
| `polynomial_ring.ts` | 🟡 Codex audit | ✅ | PolynomialRing with lagrange, vanishing, cyclotomic; Neville/divided-difference corrected to Sage's row+table semantics |
| `polynomial_element.ts` | 🟡 Codex audit | ✅ | ZZ/QQ roots: 699 bundled-source comparisons cover sparse ordering, NTL/PARI dispatch and zero errors (24.50.1). ZZ/QQ factor/irreducibility: 1,118 additional comparisons cover native routing, constants, zero errors and cached PARI-state effects (24.50.2). GF roots, distinct ZZ/QQ roots, extension factor order and PARI state have 5,796 additional comparisons (24.51.0); IntegerModRing roots and helpers have 5,126 new native comparisons (24.52.0), with full fast/affected-caller validation passing; finite irreducibility remains open. The local FLINT-style factor driver retains its separate internal tests. |
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
| `matrix_operations.ts` | ✅ 93% | ✅ 196 | `minpoly` (works over QQ now that `factor()` does; Sage's `x^3-30x^2-80x` doctest reproduced), `is_semisimple`, exact `is_positive_(semi)definite`, `right_kernel_matrix`, `density`, `eigenvalues`. **New in 0.0.12, redone in 0.0.14:** `norm(A,2)` now follows Sage's own route (`change_ring(CDF)`, `A^H·A`, SVD, `matrix2.pyx:16460-16471`) instead of an exact-rational eigenvalue isolation, so RR and CC entries work and Sage's four CC doctest values match; every previously pinned value reproduces bit for bit. Rings with no map into CDF are refused rather than mis-normed. Also `is_similar(transformation=true)` (verified against exhaustive brute force over 134 785 matrix pairs), a generic `change_ring`, and `is_diagonalizable(base_field)` |
| `matrix_integer.ts` | ✅ 99% | ✅ 860 | HNF, SNF, elementary divisors, kernel, exact integral LLL (`delta` 0.99), symplectic form. `frobenius_form` flags 0/1/2 port PARI's `RgM_Frobenius`. **0.0.16:** `p_minimal_polynomials`, `null_ideal` and `integer_valued_polynomials_generators` now solve the exact Smith congruences from `compute_J_ideal.py`; composite and negative moduli and rational generators match Sage coefficient-for-coefficient. `is_LLL_reduced` uses exact Rational Gram-Schmidt with Sage's default and validation, never IEEE-754 |
| `matrix_modn.ts` | ✅ 95% | ✅ | charpoly, determinant, echelonize, `right_kernel_matrix` (all three basis formats). **Composite modulus now delegates to parigp-ts `matkermod`** — 300 random matrices match SageMath entry for entry, 2190 brute-force cases confirm the returned rows generate the *full* kernel |
| `matrix_mod2.ts` | ✅ 95% | ✅ | GF(2) matrices; `pluq`/`ple` now use M4RI's transposition-list convention for both P and Q |
| `matrix_decompositions.ts` | ✅ 97% | ✅ 121 | RREF echelon form, LU, QR, Cholesky, Bunch-Kaufman `block_ldlt`, Smith, Hermite, LLL_gram. **New in 0.0.12:** `jordan_form(transformation=true)` (reproduces Sage's *exact* `P`, not merely a valid one) and `krylov_kernel_basis(variable=…)` (every matrix in Sage's docstring reproduced verbatim). `pivots` is now re-exported from `matrix/index.ts`. **0.0.14:** `jordan_form` honours `subdivide` (default `true`, as Sage) and a new exported `matrix_str` is a faithful port of `Matrix.str` (`matrix0.pyx:1834`), so five `jordan_form` doctests and all six of `subdivide`'s own doctests reproduce character for character. `krylov_kernel_basis`'s shifted-Popov property is now asserted with the ported `is_popov`/`is_hermite`/`popov_form`, not by byte-equality with Sage's printed output. **0.0.15:** `Polynomial.roots()` now follows Sage's factor ordering, fixing Jordan block order. **Remaining:** `Matrix.toString` in `matrix_generic.ts` is not subdivision-aware, so `jordan_form` attaches `matrix_str` per instance as a stopgap |
| `matrix_decompositions_additions.ts` | ✅ 100% | ✅ | SVD_double, QR_double, LU_double for IEEE 754 real matrices |
| `matrix_polynomial_dense.ts` | ✅ 80% | ✅ 81 | **New in 0.0.14.** Port of `matrix_polynomial_dense.pyx`: `degree_matrix`, `constant_matrix`, `coefficient_matrix`, `truncate`, `shift`, `reverse`, `row_degrees`/`column_degrees`, `leading_matrix`, `leading_positions`, `is_reduced`/`is_weak_popov`/`is_popov`/`is_hermite`, `weak_popov_form` (Mulders-Storjohann), `popov_form`, `reduced_form`, `hermite_form`, `minimal_approximant_basis` and `is_minimal_approximant_basis`. Every doctest in the ported functions passes verbatim; `is_weak_popov`/`is_popov` are additionally brute-forced over all 256 2x2 matrices over `GF(2)[x]` of degree <= 1. Exported from `matrix/index.ts` with the generic names aliased. **Not ported** (out of scope, no stubs): `inverse_series_trunc`, `solve_left/right_series_trunc`, `left/right_quo_rem`, `reduce`, `minimal_interpolant_basis`, `minimal_kernel_basis`, `minimal_relation_basis`, `basis_completion` |
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
| `ell_generic.ts` | 🟡 Codex audit | ✅ comparative | Same-base-field coordinate roots, Montgomery selection and coefficient repr: 5,822 new native comparisons (24.54.1). Generic execution 922/960 lines. Invariants, isomorphisms and division polynomials have earlier comparisons; general coordinate coercion/promotion, the finite torsion alias and remaining caller domains are still under audit. |
| `ell_finite_field.ts` | ✅ 96% | ✅ | cardinality, trace, generators, twists, torsion_basis, `abelian_group`, `set_order`/`has_order`, `frobenius_order`, Vélu. `is_j_supersingular` skips Sage's precomputed j-polynomial table (exact anyway under the default `proof=True`) |
| `ell_point.ts` | ✅ 96% | ✅ | Point arithmetic, weil/tate/ate pairings, `division_points` (Sage's `_multiple_x_numerator` and distinct-root algorithms; 211 native comparisons and a watched degree-5,100 caller), `is_divisible_by`, `point_log`. No p-adic shortcut for anomalous curves |
| `ell_curve_isogeny.ts` | 🟡 95% | ✅ 250-area oracle | Vélu, Kohel, fastElkies' BMSS, `compute_intermediate_curves`, `dual()`. **0.0.15:** ported Stark's algorithm, the quadratic `weierstrass_p` it needs and odd-degree `is_kernel_polynomial` validation; fixed even-degree `dual()` and implemented `formal()` for Vélu isogenies. Non-Vélu formal expansions still raise explicitly |
| `formal_group.ts` | ✅ 95% | ✅ 41 | `differential()`, `log()`, `inverse()`, `group_law()`, `mult_by_n()`, `sigma()`. **0.0.12 audit finding: these were already correct at 0.0.11 — the deferral note calling `differential()` hardcoded was stale.** This pass proved it by running Sage's own doctests verbatim (incl. the 35-term `w(35)`, the 16-term `mult_by_n(100,20)` for 37a, and the `# long time` GF(17) `mult_by_n(10,50)`) and by checking defining identities: the Weierstrass equation for `x,y` with a negative control, `log`/`exp` mutual inversion, `log(F(t1,t2)) = log t1 + log t2`, `log([n]t) = n·log t`, `[m+n] = F([m],[n])`, `F(t, i(t)) = 0`. Added `x_list`/`y_list` because `x()`/`y()` returned objects whose coefficients no caller could read. **Closed in 0.0.14:** the characteristic-zero branch of `mult_by_n` (`formal_group.py:644-665`) is ported line for line over the new Laurent series ring, and `group_law` computes in `MPowerSeriesRing(k, 't1,t2')`, so Sage's whole TESTS block over `GF(7)[[x,y,z]]` — including the genuine three-variable associativity — is verified. `x(10)`/`y(10)` now print exactly as Sage's doctests do. Two local workarounds remain around `power_series_ring.ts`'s `MPowerSeries.inv()` precision and `_subs_formal` truncation (performance, not fidelity) |
| `ell_torsion.ts` | ✅ 90% | ✅ | `_p_primary_torsion_basis` replaced with Sage's division-polynomial algorithm. Torsion over number fields still throws |
| `weierstrass_morphism.ts` | 🟡 Codex audit | ✅ | `order()` and native generator roots: 1,759 new comparisons (24.52.1); generic ordering, rich comparisons and parent guards covered by 2,371 additional comparisons (24.52.2); remaining scalar-root callers under audit. |
| `isogeny_class.ts` | ✅ 95% | ✅ + live oracle | IsogenyClass, IsogenyClassNumberField, IsogenyClassRational, `fill_isogeny_matrix`, and `Frobenius_filter`. **0.0.16:** `_compute` is the full BFS traversal with actual maps, prime/fill matrices, characteristic-zero kernel-polynomial discovery and QQ global-minimal normalization. Sage's complete three-curve `11a1` class and matrices match live. The Kenku exceptional dispatch is present; degrees 43/67/163 fail honestly until Sage's precomputed kernel table is ported |
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
| `weierstrass_morphism.ts` | 🟡 Codex audit | ✅ native comparisons | baseWI, WeierstrassIsomorphism, _isomorphisms, identity_morphism, negation_morphism; generator fixes covered by 1,759 comparisons; generic ordering/comparison/parent fixes covered by 2,371 additional comparisons. Broader field/caller audit continues. |

### `sage.schemes.hyperelliptic_curves`
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| `constructor.ts` | ✅ | ✅ 121 tests (module) | `HyperellipticCurve`, `_parse_multivariate_defining_equation` |
| `hyperelliptic_generic.ts` | ✅ | ✅ | `HyperellipticCurve_generic`, `HyperellipticPoint`, `genus_of`, `sage_poly_repr` |
| `hyperelliptic_finite_field.ts` | ✅ | ✅ + 7 live | Cartier-Manin, point enumeration, `frobenius_polynomial` and Jacobian cardinality; the live area covers odd/even characteristic, extensions, Cartier/Hasse-Witt, a-number and p-rank |
| `hyperelliptic_rational_field.ts` | ✅ | ✅ | |
| `hyperelliptic_g2.ts` | ✅ | ✅ | Genus-2 specialisation |
| `jacobian_generic.ts`, `jacobian_homset.ts`, `jacobian_morphism.ts`, `jacobian_g2.ts` | ✅ | ✅ | Mumford representation, Cantor composition/reduction |
| `invariants.ts` | ✅ | ✅ | Clebsch/Igusa/absolute Igusa invariants, Ueberschiebung |
| `hyperelliptic_padic_field.ts`, `monsky_washnitzer.ts`, `kummer_surface.ts`, `mestre.ts`, `jacobian_endomorphism_utils.ts`, `hypellfrob` | 🔴 | - | Not ported |

### `sage.algebras.quatalg`
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| `quaternion_algebra.ts` | ✅ (QQ) | ✅ + 8 live | `QuaternionAlgebra`, `QuaternionOrder`, `QuaternionFractionalIdeal_rational`, maximal orders, ideal classes; live Sage checks cover invariants, ramification, products, pairing, involution, trace, norm, inverse and characteristic polynomials |
| `quaternion_algebra_element.ts` | ✅ (QQ) | ✅ | `QuaternionAlgebraElement_rational_field` |
| `quaternion_algebra_cython.ts` | ✅ | ✅ | `basis_for_quaternion_lattice`, `intersection_of_row_modules_over_ZZ`, `rational_matrix_from_rational_quaternions` |
| Base rings other than `QQ` | 🔴 | - | Number-field quaternion algebras not ported |

### `sage.rings.function_field`
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| `constructor.ts`, `function_field.ts`, `function_field_rational.ts`, `element.ts`, `element_rational.ts` | ✅ | ✅ + 9 live | Rational function fields `k(x)` over an arbitrary constant field; live checks cover normalization, arithmetic, valuation, squares, factorization and places. The oracle exposed and fixed the general `GF(p)` constant element's missing `is_square`/`sqrt` surface |
| `order.ts`, `order_rational.ts`, `ideal.ts`, `ideal_rational.ts`, `place.ts`, `place_rational.ts`, `divisor.ts`, `valuation_ring.ts` | ✅ | ✅ | Both maximal orders, their ideals, places, divisors, Riemann-Roch (Hess 6.1, degree-one specialisation) |
| `function_field_polymod.ts` (finite extensions) | ⬜ | - | Not ported |
| Differentials, derivations, jacobians, valuations | ⬜ | - | Not ported |

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
| `free_module.ts` | ✅ 95% | ✅ 282 | FreeModule hierarchy matching Sage's, exact rank, echelonized bases, coordinates, kernels, saturation, discriminant, cardinality. **New in 0.0.12:** free modules over **non-ZZ PIDs** (`QQ[x]`, `GF(p)[x]`) via an exact `QQ(x)` fraction-field layer plus ports of `_echelon_form_PID`, `_generic_clear_column` and Sage's generic Smith normal form — 700 random spans match the vendored SageMath exactly; real `quotient` lift/project (Sage's `__quotient_matrices` over a field, `FGP_Module` over ZZ, both agreeing with Sage on 391 random cases); and a real embedded `tensor_product` (the previous `@see Reference` citation was **fabricated** — Sage's `free_module.py` has no such method). **Fixed in 0.0.14:** `intersection()` over `K[x]` now returns Sage's exact basis. The cause was not the echelon routine (upstream normalises nothing on that path) but `Matrix.integer_kernel`'s missing `self.denominator()` scaling — over `QQ[x]` that denominator is the lcm of the rational **coefficient** denominators, a non-trivial unit of `QQ[x]`; the port only cleared fraction-function denominators, so it never scaled. 914 randomly generated cases across `QQ[x]` and `GF(2,3,5,7,11,13)[x]` now agree with SageMath coefficient-for-coefficient on all four matrices, where 22 of 250 intersections were wrong before |
| `free_module_element.ts` | ✅ 95% | ✅ new 23 | Inner-product-matrix pairing, 7-D cross product, true p-norms, exact normalization, Python-style indexing. `free_module_element.test.ts` created |
| `free_module_integer.ts` | ✅ 90% | ✅ new 41 | Exact integral LLL, exact CVP/SVP enumeration, Voronoi relevant vectors and cell, q-ary lattices. `free_module_integer.test.ts` created. Note: this and `matrix_integer.LLL` are two independent LLL implementations where SageMath has one |
| `bkz.ts` | ✅ 100% | ✅ 34 | BKZ reduction, HKZ, Schnorr-Euchner enumeration (no audit findings) |

### `sage.stats.distributions` - Sampling Distributions
| Module | Status | Tests | Notes |
|--------|--------|-------|-------|
| `discrete_gaussian_integer.ts` | 🟡 Codex audit | ✅ | **All four** algorithms now, including `uniform+logtable` and `sigma2+logtable`: `dgs_bern.c` and `dgs_disc_gauss_sigma2p` are ported, and the acceptance test is dgs's real `mpfr_urandomb` comparison (0.0.11 used a different test, so even the two implemented algorithms consumed randomness differently). All four reproduce SageMath's **seeded sample streams** bit-for-bit, and `repr()` matches Sage's format. **Not implemented:** the `precision='dp'` mode (Sage documents its results as not reproducible) |
| `discrete_gaussian_lattice.ts` | 🟡 Codex audit | ✅ | GPV algorithm, coset sampling, exact (`Rational` basis/centre/GSO) with `sampleExact()`. **New in 0.0.12:** non-spherical Σ (matrix sigma, `sigma_basis`, Peikert's `r` by power iteration, Cholesky, offline samples, `_call_non_spherical`), `set_c`/`c()`/`sigma()`/`f()`, and `_normalisation_factor_zz` with a local `qfrep`. Every doctest in the vendored source reproduces. **Closed in 0.0.14:** the theta series delegates to parigp-ts `qfrep0`, and `_normalisation_factor_zz` honours `prec` through a `RealField(prec)` layer, so Sage's `round(prec=100)` doctest returns the full 28-digit `1558545456544038969634991553`. 21 values match a live SageMath character for character. **Remaining:** the `RealField` layer is a semantics re-implementation (MPFR is not vendored) whose transcendentals may differ in the last bit; the local LLL is not fpLLL |

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
| `misc/parser.ts` | 🟡 Codex | 566 comparative cases; focused unit coverage | Tokenizer and arithmetic grammar implemented for polynomial string construction; function calls, sequences and matrices remain stubs |
| `misc/randstate.ts` | ✅ 100% | ✅ | Centralized RNG + `set_random_seed` parity. **0.0.12: seeded streams now match SageMath exactly.** GMP 6.3.0 is not vendored under `reference/`, so the sources were obtained and `rand/randmts.c` (`mangle_seed`, `randseed_mt`), `rand/randmt.c` (incl. the 624-word `default_state` that seed 0 lands on) and `mpz/urandomm.c` were ported verbatim. CPython's `random.Random` is ported too, as `PythonRandom`, since `randstate.python_random()` is Sage's second generator. Verified against a C oracle linked to libgmp 6.3.0 and against SageMath 10.3's own doctest values |
| `quadratic_forms/quadratic_form.ts` | ✅ | ✅ 204 tests (module) | `QuadraticForm`, `DiagonalQuadraticForm`, `quadratic_form_from_invariants`, `QFEvaluateVector`/`QFEvaluateMatrix` |
| `quadratic_forms/quadratic_form__local_field_invariants.ts` | ✅ | ✅ | Rational diagonal form, signature, Hasse invariants, definiteness, `qfgaussred` (belongs in parigp-ts; see DEVIATIONS) |
| `quadratic_forms/ternary_qf.ts` | ✅ | ✅ | `TernaryQF`, Eisenstein reduction, zeros mod p, p-neighbours, level/disc search |
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
| FF_issquareall scalar roots | ✅ native profile | 4,324 | n=2 dependency and extension API; exact root/state comparisons, 117/117 new helper lines. FF_issquare and resultant/norm dependencies add 4,832 native comparisons (24.54.0). |
| ellinit | ✅ | 35 | elliptic/init.ts - Short/general Weierstrass, from j-invariant. `j`/`ellj` return an exact `Ratio` when non-integral |
| Point operations (elladd, ellmul) | ✅ | 65 | elliptic/point.ts - Jacobian coordinates for efficiency |
| ellcard, ellgroup, ellorder | ✅ | 36+ | elliptic/group.ts - faithful `Fp_ellcard_Shanks` and `gen_ellgroup`/`gen_ellgens` with the real Weil pairing. Verified against exhaustive point-count oracles (4542 curves), a brute-force group-structure oracle (476 runs, 0 wrong — was 85 wrong), and real PARI 2.15.4 via Sage 10.3 for primes up to 2^32 |
| ellordinate, random_FpE | ✅ | 32 | elliptic/points.ts - Find y from x, random point generation |
| ellgenerators, trace_of_frobenius | ✅ | - | elliptic/group.ts - 10 068-curve sweep |
| `Fp_ellcard_CM` (full CM table) | ✅ 100% | ✅ new | **New in 0.0.12.** All **thirteen** class-number-one discriminants (`FpE.c:624-666`, `:1282-1421`), delegating to `qfb.ts`'s `cornacchia2`. Verified against brute-force point counting (15 392 curves), Shanks (3744), a counting-independent `[#E]P = O` oracle (936 at 64/80/96 bits) and the **published SECG group orders** of secp160k1/192k1/224k1/256k1 |
| `Fp_ellcard_Schoof` / `ellcard_sea` | ✅ 100% | ✅ | **New in 0.0.12.** Schoof's base algorithm — the "S" of SEA. Verified exhaustively on all 121 104 curves over every prime `5 <= p <= 120` and against **PARI's own `ellsea` regression vectors** at 65, 70 and 101 bits. Superseded by `Fp_ellcard_SEA` below, but kept and still tested |
| `Fp_ellcard_SEA` (Elkies + Atkin) | ✅ 95% | ✅ 31 | **New in 0.0.14.** Full port of `ellsea.c` in `elliptic/ellsea.ts`: Elkies, Atkin, `match_and_sort`, the CM branch, `Fp_elljissupersingular` and `Fq_elldivpolmod`. The `seadata` package is replaced by `polmodular.ts`, which computes `Phi_L` on demand. **Note:** `ellcard` in `group.ts` still routes >= 96 bits to base Schoof and should be switched to this |
| `polmodular` / `polclass` | ✅ 90% | ✅ 70 | **New in 0.0.14.** Port of `polmodular.c`, `polclass.c` and `volcano.c`: `polmodular_ZM`/`polmodular_ZXX`/`Fp_polmodular_evalx`, the `polmodular_db_*` cache, the class-invariant predicates and `polclass0` (Hilbert class polynomials). Golden values match PARI's `polclass(D)` for fundamental and non-fundamental discriminants |
| `buch1` (`Buchquad`, `quadclassunit0`, `bnfinit`) | ✅ 85% | ✅ 35 | **New in 0.0.14.** Index-calculus class group and unit group of quadratic fields, with PARI's `t_REAL` Shanks distances, `ZM_hnflll`/`hnfspec`/`hnfadd`/`ZM_snf_group` and the GRH check |
| `galconj` (`galoisinit`, `galoisfixedfield`, `galoissubgroups`) | ✅ 85% | ✅ 62 | **New in 0.0.14.** Port of `galconj.c` plus the `perm.c`, `Zp.c` and `FpX.c` support it needs. Takes a *Galois* monic ZX; the Galois closure of a non-Galois field needs `nfsplitting0` (`base1.c:1413`), which is not ported |
| `qfrep` (theta series) | ✅ 95% | ✅ 38 | **New in 0.0.14.** `qfrep0`/`qfrep` from `bibli1.c` with `lllgramint` and `qfgaussred_positive`. Consumed by `stats/distributions/discrete_gaussian_lattice.ts` |
| `Z_factor` (real factoring chain) | ✅ 98% | ✅ 49 | **Rewritten in 0.0.12**, **completed in 0.0.14**: `tridiv_bound` + gcd-with-primorial trial division, then PARI's `ifac_crack` order — pure powers, SQUFOF, Pollard-Brent rho, ECM, MPQS, insisting ECM — driven by an `ifac_decomp` worklist. Plus `ispower.c`'s `Z_issquareall`, `is_357_power`, `is_kth_power`, `is_pth_power`, `Z_isanypower` and a real `isprimepower` that never factors `n` |
| `mpqs` (multiple polynomial quadratic sieve) | ✅ 95% | ✅ 27 | **New in 0.0.14.** Port of `mpqs.c`: self-initialising MPQS with the size-indexed parameter tables, the full/large-prime relation stores and GF(2) elimination. In-memory relation store instead of PARI's disk files; declines above 107 decimal digits as PARI does |
| `ffinit` (Adleman-Lenstra) | ✅ 100% | ✅ 24 | **New in 0.0.12.** `polarit3.c`'s full chain plus the supporting `FpX` layer, `FpX_composedsum`, the bivariate resultant and `polsubcyclo` for prime conductor. Reproduces PARI **coefficient for coefficient** for all of the first 60 primes × n ∈ [2,12] (660/660), each independently re-verified irreducible |
| `matkermod` / `matimagemod` / `matdetmod` / `matinvmod` | ✅ 100% | ✅ 41 | **New in 0.0.12.** `bb_hnf.c` specialised to the `Z/dZ` Hermite ring (Howell form, `gen_kernel`, `gen_matimage`, `gen_inv`, `gen_detops`). All 24 golden values decoded from PARI's own regression suite reproduced verbatim; kernels confirmed complete by exhaustive enumeration |
| `Qfb` family | ✅ 97% | ✅ 61 | **New in 0.0.12**, extended in 0.0.14. `qfbred`, `qfbredsl2`, `qfbcomp(raw)`, `qfbsqr(raw)`, `qfbpow(raw)`, `qfbsolve` (all 4 flags), `primeform`, `cornacchia`/`cornacchia2`, `Zp_sqrt`/`Z2_sqrt`/`Zn_quad_roots`, and the Schoenhage fast reduction. Verified against real PARI 2.15.4 on ~2500 golden values. **0.0.14 adds** PARI's `t_REAL` kernel, the `QfbExt` distance-carrying form and the `qfr3`/`qfr5` layer, so `flag \| qf_NOD` works |

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
- `ellcard`'s dispatch still sends >= 96 bits to base Schoof rather than to the new
  `Fp_ellcard_SEA` (`elliptic/group.ts:1318`, `:1357`)
- `ffgen` / `ffprimroot` / `charpoly` over `F_q` — blocks `irreducible_element(algorithm='ffprimroot')`
- `nf` module (`nfinit`, `nfmaxord`, `idealprimedec`, ideal HNF arithmetic, `nfrootsof1`,
  `nffactor`, `rnfequation`, `nfsplitting0`, `polredbest`) — the single largest remaining gap.
  It blocks `bnfinit` for degree > 2 (`buch.ts` stops at degree 2 for exactly this reason), the
  Galois closure of a non-Galois field, and `fixed_field`'s `polredbest`. A partial, misplaced
  version (`nfbasis`, `nfdisc`, `idealprimedec`, `nfgaloisconj`, `quadunit`) lives in
  sagemath-ts's `pari_nf.ts` and should be **moved** here
- `matfrobenius` — currently ported inside sagemath-ts's `matrix_integer.ts`
- `qfbclassno` — currently ported inside sagemath-ts. (`quadclassunit` and `qfrep` are ported
  here as of 0.0.14; the sagemath-ts call sites for `quadclassunit` have not been re-pointed yet)
- `mpqs_class_init` / `mpqs_class_rels` (`mpqs.c:1775`, `:1815`) — the class-group entry points
  of MPQS. Their only caller is `buch2.c`, so `buch.ts` uses upstream's own `imag_relations`
  fallback instead; the shared `MPQS_MODE_CLASSGROUP` branches are transcribed but **untested**
- `s4galoisgen` / `f36galoisgen` (`galconj.c:1519`, `:1698`) — need `FpX_ffisom`,
  `FpXQ_ffisom_inv`, `FpXV_ffisom`, `FpXV_chinese`, `FqC_FqV_mul`
- `find_isogenous_from_Atkin` / `find_isogenous_from_canonical` (`ellsea.c:900`, `:964`) — the
  only two routines in this repo that genuinely require the unvendored `seadata` package, and
  PARI itself never reaches them without it
- Weber / double-eta / Atkin class invariants in `polmodular`'s CM path — need `polmodular.c`'s
  ~1500 lines of double-eta tables plus `polclass.c`'s orientation machinery. SEA never asks
  for them (`ellsea.c:118-123` only ever uses `INV_J` and `INV_G2`, both complete)
- Transcendental functions (`dilog`, `incgam`) — currently ported inside sagemath-ts
- Extension field elliptic curves (t_ELL_Fq): `Fq_ellcard_SEA` with `T != NULL`, and the
  `FpXQE`/`FlxqE` group layer. This is what blocks residue degree > 1 in
  `Frobenius_filter` and higher-degree primes generally
- Polynomial operations (t_POL)
- ✅ **Barrel gap closed in 0.0.14:** `src/index.ts` now re-exports the `mpqs`,
  `polmodular`/`polclass`, `ellsea`, `qfrep`, `buch`, `galconj` and extended `qfb` surfaces (439
  exports, no duplicate-export error). Two deep imports in sagemath-ts were replaced by package
  imports. Note that `buch.ts` and `qfb.ts` each carry an independent transcription of PARI's
  `t_REAL` kernel and of the `qfr3`/`qfr5` containers — 24 clashing names — so only `qfb.ts`'s is
  re-exported at package level and buch's `qfr3_*`/`qfr5_*` types are aliased (`BuchQfr3`, …).
  They should be merged into one shared module. `isprimepower`, `Z_isanypower`, `Z_iroot`,
  `squfof`, `pollardbrent`, `ellfacteur`, `forprime`, `FactorOptions`, `Fp_ellcard_CM`,
  `Fp_ellj_get_CM`, `ec_ap_cm`, `Fp_ellcard_Schoof`, `Fp_elldivpol`, `Ratio` and `isRatio` are
  **still** unexported, so `arith/misc.ts` keeps its duplicate local `isprimepower`

### flint-ts (port of FLINT)
| Feature | Status | Tests | Notes |
|---------|--------|-------|-------|
| Audit (algorithm fidelity) | 🟡 Codex | GCD comparisons | 2026-09: dense integer/modular GCD kernels now implemented; mutable object APIs remain stubs. Other Sage-local FLINT primitives still require delegation |
| Package barrel (`src/index.ts`) | ✅ | ✅ new | **0.0.12 bug fix:** the barrel re-exported five *interfaces* (`nmod_t`, `fmpz_factor`, …) through a **value** `export` clause. `tsc --noEmit` passed (it elides them) but at runtime `import … from '@sagemath-ts/flint-ts'` threw `export 'nmod_t' not found`, i.e. **the package could not be imported at all** — pre-existing since the first commit. Fixed with `export type`. `index.test.ts` created; flint-ts had **zero** tests before, which is why it survived |
| fmpz (integers) | ⬜ | - | |
| fmpz_poly GCD array kernels | ✅ Codex | Direct FLINT C comparisons | Subresultant, heuristic and modular dispatch; Sage GCD/squarefree delegation. Existing object API remains stubbed |
| fmpz_mod_poly | ⬜ | - | |
| nmod_poly GCD/monic array kernels | ✅ Codex | Sage/FLINT comparisons | Euclid, Newton division and recursive half-GCD through degree 3,500; existing object API remains stubbed |

### ntl-ts (port of NTL)
| Feature | Status | Tests | Notes |
|---------|--------|-------|-------|
| Audit (algorithm fidelity) | 🟡 Codex | Native comparative fixtures | September audit covers integer, word/arbitrary-modulus polynomial, matrix and factorization array adapters; object APIs and downstream factorization/state paths remain open. See the current checkpoints above and `AUDIT-2026-09.md` |
| ZZ | ⬜ | - | |
| ZZ_p | ⬜ | - | |
| ZZ_pX | 🟡 Codex | 1,499 new native quotient/CRT traces | Quotient/multiplier array APIs, retained rebuild errors and CRT through 1,024 primes verified in 24.45.0; legacy coefficient-object methods and Hensel state propagation remain open |
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
| `DiscreteGaussianSampler` | 🟡 Codex audit | HIGH | Seeded integer/lattice algorithms and constructor/state behavior compared; remaining coercion and unported precision domains stay open |
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
| `TernaryQF` | ✅ | LOW | Eisenstein reduction, zeros mod p, p-neighbours and level/discriminant search |
| `BQFClassGroup` | 🟡 60% | MEDIUM | Class group of binary QFs — reduced representatives form a genuine class group (28 discriminants verified against the literature); no dedicated `BQFClassGroup` class |
| `qfbsolve` | ✅ 97% | MEDIUM | **New in 0.0.12** as `BinaryQF.solve_integer` and parigp-ts `qfbsolve` (all four PARI flags). Now backed by the complete `Z_factor` chain (MPQS included), with an optional precomputed-factorisation escape hatch |
| `qfsolve` (Simon's algorithm, rational solutions) | ⬜ | LOW | Solve quadratic equations |
| `least_quadratic_nonresidue` | ⬜ | LOW | Find smallest QNR mod p |

**Use Cases:** Class group cryptography, lattice reduction, quadratic sieve

### `sage.rings.number_field` - Algebraic Number Theory (PARTIAL)
Needed for advanced algebraic crypto constructions.

| Module | Status | Priority | Notes |
|--------|--------|----------|-------|
| Audit (algorithm fidelity) | ✅ | - | 2026-07 (`AUDIT-2026-07.md`): 17 findings, 16 fixed. **343 tests** (was 228) |
| `number_field_embeddings.ts` | ✅ 90% | ✅ 16 | **New in 0.0.14.** Certified archimedean embeddings: PARI's Uspensky/VCA real-root isolation (`rootpol.c` `usp`/`X2XP1`/`ZX_Uspensky`/`ZX_realroots_irred`) over exact dyadic arithmetic, plus SageMath's own certified `complex_roots` + `refine_root` over a dyadic complex-interval layer. Every box is **proved** to hold exactly one root by the interval-Newton inclusion test. Backs `embeddings`/`real_embeddings`/`complex_embeddings`/`places`, which previously threw. Verified against Sage's `x^3+2` (53- and 100-bit), `x^3-2` and `x^5+x+17` doctests and against mpmath at 25 decimals |
| `pari_nf.ts` | ✅ 92% | ✅ | The number-field kernel: `nfbasis`/`nfdisc` (Round 2), `idealprimedec`, `polisirreducible`. **0.0.12:** `nfgaloisconj`'s degree-8 cap is **gone** (replaced by LLL reconstruction with a two-sided proof — exact `g(beta)=0` verification and a Gram-Schmidt non-existence certificate); `idealprimedec` gained Buchmann-Lenstra **round 4** (`base2.c:2248`), so inessential discriminant divisors work; and `quadunit`/`quadunitnorm` (`quad.c:281`) were added. Belongs in parigp-ts — see DEVIATIONS.md |
| `NumberField` | ✅ 85% | LOW | Basic operations, element arithmetic, norm, trace, real field discriminants, integral bases |
| `NumberFieldElement` | ✅ 85% | LOW | Arithmetic, charpoly, minpoly, is_integral, `is_unit` (field semantics) + `is_integral_unit` |
| `CyclotomicField` | ✅ 85% | MEDIUM | Q(zeta_n) - cyclotomic polynomials, degree, exact `automorphisms()` at any degree |
| `QuadraticField` | ✅ 85% | LOW | Q(sqrt(D)) using `x^2 - D` verbatim as Sage does (it previously rewrote the polynomial to the squarefree part), discriminant |
| `RationalPolynomial` | ✅ 90% | - | Supporting polynomial arithmetic over Q |
| `class_group` | 🟡 60% | LOW | **Quadratic** class groups computed exactly (imaginary and real, incl. Sage's C38 x C2 doctest), guarded at \|D\| <= 2e6. **0.0.12:** degree > 2 answers in the one rigorous sub-case — when the Minkowski bound provably admits **no** prime ideal, `h = 1` (gives Sage's answer for `x^3-x-1`, `x^4-x-1`, `Q(zeta_5)`, `Q(zeta_7)`, `Q(zeta_12)`). **0.0.14:** the certificate also accepts a *proof of principality* for each factor-base prime (exhibit `alpha in P` with `|N(alpha)| = N(P)`), which closes `Q(2^(1/3))`, `Q(sqrt2,sqrt3)`, `Q(zeta_7)`, `Q(zeta_8)`, `x^5-x-1`, `x^6+243` and Dedekind's `x^3-x^2-2x-8`, and reproduces Sage's `[1,1,1]` Hecke doctest. Fields with `h > 1` (`x^3-19`) and fields we cannot certify (`Q(zeta_23)`) still throw — we will not return an unproved answer. Structure for `h > 1` needs PARI `bnfinit` (`Buchall_param`, `buch2.c:3946`), which `parigp-ts/src/buch.ts` deliberately stops short of at degree 2 |
| `unit_group` | 🟡 75% | LOW | **0.0.12:** real quadratic fundamental units implemented (PARI `quadunit_uv_basecase`), so `fundamental_units()`, `units()` and `regulator()` answer; `UnitGroup.log()` for rank 1 is exact (sign test + `O(log k)` doubling/binary search). `regulator()` is overflow- and cancellation-safe (the ~250-digit unit of `Q(sqrt(1000003))` gives a finite 576.646). **0.0.14:** `nfrootsof1` proves the number of roots of unity and a generator, so `unit_group()` no longer silently claimed torsion order 2 for **every** degree > 2 field (it claimed 2 for `Q(zeta_5)`, where `w = 10`); when the certificate is inconclusive `zeta_order()`/`torsion_order()` throw rather than invent a value. Sage's `x^4-8x^2+36` (C4 x Z) and `x^4-x^2+4` (`w = 6`) doctests reproduce. `fundamental_units()` and `regulator()` for degree > 2 still need `bnfinit`'s `getfu`/`makeunits` |
| `galois_group` | 🟡 70% | LOW | Built from genuine automorphisms — and **now at any degree**, since `nfgaloisconj`'s cap is gone: `is_galois()` for degree > 2 used to return `false` unconditionally. **0.0.14:** `fixed_field`/`fixed_field_data`, `decomposition_group`, `ramification_group`, `inertia_group`, `ramification_breaks`, `frobenius`/`artin_symbol`, `complex_conjugation`, `ramification_degree` and `galoissubgroups` work at **any degree**, by delegating to parigp-ts's `galconj`. Sage's `x^4+1 -> x^2+2` and `x^5-5x^2-3` `polred=False` doctests reproduce character for character, and the degree-8 `ramification_breaks == [1,3,5]` doctest holds. **Still throws** for non-Galois fields: Sage returns the Galois closure's group, whose route (`splitting_field.py:371`, or PARI `nfsplitting0`) needs `nffactor` + `rnfequation`, neither of which exists here. `fixed_field` also omits Sage's `polredbest` post-processing (`polred.c` not ported), so it returns Sage's `polred=False` model |
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

**Not Implemented (still throw, naming the missing upstream routine):**
- **Class group structure for `h > 1`**, and `regulator()` / `fundamental_units()` for degree > 2
  — all three need PARI `bnfinit` (`Buchall_param`, `buch2.c:3946`). That is not a
  "research-grade" gap but a size one: `Buchall_param`'s first 40 lines already need
  `nfinit_basic`/`nfinit_complete`, `nfmaxord` (round 4), `idealprimedec`, ideal HNF arithmetic,
  the `T2` archimedean form and `nfrootsof1`, i.e. ~20k lines of `base1.c`–`base5.c` that
  `parigp-ts` does not have at all. It is also circular: a rigorous principality test needs the
  units, and the units come out of the same relation search
- **Galois group of a non-Galois field** — Sage returns the Galois closure's group; both candidate
  routes (`splitting_field.py:371` and PARI `nfsplitting0`, `base1.c:1413`) need `nffactor`
  (factorisation over a number field, `nfsqff`/`base5.c`) plus `rnfequation`, neither of which
  exists in this repo
- **`fixed_field`'s `polredbest` post-processing** (`polred.c` not ported) — we return Sage's
  `polred=False` model, an isomorphic but uglier defining polynomial

*Resolved in 0.0.12:* fundamental units and regulators of real quadratic fields (`quadunit`),
`nfgaloisconj` at any degree (LLL), and `decomposition(p)` at inessential discriminant divisors
(Buchmann-Lenstra round 4).

*Resolved in 0.0.14:* archimedean embeddings (`embeddings`/`real_embeddings`/
`complex_embeddings`/`places`, certified — this was **not** the same gap as `bnfinit`, as the
0.0.12 note claimed); the whole Galois layer for Galois fields of any degree (`fixed_field`,
`decomposition_group`, `ramification_group`, `inertia_group`, `ramification_breaks`,
`frobenius`/`artin_symbol`, `complex_conjugation`, `galoissubgroups`) via parigp-ts `galconj`;
`NumberFieldIdeal.valuation` at primes of residue degree > 1; the number of roots of unity
(`nfrootsof1`), proved; and class number 1 whenever every factor-base prime can be proved
principal.

**Use Cases:** Ring-LWE over cyclotomic fields, class group crypto, NTRU

### `sage.rings.padics` - p-adic Numbers
| Module | Status | Priority | Notes |
|--------|--------|----------|-------|
| Audit (algorithm fidelity) | ✅ | - | 2026-07 (`AUDIT-2026-07.md`): 16 findings, all fixed. 148 tests (was 86 assertions) |
| `Zp(p, prec)` | ✅ 90% | LOW | p-adic integers with capped-relative precision. **`add()` multiplied by p^v twice for operands of unequal valuation** until 0.0.11 |
| `Qp(p, prec)` | ✅ 90% | LOW | p-adic field |
| `pAdicGenericElement` | ✅ 90% | LOW | `nth_root`, `square_root`, `is_square`, `multiplicative_order`, `additive_order`, `expansion`, `slice`, repr all ported from upstream |
| `power_series_ring.ts` | ✅ 90% | LOW | `log` (integral of the logarithmic derivative), `nth_root` (Newton via `_nth_root_series`), `pade` (returns a `PadeApproximant` — no `Frac(R[z])` type here). **0.0.14:** `__call__` now follows `power_series_poly.pyx:176` exactly (the old version silently lost precision), plus `truncate_powerseries`, `valuation_zero_part`, `is_monomial`, `V`, `shift`, `is_square`, `O`, negative shifts, the `is_one` fast path in `__invert__` and `sqrt`'s `test_exact` step. New `MPowerSeriesRing`/`MPowerSeries` (any number of variables, total-degree precision, negdeglex repr) from `multi_power_series_ring{,_element}.py`. ⚠️ `MPowerSeries.inv()` does **not** match upstream's precision (upstream inverts the background univariate series) — a known open defect, worked around locally in `formal_group.ts` |
| `laurent_series_ring.ts` | ✅ 90% | LOW | **New in 0.0.14.** Port of `laurent_series_ring.py` and `laurent_series_ring_element.pyx`: the ring (`gen`, `zero`/`one`, `change_ring`, `default_prec`, `is_field`/`is_exact`/`characteristic`, `residue_field`, `uniformizer`, `power_series_ring`, `_element_constructor_`, `_repr_`) and the element (`t^n·u` normalisation, `+ - * / ~ ^` with rational exponents, `valuation`/`valuation_zero_part`, `prec`, `degree`, `add_bigoh`/`O`, the three `truncate` forms, `shift`/`<<`/`>>`, `common_prec`/`common_valuation`, `list`/`coefficients`/`exponents`/`__getitem__`, `residue`, `lift_to_precision`, `V`, `_derivative`, `integral`, `nth_root`, `reverse`, `is_square`, `power_series`, substitution, and upstream's exact `_repr_` string-replacement chain). Every test value is a verbatim SageMath doctest. **Not ported:** ordering comparisons (`<`, `>`) — the port's `RingElement` interface has no order on coefficients; `laurent_polynomial()`, `_latex_`, `__pari__`, `_im_gens_`, slicing, `_unsafe_mutate`; ring-level `random_element`, `construction`, `fraction_field` |
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
`DEVIATIONS.md` at the project root is the single source of truth (**56 numbered entries as of
0.0.15, including its authoring template**). It is split into accepted differences and open
fidelity gaps; fixed gaps are removed and retained only in `CHANGELOG.md`.
Headlines:
- `GF()` returns `FiniteFieldPrime` for primes; `FiniteField` is aliased to `GFExtended` and accepts prime powers
- Element coercion is explicit via `__call__()` rather than automatic
- Uses TypeScript's `bigint` for arbitrary precision integers
- **The RNG now matches SageMath's seeded streams exactly** (GMP `randseed_mt` + CPython's `Random` both ported); the 0.0.11 caveat about a different seeding step is gone
- Several PARI routines are still ported *in place* rather than delegated (`nf`, `matfrobenius`, `qfbclassno`, transcendentals), because the target package lacks the module — but `ffinit`, `matkermod`, the `Qfb` family, `GF2X`, `qfrep`, `galconj`, `polmodular` and MPQS are now genuinely delegated

### Known Limitations
- Conway database covers p = 2 to n=64, 3 to 24, 5 to 18, 7 to 14, 11/13 to 12, 17/19/23/29/31 to 10. Outside it the default modulus is now **SageMath's own** `ffinit` / `BuildSparseIrred` choice, so it agrees with Sage wherever our Conway coverage matches Sage's — but Sage's table is much larger (it has entries for `37^2`, `97^2`, `2^100`, `19^21`, …), so element representations of those fields are not interoperable and the generator need not be primitive
- **`Z_factor` is complete** as of 0.0.14 (MPQS included); only inputs above MPQS's 107-decimal-digit ceiling raise
- **SEA is ported** as of 0.0.14 (`Fp_ellcard_SEA`), with `Phi_L` computed on demand by `polmodular.ts` instead of the unvendored `seadata`. `ellcard`'s own dispatch still routes >= 96 bits to base Schoof and should be switched
- **Class number / class group of a degree > 2 field**: quadratic fields can now delegate to `Buchquad`/`quadclassunit0`; degree > 2 still needs the `nf` layer
- **Polynomial factorization now has van Hoeij/LLL recombination** (0.0.14), so the 200 000-subset budget and the bounded prime scan are gone
- **Two performance gaps are open, not fidelity gaps:** `ellcard`'s dispatch (above) and `MPowerSeries.inv()`'s precision in `power_series_ring.ts` (upstream inverts the background univariate series; the port sums the geometric series, so precision grows and the loop never short-circuits). `formal_group.ts` works around the latter locally
- **`Matrix.toString` is not subdivision-aware** and pads per column rather than to Sage's single global width; `matrix_str` (the faithful port of `matrix0.pyx:1834`) exists and `jordan_form` attaches it per instance as a stopgap
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
| sage.arith | 47+ functions (incl. rational_reconstruction, CRT_basis, continued_fraction, trial_division, prime_powers, smooth_part, sum_of_squares) | Proven primality (APRCL/ECPP). MPQS landed in 0.0.14 |
| sage.rings.finite_rings | 9 modules (incl. tower_field, roots_of_unity, log in Z/nZ) | element.log() for GF(p^n); a wider Conway table; `ffprimroot` |
| sage.rings.polynomial | 10 modules (incl. FFT, Schönhage convolution, multivariate, ideals, GF2X, van Hoeij/LLL) | 9 term orders; Singular-grade Gröbner |
| sage.matrix | 10 modules (incl. HNF, SNF, mod2, modn, special, decompositions, Jordan + transformation + subdivisions, Frobenius flag 2, `matkermod`, polynomial matrices) | `rational_form`; `norm(2)` over number fields (needs a distinguished embedding); subdivision-aware `Matrix.toString`; symbolic base rings |
| sage.schemes.elliptic_curves | 11 modules (incl. pairings, torsion, twists, complete rational isogeny-class traversal, Stark/BMSS isogeny kernels, formal group over a real Laurent ring, CM, Bernardi sigma, supersingular alpha, `Frobenius_filter` over number fields) | Exceptional degree 43/67/163 isogeny kernel table; torsion over Q/number fields; non-Vélu `formal()`; **modular symbols** (gates most of `padic_lseries`); residue degree > 1 in `Frobenius_filter` |
| sage.crypto | 4 modules (lattice, lwe, boolean_function, sbox) | Seeded parity for `gen_lattice`'s `ideal`/`cyclotomic` branches |
| sage.modules | 4 modules (free_module incl. non-ZZ PIDs / quotients / tensor products, free_module_element, LLL, BKZ, CVP, SVP, Voronoi) | — (`intersection()`'s basis over `K[x]` was fixed in 0.0.14 and now matches Sage exactly) |
| sage.stats.distributions | 2 modules, all four integer algorithms + non-spherical Σ, arbitrary-precision `_normalisation_factor_zz`, `qfrep` delegated to parigp-ts | `precision='dp'` (Sage documents its results as not reproducible); full MPFR sampler arithmetic/high-precision parameters remain incomplete |
| sage.coding | 4 modules (RS, BCH, Goppa, Reed-Muller) | Brouwer-Zimmermann minimum distance |
| sage.groups.generic | ✅ 95% | `bounds`/`algorithm`/`verify`; `linear_relation`, `merge_points`, `structure_description` |
| sage.quadratic_forms | ✅ 97% BinaryQF plus `QuadraticForm`, local-field invariants and `TernaryQF` | `BQFClassGroup` class and `qfsolve`. Shanks distance forms landed in parigp-ts `qfb.ts` in 0.0.14 |
| sage.rings.number_field | ✅ 92% (archimedean embeddings, the Galois layer at any degree, valuations at any residue degree, proved `nfrootsof1`, certified `h = 1`) | `bnfinit` for degree > 2 (class group structure, regulator, fundamental units) and the `nf` layer it needs; Galois closure of a non-Galois field (`nffactor`); `polredbest` |
| sage.rings.real_mpfr / complex_mpfr | 🟡 Real audit / complex 90% implementation | Real conversion/observations now use native MPFR state; general arithmetic and complex precision remain limited (see DEVIATIONS.md) |
| sage.rings.padics | ✅ 90% | Extension fields (Zq) as a first-class module; minimal_polynomial/charpoly |

**Test Coverage (2026-07-29, 0.0.16, verbatim from the runners):**
- Exhaustive unit partition (`bun run test:fast` + `bun run test:slow`): **7122 pass, 32 skip,
  0 fail**, 3 025 438 expect() calls across 125 files. Fast: 5848 pass / 32 skip across 102 files
  in 39.59 s; slow: 1274 pass across 23 files in 793.60 s
- `bun run test:property` (live SageMath 10.3 comparison): **4673/4673 passed, 0 failed,
  0 errors** across 23 areas — `coding_crypto` 926, `matrix_extended` 913, `mpfr` 815,
  `padics_series` 462, `groups_modn` 335, `rand_stats` 274, `ec_advanced` 250, `lattices` 159,
  `quadratic_forms` 82, `finite_fields` 78, `arith` 77, `arith_special` 57, `polynomials` 45,
  `elliptic_curves` 43, `polynomial_ops` 30, `arith_extended` 26, `number_fields` 24,
  `matrix_ops` 22, `matrix` 16, `lwe` 15, `function_fields` 9, `quaternion_algebras` 8,
  `hyperelliptic` 7

**Test Coverage (2026-07-29, 0.0.15, verbatim from the runners):**
- `bun test`: **7101 pass, 32 skip, 0 fail**, 3 024 278 expect() calls across 125 files
- `bun run test:property` (live SageMath 10.3 comparison): **4643/4643 passed, 0 failed,
  0 errors** across 20 areas — `mpfr` 815, `matrix_extended` 908, `coding_crypto` 926,
  `padics_series` 462, `groups_modn` 335, `rand_stats` 274, `ec_advanced` 249, `lattices` 159,
  `quadratic_forms` 82, `finite_fields` 78, `arith` 77, `arith_special` 57, `polynomials` 45,
  `elliptic_curves` 43, `polynomial_ops` 30, `arith_extended` 26, `number_fields` 24,
  `matrix_ops` 22, `matrix` 16, `lwe` 15 (0.0.14 baseline: 433 cases)

**Test Coverage (2026-07-28, 0.0.12, verbatim from the runners):**
- `bun test`: **6216 pass, 32 skip, 0 fail**, 2 738 804 expect() calls across 106 files
  (0.0.11 baseline: 6208 pass / 33 skip / 0 fail across 105 files)
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
