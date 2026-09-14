Next after committing 24.47.0 Hensel changes: integer GCD/SquareFreeDecomp state.
No production changes for this next task have been applied yet.

Read sources:
- reference/ntl/src/ZZX1.cpp:2790 content,2805 PrimitivePart,2838 GCD.
- GCD handles zeros before FFT init; removes content; each trial calls
  zz_p::FFTInit(i) BEFORE rejecting primes dividing either leading coefficient;
  word GCD then scales by gcd of leading coefficients. Degree decrease resets
  CRT, degree increase skips; stable CRT candidate is primitive-normalized and
  certified by exact divide(f1,res) && divide(f2,res), short circuit. State must
  pass into certification divide too. Restore context and scale positive gcd.
- SquareFreeDecomp reference/ntl/src/ZZXFactoring.cpp:101 uses GCD plus exact
  divide at every stage. Existing _ZZX_kernels.divide accepts optional state.
- Current ZZX1.GCD uses private deterministic prime(i), big-modulus XGCD,
  and drops state in divide certification. prime(i) still used by divide when
  state is omitted, so do not remove that helper automatically.

Prerequisite: special FFT-prime word coefficient contexts.
- reference/ntl/src/lzz_p.cpp:88 constructor INIT_FFT uses info->q, NumPrimes=1,
  PrimeCnt=0, MaxRoot=CalcMaxRoot(q), p_info=selected prime info.
- zz_pContext(INIT_FFT,index) rejects negative with 'bad FFT prime index' BEFORE
  UseFFTPrime, whose own negative error is misspelled 'invalud FFT prime index'.
- Crossovers lzz_pX.cpp:11..18 arrays at index 0:
  mod45,mul150,newton150,div180,halfgcd90,gcd400,bermass400,trace200.
  Index1 is identical EXCEPT bermass480 and trace350.
- Current port conflates PrimeCnt with NumPrimes. Need a NumPrimes field for
  buffer allocation, CRT dimensions and transform-row counts; retain PrimeCnt
  for algorithm cutoffs. Existing ordinary contexts have both equal, 1..4.
- wordFFTPrimes must use selected FFT index's info directly in FFT mode, not
  enumerate first PrimeCnt (which is zero) or even first NumPrimes (index 0).
- Multiplier allocated-prime-count guard uses NumPrimes, not PrimeCnt.
- wordHalfGCDKernel.divrem cutoff must treat PrimeCnt 0 like 1; its half/gcd
  cutoff already based on modCrossover works. MinPolySeq needs explicit 0->400.
- Native TofftRep/FromfftRep use p_info directly. FromfftRep skips CRT and copies
  inverse FFT residues; generic static FromModularRep is NOT valid with p_info
  (its CRT coefficient tables are unset). Port transform helpers can use a
  specialized identity-mod-q reconstruction, but do not call the native static
  CRT function in a new FFT-context oracle.
- Proposal (not applied/committed): add optional fftPrime index to zz_pXOptions,
  capture index and state in F.options, initialize only that FFT prime. Require
  explicit p to match its q (document this adapter guard). In this mode use
  full CalcMaxRoot(q), as native FFTInit takes no maxroot. Nested operations must
  retain fftPrime, including normal GCD/word mul and DivRem's fresh modulus.
  An internal factory or mirrored lzz_p.ts context adapter is an alternative;
  decide cleanly and avoid inventing ambient mutable global contexts.
- No lzz_p.ts exists; ZZ_p.ts has old object/global-context stubs. Do not implement
  those unrelated big object APIs just to support this bounded dependency.

New oracle planning:
- FFT-context word oracle can derive from ntl_word_rebuild.cpp (complete matching
  bundled NTL, includes lzz_pX.cpp). Input key,index,lengths,coeffs,commands;
  init via zz_p::FFTInit(index), derive p from native context. Cover word ops
  0..25 EXCLUDING raw static-CRT op12 in this mode; add MinPolySeq around400 if
  this option is exposed publicly. Compare full outputs/errors/cache/stream.
- Integer GCD/SquareFree oracle can derive from ntl_hensel_state.cpp boilerplate,
  but call full bundled GCD and SquareFreeDecomp directly. No need extract their
  bodies unless adding diagnostic counters. Operations for GCD, squarefree and
  warm integer products permit shared-cache continuation. Native p=0/no initial
  word context is useful, with optional prewarming from existing explicit p.
- Cases need zeros, signed contents/leading coefficients, primes dividing leading
  coefficients (first several deterministic q), unlucky modular gcd degree,
  CRT degree reset/skip, stabilization and failed then successful certification,
  high-degree/crossover cases and repeated squarefree factors.

Additional oracle issue to investigate, not a confirmed production bug:
- ntl_word_rebuild TS op10 reconstructs default F with {maxroot} but omits state;
  native F=zz_pXModulus() retains ambient coefficient context/global FFT roots.
  Existing cases passed, but a default F followed by rebuild and old B with a
  partial stale FFT prefix may expose root-dependent numerical differences.
  New FFT oracle must pass the full captured options/state when clearing F.
  Probe old word oracle with nonzero stream key, oldB degree200, clearF, newF
  degree250/256 (same k/l/NumPrimes), before changing adapters or claiming bug.
- Big oracle constructor ops11/12 use C++ holder copy-assignment vs TS new-object
  rebinding. New across-context constructor tests need unique_ptr ownership on
  native side, not accidental retained allocation from assignment.
