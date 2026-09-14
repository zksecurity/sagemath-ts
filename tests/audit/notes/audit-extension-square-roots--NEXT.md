After committing 24.53.0, continue to predicate/norm dependencies, then callers.
No subagents. User still wants full coverage, each bug comparative tested.

Completed scalar sqrt behavior:
- FF_issquareall root-output-pointer adapter in parigp ff.ts delegates internal
  _finite_field_square_root.ts. n=2 specialization of actual FF_ispower route;
  no alternative Fp2_sqrt or polynomial factorization. Public extension sqrt
  added with extend=false default, all=true lists, exact errors.
- Word extension sqrt(1) MUST NOT short-circuit: native
  Flxq_sqrtn_spec_pre uses gequal1 on t_VECSMALL (always false). Large polynomial
  backend and binary/linear generic group equal1 do short-circuit. 15 initial
  root/state discrepancies caught and repaired; permanent cases cover them.
- Native C FF_issquareall roots/state oracle, Python verifies actual Sage wrapper
  option errors and unordered successful root sets/parent identity against C.
- Next predicate must NOT call randomized sqrt just to test squareness.

Predicate source:
- element_pari_ffelt.pyx:1001 -> FF_issquare (FF.c:896).
- Binary always true. Word Flxq_issquare (Flx.c:3787): zero true, otherwise
  krouu(Flxq_norm(x,T,p),p)==1. Large FpXQ_issquare (FpX.c:2516): zero/char2
  true, constant delegates Fq_issquare (even degree true else Fp_issquare),
  otherwise kronecker(FpXQ_norm(x,T,p),p)!=-1.
- Both norm functions (FpX.c:3017,Flx.c:3839) compute resultant(T,x), then
  divide by lead(T)^deg(x) if T nonmonic and x nonzero. No norm port exists.
- IMPORTANT: Sage extension norm is inherited element_base.pyx:632 and uses
  charpoly('x')[0] with degree sign, NOT direct FF_norm. Existing TS norm is a
  product of Frobenius conjugates. Do not repair that by guessing FF_norm;
  charpoly/minpoly native algorithm must be audited first. Trace uses
  __pari__().trace().lift(), so FF_trace is its dependency.
- No FpX_resultant/Flx_resultant ports. Native FpX.c:1046/Flx.c:2387 use
  Euclidean basecase below native GCD limit and half-resultant recursion above.
  Do NOT substitute generic O(n^3) determinant or Euclidean-only large inputs.
- _polynomial_gcd.ts already implements the same halfGCD recursion with no
  resultant state. Native halfres (FpX.c:635-775,989;Flx.c:1956-2101,2331)
  is the SAME algorithm with optional {res,lc,deg0,deg1,off} state. Could add
  optional accumulator to halfInternal/polynomialHalfGcdBasecase, preserving
  all old callers without it. Read exact transitions before implementation.
  Reuse matrix helpers and threshold config; port wrapper and basecase with
  native sign/nonmonic factors. Comparative tests should include threshold
  crossings, split recursion, deg drops, common factors, nonmonic and zero.
- Native scalar resultant top-level preserves zero integer output; private
  FpX_resultant_basecase returns polynomial zero only for private invalid-ish
  path but exported output may therefore be tagged depending native call.
  Verify output types against C if exposing this boundary.

Remaining curve and hyperelliptic issues as previous NEXT.md:
../audit-generic-isomorphism-order/NEXT.md
- ell_generic._square_roots still local Tonelli and misses extension
  nonresidues/QQ. _is_square uses Euler and unknown-field order assumptions.
- _poly_roots still default root+positive resort and compareFieldElements
  still strings; fix original is_x_coord/lift_x/montgomery callers source-first.
- hyperelliptic sqrt_all_unsorted now detects new scalar sqrt, so falls through
  less often, but sorts via wrong compare_elements; still misses native
  random consumption for nonsquares by doing an extra predicate first.
- hyperelliptic first distinct/default polynomial root callers and mixed-parent
  coercions remain open. Finite irreducibility notes in separate audit folder.
