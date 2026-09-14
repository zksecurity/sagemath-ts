Read-only follow-up; no production changes yet.

FLINT word is_irreducible in reference/flint/src/nmod_poly_factor/is_irreducible.c:
- length <=2 returns true at low level (Sage wrapper rejects zero/units first).
- trial x=1..p-1 only when p<=max(200,2*length), then squarefree and Shoup coarse
  distinct-degree irreducibility; do not substitute Rabin.
- Existing flint ports supply all relevant Brent-Kung vector/precomputed matrix,
  powmod, inverse-series and multiplication dependencies. factor_distinct_deg.ts
  already has the baby/giant setup; copy the native irreducibility control flow,
  which exits at first nontrivial interval gcd and never refines/removes factors.
- Sage polynomial_zmod_flint.pyx672 cached method: zero false, unit false,
  composite-characteristic nonunits raise exact NotImplementedError (constants
  included), otherwise delegates low-level nmod_poly_is_irreducible.

GF2: bundled polynomial_gf2x.pyx269 caches and calls NTL GF2X_IterIrredTest.
Existing packages/ntl-ts/src/GF2X.ts isIrreducible at555 faithfully implements
GF2XFactoring.cpp batched-GCD iterative test (not Rabin); exported IterIrredTest
already delegates to it. Use existing dependency rather than duplicating.

Large prime modular polynomial: no specialized override found in
polynomial_modn_dense_ntl.pyx. Generic polynomial_element.pyx10115 cached method
rejects zero/units, constant nonunits delegate to coefficient.is_irreducible,
otherwise asks base hook then factors. Check parent hook absence before coding.

Extensions: polynomial_zz_pex.pyx320 is uncached and has options
algorithm='fast_when_false',iter=1. It restores modulus then dispatches
NTL ZZ_pEX_IterIrredTest / DetIrredTest / ProbIrredTest, unknown->ValueError
('unknown algorithm'). Current sagemath port ignores all options and uses Rabin.
Missing dependency ZZ_pEXFactoring ports must be added there first. Existing
ZZ_pEX.ts supplies power,InvTrunc,mul,PowerMod,PowerXMod,XGCD; no modular-composition
argument/Frobenius/PowerCompose helpers currently found.
NTL IterIrredTest uses FrobeniusMap, composition table size2*sqrt(degree), and
batched GCD with batch sizes4,9,16,... . DetIrredTest uses PowerCompose then degree
factor tree/RecIrredTest. ProbIrredTest uses random degree<n polynomial, TraceMap,
rejects positive-degree traces; with iter<=0 it may deliberately report true
(e.g. odd degrees), so tests must not replace its probabilistic contract with
mathematical irreducibility. Plan native random/state adapter explicitly; the
existing public polynomial NTL contexts have documented ambient-state limits.
