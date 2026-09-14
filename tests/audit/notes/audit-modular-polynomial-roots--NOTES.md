Next audit, read-only findings. No production changes or fixtures yet.

Bundled integer_mod_ring.py:756 field(): cached GF(order()) if is_field(), else
ValueError('self must be a field'). Port IntegerModRing has no field() method.
Use the existing GF factory (FiniteFieldPrime adapter) and cache per ring instance.
Add shared modular_integers comparisons and executed LLM examples for new API.

integer_mod_ring.py:1960 _roots_univariate_polynomial default multiplicities=True:
- zero or nonfield -> NotImplementedError root finding with multiplicities for
  this polynomial not implemented (try the multiplicities=False option).
- prime/nonzero -> change_ring(f.base_ring().field()).roots(), retaining returned
  GF parent, not IntegerModRing. Returned coefficient class therefore changes.
Bundled generic polynomial_element.pyx:8895 catches NotImplementedError and falls
through. For zero/prime, bundled factor raises ArithmeticError, then generic
roots translates to the same NotImplementedError. Composite constants also reach
this error; the port's early degree-zero roots [] is wrong for these rings.

Native oracle must execute bundled generic roots with a proxy overriding word
factor using polynomial_factor_dispatch._bundled_word_factor. Installed Sage10.3
word zero.factor() returns unit0 and would incorrectly yield [] in the fallback.
K._roots_univariate_polynomial may run installed body for prime nonzero (same
source behavior), but use bundled method for zero/composite confidence. Root
outputs must compare returned parent kind/coefficient class, not only values.

Port Polynomial.roots() currently has generic Array<[C,number]> return type.
New IntegerMod-specific overload may need Polynomial<IntegerMod & RingElement>
(as monic uses Integer intersections) returning legacy FiniteFieldElement pairs.
Use actual GF factory to match root parent construction; keep other C calls typed.
Document and test the exported overload in LLM.md; never cast changed return types
away without acknowledging this existing generic API mismatch.

Other root-order callers remain open:
- hyperelliptic field_ops.compare_elements uses cmp hook then integer_representation;
  wrong for PARI extension elements. sort_roots_like_sage already negates each root,
  so prime-field zero ordering is correct there. New gen2.cmp_universal can supply
  the extension comparison after reading/replicating native sorted behavior.
- elliptic weierstrass_morphism._sortRootsLikeSage deduplicates and compares negated
  roots via numeric-or-string _elementKey; extension string order is suspect.
  Both root-without-multiplicities callers need source review (finite-field method
  computes gcd(f,x^q-x), then factors) and permanent comparative regressions.

Finite irreducibility still open:
- word Sage cached -> FLINT nmod_poly_is_irreducible; bundled is_irreducible.c uses
  small-root trial rejection and Shoup coarse DDF (not Rabin). Not yet ported.
  Existing flint nmod_poly_factor/factor_distinct_deg.ts may supply dependencies.
- GF2 Sage cached -> NTL GF2X_IterIrredTest (already exported via ntl index alias).
- large-prime polynomial_modn_dense_ntl has no own override: inspect inheritance,
  likely generic cached factor-based method.
- extension polynomial_zz_pex.pyx:320 uncached algorithm fast_when_false/fast_when_true/
  probabilistic, iter=1 -> NTL ZZ_pEX tests. Port local Rabin/no options remains open.

Follow-up source review: bundled integer_mod_ring.py:1550 has the full
_lift_residue_field_root static method. Nonzero derivative uses precision-doubling
Newton; singular roots lift one p-adic digit at a time, root-major then digit-major
ordering. Composite root lists combine factored_order prime powers with CRT_basis
and Cartesian-product order (last component fastest). Degree-one nonunit case
divides by gcd(N,a), solves over N/g, then appends increasing increments N/g.
These are output-sized loops, not enumeration of the entire composite ring.
Zero roots must enumerate because every residue is output. Implement these exact
branches, not a cap or naive search. Existing arith.CRT_basis and factor provide
dependencies; IntegerModRing currently lacks factored_order.

GF factory currently returns fresh objects (documented global factory limitation).
Cache field() per IntegerModRing instance without adding a private class member
(WeakMap avoids changing nominal diagnostics unnecessarily). Default Zmod roots
return the cached GF parent; False roots return the input IntegerModRing parent.

Prepared (scratch only): generate.py/cases.json has 3,990 root cases;
modular_polynomial_roots.py is an unexecuted oracle draft. It loads actual bundled
modular roots/lift methods and generic roots fallback, routes recursive polynomials
through proxy methods and finite-field calls through the existing bundled oracle.
Validate the proxy before accepting any outputs; framework exceptions are not
comparative expectations. _path currently uses Path.cwd() and must be made
file-relative if moved to tests/property/python.
