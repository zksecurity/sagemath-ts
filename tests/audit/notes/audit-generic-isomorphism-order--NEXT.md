After committing 24.52.2, continue the root-caller audit. No subagents.

Open source bugs:
- ell_generic._poly_roots still calls default roots and sorts positive values;
  its remaining users are is_x_coord, lift_x (char2), montgomery_model. The
  duplicate isomorphism body is now gone. Read each source caller before fixing.
- ell_generic.compareFieldElements still compares .value or strings; affects
  remaining point/model root selection. Weierstrass now has a correct private
  comparator; hyperelliptic field_ops.compare_elements still uses integer repr.
- ell_generic._square_roots searches only prime-subfield constants for a
  nonresidue; even-degree extension fields fail or run O(q). QQ also fails the
  _field_order assumptions. Source scalar hooks should be delegated, not replaced
  by polynomial factorization or another custom curve-level Tonelli loop.
- Extension FiniteFieldElement (finite_field_extension.ts) has no is_square/sqrt.
  Bundled element_pari_ffelt.pyx1001/1032 calls FF_issquare / FF_issquareall.
  Neither wrapper nor Fq square-root dependency currently exists in the port.
- IMPORTANT: FF_issquareall(x,&root) calls FF_ispower(x,2,&root), NOT FF_sqrt.
  FF_ispower (FF.c918) preserves zero, then dispatches FpXQ_sqrtn/Flxq_sqrtn/
  F2xq_sqrtn. Root choice and random state may differ from plain FF_sqrt.
  Flxq_sqrtn (Flx.c3657+) for exponent2, degree>1 selects F=2 via the factor/
  degree condition and calls Flxq_sqrtn_spec_pre -> Flxq_sqrtl_spec_pre. For large
  p, FpX.c2860+ follows analogous FpXQ_sqrtn_spec. Need read these full specialized
  algorithms and dependencies before coding. Binary F2x.c1700 uses
  gen_Shanks_sqrtn; gcd(2,2^d-1)=1 gives inverse-exponent powering, no random draw.
  Existing ellsea.Fp2_sqrt is a bespoke Tonelli loop, NOT the same exact source
  route. Do not delegate to it merely because the name mentions square roots.
- Hyperelliptic odd_degree_model chooses first DISTINCT polynomial root;
  jacobian_morphism.cantor_reduction chooses default roots()[0][0]. Both port
  callers redundantly sort via incorrect field_ops comparator. lift_x char2
  wants distinct roots then positive y sorting; mixed-parent lifting still open.
- Optional field arguments/common-parent injections for generic elliptic APIs,
  general unknown-field ordering, and global Sage factory identity remain outside
  current repaired profiles, documented in DEVIATIONS.
- Finite-field irreducibility remains open: see
  ../audit-finite-irreducibility/NOTES.md and original root-order caller NOTES.md.

- Parent equality source: finite_field_base.pyx103 explicitly compares identity.
  Current adapter tests canonical default proof/check metadata plus explicit
  modulus/name/generator and maps the port's legacy prime implementations to
  native modn. Nondefault proof/check/backend cache-key metadata are not yet
  represented by legacy parent definitions; do not claim all factory identity.
