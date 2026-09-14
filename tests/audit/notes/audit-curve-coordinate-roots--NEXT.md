After finishing/committing 24.54.1, continue source audit; no subagents.

Current profile repaired: same-base-field is_x_coord/lift_x, untwisted
Montgomery root mode/selection, scalar dependency traces and native y sorting,
extension coefficient parentheses. 5816 shared comparisons all pass, full
10196-row advanced area matches. General input coercion is the next batch.

Bundled ell_generic.py 710: is_x_coord ALWAYS calls K(x); catch only TypeError
and replace with 'x must be coercible into the base ring of the curve'. The
current port only coerces numbers/bigints. Foreign objects bypass the parent.

Bundled lift_x 913: py_scalar_to_element (sage.structure.coerce, NOT element),
L=x.parent(); K.coerce_map_from(L) first, otherwise L.coerce_map_from(K) and
E.change_ring(L). A constructor conversion is NOT a canonical coercion. QQ
coordinates may convert numerically to finite fields yet fail canonical
coercion in lift_x. If L is an extension of K, even a constant x promotes the
curve. Exact error: 'Unable to construct a point with x in {} over {}'.
General extend=True builds L.fraction_field().extension([-f,b,1],names='y').
Current TS signature lacks extend and assumes same-base-field coordinates.

Read finite_field_base.pyx 1243-1320 _coerce_map_from_: int/ZZ yes;
IntegerModRing with characteristic divisible by p yes; same finite parent yes;
same-characteristic degree-one source yes; otherwise degree divisibility AND
both parents having _prefix (algebraic closure compatibility), not arbitrary
explicit extension moduli. Prime-specific _coerce_map_from_ is at
finite_field_prime_modn.py95+. Field conversion __call__ can differ from these
maps. Existing _same_base_ring tests explicit names/moduli but native global
factory metadata is documented as open. Rational has no runtime parent property;
QQ.characteristic is a method whereas finite parents expose a bigint property.

TS base_extend/change_ring exists at ell_generic.ts1355/1382 (line numbers
before future shifts); converts coefficients via a.value ?? a, which can bypass
source conversion. Audit alongside parent promotion and preserve native errors.

Build native coordinate-parent matrix: actual bundled is_x_coord/lift_x bodies,
QQ, ZZ/Python scalars, different prime implementations, GF2 and explicit small
extensions (same/different names/moduli/degrees). Parent promotion and exact
errors must be compared. Existing local scalar roots now work; installed native
root signs cancel under source lift_x sort, but use bundled polynomial dispatcher
for consistency. Do not claim scalar-state comparison for this separate profile
unless explicitly binding bundled C root state as in curve_coordinate_roots.py.

Other pending root callers: finite torsion_points is a port-only alias, not a
native finite method; currently enumerates only characteristic p constants and
only y=0/1 in char2. Source native rational torsion_points is ell_rational_field
4126 and delegates to rational torsion subgroup. Current finite property tests
compare cardinality (adapter), not a nonexistent Sage finite torsion_points.
Hyperelliptic field_ops.compare_elements uses integer representation, differing
from PARI cmp_universal. odd_degree_model uses first DISTINCT root;
cantor_reduction uses first DEFAULT root. Both have wrong extra positive sort.
Scalar sqrt_all_unsorted still has redundant predicate/binary powering.

Field norm remains charpoly-based in Sage (element_base.pyx632/pari_ffelt982),
NOT FF_norm directly. Trace delegates FF_trace. Finite irreducibility notes in
../audit-finite-irreducibility/NOTES.md remain open.
