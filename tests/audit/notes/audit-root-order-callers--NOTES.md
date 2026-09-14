Next stage: read-only source review; no production changes yet.

hyperelliptic/field_ops.ts compare_elements at338 uses cmp hook then
integer_representation, incorrectly claiming PARI extensions compare that way.
Delegate extension comparison to committed parigp/gen2.cmp_universal; preserve
QQ cmp and prime scalar paths. sort_roots_like_sage negates elements then sorts.
sqrt_all_of sorts positive roots; lift_x/points_at_infinity/points choose first.
Native comparator same-field profiles should use impl=pari_ffelt, explicit
moduli and integer-representation input encoding, as existing root oracle.
End-to-end: native hyperelliptic_generic lift_x sorts ys, odd_degree_model takes
first distinct polynomial root. cantor_reduction default roots()[0][0] at infinity
is re-sorted incorrectly by port jacobian_morphism.ts99-109 after the now-correct
Polynomial.roots result. Inspect multiplicities before removing duplicate sort.

weierstrass_morphism.ts has bigger bugs than sort order:
- _elementKey + _sortRootsLikeSage sorts negated extension strings, not native
  PARI field order.
- All square/cube/nth-root helpers use characteristic p in place of field order q.
  Inverse Frobenius in char2/3 returns x instead of x^(p^(degree-1)).
- char2 j=0 _isomorphisms line288 solves s^2+a3E*s+c instead of the required
  quartic s^4+a3E*s+c (the preceding comment actually says quartic!).
- _roots_char2 and _cubic_roots_char3 enumerate only prime subfield constants,
  missing extension roots.
Bundled weierstrass_morphism.py317-366 constructs actual polynomials and calls
roots(multiplicities=False) for every branch. Replace local root algorithms with
that exact delegation, using current Polynomial.roots(False); remove obsolete
helpers/import discrete_log when unused. Keep function names/helper scope as
appropriate; no new arbitrary enumeration caps. Read entire _isomorphisms before
editing, including rational/infinite field support profile.

Comparative curve suggestions:
GF(25) j=0 has extension sixth roots; GF(49) j=1728 has extension fourth roots.
GF(4), GF(16) char2 j=0 exercise cubic/quartic/quadratic roots and 24 automorphisms.
GF(9), GF(27) char3 exercise quartic/cubic roots. Add transformed isomorphic pairs
with nonconstant u,r,s,t and nonisomorphic equal-j pairs, ordered raw _isomorphisms
lists plus public first-isomorphism behavior, input-parent identity and validity.
Native algorithm/oracle should execute bundled _isomorphisms with bundled finite
root method rather than relying on installed10.3 wrapper differences. Use an
external timeout on old TS runs because its wrong group-order loops can hang.
Do not run generators/compilers alongside active timed tests.

Follow-up source findings during 24.52.1 validation (NOT repaired in this batch):
- ell_generic._isomorphisms_unsorted duplicates the full generator and routes
  through _poly_roots, which still calls default roots and sorts POSITIVE values.
  This makes public isomorphism_to differ from the repaired Weierstrass constructor.
  Replace duplicate body with [..._isomorphisms(this, other)] after source review;
  existing tuple-return API is documented. _poly_roots is also used by is_x_coord,
  lift_x and montgomery_model, so audit all before changing the helper.
- ell_generic.compareFieldElements uses numeric .value or string order. QQ fractions
  and PARI extension elements need their native comparators. _compute_isomorphisms
  has the correct source sorting-key formula but the wrong scalar comparator.
  lift_x also uses it (and already explicitly sorts the positive y values).
- WeierstrassIsomorphism._comparison_impl only supports eq/ne; bundled source
  supports all six rich comparisons, including native domain/codomain comparison
  before the special (i,min(v,-v),j,v) tuple key. Wrong operand types currently
  ValueError; bundled returns NotImplemented. Read ell_generic._richcmp_ before
  implementing different-curve comparisons; preserve port operator conventions.
- Keep 24.52.1 focused on the completed generator root fixes; begin the above as
  the next caller checkpoint, then the previously noted hyperelliptic callers.
- Further ell_generic root-caller bug: _square_roots searches non-residues by
  K.__call__(integer), so it only searches the prime subfield of an extension.
  In even extension degree every nonzero prime-subfield element is a square;
  the search fails (or spends O(q) steps) instead of returning square roots.
  QQ _field_order assumptions also fail. Read actual scalar sqrt/is_square
  implementations and delegate those source hooks, with live lift_x,
  is_x_coord, Montgomery and kernel-polynomial comparisons.
- Native EllipticCurve_generic inherits WithEqualityById (fast_methods.pyx133),
  not an ainvariant lexicographic comparator. Different-curve ordering in
  morphism._comparison_impl follows richcmp_not_equal and may be TypeError.
