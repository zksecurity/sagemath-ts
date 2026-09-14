Next audit: public Polynomial.factor() and is_irreducible() over ZZ/QQ, then finite roots.
Read only / scratch probes so far. Do not claim these paths fixed.

Bundled source:
- polynomial_integer_dense_flint.pyx factor1619: c=self.content();g=self//c;
  c.factor() * g._factor_ntl() when degree<30 or>300, else g._factor_pari().
  Must factor content BEFORE polynomial factor, including RNG ordering.
- Same file _factor_ntl1571 uses native ZZX_factor and Factorization sorting.
- polynomial_element.pyx _factor_pari_helper: primitive integral factors with
  corrected content/unit, fields normalize each factor monic; final Factorization
  sorts by degree, exponent, polynomial coeffs. Constants from integer content
  must follow this sorting, not numerical prime order alone.
- rational_field.py1606 _factor_univariate_polynomial always PARI, normalize
  monic, unit=f.leading_coefficient().
- polynomial_rational_flint.pyx1613 is_irreducible: cached; length<2false,
  length2true, otherwise primitive integer numerator.is_irreducible(). Thus
  QQ factor uses PARI but QQ is_irreducible uses ZZ's degree-based NTL/PARI.
- polynomial_element.pyx10115 is_irreducible cached: zero/unit false,
  constants delegate coefficient.is_irreducible(), else self.factor count.
  Port is_irreducible currently has no cache. Const ZZ primality proof uses
  Integer.is_irreducible -> is_prime(proof=True); port uses probable prime
  currently, existing documented primality limitation.

Observed native in probe.log (Sage10.3; bundled matches source paths):
- ZZ[x](0).factor() raises ZeroDivisionError('division by zero'), unlike port's
  shared ArithmeticError. Source _floordiv_scalar checks zero at1231.
- ZZ[x](-12).factor(): unit-1 and factors2^2,3. Port earlyconstant returns[-12,1],
  leaving integer content unfactored. Const1native empty non-unit factors unit1.
- QQ[x](0).factor() ArithmeticError('factorization of 0 is not defined').
- GF(7)[x](0).factor(): empty factors, unit0; GF(2)zero ArithmeticError.
  Check bundled nmod_poly_linkage factor_helper to rule out version differences.
  These finite zero facts also affect roots(): GF(7)zero roots [] in native,
  GF(2)zero roots ArithmeticError. Port global roots zeroValueError wrong.
- For normal GF(q), finite_field_base.pyx2182 always f.factor() when
  multiplicitiesTrue, then _roots_from_factorization. Remove small-field
  enumeration only after fixing factor zero and verifying native ordering.

Current port:
- factor() degree0 earlyreturn[[self,1]] occurs BEFORE ring dispatch.
- _factorOverIntegers uses local factorIntegerPolynomial before factoringcontent.
- _factorOverRationals local factorIntegerPolynomial too.
- ZZ is_irreducible primitive-check and same localdriver;
  QQ is_irreducible uses publicfactor (wrongnativebackend fordegree<30/>300).
- local factorIntegerPolynomial plus private FLINT helpers must remain for
  _zz_factor_internal and polynomial_factorization.test.ts. Replace onlypublic
  callers / relabel docs until all privatebackendsareactuallyported.
- Current new roots helper duplicates dense factor dispatch but returns just
  linearroots; consider sharing a signed-content nativefactor helper withnew
  publicpaths, without losing content-beforefactor timing or zero semantics.
- rootcurrentNTLdefaultprivatecontext documented; lowlevelexplicitstatecompare
  exists separately. PARI live statecompare canprovepublicdriverselection.

Oracle next:
- independent Python/TS area functions for factor and irreducibility;
  normalize Sage unit into degree0 extra factor only whenunit!=1, preserving
  documentedarrayadaptation and Sage sorted factor order. Don't sortawaybug.
- cover signed constants,0,1,content,linear/nonmonic/repeated, irreducibles,
  boundaries29/30/31 and299/300/301 (use easy sparse Eisenstein/repeatedpowers).
- Include RNG state: original bundled QX_factor native C oracle op5 for
  PARI factor branch, setrand1/getrand captures; installed Sage PARI2.15 may
  differ from completebundledPARI2.18. Useappropriatecorrectversionoracle.
- Irreducibility cache should compare repeated call on same immutablepoly;
  new object must advance native state again. f.coeffs readonly property
  points to mutablearray; existing classimmutabilitycontract auditneededbefore
  addingcache blindly (RationalPolynomial separateclassalreadyhascache).

Important correction after bundled source review:
- polynomial_zmod_flint.pyx:817 NOW explicitly rejects zero with ArithmeticError,
  unlike installed Sage10.3 (GF7zeroFactorunit0). The installed GF7zero behavior
  is NOT a port factor bug; do not change factor to accept zero. Finite rootzero
  stillneeds ArithmeticError through bundled finite_field_base factor route.
- cached_method decorators verified on generic is_irreducible10114 and QQ1612.
- Polynomial.constructor clones coeffs; readonlyarray is thedeclared publicAPI.
  _new_c createsinstances viaObject.create, so anyfuture cache cannot require
  constructor-initializedprivate slots. Optionalordinaryfield/WeakMap handlesit.

Candidate implementation ordering (not yet applied):
- Dispatch ZZ factor before constant earlyreturn. ZeroZZ must throw
  ZeroDivisionError('division by zero'); QQ/finite keep bundledArithmeticError.
- ForZZfactor compute positivegcd c; factorInteger(c) BEFORE nativepolynomial
  factor; preserveinputleadingsign inprimitive passedto NTL/ZX_factor;
  signunit−1 pluscontentprimes and polynomialfactors sortwith existing comparator.
- QQfactor QX_factor(clearDenominators(this)); existingmonicunitconversionreused.
- Genericnonzero constantfield1 should [] underexistingnontrivialunit-only
  adaptation; valuesnon1return[[constant,1]]. Audit tests beforechanging.
- QQirreducibility primitiveintegerfactorization selectsNTL/PARI bydegree,
  notQQfactor'salwaysPARI route. Belowlength2 shortcutsnative.
- ZZirreducibility sourcefactors evennonprimitivepolynomials; a fastprimitive
  false shortcut changesbackend/RNG/errorbehavior. Keepcontentprimefactors
  whencheckingnativefactorizationlength; ignoreonlyunit±1,notallconstants.
- Allgenericandword is_irreducible methodsarecached upstream. Optionalfieldor
  WeakMap (includingfalse) canworkwith Object.create-based_new_c. Test cached
  repeat versus freshobjectPARIstate; nevercachethrownexceptions.
