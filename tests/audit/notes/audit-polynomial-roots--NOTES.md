Next audit candidate; source read only, no code/test/case changes yet.

polynomial_element.ts roots() at2581 calls local findIntegerRoots / findRationalRoots.
Both enumerate divisors of constant/leading coefficients using integer factorization
(getDivisorsBigInt now delegates factorInteger, so the old sqrt trial scan was fixed,
but the unnecessary integer factorization remains). Linear polynomials with difficult
semiprime coefficients are useful watched regressions. Huge content should be discarded.
No failures have been run yet; do not claim these probes as confirmed bugs.

Original sage/rings/integer_ring.pyx:1233–1460 _roots_univariate_polynomial:
- zero raises ValueError('roots of 0 are not defined').
- default dense degree<=100; sparse >100 (based on initial degree, before valuation).
- dense: strip nonunit content, factor polynomial, _roots_from_factorization.
- sparse: append zero first with valuation; shift; if now totally dense use factor
  and append to zero. Otherwise strip content, find exponent gaps greater than
  running max coefficient nbits. GCD of contiguous blocks. No gap -> factor.
  If a gap exists, compute +/-1 multiplicities via source's sparse derivative loop,
  append +1 then -1, then roots of abs>1 from g.factor. Do NOT globally sort this
  result: current sortRootsSageOrder would reorder zeros and +/-1 incorrectly.
- current port only exposes roots() with multiplicities; ring/algorithm kwargs
  remain documented omissions. Do not gratuitously expand signatures.

Original polynomial_integer_dense_flint.pyx:1619 factor(): strip content, NTL when
(degree<30 || degree>300), PARI otherwise. Current native drivers exist at
ntl-ts/ZZXFactoring.factor and parigp-ts/QX_factor.ZX_factor. Source Factorization
sorts by degree, multiplicity and factor coefficients, not raw backend order.

Original QQ rational_field.py:1606 _factor_univariate_polynomial always invokes PARI,
normalizes factors monic and puts leading coefficient in the unit. Generic
polynomial_element.pyx roots around9110 calls factor then _roots_from_factorization;
latter filters degree1 and coerces -g[0]/g[1] to base ring.
Current local factorIntegerPolynomial (~5937) is a separate FLINT-like factor driver
used by multiple exported factor/irreducibility methods, so don't replace shared
calls blindly: QQ factor and ZZ factor have distinct backend routing.
DEVIATIONS Polynomial Roots and Factorization still calls Sage's factor FLINT and
claims full correctness based on old tests; update current route statements if repaired.

Public roots() currently raises ValueError('roots of zero polynomial are not defined')
for all rings. Confirm the exact zero behavior for QQ from original source/live Sage:
generic factor() raises ArithmeticError('factorization of 0 is not defined'), which
roots catches and converts to NotImplementedError('root finding for this polynomial
not implemented'). This is not yet executed/confirmed.

Useful regression shape: ZZ x^101*(x-2)*(x+1)^2 triggers sparse dispatch and must keep
zero first; QQ uses generic factor order and may differ. Cover 99/100/101 degree,
valuation/dense-after-shift, no-gap vs gap, gaps leading to gcd1 and nonconstant gcd,
+/-1 multiplicities and nonmonic rational linear factors; source currently finds
other roots through polynomial factorization without factoring coefficient integers.

Backend nfrootsQ is NOT a replacement for this root API: it doesn't supply
multiplicities and Sage's methods use the above factor routing. The new Galois
flag-zero delegation cannot simply use flag4: full nfroots/nfsqff/get_nfsqff_data is
missing (~2000 C lines + nfinit/Trager/relative-polynomial dependencies). Keep tracked.
