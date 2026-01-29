# SageMath Tests Catalog for sagemath-ts Port

This document catalogs the test sources in SageMath that are relevant to our TypeScript port,
focusing on cryptography-related functionality.

## Summary Statistics

| Module | Files | Doctest Lines | Priority |
|--------|-------|---------------|----------|
| sage/arith/ | 9 | ~1,558 | HIGH |
| sage/rings/integer.pyx | 1 | ~1,243 | HIGH |
| sage/rings/integer_ring.pyx | 1 | ~240 | MEDIUM |
| sage/rings/rational.pyx | 1 | ~599 | MEDIUM |
| sage/rings/rational_field.py | 1 | ~219 | LOW |
| sage/rings/finite_rings/ | 23 | ~3,964 | HIGH |
| sage/rings/polynomial/ | 83 | ~22,090 | HIGH |
| sage/schemes/elliptic_curves/ | 54 | ~13,136 | CRITICAL |

**Total: ~43,049 doctest lines across 173 files**

---

## 1. Arithmetic Module (sage/arith/)

### Source Files
- `reference/sage/src/sage/arith/misc.py` (~1,235 doctests)
- `reference/sage/src/sage/arith/functions.pyx` (~42 doctests)
- `reference/sage/src/sage/arith/power.pyx` (~16 doctests)
- `reference/sage/src/sage/arith/srange.pyx` (~78 doctests)
- `reference/sage/src/sage/arith/multi_modular.pyx` (~134 doctests)
- `reference/sage/src/sage/arith/long.pxd` (~42 doctests)
- `reference/sage/src/sage/arith/rational_reconstruction.pyx` (~2 doctests)
- `reference/sage/src/sage/arith/numerical_approx.pyx` (~4 doctests)

### Key Functions Tested

#### Priority: CRITICAL (Crypto-essential)
| Function | Description | We Have? |
|----------|-------------|----------|
| `gcd(a, b)` | Greatest common divisor | YES |
| `lcm(a, b)` | Least common multiple | YES |
| `xgcd(a, b)` | Extended GCD | YES |
| `inverse_mod(a, m)` | Modular inverse | YES |
| `power_mod(a, n, m)` | Modular exponentiation | YES |
| `is_prime(n)` | Primality test | YES |
| `factor(n)` | Integer factorization | YES |
| `sqrt_mod(a, p)` | Modular square root | YES |
| `crt(r1, r2, m1, m2)` | Chinese Remainder Theorem | YES |

#### Priority: HIGH (Important)
| Function | Description | We Have? |
|----------|-------------|----------|
| `euler_phi(n)` | Euler's totient | YES |
| `next_prime(n)` | Next prime | YES |
| `previous_prime(n)` | Previous prime | YES |
| `prime_range(start, stop)` | Primes in range | YES |
| `kronecker_symbol(a, n)` | Kronecker/Legendre symbol | YES |
| `isqrt(n)` | Integer square root | YES |
| `is_square(n)` | Perfect square test | YES |
| `is_squarefree(n)` | Squarefree test | YES |
| `valuation(n, p)` | p-adic valuation | YES |
| `divisors(n)` | List divisors | YES |
| `sigma(n, k)` | Sum of divisors | YES |
| `moebius(n)` | Mobius function | YES |

#### Priority: MEDIUM (Nice to have)
| Function | Description | We Have? |
|----------|-------------|----------|
| `factorial(n)` | Factorial | NO |
| `binomial(n, k)` | Binomial coefficient | NO |
| `bernoulli(n)` | Bernoulli numbers | NO |
| `radical(n)` | Product of prime divisors | YES |
| `squarefree_part(n)` | Squarefree part | YES |
| `prime_factors(n)` | List prime factors | YES |

### Our Test Coverage
File: `packages/sagemath-ts/src/arith/misc.test.ts`
- 27 test suites
- ~150 individual test cases
- Covers all CRITICAL and HIGH priority functions

---

## 2. Integer Ring (sage/rings/integer*.pyx)

### Source Files
- `reference/sage/src/sage/rings/integer.pyx` (~1,243 doctests)
- `reference/sage/src/sage/rings/integer_ring.pyx` (~240 doctests)

### Key Operations Tested

#### Priority: CRITICAL
| Operation | Description | We Have? |
|-----------|-------------|----------|
| `__add__`, `__sub__`, `__mul__` | Basic arithmetic | YES (via BigInt) |
| `__pow__` | Exponentiation | YES |
| `__neg__` | Negation | YES |
| `__mod__` | Modular reduction | YES |
| `__floordiv__` | Integer division | YES |
| `gcd()` | GCD method | YES |
| `lcm()` | LCM method | YES |
| `inverse_mod(n)` | Modular inverse | YES |
| `sqrt()` | Square root | PARTIAL |
| `sqrtrem()` | Square root with remainder | NO |
| `factor()` | Factorization | YES |
| `is_prime()` | Primality test | YES |

#### Priority: HIGH
| Operation | Description | We Have? |
|-----------|-------------|----------|
| `digits(base)` | Base representation | NO |
| `nbits()` | Bit length | YES (via BigInt) |
| `bit_length()` | Bit length | YES |
| `trial_division()` | Trial division | PARTIAL |
| `is_perfect_power()` | Perfect power test | NO |
| `nth_root()` | nth root | NO |
| `jacobi_symbol()` | Jacobi symbol | YES |

### Our Test Coverage
- Integer operations tested indirectly via `arith/misc.test.ts`
- BigInt native operations provide most functionality
- No dedicated Integer class tests (using native BigInt)

---

## 3. Rational Numbers (sage/rings/rational*.pyx)

### Source Files
- `reference/sage/src/sage/rings/rational.pyx` (~599 doctests)
- `reference/sage/src/sage/rings/rational_field.py` (~219 doctests)

### Key Operations Tested

#### Priority: MEDIUM
| Operation | Description | We Have? |
|-----------|-------------|----------|
| `__add__`, `__sub__`, `__mul__`, `__truediv__` | Arithmetic | YES |
| `numerator()`, `denominator()` | Components | YES |
| `floor()`, `ceil()` | Rounding | NO |
| `is_integral()` | Integer test | NO |
| `height()` | Height of rational | NO |
| `valuation(p)` | p-adic valuation | NO |

### Our Test Coverage
File: `packages/sagemath-ts/src/rings/rational.test.ts`
- Basic Rational class tests
- Arithmetic operations

---

## 4. Finite Fields (sage/rings/finite_rings/)

### Source Files (23 files, ~3,964 doctests)

#### Core Files
- `finite_field_base.pyx` (~369 doctests)
- `finite_field_constructor.py` (~149 doctests)
- `finite_field_prime_modn.py` (~53 doctests)
- `integer_mod.pyx` (~615 doctests)
- `integer_mod_ring.py` (~455 doctests)
- `element_base.pyx` (~257 doctests)

#### Extension Field Files
- `element_givaro.pyx` (~248 doctests)
- `element_ntl_gf2e.pyx` (~179 doctests)
- `element_pari_ffelt.pyx` (~314 doctests)
- `finite_field_givaro.py` (~122 doctests)
- `finite_field_ntl_gf2e.py` (~61 doctests)
- `finite_field_pari_ffelt.py` (~43 doctests)

### Key Operations Tested

#### Priority: CRITICAL
| Operation | Description | We Have? |
|-----------|-------------|----------|
| `GF(p)` constructor | Prime field creation | YES |
| `GF(p^n)` constructor | Extension field creation | YES |
| `__add__`, `__sub__`, `__mul__` | Field arithmetic | YES |
| `__truediv__`, `inv()` | Division/inverse | YES |
| `__pow__` | Exponentiation | YES |
| `is_square()` | Quadratic residue test | YES |
| `sqrt()` | Square root in field | YES |
| `multiplicative_order()` | Element order | YES |
| `characteristic()` | Field characteristic | YES |
| `cardinality()` | Field size | YES |

#### Priority: HIGH
| Operation | Description | We Have? |
|-----------|-------------|----------|
| `multiplicative_generator()` | Generator | YES |
| `quadratic_non_residue()` | Non-residue | YES |
| `frobenius()` | Frobenius endomorphism | PARTIAL |
| `trace()` | Field trace | NO |
| `norm()` | Field norm | NO |
| `minimal_polynomial()` | Minimal poly | NO |
| `random_element()` | Random element | YES |
| `list()` | All elements | YES |

### Our Test Coverage
File: `packages/sagemath-ts/src/rings/finite_rings/finite_field_prime.test.ts` (~566 lines)
- IntegerMod tests
- IntegerModRing (Zmod) tests
- GF(p) prime field tests
- Field axiom verification
- Quadratic residue tests
- Large prime field tests
- Multiplicative order tests

File: `packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.test.ts`
- Extension field GF(p^n) tests
- Conway polynomial tests

File: `packages/sagemath-ts/src/rings/finite_rings/tower_field.test.ts`
- Tower field construction tests

---

## 5. Polynomial Rings (sage/rings/polynomial/)

### Source Files (83 files, ~22,090 doctests)

#### Core Files
- `polynomial_element.pyx` (~3,049 doctests) - CRITICAL
- `polynomial_ring.py` (~608 doctests)
- `polynomial_ring_constructor.py` (~180 doctests)

#### Specialized Files
- `polynomial_gf2x.pyx` (~43 doctests) - Binary polynomials
- `polynomial_zmod_flint.pyx` (~184 doctests) - Mod p polynomials
- `polynomial_integer_dense_flint.pyx` (~308 doctests)
- `polynomial_rational_flint.pyx` (~428 doctests)
- `multi_polynomial_element.py` (~593 doctests)
- `multi_polynomial_libsingular.pyx` (~1,414 doctests)

### Key Operations Tested

#### Priority: CRITICAL
| Operation | Description | We Have? |
|-----------|-------------|----------|
| Polynomial constructor | Create polynomials | YES |
| `__add__`, `__sub__`, `__mul__` | Arithmetic | YES |
| `__pow__` | Exponentiation | YES |
| `degree()` | Polynomial degree | YES |
| `__getitem__` | Coefficient access | YES |
| `leading_coefficient()` | Leading coeff | YES |
| `__call__` | Evaluation | YES |

#### Priority: HIGH
| Operation | Description | We Have? |
|-----------|-------------|----------|
| `gcd(f, g)` | Polynomial GCD | PARTIAL |
| `divmod(f, g)` | Division with remainder | PARTIAL |
| `roots()` | Find roots | NO |
| `factor()` | Factorization | NO |
| `is_irreducible()` | Irreducibility test | NO |
| `derivative()` | Differentiation | NO |
| `resultant(g)` | Resultant | NO |
| `xgcd(g)` | Extended GCD | NO |

#### Priority: MEDIUM
| Operation | Description | We Have? |
|-----------|-------------|----------|
| `monic()` | Monic polynomial | NO |
| `content()` | Content | NO |
| `primitive_part()` | Primitive part | NO |
| `reverse()` | Reverse coefficients | NO |
| `shift()` | Coefficient shift | NO |
| `truncate()` | Truncation | NO |

### Our Test Coverage
File: `packages/sagemath-ts/src/rings/polynomial/polynomial.test.ts`
- Basic polynomial operations over integers
- Polynomial over GF(2)
- Quotient rings
- ~200 lines of tests

---

## 6. Elliptic Curves (sage/schemes/elliptic_curves/)

### Source Files (54 files, ~13,136 doctests)

#### Core Files - CRITICAL
- `ell_generic.py` (~626 doctests) - Base curve class
- `ell_point.py` (~1,157 doctests) - Point operations
- `ell_finite_field.py` (~567 doctests) - Curves over Fq
- `constructor.py` (~251 doctests) - Curve construction

#### Isogeny Files - HIGH
- `ell_curve_isogeny.py` (~911 doctests)
- `isogeny_small_degree.py` (~326 doctests)
- `hom_velusqrt.py` (~266 doctests) - Velu sqrt isogeny
- `hom_composite.py` (~259 doctests)
- `hom.py` (~421 doctests)

#### Number Theory Files - MEDIUM
- `ell_rational_field.py` (~888 doctests)
- `ell_number_field.py` (~869 doctests)
- `heegner.py` (~1,143 doctests)
- `period_lattice.py` (~506 doctests)

### Key Operations Tested

#### Priority: CRITICAL (Core EC Operations)
| Operation | Description | We Have? |
|-----------|-------------|----------|
| `EllipticCurve([a, b])` | Curve construction | YES |
| Point addition | P + Q | YES |
| Point doubling | 2P | YES |
| Point negation | -P | YES |
| Scalar multiplication | nP | YES |
| `is_on_curve(P)` | Point validation | YES |
| `point(x, y)` | Create point | YES |
| `zero()` / Identity | Point at infinity | YES |
| `discriminant()` | Curve discriminant | YES |
| `j_invariant()` | j-invariant | YES |

#### Priority: CRITICAL (Finite Field Operations)
| Operation | Description | We Have? |
|-----------|-------------|----------|
| `cardinality()` | Curve order #E(Fq) | YES |
| `order()` | Point order | YES |
| `random_point()` | Random point | YES |
| `points()` | All points | YES |
| `gens()` | Group generators | YES |
| `is_x_coord(x)` | Valid x-coord test | YES |
| `lift_x(x)` | Find point with x | YES |
| `trace_of_frobenius()` | Frobenius trace | YES |
| `set_order(n)` | Set known order | YES |

#### Priority: HIGH (Crypto-relevant)
| Operation | Description | We Have? |
|-----------|-------------|----------|
| `abelian_group()` | Group structure | PARTIAL |
| `weil_pairing(P, Q, n)` | Weil pairing | NO |
| `tate_pairing(P, Q, n)` | Tate pairing | NO |
| `isogeny(kernel)` | Compute isogeny | NO |
| `isogenies_prime_degree(l)` | l-isogenies | NO |
| `division_polynomial(n)` | Division poly | NO |
| `frobenius_endomorphism()` | Frobenius | NO |

#### Priority: MEDIUM
| Operation | Description | We Have? |
|-----------|-------------|----------|
| `a_invariants()` | Weierstrass coeffs | YES |
| `b_invariants()` | b-invariants | NO |
| `c_invariants()` | c-invariants | NO |
| `short_weierstrass_model()` | Transform to y^2=x^3+ax+b | NO |
| `is_supersingular()` | Supersingular test | NO |
| `is_ordinary()` | Ordinary test | NO |
| `quadratic_twist()` | Twist curve | NO |
| `twists()` | All twists | NO |

### Our Test Coverage

File: `packages/sagemath-ts/src/schemes/elliptic_curves/ell_finite_field.test.ts` (~677 lines)
- Point operations (add, double, neg, mul)
- Point order computation
- Curve construction and validation
- Cardinality (point counting)
- Points enumeration
- Generators
- set_order
- Trace of Frobenius
- Edge cases (order 2 points, small fields)

File: `packages/sagemath-ts/src/schemes/elliptic_curves/ell_generic.test.ts`
- Generic curve tests

File: `packages/sagemath-ts/src/schemes/elliptic_curves/ell_point.test.ts`
- Point-specific tests

File: `packages/parigp-ts/src/elliptic/*.test.ts` (multiple files)
- PARI/GP compatible elliptic curve tests
- Comparison tests with PARI output

File: `tests/property/typescript/elliptic.test.ts` (~598 lines)
- Property tests comparing with PARI/GP
- 29 test cases covering all basic operations

---

## 7. Priority Ranking for Missing Tests

### CRITICAL - Must Add
1. **Polynomial GCD/division over finite fields** - Essential for crypto
2. **Weil/Tate pairings** - Essential for pairing-based crypto
3. **Isogeny computation** - Essential for isogeny-based crypto
4. **Division polynomials** - Needed for torsion point computation

### HIGH - Should Add
1. **Polynomial factorization over finite fields**
2. **Irreducibility testing**
3. **Frobenius endomorphism**
4. **Abelian group structure computation**
5. **Polynomial root finding**

### MEDIUM - Nice to Have
1. **Curve transformations (short Weierstrass, etc.)**
2. **Quadratic twists**
3. **Supersingular/ordinary classification**
4. **Field trace and norm**
5. **Binomial coefficients and factorials**

---

## 8. Test Sources We Should Replicate

### From sage/arith/misc.py
- All `gcd`, `lcm`, `xgcd` examples (DONE)
- All `is_prime` examples including Carmichael numbers (DONE)
- All `factor` examples (DONE)
- All `power_mod`, `inverse_mod` examples (DONE)
- All `sqrt_mod` examples including Tonelli-Shanks (DONE)
- All `crt` examples (DONE)

### From sage/rings/finite_rings/
- All `GF(p)` construction examples (DONE)
- All field arithmetic examples (DONE)
- All `is_square`, `sqrt` examples (DONE)
- Extension field examples (DONE)
- Conway polynomial lookup (DONE)

### From sage/rings/polynomial/polynomial_element.pyx
- Polynomial arithmetic examples (PARTIAL)
- Division/remainder examples (TODO)
- GCD examples over various rings (TODO)
- Root finding examples (TODO)
- Factorization examples (TODO)

### From sage/schemes/elliptic_curves/
- All point arithmetic examples (DONE)
- All `cardinality` examples (DONE)
- All `order` examples (DONE)
- Isogeny examples (TODO)
- Pairing examples (TODO)
- Division polynomial examples (TODO)

---

## 9. Recommendations

### Immediate Actions
1. Add polynomial GCD tests over finite fields
2. Add polynomial division tests
3. Add more edge case tests for EC point counting

### Short-term
1. Implement and test Weil pairing
2. Implement and test basic isogeny computation
3. Add polynomial factorization over GF(p)

### Long-term
1. Full isogeny support (Velu, sqrt-Velu)
2. Tate pairing implementation
3. Complex multiplication support
4. Number field support for curves

---

## Appendix: File Locations

### SageMath Reference Sources
```
reference/sage/src/sage/arith/
reference/sage/src/sage/rings/integer.pyx
reference/sage/src/sage/rings/integer_ring.pyx
reference/sage/src/sage/rings/rational.pyx
reference/sage/src/sage/rings/rational_field.py
reference/sage/src/sage/rings/finite_rings/
reference/sage/src/sage/rings/polynomial/
reference/sage/src/sage/schemes/elliptic_curves/
```

### Our Test Files
```
packages/sagemath-ts/src/arith/misc.test.ts
packages/sagemath-ts/src/rings/finite_rings/finite_field_prime.test.ts
packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.test.ts
packages/sagemath-ts/src/rings/finite_rings/tower_field.test.ts
packages/sagemath-ts/src/rings/polynomial/polynomial.test.ts
packages/sagemath-ts/src/rings/rational.test.ts
packages/sagemath-ts/src/schemes/elliptic_curves/ell_finite_field.test.ts
packages/sagemath-ts/src/schemes/elliptic_curves/ell_generic.test.ts
packages/sagemath-ts/src/schemes/elliptic_curves/ell_point.test.ts
packages/parigp-ts/src/elliptic/*.test.ts
tests/property/typescript/elliptic.test.ts
```
