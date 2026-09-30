# Deviations from SageMath

> **📋 This is a living document.** Anyone porting SageMath functionality MUST update this file
> when their implementation differs from SageMath behaviour.

**A deviation is a documented failure to be faithful to SageMath.** That framing decides where an
entry belongs and whether it belongs here at all.

## How to Use This Document

The document has two top-level parts, and the difference between them is the whole point:

- **[Part I — Accepted Deviations](#part-i--accepted-deviations)** are differences we intend to
  keep. They are limited to three kinds: adaptations forced by the language or type system,
  deliberate improvements in exactness, and upstream bugs we decline to reproduce. Each one has a
  rationale that would still be true after every dependency is ported.
- **[Part II — Open Fidelity Gaps](#part-ii--open-fidelity-gaps)** are differences that **should
  be closed by porting**. They are work items, not accepted deviations. Every entry names the
  vendored upstream file that implements the missing behaviour and an estimate of the effort to
  close it. "Out of scope for this pass" is a scheduling statement, not a rationale, so those
  entries live here.

**Routing rules (from CLAUDE.md):**

| Kind of difference | Goes in |
|--------------------|---------|
| Outputs, errors or observable shapes differ from SageMath | **this file** |
| Type mappings, naming and architectural conventions | `DESIGN.md` |
| What is implemented and what is not, per module | `SCOPE.md` |
| What changed in which version | `CHANGELOG.md` |

This file is a register of the *current* state. It is not a changelog: an entry describing
something that has been fixed must be deleted, not annotated with "resolved in 0.0.x".

**For library users:** review the relevant Part I entry before using a module, and check Part II
for whether your use case is blocked. Deviations are linked from function docstrings via
`@see Deviation:`.

**Finding functions affected by deviations:**

```bash
grep -r "@see Deviation" packages/
grep -r "SAGE_NOT_IMPLEMENTED" packages/   # the authoritative list of unimplemented paths
```

**Not recorded here:** encoding conventions internal to the property-test harness
(`tests/property/`) — flat row-major matrix transport, integer-indexed operand tables, exceptions
rendered as compared strings rather than propagated. Those are test-design choices documented in
the area modules themselves; they assert *more* than a bare comparison would and describe no
library behaviour.

## Table of Contents

### [Part I — Accepted Deviations](#part-i--accepted-deviations)

1. [Language and Type-System Adaptations](#language-and-type-system-adaptations)
2. [Return Shapes, Keyword Arguments and Signature Adaptations](#return-shapes-keyword-arguments-and-signature-adaptations)
3. [Port-Only APIs With No SageMath Counterpart](#port-only-apis-with-no-sagemath-counterpart)
4. [Infinity Representation](#infinity-representation)
5. [Exact Arithmetic Where SageMath Uses Floating Point](#exact-arithmetic-where-sagemath-uses-floating-point)
6. [No Arbitrary-Precision Floating Point](#no-arbitrary-precision-floating-point)
7. [Upstream Behaviour Deliberately Not Reproduced](#upstream-behaviour-deliberately-not-reproduced)
8. [Honest Failure Instead of Silent Approximation](#honest-failure-instead-of-silent-approximation)
9. [Bounded Search Budgets and Measured Thresholds](#bounded-search-budgets-and-measured-thresholds)
10. [Random State and Seeding](#random-state-and-seeding)
11. [Vendored SageMath 10.9.beta4 vs Installed 10.3](#vendored-sagemath-109beta4-vs-installed-103)
12. [Number Fields — Exactness-Driven Divergences](#number-fields--exactness-driven-divergences)
13. [Polynomial Roots and Factorization](#polynomial-roots-and-factorization)
14. [Finite Field Constructors and Display](#finite-field-constructors-and-display)
15. [Generic Group API and DLP](#generic-group-api-and-dlp)
16. [Matrix Module Algorithm Substitutions](#matrix-module-algorithm-substitutions)
17. [Matrix Special Constructors](#matrix-special-constructors)
18. [Lattice Algorithms — CVP, Voronoi Cells and LLL Representatives](#lattice-algorithms--cvp-voronoi-cells-and-lll-representatives)
19. [Free Module Exactness and Coordinate Types](#free-module-exactness-and-coordinate-types)
20. [Binary Quadratic Forms](#binary-quadratic-forms)
21. [Quadratic Forms (sage.quadratic_forms)](#quadratic-forms-sagequadratic_forms)
22. [Elliptic Curves and Isogenies](#elliptic-curves-and-isogenies)
23. [Hyperelliptic Curves and Jacobians](#hyperelliptic-curves-and-jacobians)
24. [Quaternion Algebras](#quaternion-algebras)
25. [Function Fields](#function-fields)
26. [Power Series, Laurent Series and Multivariate Series](#power-series-laurent-series-and-multivariate-series)
27. [Coding Theory](#coding-theory)
28. [Crypto Module](#crypto-module)
29. [Discrete Gaussian Samplers](#discrete-gaussian-samplers)
30. [ZK Sumcheck and Multilinear Extensions](#zk-sumcheck-and-multilinear-extensions)
31. [GF(2) Matrix PNG Functions](#gf2-matrix-png-functions)
32. [PARI Integer Factorization (parigp-ts)](#pari-integer-factorization-parigp-ts)
33. [PARI Elliptic Curve Algorithms (parigp-ts)](#pari-elliptic-curve-algorithms-parigp-ts)
34. [ntl-ts GF2X Representation](#ntl-ts-gf2x-representation)
35. [Newly Ported Upstream Modules — Residual Divergences](#newly-ported-upstream-modules--residual-divergences)
36. [Exact Return Types and Numeric Backends](#exact-return-types-and-numeric-backends)

### [Part II — Open Fidelity Gaps](#part-ii--open-fidelity-gaps)

37. [Number Field Class Groups, Units and Galois Closure](#number-field-class-groups-units-and-galois-closure)
38. [Number-Field Kernel Not Delegated to parigp-ts](#number-field-kernel-not-delegated-to-parigp-ts)
39. [Quadratic Class Numbers Not Delegated to Buchquad](#quadratic-class-numbers-not-delegated-to-buchquad)
40. [PARI/NTL Routines Duplicated or Ported In Place](#parintl-routines-duplicated-or-ported-in-place)
41. [parigp-ts Elliptic Curves — SEA Dispatch and Isogeny Stubs](#parigp-ts-elliptic-curves--sea-dispatch-and-isogeny-stubs)
42. [ntl-ts GF2X Factoring Stubs](#ntl-ts-gf2x-factoring-stubs)
43. [Elliptic Curves over Q and Number Fields](#elliptic-curves-over-q-and-number-fields)
44. [p-adic Precision Models, Extension Fields and L-Series](#p-adic-precision-models-extension-fields-and-l-series)
45. [Arithmetic Functions Not Delegated to PARI/FLINT](#arithmetic-functions-not-delegated-to-pariflint)
46. [Finite Fields — Conway Table and Constructor Algorithms](#finite-fields--conway-table-and-constructor-algorithms)
47. [Polynomials — Printing, Factor Shape, Term Orders and Base Rings](#polynomials--printing-factor-shape-term-orders-and-base-rings)
48. [Matrices — the matrix() Constructor](#matrices--the-matrix-constructor)
49. [Lattices — Exact SVP Rank Cap](#lattices--exact-svp-rank-cap)
50. [Parents Are Not Unique](#parents-are-not-unique)
51. [Real and Complex Precision and Rounding](#real-and-complex-precision-and-rounding)
52. [Power Series — V(0)](#power-series--v0)
53. [Coding and Crypto — Permissive Where Upstream Raises](#coding-and-crypto--permissive-where-upstream-raises)
54. [Hyperelliptic — Frobenius Polynomial Algorithms](#hyperelliptic--frobenius-polynomial-algorithms)
55. [Quaternion Algebras — Base Rings Other Than QQ](#quaternion-algebras--base-rings-other-than-qq)

56. [Polynomial String Parsing](#polynomial-string-parsing)
57. [Modular Integer Coercion and Factories](#modular-integer-coercion-and-factories)
58. [Extension Arithmetic and PARI Quotient Kernels](#extension-arithmetic-and-pari-quotient-kernels)
59. [Template for New Deviations](#template-for-new-deviations)

---

# Part I — Accepted Deviations

Differences we intend to keep. Each is forced by the language or type system, is a deliberate
improvement in exactness, or is a refusal to reproduce an upstream defect.

---

## Language and Type-System Adaptations

TypeScript has no `__call__`, no operator overloading, no keyword arguments, no duck typing and no
Python exception hierarchy. `DESIGN.md` is the normative reference for how each concept is mapped
(see its sections "Parameter Types"/"Return Types", "Options Objects", "Ring `__call__`
Pattern", "No Operator Overloading" and "Library Mapping" — named rather than line-numbered,
because line references into an edited document rot silently). Only the **observable**
consequences are recorded here.

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| JavaScript `number` as an integer argument | Python `int` is arbitrary precision, so coercion is always safe | `IntegerLike = bigint \| Integer` (`types/coercion.ts:28`); a `number` throws the native `TypeError` (`:43-49`) |
| Callable parents | `F(3)` | `F.__call__(3n)`. `GF(p)`/`Zmod(n)` are factory *functions* that build the ring; `ZZ`/`QQ` are singleton ring objects, not constructors |
| Operators | `a + b`, `a * b`, `a ** n` | `a.add(b)`, `a.mul(b)`, `a.pow(n)`; coercion is explicit |
| Keyword arguments | Python kwargs | Trailing options objects (see [Return Shapes](#return-shapes-keyword-arguments-and-signature-adaptations) for the cases that change an observable shape) |
| Generic ring elements | Duck typing | Explicit TypeScript interfaces; a value that satisfies no interface throws where Sage would coerce |
| Exception hierarchy | `ZeroDivisionError`/`OverflowError` derive from `ArithmeticError`; `NotImplementedError` from `RuntimeError` | `errors.ts` is **flat**: every class extends `Error` directly, and `TypeError` is the JS built-in re-exported (`errors.ts:14`). A Sage `except ArithmeticError` has no faithful `catch` here |
| Affected modules | All | All |

### Rationale

1. **Precision safety.** A source literal like `9007199254740993` is already corrupted to
   `…992` before any runtime check could see it, so accepting `number` cannot be made safe. For a
   library targeting cryptography, silent truncation is the one failure mode to rule out
   structurally.
2. **The other rows are not choices.** TypeScript cannot express `__call__`, operator
   overloading or Python's exception tree; there is no implementation that would be more faithful.
3. **Static types catch at compile time** what Sage's duck typing catches at run time.

### Trade-offs

- Call sites are more verbose; integer literals need the `n` suffix.
- Catching a *category* of error (all arithmetic errors) is impossible; callers must name each
  class.
- **Acknowledged exceptions to the `number` rule.** Three public integer-valued APIs still accept a
  raw `number` and are therefore vulnerable to exactly the hazard above:
  `finite_field_constructor.ts:58` `GF(order: bigint | number)` (an order above 2^53 is corrupted
  before `BigInt()` sees it), `conway_polynomials.ts:278`
  `conway_polynomial(p: number, n: number): number[]`, and `convolution.ts:129`
  `find_primitive_root(n: number, …)`. These should be widened to `IntegerLike`.

### Behavioral Impact

`TypeError` where SageMath would coerce. Error *names* match SageMath; error *classes* are JS
subclasses of `Error` with no inheritance between them. Mathematical results are unaffected.

### Bigint to Number Conversion

When internal code must convert `bigint` to `number` (floating-point math, array indexing), use
`toSafeNumber()`, which throws `RangeError` above the safe-integer range. Never use `Number(x)`.

---

## Return Shapes, Keyword Arguments and Signature Adaptations

Python keyword arguments map to options objects and Python tuples/dicts to TypeScript equivalents
(`DESIGN.md`). The cases below change an **observable shape** or argument position and are
therefore registered individually.

| Function | SageMath | sagemath-ts |
|----------|----------|-------------|
| `gcd()`, `xgcd()` | `Integer`, Python tuple | `bigint`, `[bigint, bigint, bigint]` |
| `factor()` (integers) | `Factorization` object with a separate `unit()` | `Array<[bigint, bigint]>`; see [Polynomials](#polynomials--printing-factor-shape-term-orders-and-base-rings) for the unit-handling gap |
| `power_series_ring.pade(m, n)` | `u / v` in `Frac(R[z])` | A `PadeApproximant` object with `numerator()`, `denominator()`, `power_series(prec)`, and a `toString()` matching Sage's fraction-field repr |
| `formal_group.group_law()` | Element of `PowerSeriesRing(R, 2, 't1,t2')` | A `BivariatePowerSeries` class; `toString` reproduces Sage's textual form |
| `padic square_root`, `nth_root` | `square_root(extend=True, all=False, algorithm=None)`, `nth_root(n, all=False)` | Options object `{extend, all}`, plus `square_root_all()` / `nth_root_all()` for well-typed list access. The `algorithm` keyword is not offered — there is one implementation (Sage's `'sage'` path) and PARI's returns the same value |
| `matrix indefinite_factorization` | `(L, vector(L.base_ring(), d))` | `[Matrix, R[]]` |
| `matrix krylov_basis` / `krylov_kernel_basis` | `(M, shifts, degrees, output_rows=True, var=None, basis_algorithm=None)` | Same, but the argument is named `variable` (`var` is reserved in JavaScript) and accepts a string or a `PolynomialRing`; `degrees` also accepts a single integer |
| `QQ(...)` (`RationalField.__call__`) | `rational.pyx:591-704` accepts `Rational`, `Integer`, `int`, `float`, `str`, a length-1 list, `fractions.Fraction`, anything with `_rational_()`, and `(n, d)` | A port of that method with Sage's exact `TypeError("unable to convert {!r} to a rational")` and `ValueError('denominator must not be 0')`. Two additions: the two-argument quotient spelling is kept and accepts rational operands (Sage passes `d` as a *base*, and string bases are unsupported here); single-argument pairs require integral components as in Sage, and `{numer, denom}` objects are accepted alongside Sage's `{numerator, denominator}` |
| `matrix density()` | Exact rational (`2/3`) | Exact `Rational` |
| `matrix norm(p)` | `RDF` element | JS `number` — which is what RDF is |
| `matrix is_similar(…, transformation=true)` | `(False, None)`; the returned `T` satisfies `A == T^-1·B·T` | `[false, null]`; return type widened to `boolean \| [boolean, Matrix<R> \| null]`. Sage's `T` convention is preserved |
| `matrix frobenius_form(2)` | Two elements of `MatrixSpace(QQ, n)` | `[Rational[][], Rational[][]]` — `B` is genuinely rational and there is no rational-matrix class. The empty matrix yields `[[], []]` |
| `binary_qf reduced_form(algorithm=…)` | `reduced_form(transformation=False, algorithm='default')` | `reduced_form({ transformation?, algorithm? })`, identical semantics and error messages |
| `binary_qf solve_integer(n)` | `n` may be an `Integer` **or** a `Factorization` | `solve_integer(n, { factorization? })`, mirroring PARI's own `[n, factor(n)]` input. There is no `Factorization` class |
| `Integer.multiplicative_order(n)` | No such signature — Sage's takes no argument | The two-argument form was removed; use `Mod(a, n).multiplicative_order()` |
| `Integer.log(b)` | Exact `Integer` when `b^k == self`, else a real/symbolic logarithm | Always `floor(log_b(self))` (alias for `exact_log`) — there is no symbolic ring. Test exactness with `b ** result === self` |
| `MPolynomial.monomial_coefficients()` | `dict` keyed by `ETuple` | `Map` keyed by the canonical comma-joined exponent string (`exponentToKey()`), padded to `parent.ngens_value`; `keyToExponent()` recovers the tuple. JS `Map` keys compare by identity, so array keys would never collide |
| `MPolynomial.coefficient({x: 1, y: 1})` | Dict keyed by generator objects | `Record` keyed by variable **name**. Sage's list-with-nulls and monomial forms are byte-for-byte |
| `crypto.gen_lattice(quotient=…)` | Symbolic expression or univariate polynomial | `IntegerLike[]`, the ascending coefficient list (`x^4-1` is `[-1n,0n,0n,0n,1n]`). Degree checks and Sage's messages preserved; additionally throws `TypeError('quotient must be monic')` |
| `Matrix.subs` | Transmits to the entry's `subs` | Falls back to the entry's `evaluate` (how univariate `Polynomial` exposes substitution here) before raising |
| `Matrix.denominator` | LCM of entry denominators, an element of the denominator ring | Reads `denominator` whether method or getter, takes the LCM in `bigint`, returns `ring(lcm)` — so a QQ matrix gives the *rational* 30 |
| `zk sumcheckVerify` | Round count from `len(poly.args())`; `degree_checks` defaults to `None` | `(proof, claimedSum, polyEvaluator, field, numVars, options?: {degreeCheck?})`. `numVars` is **required** and a proof with a different number of rounds is rejected — a stateless verifier cannot recover the round count without trusting the prover. `degreeCheck` defaults to no check |
| `zk binaryToInt` | No counterpart (Python ints are unbounded) | Returns `bigint`; throws `ValueError('bits must contain only 0 and 1 (got X)')` |
| `NumberFieldIdeal.norm()` | `QQ` | `Rational`. `prime_above(p, {degree})` returns one ideal, `primes_above(p)` the list, `decomposition(p)` the `[P, e]` pairs with `bigint` exponents |
| `ClassGroup.exponent()` | Invariants in decreasing order (`d_{i+1} \| d_i`), e.g. `(38, 2)` | `max(invariants)` — correct under either ordering convention |
| `convolution()` | Any commutative ring in which multiplication by two is injective | `bigint` natively, plus ring elements exposing `.parent` with `__call__` and `div`/`inv`; anything else raises `NotImplementedError`. The port's `RingElement` carries no `parent`/`zero`, so `R(x/M)` cannot be expressed generically |
| `QuotientRing` non-invertible element | `ArithmeticError` on the Singular path; `ZeroDivisionError(f'element {self} of quotient polynomial ring not invertible')` on the fallback | Always the fallback message (there is no Singular path) |
| `IntegerMod.log(b, order)` | `log(self, b=None, order=None, *, check=False)`; `order` consulted only when `check=True` | Same |
| `Matrix<Rational>`, `PolynomialRing<Rational>`, `PowerSeriesRing<Integer>` façades | n/a | `RationalMatrix`/`MatrixQQ`, `RationalPolynomial`/`PolynomialRingQQ` (`quatalg`, `quadratic_forms`) are structurally typed **views** over the repo's real `Matrix`/`PolynomialRing` classes. `Matrix<Rational>` does not typecheck under `strict` because `Rational.add` accepts `Rational \| IntegerLike` rather than exactly `this`; the runtime objects are unchanged, so all of `matrix_operations` applies |

### Rationale

1. **DESIGN.md mapping** — keyword arguments become options objects, tuples become arrays, `None`
   becomes `null`.
2. **No polymorphic parents** — TypeScript cannot express "element of the base ring", so accessors
   split by concrete type.
3. **JS `Map` semantics** — object and array keys compare by identity, so canonical string keys are
   the only correct sparse representation.
4. **Soundness** — `sumcheckVerify`'s required `numVars` closes a forgery: a short proof verified
   nothing.

### Trade-offs

- Source-breaking for callers of `Integer.multiplicative_order(n)`, `sumcheckVerify`,
  `binaryToInt`, `indefinite_factorization`, `krylov_*` and `NumberFieldIdeal.norm()`.
- Dictionary-keyed `coefficient()` calls need a string instead of a generator.
- `PadeApproximant`, `BivariatePowerSeries` and the `Rational` façades should be replaced when real
  `Frac(R[z])` / `MultiPowerSeriesRing` types and a `RingElement`-compatible `Rational` exist.

### Behavioral Impact

Values are SageMath's; only shapes and argument positions differ. Every `toString()` listed above
reproduces SageMath's printed form. (For the polynomial `toString` that does **not**, see
[Polynomials](#polynomials--printing-factor-shape-term-orders-and-base-rings).)

---

## Port-Only APIs With No SageMath Counterpart

Symbols that exist in this port and **not** in SageMath. Listed so their presence is never mistaken
for a SageMath contract, and so `@see Deviation:` docstrings have a target.

| Symbol | Module | Notes |
|--------|--------|-------|
| `polynomial_commitment.ts` (870 lines: `compute_quotient`, `batch_quotient`, `barycentric_weights`, `fri_fold`, `split_poly`, `generate_powers`, …) | `zk/` | KZG/FRI helpers. There is no `sage/rings/polynomial/polynomial_commitment.py`. `rings/polynomial/index.ts` retains a re-export block of the identical 38 values + 2 types, marked in-file as backwards compatibility only. `package.json` has **no `./zk` subpath export**, so a direct `@sagemath-ts/sagemath-ts/zk` import is not possible — the symbols are reachable from the package root (`export * as zk`) and via `./rings` |
| `src/zk/` (`sumcheck.ts`, `multilinear.ts`) | `zk/` | Ports of `reference/sage_blueprints/`, not of SageMath |
| `sparseMultilinearExtension([i])` | `zk/multilinear.ts` | Returns the selector `eq(i, x)`; the blueprint short-circuits to the constant `R(i)` — an inconsistency with every other input shape, and its branch still contains leftover `print` debugging |
| `estimateBKZBlockSize`, `bkzRootHermiteFactor`, `qaryLattice`, `qaryDualLattice` | `modules/free_module_integer.ts` | Port-invented. The BKZ estimator interpolates Gama-Nguyen/Chen values below `beta = 40` and uses the standard asymptotic formula above; it is a **heuristic** and says so in its docstring |
| `fold(codeword, challenge, domain?)`, `fold_domain(length)` | `coding/reed_solomon.ts` | FRI folding has no SageMath counterpart |
| `error_correction_capability()` | `coding/goppa_code.ts` | SageMath's `GoppaCode` registers no decoder at all. Returns the *decoder radius*: `deg(g)` in characteristic 2 with squarefree `g` (Patterson), else `floor(deg(g)/2)`. `distance_bound()` is exactly Sage's `1 + deg(g)` |
| `sampleExact()`, `samplesExact()`, `isIntegral`, `basisExact`, `cNumeric()` | `stats/distributions/discrete_gaussian_lattice.ts` | See [Discrete Gaussian Samplers](#discrete-gaussian-samplers) |
| `Integer.nth_root_mod` | `rings/integer_ring.ts` | Sage exposes this only as `Mod(a, p).nth_root(n)`; there is no `Integer` method, so there is no doctest to match. When `a == 1` and `gcd(n, q-1) > 1` we return 1 where Sage's `_nth_root_common` returns a *primitive* gcd-th root of unity (`Mod(1,11).nth_root(5)` is **9** in SageMath 10.3, deterministically across runs). Every returned value satisfies `r^n == a (mod p)` |
| `sqrt_mod(a, m)` | `arith/misc.ts` | `sage/arith/misc.py:2274` is a **commented-out stub** (`# def sqrt_mod(a, m):`), so there is no upstream function. Returns `null` for a non-residue |
| `Ti`, `GF2n(n)` | `rings/finite_rings/tower_field.ts` | Binary tower fields as used in Binius. There is no `sage/rings/finite_rings/tower_field.py`, so this module has no SageMath counterpart. `GF2n` accepts only powers of two — a property of the Binius construction, not a restriction on `GF(2^n)`, which `FiniteField(2^n)` builds for any `n` |
| `RandState.random()`, `python_random(seed)` | `misc/randstate.ts` | `random()` is a `@deprecated` alias for `c_rand_double()` retained for three call sites; `python_random`'s `seed` parameter lets a caller get a self-contained stream instead of inheriting the global state (see [Random State](#random-state-and-seeding)) |
| `PythonRandom` | `misc/randstate.ts` | A port of CPython's `random.Random`, which Sage reaches through `randstate.python_random()` rather than exposing as a class |
| `NumberFieldElement.is_integral_unit()` | `rings/number_field/` | Sage overloads `is_unit()` on the parent's type; without a distinct element class for orders the two behaviours cannot share one name |
| `NumberFieldElement.minimal_polynomial()`, `absolute_trace()`, `relative_trace()` | `rings/number_field/` | Existing convenience aliases for Sage's `minpoly()` and absolute-field `trace()`. Retained for caller compatibility; they add method names without changing the mapped values. Relative-field trace semantics are not provided by these absolute-field aliases |
| `NumberFieldElement.numerator()` outside quadratic fields | `rings/number_field/` | Extends Sage's quadratic numerator convenience to all supported absolute fields as `a * a.denominator()`. Retained to support uniform callers; general Sage elements require that explicit expression. Mapped values agree, but method availability differs |
| `QuadraticField.d` | `rings/number_field/` | The squarefree part. `.D` is Sage's `D`. Retained only so existing callers do not silently get the wrong number; **treat as deprecated** |
| `realQuadraticFundamentalUnit(K)` | `rings/number_field/unit_group.ts` | Sage exposes the fundamental unit only through `K.units()`. The quadratic case bypasses `bnfinit` entirely (as PARI's own `bnfinit` does), so it deserves a directly testable entry point |
| `from_png_data`, `to_png_data` | `matrix/matrix_mod2.ts` | The data-only substitutes for libgd file I/O; see [GF(2) Matrix PNG Functions](#gf2-matrix-png-functions) |
| `x_list(prec)`, `y_list(prec)` | `schemes/elliptic_curves/formal_group.ts` | Coefficient accessors added while `LaurentSeriesElement` had no arithmetic. That gap is closed, so these are redundant conveniences retained because callers exist; they return the coefficients from valuation −2 resp. −3 |
| `pAdicEisensteinQuadraticExtension`, `pAdicEisensteinQuadraticElement` | `schemes/elliptic_curves/padic_lseries.ts` | Sage builds this with `K.extension(f, names='alpha')` in `sage/rings/padics/`. Ours lives beside its only consumer because `pAdicExtension` is a shell with no element type |
| `Frobenius_filter` | `schemes/elliptic_curves/isogeny_class.ts` | Sage's is in `gal_reps_number_field.py`, which is not ported; `isogeny_class.ts` is its only caller |
| `FractionFieldElement`, `tensorProductVector` | `modules/free_module.ts` | Elements of `Frac(R)` for a Euclidean base ring (e.g. `QQ(x)`), and the elementary tensor |
| `change_ring(matrix, ring)`, `pivots` | `matrix/matrix_operations.ts`, `matrix/matrix_decompositions.ts` | Sage's `change_ring` is a `Matrix` method in `matrix0.pyx`; it landed in `matrix_operations.ts` for file-ownership reasons. Both are re-exported from `matrix/index.ts` |
| `PadeApproximant`, `BivariatePowerSeries` | `rings/power_series_ring.ts`, `schemes/elliptic_curves/formal_group.ts` | Stand-ins for `Frac(R[z])` and `PowerSeriesRing(R, 2)` |
| `GF2X.rep()` | `ntl-ts` | Exposes the packed bigint so `sagemath-ts` can convert cheaply |
| `monomial_coefficient(exponentTuple)` | `rings/polynomial/multi_polynomial_element.ts` | A pure **superset** of Sage's signature (which takes a monomial with the same parent, and is fully supported) |
| `Ratio`, `isRatio` | `parigp-ts/src/elliptic/init.ts` | Exact rational j-invariants; not re-exported from the package root |
| `Fp_ellcard_Shanks`, `Fp_elldivpol(l, a4, a6, p)` | `parigp-ts/src/elliptic/` | `static` in PARI / no PARI counterpart over `F_p` (PARI's SEA uses modular polynomials rather than `psi_l`). Exported so the BSGS branch and the division-polynomial recursion are testable against exhaustive oracles |
| `qfrep(Q, bound)` | `stats/distributions/discrete_gaussian_lattice.ts` | A thin adapter over `parigp-ts`'s `qfrep0`, kept only so the module's own tests can drive it; callers should import `qfrep0` directly |
| `matkermod_basis`, `zm_from_rows`, `zm_to_rows` | `parigp-ts/src/matkermod.ts` | Row-major wrappers around PARI's column-major `ZM` layout, so a caller cannot silently transpose the meaning of "the kernel basis" |
| `field_ops.ts` (`field_embedding`, `sort_roots_like_sage`, `constant_field_*`) | `schemes/hyperelliptic_curves/`, `rings/function_field/` | Sage's parents/elements make `K.characteristic()`, `a.is_square()`, `iter(K)` uniformly available; this port's field classes share no interface, so the structural dispatch is collected in one file per module rather than duck-typed at every call site |
| `ZZLattice`, `INFINITE_PLACE_QQ` | `algebras/quatalg/quaternion_algebra.ts` | Stand-ins for `sage.modules.free_module` ZZ-spans of rational vectors, and for the ring morphism `QQ -> RR` that represents the infinite place |

### Rationale

1. **ZK-specific functionality** is the project's stated focus and has no SageMath equivalent.
2. **Testability** — a few upstream-`static` symbols are exported so oracle tests can reach them.
3. **Migration safety** — a few (`.d`, `RandState.random`, `x_list`) are retained deliberately so
   that existing callers fail loudly rather than silently.

### Trade-offs

- Readers of the mirrored file layout can mistake these for SageMath APIs.
- Some are marked deprecated but not yet removed.

### Behavioral Impact

None on any SageMath-named function. Each symbol above carries an `@see Deviation:` docstring.

---

## Infinity Representation

SageMath has `sage.rings.infinity` with genuine `PlusInfinity` elements. This port has no infinity
ring, so three sentinels are in use depending on the module's return type.

| Method | SageMath returns | sagemath-ts returns |
|--------|-----------------|---------------------|
| `Integer.multiplicative_order()` (for `n` other than ±1) | `infinity` | The string literal `'Infinity'` (typed `bigint \| 'Infinity'`) |
| `Integer.valuation(0, p)` / `ord(0, p)` and negative `popcount()` / `hamming_weight()` | `+Infinity` | `'Infinity'`, typed `bigint \| 'Infinity'` |
| `Integer.exact_log(0, b)` / `log(0, b)` | `-Infinity` | `'-Infinity'`; `exact_log` returns `bigint \| '-Infinity'`, `log` returns `Integer \| '-Infinity'` |
| `Rational.valuation()`, `ord()`, `numerator_valuation()`, `val_unit()`, `RationalField.quadratic_defect()` | `infinity` | `'Infinity'` |
| `pAdicGenericElement.valuation()`, `multiplicative_order()`, `additive_order()` | `infinity` | `Number.POSITIVE_INFINITY`, typed `InfiniteOr<bigint> = bigint \| number` |
| `RealNumber`/`ComplexNumber` `multiplicative_order()`, `additive_order()` | `infinity` | `Number.POSITIVE_INFINITY` |
| `FreeModule.cardinality()`, `indexIn()` | `infinity` | `Number.POSITIVE_INFINITY` |
| `FunctionFieldElement.valuation()`, `FunctionFieldIdeal.valuation()` | `+Infinity` | `Number.POSITIVE_INFINITY` (typed `bigint \| number`) |
| `QuotientRing.cardinality()` | `Infinity` | `Number.POSITIVE_INFINITY` |
| `arith.valuation(0, p)` | `+Infinity` | `'Infinity'`, typed `bigint \| 'Infinity'`, delegated to Integer |

### Rationale

1. **`bigint` has no infinite value**, and widening every arithmetic return type to admit one would
   ripple through every caller.
2. **Numeric comparisons work** — JavaScript relational comparisons between `bigint` and `number`
   behave as expected (`Infinity > 0n`), so guards read naturally.
3. **Module-local consistency** — the string sentinel predates this register in `rings/rational*.ts`;
   changing it there is a larger breaking change than the inconsistency costs.

### Trade-offs

- Three sentinels for one concept; callers must know which module they are in.

### Mitigation

Introduce a single `Infinity` singleton (or a branded type) and migrate all sites to it.

### Behavioral Impact

Only the representation of the infinite case differs; all finite values are exactly Sage's.
`Integer(1).multiplicative_order()` is `1n` and `Integer(-1).multiplicative_order()` is `2n`, as in
Sage.

---

## Exact Arithmetic Where SageMath Uses Floating Point

CLAUDE.md forbids floating point where Sage or PARI is exact. In several places SageMath itself
uses a bounded-precision float as a *heuristic* guarded by an exact verification; we replace the
heuristic with exact integer arithmetic. Results are identical within the range where Sage's float
has enough precision, and stay correct beyond it.

| Site | SageMath | sagemath-ts |
|------|----------|-------------|
| `Z_isanypower_101` perfect-power search | `logr_abs`/`mpexp` double precision to guess `y = round(x^(1/p))` with a mod-30011 filter, then an exact `powiu(y,p) == x` check | Every prime exponent `p <= log_103(x)` tested with an exact bigint Newton k-th root (the only float left is the integral loop bound `LOG2_103`) |
| `BinaryQF._reduce_indef`, `_Rho`, `_RhoTau`, `BinaryQF_reduced_representatives` | `D.sqrt(prec=53)` and a floored real quotient | Exact `isqrt(D)`. `floor((sqrt(D)+b)/(2\|c\|)) == floor((isqrt(D)+b)/(2\|c\|))` because the numerator bound is integral, and the `\|c\| >= sqrt(D)` branch boundary yields the same `s` either way (argument spelled out at `binary_qf.ts:293-295`) |
| `Matrix_integer_dense.LLL`, `IntegerLattice` reduction | fpLLL/NTL floating-point Gram-Schmidt | Exact integral LLL (see [Matrix Module Algorithm Substitutions](#matrix-module-algorithm-substitutions)) |
| Free module `coordinates`/`echelonize`/`discriminant`/kernels | Exact over the fraction field | Exact fraction-field layer (see [Free Module Exactness](#free-module-exactness-and-coordinate-types)) |
| `Integer.real_log` and rational logarithmic heights | MPFR | Overflow-safe binary64 approximation from the integer bit length and leading bits; `Math.log` for small values. This is approximate, not exact transcendental arithmetic |
| Hyperelliptic Frobenius precision bounds | `M = 2·binomial(2g,g)·RR(q).sqrt()^g`, then `M.ceil()` and `exact_log(p)` | `M^2` computed exactly in ZZ, integer ceiling of its square root, and `p^(2B) < M^2` compared exactly. All ten doctest values over `GF(37)`, `GF(next_prime(10^9))` and `GF(11)` match |
| `qfrep` Fincke-Pohst enumeration | C `double` Cholesky data with a fudged `BOUND·(1+1e-10)` | Integral rescaling of the Cholesky data; `floor(sqrt((BOUND-y)/v) - z)` becomes `floorDiv(isqrt(…) - Z, d)`, provably the same integer |

### Rationale

1. **CLAUDE.md rule** — "Don't use floating point; use BigInt and rational arithmetic".
2. **Sage's floats are heuristics** — each is followed by an exact check upstream, so replacing the
   guess with an exact computation cannot change the verified answer.
3. **Correctness past 2^53** — the float paths fail silently on large inputs; the exact ones do not.
   This also fixes the class of bugs SageMath tracks in its issue 37635.

### Trade-offs

- The exact search can be slower than the float heuristic (negligible at the sizes involved).
- Exact LLL cannot reproduce fpLLL's choice of representative (documented separately).

### Behavioral Impact

Identical `(base, exponent)` pairs, identical reduced forms, identical logarithms — with the float
paths' silent failures removed.

---

## No Arbitrary-Precision Floating Point

JavaScript has no arbitrary-precision float type. Where SageMath uses MPFR/Arb/LAPACK this port
uses IEEE 754 doubles, **except** where the value is observable and a semantics re-implementation
was warranted (`RealNumberMP`, see [Discrete Gaussian Samplers](#discrete-gaussian-samplers)).

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Real precision | MPFR arbitrary precision | `real_mpfr.ts:478` stores a single `private readonly _value: number`; `RealField(200)(2).sqrt()` returns the 53-bit `1.4142135623730951` |
| Complex results from real operations | Promote to the complex field (`RR(-1).sqrt()` is `1.00000000000000*I`) | Return `Number.NaN` (`real_mpfr.ts:594-690`) |
| Special functions | MPFR/Arb correct rounding | Polynomial/series approximations over the host libm. `gamma(3.5)`, `zeta(2)`, `erf(1)` agree with Sage to double precision |
| SVD | `scipy.linalg.svd` (LAPACK DGESDD/DGESVD) | Jacobi one-sided SVD (`matrix_decompositions_additions.ts:62-66`), default tolerance `1e-14`. On `[[4,0],[3,-5]]` the singular values agree with Sage's `matrix(RDF,…).SVD()` to 1 ulp |
| QR | `scipy.linalg.qr` (LAPACK DGEQRF) | Householder reflections; supports full and reduced |
| LU | LAPACK | Gaussian elimination with partial pivoting; returns `P`, unit-diagonal `L`, `U` with `PA = LU` |
| Affected modules | `rings/real_mpfr.pyx`, `complex_mpfr.pyx`, `matrix_double_dense.pyx`, `matrix2.pyx` | `rings/real_mpfr.ts`, `rings/complex_mpfr.ts`, `matrix/matrix_decompositions_additions.ts` |

### Rationale

1. **Runtime limit** — there is no arbitrary-precision float in JavaScript, and no WebAssembly
   dependency is taken.
2. **Sage is itself inexact here.** `RDF` *is* doubles and `matrix2.pyx`'s `norm(A, 2)` goes through
   `change_ring(CDF)` + SVD, so following it in doubles is faithful, not a shortcut.
3. **Exact routes exist alongside** — `matrix_decompositions.ts` provides exact decompositions over
   exact ring elements; the double-precision helpers are a separate, explicitly-named surface.

### Trade-offs

- ~53-bit precision regardless of the declared `RealField(prec)`; no guaranteed correct rounding.
- Some complex-valued outputs of real operations are not representable and become `NaN`.
- Jacobi SVD converges more slowly than divide-and-conquer LAPACK for large matrices; no complex
  matrix support.
- `matrix_decompositions_additions.ts` has **no SageMath counterpart path** (CLAUDE.md requires
  mirrored paths) and its API is `SVD_double(number[][])` rather than `Matrix.SVD()`.
- Printing is JS `Number#toString`, not the declared precision — registered separately as an open
  gap under [Real and Complex Printing](#real-and-complex-printing-and-precision).

### Behavioral Impact

Numerical results are approximate and diverge from SageMath above 53 bits. Exact
rational/integer matrix operations are unaffected.

---
## Upstream Behaviour Deliberately Not Reproduced

Places where the vendored upstream is itself buggy or crashes, and we deviate on purpose.

| Site | SageMath / PARI / FLINT | sagemath-ts | Why |
|------|-------------------------|-------------|-----|
| PARI `matkermod`'s `m > 2n` shortcut | `bb_hnf.c:1049` computes the kernel as `shallowtrans(matimagemod(shallowtrans(A), d))`. For `A == 0` the image has zero columns, and a `t_MAT` with zero columns carries no row count, so the transpose collapses to 0×0 and the reported basis generates only `{0}` while the real kernel is everything. Reproduced on SageMath 10.3: `matrix(Zmod(4),3,1,[0,0,0]).right_kernel_matrix()` is `[]` | `_right_kernel_matrix_over_integer_mod_ring` calls `matkermod(A, n, wantIm=true)`; PARI's own condition is `!im && m > 2*n`, so requesting the image disables the unsound branch. The same input gives `[[1]]` | An empty basis for a non-trivial kernel is silently wrong. A 2400-case Sage sweep showed the two agree **everywhere except** the all-zero matrices with `m > 2n` (62 cases); below that threshold Sage itself returns the identity, so the fix restores consistency with PARI's own answer in the sound regime |
| PARI 2.18.1's Schoenhage fast reduction for binary quadratic forms | `qfi_red` (`Qfb.c:975-987`), `qfi_redsl2` (`:857-880`), `qfr_red_i` (`:914-940`) and `qfr_redsl2` (`:825-855`) negate `b` before `pqfbred_rec` (which requires non-negative coefficients) and never negate it back. For `b < 0` that returns the reduced form of the **conjugate** class, and `qfi_redsl2` negates only the second *row* of `U`, giving `det U = -1` | The intermediate form is conjugated back and `D·U·D` (`D = diag(1,-1)`) is used instead of a row negation, so `det U == 1` and `Q ∘ U == result`. `qfr_redsl2` additionally falls back to the base case when `a < 0` | The corrected public qfbredsl2 entry preserves the form/transformation identity. The new redimagsl2 compatibility entry used by the LLL wrapper instead preserves native behavior. Both sides of the 9000-bit threshold are compared; the corrected entry was also verified by forcing the threshold to −∞ on 80 sample forms |
| PARI's `matdetmod` documentation | `reference/pari/src/functions/linear_algebra/matdetmod` claims `matdetmod([4,2,3;4,5,6;7,8,9],27) == 9` | Returns `18` | The determinant is −9, i.e. 18 mod 27, and PARI's **own regression output** (`test/32/bbhnf` line 240) records 18. A stale doc block. Sage's `qfbsolve` doc block is stale in the same way — it predates `allsols`' `v >= 0` normalisation and the lexsort at `Qfb.c:1930` — so we follow the source and assert solution *sets* |
| `groups.generic.discrete_log_lambda` with `N = 1` | `k = 0`; step setup makes no random draws and `hash % k` raises `ZeroDivisionError` | **Restored 2026-09-12.** Preserves `k = 0`, the hash callback, native mutation calls and the exception | Exact behavior includes the original failure and random-state effects for a single-value interval |
| `ell_finite_field.twists()` `break` placement | `ell_finite_field.py:1940-1944` puts the `break` at the for-loop level, so only `twists[0]` is ever tested for isomorphism with `self` | **Replicated verbatim**, with a comment citing the line numbers | CLAUDE.md requires behavioural equivalence with the vendored Sage. Implementing the docstring's stated intent would change the returned ordering for `j = 0`/`1728` curves. The returned *set* is complete and pairwise non-isomorphic either way |
| `FunctionFieldIdeal.is_zero()` | Inherits from `Element`, whose `__bool__` falls back to `True` when the parent has no zero — and the parent is the *multiplicative* `IdealMonoid`. So `O.ideal(K(0)).is_zero()` is `False` and the zero ideal prints as `Ideal (0) of …`, never `Zero ideal of …` | Returns `false` unconditionally, with the upstream reason cited in the docstring | This **is** upstream's observable behaviour (executed on 10.3) and is visible through `_repr_` and through `divisor()` never taking its "not defined for zero ideal" branch. Recorded so nobody "fixes" it later |
| Factoring the zero element in a function field | Bundled Sage 10.9 rejects zero in its polynomial factorization backends with `ArithmeticError: factorization of 0 is not defined`; installed Sage 10.3 lacks the word-modular zero guard | Element factorization and finite-ideal factorization/divisors reach the bundled zero guard. The infinite zero ideal still has divisor `Place (1/x)`, from native degree arithmetic | The old always-nonthrowing element behavior is repaired. The comparative oracle applies the bundled polynomial zero guard before using the older runtime |
| `QuadraticForm.__setitem__` | `Q[i,j] = c` leaves `self.__det`, `self.__level` and `self._rational_diagonal_form_and_transformation` **stale** | `Q.set(i, j, c)` clears those caches. `rational_diagonal_form` still deep-copies its result, which is the behaviour the doctest at `local_field_invariants.py:150` pins | Keeping the stale cache would let a mutated form report the old determinant |
| `QuadraticForm.is_rationally_isometric` over ZZ | Raises `AttributeError: 'IntegerRing_class' object has no attribute 'real_embeddings'` | Raises `NotImplementedError` naming the missing support and suggesting `change_ring(QQ)` | There is no correct upstream behaviour to copy — it is a crash. A clear `NotImplementedError` cannot be mistaken for a computed answer |
| `SBox.differential_branch_number` for `n > m` | `sbox.pyx:1435` indexes `_S_list[b]` for `b < 2^n`, so `SBox([0,3])` raises `IndexError` | Returns the value the documented formula gives (3) | The port is right and upstream is wrong; recorded so a later fidelity pass does not reintroduce the crash |
| `hyperelliptic cantor_reduction` root extraction | `r = (x**2 + h[g1]*x - f[2*g1]).roots()[0][0]` raises `IndexError: list index out of range` when the quadratic has no rational root (reproduced on 10.3 with `f = 9x^6+6x^5+7x^3+6x^2+2x+9`, `h = 4x^3+x^2+9x+2` over GF(11)) | Raises `ValueError` naming the quadratic that has no root in the base field. The `deg a < 2g+1` / `deg b < deg a` / divisibility assertions are likewise `ValueError` rather than `AssertionError` | An unhandled edge case in an even-degree model, not intended behaviour. Code catching `IndexError` would not catch this |
| `weak_popov_form(shifts=[])` | Reaches `min([])` and raises `ValueError: min() arg is an empty sequence` (only possible for an `m × 0` matrix) | Returns the (unambiguous) zero-column weak Popov form | Reproducing an incidental Python crash would be a worse API |
| `elementary_matrix(row1 == row2)` with no scale | `special.py:1512-1516`'s four assignments collapse to `elem[r,r] = 1`, i.e. the identity; Sage raises only when a scale is *also* given | Same: identity for a self-swap, `ValueError` for the two cases Sage rejects | An earlier audit claimed Sage raises in general; the vendored source shows otherwise, so the port follows the actual code |
| pynac's `exp(-y)` rewrite in `_normalisation_factor_zz` | For `sigma > 1` pynac rewrites `exp(-y)` as `cosh(y) - sinh(y)`, so Sage evaluates that exponential with catastrophic cancellation (`RealField(53)(exp(-2.42*pi^2))` returns `0`); upstream's own source comment records the wart | Evaluated correctly | No doctest value and none of the 21 pinned oracle values changes; for roughly `1 < sigma < 1.3` our sum keeps correction terms Sage silently drops — in our favour |
| PARI `gen_ellgroup`'s `m` output | `bb_group.c:1035-1043` writes `*pm = g1` then overwrites it with the final iteration's `lcm(s,t)`, after which `gen_ellgens` can fail to terminate | Returns `m = g1` | See [PARI Elliptic Curve Algorithms](#pari-elliptic-curve-algorithms-parigp-ts). PARI 2.15.4 (shipped with Sage 10.3) does not hang, so `g1` reproduces the shipping behaviour |
| PARI's stale `qfrep` GP doc example, `qfrpow`'s double inversion, `qfr5_pow`'s per-word exponent loop, `qfrpowraw`'s distance sign, `qfr5_to_qfr`'s `mplog2(lg(d0))`, `galconj`'s `frobeniusliftall`/`testpermutation` warners | Various | Various | Each is documented in full with its `Qfb.c` / `galconj.c` line numbers under [Newly Ported Upstream Modules](#newly-ported-upstream-modules--residual-divergences) |

### Rationale

Reproducing an upstream crash or an out-of-range return value would propagate the defect into every
consumer, with no offsetting fidelity benefit — nothing can depend on a `ZeroDivisionError` or on a
basis of the wrong length. Where the upstream quirk is merely a *choice* (the `twists` break
placement, the `elementary_matrix` self-swap, the function-field zero ideal), we reproduce it
exactly.

### Trade-offs

- Divergence from the vendored 2.18-dev PARI source in `gen_ellgroup`.
- Code catching Python's incidental `IndexError`/`AttributeError` will not catch our named errors.

### Behavioral Impact

`CRT_basis` now preserves the bundled original partial-prefix behavior, including its
extra coefficient in some non-coprime cases. `discrete_log_lambda` answers where Sage raises.
`ellgenerators` always terminates.

---

## Honest Failure Instead of Silent Approximation

Where a routine cannot produce SageMath's answer, it raises rather than returning a plausible-
looking wrong one. These are listed together because they share one rationale.

| Site | SageMath | sagemath-ts |
|------|----------|-------------|
| `eigenvalues` / `eigenvectors` | `extend=True` by default, working over the algebraic closure | Same default; raises `NotImplementedError` naming the missing algebraic closure when the charpoly does not split over the base ring |
| `MPolynomialRing.__call__` from a univariate polynomial | Converts via `_mpoly_dict_recursive` | `NotImplementedError` naming that routine (it previously fell into the plain-object dictionary branch and produced nonsense) |
| Unknown multivariate term order | `ValueError` | `ValueError("unknown term order 'name'")`, Sage's own message (it previously fell back silently to degrevlex) |
| `groebner_basis` budget exhaustion / non-field base ring | Delegates to Singular | `ArithmeticError` / Sage's `TypeError('Can only reduce polynomials over fields.')` — see [Polynomials](#polynomials--printing-factor-shape-term-orders-and-base-rings) |
| `BCH minimum_distance()` | Inherits `AbstractLinearCode.minimum_distance` (GAP/Guava Brouwer-Zimmermann) | Exact enumeration, cached; `NotImplementedError` once `q^k > 2^17`, pointing at `designed_distance()` |
| `BooleanFunction.truthTable('hex')` with `n < 2` | `ValueError('negative shift count')` | The identical `ValueError` — Sage's own failure, reproduced |
| `booleanHypercube(n)` | `Tuples([0,1], n)`, lazy and unbounded | `ValueError` above `MAX_HYPERCUBE_DIM = 25`; `1 << n` previously wrapped around at n = 31/32 |
| `random_echelonizable_matrix` / `random_unimodular_matrix` `upper_bound` | Size control by rejecting row operations past the bound — **only over ZZ and QQ** | `NotImplementedError`: the port's generic constructors work over any ring with `random_element()` and have no notion of absolute value, so the option is refused rather than ignored |
| `random_unitary_matrix`, `vector_on_axis_rotation_matrix`, `ith_to_zero_rotation_matrix`, `hadamard_bound` off RDF | Implemented over RDF/CDF | `NotImplementedError` naming the requirement (`sqrt` / trigonometric functions over an inexact ring) |
| `random_diagonalizable_matrix` | — | `NotImplementedError('unexpected eigenvector layout')` on an internal invariant violation, rather than emitting a matrix that is not diagonalizable |
| `PARI Z_factor` above MPQS's ceiling | Would continue | `NotImplementedError` naming `mpqs.c` and the 107-digit decline (`mpqs.h:400`), rather than returning a composite as prime |
| Every bounded search | — | See [Bounded Search Budgets](#bounded-search-budgets-and-measured-thresholds) |

### Rationale

1. **A wrong answer is worse than no answer** in a library whose stated goal is exact behavioural
   equivalence and whose consumers are cryptographic.
2. **Grep-ability** — CLAUDE.md rule 7 requires unimplemented paths to carry
   `SAGE_NOT_IMPLEMENTED`, so the gaps are discoverable.
3. **Bounded search where the algorithm is exponential** — an explicit cap documents the real
   reachable range instead of hanging.

### Trade-offs

- Calls that "worked" in an earlier version now raise. In every case the previous answer was wrong,
  ignored an argument, or was a placeholder — but this is a **source-breaking change** for callers
  that were not checking.
- Some inputs SageMath handles still fail here; each such case is registered in
  [Part II](#part-ii--open-fidelity-gaps).

### Behavioral Impact

Errors where there used to be plausible garbage. Every message names the missing dependency or the
algorithm that would be needed.

> **One place still violates this policy** and is registered as an open gap rather than defended
> here: `gauss_sum` silently drops terms when a ring lacks an optional method
> ([Arithmetic Functions](#arithmetic-functions-not-delegated-to-pariflint)).

---

## Bounded Search Budgets and Measured Thresholds

Several upstream algorithms are unbounded loops backed by a fallback this port does not have, or
are asymptotically better than what is reachable here. Rather than hang or silently degrade, the
port imposes an explicit budget and fails with a message naming what is missing. Collected here
because they share one rationale and are easy to mistake for arbitrary magic numbers.

| Site | Upstream | Budget here | What happens at the limit |
|------|----------|-------------|---------------------------|
| PARI insisting-ECM | Loops **forever**, because MPQS backs it up | `FactorOptions.ecmRounds`, default 4 (MPQS backs it up here too, as in PARI) | Bounded work per insisting round; without the bound `Z_factor` could hang indefinitely |
| MPQS polynomial budget | Stops only on "ran out of primes for A" or Gauss failure | `MpqsOptions.maxPolys` / `FactorOptions.mpqsMaxPolys`, default 0 = **unbounded, as PARI** | Not set in any production path; it exists so the "every stage failed" branch of `Z_factor` can be tested in seconds instead of an hour |
| Integer polynomial factorization prime search | `for ( ; ; p = n_nextprime(p, 0))`, unbounded | Unbounded, plus a residual cap of `1000 + 4·len(f)·(maxbits+10)` rejected primes | `ValueError`. Can only fire on an input that violates the precondition (non-squarefree `f`, or `f(0) = 0`), i.e. it converts an infinite loop into an error |
| van Hoeij precision doubling | `while (!check_if_solved(...))`, unbounded | 32 Hensel doublings | `ArithmeticError`. Never reached: the hardest case measured doubles twice |
| Gröbner S-pair queue | Delegates to Singular; always terminates | `maxIterations`, default 10 000 | `ArithmeticError`. A truncated set is **not** a Gröbner basis: `contains()`, `reduce()` and `dimension()` would return silently wrong answers (with a cap of 3, an ideal member demonstrably failed to reduce to 0) |
| `nfgaloisconj` precision escalation | PARI escalates inside `galoisgen` | 12 attempts, squaring the p-adic precision each time | `NotImplementedError`. The Gram-Schmidt certificate means an inconclusive result is *known* to be inconclusive, never a wrong answer. `galconj.ts`'s `galoisconj4` is the unbounded route |
| `_nf_monic_cubic_has_root` (the 2-division-polynomial test in `Frobenius_filter`) | PARI `nffactor` | 200 rational primes for the irreducibility certificate; modulus `2^2048` for the root reconstruction | `NotImplementedError` naming `nffactor`. Not reached by any of the 189 cross-checked curves or the three doctests |
| `is_similar` intertwining search | Sage raises `RuntimeError` instead | 200 pseudo-random kernel combinations | `ArithmeticError` with Sage's message text. Unreachable in testing: the proportion of units in a centralizer algebra over `F_q` is at least `prod(1 - q^-i) >= 0.288` |
| `rook_vector` naive algorithm | `ButeraPernici` / `Ryser` / `Godsil` | 50 positions **and** `k > 5` | `NotImplementedError` naming the two faster algorithms |
| Class group of a degree > 2 field | `bnfinit` (subexponential) | Minkowski bound `<= 10^6`, plus the two rigorous certificate cases | `NotImplementedError`. See [Number Field Class Groups](#number-field-class-groups-units-and-galois-closure) |
| `voronoiCell` / exact SVP / BCH minimum distance / BCH field embedding | Backend-accelerated | rank 24 / rank 30 / `q^k > 2^17` / `\|E\| > 2^22` | `NotImplementedError`, except exact SVP, which **silently approximates** — registered as an open gap under [Lattices](#lattices--exact-svp-rank-cap) |
| `number_of_partitions` / `prime_pi` | FLINT / primecount | `n <= 10 000` / `n <= 10^7` | `NotImplementedError`. Registered as an open gap under [Arithmetic Functions](#arithmetic-functions-not-delegated-to-pariflint) — the upstream algorithms are vendored |
| Quadratic class group | `quadclassunit` (subexponential) | `CLASS_GROUP_DISC_BOUND = 2 000 000` | `NotImplementedError`. **Removable today** — see [Quadratic Class Numbers](#quadratic-class-numbers-not-delegated-to-buchquad) |
| `ellcard` Schoof/Shanks crossover | PARI switches to SEA at `expi(p) >= 56` | Base Schoof from `expi(p) >= 96` | Not a failure — a *measured* threshold. But the dispatch target is now the wrong one; see [parigp-ts Elliptic Curves](#parigp-ts-elliptic-curves--sea-dispatch-and-isogeny-stubs) |

### Rationale

1. **A hang is the worst failure mode** in a library, worse than an exception: it gives the caller
   nothing to act on and no signal that a dependency is missing.
2. **The budget documents the reachable range.** Where the upstream loop is unbounded only because a
   *later* stage catches the hard cases (MPQS behind ECM; van Hoeij behind Zassenhaus), removing that
   stage without adding a bound converts "slow" into "never returns".
3. **Thresholds must be measured, not copied.** PARI's 56-bit SEA crossover is correct *for PARI*;
   transplanting it into a port that only had base Schoof would be a fidelity gesture that makes the
   function unusable.

### Trade-offs

- Inputs SageMath handles can fail here, and the failure is a hard error rather than a long wait.
- The constants are tuning choices, not upstream values, and would need re-measuring on different
  hardware or after asymptotic improvements to the underlying arithmetic.
- A caller who *wants* to spend more time has an escape hatch only where one was added
  (`FactorOptions.ecmRounds`, `solve_integer`'s `factorization`, `qfbsolve`'s `fa`).

### Behavioral Impact

Every budget exhaustion raises with a message naming the missing upstream routine. No budget
silently truncates a result — with the single exception of exact SVP above rank 30, which is
registered as an open gap.

---

## Random State and Seeding

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| RNG core | `sage.misc.randstate` wraps `gmp_randstate_t` from `gmp_randinit_default` (GMP's MT19937) | `RandState`, a verbatim port of GMP 6.3.0 `rand/randmt.c` (`__gmp_mt_recalc_buffer`, `__gmp_randget_mt`, the 624-word `default_state` table) |
| Seeding | GMP `randseed_mt` (`rand/randmts.c`): seed mod 2^19937−20027, `+2`, `mangle_seed` (`r^1074888996 mod 2^19937−20023`), `mpz_export` into `mt[1..623]`, bit 19936 into bit 31 of `mt[0]`, 3 warm-up recalcs, `mti = 2000 % 624` | **Identical** — ported line for line, including GMP's non-canonical `reduce:` loop |
| Seed `0` | `randstate.pyx:556` skips `gmp_randseed` entirely when `seed` is falsy, so GMP's *unseeded* `default_state` buffer is used | Identical |
| `mpz_urandomm` | `mpz/urandomm.c`: bit length minus the power-of-two adjustment, 80-iteration rejection, then subtract `n`; `n = 1` returns 0 without drawing | Identical port |
| Second generator | `randstate.python_random()` returns a CPython `random.Random` seeded from `ZZ.random_element(1<<128)` | `PythonRandom`, a port of CPython `Modules/_randommodule.c` + `Lib/random.py` |
| Affected modules | Any `.random_element()` / random sampling | Same |

### Rationale

1. **Bit-exact parity is achievable.** GMP is not vendored under `reference/`, so the GMP 6.3.0
   sources were obtained and `randmts.c` / `randmt.c` / `urandomm.c` ported directly. Seeded streams
   are **identical to SageMath's**, not merely same-distribution.
2. **Centralization** — matches Sage's single-randstate model.

### Trade-offs

- MT19937 is not cryptographically secure (neither is Sage's).
- `python_random()` reproduces an upstream quirk that is easy to mistake for a bug: it derives its
  seed from `current_randstate()`, **not** from the receiver (`randstate.pyx:623` calls
  `ZZ.random_element`, which reads the global state). So
  `set_random_seed(0); randstate(314159).python_random().random()` is a *seed-0* value. A port-only
  `seed` parameter is offered for callers who want a self-contained stream.
- The GMP sources are not vendored under `reference/`, so the port cannot be re-diffed against them
  in-tree; the ported functions carry `gmp-6.3.0/<file>:<line>` citations instead.
- Random **matrix** constructors and `crypto/lattice.ts`'s `IntegerLattice.gen_lattice` draw in a
  different *order* from Sage's, so the values still differ for the same seed even though the
  underlying stream matches.

### Behavioral Impact

Seeded streams match SageMath exactly, verified against a C oracle linked to the installed
libgmp 6.3.0 (7 seeds × 8 measurement families) and against SageMath 10.3 itself: `c_random()` after
`set_random_seed(1207)` = `2008037228`; `c_rand_double()` after `set_random_seed(2718281828)` =
`0.22437207488974298`; seed 0's first `c_random()` = `968665204`; `python_random().random()` after
seed 314159 = `0.29929142114291285`; `random_below(1n)` consumes nothing; a 24-value
`IntegerModRing(11)` stream.

The former unconditional-draw defect in `ZZ.random_element` is resolved: it now consumes
`c_random()` before every distribution branch, as Sage does. The `rand_stats` oracle covers
the integer streams, and the September RationalField audit adds 24 shared seeded cases for
`QQ.random_element`, which now delegates numerator and denominator draws through ZZ.

---

## Vendored SageMath 10.9.beta4 vs Installed 10.3

CLAUDE.md directs the port at the vendored tree under `reference/sage`, which is 10.9.beta4. The
locally installed oracle is SageMath 10.3. Where the two disagree, the port follows the **vendored**
source, so a reader comparing against a stock 10.3 will see differences that are version drift, not
port defects.

| Site | Installed 10.3 | Vendored 10.9.beta4 / port |
|------|----------------|----------------------------|
| `Integer.radical(0)` | Raises `ArithmeticError` | The zero guard in `categories/unique_factorization_domains.py:283-284` returns 0. The integer and free-arithmetic oracles explicitly apply that vendored guard before delegating nonzero cases to installed Sage |
| `Rational.round()` default mode | `'away'`, with a `DeprecationWarning` | `rational.pyx:3403` is `def round(Rational self, mode="even")`. The port defaults to `'even'` (`rational.ts:587`). All six modes (`toward`, `away`, `up`, `down`, `even`, `odd`) exist in both and agree value for value |
| `BinaryQF.is_reduced(D=0)` | Returns `False` | Raises `ValueError('the quadratic form must be non-singular')` (Sage issue #37635 rewrote `is_reduced` after 10.3). The port implements the 10.9 behaviour; the property-test cases route `D = 0` through other methods and short-circuit `canonical` to the literal string `'singular'` |
| `discrete_gaussian_lattice.py` | Brute-force `_normalisation_factor_zz(tau=3)`; `c` a property; repr `Discrete Gaussian sampler with σ = %f`; `_c_in_lattice` ignores whether `_G == 1` | The vendored module, loaded by path in the property-test oracle. Measured differences: `nf(ZZ^8, 0.5)` is `6.81052960784091` in 10.3 vs `6.82492448921763` in 10.9; `nf(ZZ^3, 1.0)` `15.5284660320764` vs `15.7496101985309`; sampling on basis `[[1,3,0],[-2,5,1],[3,-4,2]]` takes `_call_in_lattice` in 10.3 and `_call` in 10.9. Comparing against 10.3 would have reported four false port defects |
| `groups.generic.order_from_multiple` | A bare `assert` (message-less `AssertionError`) | `ValueError(f"The order of P(={P}) does not divide {M}")` (`generic.py:1361`). Also `order_from_bounds(P, None, …)` and `IntegerMod.log(b, order, check=)` (`integer_mod.pyx:795-798`) exist only in 10.9 |
| `basis_for_quaternion_lattice` | `reverse=False` (deprecated); `maximal_order` fails with `ValueError('basis must have rank 4')` for invariants such as `(-4,-28)`, `(-292,-732)`, `(-48,-564)`, `(-436,-768)` (Sage issue 37417); no `is_definite`/`ramified_places`/`order_with_level`/`pushforward`/`pullback`/`reduced_basis`/`is_principal(certificate)` | `reverse=True` default, issue 37417 fixed, all newer methods present. Printed ideal/order bases can differ from 10.3's; the **lattices** are identical, verified by comparing echelon basis matrices |
| `discrete_gaussian_integer` `precision='qq'` message | Adds a trailing period | `discrete_gaussian_integer.pyx:400` has none; the port matches the vendored source |
| `Frac(GF(p)[x])` element comparison for small odd `p` | `FpTElement._richcmp_` compares `(numerator, denominator)` lexicographically; its own docstring says "the ordering is arbitrary" (`fraction_field_FpT.pyx:376`) | Native FpT ordering is now selected for the same small-prime domain. See [Function Fields](#function-fields) |
| Free modules over `QQ[x]` | `_echelonized_basis` lacks the `if basis.universe().coordinate_ring() == ambient.base_ring(): d = 1` guard | The guard is present. A verifier re-running the goldens against a stock 10.3 will see 37 of 250 sweep cases differ, all over `QQ[x]`, all by a rational unit |

### Rationale

1. **CLAUDE.md is explicit**: port the vendored upstream. Pinning 10.3 where 10.9 fixed a bug would
   be pinning a known upstream defect.
2. **The divergences are individually attributable** — each is traced to a specific upstream change
   or issue number, not to an unexplained mismatch.

### Trade-offs

- A user validating against a stock SageMath 10.3 will see differences that are not defects. Each row
  above names the version in which upstream changed.
- Property-test areas that need 10.9 behaviour transcribe the relevant snippet verbatim at the call
  site, marked `VENDORED` with its reference line number, rather than skipping the case.

### Behavioral Impact

Values match the vendored source. Where 10.3 and 10.9 agree, the port agrees with both.

---

## Number Fields — Exactness-Driven Divergences

`rings/number_field/` carries a **number-field kernel**, `pari_nf.ts`, which ports the PARI routines
SageMath delegates to: `nfbasis`/`nfdisc` (Pohst-Zassenhaus round 2, Cohen 6.1.8, i.e. PARI's
`maxord`), `idealprimedec` (Dedekind-Kummer **and** Buchmann-Lenstra round 4, `base2.c:2248`
`primedec_aux`, `:2150` `pradical`, `:2185` `pol_min`), `nfgaloisconj` (LLL-based, no degree cap),
`quadunit`/`quadunitnorm` (`quad.c:281` `quadunit_uv_basecase`), and `polisirreducible`. Archimedean
embeddings live in `number_field_embeddings.ts`. Where the kernel deviates from PARI it is in the
direction of **exactness**; the *architectural* problem that it lives in `sagemath-ts` at all is a
separate open gap.

| Aspect | SageMath / PARI | sagemath-ts |
|--------|-----------------|-------------|
| `fujiwara_bound`'s `log2\|c_i\|` (`rootpol.c:1628`) | Floating point | `bitLength(c_i)`, so the returned bound is *proved* rather than rounded (at most one extra bisection level) |
| `polsolve`'s Newton refinement (`rootpol.c:2139`) | Floating-point Newton | Exact bisection over dyadic endpoints, using the exact integer sign of the `ZX` at each midpoint. `O(prec)` evaluations instead of `O(log prec)`, but every endpoint is a proved bound and no `t_REAL` kernel is needed |
| `polroots` (`cleanroots`/`all_roots`, the ~1600-line Schoenhage splitting-circle method) | Certified numerically | Replaced by a Durand-Kerner estimate in doubles — sound because SageMath's own `complex_roots` (`complex_roots.py:154`, `refine_root.pyx:27`) treats the estimator as an **untrusted** black box and certifies it with interval Newton. That certification is ported verbatim, so every returned box is *proved* to contain exactly one root |
| Ramification groups | `idealramgroupstame`/`idealramgroupswild` (`base1.c:931-1038`) use PARI's uniformiser/residue-generator shortcut | `G_v(P) = {s in D(P) : v_P(s(w) − w) >= v+1 for every w in a Z-basis}` — the definition. The shortcut needs `nf_get_diff`, `zk_to_Fq_init`/`modpr_genFq` and `ZC_galoisapply`, none ported. `n` valuation tests per group element instead of one, and no differente-based cap on the filtration length. Reproduces upstream on every doctest, including the wild `p = 2 \| e = 8` case with breaks `{1, 3, 5}` |
| `NumberFieldIdeal.valuation` | Divides by PARI's anti-uniformiser `pr_get_tau` (`base4.c:3007`) | Climbs `P^k` with the exact HNF membership test, bounded by `v_p(N(x))`. `pr_get_tau` is part of the `prid` structure the port does not build |
| `primedec` return shape | PARI's 5-component prime structure with a uniformizer and anti-uniformizer | `{gens, e, f}`; `number_field.ts` then searches for a two-element representation `(p, alpha)` and **certifies** it with `N((p,alpha)) == p^f`, falling back to the full generating set. No wrong ideal can be returned; the printed form can have more generators than Sage's |
| `quadunit` | Switches to the product-tree variant `quadunit_uv` (`quad.c:429`) at `D >= 2000000` | Basecase only. Identical `[u,v]`; the product tree is a big-integer-multiplication speed optimisation (0.5 ms for the ~800-digit unit of `D = 511681`) |
| Regulator evaluation | MPFR | Double-precision embeddings (the regulator is transcendental, so no exact representation exists with the primitives here), through an overflow- and cancellation-safe `quadraticLogAbs` (log-sum-exp over bigint bit lengths). The ~250-digit unit of `Q(sqrt(1000003))` gives a finite `R = 576.646` where a naive double evaluation returns `Infinity` |
| Proof flags | GRH-conditional results are flagged | No proof flags |
| Affected modules | `sage/rings/number_field/`, `pari/src/basemath/base1.c`, `base2.c`, `base4.c`, `rootpol.c`, `quad.c` | `packages/sagemath-ts/src/rings/number_field/` |

### Rationale

1. **Every replacement is provably an upper bound or an exact evaluation**, so it can only agree
   with upstream or be more correct. Where the upstream value is *observable* rather than an internal
   heuristic, upstream's arithmetic is reproduced instead.
2. **SageMath's own architecture licenses the root-finding split** — it certifies an untrusted
   estimator, so the estimator may be anything.

### Trade-offs

- Slower: `O(prec)` bisections rather than `O(log prec)` Newton steps; `n` valuation tests per
  ramification-group element.
- The regulator is a double, so very large regulators lose relative precision.
- The printed generating set of a prime ideal can be longer than Sage's.

### Behavioral Impact

Maximal orders, integral bases, field discriminants, prime decomposition, ideal arithmetic,
automorphisms and archimedean places agree with PARI on every value tested. Dedekind's classical
inessential-discriminant example `x^3 − x^2 − 2x − 8` at `p = 2` splits into three primes with
`e = f = 1`. `nfrootsof1` **proves** the number of roots of unity, so `zeta_order()`/
`torsion_order()` either return the proved value or throw — they never return an invented one.
`QuadraticField(D)` uses `x^2 − D` verbatim as Sage does. The zero ideal is not reported prime
(Sage's `idealismaximal` does not accept it at all). `NumberFieldElement.is_unit()` implements Sage's
*field* branch; the ring-of-integers test is `is_integral_unit()`.

---
## Polynomial Roots and Factorization

Covers the supported `roots()`/`factor()` interfaces and remaining backend differences.

| Aspect | SageMath | sagemath-ts |
|---|---|---|
| Ring support | Many exact and approximate rings | Base-ring finite fields, ZZ and QQ only |
| `roots()` options | `ring`, `multiplicities`, `algorithm` | Supports `{ multiplicities?: boolean }`; ring overrides and algorithm options are unavailable |
| ZZ roots | Dense factorization through degree 100; sparse exponent gaps and block GCD above 100; factor subproblems select NTL below degree 30 or above 300, PARI otherwise | Same dense/sparse dispatch, content handling, sparse derivatives, root order and NTL/PARI calls; FLINT handles integer block GCD |
| QQ roots | PARI factorization, monic normalization and Sage factor order | Delegates to PARI QX_factor, filters linear factors and sorts exact rational roots by multiplicity then descending value |
| Public ZZ/QQ `factor()` | Dense ZZ selects NTL/PARI by degree; QQ selects PARI | Same dispatch; integer content is factored before the polynomial, and rational factors are made monic |
| `factor()` result | `Factorization` stores its unit separately | Array of factor/multiplicity pairs, including a nontrivial unit as a degree-zero factor; the constant 1 has an empty array, and ZZ constants include their prime factors |
| ZZ/QQ `is_irreducible()` | Cached; QQ factors its positive-leading primitive numerator over ZZ, including ZZ's degree dispatch | Same primitive normalization, degree dispatch and cached true/false results; unit entries are excluded from the ZZ factor count |
| Random contexts | Native process-level NTL stream/cache and PARI state | Roots, factorization and irreducibility use the NTL driver's private default context and shared PARI state; the public API does not expose NTL continuation |
| Finite-field roots | The default multiplicity-preserving method factors and extracts linear roots in factor order | Same extraction for prime and extension fields; distinct roots first use gcd(f, x^q-x). Extension factor order delegates to PARI universal comparison; zero preserves the selected factor/modular-power backend error |

### Rationale

Exact ZZ/QQ roots, factors and irreducibility use the bundled Sage algorithms and existing
dependency ports. The array factorization interface predates a separate Factorization type,
so retaining a nontrivial unit in that array preserves product reconstruction. NTL's private
context follows the constructor adapter and isolates mutable library state. The former local
FLINT-style Zassenhaus/van Hoeij driver remains available only through its test hooks; public
ZZ/QQ factor and irreducibility calls no longer use it.

### Trade-offs

Unsupported rings, ring overrides and root algorithm options remain unavailable. Callers must distinguish the extra
unit entry from genuine constant prime factors over ZZ. NTL calls do not reproduce
process-global random/cache continuation; the low-level NTL driver exposes explicit state
for testing that behavior separately. Finite-field irreducibility backends, caching and
algorithm options have separate contracts and remain under audit. Integer primality proof
limits still apply when testing constant integer polynomials.

### Behavioral Impact

ZZ/QQ root comparisons cover ordered values, multiplicities, result classes and exact zero
errors against executable bundled Sage methods. In particular, sparse ZZ roots preserve the
initial zero and then +1/-1 insertion order. Bundled QQ zero roots raise NotImplementedError;
installed Sage 10.3 instead propagates ArithmeticError, so the oracle executes the bundled
method. Factor comparisons preserve coefficients, multiplicities, unit adaptation and order.
The native PARI comparisons check content-before-polynomial effects and cached irreducibility
calls versus fresh polynomial objects. These checks do not establish whole-domain coverage
or process-global NTL-state parity. Distinct finite roots may have a different order because
their factor multiplicities are all one. Binary/extension zero distinct roots raise
`ZeroDivisionError("modulus must be nonzero")`; word-prime roots raise `ZeroDivisionError`
with an empty message; large-prime roots propagate `NTLError("ZZ_pX: division by zero")`.
The distinct-root oracle executes the bundled field method, using Sage's independent
modular-power routine where installed 10.3 predates the bundled backend wrappers. It
executes the bundled binary/extension zero guards and native word/large first reduction;
this adapter does not establish NTL state equivalence. Finite-field irreducibility and
cross-ring root overrides remain open. IntegerModRing base-ring roots follow the
[modular root contract](#modular-polynomial-roots-and-hensel-lifting), including its
documented upstream linear-case bug.

For polynomial printing and integer factor units, see
[Polynomials](#polynomials--printing-factor-shape-term-orders-and-base-rings).

---

## Finite Field Coercion and Backend Boundaries

The September audit compares 99 public members across `PrimeField`, `PrimeFieldElement`,
`FiniteFieldPrime`, its `FiniteFieldElement`, `FiniteFieldExtension`, and extension elements.
The following adapters and remaining boundaries are explicit; member dispatch coverage does
not establish full input-domain coverage.

| Aspect | SageMath | sagemath-ts |
|---|---|---|
| Scalars | Integer, bool, None, rational and string conversion through IntegerMod or PARI constructors | Integral JavaScript numbers denote integers; Integer wrappers, booleans, null/undefined, Rational values and integer strings now follow those paths. Explicit Rational conversion works when its denominator is invertible; implicit arithmetic with a Rational is rejected because QQ has no canonical map into a finite field |
| Mixed arithmetic | Canonical quotient maps, prime-subfield embeddings and named-field compatibility determine the result parent | All four element classes now share these rules. Compatible prime/modular operands promote into extension fields; both prime classes interoperate. Nonintegral numbers, null, strings and lists are rejected by ordinary arithmetic and compare unequal |
| Sequence multiplication | IntegerMod and prime elements implement __index__, so Python reflected multiplication repeats strings/lists; PARI extension elements reject it | Native modular/prime mul exposes string/list return overloads, including shallow list repetition and signed 64-bit index errors. Extension sequence errors match PARI. Host allocation limits and out-of-memory exceptions are not normalized between JavaScript and Python |
| Parent identity | Cached parents; PARI element conversion rejects a different parent | Repeated TypeScript constructors are separate objects. Compatibility uses characteristic, degree, generator name and defining modulus. Incompatible extension elements now raise TypeError, including constant elements; polynomial inputs first convert their coefficient ring |
| Extension representation | Givaro is the default for small fields; PARI is an available backend | Polynomial-basis iteration and coordinate-list construction follow PARI. Tests select `impl='pari_ffelt'` on installed Sage 10.3 and provide the matching Conway modulus; they do not assert Givaro's different iteration order |
| Integer conversion hooks | `_integer_(ZZ=None)` lifts native residues and prime-subfield extension elements; nonconstant extension elements raise ValueError | Modern prime and extension elements now expose the same hook; the optional parent is ignored. Integer constructors call it without interpreting extension coordinates as a base-p integer |
| Generator indices | Fields test Python truthiness; quotient rings and QQ test equality with zero | `gen(n: unknown = 0)` validates the port's scalar/list/finite-element representations and preserves the distinct IndexError messages. General user-defined Python truthiness and equality protocols remain outside this bounded adapter |
| Integer representation | `from_integer(n)` with `0 <= n < q`; `fetch_int` is an older spelling | Existing `fromInteger(n: bigint)` spelling is preserved. Invalid inputs now raise the same ValueError instead of wrapping |
| Seeded randomness | Prime fields use Python randrange. The bundled newer reference also uses one randrange(q) draw for an extension field with no arguments; installed 10.3 uses a vector randomizer | Both prime implementations and extensions now use the ported Python random stream. The extension oracle executes bundled finite_field_base.pyx:1058-1059 using Sage's own randrange and from_integer. Random-element keyword options remain unimplemented |
| Powers | Prime elements inherit IntegerMod's native/int64/GMP implementation; the original exponent type and magnitude select the error path | Both prime-element classes reuse IntegerMod.pow. IntegerLike, Rational, number, boolean, string and null exponents follow Integer coercion; JavaScript integral numbers denote Sage integers. Non-native exponent objects use the GMP error path even when they convert to a small integer |
| Square roots | Prime-field `sqrt(extend=True, all=False)` returns a root in the base field or GF(p^2); `all` returns base-field roots or follows the backend's exception path | Both prime-element implementations now accept `{extend, all}` with Sage defaults. Explicit `extend:false` narrows the return type and rejects nonsquares. Root arrays are sorted; the small-modulus search returns an empty array before extension while the general extension/all branch raises Sage's NotImplementedError |
| String expressions | PARI field construction parses through its polynomial ring | Arithmetic strings now use the same tokenization, grammar and coercion order. Fractions must cancel to a polynomial before extension reduction. General parser function calls, matrices and sequences remain unimplemented; see Polynomial String Parsing |
| Port conveniences | Native arithmetic, characteristic(), vector(), field iteration | Arithmetic method names, `.p`, `.coefficients()`, and `.quadratic_non_residue()` are explicit adapters. The last returns the first nonsquare in ordered prime-field iteration; GF(2) raises a port-specific ValueError |

Prime-field square-root defaults were repaired after the factory audit. Existing TypeScript
callers that require the original field should use `sqrt({extend:false})`; the default return
type now includes an extension element. Legacy elements delegate to the modern implementation
and convert base-field results back to their original parent. This shares the PARI-backed
square-root algorithm. Extension roots use the exact modulus x^2-a and generator name sqrt<a>.
The integer_mod.pyx small-search branch is reproduced without enumerating all residues:
its empty-result condition is computed directly and all square roots come from the existing
PARI kernel. The bundled source adds an explanatory NotImplementedError message to the
extension/all branch; the Sage 10.3 oracle adapter replaces only that older empty message.
The extension-element class itself still has no general sqrt method.

Direct prime parent constructors now use Integer coercion before checking primality. The
legacy check:false option retains the generic-ring positive-order check and admits order
one, as Sage does. Nonprime unchecked-parent arithmetic remains under active audit and is
not covered by the prime-field guarantees above.

The factory oracle has two explicit version adapters. Bundled Sage 10.9 validates prime-field
names, which installed 10.3 ignored; tests execute Sage's `certify_names` for that branch.
Installed 10.3 also misses small even perfect powers such as 676 (upstream issue 40846).
Tests use Sage's PARI `ispower` result, matching the corrected bundled `Integer.perfect_power`
path. The TypeScript factory now delegates to `parigp-ts/Z_isanypower`.

For a zero defining polynomial, the port preserves the executed 10.3 errors: inverse-of-zero
in characteristic two and `ValueError('leading coefficient must be invertible')` in odd
characteristic. The newer bundled FLINT `nmod_poly_make_monic` throws `FLINT_DIVZERO` for zero;
the precise Python exception mapping has not been executed on 10.9. This error-message/version
boundary remains unverified against that runtime and is not counted as proven 10.9 parity.

**Rationale:** Preserve the existing TypeScript representations while matching the selected
Sage backend and its actual coercion rules. A Rational's explicit conversion and implicit
arithmetic are different operations in Sage.

**Trade-offs:** Host memory exhaustion and allocation limits remain runtime-specific;
sequence outputs and representable index boundaries are compared exactly. Givaro iteration,
general parser function/sequence/matrix input,
random keyword arguments, and Python parent object identity are not yet reproduced. The
small Conway database and broader factory options remain tracked below. Degree-one Conway
entries for every existing characteristic were added and compared with Sage in 0.0.26.

**Behavioral impact:** Repaired scalar conversions, bounds, parent checks and seeded random
streams now match the selected Sage paths. The listed unimplemented input branches still
differ and remain part of the ongoing audit. `Fp_order` follows arith1.c:2603-2635 and
bb_group.c:670-713, including the native-word and recursive factor-splitting paths.

---

## Finite Field Constructors and Display

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `GF(q)` for a prime power | `GF` and `FiniteField` are the same factory object | The **exported** `GF` (`rings/finite_rings/index.ts:21`, re-exported from `finite_field_extension.ts:1254` as `GFExtended`) and `FiniteField` both build prime powers: `GF(9)` is `Finite Field in a of size 3^2`. Only the *module-local* `finite_field_constructor.ts:58` `GF` is narrowed to prime fields, returning the concrete `FiniteFieldPrime` |
| Extension field constructor spellings | `GF(p^n, 'a')` | `GF(p^n)` / `FiniteField(p^n)` / `GFExtended()` / `GFpn()` |
| Element display | `repr(a)`, `a.lift()` | Prime fields: `a.toString()`, `a.repr()`, `a.value`, `a.lift()`. Extension elements have `.repr()` but no `.value` |
| Affected modules | `sage/rings/finite_rings/` | `packages/sagemath-ts/src/rings/finite_rings/` |

### Rationale

1. **Type inference** — retyping the module-local `GF` to a union would break the many
   elliptic-curve and matrix call sites that depend on the concrete `FiniteFieldPrime`.
2. **Interop** — `.value` provides direct `bigint` access for prime-field elements, which is what
   downstream cryptographic code consumes.

### Trade-offs

- The two names are not interchangeable *inside* the module: the local `GF` narrows to prime
  fields, `FiniteField` does not. Consumers of the package are unaffected.
- Extension-field elements do not expose `.value`, so generic code must branch.

### Mitigation

Retype the module-local `GF` to the union once the concrete-type call sites are migrated, and export
a single factory under both names.

### Behavioral Impact

`GF(p^n)` and `FiniteField(p^n)` both construct extension fields at the package level. The
`name` and `modulus` options are now accepted by both factories — see
[Finite Fields](#finite-fields--conway-table-and-minimal-polynomials).

---

## Generic Group API and DLP

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `discrete_log` options | Supports `bounds`, `algorithm` (`bsgs`, `rho`, `lambda`), `verify`, and `ord=oo` | Implements exactly Sage's `bounds=None, algorithm='bsgs', verify=True` path — including the repair of an `ord` that is a proper multiple of the base's order, and the `<30` linear branch of `bsgs`. `ord` is optional (`generic.ts:703`) and falls back to the base element's order |
| `discrete_log_rho` options | `ord` optional; configurable `hash_function` (`generic.py:659`) | `ord` is required; both implementations require it to be prime. No custom hash function in the TypeScript API |
| Hashing in DLP algorithms | Element hashing/equality | BSGS and rho memory bucket by the existing element string key, then resolve collisions using object identity and element equality. Pollard walk partition hashes retain their documented adapters |
| `order_from_multiple` | `(P, m, plist=None, factorization=None, check=True, operation='+', …)` | `(P, m, factorization?, operation='+', identity?, inverse?, op?, options?: {plist?, check?})` — `check` defaults to **true** and `plist` is honoured, but they are passed in a trailing options object |
| Generic utilities | `linear_relation`, `merge_points`, `structure_description` | Not implemented |
| Affected modules | `sage/groups/generic.py` | `packages/sagemath-ts/src/groups/generic.ts` |

### Rationale

1. **Incremental porting** — the API surface is narrowed to core cryptographic use cases; GAP-backed
   utilities (`structure_description`) have no backend here.
2. **JS runtime limits** — there is no standard hash for custom objects, so string keys are used
   internally.
3. **Positional-argument compatibility** — `plist`/`check` are passed in a trailing options object
   rather than in Sage's positions, because `rings/finite_rings/integer_mod.ts` and other modules
   call `order_from_multiple` positionally as `(a, m, factorization, operation)`.

### Trade-offs

- Fewer algorithm choices and no bounds handling for discrete logs.
- `discrete_log_rho` requires an explicit order; both implementations reject composite orders.
- BSGS and rho memory handle non-unique string keys correctly; BSGS also distinguishes separate NaN elements. Equal elements must still produce consistent string keys; a collision-heavy key can degrade bucket lookup performance.
- Argument *positions* for `order_from_multiple` differ from Sage's even though the semantics match.

### Behavioral Impact

Calls using `bounds`, `algorithm` or `verify` on `discrete_log()` are unsupported.
`order_from_multiple()` honours Sage's `check=True` default, so
`order_from_multiple(Mod(2,7), 5, '*')` raises as Sage does rather than returning 5 for an element
of order 3. `has_order` **does** accept a `Factorization` (`generic.ts:1220` branches on
`Array.isArray(n)`). `discrete_log` in a 2^30-order subgroup runs in milliseconds because the
Pohlig-Hellman loop follows Sage's verbatim.

---

## Matrix Module Algorithm Substitutions

Where a SageMath matrix routine delegates to a backend this port lacks, an equivalent exact
algorithm is used instead. The *results* are identical unless stated.

| Function | SageMath | sagemath-ts |
|----------|----------|-------------|
| LLL | fpLLL (default) or NTL, floating-point Gram-Schmidt with precision escalation | Cohen 2.6.7 **exact integral LLL** in bigint (Gram determinants `d_i`, `lambda_ij = d_j·mu_ij`; every division exact). Default `delta` is Sage's 0.99; `eta` is accepted but advisory, since the exact algorithm always achieves `\|mu\| <= 1/2` |
| LLL on dependent rows | fpLLL/NTL run MLLL and return the zero rows first | Dependence is detected exactly, the generating set is replaced by the nonzero rows of its Hermite normal form (same lattice), and `nrows − rank` zero rows are prepended |
| BKZ | Full enumeration/pruning | Repeated LLL passes (approximate) |
| Determinant over ZZ | Multimodular / LinBox | Closed forms for `n <= 3`, fraction-free Bareiss above (`matrix_integer.ts:249-250`). Same value; different asymptotics |
| Frobenius form | PARI `matfrobenius` | **All three flags** are a verbatim port of PARI's `RgM_Frobenius` (`alglin2.c:428-720`: Storjohann's Lemmas 9.14/9.18, Ozello's theorem 4, `_frobTransL/D/S`, `_minpoly_polslice/listpolslice/dvdslice`) over exact rational arithmetic, 1-indexed to match `gcoeff(M,i,j)` line for line. `flag=2` returns `[F, B]` as `Rational[][]`; an out-of-range flag raises `ValueError('incorrect flag in matfrobenius')` — PARI's own text with the port's error class |
| `right_kernel_matrix` over Z/nZ | PARI `matkermod` for composite `n` | Prime modulus is a faithful port of `matrix_modn_dense_template.pxi:2072` (all three `basis` formats); composite `n` delegates to `parigp-ts` `matkermod`, mirroring `:2136`'s fallback |
| `matrix_modn` determinant | LinBox for prime `p > 2` | `n <= 3` uses Sage's naive formulas; `n >= 4` uses centered-lift-to-ZZ + fraction-free Bareiss |
| HNF transformation matrix | `fmpz_mat_hnf_transform` | Classical row operations; `U` differs from FLINT's for rank-deficient input (both satisfy `U·A == H`, and `H` matches Sage exactly) |
| `is_positive_definite` / `is_positive_semidefinite` | Eigenvalue signs off the 1×1 diagonal blocks of the Bunch-Kaufman `block_ldlt` factorization | After the same ring check and `is_hermitian` test, off the characteristic polynomial: with `charpoly = sum c_i x^i` the elementary symmetric functions are `e_k = (−1)^k c_{n−k}`, and a Hermitian matrix is positive (semi)definite exactly when every `e_k` is `> 0` (`>= 0`). Provably equivalent, exact, one division-free charpoly |
| `is_similar` | `A.rational_form() == B.rational_form()` | Compares, for every monic irreducible factor `h` of the charpoly, the multiset of elementary-divisor exponents recovered from the growth of `dim ker(h(A)^k)` — precisely the data determining the rational canonical form. Also detects a 6×6 counterexample that charpoly+minpoly misses |
| `is_similar(transformation=true)` | `matrix2.pyx:13052-13070`: Jordan forms over the fraction field, then over the algebraic closure, else `RuntimeError` — its own doctest at `:12918` shows two *provably similar* matrices over `GF(7^2)` for which Sage raises | Tries Sage's Jordan-form formula first (so it matches Sage exactly), then falls back to solving the intertwining equation `B X = X A` as an `n²×n²` homogeneous system. **This succeeds in cases where Sage raises.** Every candidate is verified (`B·T == T·A`, `rank(T) == n`) before being returned |
| `norm(A, 2)` | `matrix2.pyx:16466-16471`: `change_ring(CDF)`, `A^H·A`, numerical SVD, `max(S).real().sqrt()`, returning `RDF` | The same route: `_entryToCDF` mirrors `change_ring(CDF)`, then `A^H·A` and an SVD in double precision. RR and CC entries are accepted; number-field entries raise `NotImplementedError` (no distinguished complex embedding is wired in); positive characteristic raises `TypeError`, as upstream |
| `jordan_form(transformation=true)` | `matrix2.pyx:12259-12312` + `_jordan_form_vector_in_difference` (`:20895`) | Ported line for line, reproducing Sage's **exact** `P` rather than merely a valid one (Sage's `right_kernel().basis()` is echelonized, and the chains depend on which kernel vector is picked first). Eigenvalues come from `A.charpoly().roots()` (`:12228`) |
| `jordan_decomposition` | `matrix2.pyx:12383-12400`: a Newton iteration on the minimal polynomial, which succeeds **even when the eigenvalues are not in the base field** | Reads `D` and `N` off the Jordan form. Correct whenever the eigenvalues lie in the base field; otherwise propagates `jordan_form`'s `ArithmeticError` |
| `krylov_kernel_basis` | `matrix2.pyx:20343-20478` — builds the kernel directly from the Krylov basis as `relation = D·C^-1` | Ported |
| `change_ring(matrix, ring)` | `matrix0.pyx:1666-1715`, relying on the coercion framework | No coercion framework, so `_coerce_entry` asks the target ring to convert and, on failure, builds the canonical morphism only in the two cases where one **provably** exists: `QQ -> R` (`n/d ↦ R(n)·R(d)^-1`, raising if `R(d)` is not a unit) and `Z/mZ -> R` when `char(R) \| m` (or `m = 0`). Fallbacks still raise when unsupported. The previous claim that `Z/8 -> GF(7)` must be refused was wrong: explicit entry conversion succeeds in Sage and now succeeds through the corrected prime-field constructor |
| `QR` | Scales each column by `1/sqrt(<v,v>)`, so `Q` is unitary and `R` has non-negative diagonal; raises `TypeError` when the fraction field has no square roots | Delegates to `gram_schmidt_noscale`: `Q`'s columns are **orthogonal but unnormalized** and `R` differs by the diagonal factor. On `[[1,2],[3,4]]` over QQ, `Q = [[1,3/5],[3,-1/5]]`, `R = [[1,7/5],[0,1]]`, `Q·R = A` exactly, where Sage raises. `full` defaults to `true` as in Sage |
| `block_ldlt` pivot selection | Bunch-Kaufman compares `\|A_kk\|`, `omega_1`, `omega_r` against `alpha = (1+sqrt 17)/8` in C doubles | Exactly that rule when the base ring's elements expose `abs()` (QQ, RR, ZZ); over rings with no absolute value an exact rule. Pivot choice affects only numerical stability — `P^T A P == L D L^T` holds exactly |
| `principal_square_root` | Returns `False` when `check_positivity` and not positive definite; works over the algebraic closure | Skips the positivity check and diagonalizes over the base ring; raises `ArithmeticError` when not diagonalizable there. Over a finite field "the" square root of an eigenvalue is defined only up to sign, so the result is *a* square root |
| `is_permutation_of` / `permutation_normal_form` | `BipartiteGraph.is_isomorphic(…, edge_labels=True)` (bliss/nauty) | Complete backtracking pruned by the column-multiset invariant of every prefix plus row/column signature filters. Both outputs are *uniquely specified*, so any complete algorithm agrees; both were checked against exhaustive brute force. Worst case exponential, as is Sage's |
| `pluq` / `ple` | `P` and `Q` from M4RI's `mzp_t`, i.e. transposition lists applied in order | Both are transposition lists (`P[pivotRow] = foundRow`) |
| `matrix_operations.pivot_rows` | Row indices | Row indices |
| Affected modules | `sage/matrix/matrix2.pyx`, `matrix0.pyx`, `matrix_misc.py`, `matrix_integer_dense.pyx`, `matrix_modn_dense_template.pxi` | `packages/sagemath-ts/src/matrix/` |

### Rationale

1. **Exactness over reproduction** — CLAUDE.md forbids floating point where Sage is exact. The
   previous double-precision LLL stopped producing a basis of the input lattice above 2^53.
2. **The delegation targets are missing or broken** — `rational_form` is a stub, `block_ldlt` was
   itself producing invalid factorizations, and there is no graph package.
3. **Provable equivalence** — each substitution computes the same mathematical object by a different
   route, verified against upstream doctests and randomized sweeps.
4. **`norm(2)` is the exception**: upstream is explicitly inexact there, so following it faithfully
   means following it in double precision.

### Trade-offs

- An exact LLL cannot reproduce fpLLL's rounding-dependent choice of representative, so individual
  rows differ, typically by sign: `matrix(ZZ,3,range(1,10)).LLL()` row 1 is `[2,1,0]` in Sage and
  `[-2,-1,0]` here; `matrix(ZZ,[[1,2,3],[31,41,51],[101,201,301]]).LLL()` row 1 is `[-1,0,1]` in Sage
  and `[1,0,-1]` here. The result is always a `(delta, 1/2)`-reduced basis of the same lattice.
- `QR`'s `Q` is not orthonormal, so callers expecting a unitary matrix must normalize.
- `permutation_normal_form(check=true)` may return a different (equally valid) permutation;
  Bunch-Kaufman's permutation may differ from Sage's over finite fields; `principal_square_root`
  returns a non-principal root over finite fields. The returned matrix is identical in each case.
- The generic determinant paths are slower than LinBox/multimodular (values identical).
- `jordan_form` raises `ArithmeticError` where Sage raises `RuntimeError`, and so does `is_similar`
  on the (unreachable) double failure.
- Sage's `jordan_form` doctests over `PolynomialRing(QQ, 'x11,…')` and over
  `FractionField(PolynomialRing(QQ,'a'))` cannot be run: the port cannot build a `Matrix` over a
  multivariate polynomial ring or a rational function field.
- `matkermod` is called with `wantIm = true` even though the image is discarded, to disable an
  unsound PARI shortcut (see [Upstream Behaviour](#upstream-behaviour-deliberately-not-reproduced)).
  Extra Howell work on tall matrices.
- **`Matrix.toString` is not subdivision-aware and pads per column** (`matrix_generic.ts:428`).
  SageMath pads every entry to one global width and draws `|` / `---+---` separators
  (`matrix0.pyx:2180`). A faithful `matrix_str` exists in `matrix_decompositions.ts` but
  `toString` does not delegate to it; `jordan_form` attaches `matrix_str` as a per-instance
  `toString` on the subdivided `J`, which is a stopgap. `_subdivisions` is also not preserved by
  `Matrix.copy()`, whereas SageMath preserves subdivisions under `copy()`.
- `is_hermitian` is exported from `matrix_operations.ts` but **not** re-exported from
  `matrix/index.ts`.

### Mitigation

Once `rational_form` and `block_ldlt` exist, swap `is_similar` and the definiteness predicates for
the upstream paths with no visible change. Rewrite `jordan_decomposition` to Sage's
minimal-polynomial Newton iteration so it works for non-split characteristic polynomials.

### Behavioral Impact

Values match SageMath's, verified by execution: 300 random integer matrices agree with Sage on the
Frobenius form `F`, the elementary divisors and **both halves** of `flag=2`, with `B^-1 F B == A`
exact on all 300; frobenius flags 0/1/2 on `diag(1,1,2)` match character for character including
`B = [[-2,-4,-1],[1,2,1],[1,1,0]]`; 300 random composite-modulus matrices match Sage's kernel basis
entry for entry, and 2190 brute-force cases confirm the returned rows generate the *full* kernel.
The issue-12693 `jordan_form` doctest reproduces Sage's `P = [2 1 0/0 0 1/-2 0 -1]` character for
character; `is_similar`'s transformation doctest reproduces Sage's
`T = [[1,0,0],[-2/3,1/6,-5/6],[2/3,0,-1/3]]`, and the similar/not-similar verdict was checked
against **exhaustive brute force over every invertible `P`** for GF(2) 2×2 (136 pairs), GF(3) 2×2
(3321) and GF(2) 3×3 (131 328) — 0 mismatches.

---

## Matrix Special Constructors

`sage/matrix/special.py` maps to `packages/sagemath-ts/src/matrix/matrix_special.ts`. This is the
single owner for that module's divergences; no other section duplicates them.

| Function | SageMath | sagemath-ts |
|----------|----------|-------------|
| `companion_matrix`, `toeplitz`, `hankel` | See `special.py` | Sage's argument conventions (full monic coefficient list with negated border; `r` counted from the second column with `ncols = len(r)+1`) |
| `elementary_matrix(row1 == row2)` with no scale | `special.py:1512-1516` collapses to `elem[r,r] = 1`, i.e. the identity; Sage raises only when a scale is *also* given | Replicated verbatim — identity for a self-swap, `ValueError` for the two cases Sage rejects. See [Upstream Behaviour](#upstream-behaviour-deliberately-not-reproduced) |
| `block_matrix` | Flat list + `nrows`/`ncols`, ragged list, or list of lists | List of lists only; a ragged one raises Sage's own `ValueError('list of rows is not valid (rows are wrong types or lengths)')` |
| `random_echelonizable_matrix` / `random_unimodular_matrix` `upper_bound` | Size control by rejecting row operations past the bound — **only over ZZ and QQ** | `NotImplementedError`: the port's generic constructors work over any ring with `random_element()` and have no notion of absolute value |
| `random_unitary_matrix`, `vector_on_axis_rotation_matrix`, `ith_to_zero_rotation_matrix` | Implemented over RDF/CDF via QR / Haar measure / trigonometric rotations | `NotImplementedError` naming the requirement (`sqrt` and trigonometric functions over an inexact ring) |
| `hadamard_bound` | Uses `sqrt` in the base ring | `NotImplementedError` for rings without `sqrt` |
| `lehmer`, `hilbert` | Return matrices over QQ | Require a ring argument supporting `__call__` and division; `NotImplementedError` otherwise |
| `rook_vector` | `ButeraPernici` (default), `Ryser`, `Godsil` | Naive placement counting; `NotImplementedError` naming the two faster algorithms once `positions.length > 50` **and** `k > 5` |
| `is_permutation_of` / `permutation_normal_form` | `BipartiteGraph.is_isomorphic(…, edge_labels=True)` | Complete backtracking (see [Matrix Module Algorithm Substitutions](#matrix-module-algorithm-substitutions)). `permutation_normal_form(check=true)` returns 0-based index arrays with the convention `normal_form[i][j] === matrix[row_perm[i]][col_perm[j]]` instead of a pair of 1-based `PermutationGroupElement`s |
| Random matrix constructors | `sage.misc.prandom` (`randint`/`shuffle`) driven by the global randstate | `current_randstate().randint(…)` and a randstate-driven Fisher-Yates shuffle / density fraction. The underlying stream is bit-identical to GMP's, but the **draw order** differs, so values differ from Sage's for the same seed |
| `random_diagonalizable_matrix` | — | `NotImplementedError('unexpected eigenvector layout')` on an internal invariant violation |
| `matrix(...)` constructor | `constructor.pyx` accepts flat lists with `nrows`/`ncols`, dicts, callables, sparse flags and a bare `(nrows, ncols)` form | **List of lists only** — registered as an open gap under [Matrices](#matrices--j-ideals-lll-reducedness-and-the-matrix-constructor) |
| Affected modules | `sage/matrix/special.py`, `constructor.pyx` | `packages/sagemath-ts/src/matrix/matrix_special.ts`, `matrix_space.ts` |

### Rationale

1. **No inexact matrix type** — `random_unitary_matrix` and the rotation constructors are defined
   over RDF/CDF; implementing them with JS doubles inside a `Matrix<R>` would introduce floating
   point into the exact matrix hierarchy for functions nothing in this port consumes.
2. **Refusing an argument beats ignoring it** — `upper_bound` is the clearest case: SageMath's own
   implementation is ZZ/QQ-only and the port's constructors are ring-generic.
3. **Bounded naive algorithms are declared** so the reachable range is documented rather than
   discovered as a hang.

### Trade-offs

- Five constructors are unavailable where SageMath answers.
- `block_matrix`'s flat-list and ragged forms are unavailable.
- `permutation_normal_form(check=true)` may return a different (equally valid) permutation when the
  matrix has non-trivial automorphisms; the matrix returned is identical.

### Behavioral Impact

Values match SageMath's for every implemented constructor. The divergences are honest refusals or
index-base/shape adaptations — none returns a different mathematical object.

---

## Lattice Algorithms — CVP, Voronoi Cells and LLL Representatives

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `IntegerLattice.closest_vector` | Projects `t` onto the span, then Micciancio-Voulgaris over the diamond-cut Voronoi cell | Exact Fincke-Pohst enumeration seeded with Babai's nearest plane (the projection is unnecessary: the orthogonal component of `t` is a constant offset of the objective) |
| `approximate_closest_vector` | `nearest_plane` / `rounding_off` / `embedding` | Same, with Sage's `embedding` default and round-half-to-even |
| `voronoi_relevant_vectors` | Reads the defining point of each inequality of the diamond-cut `Polyhedron` | Enumerates the `2^r − 1` nonzero cosets of `L/2L` and keeps those with exactly two minimal-length vectors (Voronoi's theorem); results sorted lexicographically |
| `voronoi_cell` | Returns a `Polyhedron` | Returns an H-representation `{normals, offsets}` (gcd-normalised `normals · x <= offsets`, `offsets` exact `bigint`) |
| `isLLLReduced` | fpLLL, `delta = 0.99` | Exact bigint, `delta` default 0.99, no `1e-10` fudge (`free_module_integer.ts:1991`). Skips the leading zero rows LLL emits for rank-deficient input, and raises Sage's `'sage'`-algorithm `ValueError('linearly dependent input for module version of Gram-Schmidt')` for genuinely dependent nonzero rows |
| Affected modules | `sage/modules/free_module_integer.py` | `packages/sagemath-ts/src/modules/free_module_integer.ts` |

### Rationale

1. **Exactness** — CVP is exact here. The previous `closestVector` enumerated around the **origin**
   with coefficients in `[-3, 3]`, so for a distant target it degraded to Babai (a rank-3 example
   gave `d^2 = 125` against the true 98).
2. **No `Polyhedron` class** — the facet description is the natural representation of the same
   object; the relevant-vector set is mathematically identical and was verified sound and complete
   against brute force.
3. **Cost** — enumeration is far cheaper than building the Voronoi cell, which is exponential in the
   rank.

### Trade-offs

- When several lattice vectors are equidistant from the target, which one is returned may differ
  from SageMath's.
- `voronoiCell` above rank 24 raises `NotImplementedError` rather than exhausting memory.
- Callers wanting a polyhedron object must build it from the inequalities themselves.
- Exact SVP is capped at rank 30 and **silently approximates** above it — registered as an open gap
  under [Lattices](#lattices--exact-svp-rank-cap).
- The only remaining floating point in the module is the legacy exported `gramSchmidt()` helper
  (used by `bkz.ts` and `discreteGaussianSample`) and the heuristic estimators (`hadamardRatio`,
  `gaussianHeuristic`, `hermiteFactor`, `estimateBKZBlockSize`).

### Behavioral Impact

Results are exact closest vectors, matching SageMath's value though not necessarily its choice among
ties. All five `approximate_closest_vector` doctest values reproduce exactly (delta 0.26 ->
`(1331,1324,1349,1334)`; delta 0.99 and `nearest_plane` -> `(1326,1349,1339,1345)`; `rounding_off`
-> `(1331,1324,1349,1334)`; `(-6,5/3)` -> `(-6,2)`), as does the `voronoi_relevant_vectors` doctest
(`IntegerLattice([[3,0],[4,0]])` -> `[(-1,0),(1,0)]`).

---

## Free Module Exactness and Coordinate Types

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Linear algebra | Exact over any PID, with echelon forms over e.g. `QQ[x]` | Exact fraction-field layer for `bigint`/`number`/`Rational` entries, for fields providing `div`/`inv`, and for univariate polynomial rings over a field (`QQ[x]`, `GF(p)[x]`): a `FractionFieldElement` (exact `QQ(x)`), a port of `Matrix._echelon_form_PID` (`matrix2.pyx:17305`) and `_generic_clear_column` (`:20613`), and a port of Sage's generic Smith normal form (`:16732`, `:20730`, `:20537`). Over any *other* ring, echelon forms, coordinates, kernels and linear dependence raise `NotImplementedError('exact linear algebra is not implemented over this base ring')` |
| `tensor_product` | **No such method on free modules.** The only concrete embedded definition upstream is `free_quadratic_module_integer_symmetric.py:1343`, built on `Matrix.tensor_product` (`matrix2.pyx:9983`) | A port of that: rank `m·n` submodule of `R^(deg1·deg2)` whose user basis is the Kronecker product of the basis matrices, with the Kronecker product of the inner product matrices, plus Sage's `discard_basis=True` variant. (An earlier `@see Reference: FreeModule_generic_pid.tensor_product` citation was **fabricated** — no such method exists) |
| `quotient` over a field | `FreeModule_generic_field.__quotient_matrices` (`free_module.py:5366`) | Ported exactly (basis extension by pivot rows of `B.stack(S)`, `Q = D[:, n−m:n]`, `L = D^-1[n−m:n, :]`) |
| `quotient` over ZZ | `FGP_Module` (`fg_pid/fgp_module.py:268`) | Ported exactly, including `invariants()` with and without ones, the Smith generators and `cardinality()` |
| `quotient` over any other ring | Various | Sage's exact `NotImplementedError('quotients of modules over rings other than fields or ZZ is not fully implemented')` (`free_module.py:4472`) |
| `coordinate_vector` | Element of `FreeModule(R.fraction_field(), rank)` — always a `Rational` for a ZZ module | `bigint` when integral, `Rational` otherwise (`number` over a JS-number base ring) |
| `coordinates(check=False)` | Skips verification and can return a vector that does not reconstruct `v` | Always raises `ArithmeticError('vector is not in free module')` when `v` is outside the span |
| `indexIn()`, `cardinality()` | Base-field element or `infinity`; Sage `Integer` or `+Infinity` | `bigint`/`Rational`, or `Number.POSITIVE_INFINITY`; cardinalities routinely exceed 2^53 (`GF(2)^70`) so `bigint` is required |
| `norm(p)` for irrational results | Symbolic (`sqrt(14)`, `276^(1/5)`) | Exact `bigint`/`Rational` whenever the p-th root is rational, a double otherwise |
| `normalized()` | `v / v.norm(p)`; the base ring changes, usually to the symbolic ring | A vector over QQ, or over the double field when the norm is irrational |
| `discriminant()` | `FreeModule(R,n)` uses `det(gram)`; with an inner product matrix it is a `FreeQuadraticModule` whose discriminant is `(−1)^(rank//2)·det(gram)` | Same split, keyed on whether an inner product matrix was supplied (the port merges `free_quadratic_module.py` into `free_module.ts`, so the class distinction becomes a runtime condition) |
| `submodule(gens, check=True)` | `ArithmeticError('argument gens (= …) does not generate a submodule of self')` | Implemented; `span()` remains unchecked as in Sage |
| `Frac(QQ[x])` normalisation | Keeps unit denominators (`x/2` has numerator `x`, denominator `2`) | Divides the unit out (numerator `1/2·x`, denominator `1`). Purely representational — the two print identically |
| Affected modules | `sage/modules/free_module.py`, `free_quadratic_module.py`, `free_module_element.pyx` | `packages/sagemath-ts/src/modules/free_module.ts`, `free_module_element.ts` |

### Rationale

1. **The alternative was silent double-precision arithmetic** over arbitrary rings, which produced
   the previous defects: `span()`/`subspace()` used the *number of generators* as the rank, and the
   float RREF/kernel/determinant helpers rounded doubles back into `bigint`.
2. **Consumer expectations** — every caller of a ZZ module wants `bigint`; carrying `Rational`
   everywhere would change the entry type of the whole module.
3. **`check=False` cannot be reproduced usefully** — Sage's unchecked partial answer needs the
   rref-pivot transformation machinery and yields a meaningless vector; raising never differs when
   `v` is in the module.
4. **Delegation** — Hermite normal form goes to `matrix_integer.hermite_normal_form` and saturation
   to `matrix_integer.saturation`, exactly where SageMath delegates to `Matrix_integer_dense`.

### Trade-offs

- `intersection()` over a non-ZZ PID returns the mathematically correct module (250/250 random cases
  verified **equal in SageMath**, using Sage's own module equality) but in ~25 % of cases its
  **echelon basis** differs from Sage's by a unit of the base ring. Two upstream normalisations are
  not reproducible from the Euclidean interface: `_echelon_form_PID` is itself not canonical up to
  units, and Sage's `integer_kernel` (`matrix2.pyx:5646`) scales by `Matrix.denominator()` — for a
  `QQ[x]` matrix the lcm of the *coefficient* denominators, a notion living in `Frac(ZZ)` rather than
  in the Euclidean structure of `QQ[x]`. The module generated is always identical.
- `_echelon_form_PID` omits the reduction above the pivots (`matrix2.pyx:17419-17426`), for the same
  reason Sage omits it for `K[x]`: polynomial ideals have no `small_residue`, so Sage's own
  `except AttributeError` swallows that step.
- `FreeModuleQuotient` folds Sage's two classes (`FreeModule_ambient_field_quotient` and the
  torsion-carrying `FGP_Module`) into one, so over a field `rank = degree = dim V − dim W` and over
  ZZ `degree` is the number of Smith invariants with `rank` the free rank.
- `tensorProduct`'s degree is `deg(M)·deg(N)` (equal to `rank·rank` for ambient modules) and the
  result has a basis rather than being ambient.

### Behavioral Impact

Ranks, echelon bases, coordinates, kernels, intersections, complements, discriminants and
cardinalities are exact and match SageMath's values: 700 random `QQ[x]`/`GF(p)[x]` spans match the
vendored SageMath's echelon basis exactly; quotients were checked against the installed Sage on 200
random QQ and 191 random ZZ cases — **every** projection, lift, invariant and cardinality agrees
exactly, with `project(lift(x)) == x`, `project(W) == 0` and additivity executed on all 391; the
`IntegralLattice("D3")` tensor-product doctest reproduces number for number including the Gram
matrix and the `discard_basis` variant; 914 randomly generated `QQ[x]`/`GF(p)[x]` cases across five
sweeps agree coefficient-for-coefficient on `P`, `Q`, `P ∩ Q` and `P + Q`.

---

## Binary Quadratic Forms

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Composition / reduction backend | `BinaryQF.__mul__` calls PARI `qfbcompraw`; `reduced_form` calls PARI `qfbred`/`qfbredsl2` for non-square discriminants | **Delegated** to `parigp-ts`'s `qfb.ts` (a port of `pari/src/basemath/Qfb.c` and `quad.c`). Sage's own `_reduce_indef` is retained for **square** discriminants only |
| `algorithm` selection | `algorithm = 'sage' if self.is_reducible() else 'pari'` (`binary_qf.py:947-948`), with `'default'`/`'pari'`/`'sage'` accepted explicitly | Identical, with all three of Sage's error paths |
| Squaring dispatch | PARI's `qfb_comp` squares only when the two GEN **pointers** are identical (`if (x == y)`); Sage converts both operands separately, so that path **never fires from Sage**, not even for `Q * Q` | The `this === other -> _square()` shortcut was removed to match. Behaviourally a no-op: `qfb_sqr` and the general `qfb_comp` agreed on all 716 self-compositions tested |
| `solve_integer` | `binary_qf.py:1608-1806`: negative-definite recursion, an elementary algorithm for square discriminants, `qfbcornacchia` for prime `n` with `disc < 0`, else `qfbsolve` with `_flag` in {1,2,3} | Ported in full; Sage's `Factorization` argument becomes an optional `{ factorization }` option |
| Shanks distance forms (`qfr5_*`, `qfr5_dist`) | Present | Ported. `qfb.ts` carries a transcription of PARI's `t_REAL` kernel (`nbits2prec`, `addrr`, `mulrr`, `divrr`, `sqrtr`, `mplog2`, `logr_abs`, `shiftr`, …) and a `QfbExt` type carrying the logarithmic distance. `buch.ts` carries a **second, independent** copy of that kernel; only `qfb.ts`'s is re-exported from the package root |
| Reduction arithmetic | `D.sqrt(prec=53)` | Exact `isqrt(D)` (see [Exact Arithmetic](#exact-arithmetic-where-sagemath-uses-floating-point)) |
| Affected modules | `sage/quadratic_forms/binary_qf.py` | `packages/sagemath-ts/src/quadratic_forms/binary_qf.ts`, `packages/parigp-ts/src/qfb.ts` |

### Rationale

`qfb.ts` was verified against the **real PARI 2.15.4** (reached through the local SageMath) on
golden data: `qfbredsl2` 80/80 exact (form *and* base-change matrix), `qfbcompraw` 300/300,
`qfbpowraw` 120/120, `primeform` 60/60, `qfbcornacchia` 104/104, `qfbsolve` 1670/1800. Every one of
the remaining divergences is a **documented upstream change** between the oracle (2.15.4) and the
source we ported (2.18.1), not an error: `CHANGES-2.16 #45` ("changed `qfbred` to use standard
normalization, same as `qfbredsl2`") accounts for the 6 `qfbred` + 14 `qfbcomp` + 4 `qfbpow` cases,
and `CHANGES-2.16 #9` plus `allsols`' new `v >= 0` normalisation and the lexsort at `Qfb.c:1930`
account for the 130 `qfbsolve` cases. The 2.18 semantics were verified independently rather than
assumed.

### Trade-offs

- **`reduced_form` for indefinite non-square discriminants follows the vendored 2.18.1
  normalisation**, i.e. the same representative as `qfbredsl2`. Code validated against a PARI older
  than 2.16 may see a different (equally reduced) indefinite representative. Property tests
  therefore compare the sorted `cycle(proper=True)` class invariant plus the exactly-comparable
  `reduced_form(transformation=True)` triple, rather than the raw representative.
- `solve_integer` on non-square discriminants: PARI 2.15.4 and 2.18.1 return different (both
  correct) representations — the vendored Sage doctest itself changed sign between versions — so the
  property tests compare an exhaustive brute-force solution **set** (complete for positive definite
  forms, since `|y|` is bounded by `4an/|D|`) plus membership of PARI's answer in it. That is
  strictly stronger than pinning one pair.
- `solve_integer` for hard-to-factor `n` inherits `ifactor.ts`'s factoring chain; the optional
  `factorization` argument is the documented workaround.
- `algorithm` is an options-object field rather than a positional keyword.

### Behavioral Impact

Equivalence with the pre-delegation code was proven by execution against a side-by-side import of
the previous file: 29 944 random forms (coefficients to 10^12) with 0 differences in both the
reduced form and the SL2 base change; 59 280 compositions across all 400 valid discriminants in
`[-400,400]` with 0 differences; 1600 `BinaryQF_reduced_representatives` calls with 0 differences;
500 class-group Cayley tables with 0 differing tables. `solve_integer` reproduces every Sage doctest
exactly, and its square-discriminant branch was cross-checked against exhaustive brute force over a
241×241 box on 600 random cases with 0 wrong solutions and 0 false nulls.

---
## Quadratic Forms (sage.quadratic_forms)

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `Q.matrix()` / `Q.Gram_matrix()` | A matrix over the form's base ring (`Integer Ring` for a ZZ form) | A `RationalMatrix` — a structurally typed view of the generic `Matrix` specialised to `Rational`, always with QQ as its runtime base ring. Entries are unchanged (and integral where SageMath's are). Only `A.base_ring()` differs, which no ported code consults |
| `det()`, `Gram_det()`, `coefficients()` | `Integer` for a ZZ form, `Rational` for a QQ form | Always `Rational` (use `.numerator` for the bigint). `gcd()` and `level()` still return `bigint` because they are ZZ-only |
| One-argument matrix constructor | `QuadraticForm(M)` takes `M.base_ring()` | Infers ZZ when every entry is integral, QQ otherwise. `adjoint_primitive` depends on the ZZ inference (it calls `.primitive()`, which is ZZ-only), so this matches SageMath where it matters |
| `theta_series` / `theta_by_pari` | A power series in `ZZ[[q]]` when `var_str` is nonempty (the default), the raw vector when `var_str == ''` | Always the `bigint[]` coefficient vector, i.e. Sage's `var_str == ''` behaviour |
| `pseudorandom_primitive_zero_mod_p` | Draws `(r1, r2)` uniformly at random until it finds a zero | Scans `r1 = 0,1,…` and `r2 = 0,1,…` in order |
| `has_integral_Gram_matrix`, `level` warnings | `warnings.warn` | `console.warn` with SageMath's exact message text |
| `qfgaussred` | `self.__pari__().qfgaussred()` | A line-for-line port of PARI's `gaussred` (`alglin2.c:1650-1749`) inside `quadratic_form__local_field_invariants.ts`, because `parigp-ts` exports only `qfgaussred_positive`. Its output is pinned against real PARI in a test |
| `QuaternionOrder.quadratic_form()`, `QuaternionFractionalIdeal.quadratic_form()`, `ternary_quadratic_form()` | `sage.quadratic_forms.QuadraticForm` objects | The underlying Gram/Hessian matrix (`IntegerMatrix` for the quaternary forms — denominator-cleared and divided by the gcd exactly as upstream — and a `RationalMatrix` for the ternary form). Numerically identical to `Q.matrix()`; `theta_series`/`theta_series_vector`/`minimal_element` are provided directly on the ideal |
| Affected modules | `sage/quadratic_forms/` | `packages/sagemath-ts/src/quadratic_forms/` |

### Rationale

1. **Type-system forced** — the repo's `Matrix<R extends RingElement>` constraint is not satisfiable
   by `Rational` under `strict` (its `add` accepts `Rational | IntegerLike`, not exactly `this`), and
   a `bigint | Rational` union would infect every caller.
2. **Reproducibility** — a deterministic scan for primitive zeros makes `find_zeros_mod_p`
   reproducible. As a *set* of points of `P^2(F_p)` the output is independent of the starting zero
   (verified against Sage for 20 forms at `p = 7` and one at `p = 17`), and given the **same**
   starting vector `_find_zeros_mod_p_odd` reproduces Sage's list element-for-element (verified for
   6 cases including the `p = 1009` doctest).
3. **No `PowerSeriesRing` dependency** is worth taking for theta series when the coefficients are
   what every doctest compares.
4. **Upstream cache invalidation** — see [Upstream Behaviour](#upstream-behaviour-deliberately-not-reproduced)
   for `__setitem__`.

### Trade-offs

- `sage.quadratic_forms.QuadraticForm` is not available as a return type from the quaternion
  modules, so callers get a Gram matrix and must build the form themselves.
- Callers wanting a printable theta series must build it from the coefficient vector.
- `qfgaussred` lives in `sagemath-ts` rather than `parigp-ts` — see
  [PARI/NTL Routines](#parintl-routines-duplicated-or-ported-in-place).

### Behavioral Impact

Values are SageMath's; only wrapper types differ. `find_zeros_mod_p`'s *order* and choice of
representatives depend on the deterministic start.

---

## Elliptic Curves and Isogenies

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Weierstrass models | General form in all characteristics | Five-coefficient construction now preserves the supplied equation through EllipticCurveGeneric, including characteristics 2/3. Two-coefficient prime curves retain EllipticCurveFiniteField. The generic finite-field PARI group backend remains incomplete; see Constructor model routing below. |
| `is_j_supersingular` | Checks `supersingular_j_polynomial(p)(j) == 0` when `p` is in the precomputed table, giving an exact answer even with `proof=False` | Skips the table (`supersingular_j_polynomial` is not ported) and always falls through to the 10 random-point tests (`ell_finite_field.ts:1463`), plus the trace-of-Frobenius check when `proof` is set (the default). With `proof=True` — Sage's and our default — the answer is identical and proved |
| `montgomery_model` representative | `EllipticCurveIsogeny(GF(7) j=1728 curve, (0,0), model='montgomery')` reports `A = 1` | Returns `A = 6`, the other root of the defining cubic. Both are valid Montgomery forms; see the root-ordering row under [Polynomial Roots](#polynomial-roots-and-factorization) |
| `possible_isogeny_degrees(E)` over Q | Billerey/Larson bounds | Mazur's list `[2,3,5,7,11,13,17,19,37,43,67,163]`, optionally intersected with the degrees for which `isogenies_prime_degree` finds an isogeny. Correct as a **superset** over Q; **not valid over larger number fields** |
| `isogeny_degrees_cm(E)` | Exact | Ported including the horizontal-primes step (`isogeny_class.py:1309-1317`) and the `n/(2h)` downward-ramified test. The function's contract ("this list is not necessarily minimal") holds |
| `Frobenius_filter` good-reduction test | `E.has_good_reduction(p)` on the **minimal** model (Tate's algorithm via `local_data`) | `v_P(disc) == 0` of the *global integral* model built in place (Laska-Kraus-Connell minimisation is not ported). A non-minimal model makes a few extra primes look bad; those are skipped, which can only make the filter **weaker** (a superset), never unsound |
| `isogenies(fill=true)` | `isogeny_class.py:369-370` raises `NotImplementedError` | The same — this is **not** a deviation, it is upstream's behaviour, recorded here so it is not re-reported |
| `qf_matrix()` | `ValueError('qf_matrix only defined for isogeny classes with rational CM')` (`isogeny_class.py:329-330`) | The same `ValueError` (`isogeny_class.ts:340-345`) |
| Affected modules | `sage/schemes/elliptic_curves/` | `packages/sagemath-ts/src/schemes/elliptic_curves/` |

### Rationale

1. **Unported dependencies** — `supersingular_j_polynomial`, Laska-Kraus-Connell minimisation.
3. **Ties** — several equally valid representatives exist for the Montgomery model.
4. **Superset over unsound** — an over-reported isogeny-degree candidate set is a documented weakening
   of a filter, never a wrong answer.

### Trade-offs

- `is_j_supersingular(proof=False)` is probabilistic where Sage would be exact for small `p`.
- The Montgomery `A` can differ from Sage's printed value.
- `possible_isogeny_degrees` over a number field is not a valid bound.
- `ell_generic`'s `toString` prints `y^2 + 1*x*y = x^3 + 1` where Sage prints `y^2 + x*y = x^3 + 1`;
  `_equation_string` omits Sage's `±1` special cases and the final `s.replace('+ -', '- ')`. This is
  a **bug**, not a deviation — see [Elliptic Curves over Q](#elliptic-curves-over-q-and-number-fields).

### Behavioral Impact

Vélu's formulas, `division_points`, `multiplication_by_m`, `_isomorphisms` (all char 2/3/p branches),
`lift_x`/`is_x_coord`, the bivariate `division_polynomial`, `abelian_group`/`gens`,
`set_order`/`has_order`, `torsion_basis`, `twists`, `frobenius_order` and the whole formal group
reproduce SageMath's doctests. `Frobenius_filter` reproduces all three SageMath doctests exactly,
including the `d = −23` degree-6 verbose transcript ending `List of primes after filtering: [2, 3]`;
189 curves over six number fields (degrees 2, 2, 2, 3, 4, 6) were cross-checked against live
SageMath on **both** the filter output and the `include_2` boolean with 0 mismatches. Note the
correct Vélu accumulation is upstream's `v += vQ`, `w += uQ + xQ·vQ` (with `vQ = 2·gxQ` for
non-2-torsion); the `w += 2(uQ + xQ·gQx)` form double-counts `uQ` and yields a non-isogenous
codomain.

---

## Hyperelliptic Curves and Jacobians

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Integer polynomials | `frobenius_polynomial()` returns an element of `ZZ['x']`; `zeta_function()` a rational function in `ZZ(x)` | `ZZPoly = bigint[]` (ascending coefficients) and a `ZZRationalFunction` with numerator/denominator coefficient arrays. `zz_poly_repr()` renders either exactly the way Sage prints it, and the tests compare those strings |
| `clebsch_to_igusa`, `igusa_to_clebsch`, `ubs`, `Ueberschiebung`, `diffxy` | `clebsch_to_igusa(A, B, C, D)` — ring elements carry their parent, so constants like `-120` and `1/135000` coerce automatically | Take the base ring as an extra first argument: `clebsch_to_igusa(K, A, B, C, D)`. This port's `RingElement` has no `parent`, so there is no way to build the ring's 1 from an element alone. The high-level entry points that take a polynomial (`clebsch_invariants`, `igusa_clebsch_invariants`, `absolute_igusa_invariants_wamelen/kohel`) keep Sage's exact signature because they can read `f.parent.base_ring` |
| Invariants representation | `invariants.py` builds the differential operator symbolically in `QQ[dfdx, dfdy, dgdx, dgdy]` (`differential_operator`) and applies it with `diffsymb`; forms are `MPolynomial`s | Every object in `ubs` is a homogeneous binary form stored as a dense coefficient array of its nominal degree, and `(fx·gy − fy·gx)^k` is expanded with the binomial theorem directly. `differential_operator` and `diffsymb` have no separate counterpart. Transvectants of homogeneous forms are homogeneous, so the dense representation is exact |
| `cantor_reduction_simple` ambiguous form | Asserts `deg == genus+1`, **prints** `Returning ambiguous form of degree genus+1.` to stdout, and returns the pair (`jacobian_morphism.py:151-155`) | Same assertion (raised as a `ValueError` with the observed degree) and same return value, but nothing is written to stdout |
| Class specialisation | `constructor.py:335-368` builds the concrete class at run time by multiple inheritance from `HyperellipticCurve_g2` and the base-ring specialisation | Six explicit classes: `HyperellipticCurve_generic`, `_finite_field`, `_rational_field`, `_g2`, `_g2_FiniteField`, `_g2_RationalField`, with the genus-2 method bodies delegating to free functions in `hyperelliptic_g2.ts`. TypeScript has no multiple inheritance. The class **name** differs from Sage's `HyperellipticCurve_g2_FiniteField_with_category`; `instanceof` against every layer is asserted in `constructor.test.ts` |
| `_points_fast_sqrt` over an extension field | Iterates `GF(p^n)` in Zech-logarithm (givaro) order and uses PARI/givaro's canonical square root | Iterates in integer-representation order and uses this port's square root, so the same 7 / 31 / 122 points come out **permuted**. Neither the field's iteration order nor PARI's sqrt branch is reproducible without porting givaro. Over prime fields, where `points()` uses `_points_cache_sqrt`, the order matches Sage element-for-element including the "larger square root first" quirk |
| `field_embedding` for a non-prime base field | `Hom(K, L)[0]` (`hyperelliptic_finite_field.py:1311`) | Finds the first root of `K`'s modulus while iterating `L` and maps `sum c_i a^i -> sum c_i r^i`. Two embeddings of `GF(p^m)` into `GF(p^(mn))` differ by an automorphism of `L`, and the point count is invariant under it. Verified for `GF(9)` (`n <= 3`) and `GF(4)` (`n <= 6`) |
| Frobenius precision bounds | 53-bit `RR` | Exact in ZZ — see [Exact Arithmetic](#exact-arithmetic-where-sagemath-uses-floating-point) |
| Affected modules | `sage/schemes/hyperelliptic_curves/` | `packages/sagemath-ts/src/schemes/hyperelliptic_curves/` |

### Rationale

1. **No `parent` on `RingElement`** — the extra ring argument is the only way to build ring constants
   from a bare element; the polynomial-taking entry points do not need it and keep Sage's signature.
2. **Homogeneity is preserved by transvection**, so the dense binary-form representation is exact and
   avoids porting a 4-variable polynomial ring.
3. **A library should not print to stdout**, and this repo has no verbosity mechanism.
4. **TypeScript has no multiple inheritance**, so `dynamic_class` becomes a fixed hierarchy.

### Trade-offs

- Consumers cannot do polynomial arithmetic on `frobenius_polynomial()`'s result without converting.
- `differential_operator`/`diffsymb` have no directly callable counterpart.
- A caller cannot see that an ambiguous Cantor form was returned.
- Adding a new base-ring specialisation requires a new genus-2 combination class.
- Point **list order** over non-prime fields differs from Sage's; the point *set* is identical
  (verified for `GF(9)` twice, `GF(49)`, `GF(121)`).

### Behavioral Impact

Values match Sage: the full `ubs` dictionary over `GF(31)` and every invariant tuple over QQ match
coefficient by coefficient; all ten Frobenius precision-bound doctest values over `GF(37)`,
`GF(next_prime(10^9))` and `GF(11)` match. The cost deviation in `frobenius_polynomial`'s algorithm
selection is registered separately under
[Hyperelliptic](#hyperelliptic--frobenius-polynomial-algorithms).

---

## Quaternion Algebras

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `quadratic_form()` return type | `sage.quadratic_forms.QuadraticForm` | The underlying Gram/Hessian matrix; see [Quadratic Forms](#quadratic-forms-sagequadratic_forms) |
| `free_module()` | `sage.modules.free_module` objects (`FreeModule_submodule_with_basis` over ZZ inside `QQ^4`) supporting span/intersection/index_in/quotients | A `ZZLattice` holding the echelon (Hermite) basis, with `basis()`, `basis_matrix()`, `eq`, `is_submodule`, `contains`, `intersection`, `add`, `scale` and `index_in`. `QuaternionAlgebra_abstract.free_module()` returns `{ rank: 4, inner_product_matrix }`. The echelon basis matrices are identical to SageMath's (verified on 30+ ideals), so all containment/equality/index results agree; quotient-invariant computation `(V/W).invariants()` is not provided |
| Infinite places | `ramified_places(inf=True)` returns (finite places, `[ring morphisms QQ -> RR]`) | `[bigint[], string[]]` where the unique infinite place of QQ is the exported constant `INFINITE_PLACE_QQ = 'infinity'`. The set semantics that `is_division_algebra`/`is_matrix_ring`/`is_isomorphic` rely on are preserved exactly |
| `P1List` for prime level | `cyclic_right_subideals` uses `sage.modular.modsym.p1list.P1List(p)` | A private `p1list(p)` / `p1_normalize(p, u, v)` valid for prime `p` (the only case `cyclic_right_subideals` supports). Reproduces `P1List(3).list()` and `P1List(5).list()` and, end to end, the exact order of the returned subideals for `p = 3, 5, 7, 13` |
| `intersection_of_row_modules_over_ZZ` | `s.right_kernel_matrix(algorithm='pari', basis='computed')`, whose particular ZZ-basis appears in the doctest matrix | The repo's Smith-form kernel, so the returned 4×4 matrix can differ by a unimodular factor. The **row module** — the only thing every caller uses, via `ZZLattice.span` — is identical |
| Pickling, Magma, hashing | `unpickle_QuaternionAlgebra_v0`, `_magma_init_`, `__hash__` | Omitted; algebra identity comes from the factory cache and equality from `eq()` on canonical lattices |
| Affected modules | `sage/algebras/quatalg/` | `packages/sagemath-ts/src/algebras/quatalg/` |

### Rationale

1. **The dependencies are not ported** — `sage.quadratic_forms.QuadraticForm`,
   `sage.modules.free_module` in a form usable for ZZ-spans of rational vectors, ring morphisms and
   `sage.modular`.
2. **The lattice is what every consumer uses**, and it is identical; only the presentation differs.

### Trade-offs

- No `Factorization`, `QuadraticForm` or `FreeModule` objects to hand back.
- The raw matrix printed by `intersection_of_row_modules_over_ZZ` differs from the doctest's.
- No pickling, Magma interface or hashing.

### Behavioral Impact

Every lattice, norm, order, theta series and equivalence result agrees with Sage. Base rings other
than QQ are an open gap — see
[Quaternion Algebras — Base Rings](#quaternion-algebras--base-rings-other-than-qq).

---

## Function Fields

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Element and ideal comparison over small odd prime constant fields | `Frac(GF(p)[x])` for small odd `p` uses `FpTElement._richcmp_`, which compares `(numerator, denominator)` lexicographically and whose own docstring says "the ordering is arbitrary" (`fraction_field_FpT.pyx:376`). GF(2), GF(4), GF(65537), GF(131101) and QQ use the generic `richcmp(a.num·b.den, a.den·b.num)` (`fraction_field_element.pyx:994`) | The same backend selection: FpT for prime constant fields with `2 < p < 46341`, otherwise generic cross multiplication |
| `element.factor()` / `ideal.factor()` | A `Factorization` with a `.unit()`; `Factorization(self._factor(), cr=True)` | `{ unit, factors }` and a sorted array of `[ideal, exponent]`. Both reproduce `Factorization.sort`'s ordering: `(degree, exponent, prime)` for elements and the prime itself for ideals (`factorization.py:671`) |
| `place.residue_field()`, `valuation_ring.residue_field()`, `divisor.function_space()` | `(k, mor_from_k, mor_to_k)` as `FunctionFieldRingMorphism` objects; `(V, mor_from_V, mor_to_V)` with a `VectorSpace` | `[k, from_k, to_k]` and `[dimension, from_V, to_V]` with plain functions. Residue and function-space triples are cached and public callable maps preserve the compared native domain-conversion rules and errors; morphism introspection is unavailable. `to_V(from_V(v)) == v` is tested |
| `RationalFunctionField.field()`, `element()` | `Frac(k[x])`; the underlying `FractionFieldElement` | The actual cached `Frac(k[x])` parent and underlying fraction element; identity, normalization state, representations and arithmetic now delegate |
| `FunctionFieldMaximalOrderInfinite` | Multiple inheritance from `FunctionFieldMaximalOrder` and `FunctionFieldOrderInfinite` | Extends `FunctionFieldOrderInfinite` only and overrides `_repr_`. That is the branch `FunctionFieldPlace._richcmp_` depends on (`isinstance(prime.ring(), FunctionFieldOrderInfinite)`, `place.py:166`), so place ordering is preserved |
| `_place_class` | A class attribute set in `RationalFunctionField.__init__`, consumed by `PlaceSet` as `self.Element` | An abstract factory method `_place_class(parent, prime)`. TypeScript has no `element_class` mechanism; the constructed objects are identical |
| `ConstantField` / `ConstantFieldElement` | Any object in `Fields()`; the coercion framework does the rest | A new `constant_field.ts` (no upstream counterpart) declaring the structural interface a constant field must satisfy, plus `constant_field_characteristic/_is_finite/_cardinality/_element_list`, `compare_constants` and `divide_constants`, tolerating `characteristic` as a property or a method |
| `IndexError` | `K.gen(1)` raises `IndexError("Only one generator.")`; `O.gen(1)` raises `IndexError("there is only one generator")` | The shared `IndexError` class and the exact SageMath message |
| Valuations | `+Infinity` | `bigint \| number` with `Number.POSITIVE_INFINITY`; see [Infinity Representation](#infinity-representation) |
| Affected modules | `sage/rings/function_field/` | `packages/sagemath-ts/src/rings/function_field/` |

### Rationale

1. **Backend selection is observable.** The port now reproduces FpT's degree-first,
   low-coefficient-first numerator/denominator ordering below the native prime cutoff,
   and the generic cross-multiplication rule elsewhere. Rational constants compare
   numerically.
2. **Explicit map and container adapters** — the public function-field API retains
   callable maps, a vector-space dimension and factor arrays. The broader repository
   has vector-space classes, but this adapter does not expose those objects or
   morphism introspection.
3. **No multiple inheritance in TypeScript**; the branch that carries place ordering was kept.

### Trade-offs

- The former small-prime fraction-ordering deviation is repaired. The expanded live
  comparisons cover QQ, characteristic two, and primes on both sides of the FpT cutoff.
- No `Factorization` repr and no `.prod()`; no morphism domain/codomain introspection.
- The former `field()`/`element()` representation deviation is repaired. Native
  square-root dependency routing and broader coercion support remain under audit.
- `compare_constants` falls back to string comparison for constant fields whose elements expose no
  integer lift or a rational numeric comparator, which can order other constant fields
  differently from SageMath. Prime fields and QQ have direct numerical comparisons.
- The `_place_class` factory and `constant_field.ts` are architectural and arguably belong in
  `DESIGN.md`; they are listed here because `constant_field.ts` has no mirrored upstream path.

The original polynomial factorization zero guard is now propagated through function
field elements and finite ideals. Installed Sage 10.3 lacks this guard on word-modular
polynomials; the oracle applies the bundled `polynomial_zmod_flint.pyx:820` guard
before calling the older factorization code. QQ and extension-field zero errors also
match the installed runtime directly. Infinite-ideal factorization uses degree
arithmetic and retains the source's zero-ideal behavior.

Factory uniqueness now holds for the same constant-field parent object and variable
name. Renaming returns the native triple using plain callable maps; domain/codomain
introspection is absent. The port's distinct-but-equivalent constant-field objects
are not canonicalized by this factory.

For negative powers of zero, the source `is_nth_power` is a Cython `noexcept` method:
inversion failure emits an unraisable diagnostic and returns false. The port returns
false for this case without reproducing the runtime's stderr diagnostic. Ordinary
inversion, division and root errors retain the selected fraction backend's class and
message. Square-root lists deduplicate the characteristic-two root.

### Behavioral Impact

Places, divisors, Riemann-Roch spaces and factorization content and order are SageMath's, verified by
byte-for-byte transcript diffing including a dedicated 602-line GF(65537)/GF(131101) run.

---

## Power Series, Laurent Series and Multivariate Series

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Parent identity | Compares parents with `is` (power series rings are `UniqueRepresentation`) | `PowerSeriesRing.is_identical_to`: same base ring, variable name and default precision. There is no parent cache, so `===` would spuriously re-coerce series built from an equal-but-distinct ring object — exactly the situation in `formal_group.ts` |
| `PowerSeries.__call__` on a zero argument | Returns `self[0]`, an element of the **base ring** | Returns the constant series. TypeScript needs one return type; the value is the same |
| Laurent ordering (`<`, `>`) | `_richcmp_` implements all six comparisons | Equality only. The port's `RingElement` interface has no order on coefficients, and adding one would change an interface implemented across the repo |
| Negative left shift | `__lshift__(n)` sets the precision to `prec + n`, which for `n < 0` can be **negative** | Routed through `__rshift__`, i.e. clamped at 0. No ported upstream path uses `<<` with a negative shift; `>>` with a negative argument (which upstream does use) is exact |
| `sqrt` over `ZZ` | Computes `half = ~R(2)` up front, which for `R = ZZ` silently lands in `QQ` | `1/2` is formed only when a coefficient needs it, so an exact square with unit constant term has a square root over `ZZ` (needed for the `(t^-4).is_square()` doctest); a genuinely fractional root still raises |
| `log()` of a series with a positive non-unit constant term | `(2+x).log()` computes `QQ(2).log()` in the Symbolic Ring and then dies with a `TypeError` adding a symbolic constant to a power series (`power_series_ring_element.pyx:2685-2690`) | `ArithmeticError: constant term of power series is not 1` |
| `MPowerSeries` representation | Wraps a univariate series in an auxiliary variable over the multivariate polynomial ring (`_bg_value`), whose `T`-degree is the total degree | The exponent-vector map directly, plus a total-degree precision, with the precision rules the background ring would produce (`add`: min; `mul`: `min(p1+v2, p2+v1)`). There is no multivariate polynomial ring wired into this module, and the background representation *is* a total-degree grading of the same dictionary |
| `MPowerSeries` division by a non-unit | Falls back on `quo_rem` (documented `# needs sage.libs.singular`) | Throws. The unit case (upstream's fast path) and `ZeroDivisionError` are implemented |
| `laurent_polynomial()`, `_latex_`, `__pari__`, `_im_gens_`, slicing, `_unsafe_mutate`; ring-level `random_element`, `construction`, `polynomial_ring`, `fraction_field` | Present | Absent — they need rings the port does not have, belong in `parigp-ts`, or deliberately mutate an immutable element |
| `MPowerSeries` analytic methods (`exp`, `log`, `derivative`, `integral`, `quo_rem`, `V`, `shift`, …) | Present | Absent; the parts needed to state and check the formal group's associativity identity are complete |
| Affected modules | `sage/rings/power_series_ring*.py`, `laurent_series_ring*.py`, `multi_power_series_ring_element.py` | `packages/sagemath-ts/src/rings/` |

### Rationale

1. **No symbolic ring**, so `log()` of a non-unit constant term cannot reach upstream's crash; a
   clean rejection is strictly better than replicating a `TypeError`.
2. **No parent cache**, so structural identity is the only workable comparison.
3. **The total-degree grading is the same dictionary** the background representation produces.

### Trade-offs

- Laurent series cannot be ordered.
- `MPowerSeries` division by a non-unit and the analytic methods are unavailable.
- **Two performance workarounds are open, and are performance rather than fidelity:**
  `MPowerSeries.inv()` does not match upstream's precision (upstream inverts the *background
  univariate* series, `multi_power_series_ring_element.py:725`), which made one division in
  `group_law(50)` take 12.1 s; and `_subs_formal`'s untruncated intermediate powers made
  `mult_by_n(10, 50)` take 13.7 s. `formal_group.ts` carries a local `bivariateInverse` and a
  truncating composition that reproduce upstream's precision exactly (0.2 s and 0.3 s). The right fix
  is in `power_series_ring.ts`.
- `V(0)` diverges in a way that is **not** accepted — see
  [Power Series](#power-series--v0).

### Behavioral Impact

The formal group's `x(10)` and `y(10)` print exactly as SageMath's doctests do and satisfy the
Weierstrass relation under Laurent arithmetic; `mult_by_n`'s characteristic-zero branch
(`formal_group.py:644-665`) reproduces the 37a doctest character for character; and Sage's whole
`group_law` TESTS block over `GF(7)[[x,y,z]]` — including the three-variable associativity
`F(x, F(y,z)) == F(F(x,y), z)` — is verified.

---

## Coding Theory

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `GoppaCode` decoding | Registers **no decoder** (only a `GoppaCodeEncoder`); Sage falls back to generic syndrome decoding from `AbstractLinearCode` | Binary Goppa uses Patterson; non-binary uses the Sugiyama key equation via `_partialXGCDBalanced`, a direct port of `GRSKeyEquationSyndromeDecoder._partial_xgcd` (`grs_code.py:2145-2157`) followed by Forney `e_i = omega(L_i)/sigma'(L_i)`. Radius `floor(deg(g)/2)`, less than a generic syndrome decoder's true covering radius, but exact within it and failing loudly outside |
| `GoppaCode.distance_bound()` | `1 + deg(g)` | Same |
| `BCHCode.minimum_distance()` | Inherits `AbstractLinearCode.minimum_distance`, delegating to GAP/Guava's Brouwer-Zimmermann | Exhaustive codeword-weight enumeration, cached; `NotImplementedError` once `q^k > 2^17` (`bch_code.ts:698`). Exact wherever it answers (Golay `[23,12]` -> 7) |
| BCH field embedding | `finite_field_base.extension` uses `alpha = E.gen()^((\|E\|−1)/(\|F\|−1))` when both fields are Conway, else `self.modulus().any_root(E)` (randomized Cantor-Zassenhaus) | Tries the Conway power first and accepts it if it is a root of the base modulus; otherwise iterates over the splitting field looking for a root; `NotImplementedError` when `\|E\| > 2^22`. Deterministic, and exactly Sage's choice whenever Conway polynomials are tabulated |
| Affected modules | `sage/coding/{goppa,bch,grs,reed_muller}_code.py` | `packages/sagemath-ts/src/coding/` |

### Rationale

1. **No upstream Goppa decoder to port** — the algorithm was taken from the upstream decoder for the
   closest relative (the GRS key equation) rather than invented.
2. **No Brouwer-Zimmermann port** — the choice was between an exact-but-limited algorithm and an
   honest stub; enumeration reproduces Sage's value where it answers and fails loudly otherwise.
3. **No `any_root`/Cantor-Zassenhaus and no `is_conway` flag** — testing the Conway power directly is
   both deterministic and exactly Sage's choice in the tabulated case.

### Trade-offs

- Non-binary Goppa corrects fewer errors than a generic syndrome decoder would.
- Large BCH codes get `NotImplementedError` where Sage answers.
- Very large non-Conway splitting fields make the BCH embedding unavailable rather than randomized.
- `decode()` naming and several permissiveness differences are registered as open gaps under
  [Coding and Crypto](#coding-and-crypto--permissive-where-upstream-raises).

### Behavioral Impact

Sage's three `GRSKeyEquationSyndromeDecoder` doctests reproduce exactly; the BCH generator polynomial
divides `x^n − 1` over GF(4)/GF(8)/GF(16); Forney carries the `X_i^(l−b)` factor so `b ∈ {0,1,2,3}`
and `l ∈ {1,5,7}` all decode; Reed-Muller's recursive Plotkin decoder decodes `u` from both halves
and keeps the closer candidate (0 failures over 1.3M decodes), and its monomial order matches Sage's
`Subsets` enumeration.

---

## Crypto Module

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `gen_lattice` seeded output, `type='modular'` / `'random'` / dual | `sage.crypto.gen_lattice(m=10, seed=42)` etc. print specific matrices | **Reproduced exactly, row for row.** Sage draws the random block with `MatrixSpace(ZZ_q, m−n, n).random_element()`, i.e. `Matrix_modn_dense_template.randomize` = `rstate.c_random() % p` row-major (`matrix_modn_dense_template.pxi:2843`); we use `c_random()` rather than `mpz_urandomm`. The three doctests at `sage/crypto/lattice.py:81-105` and `:147-157` are pinned |
| `gen_lattice` seeded output, `type='ideal'` / `'cyclotomic'` | Prints a specific matrix | Does **not** reproduce it (`quotient=[-1,0,0,0,1]`, seed 42, gives the circulant block `[-2 -2 -4 1 / 1 -2 -2 -4 / …]` where Sage gives `[2 3 -5 3 / 3 2 3 -5 / …]`). Sage draws through `PolynomialRing.random_element(degree=n−1)` -> `IntegerModRing.random_element()`, a *different* consumption pattern (leading coefficient first, with a redraw loop on zero). Our `PythonRandom` stream matches Sage's `python_random().randrange(q)` exactly, but the polynomial-ring layer above it is not ported, and the expected answer is Sage-version dependent. Structural invariants (block shape, circulant rows, `\|det\| = q^n` / `q^(m−n)`, minrep range, primal/dual relation) all match |
| `IntegerLattice.gen_lattice` (`crypto/lattice.ts`) | Sage's global randstate | A local seeded xorshift, so entries differ from Sage's doctest matrices even for the same seed |
| `LWE`/`RingLWE` `repr` | Prints `None` for an unbounded `m` | Prints `None`. Every numeric field matches Sage's doctests, and the **sampler's own repr** matches too: `Discrete Gaussian sampler over the Integers with sigma = 1.915069 and c = 401.000000` |
| `SBox` LAT | Per-mask Walsh-Hadamard transform | Same (values unchanged; AES went 175 ms -> 6.5 ms) |
| Affected modules | `sage/crypto/{lattice,lwe,boolean_function,sbox}.py` | `packages/sagemath-ts/src/crypto/` |

### Rationale

1. **Two different generators upstream** — Sage's `modular`/`random` lattice block goes through GMP's
   MT19937 (`c_random() % q`) while the `ideal`/`cyclotomic` block goes through CPython's via the
   polynomial ring. Both generators are ported bit-exactly (see
   [Random State](#random-state-and-seeding)); what is missing for `ideal`/`cyclotomic` is the
   *draw order* imposed by `PolynomialRing.random_element(degree=n−1)`.
2. **The expected answer for `ideal`/`cyclotomic` is Sage-version dependent** — the vendored doctest
   and SageMath 10.x disagree on that input — so pinning either would be pinning a version rather
   than a behaviour.
3. **Verification by invariant** — those two branches are verified with seed-independent structural
   and algebraic oracles rather than by pinning random values.

### Trade-offs

- `gen_lattice(type='ideal'|'cyclotomic')` does not reproduce Sage's published matrix for a given
  seed; the other two types do.
- `crypto/lattice.ts`'s `IntegerLattice.gen_lattice` still uses a local seeded xorshift.

### Behavioral Impact

`gen_lattice`'s three `modular` / `random` / `dual` doctests reproduce exactly, verified both against
the vendored reference and against the live Sage install. Sage's doctest values also reproduce
exactly for `LindnerPeikert(20)`, `RingLindnerPeikert(16)`, `Regev(20)`, the 3-round MISTY
construction, the 8×8 LAT of `SBox(7,6,0,4,2,5,1,3)` with all three scalings, `min_degree = 2`, the
`'03'`/`'43'`/`'00ab'` hex round trips, the algebraic-immunity cases and the dual lattice.

---
## Discrete Gaussian Samplers

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `DiscreteGaussianDistributionIntegerSampler` algorithms | `uniform+table`, `uniform+online`, `uniform+logtable`, `sigma2+logtable` | **All four.** `dgs_bern.c` (`dgs_bern_uniform_t` with its 32-bit bit pool, `dgs_bern_mp_t`, `dgs_bern_exp_mp_t`) and `dgs_gauss_mp.c`'s `dgs_disc_gauss_sigma2p_t` are ported, including the `sigma <- k·sigma_2` adjustment with `sigma_2 = sqrt(1/(2 ln 2))` and `MPFR_RNDN` rounding of `k` |
| Acceptance test | `mpfr_urandomb` compared against the tabulated probability | Same |
| `repr()` | `Discrete Gaussian sampler over the Integers with sigma = %f and c = %f` (six decimals, reporting the **adjusted** sigma for `sigma2+logtable`) | Byte-identical, including `sigma2+logtable`'s adjusted `sigma = 3.397287` |
| `DiscreteGaussianDistributionLatticeSampler.__call__` | Returns a vector over the base ring of the basis matrix — integers for `ZZ^n`, rationals for a `QQ` basis | `sample()`/`samples()` keep the `bigint[]` signature and throw `ValueError('lattice basis is not integral; use sampleExact() for exact rational samples')`; `sampleExact()`/`samplesExact()` return exact `Rational[]` for any basis. Public `isIntegral`, `basisExact`, `cNumeric()` were added |
| `DiscreteGaussianDistributionPolynomialSampler` | Lives in `sage.crypto.lwe` with signature `(P, n, sigma)`; there is **no** such class in `sage.stats.distributions.discrete_gaussian_lattice` | `crypto/lwe.ts` carries the faithful `(P, n, sigma)` class; `discrete_gaussian_lattice.ts` *additionally* exports a convenience `(n, options)` form returning a coefficient array |
| Non-spherical Σ | Matrix sigma, Peikert's `r`, Cholesky, offline samples, `_call_non_spherical` | Implemented: covariance matrices (with the scaled-identity collapse and `sigma_basis` for `Σ = S Sᵀ`), `_maximal_r()` by power iteration on `Q Σ^-1`, `_precompute_data`'s Peikert branch (exact `B_inv` over QQ, `r = 0.9999·_maximal_r`, `B2 = chol(Σ − r²Q)ᵀ`), `add_offline_samples`, `_randomise`, `_call_non_spherical` |
| `_normalisation_factor_zz` theta series | PARI `Q.__pari__().qfrep(B, 0)` | **Delegated** to `parigp-ts`'s `qfrep0`; `qfrep` here is a thin adapter that builds PARI's column-major `t_MAT` and raises PARI's own `incorrect type in qfminim` for a non-integral form |
| `_normalisation_factor_zz` working precision | `RealField(prec)` (MPFR) | Honoured: the sum is accumulated in a `RealField(prec)` layer (`RealNumberMP`) built inside the module, so `prec = 100` returns the full 28-digit answer |
| `precision='dp'` constructor keyword | Routes through `dgs_gauss_dp.c` (drand48 / libc `random()`) | `'mp'` (the default) works; `'dp'` throws naming the unported `dgs_gauss_dp.c` (`discrete_gaussian_integer.ts:438-441`); any other value raises Sage's exact `ValueError("Parameter precision '…' not supported")`. `dp` results are documented by Sage itself as not reproducible, so there is no oracle |
| Affected modules | `sage/stats/distributions/`, `dgs_gauss_mp.c`, `dgs_bern.c` | `packages/sagemath-ts/src/stats/distributions/` |

### Rationale

1. **`dp` carries no oracle value** — Sage itself documents that "in the latter case results are not
   reproducible", so there is nothing to match against.
2. **No polymorphic "element of the base ring"** — TypeScript cannot express it, and every downstream
   consumer (LWE, crypto) wants `bigint`. Splitting the accessor keeps the common integral case
   statically typed while making the general case exact, instead of the previous silent
   `BigInt(Math.round(...))` corruption.
3. **Avoiding cross-module churn** — moving or removing the duplicate polynomial sampler would break
   `crypto/index.ts` and `stats/index.ts`, so it is registered instead, with an `@see Deviation:`
   docstring pointing callers at the Sage-faithful version.

### Trade-offs and remaining divergences

- **`RealField(prec)` is a semantics re-implementation, not a transcription.** MPFR is not vendored
  under `reference/`, so `RealNumberMP` implements sign/mantissa/exponent with round-to-nearest-ties-
  to-even applied to the *exact* result — which makes `+`, `−`, `*`, `/` correctly rounded by
  construction — plus `pi` (Machin), `log 2` (`2 atanh(1/3)`), `sqrt` (integer square root) and `exp`
  (argument reduction) at 96 guard bits, and `real_mpfr.pyx`'s printing (`:1897`) and `round()`
  (`:3034`) rules. The four basic operations and the printing are exact; a transcendental may differ
  from MPFR in its last bit. Verified against mpmath at 200/500 bits to the full printed width,
  against 4000 random exact-rational operations, and against V8's independent `toPrecision(15)` on
  3000 doubles.
- **Two upstream *evaluation* semantics are reproduced deliberately**, because they are what the
  doctest value encodes: pynac pulls the numeric factor out of `(sigma·sqrt(2 pi))**n`, so `sigma^n`
  is rounded at sigma's own 53 bits before `(2 pi)^(n/2)` is applied (this is why Sage's `prec = 100`
  doctest reads `…969634991553` and not the mathematically correct `…995783045323`); and for
  `sigma <= 1` the argument and `exp` are `RealNumber(53)` operations, so `prec > 53` buys no
  accuracy there. The **third** pynac artefact — rewriting `exp(-y)` as `cosh(y) − sinh(y)` — is
  deliberately *not* reproduced; see
  [Upstream Behaviour](#upstream-behaviour-deliberately-not-reproduced).
- **`RealNumberMP.str()` refuses decimal exponents beyond ±100 000** (printing them exactly would
  need a ~28 Mbit power of ten). Only printing is affected.
- **LLL inside the non-spherical branch is a local exact-rational LLL** (delta 0.99) rather than
  fpLLL. It is used only to pick the enumeration centre.
- **Non-positive-definite `Σ` raises `RuntimeError`**, matching Sage, including column padding.
- **`sigma` and `c` are methods** (`D.sigma()`, `D.c()`), matching the vendored source, rather than
  public properties. `cNumeric()` is kept as a float convenience view.
- **`BernExpMp` stores `p = 0` where dgs leaves memory uninitialised.** `dgs_bern.c:121-124` breaks
  out of the table build when `exp()` underflows but still sets `l = i+1`, so index `i` may be read
  while `p[i]`/`B[i]` were never initialised — undefined behaviour in C. We store an explicit zero,
  which deterministically rejects.
- **`sigma2+logtable` rejects `sigma < sigma_2`.** dgs computes `k = round(sigma/sigma_2)` with no
  lower guard, so `k = 0` makes `mpz_urandomm` divide by zero. We throw `ValueError` naming
  `sigma_2 = 0.849322`. An added error, not a changed result.
- The lattice constructor now defaults `sigma` to 1 and preserves the original scalar
  validation order. See [Gaussian Constructor and Center State](#gaussian-constructor-and-center-state)
  for deferred centers and the remaining input/backend limits.

### Behavioral Impact

**All four integer-sampler algorithms reproduce SageMath's *seeded sample streams* exactly**, not
merely its distributions: 14 pinned 16-sample streams across three parameter sets and the two
non-integer-centre paths, plus Sage's own `_flush_cache` doctests, plus 200 000-sample chi-squared
tests at the 0.999 level over 11 settings. Every vendored `_normalisation_factor_zz` doctest
reproduces (`15.7496101985309`, `3.16536453178580`, `6.82492448921763`, and `prec=100` giving
`1.5585454565440389696349915528e27`), as do `_maximal_r() = 0.584028653716433`, the `Σ` recovered
from a `sigma_basis`, the exact `RuntimeError` text, the two `f()` values, the three `__repr__`
forms and the five spherical normalisation values. Rounding of a non-integer centre follows dgs's
round-half-to-even (`dgs_gauss_mp.c:161-165`), so `sigma=3, c=1.5, tau=2` has support `[-4, 8]`
(`c_z = 2`), not the `floor(c)` window.

---

## ZK Sumcheck and Multilinear Extensions

These modules port `reference/sage_blueprints/`, not SageMath itself.

| Aspect | Blueprint | sagemath-ts |
|--------|-----------|-------------|
| Constant round polynomial | `sumcheck_round_prover` raises `ValueError('prover: Layer polynomial is not univariate')` when `len(res.variables()) != 1`, which includes a constant. Running the blueprint on the zero polynomial does not even reach that error — it crashes with `AttributeError: 'IntegerMod_int' object has no attribute 'variables'` | Only a round polynomial depending on a variable **other** than the free one is rejected (`'prover: Layer polynomial is not built from the correct variable'`); a constant round polynomial is returned as-is |
| Affected modules | `reference/sage_blueprints/{sumcheck,mle}.sage` | `packages/sagemath-ts/src/zk/` |

For `sumcheckVerify`'s required `numVars`, `binaryToInt`'s `bigint` return and
`sparseMultilinearExtension([i])`, see
[Return Shapes](#return-shapes-keyword-arguments-and-signature-adaptations) and
[Port-Only APIs](#port-only-apis-with-no-sagemath-counterpart). For `booleanHypercube`'s bound see
[Honest Failure](#honest-failure-instead-of-silent-approximation).

### Rationale

**The zero function is a legitimate sumcheck instance** — it has a zero round polynomial in every
round and its sumcheck is sound; the blueprint's strictness is an artefact of Sage returning a
base-ring element rather than a polynomial when all variables vanish.

### Trade-offs

- A (spurious) error signal for constant round polynomials is lost.

### Behavioral Impact

`sumcheckRun` on all-zero values succeeds here and crashes in the blueprint. The round polynomial is
built **symbolically** (exact, arbitrary degree) rather than as a hardcoded line through `p(0)`,
`p(1)`, and variables are addressed by the ring's real names instead of assumed `x0, x1, …` — so a
ring named `[a,b]` gives `-2*a + 11`, and `x0^2*x1 + x0 + 1` over GF(101) gives `x^2 + 2*x + 2`, both
matching the blueprint run under a real `sage` binary.

---

## GF(2) Matrix PNG Functions

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| PNG I/O library | libgd (C library) | No external library; data-only functions |
| `from_png(filename)` / `to_png(A, filename)` | File I/O | Retained but throw `NotImplementedError` with guidance to the data-based alternatives (`matrix_mod2.ts`) |
| Substitutes | — | `from_png_data(width, height, pixels)` and `to_png_data(A)` (`matrix_mod2.ts`), with the black = 1 / white = 0 convention SageMath's libgd path uses (`matrix_mod2_dense.pyx:2613-2700`) |
| Affected modules | `matrix/matrix_mod2_dense.pyx` | `matrix/matrix_mod2.ts` |

### Rationale

1. **No native image libraries** — JavaScript environments have no universal PNG library like libgd.
2. **Environment portability** — browser, Node.js and Deno handle PNG differently.
3. **Separation of concerns** — raw pixel data lets users choose their PNG encoder.

### Trade-offs

- No direct file I/O; users must handle PNG encoding/decoding themselves (e.g. `pngjs` on Node.js,
  the Canvas API in browsers).
- Empty-image errors now use Sage's `cannot write image with dimensions {c} x {r}`.
  Binary pixel conventions and empty errors are compared with the original libgd functions.
  Grayscale thresholding is a data-adapter extension; Sage's file API expects a 1-bit image.

### Behavioral Impact

`from_png()`/`to_png()` throw instead of performing file I/O; the data functions provide equivalent
functionality with the same pixel convention.

---

## PARI Integer Factorization (parigp-ts)

`ifactor.ts` is PARI's real factoring chain; every stage of `ifac_crack` is present.

| Aspect | SageMath (PARI/GP) | sagemath-ts (parigp-ts) |
|--------|--------------------|-------------------------|
| `Z_factor` chain | `ifac_crack` (`ifactor1.c:2786`): trial division -> pure powers -> SQUFOF -> Pollard-Brent rho -> ECM (non-insisting) -> MPQS -> ECM (insisting), driven by the `ifac_decomp` worklist | The **same chain, same order**, ported with file:line citations: `tridiv_bound` + gcd-with-primorial fast trial division (`Z_oddprimedivisors_fast`, `:3306`), SQUFOF (`:1474`, incl. `squfof_ambig`), Pollard-Brent (`:1184`/`:1361`, incl. Brent fast-forward, backtracking and multi-factor returns), Lenstra-Montgomery ECM (`ellfacteur`/`ECM_loop`, `:752`/`:1038`), MPQS (`mpqs.ts`) |
| MPQS relation store | Large-prime relations spooled to disk (`pari_unique_filename`) | In memory (a `Map`). Same relations, same combined full relations; only a memory ceiling on very large inputs |
| MPQS decline threshold | Above 107 decimal digits (`mpqs.h:400`) | The same. Only then does `Z_factor` throw |
| Pollard-Brent stage gate | `ifactor1.c:1361-1367`: one-limb inputs are declined; two-limb inputs inspect the **low** limb against 2^32, despite the source comment describing a 96-bit gate | **Restored 2026-09-12.** The former widened search is removed; inputs now follow the native gate and later stages. |
| Pollard-Brent round budget | Native integer counters and seed selection; retry exhaustion raises the native internal-error message | **Restored 2026-09-12.** Exact BigInt counters, exact seed/size addition, unchanged native budget formula, and matching PariError. |
| SQUFOF acceptance bound | 2^46 on 64-bit builds (`:1487`) | **Restored 2026-09-12.** The former 2^59 widening is removed. Larger inputs proceed to the remaining native factorization stages. |
| ECM | Up to 64 curves in parallel with Montgomery batched inversion, Montgomery's PRAC chain, and a stage 2 over a 48-entry helix of residue classes mod 210 with a 1024-entry baby-step table | **Restored 2026-09-12.** Native batching, PRAC rules, continuation, whole-modulus failure handling, prime schedule and persistent state now replace the serial/binary/additive implementation. Seeds advance with BigInt. The existing insisting-round limit remains. |
| `is_357_power`'s residue sieve | mod 211/209/61/203, then 117/31/43/71 through a 106-entry mask table | **Restored 2026-09-12.** Omitting the sieve changed the returned mask. The native table, word-path valuation filter, and distinct word/multiword search priorities are now retained. Exact Newton roots still replace verified float guesses. |
| Perfect powers / `isprimepower` | `ispower.c`: `Z_issquareall`, `is_357_power`, `is_kth_power`, `is_pth_power`, `Z_isanypower_101`, `Z_isanypower`, `isprimepower` | All ported, using **exact integer k-th roots** (Newton) instead of PARI's `sqrtnr`/`mpexp` float guesses. `isprimepower` never factors `n` |
| `pollardbrent_i` on 4 unlucky restarts | `pari_err_BUG` (`:1330`) | **Restored 2026-09-12.** Raises the native PariError and message; six direct comparisons exercise this reachable path. |
| File placement | `ispower.c` | `isprimepower`, `Z_isanypower`, `Z_isanypower_101`, `Z_issquareall`, `is_kth_power`, `is_pth_power`, `is_357_power` and `Z_iroot` live in `ifactor.ts`; they should move to `parigp-ts/src/ispower.ts`. Each carries a JSDoc naming its upstream file and line |
| Affected modules | `pari/src/basemath/ifactor1.c`, `ispower.c`, `mpqs.c` | `packages/parigp-ts/src/ifactor.ts`, `mpqs.ts` |

### Rationale

1. **No floating point** — PARI's root guesses are heuristics followed by an exact check, so
   replacing the guess with an exact integer root cannot change the verified answer.
2. **Native stage gates are preserved.** The former SQUFOF/Pollard-Brent widening
   changed direct stage results and routing. Those deviations were removed on
   2026-09-12; comparisons cover the restored boundaries.

### Trade-offs

- Memory ceiling on very large MPQS inputs, where PARI spools to disk.
- The former serial-ECM performance and output differences were removed on 2026-09-12.
  JavaScript allocation/BigInt costs remain, while native batching and continuation are retained.
- The `ispower.c` routines are in the wrong file.

### Behavioral Impact

`Z_factor` is **correct and complete** everywhere it answers: agrees with brute-force trial division
exhaustively on 1..5000 and on 4000 random `n < 1e9`; a 2000-random-semiprime sweep where the product
is restored, every factor is BPSW-prime and the factorization equals exactly `{p,q}`; the published
factorizations of `F6 = 2^64+1`, `M67 = 193707721 · 761838257287`, `M71`, `M101` and `F7 = 2^128+1`;
and a 40-digit input in 47 ms agreeing with `sage: factor(...)`. `isprimepower` does not factor its
argument: for 24-digit primes `p`, `q`, `isprimepower(p·q)` is `null` and `isprimepower(p·p)` is
`[p,2]`, both in ~0 ms. `Z_factor` and `factoru` take an optional `options?: FactorOptions`
(`{ ecmRounds, mpqsMaxPolys }`) second parameter; existing one-argument call sites are unaffected.

---

## PARI Elliptic Curve Algorithms (parigp-ts)

| Aspect | SageMath (PARI/GP) | sagemath-ts (parigp-ts) |
|--------|--------------------|-------------------------|
| `ellcard` dispatch | Naive trace enumeration for `expi(p) < 11`, `Fp_ellcard_CM`, `Fp_ellcard_Shanks` in the middle range, SEA for `expi(p) >= 56` (`FpE.c:1424-1437`) | Same naive branch below `p = 2048`, then `Fp_ellcard_CM`, then Shanks, then **base Schoof** from `expi(p) >= 96` (`group.ts:1318`, `:1357`) — a measured threshold, but the wrong dispatch target now that SEA exists. See [parigp-ts Elliptic Curves](#parigp-ts-elliptic-curves--sea-dispatch-and-isogeny-stubs) |
| SEA (Schoof-Elkies-Atkin) | `ellsea.c`, needs the `seadata` modular-polynomial package | **Ported in full** as `Fp_ellcard_SEA` (`elliptic/ellsea.ts`): Elkies, Atkin and the match-and-sort final step, plus `Fp_elljissupersingular` and the CM branch. `seadata` is replaced by `polmodular.ts`, which computes `Phi_L` on demand and caches it — which is how PARI *generates* `seadata` in the first place |
| `Fp_ellcard_CM` | Full CM table (`Fp_ellj_get_CM` + `ec_ap_cm`) | **All thirteen** class-number-one discriminants, ported line by line from `FpE.c:624-666` and `:1282-1421`, delegating to `qfb.ts`'s `cornacchia2`. Includes PARI's signed-int `(CM&3)==0 -> CM>>=2` semantics and the `case -28: ap_cm(-7, -114, …)` quirk |
| `Fp_ellcard_Schoof` `j = 0` / `j = 1728` shortcut | `ellsea.c:1990-1993` | Not taken — routing back into `ellcard` would be a recursion hazard and would remove those curves from the Schoof test oracle. `ellcard` applies the CM shortcut before ever reaching Schoof |
| `cornacchia2` failure inside the `ap_*` helpers | PARI writes `(void)cornacchia2(...)` and ignores the return value, leaving the out-parameter as `gen_0` — which would report a wrong trace | Return `null`, `Fp_ellcard_CM` returns `null`, and `ellcard` falls through to Shanks/Schoof. Unreachable in theory and never fired in ~20 000 tested CM curves, but a plausible wrong cardinality is the one outcome to avoid |
| `gen_ellgroup` `m` output | `bb_group.c:1035-1043` writes `*pm = g1` and then overwrites it with the final iteration's `lcm(s,t)` | Returns `m = g1`. When the primes of `N0` are not all settled in one iteration, the final `m` need not be a multiple of `d2`, and then `gen_ellgens` can never terminate: measured on `E/F_43: y^2 = x^3+7x+8` (group `[12,3]`), about 0.5 % of runs produce `m = 4` with `d2 = 3` — 4 hangs in 885 runs. `g1` provably satisfies `d2 \| g1 \| d1`, and PARI 2.15.4 never hangs on that curve over 4000 fresh `ellgenerators` calls, so `g1` reproduces the *shipping* PARI behaviour |
| `Fp_ellcard_Shanks` visibility | `static` in `FpE.c` | Exported, so the test suite can exercise the BSGS branch against an exhaustive point-count oracle |
| `random_FpE` | `FpE.c:369-385` returns `Fp_sqrt(rhs, p)`, the canonical smallest root | Same. `<P>` and `<-P>` are the same subgroup, so order, group-structure and pairing consumers are unaffected |
| `j` / `ellj` return type | `t_INT` or `t_FRAC` | `bigint` when `c4^3` is divisible by the discriminant, else an exact `Ratio {num, den}` (with an exported `isRatio` guard) — there is no rational type in this package |
| Advanced functions (`ellisogeny*`, `ellfrobenius`) | Fully implemented | Stubs throwing `PARI_NOT_IMPLEMENTED` — see [parigp-ts Elliptic Curves](#parigp-ts-elliptic-curves--sea-dispatch-and-isogeny-stubs) |
| Barrel exports | — | `Ratio`, `isRatio`, `Fp_ellcard_CM`, `Fp_ellj_get_CM`, `Fp_ellj_nodiv`, `ec_ap_cm`, `Fp_ellcard_Schoof` and `Fp_elldivpol` are module-level exports **not** re-exported from `packages/parigp-ts/src/index.ts`; `Fp_ellcard_SEA` and `Fp_elljissupersingular` are (`index.ts:445-446`) |
| Affected modules | `pari/src/basemath/ellsea.c`, `FpE.c`, `bb_group.c`, `ellisog.c` | `packages/parigp-ts/src/elliptic/{group,points,init,advanced}.ts`, `ellsea.ts` |

### Rationale

1. **SEA needed the modular polynomials, so they were computed rather than read.** `reference/pari`
   ships the `seadata` *reader* but not the *data* (`reference/pari/data` is empty). Porting
   `polmodular.c`/`polclass.c`/`volcano.c` gives PARI's own `seadata`-less fallback
   (`ellsea.c:118-123`).
2. **Termination over literal fidelity for `m`** — a hang is not a faithful reproduction of a value.
3. **Exactness** — returning `Ratio` is how the j-invariant stays exact.

### Trade-offs

- **Base Schoof is `O(log^5 p)` with schoolbook `FpX` arithmetic** where SEA is `O(log^4 p)`.
  Measured on this port, single random curve, Schoof vs Shanks: 56 bits 12.8 s / 0.10 s; 64 bits
  21.1 s / 0.39 s; 72 bits 82.6 s / 2.41 s; 80 bits 101.7 s / 4.85 s; 88 bits 189.6 s / 26.1 s;
  96 bits ~358 s / 296 s (9 GB rss). So `ellcard` keeps Shanks below `expi(p) = 96`, not PARI's 56.
  The value returned is unaffected.
- `m = g1` means the Weil pairing is computed at a possibly larger exponent, i.e. marginally slower.
- `Fp_elldivpol(l, a4, a6, p)` is a **new public function with no PARI counterpart over `F_p`** (PARI's
  SEA uses modular polynomials rather than `psi_l`); exported so the recursion is testable.

### Behavioral Impact

`Fp_ellcard_SEA` returns exact cardinalities at every size: verified against **PARI's own `ellsea`
regression vectors** (`reference/pari/src/test/in/ellsea` entries `v[9]`, `v[10]`, `v[11]`) at 65, 70
and 101 bits, all exact, and independently at `p = 2535301200456458802993406410683` where it returns
`2535301200456457821343807570392`, byte-identical to
`sage -c "EllipticCurve(GF(p),[3,5]).cardinality()"`. `ellcard_sea` was verified exhaustively against
brute force on all 121 104 curves over every prime `5 <= p <= 120` and against Shanks on 375 random
curves from 20 to 88 bits; the zero-divisor split path fired 5984 times over 16 308 curves with zero
errors. `Fp_ellcard_CM` covers all thirteen discriminants, verified against brute-force point
counting (15 392 curves), Shanks (3744 curves to 48-bit primes), a counting-independent `[#E]P = O`
oracle (936 curves at 64/80/96 bits), and the **published SECG group orders** of secp160k1,
secp192k1, secp224k1 and secp256k1.

---

## ntl-ts GF2X Representation

| Aspect | NTL | ntl-ts |
|--------|-----|--------|
| `GF2X` storage | A `WordVector` (`xrep`) of machine words; `normalize()` strips zero words, `SetMaxLength(n)` preallocates, and the object can be temporarily unnormalized | A single bigint whose bit `i` is the coefficient of `x^i`, so it is **always** normalized: `normalize()` is a no-op and `SetMaxLength(n)` only performs NTL's negative-length check. `SetLength(n)` still truncates coefficients `>= n` exactly as NTL does |
| `XGCD` | Half-GCD recursion above `NTL_GF2X_GCD_CROSSOVER` (`GF2X1.cpp:3625`) | Plain extended Euclid |
| `PowerMod` | Sliding-window exponentiation (`GF2X1.cpp:1743`) | Binary square-and-multiply |
| `random`, `factor`, `SquareFreeDecomp`, `DistinctDegFactor`, `EqualDegFactor`, `BerlekampFactor` | Implemented | Honest `NTL_NOT_IMPLEMENTED` stubs (`GF2X.ts:472`, `:480`, `:488`, `:497`, `:505`, `:822`) — see [ntl-ts GF2X Factoring Stubs](#ntl-ts-gf2x-factoring-stubs) |
| `BuildRandomIrred` | Implemented (`GF2XFactoring.cpp:504`) | **Absent entirely** — there is no such symbol in `ntl-ts`, not even a stub. `GF2X.random` at `:822` is what stands in for it |
| Affected modules | `ntl/src/GF2X*.cpp`, `GF2XFactoring.cpp` | `packages/ntl-ts/src/GF2X.ts`, `GF2.ts`, `GF2X_irred_tab.ts` |

### Rationale

1. **bigint already is an arbitrary-precision bit vector** with XOR and shifts; hand-rolling a word
   vector would add no fidelity and duplicate what the runtime does.
2. **Identical results** — over GF(2) the gcd is unique, and the Bezout pair with
   `deg(s) < deg(b) − deg(d)`, `deg(t) < deg(a) − deg(d)` that NTL returns is exactly the one plain
   extended Euclid produces; sliding-window and square-and-multiply compute the same power.

### Trade-offs

- `XGCD` is `O(n^2)` bit operations rather than `O(n log^2 n)` — irrelevant at the degrees Sage's
  default-modulus path uses (n up to a few thousand).
- Factoring and random generation are unavailable.

### Behavioral Impact

None on any implemented function. `IterIrredTest`, `BuildIrred` and `BuildSparseIrred` are
line-for-line ports of `GF2XFactoring.cpp` over a vendored copy of NTL's 2049-row `GF2X_irred_tab`,
verified against Sage's `polynomial_gf2x.pyx` doctests (`BuildIrred_list(2/3/4/33)`,
`BuildSparseIrred == BuildIrred` for `n ∈ [1,32]`, `BuildSparseIrred(33) = x^33 + x^10 + 1`) and
against exhaustive brute-force irreducibility over all 2046 monic polynomials of degree <= 10.
`sagemath-ts`'s `polynomial_gf2x.ts` delegates its `add`/`sub`/`neg`/`mul`/`sqr`/`leftShift`/
`rightShift`/`trunc`/`divRem`/`gcd`/`xgcd`/`powMod`/`derivative`/`reverse`/`is_irreducible`/
`buildIrred`/`buildSparseIrred` to this package.

---

## Newly Ported Upstream Modules — Residual Divergences

The residual differences inside otherwise-faithful transcriptions of MPQS, `polmodular`/`polclass`,
SEA, `buch1`, `galconj`, `qfrep`, the Shanks-distance `t_QFB`, Laurent and multivariate power
series, polynomial matrices and van Hoeij. These are *residual* differences inside a faithful
transcription, not substitutions for the algorithm.

### Shared rationale

1. **Exactness where upstream is inexact for implementation reasons, and inexactness where upstream
   is inexact by design.** Where the inexactness is an implementation artefact (`minim0_dolll`'s
   Cholesky, `fmpz_lll`'s Gram-Schmidt, `galoisborne`'s `||den·V^-1||`) the port is exact, which can
   only agree with upstream or be more correct. Where upstream is *deliberately* inexact and the
   value is observable (`fmpz_poly_CLD_bound`, MPQS's Knuth-Schroeppel score, byte-scaled logarithms
   and the target size of `A`, `matrix2.pyx`'s `norm(A, 2)`) the port reproduces the same arithmetic,
   including `Float32Array` where PARI uses a C `float`.
2. **Randomness is Las Vegas everywhere it appears here**, so a deterministic seeded stream
   (xorshift) replaces `pari_rand` without affecting any answer — which is why PARI's own golden
   outputs match despite completely different randomness. Only *which* generators or relations are
   found can differ.
3. **Unreachable upstream branches are transcribed but flagged, not claimed.**
4. **Refusing beats guessing.** Every gap throws `NotImplementedError` naming the upstream routine
   and its `file:line`.

### MPQS (`parigp-ts/src/mpqs.ts`, from `mpqs.c`)

| Aspect | PARI | Port |
|--------|------|------|
| GF(2) kernel | `F2Ms_ker` (`F2v.c:1063`): dense `F2m_ker_sp` for `nbrow <= 640`, randomized block Lanczos above | **Restored 2026-09-12.** The binary-matrix backend preserves singleton elimination, block Lanczos, native retries, basis order and RNG state. The former dense substitute differed on 420 of 766 seeded comparisons; native Lanczos can return a proper kernel subspace. |
| `Fl_sqrt` | Native Fl_sqrt_pre_i selects the smaller root, including after its generator search | **Restored 2026-09-12.** MPQS delegates to the audited native word square-root backend. The former duplicate returned the larger root in 446 of 708 direct comparisons. |
| Relation table order | Native GEN hashing, bucket chains and resizing determine kernel column order | **Restored 2026-09-12.** The port preserves the native 64-bit hash words and incremental table order. The former insertion-order Map changed 45 of 94 driver comparisons, including RNG state. |
| Sieve inner loops | `mpqs_sieve_p`/`_p1`/`_p2` are 4x/8x unrolled and interleave the two progressions | Two plain loops per factor-base entry. The multiset of byte additions is identical, so the sieve array is bit-for-bit the same |
| `mpqs_eval_sieve` bit array | `__v2di` (16 bytes) with SSE2, else `ulong` | The 8-byte scalar layout. The threshold is always `>= 128`, so both collect exactly the bytes `>= threshold`, in increasing order |
| `relaprimes` / `relp` buffers | Fixed `MAX_PE_PAIR = 60`; a candidate with more distinct factor-base divisors overruns them | Sized to the factor base / growable. Strictly a bounds fix |
| Relation exponent packing | `pi \| (ei << 20)` in a 64-bit long, so `\|ei\| < 2^43` | The same packing in 32 bits, so `\|ei\| < 2^11`. Factor-mode exponents are bounded by `log2(4·A·Q(x))` (tens); negative exponents occur only in the unreachable class-group mode |
| Factor base layout | 32-byte union on 64-byte boundaries | Parallel typed-array columns; every field keeps its type, **including** the C `float` `fbe_flogp` (`Float32Array`), whose rounding participates in `mpqs_locate_A_range` and `mpqs_si_choose_primes` |
| `MpqsOptions.maxPolys` / `.debug` | No such knobs (`MPQS_DEBUG` is a compile-time `-D`) | Added. Defaults (`0`, `false`) reproduce upstream exactly; `debug` turns upstream's own `mpqs_check_rel` (`:1069`) and post-Gauss check (`:1525`) into a permanent test oracle |
| `mpqs_class_init` / `mpqs_class_rels` (`:1775`, `:1815`) | Present | **Absent.** Their original caller is `buch2.c`, which is not ported. The shared relation combiner has native comparisons as of 2026-09-12; other class-group branches remain outside that coverage. |

*Behavioural impact:* none observed. `mpqs_increment` was compared against the upstream C function
compiled verbatim (69 999 values, 0 mismatches); the 99 parameter rows and 41 multipliers are diffed
against `mpqs.h` as a permanent test. 500 random semiprimes `>= 2^46` and 460 mixed composites split
with 0 failures.

### `polmodular` / `polclass` / `volcano` (`parigp-ts/src/polmodular.ts`)

| Aspect | PARI | Port |
|--------|------|------|
| `find_j_inv_with_given_trace` | Picks a torsion constraint `m > 1` from Sutherland's `torcosts.h` tables and draws curves from the `X_1(m)` models in `crvwtors.c` (2345 lines of model data) | The `m = 1` / `twist = 3` entry — uniformly random curves plus the two-sided filter `(p+1)P == tP`, then the faithful `test_curve_order`. `m = 1` is one of the choices the upstream tables can return, so the algorithm is the same and the returned `j` distribution is unchanged; only the constant factor grows |
| `SMOOTH_INTS` / `HURWITZ_RATIO` | 1200 hand-written entries each | Generated at module load from the GP recipes in upstream's own comments, then diffed entry-by-entry against the vendored literals (0 mismatches over all 2400) |
| Return type | `RgM_to_RgXX`, a bivariate `t_POL` | A `ZM` (column-major `M[j][i]` = coefficient of `X^i Y^j`), matching `matkermod.ts`'s convention. `vx`/`vy` are still accepted and still produce upstream's `e_PRIORITY` error |
| Machine words | 62-bit machine primes | The prime-search loops carry the candidate in `BigInt` and throw above `2^53` rather than silently losing precision. Not hit for any level tested up to `L = 71` |
| Weber / double-eta / Atkin class invariants above their internal level | Supported | **Throw**, naming `polmodular.c:500-870`, the ~1500 lines of double-eta tables at `:2457-3663`, and `polclass.c`'s orientation machinery. SEA is unaffected: `ellsea.c:118-123` only ever asks for `INV_J` or `INV_G2` |
| `polmodular0_powerup_ZM` | Reachable | Fully transcribed but **unreachable and therefore untested** |
| `quadclassnos(D)` for `\|D\| >= 500000` | Falls back to Buchquad (`buch1.c`) | Throws. Unreachable for `polmodular`, whose discriminant search is bounded by `max_max_D = 320000` |

*Behavioural impact:* verified against PARI's **own** regression oracle — the DJB-style hash in
`reference/pari/src/test/in/polmodular` reimplemented verbatim reproduces all nine golden
`modpoly_hashes` for `inv = 0` (`L = 2 … 23`), the `INV_G2` entries and all four
`check_eval_modpoly` cases with both derivatives. `polclass0` matches PARI's `polclass(D)` for 34
discriminants. **One genuine upstream-transcription bug was found and fixed here during
verification:** `common_nbr` (`volcano.c:407-427`) returns `rlen`, the count of *distinct* roots of
the degree-2 gcd, and every caller branches on it being 2; the port returned `[r0, r0]` for a double
root, so `polclass0` rejected every `j`-invariant it drew and **never terminated** for non-fundamental
discriminants such as `D = −288`. Eight non-fundamental discriminants are now pinned against PARI.

### SEA (`parigp-ts/src/elliptic/ellsea.ts`, from `ellsea.c`)

| Aspect | PARI | Port |
|--------|------|------|
| Modular equations | Caches the `seadata` table; recomputes `polmodular_ZXX` on every call when `seadata` is absent | Cached per level in a module-global `Map`. `Phi_L` over `Z` depends only on `L`, and computing it is 80–90 % of the running time (24 s for `L = 71`), so a 101-bit curve takes ~20 s cold and sub-second warm; a second 256-bit curve in the same process drops from 262 s to 13 s |
| `FpXn_inv` | Newton iteration | The `O(n^2)` coefficient recurrence. The truncated inverse is unique, so the two agree exactly; degrees in `find_kernel` are `O(ell) <= 60`. `FpXn_expint` **is** transcribed as upstream's Newton loop, because its `FpX_integXn` divides by integers with a gcd trick a naive recurrence would not reproduce |
| `grp->hash` | PARI's generic `hash_GEN` | An FNV-style hash of the x-coordinate. Any deterministic hash is correct: every match is re-verified against the actual x-coordinates (`ellsea.c:1918`) before a cardinality is recorded |
| `NULL` dereferences upstream believes cannot happen | Undefined behaviour in C | Return `false` / skip the match-and-sort attempt / `PariBugError`. None was ever hit |
| Debug traces | Global `DEBUGLEVEL` | `setSeaDebugLevel(n)`, silent by default |
| `find_isogenous_from_Atkin` (`:900`) / `find_isogenous_from_canonical` (`:964`) | Reached when the modular equation has type `'A'` or `'C'` | **Throw.** `get_modular_eqn` (`:107-123`) only sets those types from a `seadata` file; with no `seadata` it always sets `'J'`, which is the path PARI itself takes |
| `Fq_ellcard_SEA` with `T != NULL` (extension fields) | Supported | **Absent.** Every routine is transcribed in its `T = NULL` form; adding `T` means re-deriving all of them over `FpXQ` |

*Behavioural impact:* PARI's own regression file (`test/in/ellsea` + `test/32/ellsea`) reproduces byte
for byte — all 11 `ellap` values from 65 to 200 bits and all 14 `ellsea(E, smallfact)` values. NIST
P-256 and Curve25519 come out exactly right (262 s and 232 s cold, 13 s for a further 256-bit curve
once the modular equations are cached). 640 random curves agree with Shanks/Mestre and 420 with
exhaustive counting.

### `buch1` class and unit groups (`parigp-ts/src/buch.ts`)

| Aspect | PARI | Port |
|--------|------|------|
| `t_REAL` | `mp.c`'s kernel, with AGM/Newton for `logr_abs` | Shared native kernels with immutable sign/mantissa/exponent/precision records. Constructors, elementary arithmetic, comparison, conversions, square roots and transcendentals all delegate; logarithms use native binary splitting and progressive series/AGM precision. Accuracy is asserted by the algorithm itself, exactly as upstream: `get_R` (`buch1.c:996`) only accepts a regulator when `h·R·invhr` lies in `(0.8, 1.3)` |
| `ZM_pivots` | Modular rank profile (`Flm_pivots`) certified by exact linear algebra | One-step fraction-free (Bareiss) elimination with the "first unused row" rule — the same canonical row rank profile PARI certifies |
| Randomness | `pari_rand` | Seeded xorshift32 (`setBuchRandomSeed`). `no` and `cyc` are canonical (which is why PARI's golden values match); the **generators** can differ from PARI's for a given `setrand`, as PARI's own do |
| MPQS relation collection for `\|D\| >= 2^60` | `mpqs_class_init`/`mpqs_class_rels` | `use_mpqs` is permanently `false`, i.e. always upstream's own fallback `imag_relations` (`buch1.c:746`). Same relations, same class group; slower (25 digits 0.4–2.2 s, 34 digits 12 s) |
| `hnfspec`'s overflow guard | `HIGHBIT` (2^63, or 2^31 on 32-bit builds) | `2^52`, the exact-integer range of a JS number — the correct constant for this kernel, just as `2^31` is for a 32-bit PARI |
| `bnfinit` for degree 2 | Runs the general `Buchall_param` and returns a full `bnf` | Returns Buchquad's class group, regulator, torsion order and unit-norm sign. The mathematical content is identical and verified; the missing parts need the same `nf` layer as degree > 2 |
| `hnfspec_i`'s `co > 300 && co > 1.5·li` branch | Reachable | Transcribed but **unreachable from this unit**, hence untested |
| `ZM_snfall_i` for non-square-HNF input | Supported | Throws. `W` out of `hnfspec`/`hnfadd` is always a square HNF with nonzero determinant |

*Behavioural impact:* verified against PARI's own `test/32/quadclassunit` — the complete `test(10^15)`
and `test(-10^15)` tables (608 discriminants, 0 mismatches), the four `quadclassunit(±2^81+c)` values
and every bug-report case in that file. Independently: all 599 discriminants `−3 … −1200` against the
exhaustive reduced-form count, ten fundamental `D` against the exact Dirichlet formula, and 28 real
discriminants against a Pell oracle.

### `galconj` Galois groups (`parigp-ts/src/galconj.ts`)

| Aspect | PARI | Port |
|--------|------|------|
| `galoisborne` | Computes the complex roots with `QX_complex_roots` at `t_REAL` precision, forms the inverse Vandermonde numerically and takes its operator sup-norm | An **exact Hadamard/Cramer bound**: `borne = ceil(den·n^(n/2)·B^(n(n−1)/2) / floor(sqrt\|disc T\|))` with `B` Cauchy's root bound. The trace-based absolute bound is at least one, preserving native `ceil_safe` on real zero. Any *upper* bound is correct here — it only sizes the `l`-adic accuracy and a rejection threshold |
| `indexpartial` | Strict partial factorization and `ZpX_reduced_resultant_fast` | Source-based strict-factor, Sylvester/echelon and dynamic p-adic refinement; the former exponent-only approximation was removed on 2026-09-12. See PARI integral-basis denominator adapters. |
| "Combinatorics too hard" | `frobeniusliftall` and `testpermutation` print a warner and then give up or `return identity_perm(n)` — PARI can return a **wrong/partial group with a warning** | Both throw. There is no warning channel here, and the thresholds (10^15 and 10^14 tests) are never approached |
| Output certification | Trusts its p-adic bounds | An extra `certify` pass: every element must be a distinct permutation of the `l`-adic roots induced by its own polynomial, and the set must be closed under composition, else `PariBugError`. It never fired |
| Galois irreducible-factor helper | Native full `Flx_factor` machinery | `_galconj_factor_squarefree_irreducibles` delegates to the full native factor pipeline after its squarefree-input guard. `FpX_factor_squarefree` remains the separate multiplicity-indexed component API; full factorization preserves native random-state consumption. |
| Root ordering / choice of `l` | Whatever the splitting algorithm produces; cyclotomic `T` short-circuits through `galoiscyclo` | Prime-field roots now retain native word/small-degree ordering; no cyclotomic short-circuit. The permutation **labels** and `gal.roots`/`gal.p` can differ from PARI's by a relabelling (e.g. `l = 17` vs PARI's 41 for `polcyclo(8)`); everything label-independent is identical. `galoisinit(x^4-x-1)` returns `null` where live PARI returns `0` |
| `s4galoisgen` (`:1519`) / `f36galoisgen` (`:1698`) | Present | **Throw**, naming `FpX_ffisom`/`FpXQ_ffisom_inv`/`FpXV_ffisom`/`FpXV_chinese`/`FqC_FqV_mul`, none of which exists in `parigp-ts`. Falling through to the generic search would *hang* for S4 |
| `findpsi` (`:411`) | Called when `P` is not squarefree mod the current prime | Fully transcribed but **not reached by any of the 21 verified fields**, hence untested |
| `galoisgenlift_nilp` and the polycyclic layer (`:2389-2744`) | Used when `!(ga->group & ga_easy)` | Unreachable: `galoisanalysis` sets `ga_easy` for every degree `<= 104`. The guard throws by name, so it can never silently take a wrong path |
| `galoiscyclo`, `galoisinitfromaut`, `galoissplittinginit` | Present | Absent. `galoissplittinginit` needs `nfsplitting0` (`base1.c:1413`), outside `galconj.c` |

*Behavioural impact:* the exact bound is *looser* than PARI's, so the `l`-adic accuracy is higher
(`valabs` 37 vs PARI's 16 for `x^6+108`, 266 for A4, 1756 for a degree-24 field — which is why degree
~24 is slow). No effect on the answer: for 21 fields the port reproduces PARI's relative orders, its
`nfgaloisconj(T, 4)` polynomials **character for character**, its `#galoissubgroups` and its complete
`galoissubfields(G, 1)` lists.

### `qfrep` theta series (`parigp-ts/src/qfrep.ts`)

| Aspect | PARI | Port |
|--------|------|------|
| Fincke-Pohst enumeration | C `double` Cholesky data and running norms, `BOUND = borne·(1 + 1e-10)`, norms recovered by rounding | Identical enumeration order and pruning predicates, every quantity a `bigint` — see [Exact Arithmetic](#exact-arithmetic-where-sagemath-uses-floating-point) |
| Return shape | 1-indexed `t_VECSMALL` | 0-indexed `bigint[]`. cypari2's flag bit 1 is accepted and ignored, because both are the same JavaScript array |
| `pari_err_PREC` from Cholesky precision loss | Possible | Cannot occur (the arithmetic is exact). The two PREC checks that depend only on the size of `B` **are** reproduced verbatim, so the accepted range of `B` is exactly PARI's |

*Behavioural impact:* none on any oracle — PARI's own `test/in/qf` 12-dimensional form, cypari2's
doctests, exhaustive enumeration for identity forms in dimensions 1–5, A2, D4 and 3475 random
positive definite forms, an independent enumeration of E8 in its `D8+` coordinate model, and the
classical theta series of E8, D4, A2, `r_2` and `r_4`. Cost: bigint arithmetic is ~5–10x slower than
PARI's doubles in the inner loop (~6M enumerated vectors/second).

**One upstream defect is deliberately not reproduced.** PARI's GP documentation prints
`qfrep([2,1;1,3], 5, 1) = Vecsmall([1, 0, 0, 1, 0])`; the last entry counts vectors of norm 10 and
`q(-1,2) = 10`, so it must be 1. Confirmed three ways: by brute force, by reading `bibli1.c:1327`,
and by running the same call through a live PARI 2.15.4, which returns `Vecsmall([1,0,0,1,1])` — so
PARI's *code* agrees with the port and only its GP doc example is stale.

### Shanks-distance forms (`parigp-ts/src/qfb.ts`)

| Aspect | PARI | Port |
|--------|------|------|
| `t_REAL` rounding | Native per-kernel word precision, truncation and guard rounding | Addition/subtraction, multiplication/squaring and division now preserve native word algorithms. General exponent overflow remains under audit; addition uses the documented deterministic padding choice for an upstream out-of-range read |
| `qfrpow` for `n <= -2` | Inverts the form and then passes the **signed** `n` to `qfr5_pow`/`qfr3_pow`, which invert again — so PARI returns `x^\|n\|`. Verified on a live PARI 2.15.4: `qfbpow(f,-6) == qfbpow(f,6)` | Inverts once and raises to `\|n\|`, so `x^-n` is the inverse of `x^n` |
| `qfr5_pow` exponent loop | Loops over the machine **words** of `n` with `if (m == 1 && i == 2) break;` and an arithmetic `m >>= 1`; a word with leading zero bits skips squarings it owed, and a word with its top bit set never terminates. `qfbpow(f, 2^64+1)` returns `qfbpow(f,3)`; `qfbpow(f, 2^63)` overflows the PARI stack | A plain right-to-left binary chain over the whole `bigint` exponent — identical for every single-word `n`. Our `f^(10^20)` distance is congruent mod the regulator to its form's cycle distance (residual `< 1e-20` in exact rationals); PARI 2.15.4's answer has residual 1057.8 with `R = 2641.55` |
| `qfrpowraw` distance sign | Forms the distance *after* negating `n`, so `x^-k` reports `+k·d`, contradicting its own form | Uses the original signed exponent. Ten of the 360 oracle values differ from PARI by exactly this sign flip |
| `qfr5_to_qfr`'s `mplog2(lg(d0))` | Passes a **word length** where a bit precision is expected — a call site missed when PARI 2.16 changed `prec` from words to bits (`Qfb.c:488` vs `:428`) | `mplog2(precision(d0))`. The branch fires whenever `fix_expo` has (reachable: `qfbpow([f,0.], 10^8)` fires it three times); with upstream's `lg(d0)` the distance would be wrong by ~1e5 |
| `qfr_1_fill`'s `subiu(y,1)` | Reads the container `y` as an integer — a typo | `y2 − 1`, which is what `qfr_1_by_disc` computes for the same discriminant. Unreachable from every public entry point, hence untested |

*Behavioural impact:* PARI's own vendored `test/32/qfb` distances reproduce at GP's default 38 digits
(up to the last-digit rounding above); 339 of 360 values from a live PARI 2.15.4 match form-for-form
with the distance agreeing to `< 1e-35`, the 21 exceptions being the four divergences above. The
principal cycle accumulates PARI's `quadregulator(D)` on eight discriminants, and 1800+ checks
confirm that every distance the port produces is the cycle distance of its own form modulo the
regulator.

### Polynomial matrices (`sagemath-ts/src/matrix/matrix_polynomial_dense.ts`)

| Aspect | SageMath | Port |
|--------|----------|------|
| Class vs functions | `Matrix_polynomial_dense` is a Cython subclass; every operation is a method | Exported free functions taking the matrix first, exactly as `matrix_operations.ts` already does for `matrix2.pyx`. Names, defaults, error messages and outputs are unchanged |
| `degree_matrix` | `matrix(ZZ, …)` | `number[][]`, matching `row_degrees`/`column_degrees` |
| Immutability | Calls `set_immutable()` on the returned matrices | Ordinary mutable matrices — `Matrix` here has no immutability flag |
| `reverse`'s negative-degree check | `Polynomial.reverse(d)` raises for every polynomial including zero; the matrix method just forwards | Validated at the matrix level, because this repo's `Polynomial.reverse` returns early on the zero polynomial without checking. The matrix-level behaviour and message are identical to Sage's; `Polynomial.reverse(-1)` on the zero polynomial is still wrong and is flagged as a bug in a file that unit did not own |
| `hermite_form` | A method | Aliased to `polynomial_matrix_hermite_form` in `matrix/index.ts`, because `matrix_decompositions.ts` already exports a `hermite_form` for constant matrices. Same for `degree`, `truncate`, `shift`, `reverse` |
| `_hermite_form_euclidean` | Lives in `matrix2.pyx` | Implemented privately inside `matrix_polynomial_dense.ts`; architecturally it belongs in the `matrix2` port |
| `inverse_series_trunc`, `solve_left/right_series_trunc`, `left/right_quo_rem`, `reduce`, `minimal_interpolant_basis`, `minimal_kernel_basis`, `minimal_relation_basis`, `basis_completion` | Present | Not ported and not stubbed. (`minimal_approximant_basis` and `is_minimal_approximant_basis` **are** ported, at `:1721` and `:1508`, and re-exported from `matrix/index.ts:230-231`) |

*Behavioural impact:* none. Every doctest value in the ported functions' docstrings passes verbatim
(including the QQ, `GF(2^3)` and `GF(2^4)` examples and the issue #41278 regression),
`is_weak_popov`/`is_popov` are additionally brute-forced over all 256 2×2 matrices over `GF(2)[x]`
with degree `<= 1` entries, and the algebraic identities (`U·A == form`, `det(U)` a nonzero constant,
Popov idempotence, `hermite_form == popov_form(shifts)`) hold on random matrices over
`GF(2,3,5,7,11)`.

### `Frobenius_filter` over a number field (`sagemath-ts/src/schemes/elliptic_curves/isogeny_class.ts`)

| Aspect | SageMath | Port |
|--------|----------|------|
| Primes used in the walk | `for P in K.primes_above(p)`, every prime of good reduction, using `a_P` from `E.reduction(P)` | Only the primes of **residue degree one**. Residue degree `f > 1` needs `#E(F_{p^f})`, i.e. `Fq_ellcard_SEA` over an extension field, which `parigp-ts` does not have. Dropping primes can only make the filter **weaker** (a superset), never unsound |
| `division_polynomial(2).is_irreducible()` over `K` | PARI `nffactor` (`polnf.c`) | A two-sided **certificate**: a cubic is reducible over `K` iff it has a root in `K`, so irreducibility is certified by a degree-one prime at which the cubic has no root in `F_p` (one exists by Chebotarev), and reducibility by the root itself — Hensel-lifted from a completely split prime, rationally reconstructed, then **verified by exact integer arithmetic in `Z[theta]`**. Never guesses |
| `global_integral_model` | Scales until the a-invariants lie in `O_K` | Scales until they have integral coordinates in the power basis of `theta` (strictly stronger, since `Z[theta] ⊆ O_K`). Can pick a slightly larger `u`, making a few more primes look bad; those are skipped |
| Non-monic or non-integral defining polynomial | Goes through `K.pari_nf()`, which rescales to an algebraic-integer generator | Throws (`isogeny_class.ts:1370-1390`), naming that step |

*Behavioural impact:* none measured. All three SageMath doctests reproduce exactly — the `d = −23`
degree-6 case prints the verbose transcript line for line ending
`List of primes after filtering: [2, 3]`, the `Q(i)` case gives `[2, 3]` and the issue-36780 case
`[3, 5]`. 189 curves over six number fields were cross-checked against live SageMath on **both** the
filter output and the `include_2` boolean: 0 mismatches, 0 throws.

### Free modules over `K[x]` (`sagemath-ts/src/modules/free_module.ts`)

The intersection of two submodules of `K[x]^n` previously returned a basis off by a unit of `K[x]`.
The cause was not the echelon routine — upstream normalises **nothing** on this path, because
`_echelon_form_PID`'s reduction above the pivots sits inside a `try/except AttributeError` on
`Ideal.small_residue`, which only `NumberField` ideals implement. The defect was one step upstream:
SageMath's `intersection` routes the stacked basis matrix through `Matrix.integer_kernel`, which
first multiplies by `self.denominator()`, and over `QQ[x]` that is `Polynomial.denominator()` — the
lcm of the denominators of the rational **coefficients**, a non-trivial unit of `QQ[x]`. The port
computed only fraction-function denominators (always 1 for polynomial entries), so it never scaled.
Over `GF(p)[x]` the coefficients have no denominator, upstream falls back to 1, and the port was
already correct.

*Behavioural impact:* 914 randomly generated cases across five sweeps (`QQ[x]` and
`GF(2,3,5,7,11,13)[x]`, dimensions up to 5) agree with SageMath coefficient-for-coefficient on all
four matrices (`P`, `Q`, `P ∩ Q`, `P + Q`). Golden values were produced on SageMath 10.3 with the
vendored 10.9.beta4 `_echelonized_basis` patched in — see
[Vendored vs Installed](#vendored-sagemath-109beta4-vs-installed-103).

---

## Exact Return Types and Numeric Backends

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `IntegerLattice.volume()` when `rank < degree` | `gram_matrix().determinant().sqrt()`, an exact symbolic `sqrt(N)` | Returns `bigint` when the Gram determinant is a perfect square, otherwise a `SqrtInteger` (`modules/free_module_integer.ts`) carrying the radicand exactly. It prints as `sqrt(14)` and coerces to a double through `valueOf` |
| `ell_torsion.order_from_multiple` rejection | Vendored 10.9 (`groups/generic.py:1418`) raises `ValueError("The order of P(=…) does not divide …")`; installed 10.3 (`groups/generic.py:1266`) uses a bare `assert` | Raises `AssertionError('')`, matching the executable oracle. `groups/generic.ts` independently follows the vendored 10.9 behavior |
| `Rational.gamma()` at nonpositive integers | Unsigned infinity | Existing `UnsignedInfinityElement` singleton; the return type is `Rational \| UnsignedInfinityElement`. Nonintegral gamma values still raise `NotImplementedError` because the symbolic return backend is absent |
| `ComplexNumber.gamma()` return type | `ComplexNumber`, or `UnsignedInfinityRing.gen()` at a pole | Same behavior via the port-only `UnsignedInfinityElement` singleton, so callers handle a TypeScript union |
| `PowerSeriesElement.div`/`inv`/`pow(-n)` return type | An element of the Laurent-series fraction field whenever the valuation goes negative | Same behavior, declared as `PowerSeriesElement<T> \| LaurentSeriesElement<T>` |
| `ReedMullerCode.length()` / `minimum_distance()` / `parameters()` / `decoding_radius()` | Sage integers | `bigint`; JS `number` cannot exactly represent all valid code parameters |
| Bessel and error functions | MPFR, correctly rounded at the field's precision | Double-double arithmetic (`rings/real_mpfr_dd.ts`), rounded once to a double. It agrees with every pinned MPFR oracle value but is not provably correctly rounded |
| `ComplexNumber.zeta()` / `gamma()` | PARI | Borwein's Algorithm 2 for zeta and the `g = 607/128` Lanczos coefficients for gamma |

### Rationale

The return-type adaptations preserve exact values in TypeScript. The numeric backends replace
single-precision approximations but cannot provide MPFR/PARI's arbitrary-precision guarantees.
Following installed SageMath 10.3 for the one assertion keeps the executable oracle authoritative.

### Trade-offs

Several return types are unions, `SqrtInteger` and `UnsignedInfinityElement` are port-only value types,
and the special functions remain bounded by IEEE-754 output precision.

### Behavioral Impact

Values agree with the executable oracle on the covered domain. TypeScript callers must handle the
documented unions and port-only exact wrappers.

---
## Rational Comparative Adapters

The `rationals` area compares every public Rational member with Sage coercion, arithmetic,
or an explicit convenience adapter. `from(number)` uses Sage's simplest rational convention;
string construction uses GMP base-zero integer/quotient syntax. Decimal and exponential
strings are rejected as upstream. `toNumber()` transcribes `mpq_get_d_nearest` and compares
exact IEEE-754 encodings, including signed zero and subnormal/overflow rounding.

**Rationale:** the port exposes operators as methods and numerator/denominator/sign as properties.
`roundToRational`, `ndigits`, `nbits` and `factorial` are convenience methods; their oracles use
exact scaling with Sage rounding, component digit/bit counts, and integer factorial respectively.

**Trade-offs:** dispatch coverage does not cover unimplemented symbolic results. Square-root
comparisons use `extend=False`; nonintegral gamma is outside the implemented domain. Logarithmic
heights remain binary64 approximations and the oracle compares eight decimal places. The
precision arguments cannot provide Sage's arbitrary precision.

**Behavioral impact:** exact rational results and binary64 conversion agree on the tested domain.
Approximate transcendental results do not promise MPFR rounding at arbitrary precision. Infinity
uses the types documented above. The bundled Sage default rounding mode is explicitly `even`;
installed Sage 10.3 needs that argument to reproduce the bundled source.

---

## Scalar Constructor Adaptations

**SageMath:** Integer and QQ constructors accept Python integers, booleans, floats, None,
strings and ring elements. Integer strings support PEP 3127 prefixes, underscores and GMP
whitespace rules. Digit lists require a base greater than one. QQ tuple components are
coerced through Integer, and its constructor's None handling differs from a nested list.

**Port:** Integer and ZZ share a BigInt constructor path. `undefined` and `null` represent
None; booleans map to 0/1. `_integer_` hooks and `lift()` results supply IntegerLike values
through structural interfaces. Arrays represent digit lists for Integer (with a base), and
pairs/singleton lists for QQ. Rational.from accepts IntegerLike, boolean and null values.

**Rationale:** JavaScript has no Sage Element hierarchy or separate tuple type. Sharing
coercion prevents the ring and element constructors from drifting and preserves Sage's
numeric values, validation order and known-type error messages.

**Trade-offs:** only conversion protocols backed by the port's existing types are available;
Python-specific objects such as numpy/gmpy2 values are outside this API. String diagnostics
follow Python repr, and strings preserve GMP's C-string semantics, including an embedded NUL.

**Behavioral impact:** covered numeric, string, array and element conversions match Sage.
String bases use a signed C int; digit-list bases remain arbitrary precision. The surprising
`ZZ(QQ(3), base=0)` NotImplementedError is deliberately reproduced: `Q_to_Z` inherits
`Map._call_with_args` (map.pyx:854-857). `Integer(QQ(3), base=0)` succeeds through Integer.__init__.
Likewise `QQ(None)` is zero while `QQ([None])` raises TypeError. These are preserved upstream
behaviors, not new missing backends.

---

## Integer Root Backend

**SageMath:** `Integer.nth_root` delegates to GMP `mpz_root`; Rational roots delegate to
Integer roots. GMP uses precision-doubling Newton extraction in `mpn/generic/rootrem.c` and
a padded extra-limb optimization when only root exactness is needed.

**Port:** the same precision schedule, prefix/remainder update, derivative quotient, clipping,
downward corrections and padded-limb optimization use native BigInt operations. Both Rational
root methods now share Integer's backend. The reference is the official
[GMP 6.3.0 release](https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz), since GMP is not bundled here.

**Rationale:** GMP arithmetic maps to native BigInt in this project. An exact one-bit seed
replaces the nine-bit logarithm lookup seed; this avoids reproducing approximate lookup tables
while retaining precision doubling and the same correction bounds.

**Trade-offs:** the initial small-precision stages and BigInt allocation/multiplication costs
differ from GMP's machine-word implementation. No binary search over the full-size root
remains. The 64-bit padding optimization avoids the final full-size power for ordinary
non-perfect inputs, as upstream.

**Behavioral impact:** roots, truncation flags, exceptions and Cython's signed 32-bit exponent
bounds match Sage on the comparative domain. The minimum negative Rational exponent preserves
the source's signed-wraparound behavior, including its `ValueError`. Large powers, neighboring
values, negative inputs and maximal valid exponents have comparative regressions.

---

## Integer Quadratic Class Number Backend

Sage's `Integer.class_number()` calls PARI `qfbclassno` (`quad.c:610-618`), selecting
`classno`/`classno2`. The port delegates to its existing PARI `quadclassno`/`Buchquad`
(`buch1.c:1282-1288`) for either discriminant sign.

**Rationale:** this removes the incomplete 110-entry table using an already ported PARI algorithm.
**Trade-offs:** the internal PARI algorithm and proof machinery differ; no `proof` option is exposed.
**Behavioral impact:** ordinary class numbers agree on tested positive, negative and nonfundamental
discriminants, including -10000003 and 40000012. The table's input restriction is removed.

## Integer Audit Oracle Boundaries

The shared `integers` cases document their adapters in `python/areas/integers.py`:
port-only methods map to Sage free functions/operators; nth-root comparisons verify existence
and the returned root rather than requiring the intentionally different root choice. Approximate
real logarithms and heights compare to eight decimal places. General `core(t)` uses a Sage
factorization certificate; strong probable-prime checks use GMP's independent implementation.

The installed primecountpy backend crashes on this host, so prime-count cases use PARI `primepi`.
Sage `Integer.powermod` calls GMP without checking inversion (`integer.pyx:3656`) and can abort the
process for a negative exponent on a non-unit. The port keeps its safe `arith.power_mod` delegation;
those error cases compare with Sage's safe free function, explicitly in the oracle.

**Rationale:** avoid treating unavailable or crashing oracle results as successful comparisons.
**Trade-offs:** these cases establish mathematical/API equivalence on the stated domains, not
identical backend execution or identical transcendental approximations.
**Behavioral impact:** the port raises `ZeroDivisionError` instead of aborting on the unsafe
negative modular-power inputs. The other adapters normalize documented TypeScript representations.

Sage `sigma(n,k)` itself raises `AttributeError` for negative `k` and `|n| > 1`, because it calls
`divide_knowing_divisible_by` on a Rational (`arith/misc.py:1679-1680`). This port now reproduces
that class/message, while preserving the successful empty-product cases at ±1. There is no
output deviation here; this note prevents a future “repair” from changing upstream behavior.

---
## PARI Gram reduction and HNF transformation adapters

**Original:** `lllgramint` returns the image transformation from
`ZM_lll(G, 0.99, LLL_IM | LLL_GRAM)`. Singular matrices may return fewer columns;
indefinite matrices are not uniformly rejected. `ZM_hnflll` accepts an integer
remove flag: 1 removes zero HNF columns but keeps all transformation columns;
2 also removes the corresponding kernel columns from the transformation.

**Port:** `lllgramint` returns zero-indexed bigint columns, never a null sentinel.
Its pre-existing optional `{n,d}` rational delta is converted to native double
precision and passed to ZM_lll. `ZM_hnflll` uses dummy slot zero in both dimensions
and a boolean remove option, exposing native modes 0 and 1. It returns `{H,B}`;
B is null only when the caller does not request it. The original inputs are kept.

**Rationale and trade-offs:** Retain the established typed and immutable APIs while
using the native reduction algorithm. The delta option is an extension and native
remove mode 2 is not exposed by this boolean signature. Callers that want only the
image columns must select the trailing columns of B themselves.

**Behavioral impact:** Supported native transforms, kernel retention and Smith-group
representatives are compared exactly, including signs, diagonal residues, machine-word
reduction thresholds, rectangular matrices, and singular/rank-zero forms. Rational
delta controls include 1,200-bit numerator/denominator scaling. Native PARI versions
can select different valid relations; the comparative oracle uses bundled 2.18.1.

**Resource boundary:** 27 indefinite quaternion-form probes exhaust the bundled
256 MB PARI stack during Gram LLL; JavaScript instead raises RangeError when its
BigInt size limit is reached. The native stack allocation/error text is not emulated.
These failures have no mathematical result and are retained separately in
`lll_resources.native.json`, with permanent tests for the actual host errors.
They are not counted as exact transcript matches. Direct valid-input guarantees
and the earlier rank-deficient FLATTER limit remain documented under
[PARI adaptive LLL and FLATTER adapters](#pari-adaptive-lll-and-flatter-adapters).

---

## PARI rational trace and norm adapters

**Sage/PARI:** Quadratic elements use direct rational formulas. General absolute
trace/norm use the polynomial-modulus branches of gtrace/gnorm: constants multiply
by the degree or raise to that degree; polynomial traces use derivative/product/
remainder, and norms use a rational resultant divided by the leading coefficient.
They do not compute an entire characteristic polynomial for these operations.

**Port:** No-argument NumberFieldElement trace/norm follow those branches. Quadratic
formulas divide by the defining equation's leading coefficient. PARI RgX, ZX,
polarit3 and alglin2 expose the QQ polynomial dependencies; modular integer
resultants and packed integer products use the existing native ports.

**Rationale / trade-offs:** RationalPolynomialData is `[integerCoefficients,denominator]`
and RationalPair is `[numerator,denominator]`, replacing GEN rational tags. Nonzero
input denominators may be negative; outputs are canonical with positive denominators.
The gnorm/gtrace entry points here accept only a reduced RationalPolMod record,
`{value,modulus}`. They are not universal GEN dispatchers. QQ-only RgX kernels do
not cover arbitrary complex, finite-field, p-adic or nested polynomial coefficient
rings. ZX_rem/QX_ZX_rem/QXQ_mul require monic integral divisors; QXQ_norm requires a
nonzero integral modulus. Direct ZX_rem rejects nonmonic inputs instead of running
an invalid native buffer contract. Native constant-polynomial 0/1 results at a
constant RgXQ modulus are represented by the same numeric RationalPair as scalar
0/1; their GEN polynomial tag is not retained.

**Behavioral impact:** Exact trace/norm values, fractions and the compared generic
zero-divisor errors match. The port's absolute_trace/relative_trace aliases still
refer to absolute trace, as documented. The characteristic-polynomial matrix
algorithm, optional base-field arguments and complete relative extensions remain
separate fidelity work. Large-field tests disable irreducibility checking on both
sides to isolate trace/norm. Constructor degrees 30 through 300 now follow
PARI's integer factorization, and the other degrees follow NTL; see constructor
integer factorization routing.

## PARI prime-decomposition matrix adapters

**PARI vs port:** FpM image, supplementary basis and inverse retain native binary,
word and generic field dispatch, row-pivot preferences, CUP recursion, and original
independent columns. Public FpM matrices have unused slot zero on both axes; internal
row matrices omit it. F2m routines take packed BigInt columns and an explicit row
count. ZM_hnfmodprime takes zero-indexed columns. Null represents a failed inverse
or native null pivot pointer, separately from an empty matrix.

**Rationale:** Typed arrays replace GEN headers and packed C word buffers.
**Trade-offs:** malformed C buffers and composite moduli are outside the adapters'
valid-input contracts; native pointer identity is not represented. **Behavioral
impact:** compared valid matrices retain exact entries, pivot choices and native
empty-suppl errors, including unreduced inputs and rectangular matrices.

## Number-field full prime decomposition adapters

**Sage/PARI:** Prime decomposition, ideal factorization, primality and principal
ideal HNF construction use PARI's nf object. Sage caches a Factorization object and
sorts prime ideals by their HNF matrices in the chosen integral basis.

**Port:** base2 implements the Kummer, partial-Kummer and Buchmann--Lenstra paths,
including native prime generators, ramification/residue degrees, anti-uniformizers,
residue-degree limits and PARI ordering. base4 implements rational-HNF idealfactor,
idealnumden, idealismaximal, idealhnf_principal and idealadd. base3's zkmultable_inv
and zkmultable_capZ delegate to the modular ZM_gauss solver. Sage-layer ideal HNF
construction follows principal HNFs then idealadd, with native modular bounds.
NumberField.factor accepts NumberFieldIdealInput; ideal.factor returns cached
`Array<[NumberFieldIdeal,bigint]>`, sorted by Sage's row-wise HNF comparison. e/f
accessors read the cached prime data. NfIdealData has optional `index`; the factory
supplies it, with an exact basis-determinant fallback for older direct callers.

**Rationale / trade-offs:** GEN nf and factorization objects use typed records and
arrays, so Factorization methods and immutable/container semantics are not supplied.
Rational HNFs use `[integerColumns, positiveDenominator]`; idealnumden uses null for
a unit denominator and an empty numerator for zero. The inverse multiplication
vector similarly uses a common denominator. Low-level routines require valid
maximal-order data, prime moduli and certified HNFs; zkmultable_inv rejects singular
tables. These contracts do not emulate invalid C buffers. Native approximate norm
estimates use the existing exact resultant fallback. The complete nfinit reduction
and embedding pipeline remains unported; comparisons of internal prime data and
factor ordering therefore use a common nf_NOLLL basis. The original Sage
Factorization sorter and ideal comparison are executed on those same matrices.

**Behavioral impact:** Compared factors, negative exponents, cache identity,
recomposition, prime metadata, errors and common-basis order match. Public
mathematical ideal comparisons use canonical power-basis lattices. Generator
pretty-printing, strings, BNF methods and basis-dependent presentation outside that
common basis retain their existing gaps. Large-degree ideal construction no longer
uses the Sage-local unbounded HNF elimination that exhausted BigInt memory.

## Number-field ideal valuation adapters

**Sage/PARI:** `I.valuation(P)` coerces a non-ideal P through the field factory,
validates nonzero/prime/field in that order, then delegates rational ideal HNF
valuation to PARI. Zero I returns infinity. The original foreign-prime error has
an invalid two-placeholder Python format expression and raises TypeError.

**Port:** Version 24.0.0 corrects the previous reversed `P.valuation(element)` API.
`valuation(P: NumberFieldIdealInput): bigint | 'Infinity'` uses the established
infinity sentinel, preserves native error precedence (including that format error),
and delegates to `parigp-ts/base4.idealval`. Galois ramification callers use the
corrected signature. The convenience `prime_below()` method, absent from Sage's
absolute ideal class, returns the same rational prime as `pari_prime()[0]` and
`smallest_integer()`. Its old trial-division search is removed.

Native `idealHNF_val` uses the anti-uniformizer, primitive fractional content,
intersection/norm bounds, and machine-limb modular reduction schedule. `base3.ZC_nfval`
uses the native periodic content removal. `base2.primedec_end` constructs native
five-component prime metadata from certified prime HNFs, using prefix/suffix
intersections, the original uniformizer search and finite-field dependencies.
NfIdealData now also stores the defining monic integral polynomial, integral basis,
basis denominator and discriminant, all supplied by `nfmaxord_ideal_data`.

**Rationale / trade-offs:** The matrix-route `primedec_end` adapter accepts certified
prime HNFs; the complete Kummer/partial-Kummer/splitting pipeline now lives in PARI's
base2 port. Exact resultant norms replace native complex-embedding estimates, at a
potential performance cost while retaining certified integer norms. Tests compare
native generators, e/f and tau in the same integral basis (`nf_NOLLL`) and public
valuations in canonical power-basis coordinates. Complete nfinit basis reduction,
embeddings and generator pretty-printing remain open.

**Behavioral impact:** Integral/fractional/zero valuation results and supported
coercions match the comparative cases. Multi-generator ideal text inside valuation
errors is replaced by its canonical lattice in the comparison only; the remaining
message and error class are compared exactly. Raw display discrepancies remain
recorded. Strings are outside NumberFieldIdealInput and remain the previously
documented conversion gap. A zero vector is outside the low-level ZC_nfval domain:
original C loops indefinitely, while the port raises RangeError; the public zero
ideal path correctly returns infinity. That resource boundary is a separate control.
Low-level idealval also validates positive denominators and square HNF dimensions;
these checks protect the typed buffer adapter rather than emulate invalid C buffers.

## Number-field ideal class method adapters

**Sage:** Direct NumberFieldIdeal objects, including nonzero ones, lack the
fractional subclass's numerator, denominator, divides, is_coprime, factor and
residue/ramification methods. Unary inversion is unsupported on the base class.
The field factory and prime decomposition produce NumberFieldFractionalIdeal
objects for nonzero ideals. Prime-only methods preserve their distinct errors:
ramification/residue degree use pari_prime's ideal-specific message, while
residue_field raises `The ideal must be prime`.

**Port:** The existing common TypeScript base interface keeps these method names
callable on factory return values, but invoking a fractional-only method on a
direct base object raises the corresponding native AttributeError or TypeError.
The guard tests the class rather than whether its generators are zero. Prime
construction goes through the same nonzero-ideal factory as user construction.

**Rationale and trade-offs:** Factory return types include both the base zero ideal
and nonzero fractional ideals. Retaining the common static interface avoids an
unrelated TypeScript narrowing requirement. Runtime property discovery differs:
JavaScript can still see a guarded method on a base object, while Python's
hasattr reports it absent. Normal method calls and class-specific failures are
compared against Sage. Fractional ideal factorization and valuation delegate to
their independently audited PARI routines.

---

# Part II — Open Fidelity Gaps


## September 2026 Core Audit — Resolved Findings

F1–F9 from [AUDIT-2026-09.md](AUDIT-2026-09.md) are repaired. Permanent shared-case
comparisons live in `core_fidelity`; the expanded `integers` area covers all 107 implemented
Integer methods, with explicitly documented adapters for port-only conveniences. The audit
report retains the discovery evidence and records the additional repairs. These findings are
no longer open deviations. Remaining backend/API gaps below are separate work items.

**These are work items, not accepted deviations.** Each entry describes a place where the port is
less faithful than the vendored upstream allows, names the upstream file that implements the missing
behaviour, and estimates the effort to close it. None of them has a rationale that survives
"the upstream source is vendored and we have not transcribed it yet".

---

## Number Field Class Groups, Units and Galois Closure

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Class group / class number for degree > 2 when `h > 1` | PARI `bnfinit` (`Buchall_param`, `buch2.c:3946`) | `NotImplementedError` naming `Buchall_param`. Two rigorous sub-cases *are* implemented: if no prime ideal has norm `<= M_K` the class group is provably trivial, and the certificate also accepts a **proof of principality** for each factor-base prime (exhibit `alpha in P` with `\|N(alpha)\| = N(P)`, which forces `(alpha) = P`). That closes `Q(2^(1/3))`, `Q(sqrt2, sqrt3)`, `Q(zeta_7)`, `Q(zeta_8)`, `x^5−x−1`, `x^6+243` and Dedekind's `x^3−x^2−2x−8`, and reproduces Sage's `[1,1,1]` Hecke-polynomial doctest. `NumberField([-19,0,0,1]).class_group()` and `.class_number()` throw; `Q(zeta_23)` throws |
| `regulator()` and `UnitGroup.fundamental_units()` for degree > 2 | Same `bnfinit` | Throws. `NumberField(x^3−2).regulator()` throws where Sage returns `1.34737734832938`. Everything the regulator is *built from* exists — the certified archimedean embeddings that turn a unit into a row of `log\|sigma_i(u)\|`, and the proved torsion subgroup — so only the `r1+r2−1` free generators are missing |
| Galois group of a non-Galois field | Returns the group of the Galois closure (`galois_group.py:268` -> `number_field.py:9199` -> `splitting_field.py:371`), whose main loop is `nffactor` + `rnfequation`. For `x^3−2` Sage prints `Galois group 3T2 (S3) with order 6 of x^3 - 2` | `NumberField([-2,0,0,1]).galois_group()` **succeeds** and returns an object whose `toString` is `Galois group of degree 3`; the `NotImplementedError` is raised lazily from `galois_group.ts:251` only when `order()`/`elements()` is touched. A caller who only prints the object gets a silently misleading answer |
| `fixed_field`'s `polredbest` post-processing | `fixed_field(name, polred, threshold)` defaults `polred=True` for degree `<= 8` (`galois_group.py:889-899`), applying `polredbest(flag=1)` | Returns PARI's raw `galoisfixedfield` answer, i.e. exactly what Sage returns for `polred=False`. Isomorphic field, uglier model: for `K = Q[x]/(x^4+1)` and the automorphism `(1,2)(3,4)` we print `x^2 + 4` where Sage's default prints `x^2 + 1`. The other two fixed fields (`x^2−2`, `x^2+2`) are already in Sage's form |

### Why this is still open

The `nf` layer these all sit on is the single largest unported piece in the repo. It is not a design
choice; nothing about TypeScript prevents it.

### How to close

- **Class group / regulator:** port `reference/pari/src/basemath/buch2.c` `Buchall_param` (`:3946`)
  and its dependencies: `nfinit_basic`/`nfinit_complete` (`base1.c:2104`, `:2143`), `nfmaxord`
  round 4 (`base2.c:462`), `idealprimedec` (`base2.c:2386`), HNF ideal arithmetic (`base4.c`), the
  T2 form + LLL (`lll.c`), `small_norm`/`rnd_rel`/`getfu`/`makeunits` (`buch2.c:2540`, `:2860`,
  `:1126`, `:1238`). `parigp-ts/src/buch.ts:3277-3299` documents the exact missing routines, and the
  relation-matrix half (`hnfspec_i`/`hnfadd_i`/`ZM_snf_group`) is already ported there.
  **Effort: large** — several thousand lines, a multi-pass project.
- **Galois closure:** port `nffactor` (`reference/pari/src/basemath/nffactor.c`, ~2000 lines) and
  `rnfequation`/`polcompositum` (`base5.c`), plus
  `reference/sage/src/sage/rings/number_field/splitting_field.py:371`. **Effort: large.** But the
  cheap half — making the failure honest by throwing at *construction* time rather than lazily, so
  the misleading repr is impossible — is a one-line change and should be done immediately, per the
  port's own [Honest Failure](#honest-failure-instead-of-silent-approximation) policy.
- **`polredbest`:** port `reference/pari/src/basemath/base1.c:2672` `polredbest(T, flag)` (and
  `rnfpolredbest` at `:3083`). **Effort: medium** — it needs the order's T2 quadratic form and an
  LLL over it; the exact integral LLL already exists (`pari_nf.ts` imports it from
  `matrix/matrix_integer.ts`), so the missing pieces are the T2 form and the `polredabs` candidate
  loop. No `polredbest` exists anywhere in `packages/`.

### Trade-offs of leaving it open

Class numbers, class groups, regulators and fundamental units are unavailable for most number fields
of degree > 2, and `galois_group()` on a non-Galois field is actively misleading until the
construction-time throw lands.

### Behavioral Impact

`NotImplementedError` where SageMath answers — except the Galois-group repr, which is wrong rather
than absent.

---

## Number-field relative trace, norm and characteristic-polynomial routing

**Sage vs port:** Sage norm accepts a base field or embedding; generic trace accepts
an optional base (the quadratic trace override rejects extra arguments). The port
preserves QQ/null, own-field and degree-one-field values and return parents, plus
argument-count and primitive-input errors. Impossible field-degree divisibility
raises Sage's no-embedding error. Proper intermediate fields still need embedding
and relative-matrix construction; that supported-type gap raises NotImplementedError.
Morphisms are not yet accepted by the public TypeScript signature. Characteristic/
minimal polynomials still use the local matrix characteristic-polynomial routine,
without Sage's variable/algorithm options or PARI dispatch below the tuning threshold.

**Rationale for the current gap:** relative matrix/embedding conversion and
characteristic-polynomial dependency routes have not yet been ported. **Trade-offs:**
proper-subfield results, morphisms, polynomial parent presentation and native
characteristic-polynomial complexity are not preserved. **Behavioral impact:** the
covered base arguments now preserve exact values, return parents and errors, while
proper intermediate fields and large characteristic polynomials remain open.
The QQ/degree-one trace paths reuse the absolute trace kernel; Sage constructs a
matrix, whose trace has the same value. The matrix/cache interface is not yet exposed.

**How to close:** port number_field_element.pyx matrix/_matrix_over_base and
charpoly/minpoly dispatch, together with PARI polarit2/polarit3 characteristic
polynomials. RelativeFromAbsolute and relativize require exact basis conversion
and rational matrix dependencies. The constructor's PARI degree window is now
delegated, as are the NTL degree ranges.
**Effort:** medium to large for the remaining relative-field and matrix interfaces.

## Number-Field Kernel Not Delegated to parigp-ts

| Aspect | CLAUDE.md requires | sagemath-ts |
|--------|--------------------|-------------|
| Where the PARI `nf` routines live | Delegation to `parigp-ts` | `packages/sagemath-ts/src/rings/number_field/pari_nf.ts` still holds integral-basis, conjugate and unit kernels. Active prime decomposition, ideal construction, valuation, factorization and primality now delegate to PARI base1–base4; the old local primedec helper remains a separately exported legacy API |
| `nfgaloisconj` | One implementation | **Two.** `number_field.ts:31` still imports `nfgaloisconj` from `./pari_nf.js` and `number_field.ts:1939` (`automorphisms()`, the general case) still calls it, while `parigp-ts/src/galconj.ts` holds the real `galoisgen`/`galoisinit` and `galois_group.ts:15-23` already imports `galoisinit`/`galoisfixedfield`/`galoispermtopol` from `@sagemath-ts/parigp-ts`. Two independent implementations of the same PARI routine are live in one repo, only one of which is the delegation target CLAUDE.md mandates |

### Why this is still open

The fixes landed in parallel work units that did not own `parigp-ts`. That is a scheduling
constraint, not a rationale.

### How to close

- Move `pari_nf.ts` into `parigp-ts` as an `nf` module and re-point `number_field.ts` at it.
  **Effort: medium** — a file move plus import rewiring; the algorithms are already written and
  verified.
- Re-point `number_field.ts:1939` at `galconj.ts`'s `galoisconj4` and delete `pari_nf.ts:1229-1570`.
  **Effort: small** — the replacement is already written and verified character-for-character
  against PARI on 21 fields; this is a re-point plus deleting ~340 lines and its test-file section.

### Trade-offs of leaving it open

The dependency graph does not mirror SageMath's, and a divergence between the two `nfgaloisconj`
implementations would be invisible.

### Behavioral Impact

None on values today — `pari_nf.ts`'s `nfgaloisconj` reaches the same proved answer as PARI, by `n`
independent lattice reductions rather than from group generators (slower on large degrees:
`Q(zeta_25)`, degree 20, takes 3.1 s).

---

## Quadratic Class Numbers Not Delegated to Buchquad

Two remaining call sites compute quadratic class numbers by enumeration, when the real
subexponential algorithm is already ported, exported and oracle-tested in this repo.

| Site | SageMath | sagemath-ts |
|------|----------|-------------|
| `class_group.ts` | PARI `bnfinit`/`quadclassunit` | Enumerates reduced primitive forms of the field discriminant (definite for `D < 0`, rho-cycles of indefinite forms for `D > 0`), composes with Dirichlet composition (Cohen 5.4.7) and reads the elementary divisors off the `\|G[p^k]\|` counts. Guarded by `CLASS_GROUP_DISC_BOUND = 2_000_000n` (`:501`), above which `:704` throws `SAGE_NOT_IMPLEMENTED: class group of discriminant D requires PARI bnfinit/quadclassunit`. `QuadraticField(10000003).class_number()` is 2 in Sage and throws here. `class_group.ts:489` still carries a stale comment claiming `parigp-ts` has neither routine |
| `cm.ts` | `D.class_number(proof)`, i.e. PARI `qfbclassno` | Counts reduced primitive positive-definite forms (`-a < b <= a <= c`, `b >= 0` when `a == c`, `gcd(a,b,c) = 1`), memoized in a module-level `h_dict` (`:359`, `:376`). `OrderClassNumber` is used for non-fundamental discriminants exactly as in Sage. `O(\|D\|)` rather than PARI's `O(\|D\|^(1/4))` |

### Why this is still open

`parigp-ts/src/index.ts:529-531` **already exports** `Buchquad`, `quadclassunit0` and `quadclassno`
— the real McCurley-Buchmann index calculus, verified against all 608 discriminants of PARI's own
`test/32/quadclassunit` regression output and against `quadclassunit(±2^81+c)`. Executed:
`quadclassno(-10000003n)` returns `706n`, matching `QuadraticField(-10000003).class_number()`, and
`Buchquad(40000012n)` returns `no = 2`, `cyc = [2]` instantly. `class_group.ts` imports nothing from
`parigp-ts` (`:13-15`). This is a **wiring gap, not a capability gap**; the delegation was
deliberately not made in the same pass as the port.

### How to close

Import `quadclassno`/`quadclassunit0` in `class_group.ts` and `cm.ts`; delete the duplicated
enumeration and bound, and rerun the CM doctests. `Integer.class_number()` now delegates to
`quadclassno`, with the old table deleted and the ordinary/narrow convention corrected.
**Effort: about a day** — two or three call sites to rewire and ~300 lines to delete. The
`Integer.class_number()` portion is complete with 33 comparative boundary/discriminant cases.

### Trade-offs of leaving it open

The remaining class-group/CM call sites still have the bounds above;
binary-quadratic-form code is duplicated across `class_group.ts`, `cm.ts` and `binary_qf.ts` (which
already delegates to `parigp-ts`).

### Behavioral Impact

Values are correct where they answer — the reachable input range is the only difference.

---

## PARI/NTL Routines Duplicated or Ported In Place

CLAUDE.md requires that where SageMath delegates to an external library, we delegate to our port of
that library. These are the remaining violations.

| Routine | SageMath delegates to | We implement it in | How to close |
|---------|----------------------|--------------------|--------------|
| `isprimepower` / `Z_isanypower` / `Z_isanypower_101` | PARI `ispower.c`, via the integer prime-power method | Repaired in 8.11.0: both arithmetic predicates delegate to `parigp-ts/ifactor.isprimepower`; the duplicate tiny-prime/any-power/root helpers are removed | Existing exact-root and native word-path optimization differences of the dependency remain documented above |
| `matfrobenius` (`alglin2.c:428-720`) | PARI | `packages/sagemath-ts/src/matrix/matrix_integer.ts:1440` `frobenius_form_integer`, with a docstring citing `alglin2.c:617`/`:688` | Add a matrix module to `parigp-ts` and move it. **Effort: medium** — the code is written and verified; it needs a home |
| `dilog`, `incgam` | PARI | `packages/sagemath-ts/src/rings/complex_mpfr.ts` | Add a transcendental-functions module to `parigp-ts`. **Effort: medium** |
| `qfgaussred` | `self.__pari__().qfgaussred()` | `quadratic_forms/quadratic_form__local_field_invariants.ts`, a line-for-line port of `alglin2.c:1650-1749`; `parigp-ts` exports only `qfgaussred_positive` | Move it into `parigp-ts` beside `qfgaussred_positive`. **Effort: small** — the output is already pinned against real PARI in a test, so the move is safe |
| Binary polynomial factorization and standalone GF2X helpers | The generic Sage `Polynomial.factor()` uses PARI; the bundled GF2X template has no NTL factor override | `Polynomial.factor()` now delegates to the native PARI binary pipeline. The separate `polynomial_gf2x.ts` convenience helpers still use local algorithms | The standalone helper contracts and NTL stubs remain separate audit work; see [ntl-ts GF2X Factoring Stubs](#ntl-ts-gf2x-factoring-stubs) |
| `parigp-ts` error classes | — | `PariTypeError`, `PariDomainError`, `PariDimError`, `PariInvError`, `PariPrimeError`, `PariSqrtnError` and `PariFlagError` are defined once in `matkermod.ts:78-120` and imported by `ffinit.ts` and `qfb.ts`; `errors.ts` now defines the generic Gen-boundary `PariError`, while these specific classes remain in `matkermod.ts` | Consolidate the existing specific classes in `errors.ts`, preserving their exports. **Effort: trivial** |

### Why this is still open

File ownership: each fix landed in a work unit that did not own the dependency package.

### Trade-offs of leaving it open

The duplicate prime-power implementation is removed. Other dependency paths still differ from SageMath, so the
"delegate to our port" rule cannot be checked mechanically.

### Behavioral Impact

None on outputs — each is a transcription of the cited upstream source, checked against PARI/Sage
values.

---

## parigp-ts Elliptic Curves — SEA Dispatch and Isogeny Stubs

| Aspect | PARI | sagemath-ts (parigp-ts) |
|--------|------|-------------------------|
| `ellcard` above the crossover | `FpE.c:1431` calls SEA | `group.ts:1318` defines `SCHOOF_BIT_THRESHOLD = 96` and `:1357` calls `Fp_ellcard_Schoof(a4, a6, p)`. `Fp_ellcard_SEA` exists, is correct, and is exported from the package root (`index.ts:445`) — it returned the exact cardinality at 101 bits in 20.9 s cold (sub-second warm), where the measurement table for base Schoof puts 96 bits at ~300 s with 9 GB rss |
| `ellisogeny`, `ellisogenyapply`, `ellisogenycompose`, `ellfrobenius` | Implemented in `ellisog.c` | `throw new Error('PARI_NOT_IMPLEMENTED: …')` at `advanced.ts:1442`, `:1476`, `:1497`, `:1535` |

### How to close

- **Dispatch:** change `group.ts:1357` to call `Fp_ellcard_SEA`, as `FpE.c:1431` does, then
  re-measure the CM/Shanks/SEA crossover and update the threshold constant.
  **Effort: hours** — a one-line dispatch change plus a re-measured threshold and a regression sweep.
  No new code needs writing.
- **Isogenies:** transcribe `reference/pari/src/basemath/ellisog.c` (1756 lines — PARI's whole isogeny
  layer) into `parigp-ts/src/elliptic/`. **Effort: a few days** — self-contained (Vélu + isogeny
  composition over `FpXQ`), comparable in size to the completed `ellsea.ts` port.

### Trade-offs of leaving it open

`ellcard` is correct but orders of magnitude slower than it needs to be above 96 bits; four PARI
entry points throw.

### Behavioral Impact

None on values. Callers who need SEA can invoke `Fp_ellcard_SEA` directly in the meantime.

---

## ntl-ts GF2X Factoring Stubs

| Routine | NTL | ntl-ts |
|---------|-----|--------|
| `factor`, `SquareFreeDecomp`, `DistinctDegFactor`, `EqualDegFactor`, `BerlekampFactor` | `GF2XFactoring.cpp` | `NTL_NOT_IMPLEMENTED` stubs at `GF2X.ts:472`, `:480`, `:488`, `:497`, `:505` |
| `random` | Implemented | Stub at `:822` |
| `BuildRandomIrred` | `GF2XFactoring.cpp:504` | **Absent** — no symbol at all |
| Consequence | — | `sagemath-ts/src/rings/polynomial/polynomial_gf2x.ts` keeps four local factoring routines, and `irreducible_element(n, algorithm='random')` for `p = 2` uses rejection sampling from `current_randstate()` instead of `BuildRandomIrred`. (That fallback is exactly what SageMath itself takes when its NTL import fails, `polynomial_ring.py:3615-3620`, but the distribution differs from NTL's and the concrete polynomial for a given seed differs from Sage's) |

### How to close

Transcribe `reference/ntl/src/GF2XFactoring.cpp` (966 lines, vendored): `SquareFreeDecomp:66`,
`DDF:209`, `EDFSplit:279`, `EDF:312`, `SFCanZass:358`, `CanZass:410`, `BuildRandomIrred:504`. Then
delete `polynomial_gf2x.ts:552-720`. One port retires two documented deviations.
**Effort: several days** — it needs `GF2XModulus` and a `RandomStream` first (NTL's ChaCha-based
generator). A deterministic seeded stream would suffice, since EDF is Las Vegas — the same argument
`buch.ts` already makes.

### Trade-offs of leaving it open

Factoring over GF(2) does not go through the NTL port, so `sagemath-ts` and `ntl-ts` can disagree,
and `algorithm='random'` does not reproduce NTL's distribution.

### Behavioral Impact

The local routines produce correct factorizations; only the provenance and the random distribution
differ.

---

## Elliptic Curves over Q and Number Fields

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Torsion subgroup in characteristic 0 | `ell_torsion.py:173-176` is literally `G = self.__E.pari_curve().elltors()` for `K = QQ`; over a number field it uses division polynomials and height methods | `ell_torsion.ts:74-79` throws `NotImplementedError('Torsion subgroup computation over number fields requires PARI/GP elltors')`; `torsion_bound` (`:753-756`) throws for number fields. Finite fields work |
| Anomalous ECDLP | `ell_point.py:4640-4642`: when the base field is prime and `n == p`, `log()` delegates to `padic_elliptic_logarithm` (Smart/SSSA, `O(log p)`) | Falls through to the generic Pohlig-Hellman/BSGS path, `O(sqrt p)`. `padic_elliptic_logarithm` (`ell_point.ts:1478`) is a stub whose finite-field branch throws at `:1507` (and carries a dead `const order = P.order();` at `:1502`), so routing there would turn a working call into a failure |
| Rational isogenies of degree 43, 67 or 163 | Sage uses its precomputed exceptional kernel-polynomial table after the Kenku j-invariant test | The Kenku dispatch recognises the same exceptional j-invariants, but `isogenies_prime_degree` raises an explicit `SAGE_NOT_IMPLEMENTED` error naming the missing precomputed kernel table. Ordinary rational classes, including the complete three-curve `11a1` class and its degree matrix, are traversed exactly |
| `EllipticCurveIsogeny.formal()` | `hom.py` `EllipticCurveHom.formal` evaluates the rational maps for every isogeny | Vélu isogenies now return the genuine series and agree with the `ec_advanced` oracle. Other construction algorithms raise `NotImplementedError('formal expansion is only implemented for isogenies built by Velu')` rather than returning the old placeholder `t` |

### How to close

- **Torsion:** port `ellQtors` + `torsbound` into `parigp-ts`
  (`reference/pari/src/basemath/elltors.c:166-260` and `:610-624`; `ellnftors` at `:585-608`), then
  delegate. `torsbound` only needs reduction mod small primes and `parigp-ts` already has
  `ellcard`/`ellap` over `F_p`; `t2points`/`tpoint` need rational roots of division polynomials,
  which the port now has via van Hoeij. **Effort: large, ~400–600 lines**, dominated by introducing a
  `t_ELL_Q` curve type + `ellintegralmodel` in `parigp-ts` (`packages/parigp-ts/src/elliptic/` is
  entirely `F_p` today). The Q-only case (Mazur's 15 groups) is the cheap 80 %.
- **Non-Vélu `formal()`:** evaluate the generic rational maps over the Laurent-series layer rather
  than using the Vélu-specific formula data. **Effort: medium** — the formal group is complete, but
  the current y-map is exposed as an evaluator rather than coefficient data.
- **Exceptional isogeny degrees:** transcribe Sage's precomputed kernel-polynomial data for
  degrees 43, 67 and 163 and route the already-implemented Kenku j-invariant dispatch through it.
  **Effort: small-moderate** — the traversal and model normalization are already complete; the
  missing input is the exceptional table, not an algorithm.
- **Anomalous ECDLP:** transcribe
  `reference/sage/src/sage/schemes/elliptic_curves/ell_point.py:4650-4736` (a ~35-line Smart/SSSA
  attack whose only machinery is `EllipticCurve(Qp(p, 2), [ZZ(t) + k*p for t in E.a_invariants()])`,
  `lift_x(all=True)`, scalar multiplication and the formal-group ratio `-(x/y)`) and restore the
  `ell_point.py:4640-4642` dispatch. **Effort: medium** — the algorithm is trivial, but it needs
  `EllipticCurve` over `Qp` with `lift_x`; check whether `ell_generic` already accepts a `pAdicRing`
  base ring, otherwise budget a couple of days for the Qp curve layer.
### Trade-offs of leaving it open

Torsion over Q is unavailable. Rational isogeny classes now compute normally; only the three
exceptional table-backed degrees fail honestly. Non-Vélu formal expansions fail honestly.
Anomalous-curve discrete logs are `O(sqrt p)` instead of `O(log p)` — a real performance cliff on
exactly the curves an attacker would target, which for a crypto-focused port is the wrong thing to
leave documented.

### Behavioral Impact

One performance cliff and explicit `NotImplementedError` on the remaining unsupported paths.

---

## p-adic Precision Models, Extension Fields and L-Series

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Precision models | Capped-relative, capped-absolute, fixed-mod, floating-point, lattice, relaxed | **Capped-relative only.** `padic_generic.ts:546` `Zp(p, prec)` and `:554` `Qp(p, prec)` expose no `type` parameter |
| Extension fields (Zq, ramified/unramified) | Full element arithmetic | `pAdicExtension` (`padic_generic.ts:507-540`) is a shell with `base_ring`/`degree`/`absolute_degree` and **no element type**, so `nth_root`'s p-th-root extraction is written only for absolute degree 1 (`e = f = 1`), where it reproduces SageMath's result exactly |
| p-adic L-series | Full modular symbols, PARI, p-adic fields | `padic_lseries.ts`: `modular_symbol`, `measure`, `series`, `order_of_vanishing`, `frobenius`, `Dp_valued_series`/`_height`/`_regulator` and `_c_bound` all throw `NotImplementedError`. `_prec_bounds` propagates `_c_bound`'s error rather than over-reporting precision (it previously returned a literal `0`, which is not conservative: `c` is subtracted from the e-bounds) |
| Affected modules | `sage/rings/padics/`, `sage/schemes/elliptic_curves/padic_lseries.py` | `packages/sagemath-ts/src/rings/padics/`, `schemes/elliptic_curves/padic_lseries.ts` |

### Why this is still open

Sage's p-adics are Cython/`mpz`, **not** PARI, so the gap is not a missing binding — it is the
unported precision-model classes in `sage/rings/padics/`. `teichmuller` is additionally an element
method here where upstream exposes it only on the ring (`padic_generic.py:484`); the ring version
also exists, so the element method is a port-only extra.

### How to close

Port the remaining precision-model classes and the extension-element templates from
`reference/sage/src/sage/rings/padics/` (`CA_template.pxi`, `FM_template.pxi`,
`padic_ext_element.pyx`, `unramified_extension_generic.py`,
`eisenstein_extension_generic.py`). For the L-series, port `sage.modular.modsym` and Monsky-Washnitzer
cohomology. **Effort: large for both.** `pAdicEisensteinQuadraticExtension` in `padic_lseries.ts`
should move into `rings/padics/` once `pAdicExtension` has an element type.

### Trade-offs of leaving it open

Only one precision model; no `Zq`; no p-adic L-function values. `padic_elliptic_logarithm` (see
[Elliptic Curves over Q](#elliptic-curves-over-q-and-number-fields)) is blocked on the same layer.

### Behavioral Impact

Basic arithmetic matches SageMath for capped-relative elements: `Zp(7,10)` `a`, `a+a`,
`R(2).square_root()`, `teichmuller(3)`, `R(8).log()` all agree character for character, as do
`minimal_polynomial()` and `charpoly()` (`padic_generic_element.ts:1258`, `:1276`). Extension-field
operations and the L-series throw.

---

## Arithmetic Functions Not Delegated to PARI/FLINT

| Function | SageMath | sagemath-ts | How to close |
|----------|----------|-------------|--------------|
| `hilbert_symbol(a, b, p, algorithm)` | Rational square-class normalization; PARI by default, direct/all selectable | Repaired in 8.12.0: RationalLike coefficients, IntegerLike place, PARI delegation and original algorithm errors; conductor helpers normalize wrapped inputs | See [Hilbert symbol dependency domain](#hilbert-symbol-dependency-domain) for the explicit native input subset |
| `bernoulli(n, algorithm)` | Multiple backends (FLINT/Arb/PARI/bernmm) with heuristics; `B_1000` is instant | `misc.ts:2149-2187` is the **classical binomial recurrence** `B_m = -1/(m+1)·sum C(m+1,k)·B_k` over rationals with a gcd after every term — *not* the Akiyama-Tanigawa the comment at `:2135` claims, and the Von Staudt-Clausen comment at `:2113-2114` is dead text. `algorithm` is accepted and ignored ('pari'/'flint'/'arb'/'bernmm'/'gap' all return `-691/2730` for `n = 12`). The operative limit is not the safe-integer range but the `O(n^2)`-rational cost: measured `B_200` 231 ms, `B_400` 2.66 s, `B_1000` 76.1 s, so the practical ceiling is ~`n = 500`. Returns a plain `{numerator, denominator}` object, not the port's `Rational` class | Port FLINT's exact multi-modular path: `reference/flint/src/arith/bernoulli_number_denom.c` (68 lines, von Staudt-Clausen) + `bernoulli/mod_p_harvey.c` + `fmpq_ui_multi_mod.c` + CRT. Avoids the Arb/zeta path entirely and is exact integer arithmetic, so it satisfies CLAUDE.md's no-floating-point rule. **Effort: moderate, ~400 lines.** Also fix the algorithm comment and the return type |
| `dedekind_sum(p, q, algorithm)` | FLINT by default, PARI when selected, otherwise `ValueError('unknown algorithm')` | Repaired in 8.8.0: delegates to the corresponding dependency port, including negative/zero moduli, noncoprime reduction, large-word boundary behavior and IntegerLike coercion | See [Dedekind sum backend arithmetic](#dedekind-sum-backend-arithmetic) for remaining native optimization differences; 8,298 original comparisons |
| `gauss_sum(char_value, finite_field)` | Generic over `char_value.parent()`; the doctests feed it `UniversalCyclotomicField().zeta(q-1)` | `misc.ts:4793-4853` is a near-verbatim transcription of `misc.py:6428-6446` (same `resu += zq_power * zeta_p_powers[gen_power.trace().lift()]`, same `gen_power *= gen`, `zq_power *= zeta_q`) — **not** a "numeric-only fallback". The real blocker is that **no ring in this repo satisfies the `CharacterValue` interface** at `:4760-4764`: there is no `UniversalCyclotomicField` anywhere in `packages/sagemath-ts/src`, so the function has zero call sites and zero tests, violating CLAUDE.md's mandatory property-test rule. Two undocumented hazards: `:4838` `if (resu.add && zq_power.mul)` silently **skips terms** and returns `ring.zero()` when a ring lacks the optional methods, and `:4817` seeds `zq_power` from `ring.zeta(1n).powers(1)[0]` where Sage uses `ring.one()` | Make the missing-method branch **throw** (immediately — it violates the [Honest Failure](#honest-failure-instead-of-silent-approximation) policy), fix the `zq_power` seed, and either port a cyclotomic ring or delete `gauss_sum` until a usable ring exists. `reference/sage/src/sage/rings/universal_cyclotomic_field.py` is vendored (1740 lines) but is libgap-backed, so it is not a straight port. **Effort: minutes for the throw; large for the ring** |
| `algebraic_dependency(z, degree, options)` | Exact integer/rational shortcuts; RR/CC precision-controlled lattice; RDF and Python floats delegate to PARI; every branch selects an irreducible factor | Repaired in 8.10.0 for IntegerLike, Rational, number and RealNumber inputs, including height/proof dispatch and real-field method forwarding | See [PARI algebraic dependencies](#pari-algebraic-dependencies); complex/p-adic domains and native reduction/interval details remain open |
| `Integer.class_number()` | See [Quadratic Class Numbers](#quadratic-class-numbers-not-delegated-to-buchquad) | — | — |
| `number_of_partitions(n)` | FLINT, handles large `n` | `integer_ring.ts:2001` throws above `n = 10 000`; `p(10000)` matches `Partitions(10000).cardinality()` | Replace the `O(n^2)` DP with FLINT's exact pentagonal-number recurrence, `reference/flint/src/arith/number_of_partitions_vec.c` (~60 lines, no Arb needed), and delete the cap. **Effort: small.** Matching Sage for astronomically large `n` (HRR via Arb) stays out of reach and is the only part worth documenting |
| `prime_pi(n)` | primecount / FLINT `n_prime_pi` | `integer_ring.ts:2073` throws above `10^7`; `:2077-2082` is a per-integer `_is_prime` loop — exactly the naive pattern CLAUDE.md's Algorithm Fidelity section forbids (433 ms at 10^6, so ~6 s at the cap). `prime_pi(100000) = 9592` matches Sage | Replace with FLINT's sieve: `reference/flint/src/ulong_extras/prime_pi.c` + `prime_pi_bounds.c` (~200 lines). **Effort: small-moderate.** Note the module attribution: Sage has these in `combinat/partition.py:9762` (via `libs/flint/arith_sage.pyx:175`) and `functions/prime_pi.pyx:39` (primecountpy); `prime_pi` is **not** an `Integer` method in Sage at all |

### Trade-offs of leaving these open

Large Bernoulli numbers, any Gauss sum,
large partition counts and large prime counts are unavailable or wrong.
`gauss_sum` additionally violates the port's own honest-failure policy.

### Behavioral Impact

Dedekind algorithm selection and signed/zero denominators are repaired. Algebraic-dependency height misses now return `null`. `gauss_sum` can silently return
`ring.zero()`. Everything else raises. Integer reciprocal inversion is now repaired.

---

## Finite Fields — Conway Table and Constructor Algorithms

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Conway polynomial availability | The `conway_polynomials` package: every `(p, n)` in Frank Luebeck's `CPimport.txt` | A hand-pasted subset in `conway_polynomials.ts:36` `CONWAY_LOW_COEFFICIENTS`: p = 2 to n = 64, p = 3 to 24, p = 5 to 18, p = 7 to 14, p = 11/13 to 12, p = 17/19/23/29/31 to 10. `has_conway_polynomial` returns false where SageMath's does not — Sage's own `exists_conway_polynomial` is `True` for `37^2`, `97^2`, `787^5`, `2^64`, `2^65`, `2^100`, `19^21`, `3^24`, `3^25`, none of which we have |
| Default modulus for `GF(p^n)` | `irreducible_element` (`polynomial_ring.py:3560-3626` -> `:2628-2681`): `n == 1` -> `x − 1`; Conway if available; NTL `GF2X_BuildSparseIrred` for `p = 2`; else PARI `ffinit` | **Identical branch chain**, delegating to `ntl-ts` and `parigp-ts`. Outside the tabulated range the modulus is `ffinit`'s / `BuildSparseIrred`'s rather than the Conway polynomial Sage would use, so element representations are not interoperable with SageMath's and the generator need not be primitive |
| `modulus=` / algorithm keyword on the exported `GF`/`FiniteField` | `GF(p^n, 'a', modulus='minimal_weight')` | Repaired in 0.0.26: `GF(q, {name, modulus})` and `GF(q, name, {modulus})` accept algorithm strings, full coefficient lists or polynomials. Lists include the leading term, unlike the existing lower-coefficient `GFpn` shorthand. Prime fields honor custom linear roots; names and degrees are validated |
| `algorithm='ffprimroot'` | `self(pari(p).ffinit(n).ffgen().ffprimroot().charpoly())` | Throws naming `ffgen`/`ffprimroot`/`charpoly`: `parigp-ts/src/ff.ts` is `F_p`-only with no finite-field element type. It throws rather than returning some other irreducible polynomial, which would silently not be primitive |
| `algorithm='random'` over GF(2) | NTL `BuildRandomIrred` | Rejection sampling from `current_randstate()` — exactly the fallback SageMath itself takes when its NTL import fails (`polynomial_ring.py:3615-3620`). The distribution differs from NTL's and the concrete polynomial for a given seed differs from Sage's. See [ntl-ts GF2X Factoring Stubs](#ntl-ts-gf2x-factoring-stubs) |

### How to close

- **Conway table:** the data **is** vendored at
  `reference/flint/src/nmod_poly/conway_polynomial_data.c` (322 KB; `__nmod_poly_cp_primes0` covers
  all 55 primes 2..257, degrees to 409 for `p = 2`), and the decoder
  (`conway.c::conway_polynomial_lt_260`) is already ported once. Write a one-off script to regenerate
  the table for all `p <= 257` and re-run the existing irreducible/primitive/subfield-compatibility
  validator. **Effort: low-moderate** — one script plus a generated ~300 KB TS data file. Full parity
  with Sage's `conway_polynomials` package (e.g. `787^5`) still needs Luebeck's `CPimport.txt`, which
  is **not** vendored — that residue is the only honest deviation left here.
- **Factory keywords:** `name` and `modulus` are repaired. Tuple-order shorthand and broader
  keywords (`check_prime`, `check_irreducible`, `implementation`, `proof`, `names`, `prefix`,
  `repr`, `elem_cache`) remain unsupported. The default generator name remains the port
  convention `a`, instead of Sage's unnamed algebraic-closure/pseudo-Conway path.

### Trade-offs of leaving these open

`GF(p^n)` element representations diverge from SageMath's outside the tabulated range.
Unsupported factory keywords now raise instead of silently becoming a generator name.

### Behavioral Impact

Within the Conway database the polynomial is the real Conway polynomial and everything matches:
45/45 explicit-algorithm moduli match Sage exactly, and of 130 sampled default moduli there are **0
genuine mismatches** — the 74 divergences were each confirmed with Sage's `exists_conway_polynomial`
to be cases where Sage has a Conway entry our table lacks. `GF(1009^8)`, where Sage also has no
Conway entry, matches byte for byte.

> **Note for future maintainers:** `GF(2^8).modulus()` is `x^8 + x^4 + x^3 + x^2 + 1`, the **Conway**
> polynomial (`finite_field_givaro.py:69`), *not* NTL's `BuildSparseIrred(8)` value
> `x^8 + x^4 + x^3 + x + 1` (the Rijndael polynomial). Conway wins Sage's branch order. The NTL value
> is reachable via `algorithm='minimal_weight'`. Do not "fix" this.

---

## Polynomials — Printing, Factor Shape, Term Orders and Base Rings

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `Polynomial.toString()` | `x - 1`, `x^3 - 2*x + 5`, `-x` | `x + -1`, `x^3 + (-2)*x + 5`, `(-1)*x`. Negative coefficients are not folded into the sign. This matters because other sections assert that every listed `toString()` reproduces SageMath's printed form |
| `factor(n)` for `n < 0` (integers) | Keeps the unit out of the factor list: `list(factor(-360)) == [(2,3),(3,2),(5,1)]`, `unit() == -1` | Returns `[[-1,1],[2,3],[3,2],[5,1]]` — the unit is an entry. (The polynomial `factor()` does the same deliberately, to restore `prod(factors) === f`; for integers there is no such justification, since the caller can multiply the sign back in trivially) |
| Term orders | `TermOrder` supports thirteen orders (lex, invlex, deglex, degrevlex, neglex, negdegrevlex, negdeglex, degneglex, wdeglex, wdegrevlex, negwdeglex, negwdegrevlex) plus block orders and weighted gradings | `multi_polynomial_element.ts:29` declares `export type TermOrder = 'lex' \| 'deglex' \| 'degrevlex'`; `:52` throws `unknown term order '<name>'`, matching `term_order.py:796` exactly. Weighted grading and `degree(std_grading=…)` have no meaning here. Roughly 19 of Sage's `MPolynomial` methods are honest `NotImplementedError` stubs |
| `MPolynomialRing` over ZZ | Works | **Cannot be constructed at all.** `MPolynomialRingConstructor(ZZ, ['x','y'])` dies inside `gens()` with a raw JS `TypeError: coeff.isZero is not a function` (`multi_polynomial_ring.ts:154`) before any deviation check runs. The only "ZZ" exercised in the test suite is a hand-rolled mock (`multi_polynomial_ideal.test.ts:130` `makeZZ`) |
| `groebner_basis()` over Zmod(n) | `MPolynomialIdeal.groebner_basis` works over ZZ and Zmod(n) via Singular's `std` (doctests at `multi_polynomial_ideal.py:4593-4610`) | Raises Sage's `TypeError('Can only reduce polynomials over fields.')` (`multi_polynomial_ideal.ts:124`, matching `multi_polynomial_element.py:2487-2488`). Without the field guard the quotient coefficient rounds to zero, the subtrahend is zero, and the reduction loop never terminates (reproduced as a timeout) — so the guard is correct *given* the missing Singular backend, but Sage does answer here. `dimension()` over a non-field raises `NotImplementedError('implemented only over fields')` |
| Gröbner engine | Singular/FGb with F4/F5-style algorithms | Naive Buchberger with an iteration budget (see [Bounded Search Budgets](#bounded-search-budgets-and-measured-thresholds)). `dimension()` is the ported Cox-Little-O'Shea algorithm (`multi_polynomial_ideal.py:1128-1192`); all five of Sage's `dimension()` doctests reproduce (`1, -1, 1, 1, 2`) |

### How to close

- **Printing:** fold negative coefficients into the sign in `Polynomial.toString`, emitting Sage's
  `' - c*x^k'` form. **Effort: small**, and it removes a false claim elsewhere in this document.
- **Integer `factor()` unit:** reconcile this with the documented array/unit adaptation.
  `factor(0)` already raises `ArithmeticError`; the former claim of `ValueError` was stale.
- **Term orders:** transcribe `sortkey_invlex`, `neglex`, `negdegrevlex`, `negdeglex`, `degneglex`,
  `wdeglex`, `wdegrevlex`, `negwdeglex`, `negwdegrevlex` from
  `reference/sage/src/sage/rings/polynomial/term_order.py:965-1181` into `getTermOrderComparator` —
  each is 15-20 self-contained lines with no dependency beyond the exponent tuple.
  **Effort: ~150 lines, under a day.** Block orders, weighted gradings and a real `TermOrder` object
  are a separate, larger job that can stay documented.
- **`MPolynomialRing` over ZZ:** make the coefficient interface accept `bigint`/`Integer`, or wrap ZZ
  elements. **Effort: small**, and it is basic functionality that is currently documented nowhere.
- **Gröbner over ZZ/Zmod(n):** needs Singular's `std` or an equivalent; **effort: large** — this is
  the one row here that is genuinely a backend gap.

### Trade-offs of leaving these open

Printed polynomials do not round-trip against Sage's repr; nine term orders are unavailable and
silently unreachable (the rejection is correct — degrading to degrevlex would produce wrong Gröbner
bases with no error); multivariate polynomials over ZZ cannot be built.

### Behavioral Impact

Different printed forms; a different `factor()` shape and error class for integers; `TypeError`/
`NotImplementedError` where SageMath computes.

---

## Matrices — the matrix() Constructor

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `matrix(...)` constructor | `sage/matrix/constructor.pyx` accepts a flat list with `nrows`/`ncols`, dicts, callables, sparse flags and a bare `(nrows, ncols)` form | `matrix_space.ts:268` is `export const matrix = MatrixFromEntries;`, i.e. the only accepted signature is `(ring, entries: R[][])`. `matrix(QQ, 2, 2, [1,2,3,4])` — the most common SageMath spelling, used throughout `sage/matrix`'s own doctests — dies with `TypeError: undefined is not an object (evaluating 'entries[0].length')` at `matrix_space.ts:250` |

### How to close

- **`matrix()`:** implement the flat-list + `nrows`/`ncols` form from the vendored `constructor.pyx`.
  **Effort: small**, and it removes a constant source of friction when transcribing doctests.

### Trade-offs of leaving these open

Every doctest transcription using the common flat-list constructor must rewrite the call by hand.

### Behavioral Impact

One unusable constructor signature.

---

## Lattices — Exact SVP Rank Cap

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `IntegerLattice.shortest_vector` | fpylll's `SVP.shortest_vector`, or PARI `qfminim` with `algorithm='pari'` | Exact Fincke-Pohst enumeration in bigint/rational arithmetic for rank `<= EXACT_SVP_MAX_RANK = 30` (`free_module_integer.ts:1540`). Above rank 30 it **silently returns `rows[0]` of the LLL basis** (`:505-511`), a `2^((n−1)/2)` approximation, with no error — the one place in the port where a budget silently degrades a result. `algorithm='pari'` is accepted and then ignored (`:490-493`) |

### Why this is still open

The delegation target now effectively exists: `packages/parigp-ts/src/qfrep.ts` is a line-cited port
of `bibli1.c:1299-1462` `minim0_dolll` — the Fincke-Pohst core of `qfminim` — plus `minim_lll` and
`forqfvec_init_dolll`. The earlier rationale ("no external dependencies — without fpylll or PARI
bindings") is contrary to CLAUDE.md's delegation rule and no longer true.

### How to close

Export a `qfminim` from `parigp-ts` on top of the existing `minim0_dolll` (a vector-returning flag,
`min_FIRST` / `min_ALL`, plus a thin wrapper), have `shortestVector` delegate to it and honour
`algorithm='pari'`, and drop the rank-30 cliff.
**Effort: moderate, a few hundred lines** — `qfrep.ts` already contains the enumeration and the
LLL-reduced Cholesky, so no new algorithm is needed.

### Trade-offs of leaving it open

Above rank 30 the answer is an approximation presented as an answer. Below rank 30 it is exact and
matches SageMath.

### Behavioral Impact

Silently non-shortest vectors above rank 30. `algorithm='pari'` is a no-op.

---

## Parents Are Not Unique

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Parent identity | `UniqueRepresentation`: `GF(7) is GF(7)` is `True`, and Sage's coercion framework is built on `parent is parent` checks | `GF(7) === GF(7)` is `false`. Repeated GF and polynomial/series constructor calls return distinct, non-identical parents, so `is`-style comparisons that Sage relies on cannot be reproduced. `PowerSeriesRing` works around it with a structural `is_identical_to` (see [Power Series](#power-series-laurent-series-and-multivariate-series)) |
| Caching generally | Extensive (`@cached_method`, `UniqueRepresentation`) | Minimal (zero/one caching, selective caches). Zmod/Integers/Mod now share the signed-order weak factory cache, matching Sage identity while the parents are live |

### Why this is still open

This was previously recorded as a performance choice ("avoid heavy caching until profiling indicates
need"). It is not: parent identity is **observable**, and every module that needs it has had to
invent its own structural comparison.

### How to close

Add a parent cache keyed on the constructor arguments (the analogue of `UniqueRepresentation`) for at
least `GF`, `FiniteField`, `PolynomialRing`, `PowerSeriesRing`, `NumberField` and `MatrixSpace`, and
replace the ad-hoc `is_identical_to` helpers with `===`.
**Effort: moderate** — mechanical, but it touches every ring constructor and needs care with
equal-but-distinct argument objects.

### Trade-offs of leaving it open

Structural comparison must be reimplemented per module; repeated construction is slower; deep imports
and root re-exports can hand back different objects for the same ring.

### Behavioral Impact

`===` on parents is unreliable. Mathematical values are unaffected.

---

## Real and Complex Precision and Rounding

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `RealField(53, rnd='RNDD').rounding_mode()` | Returns the string `'RNDD'` | Returns a `RoundingMode` **enum member** (a number). There is no string-returning accessor |
| Directional rounding | MPFR applies the field's rounding mode to **every** primitive | Applied to `div` and `sqrt` (`real_mpfr.ts` `applyRounding`, which re-rounds the correctly rounded RNDN result after an exact comparison). The other primitives still round to nearest regardless of the field's mode |

### How to close

Add a string-returning `rounding_mode()` and route the remaining primitives (`add`, `sub`, `mul`,
`exp`, the transcendentals) through `applyRounding`. **Effort: small for the first, medium for the
second** — each primitive needs an exact error-sign oracle, which only `div` and `sqrt` have today.

### Trade-offs of leaving it open

A non-default rounding mode is honoured by division and square root but silently ignored elsewhere.

### Behavioral Impact

The stored value is a double either way — see
[No Arbitrary-Precision Floating Point](#no-arbitrary-precision-floating-point) for that.

---

## Power Series — V(0)

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `V(0)` on an exact (infinite-precision) power series | `SignError: cannot multiply infinity by zero` — an incidental failure in the precision bookkeeping | Returns `0`. **Mathematically wrong**: `V(0)` of `1 + 2x + 3x^2` is the constant `6`. So this is a real defect *and* a case that cannot be pinned against the oracle. The finite-precision case is pinned and passes |

### How to close

Compute the correct constant (the sum of the coefficients) rather than returning zero.
**Effort: trivial.**

### Trade-offs of leaving these open

`V(0)` returns a wrong value silently.

### Behavioral Impact

One wrong value.

---

## Coding and Crypto — Permissive Where Upstream Raises

The port is more permissive than SageMath in several places. Being more permissive is not
automatically wrong, but it is undocumented divergence, and one of these is a naming trap.

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `ReedMullerCode.decode` / `ReedSolomonCode.decode` | `C.decode_to_code(v)` returns a length-`n` **codeword**; `C.decode_to_message(v)` returns the length-`k` message | `decode(received)` returns the length-`k` **message** (i.e. Sage's `decode_to_message`); re-encoding it gives the codeword. A caller expecting a codeword silently gets a shorter vector |
| `SBox.to_bits(x)` without an explicit width when `m != n` | `sbox.pyx:271-303`: `n` stays `None` unless `m == n`, so `ZZ(x).digits(base=2, padto=None)` is called and SageMath raises `TypeError` | `sbox.ts:311` always pads to `this._n`, so `SBox([0,1,1,0]).to_bits(3)` returns `[1]` where SageMath raises |
| `Regev(1)` | Builds the oracle with `q = 2` and `sigma = +inf` | Throws `ValueError('sigma must be a finite number, got Infinity')`. A degenerate `n = 1` corner where SageMath produces an unusable object |
| Generalized Reed-Solomon codes with `k = n` | `parity_check_matrix()` and every syndrome-based decoder raise `ValueError`, because `dual_code()` needs a positive dimension (`grs_code.py:239` via `:476`) | `parity_column_multipliers()`, `syndrome()` (empty) and `decode()` all succeed |

### How to close

- **`decode`:** rename to `decode_to_message` and add a `decode_to_code` that re-encodes. **Effort:
  small**, and it removes the fidelity trap. (Alternatively, register it permanently in
  [Return Shapes](#return-shapes-keyword-arguments-and-signature-adaptations) — but the name is the
  problem, not the shape.)
- **The other three:** decide per case whether to reproduce upstream's raise or to keep the port's
  behaviour, and record the decision. Reproducing a Python `TypeError` from `Integer.digits` is not a
  useful contract, so `to_bits` and `Regev(1)` are probably keepers; the GRS `k = n` case should
  probably raise as upstream does, since a syndrome decoder for a code with no dual is meaningless.
  **Effort: trivial each.**

### Trade-offs of leaving these open

`decode()` is a silent shape trap. The others are only unregistered.

### Behavioral Impact

A shorter vector than a Sage-trained caller expects from `decode()`; success where SageMath raises in
four places.

---

## Hyperelliptic — Frobenius Polynomial Algorithms

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `frobenius_polynomial()` algorithm selection | `hyperelliptic_finite_field.py:616-630` picks `'matrix'` (hypellfrob) when the base field is prime, large enough, and the model is odd-degree with `h = 0`; otherwise `'pari'` (hyperellcharpoly) in odd characteristic; otherwise `'cardinalities'` | Always `'cardinalities'`. `'matrix'` and `'pari'` are `NotImplementedError` stubs naming `hypellfrob` and `parigp-ts`'s missing `hyperellcharpoly`. The same applies to `count_points()` and `cardinality()`, which always take the exhaustive path |

### Why this is still open

Neither dependency exists in this repo. The characteristic polynomial of Frobenius is uniquely
determined by the curve, so the **value** returned is identical to SageMath's — this is a cost
deviation, not a behavioural one. But the running time is `O(q^g)` instead of polynomial, so large
base fields (e.g. the `GF(3663031)` doctest in `jacobian_generic.py:429`) are out of reach.

### How to close

Port `hypellfrob` (Kedlaya's algorithm, `sage/schemes/hyperelliptic_curves/hypellfrob/`) and/or
PARI's `hyperellcharpoly` (`reference/pari/src/basemath/hyperell.c`) into `parigp-ts`, then restore
Sage's dispatch. **Effort: large.**

### Trade-offs of leaving it open

Curves over large base fields are uncomputable.

### Behavioral Impact

None on values where computable; `O(q^g)` running time.

---

## Quaternion Algebras — Base Rings Other Than QQ

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Base rings | `QuaternionAlgebra(a, b)` accepts elements of any ring in which 2 is a unit (GF(p), number fields, `Frac(QQ[x])`, Laurent polynomial rings, …), with `QuaternionAlgebraElement_generic` / `_number_field`; orders and ideals are QQ-only | **QQ only.** Anything else raises `NotImplementedError('SAGE_NOT_IMPLEMENTED: quaternion algebras over base rings other than QQ')`; `QuaternionAlgebraElement_generic` and `QuaternionAlgebraElement_number_field` are named `NotImplementedError` stubs |

### Why this is still open

The order/ideal/maximal-order machinery this module targets (Deuring correspondence, isogeny-based
crypto) is QQ-only in SageMath too. Supporting other base rings needs number-field ideals, real
embeddings and PARI `alginit`, none of which is ported.

### How to close

Port number-field ideals and real embeddings (largely the same `nf` layer as
[Number Field Class Groups](#number-field-class-groups-units-and-galois-closure)) and PARI's
`alginit` (`reference/pari/src/basemath/alglin*.c`). **Effort: large for the base rings.**

### Trade-offs of leaving it open

No GF(p) or number-field quaternion algebras, hence no `ramified_places` over number fields.

### Behavioral Impact

None for QQ, hard error (never a wrong answer) elsewhere. `minimal_element` returns a different
element of minimal norm in 4 of 254 sweep cases; the norm is always the true minimum.

---

## Polynomial String Parsing

| Aspect | SageMath | sagemath-ts |
|---|---|---|
| Grammar | `misc/parser.pyx` also parses function calls, lists, tuples and matrices | Arithmetic expressions and comparisons are ported. General function calls and sequence/matrix entry points explicitly raise `SAGE_NOT_IMPLEMENTED` |
| Parser construction | Python operators automatically dispatch on constructor results | The internal TypeScript Parser requires integer/float/name constructors and explicit operator callbacks; the first five argument positions follow Sage, and a sixth callback object supplies binary/unary operations. Python-default constructors are not supplied |
| Fraction representation | Polynomial division uses FractionFieldElement or the small-prime FpT specialization | String evaluation shares the polynomial expression parser and fraction arithmetic through a temporary function-field wrapper. Polynomial results return to the polynomial ring; fraction constructors retain the fraction result. Public fraction classes are available |
| Symbolic and non-field expressions | Sage can promote through symbolic rings or fraction fields | The implemented string path is verified over prime-field coefficient rings. General symbolic factorials, rational exponents and promotions from arbitrary coefficient rings remain open and are not claimed equivalent |
| Backend display | GF(2) polynomial rings print `using GF2X`; sufficiently large prime fields print `using NTL` | These Sage representation labels are now preserved. They do not establish that every polynomial operation delegates to that backend; the separate algorithm-delegation gaps remain |

**Rationale:** Explicit callbacks supply the operator dispatch unavailable in JavaScript.
The shared evaluator uses the existing fraction normalization and arithmetic. Its temporary
function-field wrapper still limits non-field coefficient promotions.

**Trade-offs:** This is an arithmetic parser for the implemented polynomial construction path,
not a complete symbolic parser or evaluator. Function-call, sequence and matrix stubs remain
visible, and the added keyword/constructor variants do not imply universal coefficient support.

**Behavioral impact:** Prime-field polynomial strings now match the tested Sage path, including
implicit multiplication, exponent precedence, numeric factorials, comparisons, cancellation,
QQ coercion rejection and backend-specific error messages. Outside the Conway table the
extension-string oracle supplies Sage's PARI-generated modulus explicitly, matching the port's
documented fallback rather than comparing different field presentations. Parser tokenization
uses Python's whitespace and Digit categories; Unicode-version differences outside the tested
character repertoire remain possible.

**How to close:** Port the remaining `misc/parser.pyx:573-725,976-1044` entry points and the
missing symbolic and general coefficient promotions, then compare their constructors and error paths.

---

## Modular Integer Coercion and Factories

Sage's integer-mod element constructor accepts exact integers, rationals and explicit lifts;
its arithmetic invokes the canonical coercion model. The port now handles Integer wrappers,
booleans, null/default zero, integer strings and Rational conversion, and preserves the
native-int/int64/GMP negative-power error boundaries. Binary operations use direct quotient
maps or the gcd-modulus pushout; a newly formed trivial quotient is rejected. A compatible
prime-field operand promotes the result to its field. Incompatible equality returns false.
Seeded ring randomness now uses Python randint, including the optional symmetric bound.

**Rationale and architecture:** TypeScript has named arithmetic methods and separate prime
field classes. Overloads preserve IntegerMod results for integer/modular operands and expose
the prime-element type for prime-field operands and extension types for extension operands.
Sequence multiplication now returns strings/lists as Python does; earlier tests incorrectly
called result.parent() and mistook these valid results for AttributeError. The corrected
oracle uses Sage parent(result), with permanent mixed-operand and sequence regressions.
Python operator and parent displays are
compared explicitly. The `repr()` convenience retains `Mod(value, modulus)` although Sage's
repr is just the residue; `toString()` matches Sage. Factory parents are weakly interned by
the original signed order; direct generic-ring constructors are not interned.

**Trade-offs and behavioral impact:** The shared modular cases compare construction and arithmetic,
including exact errors and result parents, but do not implement Sage's general symbolic,
p-adic, custom conversion-map or category framework. Factory category/is_field options remain
unimplemented. Zero/default Zmod orders return ZZ, negative orders normalize, and Mod with
zero returns its original input even with an explicit parent. An explicit modular or prime
parent overrides a nonzero modulus. Dynamic factory return types expose these alternatives.
Integral JavaScript numbers still denote integers; the zero-Mod oracle converts integral
Python floats to Integer to respect that documented language mapping.

**PARI delegation:** Modular multiplicative order now delegates to `parigp-ts/znorder`,
mirroring PARI's prime-power valuation and lcm algorithm. The API accepts separate bigint
residue/modulus arguments in place of a GEN t_INTMOD. Its valid-input comparisons include
modulus one, powers of two, odd prime powers, mixed composites and a 127-bit prime. Optional
order means a positive known multiple; PARI exception objects and arbitrary invalid GEN
inputs are outside this adapter's domain.

---

### How to close the remaining coercion/category gaps

Port `integer_mod_ring.py:205-230` category refinement and the general exact/symbolic conversion
paths in `integer_mod.pyx:377-404`. **Effort: large**, because category and symbolic parents are
not implemented. LWE/RingLWE continue to require a nonzero modulus even though Zmod itself
supports zero; their element/vector/sampler types still assume IntegerMod. Supporting the ZZ
case from `crypto/lwe.py` requires those consumers to handle native integer elements too.
**Effort: moderate.**

---

## Template for New Deviations

Copy this template when adding a new entry, and put it in **Part I** only if the rationale would
still hold after every dependency is ported. Otherwise it belongs in **Part II**, with an upstream
file and a "How to close".

```markdown
## [Deviation Title]

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Description | what Sage does | what we do |
| Affected modules | upstream path | port path |

### Rationale

1. **Reason 1** - Explanation
2. **Reason 2** - Explanation

### Trade-offs

- What we lose by deviating

### Behavioral Impact

Does this change outputs? Edge cases? Error messages?
```

For a Part II entry, replace *Rationale* with:

```markdown
### Why this is still open

(A cause, not a justification.)

### How to close

Port `reference/<lib>/<file>:<line>`. **Effort: trivial | small | moderate | large.**
```


## Extension Arithmetic and PARI Quotient Kernels

| Aspect | SageMath/PARI | sagemath-ts |
|---|---|---|
| Extension arithmetic | element_pari_ffelt calls FF_add/sub/neg/mul/div/inv/pow, which selects an FpXQ, Flxq or F2xq kernel | Polynomial-basis elements now delegate to parigp-ts FpX/FpXQ kernels. Local polynomial multiplication, xgcd and exponentiation are no longer the extension backend |
| Signed quotient powers | FpXQ_pow inverts before negative powers; exponent zero returns polynomial one immediately | Both entry points share the native early returns and unreduced initial bases. Word powering delegates below 2^63; separate square kernels and bb_group.c's exact exponent windows are shared |
| Inverse signature | FpXQ_inv(a,T,p) computes an inverse over F_p; ZpXQ_inv(a,T,p,e) lifts to p^e | FpXQ_inv now supports the original three arguments. Its prior four-argument (a,T,q,p) Hensel convenience remains supported for q=p^e and is compared directly with ZpXQ_inv |
| Exponent coercion | PARI extension elements compare exp with zero before Integer conversion; zero to a negative noninteger raises ZeroDivisionError first | Extension pow accepts unknown for validation and follows the same order. Null/string/list ordering errors, Rational and nonfinite-number errors are compared exactly; integral JS numbers retain the documented Integer mapping |
| Representation and performance | PARI packs small-prime/binary polynomials and selects optimized multiplication kernels | The existing dense bigint FpX representation is retained for all characteristics. FpX/Flx now use the native packing, square, division and GCD schedules described below, with documented bounds on packed intermediates; complete FF binary/extension specialization remains separate work |
| Low-level errors | PARI e_INV includes GEN rendering and native backend names | FpXQ_pow, its power alias and FpXQ_inv now preserve PariError, native backend names and polynomial payloads. Their direct and Hensel comparisons use complete bundled-native diagnostics without e_INV normalization. Other legacy error adapters remain separate work |

**Rationale:** Use the dependency port for Sage's delegated operations while preserving the
existing TypeScript field representation and public four-argument inverse convenience.
Share signed-power behavior and PARI's exponentiation strategy across the two entry points.

**Trade-offs:** Dense public coefficient arrays use bounded packed kernels internally;
native GMP behavior beyond that bound is not reproduced. Other legacy low-level
error adapters remain outside the audited quotient kernels.
These are explicit architecture/performance/diagnostic boundaries, not full backend parity.

**Behavioral impact:** Numeric quotient results, noninvertibility, zero-exponent ordering and
Sage-facing exponent-conversion errors now agree in the comparative cases. Tests call the
installed PARI C functions via install, including reducible and nonmonic quotient moduli,
rather than assuming every quotient is a field. The extension oracle's negative-infinity
fixture was also corrected to preserve the sign before comparison.

## Dense FLINT polynomial GCD kernels

| Aspect | SageMath / FLINT | TypeScript port |
|---|---|---|
| Integer GCD dispatch | Sage calls FLINT; FLINT strips powers of x and selects subresultant PRS, packed heuristic or modular reconstruction | Sage delegates to flint-ts `_fmpz_poly_gcd`, with the same dispatch, content/sign handling, heuristic verification, CRT bounds and unlucky-prime handling |
| Representation | Caller-owned fmpz/nmod buffers with explicit lengths and native limbs | Exported kernels take readonly bigint arrays, infer and normalize lengths, allocate their result and never modify inputs. Zero is `[]`; `_fmpz_poly_gcd_heuristic` returns `null` when its heuristic fails |
| Modular GCD | Word arithmetic, Euclid below tuned thresholds and recursive half-GCD above them | The same degree/word-size cutoff tables and half-GCD recursion; native BigInt replaces word reduction. The raw `_nmod_poly_gcd` result need not be monic; callers normalize it |
| Multiplication/division | Native classical, KS/KS2/KS4/FFT and division dispatch | Small classical products and larger exact Kronecker products use BigInt multiplication. Modular division uses reversed Newton inversion above the classical cutoff; integer divisibility uses divide-and-conquer with cutoff 16. Native FFT/REDC and additional packing dispatch remain unported |
| Prime selection | FLINT's modular-prime table, mod-30 wheel and native hashed Miller–Rabin witnesses | Same 64-bit modular-prime table and wheel; beyond the table a deterministic seven-base unsigned-64-bit Miller–Rabin check replaces the native witness lookup |
| Low-level errors | Native preconditions and FLINT exception/abort mechanisms | Prime moduli and valid bigint coefficient arrays are preconditions. Invalid low-level modulus/inversion inputs use RangeError; this is not a port of the C error system |

**Rationale:** Keep Sage's algorithm and dependency boundaries while using native BigInt
storage instead of requiring the still-stubbed mutable FLINT object APIs. Fast products,
recursive division and half-GCD retain the principal asymptotic algorithms.

**Trade-offs:** Native instruction-level tuning, REDC, FFT and multi-product Kronecker
variants are not reproduced. Runtime crossover performance can differ substantially;
this is not a claim of complete FLINT implementation or native performance parity.

**Behavioral impact:** Integer GCD values and the heuristic success/failure result are
compared directly against Sage's loaded FLINT C library. Modular comparisons normalize
raw remainders to the public monic result. Sage's own zero/one shortcuts remain above
FLINT and preserve their original operand sign and identity.

## Polynomial Monic Normalization

| Aspect | SageMath | TypeScript port |
|---|---|---|
| ZZ result parent | A nonmonic integer polynomial is divided by its leading coefficient in QQ; an already-monic polynomial is returned unchanged | Same coefficients, parent extension and identity, including a leading coefficient of -1 |
| Native modular backend | `polynomial_zmod_flint.monic` checks that the leading coefficient is a unit, then calls FLINT and allocates a result | Word-sized modular/prime coefficients (except GF2X) delegate to `_nmod_poly_make_monic` in flint-ts. Already-monic inputs also allocate a result |
| Zero input | QQ/ZZ and generic finite backends raise their coefficient inversion error; Sage 10.3's word-modular backend raises ValueError | Preserve those executed Sage 10.3 errors. The bundled newer FLINT zero-normalization error boundary remains unverified on a Sage 10.9 runtime, as described under Finite Field Representation above |
| Monicity predicate | Zero is not monic; the leading coefficient is compared with the base ring's one | Uses the coefficient-ring one instead of a JavaScript number literal |
| Type interface | Dynamic result parent | `monic` and its existing `_monic` alias expose an Integer-coefficient overload returning an integer or rational polynomial. Their `RingElement` intersections bridge the pre-existing polymorphic-this constraint; general categorical base extension remains unimplemented |

**Rationale:** Preserve Sage's numeric and parent behavior within the coefficient rings
implemented by this port, while retaining the existing generic coefficient protocol.

**Trade-offs:** The generic TypeScript signature cannot infer every runtime parent conversion
of arbitrary user-defined coefficient classes. This repair covers ZZ, QQ and the implemented
finite parents; it is not Sage's complete base-extension machinery.

**Behavioral impact:** Integer normalization no longer floors fractional coefficients or
returns an integer parent for rational results. Zero and nonunit errors and object identity
are compared separately. Existing tests that accepted a zero GF(2) monic result were corrected
against the Sage oracle, rather than relaxing the new comparisons.


The specialized `GF2` coefficient class now preserves the same division/inversion error as
Sage's IntegerMod backend. Its `pow(number | bigint)` delegates to the existing IntegerMod
implementation for exponent conversion and the native/GMP error cutoff. This does not
complete the specialized class's broader constructor, mixed-parent or generator audit.

## Polynomial Index Conversion and Storage

**SageMath versus the port:** `getCoeff` (the named TypeScript counterpart to scalar
`f[n]`), `shift`, `truncate` and `reverse` now preserve the distinct generic, FLINT,
GF2X and NTL conversion boundaries for the implemented coefficient parents. Generic
negative truncation uses Python slicing; QQ reversal converts `degree + 1`; native
reversal retains unsigned increment wraparound. Finite scalar arithmetic can change
reversal's padding branch. Backend-specific zero shortcuts and object identity are
compared. The implementation still stores coefficients in JavaScript arrays rather
than the native mutable FLINT/NTL objects; array copies implement these linear storage
operations. Native polynomial-template shift's construction and division by a monomial
are represented by the corresponding coefficient slice/padding operation.

**Rationale:** Preserve observable conversion, coefficient and parent behavior while
using the port's established dense storage. Signed/unsigned native bounds are checked
with BigInt before converting an array index to number.

**Trade-offs:** JavaScript arrays have smaller allocation limits than native address
spaces. Resource-exhausting valid lengths are not claimed equivalent and are not run
as comparative fixtures. The retained `shift(number | bigint): Polynomial<C>` overload
assumes finite numeric inputs; its `unknown` overload includes null, but TypeScript's
number type cannot exclude NaN. A nonzero generic QQ/ZZ polynomial shifted by NaN
returns null as Sage returns None. Unknown third-party scalar protocols and arbitrary
coefficient categories remain outside this repair.

**Behavioral impact:** Exact indexing, fractional truncation, errors, output parents
and identity match the executed fixtures. This does not establish full categorical
coercion, allocation-limit or native storage/dependency implementation coverage.

**Oracle compatibility:** Installed Sage 10.3 says "non-negative" in reversal errors;
the bundled source says "nonnegative". The comparison adapter normalizes only this
specific message prefix. Sage 10.3's compiled generic sequence-repetition error cleanup
can decrement cached coefficient references and eventually crash the test process.
For the affected failing positive shifts and generic reversal padding, the adapter
executes the same scalar operation on an empty list, preserving the rejection without
including a cached coefficient. Rational, extension-element and oversized-integer
repetition errors were checked in isolated Sage processes against the original calls.
Successful operations still call the production Sage polynomial methods.

## Polynomial Common Coefficient Parents and Representation

**SageMath versus the port:** Sage's coercion model selects a common polynomial parent
before arithmetic. The repaired `add`, `sub`, `mul` and `eq`
perform that selection for wrapped ZZ, QQ, implemented modular/prime/extension parents,
and nested univariate polynomial bases. Finite parent selection delegates to the existing
finite-element coercion implementation, including quotient-ring common divisors. It no
longer treats equality of coefficient zeros as evidence of a canonical map. Coefficients
are converted and normalized before degree comparison or multiplication's zero shortcut.

**Rationale:** A result's parent must describe its actual coefficients, including when
reduction annihilates a leading coefficient. Reuse the existing finite coercion rules
rather than duplicate their maps in polynomial arithmetic.

**Trade-offs:** This is a restricted coefficient-parent selector, not the complete Sage
category/coercion framework. Finite extension embeddings beyond the existing finite
backend, sparse/multivariate constructions, unrelated symbolic/number-field parents,
and arbitrary user-defined parent maps remain outside this repair. Scalar arithmetic
protocols and the remaining mixed-parent polynomial methods are still under audit.
The existing generic multiplication storage/algorithm limitations remain unchanged.

**Behavioral impact:** Compatible operands now produce the correct coefficient parent;
incompatible characteristics and variable names fail at the polynomial boundary.
Polynomial equality compares the normalized images in the common ring. Same-coefficient
arithmetic retains its return type; the added mixed overload exposes `Polynomial<C | D>`.
That union describes the implemented coefficient classes, not arbitrary categorical
extensions. `eq` accepts unknown inputs without a compile-time polynomial cast.

String representation now follows `polynomial_element.pyx:_repr`: signs are joined using
subtraction, coefficients 1 and -1 are omitted on nonconstant terms, and parentheses
surround compound additive coefficients. Existing unit expectations for `x + -1`,
`(-1)*x` and `1*x` were corrected against the original Sage formatting.

## Polynomial Quotient and Remainder Backends

**SageMath versus the port:** Native ZZ division uses FLINT's partial integer quotient,
not generic exact coefficient division. The new dense `_fmpz_poly_divrem` backend ports
FLINT's basecase (divisor length below 6), recursive divide-and-conquer (cutoff 16), and
unbalanced block divisions. Small-magnitude leading coefficients are skipped; other
quotients use floor division. Uncancelled high terms remain in the remainder. The exact
flag rejects nonintegral leading quotients while allowing a lower-degree remainder.
Integer GCD's exact divisibility verification now shares this backend.

Word-sized modular division delegates to `_nmod_poly_divrem`, reusing the existing
classical/Newton kernel and unsigned BigInt packing infrastructure shared with
modular GCD. Its public array adapter is in `nmod_poly/divrem.ts`; the implementation
kernel is shared with `nmod_poly/gcd.ts`. Binary polynomial division delegates to
`ntl-ts/GF2X.DivRem`, using linear string conversion to/from its packed representation.
Sage-facing parent coercion, divisor checks, zero-input identity and errors remain above
the dependency calls. Generic nested coefficient division returns a polynomial quotient
only when its coefficient remainder vanishes.

**Rationale:** Reproduce both the original dependency boundary and quotient semantics;
non-exact native integer division is a supported operation. Reuse verified fast division
and packed products rather than replace them with coefficient-by-coefficient trial loops.

**Trade-offs:** The array backends retain the existing BigInt storage adaptations and
native multiplication tuning limitations described under Dense FLINT Polynomial GCD
Kernels and ntl-ts GF2X Representation. The raw C mutable-object APIs remain separate
stubs. QQ, large-prime and extension-polynomial fast division backends are not completed
by this repair; their pre-existing generic arithmetic path remains. Arbitrary fraction
fields as coefficient rings are outside the nested exact-division branch.

**Behavioral impact:** Shared comparisons verify quotient/remainder coefficients, parents,
identity and errors across implemented coefficient parents. NTL's large-modular zero
error is represented by exported `NTLError extends RuntimeError`, matching the Sage
wrapper. Direct integer backend comparisons call the installed FLINT C symbols through
Sage's loaded extension, including the raw exactness flag. Their nested array results
are JSON-encoded on both sides because the generic transcript formatter otherwise
misclassifies a two-coefficient quotient as a factorization.

## Polynomial Scalar Equality and Coefficient Embedding

**SageMath versus the port:** Equality now uses implemented integer, rational and finite
scalar parents, and embeds a coefficient polynomial as a constant in its outer ring.
For nonintegral JavaScript numbers, comparisons over ZZ/QQ and nested polynomial bases
reproduce real-double coefficient conversion, including rounding, underflow and infinity.
Generic multiplication preserves an already-zero outer polynomial during scalar actions.

**Rationale:** Sage compares after canonical coercion; testing only the original degree
or requiring both inputs to be Polynomial incorrectly rejects equal images.

**Trade-offs:** This bounded adapter does not implement the full RDF polynomial parent
or arbitrary Python rich-comparison protocols. As elsewhere in this port, integral-valued
JavaScript numbers represent Sage integers; JavaScript cannot retain a distinction between
Python integer 1 and float 1.0. Scalar arithmetic outside polynomial coefficient embedding
and arbitrary number-field/category maps remain under audit.

**Behavioral impact:** New shared cases verify scalar equality, coefficient embedding,
result parents and object identity. Binary64 conversion can annihilate higher coefficients
or overflow a constant; equality follows that converted image. Incompatible scalar parents
and unsupported objects compare false.

## Polynomial GCD Backend Coercion and Identity

**SageMath versus the port:** GCD coerces supported parents before backend dispatch.
ZZ/QQ and word-modular GCD delegate to the FLINT port; binary GCD delegates to NTL GF2X.
The new `_fmpq_poly_gcd` takes integer numerator arrays (input denominators cancel),
removes content, invokes the integer backend, and returns numerator coefficients and a
positive common denominator. Finite template zero shortcuts retain the original operand
without making it monic; equal nonzero inputs use `monic`. Other native calls allocate.
Modular Euclid stops at a nonzero constant without attempting to invert that constant.

**Rationale:** Match the original dependency boundary and backend-specific observable
behavior. TypeScript ring constructors are not globally interned; operand-returning GCD
shortcuts preserve the original operand when its parent and the canonical result parent
have coercion maps in both directions.

**Trade-offs:** Large-prime/extension fast GCD still uses the pre-existing generic
arithmetic path after native shortcut handling; full NTL ZZ_pX/ZZ_pEX kernels remain open.
Nested polynomial GCD still lacks Sage's Singular flattening backend. The existing BigInt
storage and tuning deviations in FLINT and NTL remain. The modular array kernel reports
failed inversions with RangeError; the Sage boundary translates native calculation failure.

**Behavioral impact:** Shared comparisons check coefficients, parent, allocation identity
and errors. The mixed GCD overload returns `Polynomial<C | D>`. The bundled Sage template
has an equal-input shortcut and guarded composite GCD which installed Sage 10.3 lacks.
The Python oracle applies that exact shortcut and invokes the original native GCD inside
a compiled Cython signal guard, translating its RuntimeError as the bundled linkage does.
This prevents old FLINT aborts from terminating the oracle; it does not replace native
GCD with a test-local polynomial algorithm.

## Polynomial Extended GCD Finite Backends

**SageMath versus the port:** Word-modular XGCD now delegates to `_nmod_poly_xgcd` in
FLINT. The kernel follows the native bit-size cutoff table, classical extended Euclid,
recursive half-GCD and final cofactor division, sharing normalized arithmetic, packed
multiplication and Newton division with GCD. Its HGCD matrix maps inputs to remainders
(the inverse of the source matrix), consistently with the existing GCD port. Binary
XGCD delegates to NTL GF2X. Supported mixed parents normalize before native dispatch.

**Rationale:** Reproduce native coefficients, algorithm complexity and observable zero
conventions. Finite template zero shortcuts preserve the original nonmonic operand;
large-prime NTL XGCD of two zeros returns a unit first cofactor. QQ returns three separately
allocated zeros. Cross-parent operand shortcuts account for non-interned TypeScript parents.

**Trade-offs:** Large-prime/extension XGCD retain their pre-existing generic arithmetic
kernels; fast NTL delegation remains open. QQ and ZZ delegation was completed in 2.11.0
(see Polynomial Integer and Rational Extended GCD).
Existing BigInt packing/tuning differences remain. The shared modular kernels expose an
internal arithmetic object for reuse; it is not exported from the package root.

**Behavioral impact:** Shared tests compare all three coefficients and parents, original
operand identity, and aliasing among the three outputs. Direct FLINT Cython comparisons
bypass Sage's zero shortcuts and exercise cutoff, exact-divisor, constant-remainder,
unbalanced and long-common-factor branches. Native inversion failures map to a comparable
error category at the array boundary; the Sage wrapper preserves the bundled ValueError
message. The older installed Sage's composite-XGCD abort uses the same narrow signal guard
as GCD, translating errors according to the bundled linkage.

## FLINT Resultant Kernels

**SageMath versus the port:** The new dense-array integer resultant follows FLINT's
length-144/bit-cost dispatch between a subresultant pseudo-remainder sequence and modular
reconstruction. The modular path removes content, uses the source Euclidean norm bound,
skips word primes dividing leading coefficients and reconstructs a centered result with
a balanced CRT tree. Modular resultants use small inverse-free formulas or the existing
Euclidean/recursive half-GCD arithmetic. HGCD now optionally tracks resultant factors,
including deferred degree drops and parity offsets through truncated recursive inputs.

**Rationale:** These native dependency kernels are required for integer extended GCD and
for replacing the Sage port's generic determinant path with its original delegation.

**Trade-offs:** The norm-bound bit length is computed exactly by comparing the squared
norm product with a power-of-two boundary, avoiding materialization of its integer square
root. Reconstruction uses balanced BigInt CRT nodes instead of FLINT's mutable comb
storage. Existing BigInt packing/tuning deviations remain; Montgomery/assembly dispatch
is not reproduced. Sage-facing ZZ/QQ/word-prime resultant integration uses these kernels as of 2.13.0;
integer XGCD has used them since 2.11.0.

**Behavioral impact:** `_fmpz_poly_resultant(a,b)` and `_nmod_poly_resultant(a,b,n)` normalize
dense constant-first arrays without mutation. The modular kernel retains native same-array
alias behavior: a shared nonzero constant array returns zero, while two distinct constant
arrays return one. Zero inputs return zero. Direct Cython/native FLINT comparisons preserve
this pointer distinction and catch composite inverse failures inside a signal guard.


## Polynomial Integer and Rational Extended GCD

**SageMath versus the port:** ZZ and QQ XGCD now delegate to FLINT-port kernels. Integer
XGCD uses the length-ordered resultant, skipped word primes, centered CRT, stabilization
checks and source coefficient bounds. Rational XGCD removes content/common factors and
rescales the native integer Bezout triple. Sage's integer wrapper preserves its scalar
constant/zero returns and clears all QQ output denominators when the resultant vanishes.

**Rationale:** Reproduce the native coefficients and algorithm, including its input-order
sign convention, rather than the stronger claim in Sage's docstring that the first output
always equals the original-order resultant.

**Trade-offs:** Immutable dense arrays replace native mutable storage; existing BigInt
packing/tuning differences remain. The standalone integer array kernel rejects two
nonzero constants unless the right one is a unit: native FLINT has undefined behavior
and can fail to terminate for the other cases. The Sage wrapper handles these inputs
through Integer.xgcd before delegation. Large-prime/extension NTL XGCD remains open.

**Behavioral impact:** Valid native inputs and Sage results are checked against actual
FLINT/Sage calls, including scalar types, parents, operand identity and output aliasing.
The array-only constant guard raises RangeError outside that native kernel's supported
input domain; it does not change Sage's defined constant behavior. The Integer-specific
TypeScript overload exposes scalar components previously hidden behind an unimplemented
branch; existing field overloads retain polynomial results.


## Polynomial Sylvester Matrices and Explicit Scalar Construction

**SageMath versus the port:** Sylvester matrices now use canonical common parents before
reading degrees, accept scalar operands, reject zero polynomials and reproduce optional
variable conversion. Sage indexes the left polynomial's variable tuple after coercion;
a constant left input can therefore raise IndexError before the zero-input check. Explicit
foreign finite-element construction follows `FiniteRingElement.polynomial()` when there
is no canonical coefficient embedding. Specialized GF2 coefficient wrappers are normalized
through the prime-field constructor. QQ singleton lists retain element-constructor
conversion rather than re-entering the finite-element parent hook.

**Rationale:** Preserve source error ordering and distinguish canonical coefficient
embedding from explicit conversion through a polynomial representative. Direct scalar
constructor comparisons supplement the Sylvester comparisons so ignored variable values
cannot hide wrong polynomial coefficients.

**Trade-offs:** The existing TypeScript API returns coefficient arrays rather than a Sage
Matrix object, so matrix methods, mutability and matrix-parent identity are unavailable.
The default small-field Sage matrix implementations require Givaro's `_cache`; the PARI
finite-field parents used by this port lack it. Their observable AttributeError is retained
for odd extension orders below 256 and nonempty characteristic-two matrices of order at
most 65536. This models the full reference installation with its optional matrix backends.
The internal canonical-parent helper is shared with the constructor but is not exported
through package entry points.

**Behavioral impact:** Shared original-Sage tests check entries, dimensions, scalar types,
coefficient domains and errors. The existing resultant/discriminant determinant paths use
the extracted entry builder directly, preserving their behavior while public Sylvester
validation is corrected. Their remaining FLINT/NTL/PARI algorithm delegation is still under
audit; this change does not claim that determinant fallback matches all native backends.


## Polynomial Resultant Delegation and Real Double Coefficients

**SageMath versus the port:** ZZ, QQ and word-prime resultants delegate to the corresponding
FLINT ports. The rational kernel removes content and denominators, probes a 60-bit modular
GCD and uses exact integer GCD if the probe is inconclusive before computing the integer
resultant. Composite word rings delegate to the existing modular matrix determinant.
Canonical scalar/mixed-parent conversion, extension-method coercion ordering, proof-keyword
errors and the word-prime same-object constant result are preserved.

**Rationale:** Replace the generic determinant where Sage uses FLINT, and preserve the
observable type and error behavior around delegation. Scalar floating inputs require a
real-double coefficient parent; the existing approximate MPFR parent represents a different
Sage type. The new RDF coefficient subset follows native IEEE arithmetic and formatting.

**Trade-offs:** Immutable dense arrays replace native storage. RDF currently implements
numeric coefficient construction and basic arithmetic only; other Sage methods and input
conversion hooks remain unported. Real resultants delegate to PARI's native double conversion and Sylvester elimination
as of 2.14.0. SciPy's finite-overflow fallback remains explicitly unimplemented, as do
resultants with nonconstant coefficients in nested RDF polynomial rings. Large-prime,
binary and extension resultants retain existing generic arithmetic pending native backend
integration. This batch does not mark those paths fully audited.

**Behavioral impact:** Supported paths have direct original-Sage/native comparisons.
Integer-to-RDF overflow yields signed infinity. Native PARI `rtodbl` loses negative zero,
flushes sufficiently small subnormals and rejects the highest finite binary64 exponent;
these quirks are reproduced, including when reached through polynomial resultants. IEEE
NaN payloads are normalized to a single `NaN` marker in the RDF test transcript, because
the payload is platform dependent; all finite and signed-zero bit patterns are compared.
The Buchmann binary64 wrappers now delegate to the same native conversion functions. Unsupported real-resultant
branches raise explicit not-implemented errors instead of claiming native equivalence.


## Native Single-Word Real Elimination

**SageMath versus the port:** General RDF resultants now delegate to PARI's inexact
Sylvester determinant path. The port uses the native 2-by-2 formula and column elimination,
selects the first pivot with maximal exponent, checks approximate zero against the original
input column and preserves early-return and sign conventions. Dependency-level `Polrev`
arrays retain inexact leading zeros; Sage removes those coefficients before delegation.
Cypari2 represents binary64 zero with accuracy exponent -53, distinct from raw `dbltor`.

**Rationale:** JavaScript double arithmetic and the older approximate QFB addition adapter
do not reproduce PARI's whole-word alignment and cancellation rules. The internal
`kernel/none/add.ts:addrr` implements the native single-word addition branch;
`kernel/gmp/mp.ts:divrr` now implements both single-word and multiword native division.
Multiplication reuses the existing exact-rounded kernel.

**Trade-offs:** The separate addition entry point requires 64-bit nonzero inputs.
Division is shared with QFB/Buchmann and accepts native whole-word precision. The internal
real-only `alglin1.ts:det` uses immutable mantissas and row arrays, with null for exact
integer padding zeros. Its empty determinant is represented as a real one at this typed
boundary. It is not the general PARI matrix API. Remaining real arithmetic gaps and
SciPy's finite-overflow fallback remain separate implementation work.

**Behavioral impact:** Direct native comparisons check exact rational mantissas and zero
exponents through primitive and composed operations. Resultants compare output binary64
bits, including pivot swaps, near cancellation, repeated operands, explicit zeros and
unbalanced degrees. Native arrays and Sage polynomials have separate comparisons so
normalization in one layer cannot hide a discrepancy in the other.


## Polynomial Derivative Protocol and Native Coefficients

**SageMath versus the port:** Polynomial derivative arguments now use the source
`derivative_parse`/`multi_derivative` protocol: variable/count pairs, consecutive counts,
single already-expanded lists and zero-count identity. Added `_derivative`, `diff`,
`differentiate` and `gradient`; nested variables recurse through coefficient derivatives.
QQ delegates to FLINT's rational derivative and always creates a new polynomial. Generic
zero differentiation keeps the original object. Generic coefficients use one scalar
multiplication, matching Sage's degree-times-coefficient expression and RDF rounding.

**Rationale:** The old method ignored every argument, omitted aliases and always created
a new zero. Its repeated-addition coefficient scaling also differed from native binary64
multiplication and added a logarithmic factor to each coefficient operation.

**Trade-offs:** The TypeScript utility accepts array sequences; general Python iterable
protocols and symbolic expression variables are not ported. Integer-like JS numbers are
counts, following the existing number/Integer mapping. Alias getters expose the same
method function. FLINT uses immutable normalized arrays rather than mutable polynomial
storage; rational array denominators must be positive. These representation choices do
not alter compared values in the supported domain.

**Behavioral impact:** Shared tests check arguments, alias identity, output coefficients,
parents, zero-count identity, zero-output allocation, nested variables and RDF rounding.
The original Sage 10.3 runner uses `derivative counts must be non-negative`; the bundled
Sage 10.9 source at `misc/derivative.pyx:163` uses `derivative counts must be nonnegative`.
The oracle translates only that exact older ValueError message to the bundled wording.
Native FLINT comparisons call QQ._derivative, which invokes fmpq_poly_derivative and then
_fmpz_poly_derivative; denominator-one fixtures exercise the integer kernel directly
through that original call chain.

## Polynomial Pseudo-Division and Fraction Parents

**SageMath versus the port:** `pseudo_quo_rem` now follows the original scalar zero test,
coefficient-base membership, original-operand powers and ordered scalar arithmetic. It
retains negative final exponents for degree gaps and can return a fraction-field quotient.
The new mirrored fraction-field modules currently implement constant-denominator data
and arithmetic. Sage also supports nonconstant denominators, full fraction coercion,
subclass dispatch and native FpT arithmetic; these remain open fidelity work. Public
polynomial `pow` still lacks the complete native/negative-power protocol. Pseudo-division
uses the source's constant negative-power rules internally; nonconstant negative
coefficient powers can still reach an explicit unimplemented fraction denominator.

**Rationale:** The previous implementation replaced coefficient membership with degree
zero, required polynomial operands, lost scalar parents, ignored negative powers and used
linear repeated multiplication. Its zero quotient also always had a polynomial parent.
The fraction type makes the represented parent and numerator/denominator inspectable;
it is not a complete port of Sage's fraction-field hierarchy.

**Trade-offs:** Supported comparisons cover exact integer, rational, modular, binary,
large-prime, extension and constant nested polynomial coefficients. Native polynomial
multiplication/powering dispatch, nonconstant fraction denominators, arbitrary coefficient
hooks and general real powers need further dependency work. Unknown custom coefficient
adapters retain their supplied ring context and use generic binary exponentiation when
no power method exists. Nonconstant-denominator construction raises an explicit
`NotImplementedError` instead of returning a value with a false polynomial parent.

**Behavioral impact:** Shared original tests check coefficient values, result parents,
operand/result identities, mixed variables, scalar attributes, zero/error ordering and
negative-power inversion errors. Direct fraction comparisons verify numerator and
denominator normalization, arithmetic, field/ring access and same-parent conversion.
Zero-ring scalar equality no longer obtains a quotient-parent map by first wrapping the
scalar as a polynomial; polynomial-to-polynomial comparisons retain their own coercion.

## Native Modular Polynomial Products and Fraction Fields

**SageMath behavior:** `fraction_field.py`, `fraction_field_element.pyx` and
`fraction_field_FpT.pyx` select generic, univariate-field and small-prime native fraction
classes. Generic arithmetic cancels cross factors and preserves unreduced intermediate
representations. FpT arithmetic delegates to FLINT modular polynomial products, powers,
GCD and division; its equality compares stored numerator and denominator. Generic
numerator/denominator accessors retain identity; FpT accessors allocate polynomial copies.

**Port:** The univariate ZZ, QQ and finite-field subset now follows those class and
arithmetic protocols, including nonconstant denominators, reduction flags, negative
powers and cached parent constants. `PolynomialRing.fraction_field()` selects the class;
`FractionField(R)` currently accepts polynomial rings. The explicit native product
kernels use constant-first normalized immutable BigInt arrays with a positive unsigned
64-bit modulus; powers take an unsigned 64-bit exponent. Dispatch follows portable FLINT
classical/KS/KS2/KS4 multiplication. Packed BigInts replace native limb arrays; KS4 keeps
the source's overlapping-coefficient recovery. Native architecture-specific `fft_small`
and mutable C allocation/alias machinery are not reproduced.

**Rationale:** Preserve source algorithms and dependency boundaries in the available
univariate polynomial representation. BigInt supplies portable arbitrary-precision packed
arithmetic. Generic coefficient GCD/division still uses its existing audited backend subset;
unsupported nested coefficient operations and complete coercion maps remain open.

**Trade-offs:** No native FFT acceleration or C memory-layout compatibility; performance
can differ. This does not establish full fraction-field coverage: multivariate bases,
general scalar-ring factory inputs, arbitrary coercions, factorization/support and native
iterators remain outside this implemented subset. Unreduced native fractions deliberately
retain the source's structural equality, including `0/x != 0` until normalized.

**Behavioral impact:** Compared valid arithmetic, parent data, representations and
identities match the original. Original FpT negative powers of a zero numerator can
segfault while reading an empty denominator's leading coefficient. The port raises the
same empty-message `ZeroDivisionError` as the original's checked inverse; the comparative
oracle uses that original inverse only for this crashing input. Explicit native element
construction with denominator zero is rejected safely rather than entering the original's
unsafe normalization path. These guards do not alter valid-denominator arithmetic.

Parent-class comparisons remove Sage's runtime `_with_category` suffix to compare the
underlying mirrored class; TypeScript does not synthesize Python category subclasses.
The installed Sage 10.3 factory's trailing period in `R must be an integral domain.` is
normalized to the bundled source's message without that period. Fraction constructor
exception handling explicitly includes the port's separate `ZeroDivisionError` class where
Python catches its `ArithmeticError` base.
The parent oracle also enforces the bundled FpT constructor's explicit
`Polynomial_zmod_flint` class check, which installed Sage 10.3 lacks. It applies only
after the native characteristic-range check and reproduces the bundled source's
`TypeError('unsupported polynomial ring')`; valid native arithmetic still runs in Sage.


### Fraction conversion follow-up (2.18.0)

**SageMath behavior:** Canonical fraction embeddings and partial polynomial/scalar sections
have different acceptance rules, messages, identity and normalization effects. A generic
fraction needs a unit denominator even when a larger coefficient field could represent
its reciprocal. Native FpT sections can normalize the original element, including before
raising an error. Fraction construction retries the source's ordered numerator/denominator
resolution and clears rational polynomial coefficient denominators when mapping to ZZ.
Polynomial rings over non-domains do not expose `fraction_field` in their Sage categories.

**Port:** These paths now follow the source for the compared univariate ZZ/QQ/finite-field
parents and nested polynomial/fraction coefficient embeddings. The shared TypeScript
PolynomialRing class necessarily has the method on its prototype; calling it on the
compared non-domain parents raises Sage's AttributeError, while the separate FractionField
factory keeps its integral-domain TypeError. Presence introspection therefore still differs.
Scalar constructors' direct fraction hooks and arbitrary parent maps remain under audit.

**Rationale and trade-offs:** The port retains its documented noninterned PolynomialRing
construction and shared class representation. The nested-coefficient comparative oracle
uses Sage's original `PolynomialRing_field` or `PolynomialRing_integral_domain` constructor
directly, preserving the explicitly supplied base parent. Sage's global factory can instead
reuse an equal cached parent with a different fraction subclass; using that factory would
make unreduced representation results depend on earlier unrelated cases. This adapter
compares original arithmetic and sections for the same specified parent, but does not claim
Sage global factory caching or categorical method-presence parity. Rational coefficient
storage uses the internal exact denominator-clearing adapter described in DESIGN.md.

**Behavioral impact:** The new comparisons include converted values, messages, result
parents, numerator/denominator data, retained identity and native mutation. Nested
pseudo-division no longer truncates fractional coefficient inverses to zero; constant
polynomial sections and polynomial/fraction equality use their canonical common parent.
Independent public fraction-module imports no longer hit an initialization cycle.


### Scalar conversion follow-up (2.19.0)

**SageMath behavior:** Generic fraction scalar conversion reduces the stored numerator and
denominator before converting them and taking the denominator's inverse in the target
ring. Native FpT uses separate partial sections; QQ rejects native fractions without
normalizing them. Scalar parent constructors can use a constant-polynomial section before
the element hook, including composed coefficient embeddings. Finite field constructors
must distinguish that section from evaluating a foreign polynomial at the field generator.

**Port:** The compared ZZ, QQ and finite-field constructor and hook paths now follow this
ordering. GF2 exposes its IntegerMod integer/rational hooks and accepts object constructor
inputs through the prime-field implementation. Generic fraction reduction allocates the
original native floor-division zero result instead of retaining `quo_rem`'s zero object.
Modular unit inversion follows the generic ring's precheck and ArithmeticError, while
direct modular `inv()` retains its ZeroDivisionError. Broader arbitrary scalar parents,
inexact rings and general categories remain outside these comparisons.

**Rationale and trade-offs:** TypeScript uses explicit hooks, structural polynomial data
and the existing RingElement `inv`/`isUnit` methods instead of Python category injection.
BigInt integer results have no methods or parent object; the internal adapter supplies
integer unit inversion and comparison frames name the requested scalar parent. This does
not add arbitrary Sage category methods to every TypeScript class.

**Behavioral impact:** New shared cases compare conversions, direct hooks, exact errors,
source normalization and stored-polynomial identities. The installed Sage 10.3
`ConstantPolynomialSection` says `not a constant polynomial`; the bundled source includes
the offending polynomial before `is not a constant polynomial`. The Python oracle adds
that operand only for this exact older message, using the actual offending polynomial.
The same correction applies to the earlier nested constant-polynomial section tests.
All arithmetic and conversion decisions still execute in the original implementation.

## Polynomial Integer Powers and Portable Native Products

**SageMath behavior:** Integer polynomial powers dispatch to FLINT's small, binomial,
multinomial or binary algorithms; rational powers preserve native common-denominator
storage. Word-modular, binary and extension polynomial templates delegate to FLINT/NTL,
with backend-specific zero, inverse, overflow and object-identity rules. NTL extension
polynomial multiplication uses Kronecker substitution into ZZ_pX. Native multiplication
selects classical, Karatsuba, FFT or Schoenhage–Strassen kernels according to platform,
coefficient sizes and degree.

**Port:** The compared integer-exponent power paths now delegate to those library ports,
including signed-word bounds, negative results in fraction fields and native identities.
Full `Polynomial.mul` now uses these same dependency boundaries; generic exact products
use SageMath's Karatsuba recursion and coefficient-ring thresholds. FLINT full rational
products preserve cross-content cancellation and alias-square raw denominators separately
from the fully canonicalized truncated-product adapter.
The portable FLINT shared product uses signed packed BigInt convolution with bounded
classical leaves and an even/odd Karatsuba fallback. NTL ZZ_pX uses modular classical leaves,
packed convolution and the native balanced/unbalanced Karatsuba split when a packed product
would exceed the conservative temporary-size budget. NTL GF2X uses carryless word leaves,
Karatsuba/block splitting and exact bit-interleaving squares. ZZ_pEX uses native Kronecker
strides followed by monic coefficient reduction.

**Rationale:** JavaScript BigInt has an engine-specific temporary-size limit that can reject
packed products even though every result coefficient fits. Recursive splitting avoids this
artificial failure without introducing an unbounded quadratic fallback. BigInt multiplication
provides the portable large-integer primitive; native FFT buffer code, assembly leaves and
thread scheduling have not been ported for these product boundaries.

**Trade-offs:** Native multiplication crossover tuning and full FFT/SS asymptotic performance
are not reproduced. Above the packed budget the fallback has Karatsuba complexity. Extension
coefficient reduction currently uses monic long division; native large-degree Newton remainder
selection remains open. Memory allocation, maximum feasible outputs and timing can differ.
The stateful NTL ZZ_pX class remains stubbed; the implemented array adapters are explicit.
Fractional polynomial exponents, modular powers and arbitrary coefficient-parent behavior
remain under audit; this entry does not claim full polynomial-power coverage.

**Behavioral impact:** Permanent comparisons check complete native coefficient vectors
(or SHA-256 of their exact signed hexadecimal serialization), errors, result parents and
identity. Rational adapter powers retain noncanonical input pairs such as 2/2. Coefficient
reconstruction uses the supplied ring's integer construction and field division.
The Python oracle serializes NTL residues through their exact string form: its `__int__`
method truncates to a C int. For the signed-long-minimum GF2X power case, the original short
Sage wrapper omits the signal guard and aborts; a Cython `sig_on` wrapper executes the actual
original operation and makes its native NTL error catchable. This changes crash handling
for the comparison only; it does not substitute an expected arithmetic result.

## Polynomial Roots and Truncated Series

**SageMath behavior:** ZZ/QQ rational polynomial powers call exact roots. Polynomial
roots use coefficient roots, valuation and characteristic extraction, Newton iteration
and exact verification. Truncated products, powers and inverses dispatch to FLINT;
extension inverse series calls NTL `ZZ_pEX::InvTrunc`. Method-specific Python/Cython
argument conversion and identity shortcuts affect errors and results. Factorization sorts
by degree, multiplicity and polynomial coefficients, starting at the highest coefficient.

**Port:** The compared paths follow this dispatch through library-owned dense adapters.
FLINT integer/rational/word-modular low products and powers preserve native zero rules;
integer and rational inverses retain their basecases and high-half Newton lifts. NTL
extension inverse uses NewtonInv, Kronecker products and ZZ_pX's recursive half-GCD for
coefficient inversion. The rational inverse's reversed-division basecase is bounded at
24 coefficients. The positive single coefficient-root bridge follows Johnston, prime-power
p-adic logarithms/exponentials and CRT; it does not expose the full finite-element root API.
P-adic parents now use the existing prime predicate instead of trial division.

**Rationale:** Dense arrays and existing BigInt/parent interfaces replace mutable native
buffers and Python category injection. Low products truncate inputs before using the
bounded packed/Karatsuba kernels where native code uses FFT/SS; this preserves exact
coefficients without artificial packed-temporary failures. Root selection uses the ring's
original default category algorithm: prime finite fields and integer-modulus rings can
choose different valid roots. Numeric factor sorting repairs the downstream root choice.

**Trade-offs:** Full FFT/SS tuning and asymptotic performance remain as documented under
Polynomial Integer Powers and Portable Native Products. NTL extension coefficient reduction
still uses monic long division. The internal positive coefficient-root bridge does not add
`all`, `extend`, negative roots or categorical methods to all finite classes. Sage's dynamic
category refinement and cached multiplicative-generator state are not generally reproduced;
comparisons cover the default operation sequences recorded in the shared cases. Mixed-parent
truncated products, arbitrary inexact coefficients and modular powers remain under audit.

**Behavioral impact:** Comparisons check results, identities, precise errors, native raw
coefficients and default single-root choices. The installed Sage 10.3 generic truncated-power
error says `non-negative`; the bundled source says `nonnegative`. The oracle adapts only that
exact message. Extension inverse precision conversion runs inside Sage's signal guard;
PARI-field precision objects can produce `SystemError: calling remove_from_pari_stack()
inside sig_on()`. The port preserves that observed wrapper error without touching a native
signal stack. This is separate from valid precision arithmetic, which delegates to NTL.


### Bounded native division and modular products

The portable-product limitation above also applies to integer divide-and-conquer division
and word-modular KS/KS2/KS4 multiplication. Sage's native libraries can allocate packed
limb buffers larger than JavaScript permits for one BigInt. The port now shares the bounded
FLINT product kernel at these boundaries and reduces modular results afterward. This fixes
valid divisions and products that previously threw a temporary-size error. The trade-off is
Karatsuba performance above the packing budget instead of native FFT/SS multiplication;
exact quotient, remainder and product coefficients are preserved. Permanent comparisons
include integer coefficients above 2^16384 and modular products through degree 65,534, plus
integer coefficients above 2^32768 and both aliased and distinct modular operands.


## Truncated Multiplication Parents and Real Coefficients

**SageMath behavior:** `multiplication_trunc` performs canonical coercion before strict
precision conversion. Python real scalars can promote QQ polynomials to RDF; strings,
lists and incompatible modular scalars do not gain implicit coefficient conversions.
Nested polynomial parents promote the coefficient field when their variable towers align.
The internal typed method accepts `None` in the tested wrappers, retains native precision
errors, and produces zero for positive precision. Generic real products use ordered
schoolbook sums, term actions and a distinct aliased-square loop.

**Port:** The compared scalar and nested-parent paths now use those canonical decisions.
The TypeScript internal signature explicitly includes `null`. Real and nested-real products
follow the generic source's operation order instead of initializing every coefficient sum
with zero. Scientific exponent signs and infinity tokens do not cause extra coefficient
parentheses. This corrects both output parents and IEEE signed-zero coefficients.

**Rationale and trade-offs:** Explicit parent/type adapters replace Python's coercion model
and Cython's generated type names. Recorded real fixtures encode input coefficient bits;
results compare all finite/infinite coefficient bits and normalize NaN payloads to `nan`.
NaN payloads are platform details and are not claimed equivalent. Generic RDF methods retain
Sage's quadratic arithmetic; native exact polynomial multiplication and broader variable
pushouts remain under audit. The `null` behavior is limited to the compared internal wrapper
contract, without relying on native object-layout reinterpretation.

**Behavioral impact:** Shared comparisons include rejected inputs, error precedence,
nullable operands, common-parent identities, nested QQ/RDF promotions, signed zeros,
subnormal coefficients, infinity and aliased squares. Existing exact-polynomial kernels
retain their own dispatch; this entry does not claim a complete public multiplication audit.

## Polynomial Modular Powers

**SageMath behavior:** Three-argument polynomial powers have backend-specific validation,
coercion, overflow and reciprocal rules. ZZ/QQ reject a supplied modulus. Word-modular
powers use FLINT binary powering and reversed-modulus preinverses for large exponents;
extension powers use NTL. Binary polynomials fall back to the generic implementation
outside signed-word bounds. Generic constant and negative powers ignore the modulus.
A negative template exponent computes a polynomial power before taking its fraction-field
reciprocal; this is not generally a modular inverse. C-long negation at its minimum value
passes a negative exponent into NTL and is preserved in the port.

**Port behavior:** `Polynomial.pow(n, modulus?)` now follows those wrapper branches and
uses new FLINT modular-power/inverse-series and NTL extension-power/XGCD adapters.
The modular half-GCD/Newton helper also delegates to bounded FLINT products, avoiding
oversized BigInt temporaries for valid 4,096/8,192-coefficient inputs.
The NTL extension kernel uses classical/Newton division, reciprocal modular reduction,
sliding-window powering and half-GCD. Immutable arrays replace mutable native objects;
existing bounded multiplication and portable native tuning limitations still apply.
The low-level FLINT inverse-series adapter rejects nonpositive precision with RangeError;
native FLINT aborts instead. Low-level invalid native arguments do not reproduce C process
termination. Public NTL failures use NTLError.

RDF integer powers now delegate to a small `gsl-ts` port. It preserves GSL's integer-power
operation order, including inversion before negative powers, and Sage's large-exponent
log/exp path. Only the value-returning integer-power, log and exp subset is implemented.
JavaScript Math replaces native libm calls, so transcendental last-bit agreement on other
engines/platforms is not guaranteed. Integer powers are accepted through IntegerLike;
general real or complex exponent coercion remains outside this RDF method.

Generic real division now follows Sage's descending subtraction order and explicit
remainder truncation, including NaN/infinite coefficients. Its unchecked quotient storage
preserves trailing zeros and their original display behavior. This does not establish
coverage of every other unchecked generic-polynomial constructor path.

**Rationale:** Restore deterministic values, parent selection, errors and dispatch while
retaining portable BigInt/binary64 representations and avoiding native process aborts.
**Trade-offs:** Native memory layouts, error handlers, FFT tuning and universal libm bit
identity are not reproduced. The three-function GSL package is not a general GSL port.
**Behavioral impact:** The compared public/native domain now agrees with the original;
invalid low-level precision and possible cross-platform transcendental rounding differ
as stated above. Previously ignored moduli, incorrect RDF negative powers, invalid real
unit-inverse rejection and unbounded NaN remainder degrees are fixed.

The permanent comparison oracle uses installed Sage 10.3 plus narrow wrapper guards from
bundled Sage 10.9.beta4: ZZ/QQ modulus rejection, extension exponent/canonical coercion,
large-extension native dispatch, and template zero/unit-modulus validation. Arithmetic
still executes the original FLINT/NTL/GSL code. Guarded Cython calls isolate native NTL
build/inverse errors. Sage 10.3's overflowing generic generator list repetition can corrupt
cached GF(2) values after repeated failures, so those original operations run once per
isolated Python process. Native FLINT zero-precision aborts and the port's RangeError are
compared using an explicit invalid-precision marker; other errors remain exact.

## Polynomial Evaluation and Composition

**SageMath:** Backend wrappers choose native FLINT/NTL scalar evaluation and composition,
then fall back to canonical coercion and generic polynomial evaluation. Degrees 4 through
50,000 use a cached sparse Horner graph with binary gap powers. Floating-point results depend
on that operation order. Zero and constant polynomials still validate incompatible parents.
Generators carry a distinguished identity, except native template classes whose `is_gen`
checks equality to the indeterminate.

**Port:** The compared paths now follow these decisions, including one-level argument-list
unpacking, nested coefficient evaluation, primitive float versus RDF result parents, native
QQ composition allocation, and large-prime first-power allocation. The public `compose`
adapter means `f(g)`; it is not a claim that Sage exports a generic `compose` method.
The new compiled evaluator preserves original coefficient truth tests and arithmetic order.
FLINT integer/rational evaluation uses the native Horner/binary-splitting cutoffs; composition
uses Horner, Taylor shift or divide-and-conquer as in the bundled source. Word-modular
composition retains its monomial-plus-constant specialization. NTL scalar evaluation uses
its original Horner loop for large composite-modulus and extension backends. Large prime
fields inherit generic compiled evaluation; their native power wrapper does not supply a
native scalar-evaluation override. All discovered discrepancies have shared comparative cases.

**Rationale:** Preserve original deterministic results while mapping native storage and
Python operator dispatch onto immutable coefficient arrays, BigInt, binary64 and explicit
coercion adapters. Minimal custom coefficient parents remain supported. The bundled integer
composition wrapper's tuple-type check is retained through its reachable generic fallback.

**Trade-offs:** Native fixed-limb/block Taylor-shift arithmetic is represented by BigInt
operations with the same triangular algorithm; threaded native work scheduling is not
reproduced. Existing bounded product and native tuning limitations still apply. Invalid
low-level moduli/denominators raise RangeError instead of C aborts or undefined behavior.
Globally interned polynomial parents, arbitrary symbolic/complex inputs, specialized matrix
evaluation classes and general multivariate argument dispatch remain open fidelity work. The mirrored
compiled class exposes its callable/eval/string interface, not native instruction-node
classes or Python pickling protocols. Its direct constant-list behavior intentionally
matches the original even though public evaluation bypasses compilation for constants.

**Behavioral impact:** Compared exact outputs, IEEE result bits, parents, identities and
errors agree. NaN payloads are normalized in comparisons. Portable storage, invalid native
arguments and unsupported input domains differ as described; full behavioral coverage is
not inferred from line coverage. Source evidence is in `polynomial_element.pyx`,
`polynomial_compiled.pyx`, the native polynomial subclasses, FLINT's `fmpz_poly`,
`fmpq_poly`, `nmod_poly` and `gr_poly` evaluation/composition/Taylor routines, and
NTL's `ZZ_pX1.cpp` and `ZZ_pEX.cpp`. Remaining dispatch domains require those original
coercion/action protocols rather than ad-hoc scalar conversions.

**Oracle version adapter:** Installed Sage 10.3's FLINT canonicalizes some internal rational
evaluation results, while the bundled FLINT preserves unreduced denominators. The Cython
oracle template compiles the exact bundled Horner, divide-and-conquer and wrapper C function
bodies against installed native FLINT primitives. Only function-name prefixes and the
installed header's `slong` alias are adapted. This independently verifies unreduced pairs;
it does not calculate expected answers with the TypeScript algorithm. Compilation failures
are cached and surfaced as test failures. Whole PARI-backed polynomial calls are not placed
inside Cython signal guards.

Installed Sage 10.3 also lacks bundled `Polynomial_dense_mod_p.__pow__`: it inherits generic
powering and returns an unflagged input unchanged at exponent one. The oracle's narrow
`bundled_plain_power` adapter calls the original NTL power and reconstructs the polynomial
for that case, as the bundled wrapper does. It does not patch the expected identity bit.

The installed large-composite NTL scalar wrapper also constructs its answer in the
polynomial parent; the bundled wrapper constructs it in the coefficient base. The oracle
executes native evaluation then applies that original base constructor only on its canonical
scalar branch. NTL's large-composite `is_gen` returns an integer, represented as `0n`/`1n`;
other implemented polynomial backends return booleans. Public integer evaluation accepts
arbitrary-size points and never applies the word-sized polynomial index conversion.

Explicit generator-flag construction also mirrors the backend constructors: QQ and binary
backends initialize the indeterminate, while other list-based constructors keep the supplied
coefficients. The integer backend's native representation ignores the distinguished flag.
These cases are independently compared through the original element-class constructors.

## Polynomial Matrix Evaluation Actions

**SageMath:** Generic polynomial evaluation resolves canonical matrix/scalar parents and keeps
fast scalar actions. Matrix-space pushout can change the implementation when the coefficient
base changes. Integer-zero canonical coercion admits rectangular matrices, but nonzero scalar
addition and incompatible matrix products raise. Integer-to-real matrix conversion calls
`fmpz_get_d`, which truncates rather than using scalar integer rounding.

**Port:** Generic, integer and modular matrices now follow those evaluation actions and
validation paths. Sparse compiled powers retain their original base until promotion is
needed. Matrix parent metadata records default versus explicitly generic implementations.
The new `integer_to_real_double_dense` delegates to the FLINT bigint conversion adapter.
Default QQ, RDF and binary extension matrix values use existing generic matrix storage;
existing native integer/modular classes remain concrete where available.

**Rationale:** Restore original evaluation values, rounding, parents and exceptions while
reusing the port's available matrix arithmetic and portable storage.
**Trade-offs:** This does not port NumPy/BLAS, M4RIE or every native matrix backend. Existing
matrix arithmetic complexity/tuning limits remain. Platform-specific BLAS last-bit behavior
outside the compared domain is not guaranteed. Additional polynomial coefficient domains,
arbitrary symbolic/complex actions and global matrix parent interning remain open.
**Behavioral impact:** Compared exact values, IEEE bits, result coefficient parents and
errors agree. Concrete native object layouts differ. The original small odd-extension
MeatAxe constructor's missing PARI `_cache` error is preserved because evaluation itself
raises; valid binary extension results are represented by the existing generic matrix.

**Oracle adapters:** Installed Sage 10.3 lacks bundled `MatrixSpace.is_exact`; the oracle
executes that exact bundled method, extracted from its AST, for the duration of evaluation.
The bundled method delegates to the base ring instead of inheriting unconditional exactness.
When original binary-extension evaluation returns a native M4RIE matrix whose public getter
fails on a PARI coefficient parent, a separate Cython reader serializes the actual stored
words and constructs coefficients from them. It does not replace matrix evaluation or hide
its exceptions. A separate module avoids conflicting FLINT/M4RIE header declarations.
The FLINT conversion oracle compiles the exact bundled `fmpz_get_d` body with its 53-bit
constants and installed header typedef aliases. Native matrix conversion is also compared
directly; no expected values are computed with the TypeScript implementation.

## Specialized Binary Matrix Multiplication

**SageMath:** Binary matrices delegate products to M4RI. The default uses Strassen-Winograd
with a distinct Bodrato squaring schedule, M4RM Gray-code tables and native classical
fallbacks. Explicit algorithm parameters are signed C ints. Some small cutoff geometries
produce invalid zero-width native windows, causing an abort or segmentation fault that
Sage's signal guard turns into RuntimeError.

**Port:** Binary matrix products now delegate to `m4ri-ts`, preserving the packed algorithms,
squaring identity, word-aligned splitting and rectangular tails. Scalar bit-array cubic
multiplication is removed from the public/default path. `_multiply_m4rm` and
`_multiply_strassen` accept IntegerLike parameters; their signed C-int overflow, empty
shortcuts and dimension errors are compared. The native cutoff faults are represented by
RuntimeError without terminating JavaScript. Polynomial evaluation supports this matrix
class, and its overloads retain matrix results for mixed and argument-list inputs.

**Rationale:** Restore delegation and algorithm complexity while keeping browser-compatible
bit storage and the existing public matrix interface.
**Trade-offs:** Native output pointers, writable shared windows, SSE/OpenMP scheduling,
allocation caches and hardware cache detection are replaced by owned immutable packed rows
and fixed tuning. Packed adapters omit native destination buffers and `mzd_make_table` uses
column offset zero. The low-level adapters validate shape/count/window/table preconditions
with RangeError instead of invoking undefined C behavior. They define zero-dimensional
products through the same zero-result convention as Sage's public wrappers. Zero-column
Gray-code tables with positive k are outside the native comparison domain. The existing
classical public guard still rejects incompatible nonempty dimensions safely; native C has
no corresponding safe contract, so comparisons use an explicit invalid-native-input marker.
**Behavioral impact:** Compared exact matrix values, errors and allocation rules agree.
Native storage aliasing, process termination and timing differ as stated. Ordinary native
matrix arithmetic, other binary methods and broader coefficient/action domains remain
separate audit work; these tests do not establish whole-repository behavioral coverage.

**Oracle:** Most valid kernels run installed M4RI through Cython. The comparator reads and
writes its actual native row words; zero-width rows are never dereferenced. Small Strassen
boundaries run in separate processes linked to M4RI 20251207, selected by the bundled Sage
package metadata and verified against SHA-256
`7b195a4d88fa827b9ec6d087c3cac739ab6e5100da05faaed3a1d2c20ca3a930`.
The build cache downloads that archive on first use, compiles the unmodified library and
links a small fixture/serialization program. Both the installed 2020 library and the pinned
2025 source exhibit the compared zero-window faults. Isolated native signals become the
same RuntimeError names/messages seen under Sage's signal guard. The public oracle uses
this native delegate plus original parent reconstruction for bounded small-cutoff cases;
it does not calculate expected products using TypeScript or a replacement algorithm.

## Binary Matrix Indices and Slices

**SageMath:** Matrix subscription checks integer-index protocols before assigning to a C
index. Reads use a signed C int and rational `__int__` truncation; assignment, row access
and submatrix arguments use signed `Py_ssize_t` and exact rational `__index__` conversion.
Entry/row indices may be negative. Any negative submatrix size selects the remaining
extent. Nonempty rows and slices delegate extraction to M4RI.

**Port:** `get`, `set`, `row`, `submatrix`, `swap_rows` and `swap_columns` preserve these values, conversion/bounds error
classes, messages and validation order for supported scalar indices. They accept bigint,
Integer, Rational and boolean, plus the existing integer-valued JavaScript number indices.
Slices delegate to `m4ri-ts/mzd_submatrix`; only the intersecting source words are packed,
so extracting a small slice does not scan the entire input matrix.

**Rationale:** Restore original indexing and extraction while preserving the established
array-index interface and portable bit-array storage.
**Trade-offs:** Integer-valued JS numbers represent integer indices here, even though
Python float indices are rejected; non-integral numbers and nonfinite numbers are rejected.
Arbitrary Python index-protocol objects and native vector/cache identity are not exposed.
Returned row arrays are owned copies. Native destination buffers and writable aliasing are
omitted. The low-level M4RI adapter validates bounds with RangeError and uses Sage's empty
slice shortcut, avoiding invalid zero-width C row accesses. Extremely large positive slice
extents that would overflow Sage's internal C-int endpoint arithmetic are rejected against
the actual matrix dimensions; undefined native memory accesses are outside the comparison.
**Behavioral impact:** The compared entries, row/slice values, allocation identity, integer
conversion and exceptions agree. Representation, aliasing, unsupported protocols and
undefined native inputs differ as stated. Algorithm selection, cached matrix results,
density and other binary matrix methods remain in the active audit.

## Native Real Literal Conversion

**SageMath:** `RealLiteral` retains its original text and uses MPFR for conversion at the
requested precision, printing and exact rational extraction. The bundled version reparses
at 53 bits when converting a literal to a Python float.

**Port:** `RealLiteral`, `create_RealNumber` and precision-converted results delegate decimal
conversion to `mpfr-ts`. The supported domain is nearest-even rounding, precision 1–4096
bits, at most 4096 significant input digits, decimal scale from −4096 through 4096, and at
most 4096 output digits. Other bases, rounding modes and larger working domains raise
`NotImplementedError`. MPFR's process-global flags and allocated buffer interfaces are
not exposed. Invalid precision at the Sage layer now raises the original `ValueError`.

**Rationale:** Exact bounded BigInt base powers and quotient/remainder rounding provide the
original conversion results without a native MPFR runtime. This replaces MPFR's adaptive
interval/Ziv working-precision loop within an explicitly bounded domain.

**Trade-offs:** This is a conversion subset, not a full arbitrary-precision real arithmetic
port. Inherited `RealNumber` arithmetic and transcendentals still use the existing binary64
implementation. General field coercion and nondecimal literals remain outside this repair.

**Behavioral impact:** Literal precision, decimal formatting, exact rationals, signed zeros,
IEEE conversion and mutable low-level parsing match the original within the compared
conversion domain. Precision changes reparse the text; absolute value and negation preserve
native state. Outside the supported domain, explicit failures replace native capabilities.

**Oracle:** Native conversion tests execute MPFR 4.2.1, the release pinned by the bundled
Sage source (verified archive SHA256
`b9df93635b20e4089c29623b19420c4ac848a1b29df1cfd59f26cab0d2666aa0`). Installed Sage 10.3 lacks
the bundled `RealLiteral.__float__` override, so that observation executes the bundled
expression `float(x.numerical_approx(53))` explicitly. Other observations call the originals.

## Binary Matrix Density

**SageMath:** Default density is an exact rational, except empty matrices return Python
integer zero. Approximate density calls M4RI with resolution one and wraps its double in a
`RealLiteral`. Native sampling omits the final full word of wide matrices whose column count
is a multiple of 64; that behavior is retained.

**Port:** Default density now returns `Rational | bigint`; approximate density returns
`RealLiteral` and delegates to `mzd_density`. The dependency also exposes the original
resolution parameter, including its automatic sampling stride.

**Rationale:** Preserve both native sampling and the original exact/approximate type and
precision distinction instead of collapsing all results to JavaScript numbers.

**Trade-offs:** This corrects a public return type and requires a major version change.
For positive row count and zero columns, native approximate density dereferences a null
row and faults. The port raises `RuntimeError` instead of terminating the host process.

**Behavioral impact:** Exact fractions and approximate literal observations match compared
originals, including zero-row NaN and native word-boundary sampling. The safe exception for
the native fault and the MPFR conversion limits above remain explicit differences. Fault
comparisons run the original kernel in isolated child processes.

## Binary Elimination and Cached Forms

**SageMath:** Binary rank and elimination select M4RI algorithms, preserve algorithm-specific
unreduced forms, cache results, and return immutable cached echelon forms. Classical
elimination always reduces. Empty `echelonize` returns self; other successful calls return
None. Cached values bypass later algorithm validation. Lexical ordering uses unsafe native
swaps and retains its cache, including potentially stale cached forms.

**Port:** These algorithms, cache transitions, return values and mutability checks are now
preserved. `echelonize` returns `void | this`; its optional fourth options object carries the
original M4RI k keyword. Pivots remain copied number arrays rather than Python tuples.
The dependency returns owned packed matrices, ranks and transposition lists instead of
mutating native pointers. Triangular solvers also return owned results.

**Rationale:** Preserve the original method dispatch, algorithm complexity and stateful
observations while representing native buffers with portable BigInt rows and owned values.
The explicit classical loop remains only for the original classical algorithm.

**Trade-offs:** Native scheduling, vector instructions and cache discovery are replaced by
a portable 4 MiB tuning model. BigInt interval operations have different memory costs than
in-place word copies. Rank/elimination dependency adapters return immediately on zero dimensions,
matching the Sage wrapper instead of exposing undefined native empty-buffer accesses.
Invalid triangular adapter dimensions and invalid low-level table ranges raise `RangeError`
instead of invoking undefined native memory access or terminating the process.
Standalone Sage `ple`/`pluq` wrappers are repaired in the follow-up below; other binary
methods remain under audit. This entry does not claim the entire matrix module is repaired.

**Behavioral impact:** Compared ranks, entries, factorization permutations, triangular
solutions, validation order, cache identity and immutability match originals. The return
type and immutable cached forms correct observable public behavior in 5.0.0. M4RI k values
11–16 can create panels wider than a native word: the pinned release disables assertions
and its 64-bit shift behavior is reproduced, even when it produces a noncanonical result.
That behavior is tied to the compared release/build; assertion-enabled builds may abort.

**Oracle:** Algorithms execute the checksum-pinned bundled M4RI release. Installed older
M4RI faults on empty column swaps; these cases execute the pinned guarded kernel with the
original Sage cache invalidation. Installed Sage 10.3 lacks `doubly_lexical_ordering`, so its
bundled method is compiled unchanged apart from dedentation and a typed self parameter.
This preserves native swaps and cache effects rather than substituting a TS-derived oracle.

## Standalone Binary Factorizations

**SageMath:** Standalone `ple`/`pluq` convert their parameter to a C int before checking the
algorithm, copy the matrix and delegate to the selected standard, Russian or naive M4RI
kernel. Fractions and finite floats truncate; nonfinite values and nonnumeric arguments
raise the original conversion errors.

**Port:** These functions now delegate every algorithm and accept bigint/Integer, the
existing number inputs, Rational and boolean parameters. Numeric conversion, C bounds,
error messages, mutable copied results and P/Q transposition lists are preserved.
The Russian backend uses the pinned seven-table schedule (correcting the earlier eight-table
port), with native 64-bit reads for oversized panels. Explicit k is supported through 16;
other sizes raise `RangeError` for nonempty matrices, even where a particular native input
might finish.

**Rationale:** Preserve the original wrapper and backend algorithms without exposing native
pointers, process termination or unbounded exponential table allocations. The native header
specifies table k through 8; the tested larger domain extends that through 16.

**Trade-offs:** Arbitrary user-defined numeric conversion protocols are not implemented.
Guarded invalid table sizes do not reproduce native allocation failures or undefined memory
access. Build-dependent oversized-panel behavior remains tied to the pinned oracle.

**Behavioral impact:** Compared values, factors, permutations and errors match. Standard
PLE/PLUQ on a positive-row, zero-column matrix dereferences invalid native storage; the port
now raises `RuntimeError('Segmentation fault')`. Explicit Russian and naive empty kernels
retain their native zero-rank results. Other standalone parameter and matrix observations
are exact in the compared domain.

**Oracle:** Normal wrapper calls use Sage directly. Empty/high-table cases project the
bundled wrapper's copied result from the pinned native kernel. The PLE kernel oracle's old
zero-dimension shortcut has been removed: comparisons now execute and isolate actual native
faults. Elimination's separate Sage empty shortcut remains intact. High-parameter regressions
also verify the seven-table correction and its effect on factors and rank.


## Binary Matrix Permutation Image Lists

**SageMath:** `permute_rows` and `permute_columns` accept a permutation group element,
check mutability, then apply reversed nontrivial cycles in order. Fixed points outside
the matrix are ignored; a later cycle bounds error retains earlier completed mutations.

**Port:** The existing 1-based number array API represents the image list of an element
of `SymmetricGroup(perm.length)`. Validated lists use the same cycle order and mutations.
Invalid image lists fail before matrix mutation, as constructing the original group
element before calling the method does. Arbitrary permutation-group objects are absent.

**Rationale:** Retain the established array API while restoring the underlying group action.
**Trade-offs:** No parent/group metadata or user-defined cycle protocol is exposed.
**Behavioral impact:** Compared actions, degree boundaries, exceptions and cache effects
match under the explicit image-list adapter. Raw Python lists are not passed directly to
Sage matrix methods. Empty swaps execute the pinned native guard where installed older
M4RI would crash, as documented for binary cache comparisons.

## Native Binary Inverse Adapter

**SageMath/M4RI:** Sage checks square shape and rank before `mzd_inv_m4ri`. The native
kernel places identity at the next 64-bit boundary in a twice-word-aligned augmented
matrix, fully reduces it, and copies that window. It ignores k and does not reject singular
matrices when called directly.

**Port:** Sage follows the same checks and delegates. `m4ri-ts/mzd_inv_m4ri(A, k = 0)`
returns an owned packed matrix, retaining singular raw-kernel outputs and ignored k.
**Rationale:** Preserve the original dependency algorithm and window layout.
**Trade-offs:** Native destination buffers and aliasing are omitted; nonsquare native
inputs raise RangeError instead of invoking assertion-disabled undefined memory behavior.
**Behavioral impact:** Compared valid native inputs and Sage values, errors, copied
mutability and rank-cache effects agree. Unsupported destination/undefined domains differ.


## Binary Linear Solve and Kernel Adapters

**SageMath:** Nonsingular square systems use augmented elimination; rectangular and
singular systems delegate to M4RI with zero-padded copies. Consistency checking defaults
to true. Kernels use PLUQ unless the generic algorithm is requested, then transform the
basis as requested. A zero-row matrix's kernel is an immutable identity matrix.

**Port:** `solve_right(B, check = true)` follows that dispatch for binary matrix operands.
`right_kernel_matrix` supports the original algorithm, basis and proof keyword validation.
`m4ri-ts/mzd_solve_left` returns `{matrix, rhs, status}` for the mutated input pair and
native return code. `mzd_kernel_left_pluq` returns `{matrix, kernel}`, with a null kernel
for full column rank. These native adapters preserve the original padding-check offset,
unchecked inconsistent outputs and tested empty-buffer faults as RuntimeError exceptions.

**Rationale:** Restore the original dependency algorithms, outputs and state changes.
**Trade-offs:** Native destination aliasing is omitted; invalid raw solver dimensions
raise RangeError instead of aborting the process. Sage vector operands and cross-ring
coercion are not yet exposed by this specialized matrix-only API. Native empty triangular
buffers inside the solver retain their original fault behavior; Sage's empty RHS/column
shortcuts remain safe. Standalone triangular adapters retain their documented empty shortcut.
**Behavioral impact:** Compared matrix values, kernel bases, options, errors, mutability
and cache effects agree. Unsupported object protocols and unsafe native dimensions differ.

The installed Sage 10.3 has no binary `_solve_right_general` override. Its oracle therefore
compiles the bundled method, changing only indentation, typed self and a missing old .pxd
function declaration; public shape/rank dispatch follows bundled `matrix2.pyx`. The native
solve/kernel outputs are separately compared with the pinned M4RI release. Installed Sage
also lacks the newer kernel `linbox` option: its bundled GF(2) rejection is applied before
other kernel options. These are version adapters, not port-output substitutions.


## Binary Basic Operations and Column Cache

**SageMath:** Binary addition, transpose and augmentation delegate to M4RI. Column
extraction transposes first, caches immutable column vectors, and returns a shallow copy
of the outer list unless `copy=False`. Matrix mutations invalidate that cache. Shape
errors follow the original parent/operator and specialized method checks.

**Port:** These operations delegate to `m4ri-ts`. Columns are frozen number arrays with
an assignment trap preserving the original immutable-vector ValueError. The outer cache
list is shared for `columns(false)` and copied otherwise; deliberate edits to that shared
list retain the original cache effects. Mutating the matrix discards the cached list.
Native transpose/concatenation return owned packed matrices. Raw empty rectangular native
transpose calls retain their original abort as RuntimeError; Sage transpose skips them.

**Rationale:** Restore dependency algorithms, column aliasing, mutability and exact errors.
**Trade-offs:** Column arrays replace Sage vector objects, without vector-parent metadata
or arbitrary Python vector/list protocols. Native destination aliasing and SIMD intrinsics
are omitted. Transpose retains the recursive 512-bit block partition and six native
bit-swap masks; small padded word blocks replace the native unrolled packing variants.
Concatenation shifts a packed bit interval instead of individually appending native bits.
Invalid native adapter dimensions use RangeError where documented; native copy aborts
cannot terminate the host process.
**Behavioral impact:** Compared values, errors, cache identities and mutations agree.
Column-entry assignment now rejects mutation, requiring callers to copy a vector first.
Unsupported protocols, low-level aliasing and instruction-level tuning differ. Subdivision
metadata and extended formatting remain separate open work.


## Binary Matrix Formatting and Subdivisions

**SageMath:** Binary custom mappings delegate to `matrix0.str`. Keyword zero/one strings
override and mutate dictionary entries; callable mappings override those keywords. Width
is global and counts Unicode code points. Binary minus_one is ignored. Subdivision changes
clear caches before sorting and ZZ coercion; duplicate lines are retained. Copy, negation
and transpose preserve metadata, while augmentation follows its original boundary/empty
shortcuts and compatible-row rules.

**Port:** Binary `str` delegates to a matrix0 renderer adapter with a lazy entry snapshot.
Dictionary keys 0 and 1 represent the field elements; callbacks receive number bits.
`subdivide`, `subdivisions` and `get_subdivisions` expose the original metadata operations,
with IntegerLike, numeric, rational, boolean, string and null line adapters. Returned line
arrays contain bigints and are owned copies. Pass row/column lists as separate arguments,
or spread the pair returned by `subdivisions()`.

**Rationale:** Restore existing mapping/subdivide behavior without duplicating a whole
Sage matrix inheritance hierarchy or Python container runtime.
**Trade-offs:** Python tuple-as-single-argument subdivision syntax is not exposed because
JavaScript arrays already represent line lists. The renderer covers the existing four
positional formatting arguments and ASCII subdivision lines; Unicode bracket options,
label borders and CharacterArt objects are not exposed. Arbitrary Python mapping/comparison
protocols are absent. Original clamped row insertion and signed-size errors are preserved;
unsafe column positions use checked generic rendering;
undefined native fast-path byte writes outside matrix bounds are not reproduced.
**Behavioral impact:** Compared strings, mapping mutations, errors, metadata and cache
states agree in this domain. Unsupported rendering/container protocols and unsafe invalid
layout positions differ. Installed Sage 10.3 drops transpose metadata; its comparison
applies the explicit metadata step from the bundled transpose method before observation.
Its older binary formatter uses narrower C integers, so subdivision string observations
compile the bundled method with only dedentation and a typed self argument.


## Binary Constructor and Row Ownership

**SageMath:** MatrixSpace converts and validates dimensions before binary native allocation.
MatrixArgs accepts scalar, flat and nested inputs, validates sized iterators in traversal
order, and converts entries through GF(2). A nonzero scalar requires a square matrix before
field conversion. `row(i, True)` returns the immutable cached row; default row extraction
returns a new mutable vector. Matrix mutations invalidate cached vectors.

**Port:** Dimensions accept IntegerLike, number, Rational and boolean through a MatrixSpace
adapter. Scalar and flat/nested array entries and `set` values use the prime-field backend.
`matrix/args.ts` retains sequence classification and the original sized-iterator validation
order, including Rational's single-coefficient sequence protocol. Cached rows use the same
immutable number-array adapter as columns; copied rows remain mutable.

**Rationale:** Preserve these existing constructor/setter/row behaviors without recreating
Python's complete matrix parent and vector hierarchy.
**Trade-offs:** Native pointer allocation becomes owned JavaScript arrays; practical memory
limits differ. Dictionary/callable/matrix/NumPy input protocols are not exposed by this
constructor adapter. Row arrays do not carry Sage vector parents or arithmetic methods.
Integer-valued JavaScript numbers represent Python integers; fractional numbers represent
Python floats. Separate bigint and Integer inputs preserve exact large values.
**Behavioral impact:** Compared values, shape/coercion errors, identity, immutability and
cache state agree. The listed container/vector protocols and allocation limits remain
outside this domain. Installed Sage 10.3 predates bundled MatrixArgs string classification
(issue 34821); numeric-string fixture adapters apply the bundled scalar-shape branch and
convert valid flat strings before calling the old constructor. No invalid numeric strings
are covered by that version adapter. Default row copies remain independently mutable;
assigning to `row(i, true)` now raises the original immutable-vector error.


## Binary Seeded Randomization and Factories

**SageMath:** Binary randomization uses the current GMP state. Full-density uniform rows
consume two 32-bit draws per 64-bit word, including padding, reverse the bits for M4RI's
word layout and mask unused columns. Sparse uniform updates select a fixed truncated
number of column positions with replacement. Nonzero updates draw a two-part uniform
floating value for every entry. Explicit random-factory density selects nonzero updates;
omitted density selects uniform bits. Dimension conversion precedes identity creation.

**Port:** Binary matrices delegate random draws to the existing GMP randstate port, using
the same draws, bit placement, update counts and cache/error order. Convenience factories
bind Sage's base-ring argument to GF(2). The entry factory infers dimensions and validates
all row lengths before field coercion, using the original error messages.

**Rationale:** Preserve seeded values and generator advancement while retaining the public
bit-array and GF(2)-specific convenience APIs.
**Trade-offs:** Packed native writes become assignments into the owned bit arrays; raw
native memory layout is not exposed. Density conversion uses the Python/Sage adapter described below; arbitrary Python object
protocols remain outside the exposed numeric/text/byte-buffer API.
**Behavioral impact:** Seeded numeric-density outputs, subsequent random draws, caches,
mutability and compared factory coercion/error behavior agree. The broader density-input
protocols and native memory representation are not covered by this adapter.


## Binary Density Input Conversion

**SageMath:** Empty matrices return before converting density. Otherwise, Python `float`
conversion precedes the nonpositive-density shortcut and mutability/cache checks. Built-in
integers reject float overflow, while Sage integers/rationals permit infinity. Text conversion
accepts decimal underscores, Unicode decimal digits/whitespace and case-insensitive infinity
and NaN; byte buffers use ASCII conversion. Invalid conversion leaves matrix caches intact.

**Port:** `randomize` and the random factory accept number, bigint, Integer, Rational,
RealNumber/RealLiteral, boolean, prime/modular field elements, strings and Uint8Array buffers.
The internal `types/python_float` adapter uses native binary64 conversion with CPython's
accepted grammar and pinned Unicode 15.0.0 digit/printable tables. Bigint maps to Python int;
Integer preserves Sage's distinct overflow behavior. Direct null density raises TypeError;
the random factory retains Sage's separate None/default dispatch.

**Rationale:** Preserve source conversion and error precedence, with explicit equivalents
for TypeScript coefficient values and byte buffers.
**Trade-offs:** Custom Python `__float__`/`__index__` methods and arbitrary buffer objects
are not exposed. Uint8Array represents Python bytes. JavaScript owns string/array allocation;
its memory limits differ. Unicode data is pinned instead of inherited from the JS engine.
**Behavioral impact:** Compared density outputs, random-state consumption, cache effects,
conversion errors and binary64 text results agree, including signed zero and parsed NaN
signs in Bun. Native JavaScript operations may canonicalize NaN payloads after conversion;
randomization does not observe them. Pixel-data empty errors now match bundled Sage wording.
The Gamma/PNG oracle uses original compiled/libgd functions; the older installed PNG method's
capitalized empty error is adapted to the bundled source's lowercase message.

The text oracle uses CPython 3.12 with Unicode 15, matching the bundled Python dependency;
Sage 10.3's older Python cannot serve as the Unicode oracle. Sources:
[CPython float conversion](https://raw.githubusercontent.com/python/cpython/v3.12.5/Objects/floatobject.c),
[Unicode conversion](https://raw.githubusercontent.com/python/cpython/v3.12.5/Objects/unicodeobject.c),
[UnicodeData 15.0.0](https://www.unicode.org/Public/15.0.0/ucd/UnicodeData.txt).


## Gaussian Constructor and Center State

**SageMath:** `DiscreteGaussianDistributionLatticeSampler(B, sigma=1, c=0, ...)`
converts sigma before the basis, accepts empty and dependent bases, and validates scalar
positivity when an integer sampler is needed. `set_c(None)` stores None and returns without
precomputation. Other center changes preserve the original cached trivial-lattice flag and
integer samplers unless a branch replaces them. This includes upstream's surprising reuse
of the simple sampler after changing an integral center to a fractional center.

**Port:** The options object and convenience factory default sigma to 1. Zero centers use
`0`, `0n` or a zero vector; explicit null defers precomputation. `c()` and `cNumeric()` now
return null for a deferred center. The low-level `_call()` returns literal zero on a basis
with no rows, as upstream does; `sampleExact()` preserves the corresponding upstream
attribute error when that result reaches the public sampling wrapper. Empty/dependent
basis state, row-space consistency, delayed assertion/index failures, signed-zero errors
and center cache transitions are compared against the bundled Python module.

**Rationale and trade-offs:** Dense arrays and exact Rational entries remain the port's
matrix/vector representation. Lattice row-space solves now use the shared matrix solver
instead of requiring an inverse; its existing backend limits remain. `sampleExact()` is
the Sage comparison surface when samples can be rational; `sample()` and `call()` retain
the documented bigint-only convenience restriction. The extra `tau` option and finite
number-based scalar API remain. Arbitrary Python matrix/coercion protocols and full
high-precision scalar/covariance propagation are not supplied by this change. Positive
infinity and NaN can be stored when precomputation does not construct an integer sampler;
the existing integer backend still rejects those parameters when it is reached.

**Behavioral impact:** Omitted sigma, null centers, degenerate bases and subsequent cache
transitions now agree in the compared domains. Callers must handle nullable center
accessors, hence version 8.0.0. Nonpositive integer-sampler errors use CPython's `%f`
nearest-even rounding at six decimal places, including large finite values and negative
zero/infinity. The sampler repr's distinct RealNumber/Decimal formatting remains separate.

**Oracle and harness:** The shared comparator now requires equal exception messages as
well as classes. Installed Sage 10.3 capitalizes the independent-basis error and adds a
period; bundled `modules/free_module.py:6738` uses lowercase without punctuation. The
lattice oracle adapts only that exact older message for `span_of_basis`/IntegerLattice
construction. The TypeScript implementation already matches the bundled spelling.


## Gaussian Parameter Representation

**SageMath:** Lattice parameter repr uses `RealNumber.str(truncate=True)` and MPFR's
nearest-even decimal digits. Integer-sampler repr uses `RealNumber.__format__`, which
formats `Decimal(repr(value))` to six places. This differs from both JavaScript `toFixed`
and the integer sampler's separate CPython `%f` validation message.

**Port:** Both sampler repr paths now use a shared 53-bit formatting adapter. It delegates
binary64 extraction/normalization/precision rounding to `mpfr-ts/mpfr_set_d` and decimal
digits to `mpfr_get_str`. Integer repr then quantizes those decimal digits with exact
nearest-even arithmetic. Scientific notation follows Sage's exponent threshold and spelling;
fixed notation preserves signed zero and expands large finite values.

**Rationale and trade-offs:** The scalar sampler API remains number-based, so this adapter
models `RealField(53)`. It does not supply arbitrary-format Python Decimal APIs or full
high-precision sampler parameters. MPFR setters retain the existing precision range
1..4096 and nearest-even-only port scope; global MPFR exception flags and custom exponent
limits are not exposed. The double setter returns the native inexact direction and preserves
the old destination sign when assigning NaN, which does not propagate the input NaN sign.

**Behavioral impact:** Parameter strings now agree in 1,385 original comparisons, including
798 repaired integer-sampler and 12 lattice-sampler representations. Sampling arithmetic
and streams are unchanged. Another 6,720 cases compare native MPFR 4.2.1 setter state,
round-trip doubles, decimal output and destination reuse, including subnormal and singular values.
Another 160 direct formatting comparisons cover both adapters across finite and singular inputs.


## Gaussian Constructor Validation and Integer Centers

**SageMath:** The integer sampler checks sigma, tau, algorithm, integral-center requirements
and precision in that order, then converts its center to the sigma's real field. Missing
sigma raises the constructor arity error. Integer centers round to the real field precision;
Python float tail cutoffs must also survive `Integer(tau)`, which rejects fractional values.

**Port:** The supported number/IntegerLike constructor inputs follow that order and preserve
exception classes/messages, including null centers. Integer centers round to binary64;
exposed integer bounds add/subtract their half-width with BigInt arithmetic. Fractional
number tail cutoffs now raise the original conversion error.

**Rationale and trade-offs:** This repairs observable constructor and sample behavior within
the existing 53-bit sampler API. Arbitrary-precision centers beyond binary64's exponent range
remain unsupported and raise an explicit NotImplementedError. The `dp` backend remains
unported. Center and tail coercion outside the declared input types remains under audit.
The fractional-tail rejection is raised after frontend/center validation, avoiding the
native temporary sampler allocation that precedes Sage's final `Integer(tau)` check.

**Behavioral impact:** 806 shared comparisons cover validation precedence, missing/null
inputs, integral/fractional tail cutoffs, large integer centers, algorithm selection,
exposed bounds and seeded streams. The bounds are port convenience accessors; the original
oracle derives them from the rounded center and native half-width expression.

**Oracle:** Installed Sage 10.3 adds a period to the unsupported-precision ValueError.
The bundled `discrete_gaussian_integer.pyx:403` omits it. The comparative dispatcher removes
only that exact older message's period before comparison.


## Native Real Numeric Construction

**SageMath:** `RealNumber._set` delegates Python floats to `mpfr_set_d`, integers to
`mpfr_set_z`, and strings to MPFR parsing. Every conversion uses the parent's precision.
Exact rational values, representations and precision changes observe the stored MPFR value.

**Port:** Nearest-even `RealNumber` construction now delegates doubles and BigInts to those
MPFR setters and strings to the existing native parsing adapter. Ordinary real values
retain MPFR state for representations, exact rationals, sign/mantissa/exponent, negation
and absolute value. The new `mpfr_set_z` follows `set_z_2exp.c` with exponent zero, including
nearest-even ties, mantissa carry, positive integer zero and native inexact direction.

**Rationale and trade-offs:** Reuse the dependency conversion layer instead of losing input
precision through `Number`/`parseFloat`. The existing nearest-even precision range 1..4096
applies. Directed-rounding construction and ordinary arithmetic retain their existing
approximate model; this change does not provide full MPFR arithmetic. Global MPFR exception
flags/custom exponent limits remain outside the exposed API. Native values beyond the
binary64 range retain their exact conversion state, but methods still using `_value` are
part of the continuing arithmetic/predicate audit.

**Behavioral impact:** 1,995 shared constructor/negation/absolute-value observations repair
1,661 old mismatches, including invalid string rejection, subnormals, signed/singular
values, integers through 4097 bits and precisions 2..4096. Another 4,272 direct native
MPFR 4.2.1 comparisons cover integer conversion at precisions 1..4096, decimal output,
round-trip double bits and exact/inexact status.


## Native Real Predicates and Integer Rounding

**SageMath:** Predicates inspect MPFR state, and real integer rounding delegates to
`mpfr_round/floor/ceil/trunc` followed by integer extraction. Finite values outside the
double exponent range remain finite; underflow in a double projection does not erase sign.
NaN and infinity integer conversions raise the corresponding Sage ValueError messages.

**Port:** Native-backed reals now follow these paths through the MPFR port. Classification,
sign, integer tests and comparisons with ±1 use the native significand and exponent.
Rounding operates directly on the integer/target-precision grid and supports destination
aliasing. `mpfr_get_z` returns `[bigint, inexactStatus]`, replacing the C output parameter.
Its temporary precision holds all integer bits, even beyond the public real precision cap.

**Rationale and trade-offs:** Preserve conversion state instead of treating its binary64
projection as the real value. The new integer-rounding routines support nearest-even,
toward zero, toward either infinity, away from zero and the round wrapper's nearest-away
mode. Faithful `RNDF` rounding and global MPFR exception/exponent-limit state remain
unsupported. Other MPFR setters/formatters retain their documented nearest-even scope.
This repairs predicates and integer rounding; general real arithmetic and directed-mode
construction still retain their existing approximation limits.

**Behavioral impact:** 3,600 Sage comparisons cover twelve public observations at precisions
1..4096, signed/singular inputs, values near ±1 and decimal exponents through ±4096.
The initial 576 probes exposed 100 discrepancies, all repaired. Another 18,300 native
MPFR 4.2.1 cases cover classification, signed-machine-integer comparison, integer extraction,
all supported rounding modes, narrow destinations, in-place operations and native inexact
status (magnitude 2 for fractional inputs, 1 for integers losing destination precision).
Another 1,584 native cases use binary exponents through ±100,000,000: rounding work stays
bounded by the input/output precision; only integer extraction allocates its integer result.

**Oracle bindings:** Installed Sage 10.3's Cython declarations omit `mpfr_roundeven` and
incorrectly discard `mpfr_get_z`'s return value. The oracle declares both directly from
`mpfr.h` and observes their actual native results without adapting their outputs.


## Native Real Fractional Parts and Comparisons

**SageMath:** Fractional extraction uses `mpfr_frac`, preserving signed zero and mapping
infinities to signed zero. Real comparisons use canonical common-parent coercion: two
MPFR fields meet at the lower precision, reparsing RealLiteral text when coercing it.
Python floats coerce into MPFR parents through 53 bits; higher-precision real operands
meet them in RDF, whose conversion invokes the value's `__float__` method.

**Port:** Native-backed fractional extraction and comparisons now use the MPFR port.
`mpfr_set` copies the native sign, even for NaN; `mpfr_frac` isolates fractional bits and
normalizes them; `mpfr_cmp` compares sign, exponent and aligned significands. Real comparison
helpers reproduce the precision/RDF conversion paths. `equals` guards NaNs, while the
port-only three-way `cmp` retains zero for unordered comparisons, matching its existing
sign-of-relations adapter.

**Rationale and trade-offs:** Preserve the stored value and canonical coercion behavior.
Copy/fraction conversion retains the existing nearest-even-only scope; global MPFR flags
and custom exponent bounds remain unexposed. Fractional halfway status follows MPFR 4.2.1's
64-bit limb implementation: the direct raw-rounding branch returns ±2 on exact ties, while
its temporary-buffer/set branch returns ±1. Ordinary arithmetic and directed-mode
construction remain under their existing approximation limits.

**Behavioral impact:** 31,500 shared comparisons cover 480 fractional-value observations,
26,712 field/literal/Python-float comparisons, 2,400 native copy/fraction/comparison cases,
1,728 native limb-boundary halfway cases and 180 double-rounding controls. The initial
912 probes exposed 54 mismatches; expanded float cases exposed 198 additional coercion
mismatches in the first fix, and two native cases exposed inexact-status differences.

**Oracle:** Installed Sage 10.3 lacks the bundled RealLiteral.__float__ override. For the
literal-to-RDF path the oracle evaluates its bundled `float(x.numerical_approx(53))`
expression before comparing. Native real-to-real literal coercion still uses RRtoRR.
All native copy/fraction/status observations call MPFR directly without output adapters.


## Gaussian Copy Helper

**SageMath:** There is no `withOptions` method. Constructing a fresh integer sampler with
an existing sampler's sigma/center/tau/algorithm and selected overrides follows the ordinary
constructor validation and initializes independent sampler caches.

**Port:** `withOptions` is a convenience adapter for that construction. Undefined overrides
retain the existing parameter; other supplied values, including precision, pass directly
to the constructor. It preserves real/large integer centers and original validation order.
Public center/tail option types now accept `IntegerLike | number`, matching existing
constructor behavior and allowing mixed numeric/integer option objects.

**Rationale and trade-offs:** Keep the existing convenience API while making its output
and seeded streams equivalent to explicit original construction. Invalid runtime nulls
receive constructor errors, and a null algorithm selects the original default. The `dp`
backend and existing constructor-domain limits remain unchanged.

**Behavioral impact:** 432 shared constructor/paired-stream cases expose and repair 225
previous discrepancies. The oracle constructs two original samplers and compares both
streams, checking independent caches and that copying leaves the original parameters intact.


### Successive differences recursion limit

Sage's `arith.misc.differences` recursively applies each requested difference, including
when the list is already empty. The port iterates and stops when empty. This avoids a
runtime-specific recursion limit and bounds work by the list length. The trade-off is
that extremely high positive orders return `[]` where a particular Python runtime could
raise `RecursionError`; positive orders within its recursion limit have identical outputs.
Nonpositive orders retain Sage's exact `ValueError` and message. Comparative cases cover
empty/singleton inputs, polynomial sequences, and orders below/equal to/above list length.


### PARI prime lookup sizing and cache

Sage's `nth_prime` delegates to PARI `prime(long)`, which uses a static prime/ordinal
table, a process-global prime cache, and a sieve or next-prime successor search. The port
now delegates to `parigp-ts/prime` and uses the original 64-bit static table and the same
search family. It does not reproduce the process-global PARI cache.

For interval sizing it replaces the floating-point logarithm estimate with an integer
bit-length upper bound. The sieve cost decision uses an integer square-root bound and
falls back to BPSW search for short intervals at large offsets. This keeps endpoints exact
and avoids runtime-dependent floating calculations. The trade-off is different cache,
interval and strategy performance; prime values and valid-index results do not change.
The existing probable-prime/proof limitations of the dependency still apply. Input C-long
conversion is represented by a bigint range check, matching the 64-bit oracle's errors.
Tests include every static table anchor and two successors through ordinal 100,000,000,000.


### Dedekind sum backend arithmetic

Sage's default/flint path calls `fmpq_dedekind_sum`; its pari option calls `sumdedekind`.
The port now delegates accordingly. FLINT's continued-fraction matrix and alternating
quotient recurrence is performed with BigInt. It includes the native `k <= 2`/zero-numerator
shortcuts and final rational normalization. PARI uses its Knuth recurrence with the original
64-bit word/generic threshold, preserving their different negative-denominator conventions.

**Rationale and trade-offs:** BigInt replaces native limb packing and signed two-word
products. The FLINT port currently processes exact Euclidean quotients individually;
Lehmer batching and recursive half-size quotient extraction from `get_cfrac_helpers.c`
remain performance gaps for huge inputs. It retains a logarithmic number of Euclidean
steps and does not enumerate denominator residues. Native quotient-batching performance
is not claimed. PARI's machine-word recurrence is represented in exact BigInt arithmetic.

**Behavioral impact:** numerator/denominator outputs and backend error classes/messages
match the original cases, including negative moduli at both sides of the word threshold.
The Sage wrapper retains its existing plain rational-pair object; dependencies return a
reduced numerator/denominator tuple. Comparative coverage includes bundled FLINT vectors,
noncoprime pairs, zero errors, word boundaries, Fibonacci quotient chains and 32,769-bit
inputs. Installed Sage 10.3's PARI subprocess wrapper is bypassed in the oracle in favor of
the direct cypari2 call used by bundled 10.9, preserving native `PariError` messages.


### Native MPFR binary arithmetic

**Original:** MPFR add/mul round exact mantissa sums/products to the destination precision,
including aliasing, singular signs and implementation-specific inexact magnitudes.

**Port:** `mpfr_add` and `mpfr_mul` now implement nearest-even arithmetic over native state.
BigInt replaces GMP limb storage; addition aligns significands with a sticky remainder
when exponent gaps exceed working precision, and multiplication retains the exact product.
The subtraction small-operand branch preserves MPFR 4.2.1's ±2 midpoint status.

**Rationale and trade-offs:** Keep storage bounded by operand/destination precision rather
than exponent distance. Native GMP's limb-specific multiplication optimizations are left
to BigInt. Directed modes, global flags and custom exponent ranges remain unsupported.
Sage configures MPFR exponents to ±(2^62−1); the port's number exponent representation
supports safe-integer exponents within that range and rejects unsafe results. General
RealNumber arithmetic retains its existing approximation limits; these primitives are
used by the real-field algebraic-dependency path.

**Behavioral impact:** 23,256 original native comparisons cover mixed precision, both
operand aliases, signed zero/NaN/infinity, wide exponents, cancellation, and midpoint
status across 64-bit limb boundaries. The initial binary cases exposed two status
mismatches, now repaired; 5,184 additional midpoint cases also pass unchanged.


### PARI algebraic dependencies

**Original:** Python floats become RDF, which takes PARI's algdep path, not Sage's RR
branch. PARI converts a double to a 64-bit real, constructs its relation lattice with
51 retained bits, and reduces at delta=.99. RealField inputs instead use Sage's lattice
at delta=.75 with their stored precision and the original precision-hint precedence.
Integer/rational values return exact linear relations. Height bounds/proof apply to
exact and real inputs; RDF rejects them. Every approximate path selects the best-fitting
irreducible factor.

**Port:** Number inputs now delegate to `parigp-ts/bibli1.algdep`; RealNumber callers pass
their full native value and use MPFR powers/factor evaluation. Exact shortcuts, null
height misses, proof errors and irreducible-factor selection are restored. The source's
symbolic/interval equality failure in a uniqueness check is preserved as an empty-message
TypeError. Polynomial coefficients remain ascending bigint arrays, the established port
representation.

**Rationale and trade-offs:** Double inputs delegate to the bundled PARI relation-basis
LLL algorithm. The separate MPFR/real proof path still uses the port's exact row
reduction rather than fpLLL. Its proof inequalities compare squared integer norms;
native RIF outward-rounding overlap near a strict proof boundary remains an audit gap.
Only the source's exact symbolic equality failure is currently reproduced. These details
can affect the selected relation or proof outcome in untested boundary cases and are not
claimed equivalent. ComplexNumber retains its separate precision-limited path; p-adic
inputs remain unsupported. Large-degree polynomial factorization inherits its existing
backend limits.

**Behavioral impact:** 2,028 float-dispatch cases exposed 1,724 prior mismatches. Another
7,098 exact/real comparisons cover degrees, precision hints, bounds, proofs and method
aliases; 1,260 direct PARI comparisons cover raw factors, random doubles, subnormals and
C-long argument boundaries. This corrects the earlier deviation entry's erroneous claim
that JavaScript numbers already took Sage's real-field path.


### PARI prime successor cache

**Original:** Integer next/previous prime methods delegate to PARI `nextprime`/`precprime`;
PARI uses cached small primes and a 210-residue wheel for larger values, then native
word-primality/BPSW testing. Sage's arithmetic `previous_prime` function instead walks
odd candidates and has a different lower-bound error from the Integer method.

**Port:** Integer successor/predecessor operations now delegate to PARI ports using the
original residue gaps modulo 210 and existing primality backend. The arithmetic previous
function retains its own odd-candidate loop and error. Prime-power predicates delegate to
the existing dependency implementation, removing duplicate extraction code. Prime-power
successor/predecessor functions now test odd candidates between adjacent powers of two,
returning the skipped power of two at the boundary. Eratosthenes uses its original odd sieve.

**Rationale and trade-offs:** Native process-global prime-table caching is not reproduced;
small cache-eligible values also take the wheel. BigInt arithmetic spans the native word
boundary directly. This affects performance, not the returned primes. Existing BPSW/proof
limitations and the dependency's exact-root/native word-path optimization differences remain.
Eratosthenes uses an exact integer square-root bound and requires safely indexable arrays.

**Behavioral impact:** 9,870 new original comparisons cover scalar/wrapped arguments,
prime-power bases and exponents, range endpoints, all three 210-residue cycles through
630, native 64-bit boundaries and larger inputs. They expose 32 prior mismatches: wrapped
large pseudoprime powers and the arithmetic previous-prime error. Public prime/power range
and traversal inputs now accept IntegerLike consistently with their runtime coercion.


### Hilbert symbol dependency domain

**Original:** Sage normalizes rational coefficients to numerator times denominator,
validates the prime place, and delegates to PARI by default. The direct algorithm and
cross-checking all option are separate paths. Unknown algorithms raise ValueError even
for zero coefficients or the real place. Conductor helpers coerce their inputs to ZZ.

**Port:** Those paths are now preserved. `parigp-ts/arith1.hilbert` implements the original
integer valuation-parity and Kronecker formulas. Its default/zero place denotes the real
place. The source's native rejection of negative/unit moduli precedes zero shortcuts;
composite-modulus behavior is also retained in direct dependency comparisons.

**Rationale and trade-offs:** The new dependency entry point currently accepts BigInt
coefficients and place only. PARI's real, rational, modular and p-adic native dispatch
remains outside that entry point's typed domain; Sage-level rational inputs are supported
through their exact source-defined integer square-class conversion. BigInt replaces GMP
limbs; the existing prime-validation/proof limitations remain. No output difference is
introduced for the supported integer/rational Sage inputs.

**Behavioral impact:** 22,724 new original comparisons cover algorithm selection, rational
normalization, wrapped/zero coefficients, signed valuation parity, native modulus errors,
and conductor/inverse validation and canonical pairs. A safe subset of 16,010 old-code
cases exposed 4,231 mismatches. Wrapped zero Hilbert inputs and the wrapped unit conductor
inverse could loop indefinitely; permanent cases now exercise both repaired paths.


The Hilbert backend uses `gen2.Z_pvalrem`: the native 16-division threshold switches to
recursive squared-divisor extraction for word-size bases. The two-adic path uses BigInt's
isolated low bit; wider bases retain the original repeated-division path. Another 336
comparisons verify signed units, composite bases and both native thresholds. This avoids
a per-exponent division loop where PARI uses divide-and-conquer. The internal helper's
source precondition is nonzero n and p>1; the port raises RangeError outside that domain
instead of invoking the original primitive with invalid arguments.

### Bounded native arithmetic and Cython error sentinels

**SageMath:** `get_gcd` and `get_inverse_mod` select `arith_int` through 46340,
`arith_llong` through 2147483647, and the generic function above that. Cython checks
32/64-bit arguments; the inverse implementations have different shortcuts and use C
signed remainder. A returned -1 collides with their `except -1` sentinel and raises
SystemError. CPython's diagnostic can contain a process-specific object address or read
`error return without exception set` after call-site specialization.

**Port:** The mirrored internal `rings/fast_arith.ts` implements the selected gcd and
inverse methods with the original Euclidean recurrences, conversion bounds, shortcuts,
errors and fresh bound-function identity. SystemError uses the stable shorter diagnostic;
the comparative oracle normalizes only that diagnostic. Inputs and outputs use the
project's IntegerLike/bigint mapping. The module is reached through the public arithmetic
factories; rational reconstruction and the separate public xgcd wrapper are not ported here.

**Rationale:** JavaScript cannot reproduce CPython object addresses or C undefined
behavior portably. Exact BigInt implements the defined native arithmetic safely.

**Trade-offs:** For remainder by zero, we reproduce the observed 64-bit ARM Sage oracle
(dividend retained), including its ensuing error-sentinel behavior. Other native targets
may trap or behave differently. Signed native intermediate overflow, including taking
the absolute value of the minimum signed integer, has no portable behavioral guarantee;
BigInt retains the exact result. These cases are outside the factories' advertised size
bounds and are not claimed to reproduce platform-specific undefined behavior.

**Behavioral impact:** Defined bounded inputs, conversion overflows, negative moduli and
function selection now match the original. The unstable SystemError text and native
undefined arithmetic retain the explicit platform limits above.

### CRT source version and integer-domain representation

**SageMath:** The bundled 10.9 bodies implement list CRT with a balanced binary tree,
a modular-element overload and a singleton parent coercion. CRT_basis supports a
non-coprime option and preserves any partial coprime prefix before its fallback.
The installed Sage 10.3 predates these interfaces and the singleton change.

**Port:** The integer/modular-element domain now follows the bundled algorithms,
errors and fallback output, including the original extra prefix coefficients.
The comparative oracle executes the four bundled function bodies unchanged using
Sage's arithmetic globals, linking CRT to the loaded crt body. Arrays represent Python
lists/tuples; IntegerMod and the two prime-field element backends represent Sage's
IntegerMod_abstract subclasses. A singleton modular input retains object identity.

**Rationale:** Version differences should not be mistaken for port discrepancies.
The prior intentional basis-prefix deviation is removed to restore source fidelity.

**Trade-offs:** Arbitrary Euclidean-domain and polynomial CRT coercion remains outside
these exported integer signatures. Parent-specific Python numeric subclasses map to
IntegerLike/bigint. The original non-coprime basis quirk can also cause CRT_vectors to
raise where a mathematically valid solution exists; this behavior is preserved.

**Behavioral impact:** Supported list forms, signed moduli, errors and source-specific
non-coprime output now agree with the bundled implementation. No substitute expected
arithmetic is computed in the oracle.

### Product trees and factor-base nontermination

**SageMath:** `rings/generic.py` implements generic-ring ProductTree layers, remainder
propagation and cached CRT interpolation, plus streaming product/derivative evaluation.
The bundled `smooth_part` computes remainders of the original input, then strips factors
from a running floor quotient. Factorization sorts and combines equal factors. With
repeated or overlapping bases this can produce a smooth product that does not divide
the input; coprime_part nevertheless returns its floor quotient. Some inputs never
terminate: units in the base or a zero running quotient cause an infinite extraction loop.

**Port:** ProductTree and prod_with_derivative implement the integer domain, with
IntegerLike inputs, bigint outputs, frozen arrays for tuple layers, __len__/__iter__
and JavaScript iteration. Remainder/interpolation trees and the streaming balanced
product stack follow the source. Smooth/coprime extraction accepts iterable or cached-tree
bases and preserves terminating source quirks, signs, sorting, repeated factors and
zero-base errors. Source-proven nontermination raises NotImplementedError with
`SAGE_NOT_IMPLEMENTED: smooth_part: original does not terminate for this factor base`.

**Rationale:** Cached trees restore the original algorithm and supported integer API.
A call that cannot produce an original result should not silently return a plausible
factorization or keep the JavaScript event loop blocked.

**Trade-offs:** General polynomial/Euclidean-domain trees remain outside the integer
signatures. The explicit error for a source infinite loop differs from Sage's absence of
a result. Layers map immutable Python tuples to frozen arrays, without Python parent
objects. General Factorization objects continue to use the project's array representation.

**Behavioral impact:** All tested terminating integer behavior follows the original,
including partial interpolation caches after errors and overlapping-factor quirks.
Eight fixed source-proven infinite-loop comparisons run the actual original with a
100 ms alarm and map only its timeout to the documented unsupported-case error.
This adapter does not apply to any terminating-output comparisons. Installed 10.3 lacks
smooth_part/coprime_part, so their bundled bodies execute unchanged in Sage's globals;
ProductTree/prod_with_derivative match the bundled AST after removing docstrings.

### Native integer totients, divisor counts and valuation routing

**SageMath/PARI:** Euler_Phi delegates positive values above two to PARI eulerphi;
number_of_divisors checks zero, then delegates to numdiv. PARI arith2.c uses native
factorization and products of prime-power contributions. Its own eulerphi(0) is 2,
while Sage's wrapper returns 0; native numdiv(0) raises PariError. Factorization and
quadratic roots use gen2.c's word/divide-and-conquer valuation routines.

**Port:** The Sage wrappers now delegate to parigp-ts/arith2.ts and preserve their
own boundary behavior. Native integer entry points retain signs, zero conventions,
factor contributions and balanced large-vector multiplication. ifactor.ts and qfb.ts
share gen2.ts valuation/unit extraction instead of separate repeated-division loops.

**Rationale:** Native arithmetic belongs in the dependency layer, and a shared native
valuation routine preserves the source's asymptotic behavior for large prime powers.

**Trade-offs:** The new entry points accept bigint, without PARI factorization-matrix
or `[integer, factorization]` argument forms. Native word packing and assembly optimizations
use BigInt and the existing Z_factor backend, with its previously documented primality
and factorization limits. Fixed-width multiplication tricks in small-vector products
are represented with exact integer products.

**Behavioral impact:** The tested integer values and errors agree with native PARI;
Sage wrapper conventions remain separate. Direct quadratic-root comparisons bind the
installed Zn_quad_roots C symbol under a test-local GP name and exercise only provably
solvable equations, since the raw function's NULL failure return is not a GP value.
Unsolvable cases remain covered through existing qfbsolve comparisons and unit tests.


### Square decomposition input mapping and native bounds

**SageMath:** Arithmetic two/three/four/k-square functions coerce n through ZZ and k
through Python int. Impossible decompositions raise ValueError. Below 2^32, arithmetic
delegates to rings/sum_of_squares.pyx; its public methods enforce Cython uint32 bounds.
Large three/four-square decomposition uses Integer.valuation(2)'s GMP low-bit scan.

**Port:** These signatures accept IntegerLike n and return bigint arrays for tuples;
sum_of_k_squares also retains numeric k and applies Python's finite-float truncation and
nonfinite errors. Native uint32 methods now live in the mirrored rings module and retain
original decomposition order and conversion errors on the 64-bit target. BigInt square
roots replace libc sqrt for this exactly representable small domain. Large power-of-four
extraction uses an isolated low bit and one shift. The former undocumented null failures
are removed; callers must catch ValueError, reflected in the major version change.

**Rationale:** Normalize before identity checks to prevent wrapped-zero loops and missed
perfect-square shortcuts, while restoring the source's errors, algorithm dispatch and
native API. Count normalization preserves existing numeric callers and adds exact counts.

**Trade-offs:** Generic ZZ-convertible strings, fractions and foreign scalar protocols
remain outside these IntegerLike signatures. Native unsigned-long overflow text follows
the 64-bit Sage target. BigInt arrays replace Python tuples and Sage Integer elements.
The preexisting warning for a conjecturally unreachable large-N brute-force fallback is
not emitted by the port; outputs and errors on the compared domain are unaffected.

**Behavioral impact:** Successful tuples, impossible-decomposition errors, native bounds
and count coercion agree in all added comparisons. The oracle loads bundled arithmetic
bodies unchanged because Sage 10.3's negative-k error has an older spelling; native Cython
methods use the installed original. Tuple-to-array mapping is explicit in the test adapter.


### Maximal-quotient reconstruction and exact source division

**SageMath:** mqrr_rational_reconstruction uses `/` on Sage integer elements in a pure
Python body. This produces an exact rational quotient despite the stale comment saying
C division implicitly floors. For integer inputs the first quotient makes the remainder
zero, so the function either selects `(u, 1)` or returns None, with a separate u=0 branch.

**Port:** The recurrence now uses the Rational implementation for quotient, threshold,
remainders and coefficient updates, preserving the original behavior and sign handling.
All three inputs normalize IntegerLike before zero/threshold checks; results map the
original integer pair to bigint[] and None to null.

**Rationale:** The audit targets executable original behavior, including source quirks;
implementing the described floor-division algorithm produces different reconstructions.

**Trade-offs:** Generic noninteger Sage scalar inputs are outside the IntegerLike
signature. This function retains the original limited usefulness; it does not substitute
a mathematically improved algorithm. The regular rational_reconstruction API is separate.

**Behavioral impact:** 44,332 original comparisons cover signed/zero, exact threshold
boundaries and integers through 1024 bits. They expose 18,629 previous mismatches, 4,973
inside the documented input domain. Installed and bundled original function bodies have
identical ASTs after removing docstrings, so this oracle uses installed Sage directly.


### Complex display ordering and Python sorting

**SageMath:** arith/misc.py:_key_complex_for_display classifies only exactly zero
imaginary parts as real. Non-real keys near zero use zero; other keys use
ar.n(digits=9), which numerical_approx.pxd converts to 34 binary bits. Tuple inputs sort
by their first entry, metadata and object identity survive, and empty input is returned
unchanged. The built-in Python stable sort determines the order of unordered NaN keys.

**Port:** The exported function accepts binary64 `{re, im}` values and tuple records
beginning with one. MPFR keys retain the native rounded value, including finite keys
above binary64's maximum after rounding. A CPython powersort adapter preserves natural
runs, binary insertion, merge scheduling, galloping and NaN comparison order. Nonempty
results are copies; empty results and all record identities are preserved.

**Rationale:** Decimal toPrecision and JavaScript's stable sorting algorithm produce
observable differences. Keeping native binary keys avoids collapsing rounded finite
values into infinity. The language adapter restores the source behavior in the supported
binary64 domain without a quadratic general sorting algorithm.

**Trade-offs:** General ComplexField/interval objects and their arbitrary precision
parents remain outside this numeric-pair signature. Unordered comparison results are
runtime-version dependent: the adapter follows bundled Sage's configured CPython 3.12.5.
Its eleven sort routines match installed Sage's CPython 3.11.8 after comments/whitespace
are removed, permitting a direct native oracle. A future Python algorithm change may
alter NaN ordering. The CPython license is retained beside the adapter.

**Behavioral impact:** 3,120 original display cases expose 2,969 old mismatches, including
tuple handling, rounding, empty identity and nonfinite inputs. Another 3,248 direct
Python-sort comparisons check results and comparison-trace hashes through 16,384 entries,
covering natural runs, minrun boundaries, both merge directions and galloping.


### Valuation dispatch and native GMP factor removal

**SageMath:** The free valuation function forwards to m.valuation and falls back to ZZ(m)
on AttributeError, including an AttributeError raised inside the method. Integer converts
the base first, then returns infinity for zero before testing p >= 2. Its finite path uses
GMP mpz_remove. Rational valuations delegate to numerator/denominator integers; val_unit
uses native removal and skips denominator extraction when the numerator has a factor.
Integer nbits delegates to bit_length, which reads GMP's stored magnitude-size metadata.

**Port:** Free integer valuations now return the existing 'Infinity' sentinel, with a
generic overload preserving an object's valuation result type. Rational and protocol
objects dispatch as in the source, including absent-method integer-coercion fallback,
single invocation of failing conversion hooks and exact error text. The internal positive-factor GMP adapter uses native BigInt, a low-bit scan for two,
and power-doubling extraction followed by descending powers for other bases. Integer and
Rational paths delegate accordingly; the existing Rational convenience aliases now accept
IntegerLike bases. Integer stores bit-length metadata at construction and reuses it when
copying an Integer wrapper, making bit_length and nbits constant-time queries.

**Rationale:** Restore executable method routing, zero/error precedence and native
algorithm complexity. The old free function's rejection of infinity and Rational values
was inconsistent with the original and with the port's scalar classes.

**Trade-offs:** Infinity remains a string, as in Integer/Rational. Native BigInt does not
expose limb counts, so wrapper construction computes magnitude metadata with a binary
string once; existing-Integer copies reuse the cached value. This assumes the declared
immutable value contract. The GMP adapter has the caller's positive-factor domain and
raises the Integer base error for invalid factors; it does not expose raw GMP's signed
factor/unit conventions or division-by-zero trap. Native limb division and word-size
shortcuts map to exact BigInt operations. padic/component valuation names remain convenience
aliases for the corresponding original operations, rather than extra Sage methods.

**Behavioral impact:** 32,307 original comparisons cover integer/rational wrappers, method
and AttributeError fallback, signed units, large prime powers and bit metadata through
262,144 bits. The initial 22,764 scalar cases expose 3,124 old mismatches; all 240 method-
protocol cases (including absent-method coercion) also differ, for 3,364 demonstrated old discrepancies. Direct native
factor-removal comparisons use gmpy2 linked to GMP 6.3.0, the bundled Sage dependency version.
The major version change reflects the free function's widened infinity return type.


### GMP factorial kernels and native resource bounds

**SageMath:** Integer.factorial checks sign and unsigned-long conversion, then delegates
to GMP mpz_fac_ui. GMP removes powers of two and uses odd-factorial product splitting
and prime-swing recursion, with a blocked sieve and platform-specific thresholds.

**Port:** Integer.factorial now uses the GMP adapter in types/gmp_factorial.ts and preserves
the bundled source's negative-input and 64-bit overflow errors. The adapter retains the
64-bit Apple M1 tuning (odd threshold 0, divide/swing threshold 252, product threshold 26),
small tables, packed factor products, both odd-factorial flags, 1,792 presieved bits and
2,048-limb sieve blocks. Two Uint32 words represent each original 64-bit sieve limb; masks
for 5/11 and 7/13 preserve the native bits and candidate order.

**Rationale:** Restore native dependency routing and asymptotic behavior instead of a
sequential growing-integer multiplication loop. BigInt represents the exact arithmetic;
Sage's configured 64-bit dependency determines conversion bounds.

**Trade-offs:** Native temporary-buffer reuse and mutation become immutable BigInt values.
The internal sieve returns its bit array instead of filling an output pointer and returning
a count. Internal unsigned-input/flag preconditions raise RangeError instead of exposing
GMP assertions or undefined C behavior. Runtime BigInt/array allocation limits remain: this
Bun runtime rejects a 100,000! result that native GMP can allocate. No invented mathematical
cutoff is imposed, and other runtimes may support larger results. Small constant tables
are generated once from the original products; the 64-bit popcount replaces its lookup.

**Behavioral impact:** 3,098 comparisons cover factorial/odd-part values through 60,000!,
both odd-factorial flags, exact sign/overflow errors and native sieve bytes through five
million. Full integer bytes are compared by bit length and SHA-256. Direct sieve tests call
the loaded GMP 6.3.0 internal symbol; factorial and double-factorial tests use gmpy2. The
Integer oracle changes only installed Sage 10.3's negative-error spelling to bundled 10.9.
Of 548 safe old Integer calls, 33 have the wrong error message. Another 436 safe free-
function calls (both wrapper modes) expose 164 ignored-selector mismatches, for 197 total.
Inputs >= 2^64 are excluded
only from the old nonterminating implementation, and included in permanent regressions.

### PARI factorial real representation and transcendental dependencies

**SageMath:** `arith/misc.py:factorial(..., 'pari')` returns cypari2's real Gen,
backed by `trans2.c:mpfactr`. Below 410 it rounds an exact PARI factorial; medium
inputs use precision-aware interval products; larger inputs use the gamma branch.
The free function uses PARI's global default precision, normally 64 bits.

**Port:** The free PARI selector now delegates to `parigp-ts:mpfactr`, returning an
immutable `MpReal<bigint>` record `{s,e,m,p}` representing `s*m*2^(e+1-p)`. The exponent
and mantissa are bigint so values beyond the safe Number exponent range remain
exactly representable. The GMP/default selector still returns bigint. The native
backend supports explicit precision in positive multiples of 64 bits. The free
selector uses 64 bits; there is no mutable cypari2 global-precision context or full
Gen method interface. Signed-C-long conversion and native exponent-overflow errors
are preserved at the free-function boundary.

**Rationale:** A rounded real is observably different from the old exact integer
return. The native representation preserves all result bits without constructing
an impossibly large exact integer. Existing quadratic-form MpReal users keep number
exponents through the generic type's default parameter.

**Trade-offs:** Gen methods and global precision settings remain outside this API.
The native helpers require their documented C input domains; invalid machine-word,
precision or nonempty powering domains use RangeError rather than undefined C
behavior. General qfb exponent-range limitations remain
as documented. New exponential, reciprocal, gamma, Bernoulli and AGM algorithms use
the original branches and operation order; comparisons cover their full output bits,
but do not prove all possible rounding boundaries. Native zero records retain allocation precision and their accuracy exponent; their
unspecified payload words are represented by a zero mantissa.
The internal `exp1r_abs` helper now retains a BigInt scale during repeated doubling,
including direct huge inputs and native exponent overflow. PARI 2.15's older routine
can exhaust its stack on inputs accepted by the bundled 2.18.1 routine. Comparisons
compile the bundled routine with symbol renaming and installed precision macros,
and are cross-checked against a full isolated build of the bundled PARI.
The bundled routine can also request initial series precision `l1` greater than its
allocated `L` before `setprec(X,l1)`, reading beyond the real allocation. Ten of 792
initial probes faulted in the full native build. The typed record has no adjacent
allocation to read and extends with zero bits. This preserves deterministic behavior
without emulating invalid memory access, but no equivalence claim is made for those
allocation-invalid paths; permanent comparisons constrain `l1 <= L` (or skip the
series branch). Generic wide-exponent arithmetic remains a separate audit gap.
The supplied Bernoulli constants are used lazily before filling the native cache;
this changes allocation timing without changing mathematical values.

**Behavioral impact:** The old free PARI path disagrees in all 1,012 nonnegative calls
of the 1,032 safe original-selector controls. New comparisons include small/medium/
gamma thresholds, native errors, exponents beyond 2^53, precision through 8,192 bits,
Bernoulli cache values and exact product/power operation traces. The installed PARI
uses word-length precision and older tuning constants: the Python oracle temporarily
sets its exported exponential/logarithm/reciprocal thresholds to the bundled 64-bit
GMP values (4224/384/4800 bits), then restores them. Internal real kernels are called
directly where the public wrapper would select a different function. No output
mantissas are changed to force agreement.

### PARI real-addition padding

**Original:** `kernel/none/add.c:addrr_sign` aligns mantissas in 64-bit words,
truncates carries and normalizes cancellation using its guard word. At a nonzero
whole-word exponent gap, when the smaller operand has fewer remaining words, the
`extend` branch reads `x[lx]` beyond that operand's logical length. The shift branch
that would initialize this guard only runs for nonzero residual bit shifts.

**Port:** `qfb.ts:addrr` now follows that word algorithm, including output precision,
zero accuracy and guard rounding. The unavailable whole-word guard is zero.
`addir` first chooses the original integer-conversion precision and delegates.

**Rationale:** BigInt has no adjacent allocation to read. Reproducing an out-of-range
native read would make a mathematical result depend on unrelated memory contents.

**Trade-offs:** This narrow upstream case has no unique native result to reproduce.
A direct C probe with a 64-bit smaller mantissa and a 192-bit larger mantissa at
exponent gap 64 returns mantissa `170141183460469231737836218407120622932` at
exponent 65 and precision 128 with zero padding; padding 2 changes its last bit to
`170141183460469231737836218407120622933`. The port chooses the former.

**Behavioral impact:** 57,825 exact native comparisons cover addition/subtraction,
integer and machine-word operands, signs, carries, cancellation and zero accuracy
at 64–256-bit precision. The Python oracle invokes the loaded PARI C kernels with
one explicitly allocated zero padding word after input reals and restores the
PARI stack after each call. It does not modify output frames. All cases match;
23,993 failed in the initial 57,820 cases before repair. Existing multiplication/division/transcendental
rounding limitations remain as documented separately.


## Native PARI real division kernels

**SageMath versus the port:** PARI divides native real mantissas using a single-word
branch, a short Knuth quotient algorithm or GMP division, depending on operand precision.
The port now preserves those paths, including their original guard-word rounding (which
is sometimes deliberately less accurate than rounding the exact rational quotient).
Real/integer division discards low denominator words before normalization, and
integer/real division uses an extra conversion word or the native reciprocal iteration.

**Rationale:** Correctly rounding an exact BigInt ratio changed native mantissa bits.
The port follows `kernel/gmp/mp.c:divrr/divri` and
`kernel/none/mp_indep.c:divir/divur/divru`, using the bundled GMP thresholds of 256 bits
for real division and 4,800 bits for reciprocal dispatch. BigInt quotient/remainder
operations replace GMP limb division while the native truncation and correction order
remain explicit. The unsigned-word API accepts bigint to preserve all 64 input bits.

**Trade-offs:** These C entry points require normalized whole-word real operands and
an unsigned divisor in `1..2^64-1`; `divru` rejects values outside that native domain.
The immutable record uses number exponents, so general exponent overflow remains
separate audit work. Zero conversions now retain their native allocation metadata. Existing inverse error objects retain their
documented missing-GEN rendering limitation. General exponent-range fidelity and direct huge exponential/resource behavior
are also still under audit.

**Behavioral impact:** The repaired kernels match 11,374 direct native comparisons,
including unequal precision, signs, zero numerators, full unsigned words, reciprocal
thresholds and quotient correction. The first 9,784 inputs exposed 943 old mismatches.
The oracle temporarily selects the bundled tuning thresholds in the installed library
and compares exact sign/exponent/mantissa/precision frames. No claim of exhaustive
input or branch coverage follows from that count.


## Native PARI allocated zero records

**SageMath versus the port:** PARI's `real_0_bit` has no mantissa allocation, whereas
`itor(0,p)`, `rtor` and `rdivii` can produce zeros with allocated precision. The port now
preserves that distinction. Real conversion sets zero accuracy to `min(old_exponent,-p)`;
integer conversion sets it to `-p`. Integer multiplication and division consequently
use the original allocated precision when deriving zero accuracy. Inverse validation
checks the native allocation before dispatch, so some failures originate in `divrr`
and zero divided by an allocated zero returns the original inexact zero.

**Rationale:** Allocation is observable through `realprec` and later arithmetic;
canonicalizing all zeros changed outputs. The immutable record stores the meaningful
allocation/exponent fields and uses `m=0n` for unspecified native payload words.

**Trade-offs:** A large allocated zero passed to the reciprocal Newton branch has no
normalized, initialized native mantissa. The port retains a deterministic inverse error
for that input instead of reading arbitrary storage. General C pointer identity and
uninitialized payload bytes are not modeled. Inverse comparisons retain the existing
error-class mapping and omitted GEN rendering, checking the actual originating kernel;
all successful results compare sign, exponent, mantissa and allocation exactly.

**Behavioral impact:** 2,052 new native comparisons expose 1,089 old mismatches across
conversion, integer operations, inverse dispatch and zero rational conversion. The native
oracle runs under cypari2's signal handler, allowing errors to be observed safely.
An older serializer incorrectly erased allocated precision from zero outputs; fixing it
revealed and repaired 90 additional existing `rdivii` comparisons. Its expected frames
were regenerated from PARI, not adjusted to the TypeScript outputs.


## Native PARI real multiplication kernels

**SageMath versus the port:** PARI uses truncated limb products for small real
multiplication and squaring, with different thresholds for the full integer kernels.
The port now keeps the original guard accumulation, normalization, carry rounding,
extra word for unequal precision, integer conversion precision and word shortcuts.
Pointer identity selects squaring before general multiplication, as in the original.

**Rationale:** Rounding an exact product can differ from rounding PARI's truncated
intermediate product. Native comparisons near half-ulp boundaries expose these
one-bit differences even when ordinary inputs all agree. The port follows
`kernel/none/mp_indep.c:139-449`, using bundled GMP thresholds of 3,520 bits for
multiplication and 768 bits for squaring. The shared small-product code reproduces
the arithmetic of both source loops; BigInt handles full integer products and limbs.

**Trade-offs:** Native whole-word normalized inputs are required. Immutable real
records retain pointer-identity dispatch but do not expose mutable C buffers.
BigInt supplies the integer multiplication implementation, so its machine-level
tuning differs from GMP. General exponent overflow remains separate fidelity work. The documented addition-padding choice still applies to dependents.

**Behavioral impact:** 21,466 direct native comparisons cover ordinary and constructed
rounding boundaries, signed aliases, mixed integers, square/product identity and
rounding carry. The old implementation mismatches 3,404 of these cases. For example,
at 832 bits the square of `(5*2^829-3)*2^-831` differs by one mantissa bit depending
on whether the same operand object or two distinct equal objects are passed to
`mulrr`; this matches native kernel selection rather than imposing a common result.


## Native PARI square-root results and GMP dependency

**SageMath versus the port:** The exported PARI real square root now follows
`kernel/gmp/mp.c:sqrtr_abs` and `pariinl.h:sqrtr`: it keeps the odd/even exponent
remainder tests, returns a complex value for negative input, and halves the full
zero exponent with signed floor division. The complex result has an exact bigint
zero real component. Both low-level `sqrti` copies now read integer magnitude, as
PARI's unchecked integer kernel does, rather than rejecting negative signs.

**Rationale:** Correctly rounding a mathematical square root is not always the native
result. For example, at 64 bits the even-exponent guard correction rounds
`sqrt(1+2^-63)` upward. The port preserves this result and delegates its integer
root/remainder to the GMP 6.3.0 Karatsuba algorithm, including normalized reciprocal
table seeds and quotient correction. It avoids a full-precision Newton iteration
at every step. Architectural component and output-pointer mappings are in DESIGN.md.

**Trade-offs:** Native real precision is a positive multiple of 64 bits and stored
square-root exponents remain numbers. General exponent overflow remains audit work;
the Buchmann adapter now shares this kernel. `sqrtr_abs(0)` retains the safe typed zero extension even
though the original low-level helper requires nonzero input. GMP helpers use packed
BigInt carries rather than mutable limb buffers. Root-only wrappers compute and discard
a remainder using the same divide-and-conquer family, with extra work relative to GMP's
root-only optimization; machine-specific tuning and assembly are not reproduced.

**Behavioral impact:** 6,276 new shared cases include 66,048 integer root/remainder
values compared through hashes. The first 2,932 real cases expose 864 old mismatches:
720 negative returns, 120 rounding results and 24 wide zero exponents. Ten additional
negative integer-root comparisons expose the two old sign checks. Exact components,
precision, exponent and remainder are retained. A subsequent direct Buchmann audit
also replaced its independent real-root helper with this shared kernel.


## Native real-to-integer error exponents

**SageMath versus the port:** `gen3.c:gcvtoi` first tests the input exponent, then
computes the integral precision from the stored allocation. When subtraction is
needed, its error estimate is the exponent of the native `subri` result, including
finite-accuracy cancellation to zero. The qfb helper now follows this order; the
Buchmann helper preserves the corresponding cancellation and zero metadata instead
of substituting a negative sentinel or returning the input zero exponent too early.

**Rationale:** A real record containing an exact integer still has finite declared
accuracy. The old sentinel claimed much stronger accuracy than PARI, while some zero
inputs followed the wrong branch. These values guide subsequent precision decisions.

**Trade-offs:** The qfb API returns `[integer,errorExponent]`, whereas the existing
Buchmann adapter returns `{z,e}`. Both retain number exponents and their documented
input representations. Buchmann arithmetic now delegates to shared kernels; its
requested-bit conveniences round up to whole-word precision.

**Behavioral impact:** 4,012 direct native comparisons cover signs, fractional/exact
values, precision boundaries, both kinds of zero allocation and wide zero exponents.
They expose 798 old mismatches, split evenly between the two helpers. For example,
`gcvtoi(itor(1,64))` returns integer 1 with error exponent -63, while a canonical zero
with exponent 0 returns error exponent 1. Public and focused regressions preserve both.


## Buchmann binary64 working precision

**SageMath versus the port:** Native PARI `dbltor` accepts one double and creates one
mantissa word, or a canonical zero with accuracy exponent -1023. The Buchmann wrapper
retains its optional working-precision argument. When supplied, it means
`rtor(dbltor(d), nbits2prec(p))`: requested bits round up to complete 64-bit words.
Both binary64 wrappers now share the native conversion kernel, including rounding,
subnormal behavior and `PariError` overflow exceptions. `cmprr` also delegates to the
native comparison kernel, including approximate-zero comparisons.

**Rationale:** The precision argument is used by existing Buchmann callers. Composing
native conversions preserves that convenience without another rounding implementation.

**Trade-offs:** This extra argument has no direct native `dbltor` counterpart; requests
such as 200 bits allocate 256 bits. Buchmann arithmetic, square roots and
transcendental helpers now delegate to the shared kernels.

**Behavioral impact:** Omitted precision now retains the original zero allocation and
accuracy. Explicit precision preserves the original zero accuracy when possible, and
nonfinite inputs/overflow report `PariError`. Real comparison follows native zero accuracy
and compares huge exponent gaps without allocating an exponent-sized mantissa. Shared
native comparisons exercise both the raw conversion and precision-composition paths.


## Buchmann real working precision

**SageMath versus the port:** Buchmann real constructors and elementary arithmetic now
use the shared native PARI kernels. `real_0_bit` has the native one-argument signature
and canonical allocation. Zero integer conversion retains allocated precision. The
existing constructor defaults and immutable `setprec` helper are TypeScript conveniences:
explicit requested bits round up to a multiple of 64, and `setprec` composes native
`rtor` rather than changing a C buffer's header in place.

**Rationale:** Shared kernels preserve native truncation, rounding and allocation while
removing a second arithmetic implementation. Whole-word working precision supports
existing callers requesting, for example, 200 or 260 bits.

**Trade-offs:** Requested bits can round up; `real_0_bit` no longer accepts the old extra
allocation argument. This signature correction is released as version 13.0.0.
The former independent square-root/logarithm/exponential helpers now delegate to the shared native kernels. General wide-exponent and huge-exponential resource behavior remains under audit.
Native inverse errors retain the previously documented exception-subclass/GEN-rendering
mapping; comparative zero tests compare the preserved inverse origin.

**Behavioral impact:** Native mixed precision, cancellation accuracy, truncated products,
quotient rounding and unsigned-word operands replace the older approximate behavior.
`truncr` now raises native precision errors for a nonzero real whose integer part
exceeds its mantissa precision. `gcvtoi` and `trunc2nr` preserve their distinct native
precision-loss semantics. Algebraic dependencies use `trunc2nr` as upstream does.


## Shared native logarithms and Buchmann transcendental results

**SageMath versus the port:** Buchmann square roots and transcendentals now use the
same native kernels as qfb. Negative real roots return a complex record. Public
`mpexp` (Buchmann `expr`) returns `MpReal<bigint>` so the native signed exponent is
preserved; `rtodbl` accepts both exponent types. The shared logarithm implementation
uses PARI's progressive series precision, native double estimate, integer loop-bound
arithmetic and final precision clamping. Its log(2) cache uses native `atanhuu` and
Haible/Papanikolaou binary splitting with one guard word.

**Rationale:** Remove independent approximations and restore the original algorithms,
intermediate precision and result types. BigInt avoids loss of native exponent bits.

**Trade-offs:** The Buchmann square-root union and public exponential exponent type
are breaking corrections, released as 14.0.0. Other arithmetic entry points still
use numeric exponents; conversion back to that domain requires an actually safe
integer exponent. Native stack-buffer mutations are represented by immutable records,
with explicit truncation or precision conversion as appropriate.

**Behavioral impact:** Complex roots, wide zero accuracy, logarithm rounding and
exponential precision growth now match direct native probes. General huge-exponential
resource behavior remains unresolved. `logr_abs` assumes nonzero input in C; the typed
port retains a deterministic `PariDomainError` guard for zero rather than reading an
absent/unspecified mantissa. That guard is a defined extension outside the native
helper's input precondition and has no native result to compare.


The binary-splitting oracle compiles the bundled `trans2.c:get_nmax` and `atanhuu`
functions with only symbol renaming and the installed precision ABI. Installed PARI
2.15.4 reports stack allocation errors on ten rounded-equal word ratios where the
bundled 2.18.1 source reports `overflow in atanhuu`. All 410 other split/atanh cases
also match the installed native functions directly. Generic split results call native
`abpq_init`/`abpq_sum` through explicit C struct layouts and preserve complete integers.


### Fraction string fallback and source version (17.0.0)

**SageMath behavior:** Fraction construction first attempts polynomial conversion, then
uses `sage_eval` with fraction generators. Compilation precedes unknown-name lookup;
only `NameError` becomes `unable to evaluate ... in ...`. The bundled preparser no
longer implements Sage 10.3's backslash operator.

**Port behavior:** The shared arithmetic parser supports this fallback, numerator and
explicit denominator expressions, Boolean results and the compared syntax diagnostics.
It performs a syntax pass before evaluation. Full Python evaluation, Sage globals,
string-literal expressions, symbolic operations and general function calls are not ported.

**Rationale and trade-offs:** Explicit arithmetic callbacks replace Python operator
and evaluation dispatch. This preserves the tested arithmetic subset without executing
arbitrary code, but does not provide a general `sage_eval` implementation.

**Behavioral impact:** Compared rational expressions over QQ and the tested prime fields
match. The oracle executes the bundled `repl/preparse.py` during native string fallback;
46 backslash cases differ from the installed Sage 10.3 preparser. Explicit FpT zero
denominators remain governed by the existing safe-construction deviation; 93 such
probe rows are excluded from the live suite after the native process crashed.


### Infinite order basis and Sage runtime compatibility (18.0.0)

**SageMath behavior:** Bundled `order_rational.py:516` returns the Python integer tuple
`(1,)` for the infinite maximal order's basis. Installed Sage 10.3 instead returns the
function-field element `1/x`, contrary to its own docstring.

**Port behavior:** The public return type is now `[bigint]`, returning `[1n]`. The
finite order retains a basis of function-field elements. The oracle explicitly uses
the bundled method's literal tuple for the infinite order.

**Rationale and trade-offs:** Follow the bundled original and the project's native
integer mapping. Callers that previously received function-field elements from the
infinite basis must explicitly coerce its integer into their field when needed.

**Behavioral impact:** The return container, integer entry and tested basis operations
match the bundled source. This adapter intentionally does not reproduce the older
runtime's incorrect reciprocal. Ideal monoids now preserve their native unique-parent
identity, and finite/infinite orders accept each other's ideals by their generators.


## FLINT word square roots and arithmetic

**SageMath/FLINT behavior:** `ulong_extras/sqrtmod.c` selects a bounded small-modulus
search, congruence shortcuts or Tonelli iteration. Its dependencies use binary word
Jacobi arithmetic, normalized reciprocal-based products, and square-residue filters
followed by a binary64 square-root candidate.

**Port behavior:** Root selection, all algorithm branches and the native Tonelli
iteration limit are preserved. Word operands use bounded BigInt. Jacobi follows the
original binary subtraction algorithm; limb preinversion uses the native-division
formula. Modular powering follows the original binary schedule with exact BigInt
modular products, without consuming the reciprocal argument. The square checker
retains the original residue filters and computes its candidate by integer Newton
iteration instead of converting to binary64.

**Rationale:** Keep integer arithmetic exact and portable in JavaScript, without native
limb or floating-point dependence. BigInt supplies the underlying division/product
operations; the library-level algorithms and root representative choices remain intact.

**Trade-offs:** Hardware reciprocal optimization and constant-time word square-root
approximation are not reproduced. Incorrect precomputed reciprocals are outside the
supported native contract. Explicit RangeError guards replace undefined/unchecked
native behavior outside word bounds, positive moduli or Jacobi denominator conditions.

**Behavioral impact:** All 69,725 direct native comparisons match, covering every new
kernel and 64-bit boundaries, including source-defined composite-modulus guard paths.
This is not a claim that `n_sqrtmod` solves arbitrary composite-modulus equations; its
contract requires a prime modulus. Polynomial/function-field root routing follows in
the next audit batch.


## FLINT polynomial square-root kernels

**SageMath/FLINT behavior:** Small-prime rational-function roots call FLINT polynomial
square roots. FLINT removes even valuation, computes a truncated series root, and
certifies the remaining high coefficients. Generic-ring series kernels use basecase
recurrences or Newton inverse roots with the final Karp–Markstein correction.

**Port behavior:** Dense BigInt arrays implement these algorithms, including the
portable `gr/nmod.c` cutoff table and classical/packed high-product dispatch. The
`gr_poly` functions currently instantiate only prime modular coefficient rings; they
are not general context-based GR APIs. Existing BigInt product kernels replace
machine-limb operations. Zero denotes the zero polynomial as `[]`; full-root failure
is `null`. Native GR status codes are represented by exceptions outside the supported
series preconditions (an invertible square constant and invertible 2).

**Rationale and trade-offs:** Reuse the port's exact dense representation and native
algorithm structure without C contexts, output buffers or aborts. The complete
`nmod_poly` class facade remains unimplemented; these are callable dense kernels.
Sage's general polynomial squarefree dispatch outside the repaired integer path
still has previously existing backend-routing gaps and remains under audit.

**Behavioral impact:** 18,285 native comparisons cover complete roots, series,
basecase/Newton transitions and high products, including native low-coefficient
behavior. The oracle calls installed FLINT through Cython. For length-zero inverse
series it uses the bundled GR basecase's explicit empty result; the installed public
nmod wrapper instead aborts. Invalid native preconditions are not compared as a
portable exception contract.


## Fraction field square roots

**SageMath behavior:** Generic fractions test the product of numerator and denominator;
FpT instead checks both raw polynomials using FLINT, makes the root denominator monic,
and chooses the smaller leading numerator coefficient. FpT retains raw zero identity.
Function-field square roots always wrap the returned fractions in new elements.

**Port behavior:** These paths, root lists and unnamed-extension errors now match.
`Polynomial.is_square(root)` supplies the generic fraction's squarefree-based test;
its optional tuple uses `null` for no root. Generic fraction `sqrt(extend, all, name)`
accepts the native arguments, but a named nonsquare extension raises an explicit
`NotImplementedError`; quotient extensions over fraction fields are not implemented.

**Rationale and trade-offs:** Restore the actual underlying fraction dispatch and
representations. This does not introduce symbolic or algebraic extensions, so callers
requesting a named nonsquare root cannot obtain Sage's quotient-field element yet.

**Behavioral impact:** 25,488 fraction comparisons exercise both normalized and raw
representations, Boolean results, exact roots, root ordering, identity, and exception
class/message. A further 5,580 polynomial comparisons include root parent identity.
The integer polynomial cases use the already documented Integer-object coefficient
adapter: the public ZZ parent returns raw BigInt, which is not a RingElement.


## NTL integer polynomial squarefree decomposition

**SageMath/NTL behavior:** Integer polynomial squarefree decomposition removes signed
content and delegates to NTL `SquareFreeDecomp`. NTL computes integer GCDs using
modular GCD, balanced CRT and exact divisibility, and selects plain or homomorphic
exact division according to polynomial degrees.

**Port behavior:** The dense routines retain signed content, modular word GCD,
balanced CRT and exact-division certification. `GCD(a,b,state?)` and
`SquareFreeDecomp(f,state?)` accept a shared PolynomialProductState. Each modular
iteration initializes its FFT prime before rejecting leading coefficients, then
uses the native FFT-prime word context. Both certification divisions and all
squarefree stages retain that state. Zero shortcuts do not initialize primes.
SquareFreeDecomp requires primitive input with positive leading coefficient, as
specified in NTL's ZZXFactoring documentation.

**Rationale:** Explicit moduli and state replace native ambient contexts without
changing the modular reconstruction or exact-division algorithms.

**Trade-offs:** Omitting state creates a private zero-key stream/cache; it cannot
continue a caller's native ambient stream. The native profile is the bundled
60-bit, nonthreaded build; other FFT-prime configurations are not modeled.

**Behavioral impact:** Coefficients, multiplicities, cache contents and stream
continuation match the native comparisons, including skipped primes, modular
degree changes, stable CRT candidates rejected by either certification division,
and repeated calls. Native context backup/restore becomes explicit state passing.


## Polynomial squarefree backend dispatch

**SageMath behavior:** Integer polynomials use NTL; word modular polynomials use
FLINT's squarefree factorization. Generic characteristic-zero fields preserve the
nonmonic first factor produced by the category's successive-GCD algorithm. Finite
fields use the characteristic-aware algorithm and inverse Frobenius on coefficients.
IntegerModRing parents without the field-base helper report NotImplementedError on
the generic binary/large-modulus path, even though ordinary factorization can work.

**Port behavior:** These squarefree paths and their zero/composite guards now match.
The new dense `nmod_poly_factor_squarefree(f, p)` returns factor/multiplicity pairs
in native order for a prime word modulus. Category/finite-field helpers live in their
mirrored source modules; explicit polynomial and characteristic/order arguments replace
the native category context. The existing public array representation carries a unit
as a constant polynomial pair when needed; the ZZ zero decomposition is `[[zero, 1]]`.

**Rationale and trade-offs:** Retain native algorithms, units and coefficient parents
without implementing Sage's dynamic category machinery or Factorization class. Generic
IntegerModRing factorization keeps its pre-existing finite-field factorization path
while its separate PARI backend is audited; that remaining routing gap is not closed
by this squarefree change. Coefficient-ring category mutation from an earlier explicit
field certification is not modeled by the current TypeScript parent classes.

**Behavioral impact:** Shared comparisons cover 2,128 integer/rational/modular inputs,
890 direct FLINT inputs and 2,187 actual finite-field inputs including extension
coefficients. Installed Sage 10.3 omits the bundled word-polynomial zero check; the
oracle explicitly executes the bundled ArithmeticError guard before calling the native
method. Integer zero decompositions and generic ValueError paths execute the original
runtime directly. Backend field guards remain ahead of constant shortcuts.


## Binary field coefficient roots and caching

**SageMath behavior:** GF(2) elements use IntegerMod_int's modulus table. Arithmetic and
square roots return its cached entries; negation preserves an uncached zero. Parent
conversion of an existing element preserves that object. Both residues are squares.

**Port behavior:** The dedicated GF2 coefficient class now exposes `is_square()` and
`sqrt({extend, all})`, and its parent retains the two canonical coefficients. Arithmetic,
roots, distinguished constants and iteration use those entries, except native zero
negation and existing-element conversion. Direct GF2Element construction remains fresh.

**Rationale and trade-offs:** A two-element parent-owned tuple replaces the native
modulus table. This preserves the compared identities with no native allocation API.
JavaScript integral numbers and BigInt denote Sage integers under the existing input
mapping; Python-int-specific allocation is not a separate TypeScript input domain.

**Behavioral impact:** All 1,582 native comparisons cover roots, polynomial root
reconstruction, errors, arithmetic and identities for canonical and directly constructed
coefficients. The oracle uses the actual IntegerMod_int class, not the unrelated
category-generated `parent.element_class`, and explicit Sage Integer constants.


### PARI squarefree zero-input safety

Native `FpX_factor_squarefree` requires a valid dense finite-field polynomial and a
prime modulus. Word-prime zero input raises `PariError` with the native inverse-error
message; the port matches it. Large-prime zero input reaches invalid native polynomial
division and can segfault. The port raises `RangeError` before entering that branch.

**Rationale:** an unsafe C input has no useful polynomial output to preserve.
**Trade-off:** consumers cannot reproduce the native crash. **Behavioral impact:**
only invalid large-prime zero input differs; valid constants, units, squarefree and
repeated factors match direct PARI calls, including the differing constant results
of the word and large-prime branches.

## Negative Polynomial Monomial Degrees

**SageMath:** `PolynomialRing.monomial(n)` documents nonnegative exponents, but its
unchecked dictionary constructor passes negative values into native backends.
QQ's FLINT coefficient setter can segfault; word-modular FLINT powering can abort,
segfault or request an enormous allocation. NTL GF2X/ZZ_pEX returns zero, while the
generic dense constructor raises `IndexError('list assignment index out of range')`.

**Port:** Preserve the terminating NTL/generic outcomes, including dictionary
construction through `PolynomialRing.__call__`. For negative QQ or odd
word-modular FLINT exponents, raise
`RangeError('negative monomial degree is outside the native FLINT contract')`.

**Rationale:** Native invalid memory access and process aborts have no deterministic
mathematical result to reproduce. **Trade-off:** failure type differs outside the
specified nonnegative-input domain. **Behavioral impact:** valid monomials are
unchanged; invalid FLINT inputs terminate predictably. Shared oracle cases explicitly
substitute this guard for those unsafe native calls; all other cases execute Sage.

## FLINT Random Bit-Count Bounds

**SageMath/FLINT:** The word random helpers assume bit counts from 0 through the
native word width; larger or negative counts reach nonportable C shifts.
**Port:** The 64-bit port raises `RangeError('bits must be between 0 and 64')` for
invalid or nonintegral bit counts, before consuming random state.
**Rationale:** Undefined shifts do not supply a deterministic result to reproduce.
**Trade-off:** These invalid calls terminate with a portable exception.
**Behavioral impact:** All valid bit counts and state transitions match the bundled
64-bit implementation; the comparative oracle substitutes this guard only outside
that domain.

## Native FLINT factor array kernels

The native `nmod_poly_factor` word-prime path now delegates to `flint-ts`, including
FLINT's deflation, squarefree, Cantor–Zassenhaus and Kaltofen–Shoup dispatch.
The latter retains baby/giant distinct-degree splitting, vector Brent–Kung
composition and modular matrix products. Matrix multiplication uses the portable
single-thread classical/Strassen cutoffs and Bodrato schedule; row evaluation uses
FLINT's Horner/rectangular split. BLAS and native threading are not available in
this JavaScript backend. The asymptotic portable algorithms are preserved, but
native tuning and parallel throughput are lost. Valid finite-field outputs and
native insertion order are compared against the bundled FLINT C library.

The C mutable polynomial/factor/matrix outputs are represented by dense arrays and
returned tuples; `nmod_poly_remove` returns `[remaining, multiplicity]`, factor
routines return factor/exponent pairs, and the three complete factor entry points
return `[leadingCoefficient, factors]`. Distinct-degree pairs carry the common
irreducible degree in their second entry. Probabilistic equal-degree splitting
returns `null` on failure instead of exposing the C scratch candidate; the supplied
random state still consumes exactly the native stream. These adapters make native
ownership explicit in TypeScript; caller aliasing and unused C output storage
cannot be observed. Valid mathematical results are unchanged.

Invalid lengths, strides, dimensions and equal-degree inputs outside native
preconditions have `RangeError` guards where the C code aborts, accesses invalid
memory or cannot terminate. In particular deflation must be a positive safe integer,
inflation nonnegative, and removal requires a positive-degree divisor. Equal-degree
probing requires degree at least two and positive `d`; recursive splitting requires
positive `d` dividing the input degree. It retains the native caller precondition
that the input is squarefree and all irreducible factors have degree `d`, rather
than adding a second factorization to validate it. Matrix and composition guards
reject incompatible dimensions or outer degrees at least the modulus degree.
This avoids undefined native behavior; error classes/messages outside the native
contract differ. Comparative tests call C only for valid inputs and explicitly
model the documented guards otherwise.

These changes repair the word-prime branch (odd primes below 2^63) selected by
Sage's polynomial backend. Binary and larger-prime factor routes now delegate to
PARI; extension factor routing and
prime-power/p-adic factorization, are separate audit work and remain open.

## PARI word random-state adapter

`parigp-ts` now provides the bundled PARI 64-bit XORGEN generator, global
`pari_init_rand`/`setrand`/`getrand` state, word and arbitrary-integer rejection
sampling, and polynomial/vector sampling. Native 64-bit wraparound, state-word
ordering, seed validation, discarded initialization blocks and random consumption
are preserved. Saved states are native positive integers and can be round-tripped;
`random_bits(64)` and `random_zv` preserve PARI's signed machine-word results.

PARI polynomial variable identifiers belong to its GEN representation. The array
APIs `random_Flx(length,p)` and `random_FpX(length,p)` omit that identifier and return
ascending reduced coefficients; `random_F2x(length)` uses a bigint with coefficient
`i` in bit `i`. This matches the port's existing coefficient adapters without
introducing GEN ownership. Variable-name metadata and native object aliasing are
not represented; the sampled coefficients and saved random state are unchanged.

Native random routines assume positive limits and valid machine-word lengths.
Out-of-contract limits/lengths receive pre-consumption `RangeError` guards in the
port instead of undefined shifts, invalid allocation or nontermination. In
particular `random_bits(0)` shifts a native 64-bit value by 64 and has no portable C
result; the port accepts widths 1 through 64 and guards zero. Empty vectors and
polynomials consume no randomness, and zero-length polynomial sampling does not
validate an unused modulus, as in PARI. Native `setrand` domain errors retain their
messages and the `PariError` adapter. Tests execute original C for defined inputs
and model only these documented guard boundaries.

The global state is local to one JavaScript module instance rather than native
thread-local storage. Separate worker/module instances are independent; multiple
callers in one instance share it. Existing local Las Vegas generators elsewhere
in the port have not yet been replaced and remain separately documented audit
work. This dependency alone does not establish generic polynomial-factor routing.

## PARI packed binary-polynomial kernels

Binary `Polynomial.factor()` now follows the generic Sage factor method's PARI
route. Its backend retains PARI's direct degree-at-most-two cases, Cantor pipeline
through degree 20, and Berlekamp pipeline above degree 20. Squarefree/DDF splitting,
EDF traces, native kernel basis order and random draws are compared to the bundled
C library, including complete RNG state. The Sage-facing results retain Sage's
factor ordering and existing unit/zero adapters. Large-prime factorization also
delegates to PARI. Extension factor routes and the standalone GF2X convenience
helpers remain separate audit work.

`F2x` values use nonnegative bigint coefficient bits rather than PARI GEN word
vectors. Variable identifiers, unused trailing storage and object identity are
omitted. Native arithmetic is preserved: 64-bit GMP multiplication cutoffs at 11
and 41 words select basecase, Karatsuba and integer packing; square/root operations
interleave/deinterleave bits; division and GCD follow the original binary kernels.
This keeps the portable algorithm complexity while relying on native BigInt for
integer multiplication instead of GMP. The results have no variable metadata.

`F2m_ker_sp(columns, rows, deplin)` receives packed column values and mutates their
array slots. Zero flag returns a kernel-basis array; nonzero flag returns the first
dependency vector or null. The native empty-column early return remains an empty
array even with a nonzero flag. `F2m_ker` copies its columns first. Native shared
storage between different vector objects has no representation in bigint slots;
valid independently owned matrix columns retain the original mutation and basis
order. The C oracle allocates independent columns before testing the in-place API.

Negative polynomial bits, invalid dimensions/counts and nonzero padding above a
matrix's row count raise `RangeError`. The native square-root kernel assumes a
square and can read beyond its lookup table for nonsquares; these are guarded.
The native remainder kernel requires a nonzero divisor, and Frobenius matrix
construction requires positive degree; the port guards undefined calls rather than
running invalid native memory/shift loops. Defined native division-by-zero errors
retain the `PariError` message. These choices provide deterministic errors outside
native preconditions; they do not change valid mathematical results.

`F2x_valrem` returns `[valuation, remainder]` with both entries bigint, including
native `LONG_MAX` for zero. Factor results are `[packedFactor, multiplicity]` pairs;
DDF pairs contain common irreducible degree. Squarefree arrays are indexed by
multiplicity minus one and retain constant-one placeholders. Direct PARI factoring
of zero returns `[[0n, 1]]`; Sage's public polynomial method rejects zero before the
backend call, as its source requires.


### PARI quotient power count adapter

- Sage/PARI takes signed machine counts for FpXQ_powers/FpXQ_autpowers and an
  unsigned machine exponent for FpXQ_autpow. Negative vector counts are outside
  the native allocation contract. The existing TypeScript signatures use number;
  they now reject negative, fractional or unsafe counts with RangeError instead
  of silently returning an empty/initial result or losing exponent bits.
- Rationale: preserve exact binary exponents within the established number API.
  Trade-off: native unsigned exponents above Number.MAX_SAFE_INTEGER are not
  representable by this signature. This guard affects invalid/unrepresentable
  inputs only. Valid results, including unreduced initial vector entries, are
  compared directly with the bundled PARI library. Numeric trailing zero storage
  is canonicalized before measuring degrees or copying the first power; actual
  coefficient values in that entry remain unreduced. Native reciprocal
  preparation occurs before the zero/unit table shortcuts, so nonunit leading
  coefficients can fail at the generic/word cache thresholds even for count zero.
- FpXQ_autpowers retains the existing 1-indexed vector adapter (slot zero unused).
  Its initial automorphism and FpXQ_auttrace's unit-count pair now canonicalize
  storage without reducing copied coefficients. Both prepare the native generic
  reciprocal before their initial-result shortcuts; FpXQ_autpow prepares it only
  after its count-zero/count-one remainder paths. Direct word composition skips
  the word reciprocal if coefficient reduction makes the input zero. These are
  fidelity repairs with no additional input restrictions or output deviations.
  Binary powering now follows PARI's composition order and reuses native
  Brent–Kung tables. The larger-prime minimal-polynomial/factor routes and
  underlying quotient-reduction algorithms remain under audit.


### PARI matrix product representation

- PARI stores ZM/FpM/Flm matrices as native columns. ZM_mul retains the existing
  TypeScript column arrays with unused slot zero in both dimensions; the new
  FpM_mul/Flm_mul use that same adapter. F2m_mul accepts packed bigint columns and
  a separate row count, consistent with the existing F2m kernel adapter.
- Rationale: preserve the established Buchmann API while sharing the original
  backend dispatch with modular composition. Trade-off: callers must supply
  rectangular, dimension-compatible matrices, just as the native C routines
  require; the unused slots and packed row count are explicit in TypeScript.
- Flm_mul rejects nonpositive or non-word moduli; F2m_mul rejects invalid row
  counts and column bits outside their declared dimensions. FpM_mul now follows
  native signed-modulus and empty/zero-product boundaries; its former blanket
  nonpositive guard rejected defined native behavior (see the signed composition
  and matrix entry below). Modulus one returns zero matrices.
- Valid outputs and algorithm cutoffs match native ZV.c/FpV.c/F2v.c. Integer
  products use size-dependent Winograd thresholds, the native 2x2 kernel cutoff,
  and balanced-input modular reconstruction at dimension 70. Word products use
  the 64-bit thresholds 140/40/70. CRT uses the native prime range and bit bound;
  its balanced merge tree runs synchronously instead of native worker scheduling.
  This execution scheduling changes no mathematical outputs or random state.


### PARI modular composition input contracts

- FpX/Flx polynomial arrays have ascending coefficients. Their precomputed power
  vectors are zero-indexed, matching the existing FpXQ_powers adapter. Native
  table construction, optimal parameter selection, blocked matrix evaluation and
  giant-step folding are preserved; automorphism callers reuse their tables.
  Polynomial and table entries remove numeric trailing zeros before determining
  degrees. After coefficient-matrix multiplication, evaluation prepares the
  native reciprocal even for one block. Canonical zero Q returns before that
  preparation; prepared internal contexts reuse their supplied reducer.
- PARI's FpX vector conversion truncates polynomial coefficients to the modulus
  degree. That behavior is preserved, including zero results for constant moduli.
  Its Flx vector conversion instead writes past the allocated buffer if a copied
  power has excess coefficients. The port guards that undefined native contract
  with RangeError('power polynomial exceeds modulus degree').
- Empty/insufficient power tables and zero moduli are guarded when evaluation
  requires them. A zero evaluation polynomial returns zero before these checks,
  as the original does. Word modulus and safe table-count guards reuse the existing
  adapters. FpXQ_auttrace takes a nonzero unsigned-word bigint exponent, matching
  the generic powering schedule, including its sliding-window boundaries.
- Rationale: avoid native memory corruption or invalid allocations while matching
  all defined results. Trade-off: unsafe C input combinations receive explicit
  errors. Valid native coefficients, table sizes, general trace pairs
  and signed scheduler weights are compared directly with bundled PARI.
- Scalar and polynomial multiplication now reduce signed coefficients canonically,
  as FpX.c does. Multiplication now uses the shared native packing and bounded
  Karatsuba adapters described below. Quotient reduction uses the native schedules
  described below; minimal-polynomial and factor routes remain separate audit work.


### PARI polynomial multiplication adapters

- Native FpX multiplication/squaring delegates word moduli to Flx and larger
  moduli to ZX/GMP. The port now preserves those dependency boundaries and the
  separate square kernels. ZX uses signed Kronecker packing, scalar/linear
  shortcuts and the bundled 64-bit GMP square-tuning tables. Flx uses the native
  small-word threshold 3037000493 and packing cutoffs 30/8 for products and 37/14
  for squares. Under this fixed tuning the packing cutoff precedes Karatsuba.
- Polynomial coefficients remain ascending BigInt arrays. Direct Flx products
  require reduced coefficients and a positive 64-bit modulus, with explicit
  RangeError guards outside that native contract. FpX accepts signed integer
  coefficients and performs the original conversion/reduction. Zero tails are
  normalized as native polynomial construction does.
- JavaScript rejects sufficiently large temporary BigInts even when every result
  coefficient fits. Packed products are therefore bounded at 2^16 temporary bits.
  Above that budget, the port uses the native balanced/unbalanced Karatsuba split
  and the characteristic-two square shortcut until packed leaves fit. Signed
  integer products use the same algebraic split over Z. Earlier recursive splitting
  also limits the cost of large packed multiplications on Bun.
- Rationale: avoid artificial allocation failures without an unbounded quadratic
  fallback, consistent with the existing FLINT/NTL portable product policy.
  Trade-off: native GMP FFT/SS crossover tuning and asymptotic performance above
  the budget are not reproduced; the fallback has Karatsuba complexity. Outputs
  are exact and compared with bundled PARI, including the three previously failing
  temporary-limit cases. Engine limits on individual BigInt results still apply.
  Quotient reduction, half-GCD and minimal polynomials now use the native schedules
  described below; complete factor/extension routes remain open.

### PARI polynomial division adapters

- FpX/Flx division and remainder now preserve the bundled 64-bit GMP thresholds,
  reciprocal basecases, Newton precision schedule and descending Barrett blocks.
  Multiplication delegates to the corresponding PARI product port. Public quotient
  results use `[quotient, remainder]`; reciprocal inverses use ascending arrays,
  truncated at degree(T)-1 as in PARI. Trailing zero inputs are normalized by the
  array boundary, as native polynomial construction does.
- Native FpX reduces smaller dividends, skips inversion for constant remainders,
  and retains an unreduced leading quotient coefficient in the large-modulus monic
  basecase. These outputs and native PariError messages are preserved, including
  the differing Fp_inv (gcd) and Fl_inv (input) error payloads. Native Newton
  inversion preserves its logical `[0]` result at modulus one; the adapter does
  not trim that result to the empty polynomial.
- Direct Flx inputs require reduced coefficients and a positive unsigned 64-bit
  modulus. The same RangeError guards as the multiplication adapter enforce that
  contract. For a zero divisor, the native word remainder basecase reads the zero
  variable header as its leading coefficient: ordinary moduli raise Fl_inv's
  inverse error, while modulus one continues into invalid indexing. The port
  preserves the former error and raises `RangeError('remainder divisor must be
  nonzero')` only for the latter unsafe path, including FpX word conversions.
- Rationale: preserve defined native outputs while preventing undefined native
  memory operations at the JavaScript boundary. Trade-offs: flat arrays do not
  carry native precomputed `[inverse, polynomial]` reduction objects; callers
  currently recompute the inverse. BigInt replaces machine-word overflow/reduction
  instructions. Behavioral impact: the documented guard affects only the unsafe
  modulus-one zero-remainder path; valid-domain results are compared to bundled C.
  Half-GCD, minimal-polynomial and factor routes remain under audit.

### PARI polynomial GCD adapters

- FpX/Flx gcd and extended gcd return PARI's unscaled final remainder and Bézout
  coefficients. GCD, extended-GCD and half-GCD use the original 64-bit GMP cutoffs,
  recursive high/low polynomial splits and seven-product polynomial matrix kernel.
  The former unconditional monic normalization has been removed. Inverse, Hensel
  and distinct-degree callers normalize only where their own native algorithm
  requires it. The split-part constant/linear shortcut preserves unreduced input.
- `FpX_halfgcd_all` / `Flx_halfgcd_all` return `[M, a, b]`; their matrix-only variants
  return `M`. Here M is a zero-indexed two-by-two array of **rows**, with ascending
  polynomial coefficient arrays in each entry. This differs from the native GEN
  column storage and from the existing one-indexed integer-matrix adapters.
  It satisfies `[a,b]^t = M*[x,y]^t`, matching the native transformation.
- Rationale: make the small polynomial transformation explicit without exposing
  PARI output pointers or memory headers. Trade-off: callers must use the stated
  row convention rather than native GEN indexing. Behavioral impact: corresponding
  matrix entries, transformed polynomials and coefficient scales are compared
  exactly with bundled PARI. Direct Flx inputs retain the reduced-word guards.
- The full big-prime extended-GCD basecase requests both cofactors and can raise
  native division-by-zero for x=0, y!=0; the half-GCD route has different zero
  behavior. Both native branches are preserved. Inversion requests only the
  second cofactor in native `(T,x)` order and rescales by the gcd constant.
  The existing compact PariInvError quotient-inverse adapter remains documented
  under Extension Arithmetic and PARI Quotient Kernels; GCD errors themselves
  retain their native class/message payloads.
- Quotient powers preserve unreduced initial bases for exponents beyond +/-1,
  delegate to word arithmetic only below 2^63 (the native is_bigint boundary),
  and use separate square kernels. `gen_pow_i` implements PARI's nonzero exponent
  magnitude schedule, including every arbitrary-precision sliding-window cutoff.
  Its zero-exponent RangeError guards a native helper precondition. Barrett inverse
  contexts, minimal polynomials and full prime-field factors now use the native
  algorithms. Extension factors remain under audit; prime-field roots now use the native splitter.


### PARI minimal-polynomial input contracts

**SageMath/PARI:** FpXQ_minpoly delegates word primes to Flxq_minpoly; both use
Shoup's randomized transposed-power projections, half-GCD reconstruction and
reusable Barrett inverses. They consume PARI's global random state. A positive
prime modulus and positive-degree polynomial modulus are native preconditions;
a word polynomial with degree at least the modulus degree reaches an unchecked
Flx_to_Flv vector copy. Modulus one can loop forever.

**This port:** Both routes now follow those algorithms and consume the exact
native random stream. Word dispatch uses the unsigned 64-bit boundary; only the
large-prime route reduces the input modulo the polynomial modulus first.
Nonpositive/one moduli, constant/zero polynomial moduli and oversized word input
raise explicit RangeError messages. Direct Flx input also uses the existing
positive-word-modulus and reduced-coefficient guards.

**Rationale:** Preserve defined native results, random state and algorithmic
complexity while rejecting native memory-unsafe or nonterminating contracts.
**Trade-offs:** These guards replace undefined native behavior; they do not
validate primality or extend these field algorithms to composite moduli.
**Behavioral impact:** Defined finite-field cases compare both the returned monic
polynomial and full getrand() state against the bundled native executable.
Cached reduction contexts remain internal; callers pass ordinary coefficient
arrays rather than PARI reduction-object GENs.


### PARI distinct-degree factor adapters

**SageMath/PARI:** FpX_ddf converts and normalizes positive-degree input, selects
F2x for modulus two, Flx for other word primes, and FpX otherwise. The latter two
use Shoup's baby-step/giant-step splitting. Word refinement retains the unscaled
GCD factors; the large-prime refinement normalizes them. Native results contain
parallel factor and common-degree vectors. Flx_nbfact_by_degree returns a native
vector and an out-parameter total; its documented input is squarefree.

**This port:** The mirrored kernels use factor/degree pairs, matching F2x_ddf's
existing representation. The existing public FpX_ddf still returns a Map from
common degree to factor. Flx_ddf returns `[factor, degree][]`. Factors retain the
native scale. Count results use `{D, nb}`, with unused `D[0] = 0`; the existing
FpX_nbfact_by_degree convenience delegates word inputs to Flx and uses the FpX
Shoup components for large primes. Native integer division replaces fractional
counts on repeated-factor controls. Such controls do not extend the mathematical
squarefree contract.

**Rationale:** Preserve existing consumer signatures while matching native
outputs, grouping, backend selection and algorithmic complexity.
**Trade-offs:** GEN vectors and out-parameters use TypeScript containers. Moduli
less than two receive an explicit RangeError, replacing an invalid native field
contract. Direct word kernels require a positive word modulus and reduced
coefficients. Modulus-two zero-polynomial DDF retains the existing binary
RangeError adapter; the native routine does not terminate. Defined odd-prime
zero-polynomial errors, including count-vector allocation overflow, are preserved.
**Behavioral impact:** Shared comparisons check exact coefficient scales, order,
counts and exceptions. Factor callers explicitly normalize DDF components where
native factor_Shoup does; full prime-field factors now use native equal-degree
splitting and random state. Prime-field roots now use the native splitter;
extension factors remain open.

### PARI modular square-root adapters

| Aspect | Bundled PARI | TypeScript adapter |
|---|---|---|
| Word nonresidue | `Fl_sqrt` returns `ULONG_MAX` | Returns `null`, like `Fp_sqrt`, avoiding confusion with a field element |
| Word input contract | Requires reduced unsigned-word arguments and a prime modulus | Rejects a nonpositive/oversized modulus or an unreduced argument with `RangeError`; preserves native zero/composite results inside the representable contract |
| Fused exponent | `gen_pow_fold` assumes a nonzero integer and uses its magnitude | Rejects exponent zero with `RangeError`; generic callbacks and returned values follow the existing immutable-value powering convention |
| Prime search for large-modulus shortcuts | `u_forprime_init/next` traverses primes through `ULONG_MAX` | Uses the existing segmented sieve in 2^16-entry blocks through 2^20+2, then the existing arbitrary-precision next-prime search |

Rationale: `null` is the established square-root failure type; explicit guards
replace invalid native machine-word inputs. Prime traversal reuses the audited
sieve/search dependency without allocating a sieve for the entire word range.
Trade-offs: raw C sentinel values and unsafe input behavior are not exposed;
prime iterator storage and switching thresholds are adapted to JavaScript.
Behavioral impact: prime-field roots, callback schedules and the sequence of
candidate primes are preserved. Native composite controls are regression tests,
not an extension of the original prime-modulus precondition. The native word
route, signed-small and Gauss-sum shortcuts, Atkin formula, Cipolla cost test and
Tonelli-Shanks schedule are retained. Square-root searches consume no random state.


### PARI polynomial normalization adapters

**SageMath/PARI:** FpX_normalize preserves zero and monic input, otherwise scales
lower coefficients and explicitly writes a leading one. Flx_normalize attempts
word inversion even for zero, so zero raises the native Fl_inv error for p > 1.
The two inverse routes have different native error payloads.

**This port:** Ascending bigint arrays replace GEN polynomials; trailing zeros are
removed as by the oracle's Polrev construction. The generic zero result is `[]`;
word zero preserves the native error (and the modulus-one empty result). Forced
leading ones at modulus one and native Fp_inv/Fl_inv errors are preserved. Direct
word entries reject moduli outside 1 <= p < 2^64 and unreduced coefficients.

**Rationale:** Use canonical coefficient arrays and reject invalid machine-word
contracts while retaining defined native behavior. **Trade-offs:** Native headers,
variable metadata and object identity are not represented. **Behavioral impact:**
Shared tests compare exact coefficients, empty values and exception classes/messages,
including safe zero, negative-modulus and composite controls. Such controls do not
extend the finite-field preconditions of callers.

### PARI full polynomial factorization adapters

**SageMath/PARI:** FpX_factor and Flx_factor return a GEN matrix of irreducible
factors and multiplicities, discard the leading unit, and consume global PARI
random state. Generic factorization of zero returns the zero factor with exponent
one; direct Flx factorization normalizes first and raises an inverse error.

**This port:** `PolynomialFactor` is `[bigint[], number]`, and full factorization
returns `PolynomialFactor[]`. Factors retain native degree/coefficient order and
multiplicities. Both entry points reject p < 2; direct Flx also requires a positive
64-bit modulus and reduced coefficients. The Galois support helper preserves its
constant-empty and squarefree validation adapter before delegation. Sage's polynomial
wrapper restores its leading unit and Sage ordering; it retains its own zero error.

**Rationale:** Preserve native algorithms, values, random draws and backend routing
without exposing GEN matrices or unsafe field contracts. **Trade-offs:** Container
shape and multiplicity storage differ; no primality validation is added. **Behavioral
impact:** Defined prime-field outputs and full random states match the bundled C
oracle. Word EDF collision restarts and EDF_simple recursion retain native schedules.
Large-prime Sage factorization now delegates to this backend. Extension-field factor
routing remains under audit; the separate prime-field root splitter is now native.

Packed polynomial products use hexadecimal strings when coefficient-block widths
are divisible by four, and binary strings otherwise. This reduces conversion costs
without changing the represented packed integer, recursive product budget, native
algorithm cutoffs or output coefficients. The trade-off is JavaScript string storage
instead of native limbs; aligned/unaligned and signed boundary comparisons cover it.


### PARI polynomial root adapters

**SageMath/PARI:** FpX_roots reduces its input before zero/degree checks and
returns distinct roots. Flx_roots accepts reduced word coefficients. Both use
native small-degree shortcuts, word square/nonsquare cutouts and deterministic
shifts; neither consumes global random state. General word roots use signed-long
vecsmall_sort, while the quadratic shortcut compares unsigned roots. Larger-prime
roots use ordinary integer sorting. These ordering differences are preserved.

**This port:** Ascending bigint coefficient arrays replace GEN polynomials;
root arrays replace GEN columns/word vectors. Trailing zero coefficients are
canonicalized. Roots reject p < 2, and all direct Flx entries require
1 <= p < 2^64 with reduced coefficients. FpX_is_totally_split retains raw-degree
constant/linear returns and the unsigned-degree bound before rejecting remaining
p < 2 inputs with RangeError. The native zero-predicate difference survives:
zero is false below the word limit, but larger prime moduli reach an inverse error.
Counts retain -1 for zero, and direct word counts use the quadratic discriminant
and binary square shortcuts.

**Rationale:** Keep defined prime-field outputs, errors, ordering and algorithm
complexity without exposing GEN metadata or unsafe machine-word inputs.
**Trade-offs:** Explicit guards replace invalid native field contracts, including
negative/nonprime-modulus behavior outside the tested controls. There is no
primality validation added to every entry. **Behavioral impact:** Shared tests
compare native values and exception classes/messages, including selected safe
composite errors. Those controls do not extend the original prime-modulus
precondition. The remaining Galois prime-selection and cyclotomic shortcuts are
separate deviations; matching root splitting alone does not close them.

### PARI extension-polynomial coefficient adapters

**SageMath/PARI:** FpXQX coefficients may be GEN integers or GEN polynomials;
FlxqX uses word coefficient polynomials and F2xqX uses packed binary polynomials.
Coefficient tags affect dispatch and results. In particular, FpXQX_sqr on integer
coefficients returns ZX_sqr without reduction. Generic monic normalization retains
raw lower coefficients, while word/binary normalization inverts even a leading one.
These native behaviors, including constant-polynomial tags, are preserved.

**This port:** Ascending arrays replace GEN polynomial containers. Generic
`ExtensionCoefficient` is bigint or bigint[]; word outer polynomials are bigint[][],
and binary outer polynomials are bigint[] with packed coefficients. GEN variable
metadata and cached reduction objects are not public inputs. Inverse errors render
the inner variable as `y`, matching the canonical native oracle variable order.
Existing standalone quotient-inverse adapters retain their separately documented
error convention; these new normalization errors compare exact classes/messages.

All inputs require a positive-degree extension modulus. Generic/word entries
require p >= 2; word entries additionally require p < 2^64 and reduced coefficients,
including those of T. Negative binary coefficients raise RangeError. Binary products
retain native 2*deg(T)+1 XOR packing, including raw coefficients that overlap
blocks. For n nonzero-storage coefficients and d=deg(T), native allocation has
64*ceil(((n-1)*(2*d+1)+d+1)/64) bits. A coefficient shifted to its block must fit
that allocation; otherwise the port raises RangeError instead of permitting a
native out-of-bounds write. Binary squares/reduction/normalization accept
unreduced coefficients. Generic/word Kronecker degree failures retain native BUG
errors, and mixed generic scalar products retain the native transpose shortcut.
Trailing zero storage is canonicalized. No primality or irreducibility check is added.

**Rationale:** Preserve coefficient-dependent algorithms, asymptotic complexity,
outputs and native errors while omitting GEN memory/variable metadata and unsafe
machine-word inputs. **Trade-offs:** Callers must choose the coefficient container
explicitly; invalid native field contracts gain guards, and arbitrary variable names
are not retained in errors. **Behavioral impact:** Defined field results and selected
safe reducible-modulus/error controls match the bundled C oracle, including unreduced
integer squares. These controls do not extend native field preconditions. The new
arithmetic dependencies do not yet replace Sage's extension-factor backend.

### PARI extension-polynomial division adapters

**SageMath/PARI:** Extension division has separate quotient, remainder, pair and
exact-divisibility output-pointer modes, with distinct basecase/Barrett cutoffs.
Reciprocal inversion uses its own basecase/Newton cutoff. A GEN pair can carry a
precomputed outer reciprocal. Generic division below 2^64 converts to word
polynomials and collapses constant polynomial coefficients to integers; the early
remainder return precedes this conversion. Word/binary basecase reciprocals
reduce coefficient products but retain the raw initial coefficient term; the
generic basecase reduces the complete sum. These native differences are preserved.

**This port:** Dense tagged arrays retain the arithmetic representations described
above. Public divrem returns a pair, div/rem return polynomials, and the private
kernel additionally preserves ONLY_DIVIDES's polynomial-or-null result.
ExtensionReduction stores polynomial and inverse fields; ExtensionModulus accepts
either this object or a plain polynomial. Use the matching get_red constructor
for valid cached data. Existing reduction objects are retained. The inner T
remains a plain polynomial/packed value; public inner-cache metadata is not added.

All entries retain the shared guards: positive-degree T, p >= 2 for generic/word
arithmetic, p < 2^64 and reduced numeric coefficients for word inputs, and
nonnegative packed binary coefficients. Division preserves native raw inner-coefficient behavior. Generic/word Kronecker
packing requires degree below deg(T); binary packing preserves raw XOR overlaps within its native word allocation. Trailing input zero storage is canonicalized, but native
raw quotient coefficients and output tags are preserved. Native inverse and zero
divisor errors use PariError and the canonical inner variable y. Primality and
irreducibility are caller contracts, not additional validation.

**Native binary quirk retained:** F2xqX_divrem_Barrett's remainder-only branch uses
lg(r), rather than the reduced logical length, when its final short block is
already below the divisor degree. Stale high coefficients can survive even over a
valid field: rem(X^100+1, X²) with inner modulus y²+y+1 returns X^100+1.
The port preserves this source behavior, including overwritten intermediate
coefficients. The paired divrem operation follows its separate native branch and
returns remainder one in that example. This is an upstream behavior, not a
mathematically canonical remainder convention.

**Rationale:** Preserve native values, coefficient tags, error timing, reciprocal
reuse and algorithmic complexity while replacing GEN/pointer storage with typed
containers. **Trade-offs:** Cached objects must contain a valid matching reciprocal;
GEN variable metadata and inner reciprocal objects are not accepted. Guards replace
unsafe native field contracts. **Behavioral impact:** Shared native comparisons
cover all output modes, cached and uncached cutoffs, selected safe reducible and
composite controls, raw leading zeros and the binary source quirk. Barrett traversal
retains a single backing array; its linear copying bound is checked independently.
The complete extension GCD/factor backend remains unfinished.

### PARI extension-polynomial GCD adapters

**SageMath/PARI:** GCD returns its unscaled last remainder, and extended GCD
returns the corresponding exact Bézout coefficients. Extension basecases track V
and derive U by quotient division; requesting U with a zero first polynomial
therefore raises native division-by-zero, while the V-only mode can succeed.
Half-GCD uses recursive seven-product matrix multiplication and native thresholds.
Its public matrix operation does not reduce inputs first. These details differ
from some base-field routines and are preserved.

**This port:** Public gcd/extgcd/halfgcd entry points use the same generic, word
and packed-binary coefficient containers as extension arithmetic. Extgcd returns
[gcd,U,V]; ExtensionMatrix stores the native matrix in row-major order. A private
kernel additionally supports the native V-only output mode. Generic primes below
2^64 delegate to word routines and collapse constant polynomial coefficients to
integers. Shared field guards and canonical input storage still apply.

**Rationale:** Preserve source algorithms, complexity, native scaling, coefficient
tags and output-pointer behavior using TypeScript containers. **Trade-offs:** GEN
variable metadata, inner caches and native output pointers are not public inputs.
Unsafe allocation contracts retain documented guards. **Behavioral impact:** Shared
comparisons cover values, matrix entries, both Bézout output modes and exact
exceptions. Half-GCD matrices follow the bundled source's initial
swap/remainder conventions. These dependencies do not complete extension factoring.

The retained binary remainder quirk also makes native gcd(X^100+1,X²) over T=7
cycle: successive remainder states repeat. The GCD port preserves that source
behavior. Permanent comparisons run these cases in isolated processes with a
two-second deadline and pair them with nearby terminating degrees and extgcd
controls, whose paired-division path avoids the cycle. This is an upstream
nontermination case, not a new canonical GCD convention.

### PARI extension quotient adapters

**SageMath/PARI:** Extension quotient operations combine extension products,
remainders and V-only Bézout inversion. Inner and outer reciprocal caches are
prepared at distinct entry points; zero/unit exponents precede cache construction
in several powering routines. Error payloads use native GEN representations.

**This port:** Private extension kernels now accept an optional inner reciprocal
array while retaining the raw modulus separately. A supplied reciprocal is reused
across coefficient reductions and converted along with generic-to-word dispatch.
This replaces GEN cache metadata without recomputing an inner inverse per operation.

Native error display is also retained: payloads shorter than 1,600 characters are
printed in full; larger payloads use the native omitted-object message, including
its newline and GEN type. This follows the initialized native breakloop setting.
Word/binary outer coefficients display the canonical variable header and signed
64-bit Vecsmall words. Generic coefficients use ordinary polynomial formatting.
Native test oracles JSON-encode error text so multiline payloads remain one result
record; decoded messages are compared without normalization.

**Rationale:** Preserve reciprocal reuse, coefficient representations, error timing
and the initialized native display convention while using TypeScript arrays.
**Trade-offs:** Arbitrary GEN variable metadata and GP breakloop configuration are
not exposed. **Behavioral impact:** The error cutoff changes previously overlong
port messages to their native text. Public quotient arithmetic and power APIs retain these representations.
Power tables take a nonnegative number below 2^32-1; unsigned word exponents must
fit 64 bits. These guards replace unsafe native allocation/machine-word inputs.

Generic powering at primes below 2^64 calls ZXX_to_FlxX on S rather than its cached
container variant. A cached scalar-coefficient S is therefore read as one inner
coefficient; the port preserves that defined output quirk. A cached S with GEN
polynomial coefficients would feed invalid types to the native integer reader.
The port guards that conversion with RangeError after cache construction, while
preserving the 0 and ±1 shortcuts. This guard excludes unsafe GEN memory reads;
word powering and generic powering with a plain S retain their normal paths.
The complete extension-factor backend remains unfinished.


## PARI extension composition adapters

- **SageMath/PARI:** extension-polynomial evaluation uses generic Brent–Kung
  blocks or word-extension matrix packing. Matrix inputs are native GEN columns;
  powers, moduli and entries must have valid low-level storage and dimensions.
- **TypeScript:** the six direct/table evaluation APIs use the existing ascending
  coefficient arrays and extension-modulus cache objects. The new FlxqM_mul API
  uses columns with unused column/row slot zero, as the existing PARI matrix APIs
  do. Inner scalar polynomials remain bigint arrays. Native packed multiplication
  delegates to the existing integer-matrix backend. Plain and cached inner
  reductions follow the original construction schedule.
- **Rationale:** arrays and explicit cache objects replace native GEN pointers;
  validation makes native unsafe allocation and division domains deterministic.
- **Trade-offs:** matrix inputs must be rectangular and dimension-compatible;
  p must satisfy 2 <= p < 2^64 and word coefficients must be reduced. Positive
  inner modulus degree and existing generic/binary coefficient guards apply.
  Word evaluation rejects a nonzero Q with an empty/insufficient power table,
  or a zero outer modulus when it reaches matrix conversion. Earlier native
  cache/power errors and zero-polynomial shortcuts keep their source order.
- **Behavioral impact:** accepted values preserve native output tags, outer
  coefficient truncation in word matrix conversion, inner-only block reduction,
  and exact errors. Generic/binary short tables retain the native domain error;
  unsafe word cases raise RangeError. These representation/domain adaptations
  are explicit comparative-oracle preconditions, not claims about native
  behavior outside its valid storage domain. Frobenius, minimal-polynomial and
  full extension-factor dependencies remain under audit.

### PARI coefficient-substitution adapters

Sage/PARI's `FpXY_FpXQ_evalx` / `FpXY_FpXQV_evalx`, corresponding FlxY and
F2xY routines, and binary scalar composition operate on GEN polynomials and
native power vectors. The port uses ascending coefficient arrays, existing
integer-versus-polynomial coefficient tags, and packed nonnegative bigint bits
for F2x. Supplied power tables have no unused slot zero. Variable metadata and
native object aliasing are omitted; returned coefficients retain native values.
This follows the existing polynomial representation, at the cost of accepting
neither native GEN pointers nor their variable identifiers.

Generic integer coefficients are copied without modular reduction; polynomial
coefficients are evaluated in the quotient. Canonical trailing zero storage is
removed without changing those tags. An all-scalar supplied-table call performs no
arithmetic with the field modulus or inner polynomial. No extra positive-degree restriction
is imposed on an unused inner modulus. Direct calls still construct their power
table before visiting coefficients, including empty inputs, and preserve native
reciprocal errors at that stage. Binary scalar composition uses native Brent–Kung
blocks and giant steps rather than a linear evaluation loop.

Word inputs retain the existing reduced-coefficient, positive 64-bit modulus
contract. Out-of-word moduli receive `RangeError`; negative binary polynomial
bits are outside the packed representation, with the existing arithmetic guards
retained. Generic/word supplied tables reuse the previously documented short-
table, zero-modulus and word-column-capacity guards. Native word columns beyond
the modulus degree can overwrite memory, so the oracle executes native power
preparation and the safe coefficient prefix before modeling that guard. Earlier
native inverse errors remain observable. These guards affect invalid native
inputs; defined results and exact native errors are compared directly.

Public signatures take the raw inner polynomial. Internal callers can retain a
prepared reciprocal through the private coefficient context. Ordinary raw calls
still prepare a reciprocal separately for each polynomial coefficient, as the
original implementation does. This preserves cache-dependent errors without
exposing native cache containers in the array API.


### PARI extension automorphism adapters

- **Sage/PARI:** the eight FpXQXQ, FlxqXQ and F2xqXQ automorphism routines take
  native tuples, GEN polynomials and machine-word exponents. Seven entries take
  signed long and pass its unsigned bit pattern to the powering kernel. Word
  auttrace takes unsigned long. The kernel requires nonzero n; zero is undefined.
- **TypeScript:** ascending arrays and existing extension-cache objects replace
  native GEN storage. Tuples omit the native unused slot zero. Binary inner
  polynomials are packed bigint bits. Exponents are bigint restricted to the
  corresponding nonzero 64-bit domain; signed negatives retain the native cast,
  including LONG_MIN. Invalid counts raise RangeError before context preparation.
- **Rationale:** explicit arrays and bigint make native layouts and integer widths
  deterministic without exposing unsafe pointers or undefined zero-bit scans.
- **Trade-offs:** variable identifiers and native aliasing are omitted; returned
  tuples copy their coefficient arrays. Existing extension-field guards require
  p >= 2 and positive-degree T; word inputs require reduced 64-bit coefficients,
  and packed binary inputs must be nonnegative. Existing evaluator guards stop
  insufficient tables and memory-unsafe native word-column copies when reached.
  These are domain adaptations, including for n=1, not native error claims.
- **Behavioral impact:** accepted inputs retain native tuple layouts, additive
  trace versus multiplicative autsum behavior, signed-to-unsigned count casting,
  coefficient tags and raw values. Numeric trailing zero storage is canonical.
  Inner then outer reciprocals are prepared even for n=1, retaining native inverse
  errors. Prepared caches are reused through the native window and Brent–Kung
  schedules. Comparative oracles call the original functions for defined inputs
  and model only the explicit representation/count guards outside that domain.


### PARI extension projection adapters

- **Sage/PARI:** extension random projections, dot products and truncated products
  use GEN/Flx arrays, native variable identifiers and signed-long lengths. Dot
  products defer inner reduction until after summing. Generic truncated products
  pack and truncate before reducing coefficients; word routines reduce the full
  product first. Generic integer-only squares do not reduce modulo p.
- **TypeScript:** the eight APIs use ascending coefficient arrays and existing
  scalar/polynomial tags. Lengths are nonnegative safe integers. Truncated
  arithmetic and dot products inherit the extension-field guards (p >= 2,
  positive-degree T, and reduced word coefficients). Native packing, backend
  multiplication and reduction schedules are retained. Internal callers may
  supply a prepared inner reciprocal without exposing its native container.
- **Rationale:** this follows the established array representation and prevents
  undefined native allocation lengths and invalid word storage.
- **Trade-offs:** native variable metadata, pointer aliasing and negative/native
  oversized allocation domains are omitted. Representation guards can reject
  invalid field inputs even when arithmetic would otherwise be skipped.
- **Behavioral impact:** accepted values preserve tags, exact generic scalar
  squares, delayed dot reduction and the different truncation/error ordering.
  Word dot products are mirrored in Flx.ts, where PARI defines them; other entries
  reside in FpXX.ts or FlxX.ts. Dot products ignore coefficients beyond the shorter
  canonical outer length, except for the common representation validation.

Sampling only needs the degree of T and preserves its raw coefficient values.
Word T coefficients must be unsigned 64-bit values, but do not need reduction
modulo p when unused for arithmetic. The oracle constructs this raw Flx storage
instead of reducing it during conversion, which would alter the degree at p=1.
Zero outer length consumes no randomness and permits zero T or an unused p=0.
Positive outer length rejects zero T because its coefficient length would be -1;
constant T gives zero coefficients without consulting p. Word p must fit unsigned
64-bit storage, and a nonempty scalar sampling request retains the existing
positive-bound guard. Full random state is compared with the native implementation.


### PARI extension minimal-polynomial adapters

- **Sage/PARI:** FpXQXQ_minpoly and FlxqXQ_minpoly use randomized Shoup projection,
  transposed products, half-GCD and composition to find the minimal polynomial
  over the coefficient field. Low-level callers supply valid GEN/Flx polynomials
  and quotient representatives. Reciprocal helpers write directly into a native
  allocation based on the modulus degree.
- **TypeScript:** ascending arrays and existing modulus-cache objects replace
  native storage. The algorithm, random-state consumption, coefficient tags,
  projection order and error phases are preserved. Numeric trailing zeros are
  canonicalized. The adapter requires positive outer modulus degree and
  deg(x) < deg(S), in addition to the existing extension-field/word guards.
- **Rationale:** zero-degree moduli and oversized reciprocal inputs can make the
  native helper write before its allocation. Explicit quotient-representative
  guards prevent that undefined memory behavior without emulating native pointers.
- **Trade-offs:** the adapter also rejects oversized raw representatives for which
  a particular native path might avoid an invalid write; callers must reduce
  them modulo S first. Native variable identifiers and pointer aliasing are
  omitted. These are documented domain restrictions, not native exception claims.
- **Behavioral impact:** supported inputs preserve monic results, repeated-factor
  multiplicity, mixed generic coefficient tags and the complete random stream.
  Nonpositive/constant S and oversized representatives raise RangeError. Native
  field assumptions (prime p, irreducible T) are not rechecked; composite/nonfield
  inputs can still expose native inverse errors or fail to terminate. Comparative
  tests exercise valid coefficient fields and defined early-error paths.

The outer reciprocal is prepared with the original T before the initial power
call prepares its own inner reciprocal. Subsequent operations retain the original
inner context; an explicitly supplied private reciprocal remains reusable. This
preserves the native preparation schedule, including its distinct Fp_inv/Fl_inv
errors. The source's defensive projection restart remains in the algorithm.


### PARI prime quotient cache and inverse errors

- **Sage/PARI:** signed FpXQ powering returns immediately for exponent zero and
  handles +/-1 before dispatch. Other negative exponents invert first, then
  prepare FpX/Flx reciprocal caches before powering. Generic-to-word dispatch uses
  the signed 2^63 boundary. FpXQ_invsafe uses a safe scalar inverse; Flxq_invsafe
  uses Fl_inv, whose nonunit error can precede the quotient-inverse error.
- **TypeScript:** the existing array APIs now share these exact phases. Reciprocals
  are prepared even for exponent two and are reused through the exponent window.
  Inverse failures raise PariError with the native backend name, original or
  converted polynomial payload and the native large-display cutoff.
- **Rationale:** repairing the skipped preparation and incomplete exceptions makes
  the existing deterministic APIs agree with the original implementation.
- **Trade-offs:** arrays still omit native variable identifiers and GEN containers;
  error polynomials use the existing x variable. The public three-argument inverse
  and an explicit undefined fourth argument select FpXQ_inv. The existing explicit
  four-argument convenience selects a Hensel inverse modulo q=p^e, initialized via
  the unsigned-word backend for p < 2^64, as ZpXQ_inv does.
- **Behavioral impact:** calls that previously returned values after skipping a
  failing reciprocal now raise the native inverse error. Error classes and messages
  that previously used PariInvError without payloads now match PARI. Exponent
  shortcuts, valid inverse/power results and the documented Hensel convenience
  remain supported. The older comparative area now checks full diagnostics through
  the bundled C oracle instead of normalizing inverse errors to a common label.


### PARI extension Frobenius adapters

- **Sage/PARI:** prime Frobenius uses the native signed or unsigned polynomial
  power entry. Extension Frobenius computes the p-th powers first, then selects
  automorphism composition or a direct p^deg(T) power using the original integer
  degree/bit-size criterion. Half Frobenius prepares inner then outer caches and
  aggregates conjugates of a^floor(p/2). Its generic entry uses word arithmetic
  for every unsigned-word p, converting both parts of a supplied outer cache.
- **TypeScript:** six mirrored exports reuse the corresponding power, division
  and automorphism kernels. Generic integer/polynomial tags, the characteristic-
  two square shortcut and cache-preparation phases are preserved. Arrays and
  existing reduction objects replace GEN storage; inner caches remain private.
- **Rationale:** retaining the native algorithms and shared caches avoids linear
  enumeration and redundant reciprocal construction. Native variable metadata
  is omitted; standalone Flx inputs correspond to variable zero, while nested
  coefficient polynomials use the established extension representation.
- **Trade-offs:** existing word bounds, reduced-coefficient guards and positive
  extension-degree guards remain. Generic extension inputs must also retain
  positive inner degree after word conversion. A leading multiple of p can
  otherwise produce invalid native extension contexts, stack exhaustion or
  memory-dependent behavior. This guard precedes Frobenius operations; shared
  fixtures apply it before constructing optional native/TS caches. Primality and
  irreducibility remain caller assumptions. Generic full Frobenius inherits
  FpXQXQ_pow's cached-polynomial safety guard, including preceding cache/prime-
  Frobenius errors; half Frobenius uses the safe native recursive cache conversion.
- **Behavioral impact:** supported outputs and native diagnostics match. Invalid
  converted coefficient contexts raise deterministic RangeError instead of
  following undefined native assumptions. Standalone Flx_Frobenius preserves
  constant-modulus behavior and p=1's return-before-reduction shortcut. For p=2,
  half Frobenius preserves PARI's exponent 2^deg(T)-1, rather than an odd-prime
  half exponent. These kernels consume no random state.


### PARI prime-polynomial linear boundaries

- **Sage/PARI:** FpX_red, FpX_add, FpX_sub, FpX_neg and FpX_Fp_mul use native
  integer modular primitives. Nonzero signed moduli produce the corresponding
  nonnegative residues. Cancellation in Fp_add/Fp_sub and zero negation return
  before division. Scalar multiplication tests the original scalar for zero,
  then reduces individual products; an empty polynomial never reaches reduction.
- **TypeScript:** the existing exports now preserve these branches, reduce with
  the absolute modulus, and raise PariError with the exact dvmdii diagnostic
  when division by zero is reached. Trailing numeric zeros are removed at the
  array boundary before arithmetic, matching construction of the native GEN.
- **Rationale:** the old signed remainder adjustment added a negative modulus,
  producing negative values outside the native residue range. Unconditional
  remainder calls also missed successful early returns and exposed a JavaScript
  RangeError instead of the native exception.
- **Trade-offs:** arrays still omit GEN variable/pointer metadata. Negative and
  zero modulus comparisons describe the observed low-level integer behavior;
  they do not make these values valid finite-field characteristics. Native field
  algorithms retain their existing positive-prime assumptions and adapter guards.
- **Behavioral impact:** negative-modulus outputs and reached zero-modulus errors
  now match. Zero/cancellation shortcuts return the native zero polynomial.
  Positive-modulus numeric results remain unchanged. Operations remain linear
  in the number of polynomial coefficients and preserve input arrays.


### PARI extension root-count and derivative adapters

- **Sage/PARI:** extension split parts reduce f, compute X^q-X modulo f and take
  an unscaled GCD. Root counts return its degree, including -1 for zero. Native
  FpXQX split parts convert to the word backend before checking degree whenever
  p occupies one limb. Larger generic p preserves raw degree-zero/one input.
  Squarefree predicates compute the outer derivative before GCD and test degree
  zero. Derivatives preserve generic integer/polynomial coefficient tags.
- **TypeScript:** nine mirrored exports reuse native Frobenius and recursive GCD
  dependencies. Generic split parts use the magnitude of a one-limb modulus when
  converting, as the GEN limb dispatch does. Degree shortcuts preserve unused T
  and raw coefficient storage. Root counts return number, and native 0/1
  squarefree results map to boolean. No field enumeration or primality test is
  introduced. Derivatives retain their original scalar/word multiplication steps.
- **Rationale:** using the original algorithms preserves both complexity and the
  unscaled split polynomial. Distinguishing coefficient tags also preserves
  zero-modulus error phases: zero integers skip reduction, while FpX_mulu reaches
  umodui before inspecting even a zero polynomial coefficient.
- **Trade-offs:** arrays omit GEN variables and pointer aliasing. The existing
  reduced-word/positive-word bounds apply before direct word operations. Binary
  bits must be nonnegative. For nonlinear split parts and squarefree GCDs, the
  existing extension-modulus/field guards remain; some native low-degree shortcuts
  accept an unused constant or zero T. With T=null, FqX_nbroots requires integer
  coefficients and rejects polynomial-tagged entries with an explicit TypeError,
  avoiding undefined native integer-kernel reads. Its prime path delegates to
  the existing FpX_nbroots adapter, including that adapter's limits.
- **Behavioral impact:** supported counts, repeated-factor detection, mixed tags,
  unscaled GCDs and defined native errors match. Raw constant multiples of p can
  become zero before the word split shortcut (-1 roots), while the large generic
  shortcut keeps their original degree (zero roots). These are low-level source
  semantics, not a claim that invalid moduli define finite fields.

### PARI signed scalar inverse and division residues

- **Sage/PARI:** invmod computes inverses modulo the magnitude of a nonzero
  modulus, returns a nonnegative GCD on failure, and rejects modulus zero with
  abs(a) as the inverse-error payload. modii returns zero immediately for a zero
  numerator, even at modulus zero; other zero-modulus reductions raise dvmdii.
- **TypeScript:** the shared polynomial division helper now preserves these
  signed-modulus, zero and error-payload branches. Generic Fp_inv and word Fl_inv
  retain distinct payload conventions. The ffinit scalar reducer shares the
  zero-numerator shortcut. Existing public root-count callers inherit the repair.
- **Rationale:** using a signed Euclidean modulus had produced negative residues
  or a negative GCD, incorrectly rejecting invertible values. Unconditional
  remainder at zero also exposed JavaScript exceptions instead of native results.
- **Trade-offs:** BigInt and the existing Euclidean implementation replace native
  GMP limb storage. Variable metadata and native allocation details remain absent;
  huge diagnostic payloads retain the existing native display cutoff.
- **Behavioral impact:** signed inverses and residues now agree with the native
  kernel. All inverses modulo zero fail with Fp_inv and abs(a), including units;
  zero reduction itself returns zero. Ordinary positive-modulus results remain
  unchanged. These low-level boundary results do not relax field preconditions.


### PARI signed polynomial word dispatch

- **Sage/PARI:** FpX.c tests `lgefint(p) == 3` and reads the unsigned magnitude
  limb in `to_Flx`/`to_Flxq`. Multiplication, squaring, division, GCD and half-GCD
  therefore select Flx for either sign when 0 < abs(p) < 2^64. Basecase division
  performs zero-divisor, shorter-dividend and constant-divisor shortcuts first.
- **TypeScript:** the corresponding existing FpX exports now use this magnitude
  test and pass the positive word to Flx. Coefficients are reduced before word
  calls. Generic and constant shortcut paths retain the original signed p.
- **Rationale:** a positive-only test selected a different algorithm, missed word
  coefficient conversion and reported generic inverse errors on negative words.
- **Trade-offs:** negative characteristics are not valid finite fields. These
  comparisons cover deterministic low-level native integer/modular boundaries;
  other algorithms keep their documented domain guards. Word callers still use
  the existing explicit contracts where native input would be unsafe.
- **Behavioral impact:** signed-word division, GCD and half-GCD results and errors
  now match; 1,902 baseline records were wrong. Multiplication and squaring had
  matching values but now delegate to the native word backend too. Large-modulus
  and positive controls match. No exported signatures or complexity are changed.
  FpX_invBarrett retains its generic reciprocal schedule and delegates its product
  steps according to the same corrected multiplication cutoff.


### PARI signed composition and matrix boundaries

- **Sage/PARI:** FpXQ_powers switches to the word backend only for counts greater
  than two and a one-limb modulus magnitude. FpX_FpXQ_eval uses the same magnitude
  conversion after its raw-zero shortcut. Cached powers convert both the modulus
  polynomial and any supplied reciprocal. FpM_mul returns empty dimensions before
  using p, dispatches abs(p)=2 to F2m and other word magnitudes to Flm, and otherwise
  reduces the integer matrix product with modii. At p=0, zero product coefficients
  return zero; a reached nonzero coefficient raises the native dvmdii error.
- **TypeScript:** these existing APIs and their trace/automorphism contexts now
  follow those phases and preserve the magnitude cutoff. Signed matrix reduction
  uses the shared exact scalar residue primitive. Counts zero through two retain
  copied raw initial power coefficients; longer word tables convert them first.
- **Rationale:** positive-only word tests missed initial coefficient conversion.
  The matrix guard rejected defined low-level native results and blocked generic
  table composition. Cached contexts also passed negative/oversized values into
  word routines. The old matrix guard is removed based on native source and tests.
- **Trade-offs:** negative and zero moduli are still not finite-field parameters.
  This support describes the observed deterministic low-level arithmetic paths.
  Rectangular compatible matrices, safe counts and native word-vector packing
  limits remain required. Minimal-polynomial APIs retain their separate documented
  domain guards; this change does not alter those preconditions.
- **Behavioral impact:** repaired 426 power-table mismatches, 12,320 composition
  guard rejections and 891 direct matrix guard discrepancies in initial probes.
  A further 180 cached trace failures were repaired. The 9,120 automorphism alias
  controls already produced native values after the dependency repairs; their
  direct composition now also delegates using native signed-word routing. The
  two old generic empty-matrix guard expectations now execute native C instead.
  No exported signature changes. Positive controls retain matching results.


### PARI polynomial observation boundaries

- **Sage/PARI:** FpX_deriv reduces the integer derivative. FpX_eval uses sparse
  Horner evaluation with Fp_powu for zero gaps; constants and zero evaluation
  arguments return before unnecessary arithmetic. FpX_center compares coefficient
  and half-modulus magnitudes without first reducing coefficients, and preserves
  its output storage. FpX_div_by_X_x computes a remainder only when its output
  pointer is supplied; the existing TypeScript export requests only a quotient.
- **TypeScript:** the four existing exports now preserve those operations and
  shortcut phases. Input trailing zeros are removed as at native GEN construction.
  Centering deliberately does not remove zeros created in the result. Evaluation
  uses the shared logarithmic gen_powu_i schedule and exact scalar residues for
  sparse gaps, including the raw-square branch for generic exponent two.
- **Rationale:** pre-reducing the evaluation argument and computing an unused
  linear remainder caused premature failures. A direct JavaScript remainder in
  differentiation exposed RangeError. Signed remainder adjustment produced wrong
  values, and centering incorrectly performed modular reduction and normalization.
- **Trade-offs:** native centering assumes reduced nonnegative coefficients and a
  half-modulus threshold; signed/unreduced tests record the defined low-level C
  output outside ordinary finite-field inputs. Nonpositive moduli are not finite
  fields. Generic scalar powering uses BigInt residue callbacks rather than GEN
  Montgomery storage; its logarithmic schedule, outputs and reached errors match.
  The separate Galois helper still reduces before centering its integer dot product.
- **Behavioral impact:** 3,379 baseline records are repaired: 19 derivative errors,
  506 evaluation results/errors, 2,314 centering values/storage/errors and 540
  linear-quotient values/errors. All squarefree controls matched before the repair.
  Valid positive-modulus controls retain their results. No exported signatures
  change, and centering/derivation/linear division remain linear in stored degree.


### PARI exported scalar arithmetic boundaries

- **Sage/PARI:** scalar Fp operations preserve distinct modii/remii shortcuts.
  modii(0,0) returns zero, while remii(0,0) raises dvmdii; Fp_sqr uses the latter.
  Fp_inv uses invmod with native GCD error payloads. Fp_div dispatches one-limb
  divisors through Fp_divu: word moduli return a reduced zero numerator before
  inversion, while larger moduli invert p modulo the divisor then divide the
  corrected integer exactly. Fp_center compares magnitudes without reducing;
  Fp_double performs one conditional subtraction without general reduction.
- **TypeScript:** the existing exports now share exact residue/inverse primitives,
  preserve reached errors and signed-word division dispatch, and implement the
  native exact-division branch. Fp_sqr retains its separate zero-divisor check.
  Addition, subtraction, negation and multiply-add keep their native early returns.
- **Rationale:** signed products returned negative remainders, a negative Euclidean
  GCD produced wrong inverses, and generic Error/RangeError replaced native errors.
  Uniform inverse-then-multiply division missed native shortcuts and diagnostics.
  General modular addition also did not match native single-subtraction doubling.
- **Trade-offs:** centering/doubling normally assume reduced field coefficients;
  raw and nonpositive-modulus cases describe observed low-level native behavior.
  Word-divisor division can return a signed integer for a negative large modulus;
  it is deliberately not normalized afterward. Fp_mulu retains its number-valued
  unsigned-word argument and its existing representation limits. These cases do
  not define finite fields of nonpositive characteristic.
- **Behavioral impact:** 6,460 baseline scalar records are repaired, including
  Fp_inv(-1,17)=16 and Fp_mul(-1,1,17)=16. Native error classes, full payloads,
  cancellation, zero/unit and word-boundary phases are preserved. All halve and
  equality controls already matched. No exported signatures change. General
  scalar powering and higher field algorithms remain under separate audit.


### PARI scalar exponentiation boundaries

- **Sage/PARI:** Fp_pow returns zero for exponent zero when the modulus divides
  the base (including zero base and zero modulus). Nonzero word exponents and
  word moduli use Fl_inv/Fl_powu; other negative exponents use Fp_inv. Larger
  exponents use native windows or fused base-two powering, selecting Montgomery,
  modified Barrett or signed remainder reduction at the bundled GMP cutoffs.
- **TypeScript:** mirrors these phases with BigInt arithmetic and the shared
  gen_pow kernels. Montgomery retains magnitude-limb arithmetic; the final
  subtraction and the large-exponent sign adjustment still use signed N.
- **Rationale:** a uniform square-and-multiply loop lost zero shortcuts, inverse
  error classes/payloads and observable signed behavior. Preserving the reduction
  algorithms also preserves the original asymptotic strategy for large operands.
- **Trade-offs:** native word multiplication/preinverse storage is represented
  with exact BigInt products/remainders. Thresholds follow bundled 64-bit GMP
  tuning (127 GEN words for Barrett, 17 for Montgomery), not host-dependent tuning.
  For an invalid negative modulus, _mul2_montred may repeatedly subtract a
  negative number until the PARI stack overflows. We detect that exact nonprogress
  condition and raise RangeError instead. Comparative adapters normalize only a
  reached native negative-modulus stack-overflow error to that checked boundary.
  Defined negative-modulus shortcuts and results remain preserved; this does not
  define a finite field of negative characteristic.
- **Behavioral impact:** 479 baseline records are repaired across 12,128 exact
  native comparisons. Nine additional comparisons include five native stack
  exhaustion boundaries and four defined controls. The Fp_pow signature and
  package exports are unchanged; no other scalar powering API is added.


### PARI scalar predicates and order boundaries

- **Sage/PARI:** Fp_issquare checks characteristic magnitude two or a Kronecker
  symbol other than -1, without reducing first. bezout(0,0) chooses both
  coefficients zero. Fp_order at word moduli uses the magnitude of a nonzero
  word order, otherwise defaults to |p|-1. Generic order requires a positive
  bound before identity checks and keeps the raw base. znorder preserves Mod
  construction errors and gdisplay-formatted coprimality errors.
- **TypeScript:** matches these branches. Word factorization delegates to factoru;
  generic bounds use Z_factor and the recursive factor-splitting schedule.
  znorder's separate bigint arguments model PARI gmodulo construction followed
  by znorder, retaining reduction and sign normalization before coprimality.
- **Rationale:** Euler's criterion substituted a different predicate at composites
  and zero/negative moduli. Uniform bound handling lost word defaults and generic
  validation. Early reduction of the generic base changed its identity checks.
  Generic Error replaced the original PariError class and complete payload.
- **Trade-offs:** arbitrary GEN/factored order inputs still have no corresponding
  exported TypeScript parameter. Word bounds and unreduced/invalid inputs retain
  observed low-level results, which are not promises of group order without the
  native a^N=1 precondition. factoru retains its existing BigInt-backed factorizer.
  This extends the previously documented valid-input znorder adapter to the
  scalar invalid-input boundaries tested here; host resource limits remain distinct.
- **Behavioral impact:** 5,872 baseline records repaired across 28,782 permanent
  native comparisons. Kronecker and plain GCD controls already matched. Native
  errors include the 1,600-character gdisplay cutoff and 2,100-digit moduli. All
  exported signatures remain unchanged.


### PARI rational Galois polynomial boundaries

- **Sage/PARI:** rational polynomial construction divides the integer polynomial
  by its denominator; Q_remove_denom gives a positive common denominator, with
  denominator one for zero. RgX_to_FpX converts each canonical rational coefficient
  in order through Rg_to_Fl/Rg_to_Fp, canceling its fraction before inversion and
  returning zero before inversion when its reduced numerator is zero. Galois
  vectopol centers signed remainders with centermodii, retaining negative half
  ties and comparing the magnitude of a supplied threshold.
- **TypeScript:** QPoly helpers now model that canonical polynomial. They trim
  trailing zeros, normalize signs/content, reject denominator zero, and convert
  coefficients through shared native residue/inverse primitives. permtopol uses
  centermodii's signed remainder and exact native permutation/error phases.
- **Rationale:** the old helpers retained negative denominators and noncanonical
  zeros; one global modular inverse failed even when coefficient cancellation
  made the result defined. Euclidean centering changed signed half ties and
  ignored the original magnitude comparison. Generic errors lost native payloads.
- **Trade-offs:** QPoly remains a TypeScript common-denominator representation,
  not a GEN. Its canonical model is gdiv(integer polynomial, denominator), followed
  by Q_remove_denom or coefficient extraction. This defines a zero-denominator
  error even for an empty numerator. Original vectopol instead divides individual
  coefficients: an empty vector performs no divisions, while a zero denominator
  reaches gdiv for numerator one and dvmdii otherwise. These distinct phases are
  preserved. Invalid memory-index permutations/matrix shapes remain outside the
  native oracle; valid permutations and checked length mismatches are compared.
- **Behavioral impact:** 3,771 preserved records repaired across 7,003 permanent
  comparisons. The oracle compiles the original galconj.c, including its static
  permtopol/vectopol, and uses native rational construction/conversion. Canonical
  valid Galois outputs remain unchanged. No exported signatures or names change.


### PARI Galois integer helper boundaries

- **Sage/PARI:** word GCD/LCM use integer magnitudes and unsigned 64-bit arithmetic;
  LCM silently wraps on overflow. factoru uses the existing PARI factorizer and
  represents zero as 0^1. Integer logarithms use word/base-two shortcuts or binary
  splitting; the checked logint0 entry point validates the base before the operand.
  u_forprime uses a cached table, segmented sieving, then word prime successors.
- **TypeScript:** the existing number-valued helper facade converts integer inputs
  as native itou magnitudes, checks word overflow, and rounds native word outputs
  back to number. LCM wraps before that rounding; factoru_small delegates to factoru.
  The Galois eulerphiu facade also converts word magnitudes, delegates to the
  exact arith2 word-factorization backend, and preserves PARI's value 2 at zero.
  The logint wrapper calls the source-based ispower implementation. Forprime keeps
  BigInt state, a cached table and segmented/large-successor paths, returning zero
  after word exhaustion. No package-root API signature changes.
- **Rationale:** floating multiplication lost word overflow, trial division took
  too long on large word primes, and number increments stopped advancing beyond
  the exact-number boundary. Unchecked logarithm bases could loop indefinitely.
  The native checked logarithm entry point supplies deterministic domain errors.
- **Trade-offs:** the Number facade cannot distinguish all native words. Tests
  project native words through binary64 and compare the exact represented integer;
  inputs must be integer numbers, and native itou rejects magnitudes >= 2^64.
  Prime state remains exact even when successive public results round identically.
  The iterator fixes the native initialized prime-table limit at 500000; it preserves
  that initialization's observed transition which skips 499979 when starting
  exactly there. Process-global changes to the native prime table are not modeled.
  Segment allocation uses the existing TypeScript sieve, not native CPU-cache tuning.
  The earlier prime-successor-cache deviation still describes standalone nextprime/
  precprime; this iterator additionally preserves its selected cached-table behavior.
- **Behavioral impact:** 242 preserved baseline failures, including 27 two-second
  time-limit cases, are repaired. Of those limits, 17 were prime iteration, nine
  invalid logarithm bases and one slow trial factorization. All 1291 permanent
  native comparisons match, including 30 additional sieve/log controls. Invalid
  integer-log inputs now use logint0's class, messages and validation order rather
  than invoking logintall outside its documented preconditions.


### PARI integral-basis denominator adapters

- **Sage/PARI:** ZX negation and scalar products return normalized polynomials;
  nonzero constants are squarefree. Squarefreeness deflates polynomial exponents
  before a modular integer GCD. indexpartial uses strict partial factorization,
  refines small prime powers through Sylvester column echelon reduction modulo
  p^m, and retains a composite perfect-power cofactor with a rounded-up exponent.
- **TypeScript:** those outputs and algorithms now follow the bundled source.
  The strict-factor helper specializes absZ_factor_limit_strict(n,0,&U) to the
  initialized factor limit 500000 and no user-added special primes. Modular GCD
  uses content/valuation extraction, primes above 2^63, finite-field half-GCD,
  centered CRT and exact divisibility certification. P-adic routines preserve
  word versus generic column arithmetic and native precision doubling.
- **Rationale:** the prior denominator bound omitted native refinement and fully
  factored unresolved cofactors, changing deterministic results. The former
  discriminant test also rejected every nonzero constant as non-squarefree.
- **Trade-offs:** process-global native factor-limit/special-prime mutations are
  not modeled. Modular GCD processes prime images sequentially instead of native
  worker batches; its reconstruction and certification remain exact. Native
  Sylvester allocation requires positive-degree input; the port raises RangeError
  after the remainder step for invalid constant input that would reach unsafe
  allocation. Prime-power and matrix helpers otherwise retain native preconditions
  (prime p >= 2, positive power modulus, rectangular column matrices).
- **Behavioral impact:** 855 preserved baseline failures are repaired and 8161
  shared native comparisons match. Zero discriminants now return zero. For x^4-2,
  indexpartial returns 8 instead of 32. The former documented indexpartial
  approximation is closed; no existing package-root signature changed.


### PARI Hensel lifting adapters

- **Sage/PARI:** scalar root lifting updates the root and inverse derivative using
  quadratic_prec_mask. Factor lifting normalizes the input polynomial (preserving
  already monic representatives), builds a stable tree by factor degree and doubles
  precision. Bezout coefficients propagate through that tree. Complete root lists
  use factor lifting. Quotient-root lifting preserves S at precision one, then uses
  shared evaluation powers and a lifted reciprocal derivative.
- **TypeScript:** galconj delegates to the source-based Zp dependency. Its existing
  one-indexed factor/root-vector facade remains; dependency vectors are zero-indexed.
  Scalar polynomial evaluation/derivatives and split-part operations now live in
  their native dependency modules, retaining their prior algorithms and outputs.
- **Rationale:** the old linear-precision two-factor loop and repeated root/inverse
  calculations bypassed the native schedules. Missing normalization broke nonmonic
  factors and Bezout coefficients; single-factor reduction changed raw monic
  representatives. Inverse errors came from different operations/precisions.
- **Trade-offs:** native exact polynomial divisions are unchecked and require an
  actual initial root/factorization. The port raises RangeError when that divisibility
  precondition fails. For multiple factors, native bezout_lift_fact at precision one
  reads an uninitialized tree structure; the port raises RangeError at that point.
  No package-root signature changed. The scalar lift requires p > 1, e >= 1 and
  a simple root; full factor lifting requires pairwise-coprime monic modular factors.
- **Behavioral impact:** 2234 preserved failing records within defined native
  execution are repaired. The shared suite has 13016 direct native comparisons and
  42 explicitly labeled guard comparisons. Those 42 evaluate the native initial
  inverse/evaluation phases and stop before an inexact unchecked division; they do
  not claim PARI itself raises RangeError. Original unchecked outputs are retained
  only as diagnostic evidence, excluded from the bug count. All comparisons match.


### PARI Vandermonde interpolation adapters

- **Sage/PARI:** inverse Vandermonde construction uses the native degree-one/two
  product-tree scheme, a remainder tree for evaluating the derivative, Montgomery
  batch inversion and only then multiplication by the supplied denominator.
  Errors therefore depend on the product of derivative values and the phase
  reached before inversion. FpV_inv assumes its input vector is nonempty.
- **TypeScript:** galconj delegates to the source-based FpX dependency and adapts
  its native column-major matrix to the existing one-indexed row-major facade.
  Batch inverses use one scalar inverse. The shared producttree_scheme also
  supplies the existing balanced-product helper without changing callback order.
- **Rationale:** separate inverses produced different classes/messages and used
  different failure phases; denominator reduction could fail prematurely. The
  native trees and batch inverse preserve the algorithm and its operation order.
- **Trade-offs:** the existing galconj empty-matrix adapter remains: [0n] maps to
  [[0n]], and [] maps to []. Native FpV_inv/Vandermonde access a first element and
  require nonempty vectors; direct dependency calls reject empty vectors with
  RangeError. These empty adapters are tested separately from native comparisons.
- **Behavioral impact:** 1792 preserved error discrepancies are repaired. All
  4227 shared native comparisons match, including 240 valid-matrix controls that
  also matched the previous implementation, direct batch-inverse cases and tree
  sizes through 4097. No existing package-root signature changed.


### PARI permutation and word-vector adapters

- **Sage/PARI:** permutation powering takes an unsigned word; cycle powering takes
  a signed word. The checked integer conversions use itou magnitudes and itos's
  signed bound (including its rejection of -2^63). Permutation order repeatedly
  applies unsigned-word LCM; group_order calls zv_prod, which explicitly assumes
  that signed products do not overflow. vecsmall_lexcmp returns -1/0/1, and uniq
  uses counting when all entries are nonnegative and below the vector length.
- **TypeScript:** checked exponent conversions precede even empty-vector shortcuts.
  Permutation-order intermediates stay exact and wrap as unsigned 64-bit words;
  only the final result converts to Number. group_order delegates to a zv_prod
  dependency with exact intermediates and a single final Number conversion.
  Prefix comparisons return a sign, and uniq follows the native counting cutoff.
- **Rationale:** negative unsigned exponents produced malformed permutations;
  missing conversion checks changed errors. Intermediate Number rounding changed
  valid word-order outputs. Prefix comparisons returned a length difference.
- **Trade-offs:** existing Number-valued interfaces remain and cannot distinguish
  every native word. Native transcripts project final words through binary64 and
  compare the exact represented integer. zv_prod rejects a prefix product outside
  the signed-word range with RangeError; the native oracle certifies the same
  precondition using arbitrary-precision products before calling the unchecked
  routine. These eight guard records do not claim PARI itself raises RangeError.
  Comparison sorting uses JavaScript's sorter with the same O(n log n) bound;
  the native O(n) counting path and its cutoff are preserved. No root API changed.
- **Behavioral impact:** 1926 defined-execution baseline failures are repaired.
  All 10704 permanent comparisons match: 10696 direct native records and eight
  overflow-guard cases. Five additional valid signed-product failures demonstrate
  intermediate rounding independently of unsupported native overflow observations.


### PARI unit-subgroup dependency boundaries

- **Sage/PARI:** listznstarelts delegates to flag-zero znstar, exact-index
  subgrouplist, modular HNF and znstar_hnf_elts. General subgrouplist also accepts
  a scalar index bound or an omitted bound.
- **TypeScript:** the Galois helper now follows those dependencies and their
  subgroup/generator ordering. Exact-index, scalar-bound and unbounded modes are
  implemented for cyclic-factor vectors. Scalar bounds delegate to strict partial
  factorization before the same Birkhoff and HNF engine.
- **Rationale:** closure-of-subsets enumeration returned different ordering and
  used a different algorithm. Porting the dependency route restores the native
  behavior. The default and explicit strict-factor limits share one implementation.
- **Trade-offs:** these new internal modules use zero-indexed arrays and bigint
  arithmetic, with the old one-indexed Number facade retained. They are not added
  to root exports. Modular HNF's coefficient reductions are eager, preserving
  its algorithm and canonical outputs without copying allocation thresholds.
- **Behavioral impact:** 555 existing-port mismatches are repaired: 543 subgroup
  ordering differences and 12 incorrect selections from negative order arguments.
  The 3505 permanent native records cover the implemented domain. PARI's packed
  multiword subgroup-cell behavior is preserved, including its removal of low
  zero limbs; the comparisons do not silently replace that native output with
  the mathematically expected matrix. A further 2781 comparisons cover all three
  bound modes, invalid bounds and strict partial factorization. GEN container
  variants (diagonal-matrix or znstar-structure inputs) are represented by their
  cyclic-factor vector in this internal interface.


### PARI strict partial-factorization state

- **Sage/PARI:** strict partial factorization returns known factors and an output
  parameter for an unresolved base/exponent. The configured factor limit, prime
  table, primorial buckets and special-prime list can affect that partition.
- **TypeScript:** ifactor1.ts returns `[factors, unresolvedOrNull]`. Limit zero
  uses 500000, the initialized table ends at 499979, and primorial buckets end at
  the same native primes. No user-added special primes are configured.
- **Rationale:** subgroup bounds need the native known-factor partition, not a
  full factorization or a list filtered only by prime size. Word-sized integers
  use restricted 2/3/5/7 power detection and primality testing; larger integers
  use the native general-power and conditional primality branches.
- **Trade-offs:** dynamic PARI session settings and special primes have no setter
  in this port. The shared word-prime table initializes lazily to support the
  dependency cycle without changing iterator output. Prime-root candidates are
  certified with the existing exact integer root dependency.
- **Behavioral impact:** 2781 native records match for the initialized state,
  including signs, overflow errors, cache edges and integers through 4096 bits.
  An unresolved base is not necessarily composite; native multiword branches may
  leave a prime base unresolved at a small limit. The default helper delegates
  to this same engine.


## Restored perfect-power signs and failure sentinel — 2026-09-12

The former `Z_isanypower` implementation exposed the magnitude-only internal
search under the public PARI name. Native `Z_isanypower` removes powers of two
from a negative input's exponent and returns a negative base; if no odd exponent
remains, it reports failure. The port now follows that wrapper and retains the
original signed input in its `[0, x]` failure tuple. The returned-base-on-failure
adapter is explicit because native PARI leaves its output pointer untouched.
The source-subpath `is_pth_power` failure exponent is also restored from one to
native zero; both internal callers now consume that sentinel.

Rationale: these are externally observable results and must match the original.
Trade-off: consumers relying on the former magnitude-only or exponent-one result
must adopt the corrected behavior. Behavioral impact: negative perfect-power
results and helper masks/sentinels change; positive final powers and factorization
remain unchanged in the comparative controls. The existing exact-root algorithm
substitution remains documented above.


## Restored native ECM stages — 2026-09-12

The former serial-curve, binary-ladder and additive-continuation substitutions
changed which factor was found and whether a round succeeded. These were
behavioral differences, despite the earlier equivalence claim. They are removed:
the port now follows native batched inverse failures, PRAC transformations,
helix/baby-step tables, giant-step gcd timing and persistent pointer aliases.
The seed counter uses BigInt to preserve increments beyond 2^53.

Rationale: direct deterministic stage outputs must agree, and native batching
and continuation also preserve the original complexity. Trade-off: the internal
port retains PARI's detailed storage and control flow, increasing implementation
size. Behavioral impact: native factor selection and round failures are restored;
final complete factorizations remain sorted. The existing configured insisting
round limit remains a documented deviation. Native-reference tests implement that
limit as a prefix of the original driver; they do not change its operations.

The ECM prime table is initialized through 500,000, as in the other factorization
ports and native oracle. Above its last entry, the B2 wheel uses native base-two
strong probable-prime tests. PRAC retains the original binary64 golden-ratio
split solely to select an addition chain; curve arithmetic and inversions remain
exact BigInt. Its square-root bound is computed exactly.


## Restored signed factorization-helper behavior — 2026-09-12

`is_kth_power` now preserves the native sqrtnr unsupported-negative-input error
when its modular filters succeed. Negative inputs rejected by a filter still
return null. The exact-root helper also checks the bit-length bound before
Newton iteration, so tiny inputs with huge exponents do not allocate a huge
power merely to establish that the integer root is one.

MPQS's implemented class-group relation combiner now inverts modulo the magnitude
of a signed discriminant, matching native invmod. It previously returned -1
instead of a relation. The class-group entry points remain unported; this repair
covers their existing shared relation helper. The MPQS square-root duplicate
also now delegates to the native canonical-root backend as recorded above.

Rationale: errors, factor-base representatives and relation vectors are observable
results. Trade-off: consumers of the internal helpers see the corrected error and
representatives. Behavioral impact: exact errors/relations/roots now match the
native comparisons; the large-exponent boundary returns null without a BigInt
allocation error. The existing exact-root algorithm substitution and MPQS sparse
kernel substitution remain documented separately.


## PARI sparse binary kernel representation

SageMath's delegated PARI sparse kernel uses a GEN vector of one-based sparse
row-index vectors and returns packed binary columns. The port uses zero-based
outer arrays containing one-based row indices. F2Ms_colelim returns one-based
original column indices; F2Ms_ker returns BigInt bitsets with bit i representing
column i of the zero-based outer array. MPQS adapts those results to its existing
one-based Uint32Array bitsets, preserving the unused bit zero.

Rationale: preserve native matrix and word operations while retaining the
existing TypeScript packed-column representation. Trade-off: callers must use
the documented array/index representation rather than native GEN headers.
Behavioral impact: no difference on valid sparse matrices. Above 640 rows the
backend now consumes the native random stream, retries failed Lanczos starts,
and returns the same ordered kernel vectors, which need not span the entire
kernel. The original dense substitute returned a full basis and did not consume
randomness. Its claims of identical outputs and unreachable matrix sizes were
incorrect. The 766 comparisons cover results, column elimination and complete
final RNG state, including native retry cases.


## Restored MPQS symbol and debug relation behavior — 2026-09-12

The internal MPQS krouu/kroiu facades now delegate to the exact Kronecker backend.
The old signed remainder mishandled negative numerators, and Number bitwise shifts
lost high word bits or failed to terminate when the low 32 bits of a positive
power-of-two denominator were zero. The signatures still use Number for native
word inputs; values must be exactly represented integers. The signed numerator
of kroiu remains BigInt. The native helpers require a positive denominator.

The debug relation reconstruction helper now delegates modular negation,
multiplication and signed powering to the native scalar backend. It previously
ignored negative exponents and left negative residues at signed discriminants.
The relation checker retains native PariError messages, including failed inverses
and the full/large-prime relation bug diagnostics. Its class-group early returns
are preserved. The separate post-Gauss diagnostic is covered by the later restoration below;
this section covers mpqs_factorback and mpqs_check_rel.

Rationale: preserve native symbols, modular products and errors instead of local
approximate/unsigned-only substitutes. Trade-off: consumers relying on incorrect
values or Error rather than PariError see the corrected behavior. Behavioral
impact: 4,802 native comparisons match, repairing 1,427 baseline discrepancies,
including 48 isolated nonterminating inputs. Class-group entry points remain
unported; the shared helpers are directly comparable to the bundled source.


## MPQS word inverse and source controls — 2026-09-12

The MPQS Fl_inv facade now uses the audited word-inverse dispatch through Fp_pow
with exponent -1. Native gcdll.c explicitly requires 0 <= a < p; p is positive
and both Number arguments must represent exact integers. The former local
Euclidean helper had correct values on the reduced controls but threw a generic
Error with different text for nonunits. Its 200 failing-input records now retain
the native PariError diagnostics. Inputs violating the native reduction
precondition were excluded from behavioral claims.

Rationale: preserve delegated native word inversion, including error behavior.
Trade-off: callers catching the old generic error must handle PariError.
Behavioral impact: 664 reduced inverse comparisons match, including the valid
modulus-one case. Another 429 initialization/factor-base/candidate controls
match before and after; they expand coverage without production behavior changes.
The class-group entry points themselves remain unported.

The candidate-scanner oracle supplies 2,016 long slots, as the TypeScript port
already does, so finishing an eight-byte block near the nominal 2,000-candidate
limit cannot overwrite adjacent memory. The original scanner code is unchanged;
its returned count, candidate order and terminator are compared. Initialization
comparisons never exceed 252 candidates; the dedicated capacity controls reach
2,002. This tests the scanner with adequate storage, not native out-of-bounds
allocator behavior.


## MPQS relation hash representation

MPQS relations use the native GEN shape [Y, Vecsmall(relp)]. The internal
_mpqs_hash.ts backend represents that shape as {Y: bigint, relp: number[]} and
reconstructs native type, length, sign and 64-bit magnitude words when hashing.
Word entries must be exactly represented Number integers. RelationTable follows
hash_init, hash_search2, hash_insert2 and hash_keys_GEN, including the native
prime-size sequence, 65% resize threshold, head insertion, collision-chain
relinking and duplicate suppression. Values must remain unchanged while stored,
as with native hash keys. This is the subset required by MPQS, not a general
implementation of all GEN types or every PARI hash-table operation.

Rationale: relation order determines the matrix columns and therefore the ordered
factors and Lanczos random-state consumption. Trade-off: consumers of the internal
MPQS API see the corrected factor order. Behavioral impact: 94 seeded MPQS driver
comparisons now match, repairing 45 records, including three final-state changes.
Another 136 native hash/table controls match, covering duplicates, deliberate
collisions, growth thresholds, signed/multiword integers and size-overflow errors.
The MPQS pair-array representation is retained: an early native integer factor
is adapted to that factor plus its residual cofactor, and native ifac class flags
are omitted, as in the existing TypeScript interface.


## Restored MPQS nonfatal post-Gauss diagnostic — 2026-09-12

With debugging enabled, native MPQS warns when a kernel relation fails
X^2 = Y^2 modulo N, then continues computing gcds. The TypeScript port previously
threw, suppressing the native result. It now emits the same core warning text
through console.warn and continues. The earlier per-relation mpqs_check_rel
failure remains a native PariError; these are distinct source paths.

Rationale: preserve the original nonfatal diagnostic and returned factors/null.
Trade-off: debug consumers must observe the warning instead of catching the old
exception. Behavioral impact: 216 of 240 diagnostic records are corrected.
The native harness captures the warning call while forwarding it to pari_warn;
comparisons verify the core message and count. Native terminal prefixes/styling
are represented by the host console's warning output. This is a logging adapter,
not a change to the mathematical diagnostic. Default debugging remains off.

Another 232 candidate-evaluation controls match, including 7,642 class-group
relations, and eleven factor-base early-return controls match. The class-group
entry points remain unported; the shared candidate evaluator itself now has
native coverage through the bad-prime-at-two branches with debug checks enabled.

## Restored Pollard lambda signed hashes and seeded steps — 2026-09-12

Sage `groups/generic.py:1214` draws lambda step sizes from `sage.misc.prandom`,
which uses the current state's CPython generator. The port used the GMP stream.
It now uses `python_random().randint(1, N-1)`, equivalent to native `randrange(1,N)`.
Python's nonnegative hash remainder is also restored: JavaScript's signed `%`
previously missed the step-table entry and silently substituted a unit step.
Invalid-bound errors now preserve the source's exact spacing.

Rationale: reproduce native walks, error text and the chosen log in intervals
containing multiple solutions. Trade-off: seeded callers may receive a different
valid logarithm or observe a different subsequent random state than with the old
port. Behavioral impact: 1,038 records differed in the 1,340 added controls;
results, errors and complete walk traces now match for the explicit shared hash.
The default string hash and the width-zero workaround remain documented adapters;
this repair does not claim native default hashing for arbitrary group objects.


## Restored Pollard rho draws, shortcuts and validation — 2026-09-12

Bundled Sage `groups/generic.py:774–797` samples all 20 m coefficients, then all
20 n coefficients and the initial exponent through IntegerModRing.random_element.
The port interleaved m/n and drew from GMP. It now delegates to the existing ring
method, which uses Sage's CPython stream. Identity and target-membership shortcuts
that skipped native work are removed; the small-order BSGS interval includes the
supplied order. Nonprime orders and exhausted retries retain native error text.

Rationale: preserve the selected result/error and random-state side effects.
Trade-off: identity targets and invalid subgroup requests now perform native
walk/retry work instead of returning early. Behavioral impact: 887 of 1,304 new
records differed before repair. All new comparisons now match result, error and
the next 64 CPython random bits. Tests execute the bundled original Python module
under the installed Sage runtime; this avoids a 10.3-only BSGS identity error.

Hashing remains the documented adapter: the oracle explicitly receives the same
DJB2 string hash the port uses. These tests do not claim native default object
hashing, arbitrary nonunique string-key equality, or the unported optional-order
and custom-hash arguments. A stale earlier statement that Sage accepted composite
rho orders is corrected: the source explicitly rejects them too.

## Restored generic group parsing and operation schedules — 2026-09-12

Bundled group parsing rejects incomplete custom identity/inverse/operation sets
before validating bounds or declared orders. The port now shares that validation
and preserves the exact native message, including for arbitrary nonstandard
operation strings accepted at runtime. The public OperationType union is unchanged.

Custom `multiple` now follows `generic.py:342–379`: detect idempotence, handle
powers two through four, then use the native binary schedule. Standard-operation
element fast paths remain. `order_from_bounds` divides bounds exactly with BigInt
ceil/floor semantics, including negative bounds, and passes `check: false` to
`order_from_multiple` after BSGS, matching the original dependency call.

Rationale: preserve validation precedence, exact bounds and native callback work.
Trade-offs: callers observing callback count/order or catching the old custom
error text see the native behavior instead. Behavioral impact: 336 parsing
records, 154 bound/trace records (27 result/error differences), and 563 standalone
multiple traces differed before repair. All 1,620 new comparisons now match the
bundled original. A min operation with exponent 2^1000 takes one callback instead
of 1,002; a lower bound of -2^54 divided by 2^54 stays -1 rather than becoming 0.
This does not assert identical callback schedules for standard element fast paths.

## Restored order-from-multiple factor-list semantics — 2026-09-12

Sage treats empty factorization and prime lists as absent. The port treated empty
arrays as present and recursed into an empty factor tree until stack overflow.
Both optional lists now use their lengths, retaining explicit-factor precedence.
Prime-list valuations delegate to the existing Integer valuation/GMP removal path;
the former division loop hung at base 1 and rejected other invalid bases incorrectly.

At a zero multiple, the source valuation is +Infinity. The single-base helper is
preserved, including terminating prime-power cases. Multi-factor infinite costs
retain the native phase/sign/power errors. SignError is represented by an
ArithmeticError with its native name. This implements these specific source paths,
not a general symbolic infinity ring. Prime-list entries in the oracle are Sage
Integers, matching the IntegerLike adapter; Python built-in integer entries can
have a different infinity-power error text. Source-nonterminating cases remain
nonterminating and are excluded from permanent execution cases.

Rationale: preserve empty-list defaults, backend delegation, results/errors and
callback traces. Trade-offs: consumers catching stack overflow or division errors
now receive a result or the native validation error. Behavioral impact: 96 of 384
list controls and 77 of 96 valuation controls differed, including 63 former port
hangs captured in isolated processes. All 480 permanent comparisons match.

## Generic group iterator and parent adapters

### Restored behavior (2026-09-12)

Sage's multiples class validates and copies the step/initial value at construction.
Each next call increments the index and computes the following value before
returning; an operation error leaves the incremented index and previous value
available for the next call. The former TypeScript generator deferred validation
and copying, performed the operation after yielding, and closed on any exception.
These source behaviors are now restored. Standard identities come from the parent
when available. Parent lookup accepts the port's stored-parent and parent-method
representations, fixing parseGroupOps on number-field elements.

### Host adapters and limits

Python copy.copy maps to a callable __copy__ hook, a callable copy method, or a
shallow copy of own property descriptors with the prototype retained. Objects
with JavaScript private/internal state need a copy hook. This is not a port of
Python's complete pickle/reduction protocol. An actual generator object retains
the public TypeScript Generator interface and inherited iterator helpers, while
its next method runs the native state machine. JavaScript return/throw and for-of
closure behavior are retained; Sage's iterator does not expose these methods.

Rationale: reproduce construction, callback timing, copies and recovery while
preserving the host interface. Trade-offs: constructor errors occur earlier and
callbacks run before their values are returned; custom objects must follow the
copy adapter. Behavioral impact: 576 of 608 new state/parent records differed
before repair. All now match, including initial-value copy failures and 32 direct
quadratic number-field parent comparisons. Arbitrary host-object internals and
Python reduction protocols remain outside the adapter.


### Number-field reduction and native array adapters

**Restored 2026-09-12:** Number-field construction reduces all input coefficients
modulo the defining polynomial through FLINT's rational polynomial remainder.
Degree-one generators represent the root, and defining-polynomial getters and
field display preserve the original polynomial's scaling. Integer powers use
Sage's generic_power schedule and return the original element for exponent one.
Other number-field backend architecture gaps remain open.

The new FLINT pseudo-division and rational-remainder kernels use immutable input
arrays and newly allocated outputs instead of native mutable buffers. The basecase
retains fixed quotient/remainder slots; the public divide-and-conquer adapter
normalizes trailing zeros. Comparative tests trim only unused output slots when
comparing these arrays with native public polynomial objects. Invalid zero-divisor
or nonpositive-denominator buffers raise RangeError instead of entering a native
abort or violating its preconditions. Rationale: represent safe native polynomial
storage in TypeScript without pointer aliasing. Trade-off: native output-buffer
aliasing and invalid-memory behavior are unavailable. Behavioral impact: valid
canonical polynomial results and pseudo-division exponents are unchanged. The
short divide-and-conquer branch expresses native buffer moves as the same low-part
polynomial subtraction using the existing FLINT multiplication port.

The internal arith/power helper implements the Sage-element protocol supplied by
number-field callers (`mul`, `inv`, `parent().one()`); Python's arbitrary-object
constructor and operator fallback is outside this helper's type domain. It is not
exported from the public arithmetic namespace. This protocol mapping avoids
pretending that JavaScript objects support Python operators; supported element
powers retain the native value, operation order, inverse errors and aliasing.


Scaled quadratic-unit comparisons use native PARI quadunit expressed at the
larger real root, matching the port's documented generator normalization;
Sage bnfinit can choose another fundamental generator. Unit log/exp cases compare
exponent round trips, roots of unity compare the complete sorted coefficient set,
and regulator controls compare integers rounded after multiplication by 10^8
(the existing approximate-regulator deviation remains). The torsion_order port
alias maps to Sage's zeta_order. These adapters do not assert bitwise equality of
regulators or a shared bnfinit basis.

The scaled equation-order different controls in this batch have a monogenic
maximal order. The existing derivative-based different calculation is not valid
for arbitrary maximal orders; that older gap requires a separate backend repair.


### Number-field ideal backend adapters

**Sage/PARI vs port:** Sage sends ideal inversion, products and the field different
through PARI's number-field object. The port now delegates these operations to
`parigp-ts/base1.ts` and `base4.ts`, using the existing integral basis and multiplication
table. Ideal comparisons normalize rational power-basis lattices, since the chosen
integral basis and displayed generators can differ. Product generators are bounded
by the field degree rather than retaining every pairwise product.

For degree at least six, native `mat_ideal_two_elt` can factor the part of the ideal's
intersection with Z supported on primes at most 47 and use `idealapprfact`. The port
currently uses PARI's certified `get_random_a` route for this branch too. **Rationale:**
the smooth-part ideal-approximation dependency is still unported at the PARI layer.
**Trade-offs:** that specialization's performance and random-state consumption are
not preserved; exact branch/RNG fidelity remains open. **Behavioral impact:** the
returned ideal HNF is certified and matches the comparative cases through degree 16;
the selected two-element generator and subsequent random outputs can differ.

The internal integer-matrix inverse uses an exact Hadamard iteration bound in place
of native low-precision real arithmetic. Half-GCD expands the native basecase Lehmer
batches into exact Euclidean steps, retaining the 66-word recursive splits and
fixups. These adapters avoid dependence on native real/GMP buffers; they can change
iteration counts and constant factors. Direct bundled-native comparisons include
matrices with 4,096-bit entries and half-GCD inputs through 16,384 bits. No result
difference was observed. Native `ZM_inv` can return a zero matrix with denominator
zero for singular matrices of dimension at least three; this low-level result is
preserved, while the trace-pairing caller rejects it. Array shape/domain guards are
port preconditions, not claims about native C behavior on invalid buffers.

`AbsoluteOrder.different()` and `codifferent()` are port convenience methods absent
from Sage's Order API; they refer to the ambient field's maximal order. They do not
implement a separate different for a nonmaximal supplied order basis. This preserves
the existing convenience API while fixing its derivative shortcut; support for
nonmaximal order objects remains incomplete. Ideal `pow(IntegerLike)` follows the
Sage integer-power schedule and object reuse. Bigint versus Integer wrappers also
preserve the distinct Python-int and Sage-Integer negative-zero-ideal error paths.


### Ideal basis and generator representations

**Sage vs port:** `integral_basis()` now returns the HNF integral-basis elements,
including fractional and nonprincipal ideals. Both sides have rank zero for the
zero ideal. The port's older `zk_basis()` name is a convenience alias for the same
basis construction; Sage instead exposes `basis()` and `integral_basis()`. The
chosen integral basis may differ, so comparative fixtures compare the exact
rational power-basis lattices and ranks, not arbitrary displayed generators.

`free_module()` retains the existing `{ basis, rank }` record, whose basis entries
are field elements, rather than Sage's full cached FreeModule object with rational
coordinate vectors. The record is now cached too. **Rationale:** preserve the
port's existing consumer interface while fixing the represented module.
**Trade-offs:** module operations and immutable basis-container behavior are not
provided by this record. **Behavioral impact:** lattice, rank and repeated-result
identity match the comparisons; object methods and container mutation differ.

`gens_two()` now returns two field elements and caches the pair, matching Sage's
parents, fractional QQ intersection generator and zero second generator for
rational ideals. It uses PARI's two-element reduction, with the already documented
higher-degree smooth-part specialization gap. The cached TypeScript tuple remains
an ordinary array under the project's tuple convention, unlike Python's immutable
tuple; callers must not mutate it. The old bigint first-entry type could not
represent fractional intersections, so its correction requires version 22.0.0.


### Fractional ideal arithmetic and remaining intersection routing

**Sage vs port:** numerator and denominator now mean coprime integral ideals, as in
Sage, and are cached. The older integer denominator was an incorrect API mapping;
version 23.0.0 changes its return type. Fractional coprimality follows Sage's four
numerator/denominator ideal-sum tests. The zero ideal has a different class in Sage,
so numerator, denominator, divides and is_coprime are absent there. The port's shared
class keeps the method names but raises the corresponding AttributeError for zero.
**Rationale:** retain the current class structure while matching method-call behavior.
**Trade-off/impact:** introspection can still see methods that Sage's zero ideal lacks;
these call results and error text match the comparative fixtures.

Sage's single-generator division uses field-element division; other fractional
ideal division delegates to PARI. Zero left operands use the monoid inverse/product
path. The port now follows that dispatch and PARI's inversion-before-zero-dividend
order. Native empty-matrix-to-ideal conversion can raise TypeError for zero products
with multiple generators; that behavior is retained. Empty/all-zero field-factory
generators normalize to a fresh single-generator zero ideal.

**Intersection routing closed in 23.0.10:** intersection now delegates to the PARI
`idealintersect` port and its LLL kernel of concatenated rational ideal HNFs. See
`Number-field ideal intersection and construction adapters` for the typed boundary
and remaining coercion/display limitations.

## PARI fast LLL floating-point profile

**Sage/PARI:** `lll.c` uses C-double Gram-Schmidt in its initial fast stage. Compiler
contraction and ABI behavior can change its partial basis, transformation and failure
status. The public LLL wrapper subsequently uses heuristic and certified stages.

**Port:** The new internal `fplll_fast` helper reproduces the audited AArch64 Apple
Clang 21 `-O3` profile using explicit binary64 FMA and separate-product boundaries.
It copies nonempty rectangular column matrices and returns partial B/U with status;
U starts at the identity or is absent. Invalid empty/ragged adapter inputs raise
RangeError instead of invoking a native routine whose shape preconditions are unmet.
The helper is not connected to ideal intersection until the native wrapper is ported.

**Rationale:** Separate JavaScript multiplication and addition disagree with the
normal native build in 120 completed probes. Exact fused rounding and the native
loop boundaries reproduce those results without weakening the oracle or disabling
its normal optimization. Explicit AArch64 signed casts/shifts also preserve partial
failure states where the upstream C expression has undefined behavior.

**Trade-offs:** This fixes one concrete native floating profile; it does not promise
identical intermediate transforms for other C compilers, architectures or FP modes.
The live stage oracle requires AArch64 Clang. Software FMA costs BigInt work. The
fast stage alone is uncertified and may not terminate on dependent keep-first bases.

**Behavioral impact:** 1,970 completed stage comparisons cover optional transforms,
keep-first, rank deficiency, zero columns, four reduction parameter pairs, dimensions
through 17 and coefficients through 4,096 bits. Another 3,375 native binary64 cases
check FMA, scaling and decomposition, including special values and exact bit patterns.
A separate inventory retains 94 native probes exceeding a two-second exploratory
bound; 86 of those were also checked to exceed the same TypeScript bound. Timeouts
are not counted as successful comparisons. Cross-profile equivalence remains open. The adaptive wrapper and
ideal-intersection routing now use the ported stages.

## PARI DPE LLL stage and resource boundary

**Sage/PARI:** `lll.c`'s DPE stage uses exact integer Gram updates and normalized
binary64 significands with signed-word exponents. Its caller may request incremental
Gram construction, a supplied Gram, no basis, optional U and output squared norms.
Dependent keep-first inputs can fail, loop, produce a dbltor overflow while converting
norms, or request a temporary with an exponent near LONG_MAX and exhaust the PARI stack.

**Port:** The internal `fplll_dpe` stage follows these updates and returns copied
`[status, G, B, U, norms]`, retaining partial failure state and native norm records.
It rejects empty/incompatible adapter matrices with RangeError. A requested shift
above the exact JavaScript integer range raises RangeError with
`fplll_dpe: Gram-Schmidt coefficient requires an unrepresentable shift`, before
allocation. It does not emulate a configurable PARI stack or its allocation diagnostic.

**Rationale:** Near-LONG_MAX shifts produced by nonfinite Gram-Schmidt coefficients
are impossible to materialize. Native temporary reservation fails even when the
rounded multiplier is zero; skipping that reservation had incorrectly returned a
normal stage-failure status for six probes. An explicit guard preserves a failure
boundary without attempting an enormous allocation or inventing a native stack size.

**Trade-offs:** The error class/message at this resource boundary intentionally
differ. Nine reached native stack failures have a separate comparative dispatcher
which matches the explicit port guard; those records are not exact error comparisons.
Normal native errors use cypari2's punctuation convention, verified directly for
`dbltor(NaN)`. The adaptive wrapper now composes all four stages and FLATTER.

**Behavioral impact:** 3,588 exact comparisons check Gram/basis/transform/norm state
or ordinary errors, and nine additional records check the documented resource boundary.
They cover all three Gram/B modes, optional U/norms, four reduction parameter pairs,
dimensions through 17 and coefficients through 4,096 bits. A separate inventory keeps
75 native probes exceeding the two-second exploratory bound; these are not passing
comparisons. The adaptive wrapper is also used by the ideal-intersection LLL-kernel route.

## PARI real LLL stage adapters

**Sage/PARI:** Private `fplll_heuristic` and `fplll` stages update pointer arguments
and use PARI real arithmetic with precision escalation performed by their caller.
The native adaptive caller uses positive word-aligned precision, beginning at 64 bits
and increasing by 64. The heuristic stage can allocate its Gram approximation at a
different precision. The proved stage can accept an exact Gram with no basis and
return squared norms. Intermediate failure is distinct from complete LLL reduction.

**Port:** Pure column adapters copy inputs, start optional U at the identity, and
return native status/partial state. Precision is in bits and must be a positive
multiple of 64; empty or incompatible matrix shapes raise RangeError. Native callers'
word-aligned precision sequences are supported; direct calls with non-word-aligned
precision are outside this adapter contract. Zero-pivot division errors preserve
the native Gen-boundary error name and payload, with cypari2's punctuation convention.

**Rationale:** Matching the native arithmetic kernels, destination precision and
integer-rounding branches is needed to reproduce reduction representatives and
failure states. Keeping a checked word-aligned contract matches the actual native
wrapper and the bit/word assumptions of its mantissa-truncation step.

**Trade-offs:** These are individual stages with identity-started transformations,
not the complete adaptive/FLATTER wrapper. Native garbage collection, configurable
stack limits and direct non-word-aligned precision calls are not emulated. The
existing real-kernel padding/precision boundaries remain documented separately.

**Behavioral impact:** 3,024 heuristic stage cases, 4,908 proved stage cases and
5,406 direct scalar cases compare exact outputs or errors against bundled PARI.
They cover bit precisions through 512 in the stages, scalar precision through 1,024,
coefficients through 4,096 bits, optional B/G/U/norms, distinct Gram/GSO precision,
half-integer ties, word-rounding wrap and trailing mantissa words. Twelve native and
TypeScript probes exceed the same two-second exploratory bound; their inventory is
separate and they are not passing comparisons. The full wrapper and ideal-intersection
routing remain open.


## PARI Householder QR adapters

**Sage/PARI:** Internal Householder routines consume t_INT/t_REAL matrices and output
GEN pointers. QR_init accepts square or tall inputs used by lattice reduction;
gaussred_from_QR requires a square full-rank input. Precision failures return zero
without assigning output pointers. The native Q vector declares n-1 public entries;
for tall inputs its last scratch reflector overlaps B's allocation header.

**Port:** Pure adapters use zero-indexed columns with bigint/MpReal cells and return
copied numerical state. QR_init returns `[success, B, Q, L]`, with null outputs on
precision failure and n-1 public reflectors on success. R_from_QR and
gaussred_from_QR return matrices or null. Nonempty compatible shapes and positive
64-bit-aligned working precision are checked with RangeError.

**Rationale:** These are dependencies for native FLATTER, where integer zeros, real
accuracy, rounding and representative selection matter. Separate managed arrays
preserve all numerical entries without emulating native allocation metadata.

**Trade-offs:** Generic complex/rational coefficients, invalid native pointer/shape
contracts, garbage collection and PARI stack exhaustion are not emulated. This is
not the generic public matqr API. Native scratch-header aliasing has no TypeScript
analogue; the native oracle reads B by its known size to avoid interpreting that
header as a vector length.

**Behavioral impact:** 4,752 exact native comparisons cover QR state, upper R and
Gaussian reduction, all-integer/all-real/mixed cells, square/tall/degenerate matrices,
working versus input precision through 8,192 bits, large coefficients and shifts,
and success/precision-failure branches. Full LLL/FLATTER routing remains open.


## PARI real Cholesky and triangular inverse adapters

**Sage/PARI:** Generic qfgaussred_positive, RgM_Cholesky and RgM_inv_upper support
several numeric coefficient types. FLATTER converts its matrices to real arithmetic
before calling these dependencies. Gaussian reduction returns null for a nonpositive
pivot and ignores the lower triangle. Real square roots keep operand precision.
The generic numeric boundary distinguishes exact zero from finite-accuracy real zero.

**Port:** These adapters support MpReal and exact 0n cells in zero-indexed columns,
preserve inputs and execute the native update order. qfgaussred_positive and
RgM_Cholesky return null for nonpositive pivots; RgM_inv_upper uses native back
substitution. Empty matrices work. Non-square Gaussian inputs retain the native
PariError; triangular shape and nonpositive/non-word-aligned Cholesky precision are
checked with RangeError. `_real_matrix.ts` supplies only the MpReal/0n numeric branches.

**Rationale:** This is the coefficient domain needed by the real FLATTER route.
Retaining exact zero and real allocation/accuracy is necessary to match rounding,
precision failure, transformed representatives and native error dispatch.

**Trade-offs:** These internal adapters are not generic rational/complex matrix
arithmetic and do not yet implement those input domains. Native garbage collection
and configurable stack exhaustion are not modeled. Large allocated-zero reciprocals
retain the deterministic rejection documented under Native PARI allocated zero
records; native Newton inversion reads unspecified zero payload words.

**Behavioral impact:** 10,445 exact native comparisons include matrix state and direct
numeric branches, bit precision through 8,192, square dimensions through nine,
nonpositive pivots, exact/inexact zeros, disparate input/working precision and error
names/payloads. Nineteen large allocated-zero matrix inverses are retained separately:
four native probes returned zero records and fifteen raised errors. They are not
passing comparisons. Port-only regressions verify the deterministic guard. Generic
coefficient domains, FLATTER/full LLL and ideal-intersection routing remain open.


## PARI integer and real matrix products and rescaling

**Sage/PARI:** Generic RgV dot/Gram/product routines support many coefficient types
and delegate integer products to ZM_mul. RgM_rescale_to_int supports integer, rational
and real coefficients. In mixed inexact input it chooses a binary scaling exponent
then rounds all entries; integer entries can lose exact ratios. All-zero real input
with negative accuracy overflows when shifted by the untouched HIGHEXPOBIT sentinel.

**Port:** Internal RgV.ts adapters support bigint/MpReal coefficients in zero-indexed
columns, preserve native arithmetic order and delegate integer products to the existing
ZM_mul port. RgM_rescale_to_int supports the same integer/real domain, retains native
rounding and the all-zero expo-overflow error, and preserves all-integer inputs.
Malformed matrix/vector dimensions raise RangeError at the typed boundary.

**Rationale:** These coefficient domains and native rounding choices are required
by FLATTER's QR, Gram, size-reduction and real-to-integer steps. Keeping the integer
backend avoids replacing its multiplication algorithms with generic scalar loops.

**Trade-offs:** Generic rational, modular, polynomial and complex coefficient types
are not supported by these internal adapters. Native scratch allocation, garbage
collection and configurable stack limits are not modeled. The existing ZM_mul
representation and arithmetic-kernel boundaries continue to apply.

**Behavioral impact:** 5,903 exact native comparisons cover dot/square/Gram/products
and rescaling, empty and rectangular matrices, mixed exact/inexact zeros, precision
through 8,192, coefficients through 4,102 bits, negative half-integer ties, all-zero
accuracy overflow and integer products through 72 by 72 crossing backend cutoffs.
The typed FLATTER/adaptive LLL wrapper is now implemented; ideal-intersection routing remains open.


## PARI adaptive Gram-Schmidt dependencies

**Sage/PARI:** FLATTER uses private exponent statistics and adaptive QR/Cholesky
routines. A full-column-rank integer basis or positive-definite integer Gram is a
caller precondition. Precision begins at n+31 bits, doubles on failure and otherwise
grows by the maximum of 4/3 and the condition/spread requirement. Allocations round
to whole 64-bit words. Upper-triangular bases have a direct conversion shortcut.

**Port:** lll.ts re-exports the corresponding pure dependency adapters from
_lll_gso.ts. Integer/real exponent statistics return bigint; integer zeros retain the
native HIGHEXPOBIT sentinel. Matrices contain zero-indexed columns. Empty/malformed
shapes and mismatched diagonal types in drop raise RangeError. A computed allocation
precision outside the exact JavaScript integer range raises RangeError.

**Rationale:** Precision selection changes the actual real records and later FLATTER
representatives. Preserving the native retry schedule and exponent sentinels is needed
in addition to matching the QR/Cholesky arithmetic itself.

**Trade-offs:** The full-rank/positive-definite preconditions remain the caller's
responsibility; invalid inputs can retry until an allocation error, as upstream.
The port does not reproduce native stack limits, unchecked signed-long overflow or
undefined pointer/shape contracts. Its unrepresentable-allocation guard is explicit,
not an exact emulation of PARI's configurable resource error.

**Behavioral impact:** 1,230 exact native comparisons cover returned QR/Cholesky
records and all seven statistics, integer/real coefficients, square/tall bases,
upper/lower/reversed and signed diagonals, dimensions through 17, basis coefficients
through 8,198 bits, real result precision through 16,512 bits and statistic exponent
shifts through +/-10^9. Unit tests separately cover the typed and resource guards.
Rank and the typed FLATTER/adaptive LLL wrapper are now implemented; ideal-intersection routing remains open.


## PARI word matrix pivot and solve adapters

**Sage/PARI:** Flm_pivots selects native pivot rows and nullity; Flm_gauss solves over
a word-prime field and returns null on rank failure. The library uses small Gaussian
elimination or its cutoff-eight CUP decomposition and recursive triangular solves.
Overdetermined solves select independent rows and do not verify remaining RHS rows.

**Port:** Public dependency adapters retain those algorithms and values with unused
slot zero on both matrix axes and on nonnull pivot vectors. Inputs are preserved.
The shared row kernels return pivot values in one-based indexing without a dummy
vector entry. Word prime range and malformed adapter dimensions raise RangeError;
primality remains a caller precondition.

**Rationale:** The representation matches existing word-matrix dependencies while
sharing the native CUP and solve kernels needed by integer rank certification.

**Trade-offs:** Native destructive scratch updates, pointer layouts, lazy unsigned
intermediates and garbage collection are not exposed. Arithmetic intermediates are
reduced with BigInt; the prime-field results and source algorithm branches are kept.
Unsupported composite moduli and invalid pointer contracts are not emulated.

**Behavioral impact:** 2,802 exact native comparisons cover pivots, rank failures and
solve representatives, word primes from 2 through 2^64-59, empty/square/tall/wide
matrices, sparse/dependent rows and columns, dimensions through 65 and both sides of
the cutoff-eight boundary. Integer rank and solve callers use these same kernels.

## PARI modular integer rank and solve adapters

**Sage/PARI:** ZM_pivots and ZM_rank use word-prime rank profiles and exact kernel
certification. ZM_gauss selects independent rows, solves modulo big word primes,
combines growing CRT batches, reconstructs rational coefficients and verifies the
candidate. Overdetermined systems use selected independent rows. In a native edge
case, primes dividing all entries can leave a zero-rank candidate and cause an
`inconsistent dimensions in gauss` error during certification.

**Port:** These algorithms share row kernels with the word adapters. Public matrices
and pivot-vector storage use unused slot zero. ZM_gauss returns `[N, D]`, where N is
an integer column matrix and D is the positive common denominator, or null on rank
failure. The observed native rank-certificate error is preserved. Matrix inputs are
copied; malformed dimensions and a solve with fewer rows than columns raise RangeError.

**Rationale:** Exact certification is required after a low rank modulo trial primes;
returning that provisional rank can be wrong. A common denominator maps native
rational matrices into the dependency's existing integer-matrix representation.

**Trade-offs:** C worker execution is sequential in TypeScript, with a balanced CRT
tree and native batch-growth schedule retained. Native GC, configurable stack limits,
unchecked prime-list exhaustion and pointer mutation are not emulated. Solves accept
matrix RHS values only; implicit-identity/column GEN overloads are outside this adapter.
The native full-column-rank/selected-row behavior is preserved, including errors.

**Behavioral impact:** 2,616 exact native comparisons cover pivot profiles, rank,
solutions and errors, dependent/sparse/rectangular matrices, misleading small primes,
failed certificates, later better ranks, bad big CRT primes, denominator reuse and
nonuniform coefficient sizes. No trial-only rank or inverse-based solve substitutes
for the native algorithm. The typed FLATTER/full LLL wrapper and ideal-intersection routing are implemented.

## PARI adaptive LLL and FLATTER adapters

**Sage/PARI:** lll.c selects quadratic-form reduction for the two-column image
shortcut, computes modular rank, chooses FLATTER using its fixed tuning tables,
then runs double, heuristic, DPE and arbitrary-precision stages. FLATTER recursively
reduces overlapping blocks and retains the original stopping statistics. Native
functions mutate basis/transformation state and accept generic GEN matrices.

**Port:** Internal `src/lll.js` exports `ZM_lll`, `lllfp`, `flat`, `ZM_flatter`,
`ZM_flatter_rank`, `ZM_flattergram`, the native flag constants and `LllResult`.
Matrices are zero-indexed columns, and inputs are preserved. Stage transformations
are composed exactly between calls. ZM_lll accepts bigint coefficients; lllfp also
accepts MpReal coefficients. Inexact Gram input currently requires MpReal or exact
integer zero cells; nonzero exact integer cells mixed with reals raise RangeError.
Full column rank is a caller precondition for direct flat/ZM_flatter, positive
definiteness for direct Gram compression, and the direct routines require at least
two columns. Rational/complex/generic GEN input is not exposed. The norms-returning wrapper
is available through the typed adapter described below. Native debug/memory-accounting machinery is not reproduced.

**Rationale:** These adapters let the port use the original lattice algorithms and
backend dependencies while retaining explicit TypeScript numeric and ownership
conventions. The Gram coefficient restriction is the domain of the audited real
Cholesky dependency. Native ZM_lll_norms crashes in the two-column image shortcut
when norms are requested (confirmed by the isolated bundled-C probe); ZM_lll itself
passes no norms pointer and avoids that source bug.

**Resource limit:** The original rank-deficient Gram FLATTER augmentation loop can
run indefinitely: its inner routine returns a nonnull identity transformation,
while the outer routine requires null to stop. Three eight-column probes timed out
at five seconds. The port follows the same arithmetic and preserves exceptions,
then throws `RangeError('rank-deficient Gram FLATTER exceeded 64 augmentation attempts')`
after 64 attempts. These three controls are recorded separately and are not exact
comparative successes. The return after a null transformation remains in source
and coverage despite being unreachable for the nonsingular transformations in this
domain. The native rank-deficient INPLACE path's use of augmented RHS columns is
preserved, including ZM_mul's truncation to the left operand's column count.

**Trade-offs / behavioral impact:** Supported deterministic results, flags and
arithmetic errors are compared exactly. Mutability, restricted coefficient domains
and resource failure presentation differ as described above. The native two-column
large-input sign behavior is preserved through `qfb.redimagsl2`, whose transformation
is row-oriented like the rest of qfb.ts and transposed at the LLL boundary. This
new compatibility entry retains the original negative-b Schoenhage result, including
determinant -1. The existing public qfbredsl2 correction remains unchanged; the old
claim that its >9,000-bit path is unreachable does not apply to the new LLL wrapper.

The four existing LLL stages now multiply by their word/truncated mantissa before
shifting, as native mului/submulshift do. Expanding the full coefficient first had
preserved values while worsening multiplication cost; permanent large-exponent
comparisons and a generously bounded end-to-end regression cover that correction.

## PARI adaptive LLL norm output

**Sage/PARI:** ZM_lll_norms optionally writes squared Gram–Schmidt norms through
pN while reducing a lattice or Gram matrix. Trivial and successful NOCERTIFY
shortcuts leave that pointer unset. The bundled two-column positive-definite image
shortcut passes a t_QFB to a matrix routine and raises a caught segmentation-fault
PariError when pN is requested. Ill-conditioned KEEP_FIRST inputs can yield
zero-precision norm records, exhaust the PARI stack, or fail to terminate.

**Port:** ZM_lll_norms returns [result, norms], using null for an unset pN and
MpReal records for supplied norms. Inputs are copied. The same adaptive stages
produce all native precision/exponent/mantissa records, including low-accuracy
zeros. The two-column source bug raises the same caught error explicitly, without
an unsafe memory access. Existing ZM_lll requests no norms and keeps its result shape.

**Rationale:** The tuple replaces native output pointers, and explicit errors allow
safe reproduction of observed native failures. The existing DPE enormous-shift
guard remains necessary; it does not emulate PARI's configurable memory stack.

**Trade-offs / behavioral impact:** Exact comparisons cover normal results and
ordinary errors, including the caught two-column error. Separate resource probes
compare native stack exhaustion to the documented RangeError guard; their error
messages intentionally differ. Two-second bounded probes retain native and port
nontermination observations separately from completed values. These are finite
observations, not proofs that every KEEP_FIRST input terminates or hangs. Some
source probes combine GRAM and INPLACE, flags the native manual declares incompatible;
matching those observations does not extend the advertised input contract.

## Number-field ideal intersection and construction adapters

**Sage/PARI:** NumberField.ideal delegates to fractional_ideal and maps ValueError
back to the zero ideal. Fractional construction accepts a generator list or an
existing fractional ideal, preserves same-field ideal identity and rejects zero.
Nonzero factory results are NumberFieldFractionalIdeal objects. Intersection coerces
its other operand to an ideal, then delegates to PARI's rational HNF/LLL-kernel
intersection. Rational constants from another abstract number field coerce exactly;
nonrational elements without compatible natural embeddings raise TypeError.

**Port:** The constructors and factories now preserve these supported behaviors,
including nested singleton generator lists at the native constructor boundary.
`NumberFieldIdealGenerator` is RationalLike | number | NumberFieldElement;
`NumberFieldIdealInput` additionally accepts a generator array or existing ideal.
Direct constructors take generator arrays; the optional Sage coerce keyword and
PARI GEN constructors are not exposed. Cross-field rational constants are accepted.
General compatible embeddings between distinct fields remain unimplemented.

Internal base4.idealintersect accepts the existing NfIdealData metadata and two
canonical rational ideal HNFs as zero-indexed columns with positive denominators.
It removes denominators, delegates to ZM_lll in kernel mode, multiplies the truncated
kernel by the first ideal, and delegates modular HNF normalization to hnf_snf.
Results are [integerColumns, positiveDenominator], with [[],1n] for zero. Inputs are
preserved; shape and denominator misuse raises RangeError. This replaces the old
I*J/(I+J) shortcut at the Sage layer.

**Rationale:** The constructor/class and delegation rules are observable behavior,
and the original intersection avoids the multiplication/inversion route and its
two-generator reduction work. The tested old and corrected Sage-layer calls both
preserve the observed post-call random state; no random-state defect is claimed. Explicit HNF and numeric adapters retain the port's
existing ownership and coercion conventions.

**Trade-offs / behavioral impact:** Supported ideal construction, coercion, identity,
zero validation and lattice results are compared with Sage. Generic GEN conversion,
relative-field ideal classes, the coerce keyword and general field embeddings
remain separate implementation gaps. Integral-basis conventions can differ between
the installed Sage oracle and the bundled source; ideal equality comparisons use
canonical power-basis lattices, while backend comparisons use identical HNF inputs.

Class-aware string prefixes now distinguish base ideals from fractional ideals.
Generator display still uses the supplied generators rather than Sage's full
`_gens_repr` reduction. Multi-generator/noncanonical pretty-printing remains an
observable gap; the string comparisons here cover zero and base principal ideals.


## Number-field ideal coercion and centered integral bases

**Sage/PARI:** Membership first coerces through the ambient field and catches only
TypeError. Coefficient vectors must have exactly the relative degree and coerce
through QQ. Addition and multiplication coerce non-ideal operands through the field's
ideal factory; addition then coerces the combined generators. Divisibility coerces
only non-ideal operands, while coprimality always coerces its argument to an ideal.
Division uses Sage's monoid coercion model, which rejects lists, foreign-field
parents and nonzero host floats. Numeric zero has a special canonical-coercion path;
Python-native scalar fallbacks can mask TypeError on a base zero ideal.

**Port:** These supported scalar, list, field-element and ideal routes now follow
Sage. `NumberFieldInput` is the exported scalar-or-coefficient-vector type. The
membership method also accepts an ideal object so its failed field conversion can
return false. Vector lengths and nonfinite coefficient errors remain ValueError.
Relative fields, general compatible embeddings, string/symbolic/PARI GEN conversion
and the full dynamic coercion model remain unimplemented. Directly created base-class ideals now preserve fractional-only invocation errors;
the remaining interface visibility adapter is documented under
Number-field ideal class method adapters.

The native final ZM_hnfcenter normalization is now applied to the maximal integral
basis, replacing the port's nonnegative residues. It uses diviiround with half ties
toward positive infinity. This closes the observed centered-basis discrepancy,
including a degree-eight case whose column reductions produce a negative half tie.
The native idealmul matrix route now checks dimensions, removes primitive content
before two-generator reduction and restores rational content afterwards. This also
preserves Sage/PARI's unusual interpretation of foreign ideal HNF coordinates when
its multi-generator multiplication dispatch skips field coercion.

**Rationale:** Coercion/error dispatch, exact basis representatives and primitive
content affect observable results or the native work schedule. Dependency routines
live in gen3.ts, hnf_snf.ts and base4.ts; the existing Sage-local maximal-order
algorithm delegates its final normalization to the PARI port.

**Trade-offs / behavioral impact:** The HNF adapters use zero-indexed columns and
copy inputs, unlike native routines that may reuse allocations. ZM_hnfcenter rejects
malformed shapes/nonpositive diagonals with RangeError; idealmul requires positive
denominators and canonical square HNFs (or [] for zero). Broader GEN domains are
excluded. Native direct comparisons check quotient ties, centered matrices, exact
integral bases and products including zero/dimension errors and random state.
The previously documented degree-at-least-six ideal-approximation specialization
and the remaining Sage-local maximal-order implementation remain open.


The existing partial class-number certificate searches bounded generator boxes.
To keep this certificate usable after changing the ambient integral basis, its
search now uses both nonnegative and centered HNFs in fixed power-basis coordinates.
Each box retains the existing candidate bound and exact norm/membership proof;
only certified generators can establish class number one. This adds a constant
factor to an already partial search, not a replacement for the missing native
bnfinit relation search. Forty-eight class-number/class-group comparisons cover
both signs of Dedekind's example and the existing certified fields under rational
rescaling of their equations.


## PARI permuted Hermite and knapsack adapters

**Sage/PARI:** ZM_hnfperm copies/reduces integer matrix columns while optionally
writing a full transformation and row permutation. ZM_hnf_knapsack tests that
HNF before restoring original row order. The bundled source permits ±1 entries
and empty rows, despite a narrower description in its manual; a zero-rank HNF
causes a caught segmentation-fault error when the source reads its first column.

**Port:** Matrices are zero-indexed columns; ZM_hnfperm returns [H,U,permutation]
with null for unrequested outputs. Permutation entries remain one-based. Inputs
are preserved, and knapsack column order, null rejection and caught source errors
are retained. Ragged input raises an adapter RangeError.

**Rationale:** Tuples replace output pointers and copied arrays replace native
ownership. Explicitly raising the caught native error avoids unsafe memory access.
**Trade-offs:** Native allocation reuse and invalid-buffer behavior are not exposed.
**Behavioral impact:** Covered matrix/transformation/permutation values and ordinary
errors match; the typed input contract requires rectangular integer columns.

## PARI factor recombination progress adapters

**Sage/PARI:** LLL_check_progress requests adaptive LLL norms, scans the native real
bound, truncates retained columns and returns null on its irreducibility decision.
It optionally accumulates debugging time. Missing norms or too few reduced columns
can make the bundled source dereference invalid memory; retaining zero columns
raises the native no-factor bug error.

**Port:** The same flags, precision selection and real comparisons determine the
returned columns. A copied matrix or null replaces GEN; the debug timing pointer
is omitted. The port safely raises the corresponding caught errors. n0 must be a
nonnegative integer no larger than the source row count.

**Rationale:** This preserves algorithm decisions without native pointers or memory
faults. **Trade-offs:** Memory/debug instrumentation and invalid vector-length writes
are outside the typed interface. **Behavioral impact:** All compared retained columns,
null decisions and caught errors match. The integer-polynomial factorization
backend uses this helper; constructor routing outside PARI's degree window uses
the NTL integer driver described below.

## PARI factorization scalar and bound adapters

**Sage/PARI:** QX_factor uses native real coefficient norms, Mignotte/Beauzamy bounds,
ceil_safe and a Fujiwara estimate followed by exact root-bound refinement. Mixed
integer/real comparison first rounds the integer to real precision. ceil_safe
accounts for uncertainty, including finite-accuracy zeros; it can exceed ordinary
mathematical ceiling. Odd unsigned half powers of negative reals return complex.

**Port:** The factor-bound helpers accept normalized positive-degree bigint
coefficient arrays. Scalar ceil_safe accepts bigint, rational pairs or MpReal;
rational pairs normalize through the native Qdivii-compatible construction adapter.
Half powers accept a real and an unsigned 64-bit exponent, returning MpReal or
MpComplex. Mixed comparisons return numeric -1/0/1 with native zero behavior.
vecbinomial uses native half-vector construction and reflection; sparse evaluation
uses integer powers across coefficient gaps. Inputs are preserved.

**Rationale:** These typed numeric records replace generic GEN while preserving
operation order, native rounding, pointer-free results and sparse arithmetic.
**Trade-offs:** Other GEN scalar/container types, native stack accounting and
malformed real records are not exposed by these signatures. JavaScript log2 supplies
the native binary64 Fujiwara estimate; general cross-platform libm bit equivalence
is not claimed. **Behavioral impact:** Native integer, real/complex record, bound,
comparison and error outputs match the registered cases. Exact refinement retains
the original conservative root bound. The NTL constructor ranges use the complete
integer driver described below.


## PARI bounded factor recombination adapters

**Sage/PARI:** QX_factor.c performs exact trial division, chooses lifting precision
with a 31-bit trace window and combines small sets of modular factors, returning
factor groups plus maximum-subset and completion output pointers. Larger unresolved
groups are left for LLL_cmbf. polarit2.c retains signed ties in centered remainders;
FpXV_prod uses balanced multiplication and returns integer 1 for an empty vector.

**Port:** Native private helpers are available as typed module functions. Ascending
bigint coefficient arrays replace integer GEN polynomials; divisors must be
normalized and nonconstant. ZX_divides_i's optional bound deliberately skips the
leading quotient coefficient, as in the source. cmbf_precs returns
[flag, a, b, q^a, q^b], with positive bounds and 2 <= q < 2^31. cmbf returns
[factors, modularGroups, maxK, done] and preserves caller input arrays. It requires
the original primitive/squarefree, modular-congruence and trace-bound contracts;
it does not certify arbitrary supplied factors. Empty and singleton modular products
preserve their native scalar/polynomial distinction. centermodii accepts a null
half-modulus for the original uncentered route.

**Rationale:** Tuples and owned arrays replace native output pointers and mutation
while retaining operation order, stopping decisions and values. **Trade-offs:**
Malformed GEN buffers, native stack ownership and debug timing are not reproduced.
The precision helper is restricted to the small-base range needed by recombination.
**Behavioral impact:** All registered exact quotients, remainders, products, bounds,
factor groups, completion flags and ordinary errors match bundled PARI. These
helpers leave native unresolved groups intact for the native trace-knapsack
phase. Sage's NTL constructor degree ranges use the complete integer driver below.


## PARI general Hermite form adapters

**Sage/PARI:** ZM_hnf dispatches to a private small-column routine for at most
seven columns and to ZM_hnfall for larger matrices. ZM_hnfall_i maintains native
pivot-height state, optionally updates a full unimodular transformation, and can
retain zero columns or remove them from one or both outputs. hnfall requests the
full transformation while removing zero columns only from the Hermite form.

**Port:** Rectangular zero-indexed bigint columns replace GEN matrices. ZM_hnf
returns the Hermite form; ZM_hnfall_i/ZM_hnfall return [H, U | null], selected by
withU. Their remove argument is 0, 1 or 2 with the original meaning. hnfall returns
[H, U] with the full transformation. The original column-count dispatch and
Bézout/sign/remainder schedules are retained; the permuted-Hermite entry shares
those same column operations. Inputs and returned storage are independent.

**Rationale:** Owned arrays and tuples replace native stack ownership and output
pointers. **Trade-offs:** Matrix allocation reuse, debug instrumentation and
malformed/ragged native buffers are not exposed; ragged arrays raise RangeError.
**Behavioral impact:** Compared Hermite forms and exact transformations match.
The private small-column routine is not exposed as a public TypeScript function;
forcing it directly on wide matrices can exhaust the original stack, while the
public native dispatch handles those same matrices. The matching public routes
are permanent comparative regressions. The integer-factor trace lattice uses
these general Hermite forms; the NTL factor backend remains open.


## PARI integer polynomial Newton sum adapters

**Sage/PARI:** polsym_gen computes Newton sums in characteristic zero, modulo N,
or over a coefficient quotient specified by T. It can copy an existing prefix,
including its degree/zeroth entry, and continues without verifying that prefix.
Modular sums are negated after reduction, so their stored residues are nonpositive.
Leading-coefficient inversion happens even when no new sums are requested.

**Port:** This integer-polynomial adapter takes ascending bigint coefficients and
T=null. N=null returns normalized rational pairs; a positive integer N returns
bigints with the native signed residues and unreduced zeroth entry. polsym is the
uncached characteristic-zero wrapper. Cached rational pairs are normalized as on
native construction and copied. Oversized prefixes raise an adapter RangeError
instead of writing beyond a native vector allocation. Zero polynomials, negative
counts and the covered modular inverse failures retain the original errors.

**Rationale:** Rational pairs replace scalar GEN tags, and owned arrays replace
stack-backed prefix storage. **Trade-offs:** Arbitrary GEN coefficients and the
non-null-T coefficient quotient variants are not exposed by these overloads;
those remain separate dependency work. **Behavioral impact:** All compared values,
prefix effects, signed residues and ordinary errors match. Integer factorization
uses the implemented modular route directly; exact zero remains [0,1] in the
characteristic-zero pair representation.

## PARI integer polynomial factorization adapters

**Sage/PARI:** QX_factor.c combines squarefree decomposition, deflation/inflation,
seven-good-prime selection, native finite-field splitting, Hensel lifting, bounded
recombination and van Hoeij trace lattices. ZX_factor/QX_factor omit polynomial
content from their factor matrices and retain exact factor order. ZX_gcd_all also
returns a scaled primitive quotient, rather than simply dividing the raw input by
the returned content-bearing GCD. Low-level stage helpers assume the inputs produced
by their callers; LLL_cmbf receives the unresolved group after bounded recombination.

**Port:** The same stages delegate to the existing PARI finite-field, Hensel,
adaptive LLL/norm and general HNF ports. Matrix boundaries retain their documented
dummy slots. Full factors are [polynomial, numericMultiplicity] pairs, matching the
other typed polynomial factor backends; QX_factor takes [integerCoefficients, denominator].
ZX_squff returns [squarefreePolynomials, multiplicities]. Characteristic-zero GCD
quotient output is reconstructed from the existing modular integer GCD and native
primitive/leading-coefficient normalization. chk_factors_get exposes the T=null
integer-coefficient route. Other source-private helpers are typed module exports
with their original caller preconditions. chk_factors requires at least two retained
factor columns. Inputs are preserved, and the global PARI random state advances.

**Rationale:** Typed arrays, tuples and shared kernels replace GEN pointer tags,
parallel worker storage and native allocation lifetimes. **Trade-offs:** Non-null-T
coefficient-field recombination, malformed buffers, native stack accounting and
debug timing are not exposed by these adapters. Characteristic-zero GCD keeps the
existing port's modular CRT schedule, without native worker parallelism. **Behavioral
impact:** Registered factors, multiplicities, stage order, bounds, scalar empty
products and random states match bundled PARI. Public factorization uses the original
stopping rules and lattice phase; it has no one-large-prime or unlimited-subset fallback.
Exploratory rec=0 calls to the private lattice stage can exhaust the original stack
and are recorded separately; they are not used by the native combine_factors caller
or counted as successful factorization comparisons.

## Number-field constructor integer factorization routing

**Sage:** Rational-polynomial irreducibility clears denominators, takes the primitive
integer polynomial and factors it. Degrees 30 through 300 select PARI; degrees below
30 or above 300 select NTL. The rational-polynomial method caches returned booleans,
including false, so repeat checks do not consume additional backend randomness.

**Port:** RationalPolynomial.isIrreducible delegates the inclusive 30..300 window
to PARI ZX_factor and the other degrees to NTL ZZXFactoring.factor. A single
primitive factor with multiplicity one means irreducible. Successful boolean
results are cached on the immutable polynomial; exceptions are not cached.

**Rationale:** Preserve Sage's backend degree dispatch and native algorithms.

**Trade-offs:** The constructor needs only the deterministic irreducibility result;
its NTL call uses the driver's private default stream/cache. The separately
exported pari_nf factor helpers delegate to PARI and therefore use PARI's shared
random state; they are not the constructor's degree-dispatch interface.

**Behavioral impact:** Constructors, irreducibility booleans, exact rejection text
and cached repeat calls match the Sage comparisons, including degrees above 300.
Full NTL factor order and explicit stream/cache continuation are covered by the
complete integer-driver adapter below.


## Legacy integer polynomial factor helpers

**PARI:** ZX_factor returns a matrix of primitive factors and multiplicities in
canonical polynomial order. Integer polisirreducible uses full factorization,
including for repeated polynomials, and compares the first factor's degree.

**Port:** zpFactorSquarefree normalizes signed content, delegates to ZX_factor
and returns the factor arrays. zpIsIrreducibleOverQ follows polisirreducible's
full-factorization route and advances the shared PARI random state identically.

**Rationale:** Retain the existing convenience signatures while using the native
factorization algorithm, ordering and state transitions.

**Trade-offs:** The squarefree helper omits content and multiplicities; callers
must supply a squarefree polynomial. Its zero/constant result is an empty list,
where native ZX_factor represents the zero polynomial as a zero factor.

**Behavioral impact:** Primitive factors and irreducibility results match the
native comparisons after this explicit container adaptation. There is no local
coefficient-size limit or separate local prime/subset search.


## NTL integer factorization trace adapters

**NTL behavior:** ZZXFactoring.cpp computes modular Newton traces from a supplied
prefix, chops or combines those traces, advances prime-power precision bounds
and constructs the lattice used by van Hoeij recombination. Native vectors and
matrices are mutable output arguments; vec_ZZ(i) and matrix(i,j) use one-based
indexing. Internal helpers assume valid prime powers, dimensions and reduced
trace residues. Argument errors use NTL LogicError. Verbose lattice construction
can print diagnostic information.

**Port behavior:** The corresponding module exports ComputeTrace, ChopTraces,
DenseChopTraces, Compute_pb, Compute_pdelta, BuildReductionMatrix, Compute_pb_eff
and d1_val. Vectors are zero-indexed arrays; matrices use rows, without dummy
entries. Mutated outputs are returned as copies or tuples, with unused vector
tails retained when native code retains them. BigInt implements exact arithmetic;
safe JavaScript integers represent native long dimensions and exponents. The two
native floating-point heuristics retain their original double formulas.

**Rationale and trade-offs:** Explicit values replace output references and native
memory ownership. Compute_pb returns [b,pb], Compute_pdelta [delta,pdelta],
BuildReductionMatrix [M,C], and Compute_pb_eff [b_eff,pb_eff]. Native argument
messages map to the existing JavaScript Error convention of the NTL port. Native
preconditions remain required: p is a prime at least two, moduli and chopping
powers are positive, trace residues are reduced, and matrix/vector dimensions
match their arguments. Internal helpers do not gain arbitrary invalid-input
semantics or verbose console output. Full NTL factorization is still incomplete.

**Behavioral impact:** The numeric values and native validation order match 4,289
shared comparisons, with all eight helper entries represented. Inputs remain
independent of returned vectors/matrices. The oracle compiles the eight bundled
NTL 11.6.0 bodies verbatim against Sage's installed NTL 11.4.3 primitive types and
arithmetic; it does not substitute installed factorization bodies. An NTL error
callback captures the original messages instead of allowing the installed
non-exception build to abort the process. The callback/JavaScript exception
mapping is an explicit adapter, not a claim of C++ exception-class identity.


## NTL multifactor Hensel lifting and integer products

**NTL behavior:** MultiLift normalizes word-prime factors in the active context,
requires at least two monic factors and a monic target, and returns the original
factors when the exponent is one. Otherwise it builds the same minimum-degree
factor tree, computes coprimality witnesses, lifts factors/inverses through the
ceiling-halved exponent sequence, skips the final inverse lift and restores the
input factor order. Native NTL_OVERFLOW(e,1,0) rejects exponents at 2^60 on the
64-bit build, before coprimality checks. Integer polynomial products dispatch to
classical, Karatsuba, homomorphic FFT or Schoenhage–Strassen kernels.

**Port behavior:** ZZXFactoring.MultiLift(a,f,e,p,options?) returns independent coefficient
arrays and uses an explicit word-prime modulus p. Native error messages and
validation order are preserved, with the existing NTL JavaScript Error adapter.
ZZX1.mul(a,b) and sqr(a) supply its exact integer products through the already
documented portable NTL polynomial multiplication boundary. Signed coefficients
are reduced modulo a power of two strictly larger than twice the exact product
coefficient bound, then recovered by centered representatives. The existing
bounded classical leaves, packed BigInt products and Karatsuba fallback handle
the convolution. No integer coefficient is converted to floating point.

**Rationale and trade-offs:** Arrays replace mutable native objects and output
references; p replaces the active zz_p context. e must be an exactly represented
integer number and p a valid native word prime. Explicit context settings and
shared state follow the stateful Hensel entry below; native buffer layouts,
process-global contexts and thread schedules are not reproduced. The portable product retains
the existing documented Karatsuba fallback complexity above its packing budget;
it does not claim native FFT/SS crossover or asymptotic performance. Native
validation is preserved, while resource ceilings remain JavaScript-dependent.
The factor tree and Hensel correction/inverse schedules are unchanged.

**Behavioral impact:** 4,058 shared comparisons cover full coefficient vectors,
errors, factor order, large word primes, normalization, staged exponents,
32-factor trees, signed products and packing boundaries. MultiLift's six native
bodies are compiled verbatim from bundled NTL 11.6.0 against Sage's installed
NTL 11.4.3 arithmetic types/primitives. Integer-product comparisons call those
installed native primitives; both bundled product dispatch and the portable
boundary have been read explicitly. Error callbacks capture native messages.
The complete NTL factorization caller is described below.


## NTL exact integer lattice adapters

**NTL behavior:** LLL.cpp performs exact incremental Gram–Schmidt computation,
size reduction with nearest-integer ties toward zero, independent/dependent row
swaps and optional transformation tracking. LLL_plus returns the ordered integer
Gram determinants; LLL/image return the final determinant and rank. Dependent
zero rows stay at the front of the full output matrix. LatticeSolve preserves
its output vector on failure and offers no reduction, image reduction or
LLL/image reduction of a successful solution. lip.cpp's extended GCD gives
coefficient one for the first input when both inputs are zero.

**Port behavior:** LLL, LLL_plus, image and LatticeSolve are exported from the
mirrored LLL module. Matrices use rows without dummy entries, and returned
arrays replace native output references. Internal reduction retains the native
one-based rank/index bookkeeping. LLL returns [rank,det,B,U], LLL_plus returns
[rank,D,B,U], image returns [rank,det,B,U], with U=null when unrequested.
LatticeSolve returns [success,x], including a copy of the supplied initial x on
failure. lip._ntl_gexteucl returns [s,t,d] with native signs, zero and coefficient
conventions. BigInt supplies the integer arithmetic and Euclidean primitive.

**Rationale and trade-offs:** Immutable output values replace C++ reference and
aliasing behavior. Dimensions and a/b/reduce parameters are integer numbers;
matrices must be rectangular. LatticeSolve's optional columns argument supplies
the dimension of an empty row matrix and must agree with existing rows. Native
verbose diagnostics and allocator/thread state are omitted. BigInt Euclidean
arithmetic does not reproduce GMP's limb-level half-GCD scheduling. Native error
messages map to JavaScript Error through the existing NTL adapter.

**Behavioral impact:** All 5,576 shared records match, including exact bases,
transformations, Gram determinants, image rows, lattice-solver failure outputs
and extended-GCD coefficients. Tests span zero dimensions, dependent rows,
several reduction parameters, through 16-by-16 dense bases and 4,096-bit GCD
inputs. The native oracle compiles the entire bundled NTL 11.6.0 LLL.cpp against
Sage's NTL 11.4.3 primitive arithmetic/types; direct extended-GCD comparisons use
those installed primitives, with the bundled lip.cpp conventions read explicitly.
The integer factor driver uses this exact lattice backend. Floating-point NTL
lattice variants remain separate.


## NTL factor recovery bounds and lifting adapters

**NTL behavior:** PolyEval performs exact Horner evaluation but rejects the zero
polynomial. RootBound makes all nonleading coefficients nonpositive, doubles an
integer upper bound, binary-searches the first nonnegative evaluation and scales
by the absolute leading coefficient. CutAway removes rows using ordered Gram
determinants, divides the retained leading block by C and computes its integer
image; native inputs D/M are discarded. AdditionalLifting selects the next
exponent, normalizes a nonmonic target and calls MultiLift after reducing its
old factors modulo p. The old P1 parameter is an output slot, not a checked input.

**Port behavior:** The four mirrored functions preserve the native algorithms,
returned coefficient vectors and error order. Arrays replace mutable output
references. CutAway returns its basis while retaining independent D/M inputs;
AdditionalLifting returns [new_P1,new_e1,new_factors]. Explicit word-prime p
replaces the native global context. Integer dimensions/exponents use number;
all arithmetic coefficients and bounds use bigint. AdditionalLifting retains
native word-power overflow checks and the original error strings, including
NTL's `_ntl_zexps` spelling for negative word-base exponents.

**Rationale and trade-offs:** Returning independent values fits the existing NTL
array adapters and avoids destructive vector clearing. Native verbosity and
allocator limits are omitted. Explicit context and product state follow the
stateful Hensel entry below. Caller preconditions require rectangular
matrices, enough Gram determinants/columns, valid word-prime p and exactly
represented integer parameters. The JavaScript BigInt resource limit can still
differ from native allocation limits; no infeasible-output guarantee is made.

**Behavioral impact:** All 3,848 shared comparisons match, including signed and
constant root bounds, actual factor-lattice caller inputs, signed division,
ignored old modulus values, target normalization and native early errors.
The oracle compiles bundled NTL 11.6.0 recovery/lifting bodies together with the
entire bundled LLL.cpp, using installed NTL 11.4.3 primitive arithmetic/types.
MultiLift's symbol is renamed only for linkage isolation; its body and error
strings are preserved. Error callbacks provide the existing explicit JavaScript
Error mapping. The complete integer factor driver and constructor routing use
these helpers as described below.


## NTL factor selection adapters

**NTL behavior:** ZZXFactoring.cpp uses degree-count patterns and bitsets to
prune subset reconstruction. Its subset-degree overload returns one slot per
factor, including trailing slots left from lower cardinalities. RecordPattern
assigns the last entry for each degree. ConstTermTest caches modular prefix
products, centers the final coefficient and uses NTL integer divisibility
(including 0 dividing 0). Polynomial products repeatedly combine the two factors
of lowest degree. Selection/removal assumes sorted, distinct valid indices.

**Port behavior:** Eleven mirrored names implement these helpers, including
both CalcPossibleDegrees overloads and both modular-vector mul variants.
Arrays replace native polynomial/vector containers and mutable output slots.
CalcPossibleDegrees(pat) returns a bigint; CalcPossibleDegrees(fac,k,p) returns
bigint[] with all native suffix slots. ConstTermTest returns
[integer_result,updated_prod,updated_ProdLen]. inplace_rev, RecordPattern,
SubPattern and RemoveFactors return independent arrays. mul(W,p,I?) selects
all factors when I is omitted and the supplied subset otherwise; its scratch
factor list is copied. p replaces the native modular context. Other overloads
of the native mul name are not exposed by this adapter.

**Rationale and trade-offs:** The explicit context and returned values follow
the existing NTL adapters. Caller-owned arrays are preserved even when the C++
helper mutates or destroys them. Native partial output mutations before errors
are not exposed. Patterns, indices, dimensions and counts must be exactly
represented integer numbers and their sums must fit that representation.
RecordPattern requires positive degree labels within the supplied pattern size;
ConstTermTest requires a nonempty valid index list, sufficient cache length and
a valid prefix length. p must exceed one; RecordPattern additionally requires
a valid word-prime modulus. These private helper preconditions remain explicit.
The product scheduler follows native NTL; coefficient multiplication uses the
existing documented portable NTL polynomial kernel instead of native FFT tables.

**Behavioral impact:** All 7,291 shared comparisons match, including zero and
constant polynomials, negative bit inputs, degree collisions, empty factor lists,
invalid subset cardinalities, cache reuse, positive midpoint centering and
noncanonical coefficients. The oracle extracts the bundled NTL 11.6.0 helper
bodies (including both degree and product overloads) and links installed NTL
11.4.3 primitive arithmetic. This is dependency coverage for factor recovery;
complete constructor factor routing is described below.


## NTL word polynomial quotient adapters

**NTL behavior:** lzz_pX.cpp provides prebuilt polynomial moduli, reciprocal-based
block reduction, specialized multiplication by X, reduced products/squares and
left-to-right binary powering. Negative exponents invert the final result.
InvMod/InvModStatus use the half-GCD algorithm in lzz_pX1.cpp, returning either
an inverse or (for the status form) a monic gcd. Modulus construction computes
its reciprocal only beyond the context-dependent native crossover. This affects
error timing when a composite coefficient modulus has a nonunit leading term.

**Port behavior:** lzz_pX.ts exports zz_pXModulus, build, rem, MulMod, SqrMod,
MulByXMod, InvMod, InvModStatus, PowerXMod, PowerXPlusAMod and PowerMod. Coefficient
arrays and explicit bigint p replace native zz_pX/global context state. The
constructor takes (f,p), with null f representing an uninitialized modulus;
build(F,f) replaces its copied polynomial and caches. rem/MulMod/SqrMod and
power functions use F; MulByXMod and inverse functions take raw f and p.
InvModStatus returns [status,polynomial], with status 0 for an inverse and 1 for
a gcd. Other native overloads and FFT-prime context constructors are not exposed.

**Rationale and trade-offs:** The existing NTL bigint coefficient arithmetic is
shared across word and arbitrary-precision modulus adapters. Word inverses delegate
to the word-specific extended half-GCD schedule, preserving native 90/180/350
cutoffs, matrix-transform and ordinary-product FFT checks. Ordinary `mul`/`sqr`
retain their prime-count-dependent 150/300/500 cutoffs. Final optional context
settings are accepted by products and inverses and retained by negative powering.
Exact portable products and reciprocal division replace native packed storage; native FFT storage, preconditioned scalar multipliers,
SIMD and threading are omitted. The constructor preserves the exact native
45/90/180 crossover profile for 60-bit words, FFTMaxRoot 25 and FFTFudge 4.
The first two deterministic FFT primes from bundled FFT.cpp define the context
thresholds, with native probes around both transitions. Coefficients and p use
bigint; caller inputs must fit native signed-word p, and values p <= 1 or
p >= 2^60 preserve native context errors. Composite moduli are accepted, but
operations requiring a nonunit inverse fail as in NTL. Internal arithmetic,
reciprocal and crossover fields are implementation details. Direct mutation of
F's fields between operations is outside the contract; use build to rebuild it.
Native allocator/FFT size limits are not represented as feasible-output promises.

**Behavioral impact:** All 10,880 shared comparisons match complete bundled
NTL 11.6.0 lzz_pX.cpp and lzz_pX1.cpp compiled against installed NTL 11.4.3
primitive types, scalar arithmetic and FFT support. The initial uncommitted
adapter computed its reciprocal too early, causing 135 mismatches on composite
moduli. Matching the native construction threshold repairs all of them; the
native cases and explicitly labeled prototype baseline are permanent. Tests also
cover 28 rebuild sequences, 70 exact context-transition controls, negative/zero
powers, inverse gcd results and polynomials through degree 1,024. Other standalone modular APIs remain separate; complete integer factorization
and constructor routing are described below.


## NTL word matrix multiplication adapter

**NTL behavior:** mat_lzz_p.cpp uses blocked scalar dot products at the base
case and Strassen–Winograd above a platform/thread-dependent crossover. Its
seven-product schedule handles odd dimensions using additional last-row,
last-column and inner-dimension products. Input/output aliasing is supported
through a temporary native matrix. Dimensions and coefficient context determine
native early errors.

**Port behavior:** mat_lzz_p.ts exports mul(A,B,p,columns?,inner?), with arrays
of bigint rows and an explicit word modulus. The port preserves the seven-product
schedule and all odd-dimension corrections. It uses the single-thread non-AVX
crossover 448. Base-case dot products pack canonical coefficients into integer
slots, extract the exact middle coefficient of each row/column product and reduce
modulo p. Each slot can hold every possible dot-product coefficient, so no
carry can enter an adjacent slot. Inputs and output are independent arrays.

**Rationale and trade-offs:** BigInt packing replaces native scalar/SIMD storage
while retaining exact dot products and the asymptotic Strassen recursion. Native
SIMD choices, parallel scheduling and scratch-window reuse are omitted. `columns`
defaults to the first B row's length, or zero. `inner` defaults to the first A
row's length, or B.length; explicit dimensions preserve zero-row/zero-inner
contracts. Inputs must be rectangular, with exactly represented integer dimension
arguments. A zero-row result is [], which cannot itself store its column count.
Word-modulus limits match the existing 60-bit NTL profile. Allocator limits remain
host-specific, and native speed or memory usage is not promised.

**Behavioral impact:** All 1,246 shared comparisons match the entire bundled NTL
11.6.0 mat_lzz_p.cpp compiled with installed NTL 11.4.3 types/scalar primitives.
Tests force the oracle to one thread and include empty dimensions, noncanonical
signed coefficients, composite moduli, maximum coefficient carry bounds,
rectangular products, all three sides of the 448 cutoff, even/odd Strassen splits
and a recursive 899-by-899 product. Large generated inputs use identical specified
64-bit recurrences with high/low-bit mixing on both sides. Every output entry is
compared using reversible fixed-width hexadecimal encoding, not a hash. No
production discrepancy was observed in this batch. Modular composition and the
complete integer-factorization route are described in their separate entries.


## NTL word polynomial composition adapters

**NTL behavior:** lzz_pX1.cpp builds matrices of powers and a next-block power
for modular composition. CompMod evaluates coefficient blocks by native matrix
multiplication, then combines them with a Horner step over blocks. Comp2Mod and
Comp3Mod share one argument. reduce rebuilds that argument from its original
base after reducing the modulus. GCD uses the native half-GCD/Euclidean schedule.
lzz_pXFactoring.cpp TraceMap and PowerCompose use binary composition, retaining
zero-accumulator and identity shortcuts. PowerCompose allocates scratch storage
before the zero-exponent shortcut, so an uninitialized modulus raises the native
negative-vector-length error even when the requested iteration count is zero.

**Port behavior:** lzz_pX1.ts exports zz_pXNewArgument, build, CompMod, Comp2Mod,
Comp3Mod, reduce, GCD and XGCD; lzz_pXFactoring.ts exports TraceMap and PowerCompose.
The composition argument stores bigint matrix rows in mat and a bigint polynomial
in poly. build/reduce update it only after success. Output references become
independent arrays or tuples. All modular operations use the existing explicit
zz_pXModulus context; GCD and XGCD take an explicit bigint p. XGCD returns
[d,s,t] with the native monic gcd and Bezout coefficients, including [[],[1n],[]]
for two zero inputs. Its extended recursion mutates local coefficient arrays only. CompMod accepts a base
polynomial or a prebuilt argument. Matrix products delegate to mat_lzz_p.ts.
Native word half-GCD cutoffs 90/180/350 and GCD cutoffs 400/800/1400 are retained
for the existing default 60-bit context profile.

**Rationale and trade-offs:** Row arrays and returned values follow the existing
NTL adapters. Portable polynomial products/division replace native FFT storage
and shared transform scheduling; the original block composition and half-GCD
algorithm families are preserved. Integer m/d/q parameters use exactly represented
numbers and binary loops avoid 32-bit JavaScript truncation. Native matrix-argument
width reduction, overwrite timing and error precedence remain observable.
The older zz_pXArgument representation, global argument memory bound, verbose
progress and alternate FFT context constructors are not exposed. Cache fields
must remain consistent with their coefficient context; use build/reduce to update
them. Native allocation and FFT resource ceilings remain host-specific.

**Behavioral impact:** All 10,559 shared comparisons match complete bundled NTL
11.6.0 lzz_pX.cpp, lzz_pX1.cpp, lzz_pXFactoring.cpp and mat_lzz_p.cpp compiled
with installed NTL 11.4.3 primitive types/scalar/FFT support. The initial adapter
missed PowerCompose's scratch-allocation error, producing 120 mismatches for an
uninitialized modulus. The corrected check, all reproductions, an explicitly
labeled uncommitted-prototype baseline and 45 extra precedence controls are
permanent. Other comparisons cover cache widths, cache reduction, direct and
batched composition, large iteration counts, scalar inverse errors and GCD
inputs through degree 2,803. Complete integer factorization and constructor
routing are described below.


## NTL word polynomial distinct-degree factorization adapters

**NTL behavior:** lzz_pXFactoring.cpp NewDDF uses baby/giant Frobenius steps,
batches four products per GCD, shrinks composition caches as factors are removed,
and refines giant intervals into distinct-degree groups. Its caller supplies a
monic squarefree polynomial over a prime field and X^p modulo that polynomial.
SFCanZass1 computes X^p and invokes NewDDF, rejecting constant inputs first.
Native output references hold factor/degree pairs and the Frobenius polynomial.

**Port behavior:** NewDDF(f, h, p) returns [bigint[], number][] in native order;
SFCanZass1(f, p) returns [factors, h]. Polynomial arrays are copied/normalized in
an explicit bigint coefficient context. The native baby/giant schedule, degree
threshold 1000 for tuning the baby-step count, power/composition selection and
four-entry GCD table are retained. Existing word polynomial and matrix ports
provide the arithmetic. Native thread-local scratch vectors become local arrays.

**Rationale and trade-offs:** Returned tuples and explicit p follow existing
NTL adapters and avoid global coefficient state. Portable polynomial arithmetic
replaces native FFT representations but retains the algorithm family and degree
cutovers. The default GCD table size is fixed at four; verbose progress and tuning
globals are not exposed. Native floating-point arithmetic is retained only for
its bounded baby-step tuning factor; polynomial arithmetic remains exact. The
native squarefree/prime-field caller contract is not broadened into an additional
validation layer. NewDDF groups equal-degree irreducible factors; splitting those
groups is a separate, still outstanding stage.

**Behavioral impact:** All 2,481 shared comparisons match complete bundled NTL
11.6.0 lzz_pX.cpp, lzz_pX1.cpp, lzz_pXFactoring.cpp and mat_lzz_p.cpp compiled
with installed NTL 11.4.3 primitive types/scalar/FFT support. They cover degree
1024, both baby-step strategies, tuning boundaries, ordered groups, early returns,
unreduced h and invalid coefficient contexts. No mismatch was observed in this
batch; two ownership checks additionally verify the TypeScript array contract.


## NTL word polynomial sequence adapters

**NTL behavior:** lzz_pX1.cpp MinPolySeq reconstructs the minimum polynomial of
a linearly generated sequence using its first 2*m entries; m is a degree bound.
It selects Berlekamp–Massey below the native context-dependent crossover and a
half-GCD transformation above it. For the default word profile, the crossovers
are 480, 900 and 1600. The native function checks m and sequence length before
reconstruction and preserves scalar inverse errors.

**Port behavior:** MinPolySeq(a: readonly bigint[], m: number, p: bigint) returns
bigint polynomial coefficients. It uses exactly represented integer m, explicit
coefficient context p, the native discrepancy updates and the native reversed
sequence/half-GCD construction. The existing GCD transformation is shared without
changing its schedule. Outputs are independent of a and ignore entries after 2*m.

**Rationale and trade-offs:** Explicit contexts and returned arrays follow the
existing NTL adapters. Portable polynomial arithmetic replaces FFT storage while
retaining the two original algorithm families and their crossovers. The native
recurrence-degree bound is a caller precondition; the function does not validate
that an arbitrary prefix extends to a recurrence of degree at most m. Exploratory
prefixes violating that assumption retain the selected native algorithm's output,
which can differ between algorithms. This is not a general recurrence-fitting API.

**Behavioral impact:** All 3,245 shared comparisons match complete bundled NTL
11.6.0 modules compiled with installed NTL 11.4.3 primitive types/scalar/FFT support.
They cover all crossovers, degree bounds through 3201, generated recurrences,
impulses, zero sequences, unused suffixes and 317 matching native errors. No
mismatch was observed in this batch. Existing composition/GCD comparisons remain
in place for the shared transformation, and two ownership checks cover arrays.


## NTL prime sequence adapters

**NTL behavior:** ZZ.cpp PrimeSeq generates successive small primes with a shared
base sieve and a per-instance moving sieve. reset(b) makes the next result the
smallest prime at least b. The native 60-bit profile uses NTL_PRIME_BND=16383,
ends at 32767^2=1073676289 and returns zero after exhaustion. Resetting to at most
two restarts at two, including after exhaustion. Native instances are noncopyable.

**Port behavior:** ZZ.ts exports PrimeSeq with next(): bigint and
reset(b: bigint): void. The native segmented sieve, block shifts and reset state
transitions are retained. Uint8Array replaces native character vectors; each
instance owns its moving buffer and shares the module's base sieve. Internal
indices remain exactly represented numbers below 2^31; input/output integers use
bigint. The initial sieve's integer square root is the constant 181 for this
explicit word profile.

**Rationale and trade-offs:** Explicit bigint results match the dependency's
integer adapters. A module-local base sieve replaces native process-global lazy
storage; worker isolates have independent caches. The sequence retains NTL's
finite prime bound rather than silently switching to another prime generator.
The older ZZ wrapper's unimplemented operations remain separate from this class.

**Behavioral impact:** All 1,031 shared operation traces match the five original
PrimeSeq method bodies and their base-sieve definition extracted directly from
bundled NTL 11.6.0 ZZ.cpp, compiled with installed NTL 11.4.3 types/primitives.
Comparisons include 40,000 consecutive primes, block boundaries, repeated resets,
two interleaved instances, exhaustion and restart. Every returned prime is
compared; no digest is used. No mismatch was observed in this batch.


## NTL word polynomial projection adapters

**NTL behavior:** lzz_pX.cpp caches a reduced multiplier's polynomial and FFT
representations. build writes its polynomial before rejecting an excessive degree,
leaving previous FFT caches intact on that failure. MulMod selects plain or cached
multiplication by native degree thresholds. lzz_pX1.cpp UpdateMap computes the
transpose of modular multiplication, using either shifted inner products or FFT
correlations. ProjectPowers transposes block composition through matrix products
and repeated UpdateMap calls. Vector lengths are significant: ProjectPowers checks
a's original length, while UpdateMap strips trailing zero coefficients first.

**Port behavior:** zz_pXMultiplier is exported from lzz_pX.ts, with empty and
(b, F) constructors, b, UseFFT and val(). build and MulMod accept this cached form
in addition to their existing signatures. UpdateMap(a, B, F) and
ProjectPowers(a, count, h, F) are exported from lzz_pX1.ts; h is a polynomial or a
zz_pXNewArgument. Returned vectors use bigint and retain ProjectPowers' requested
length. val returns a copy instead of the native const reference. Polynomial
composition now also uses the cached multiplier adapter at native call sites.

**Rationale and trade-offs:** Explicit contexts, returned arrays and independent
val copies follow existing NTL adapters. Portable polynomial convolutions implement
the native transposed-remainder calculation in the FFT branch, retaining its fast
algorithm family; native FFT buffers become cached multiplier and reversed-correction polynomials.
The correction is built once with the multiplier, matching native B1 setup. The
plain branch retains shifted inner products. The same cached/plain thresholds,
failed-build overwrite timing and vector validation order remain observable.
Cache fields must stay consistent with their coefficient context; use build to
update them. Native allocation/FFT resource limits remain host-specific.

**Undefined native zero-count boundary:** With a valid prepared argument and a
valid vector, ProjectPowers(count=0) indexes row zero of an empty matrix. An
isolated unmodified native process returned SIGSEGV (-11); count=1 returned the
expected vector. The TypeScript adapter instead throws Error with message
`ProjectPowers: prepared argument requires a positive count`. The direct polynomial
form still returns [] for count=0. Invalid vector lengths and empty argument caches
retain their earlier native errors. This guard intentionally replaces undefined
native behavior, rather than claiming identical output. The initial native probe
is archived in lzz_pX_projection.baseline.json. The Python oracle rebuilds an
unguarded executable and rechecks the crash/positive control once per source,
compiler and NTL header-profile cache key; a changed boundary fails verification.

**Behavioral impact:** 8,432 shared records match complete bundled NTL 11.6.0
modules compiled with installed NTL 11.4.3 primitive types/scalar/FFT support.
Another 423 records check the explicit zero-count guard in both adapters, for
8,855 total shared records. Coverage includes nonmonic moduli, composites, every
cached/plain cutoff, degree 512, projection counts through 1027 and 300 multiplier
rebuild traces. No TypeScript/native mismatch was observed within the defined
native input domain. Array ownership and the guarded native boundary have unit
checks in addition to the shared comparisons.


## NTL deterministic byte stream adapters

**NTL behavior:** Bundled ZZ.cpp supplies SHA-256, HMAC-SHA256, DeriveKey and a
RandomStream implementation selected by the native build profile. The installed
comparison profile uses the generic ChaCha20 implementation (64-bit long, 32-bit
state words, no AES/AVX2/SSSE3 backend), with a 64-bit counter and nonce. Historical
private kernel names salsa20_init/core/apply actually implement ChaCha20. Native
pointer inputs require valid storage: stream keys contain at least 32 bytes and
kernel states at least 16 words. SHA/HMAC output lengths are clamped to 0..32;
negative DeriveKey/get counts throw. HMAC passes its data length through a 32-bit
sha256_update argument, whereas standalone SHA processes long inputs in chunks.

**Port behavior:** The mirrored ZZ.ts module exposes these helper names plus
RandomStream. Uint8Array/Uint32Array replace pointer/length pairs; lengths and get
counts are exactly represented integer numbers. Output arrays contain the written
bytes rather than requiring caller allocation. sha256/hmac_sha256 default to 32
output bytes. DeriveKey uses the native HMAC counter construction. Stream keys use
the first 32 bytes, including nonzero-offset views; all returned byte arrays own
their storage. salsa20_core mutates the first 16 words; salsa20_apply returns its
16-word block and updates only the input counter. Extra input words are preserved.
Nonce bigint values wrap modulo 2^64 as native ZZ-to-unsigned-long conversion does.
Copy construction and assign preserve native stream state and alias visibility.

Short keys and short kernel states throw explicit adapter errors instead of
allowing native out-of-bounds pointer reads. These are caller-precondition guards,
not claimed native exception equivalence. Negative stream reads retain their
native error without advancing state. HMAC retains the low 32 bits of data length;
standalone SHA retains the full input length and modulo-2^64 bit count.

**Rationale and trade-offs:** Typed arrays and explicit copy assignment preserve
the native bytes and mutable state in JavaScript. Public deep-module kernel
helpers permit direct comparison of counter carry/wrap without generating hundreds
of gigabytes. The generic native algorithm is retained with fixed-size hashing
workspace. Other native acceleration profiles, unseeded entropy, process-global
SetSeed and global FFT-table random consumption are outside this batch. Native
allocation limits differ from JavaScript typed-array limits; huge allocations are
not a supported cross-runtime exception boundary.

**Behavioral impact:** All 4,038 defined native records match the exact extracted
bundled NTL 11.6.0 implementations compiled with installed NTL 11.4.3 primitive
types. The public RandomStream class is also extracted from the bundled header.
Another 28 records cover explicit short-buffer guards. The comparison includes
padding/block boundaries, HMAC key hashing, derivation counter carry, every stream
byte, chunked reads, copies, aliases, assignment, nonces and kernel counter wrap.
The build rejects incompatible acceleration/word profiles and hashes source,
compiler and installed version/configuration/machine headers. The 216 HMAC length
wrap probes supply short backing storage with an explicit native long length;
only its low 32-bit prefix is accessed. The TypeScript test adapter supplies the
same reported length on a short typed-array view, avoiding a multi-gigabyte
allocation while exercising the actual implementation's conversion.


## NTL integer sampling stream contexts

**NTL behavior:** ZZ.cpp random integer routines use a thread-local current
RandomStream. RandomWord draws eight bytes in the 64-bit profile. RandomBits_long
rejects lengths at least 64; RandomBits_ulong permits 64. RandomLen_long draws
ceil((l-1)/8) bytes and sets the top bit, while arbitrary-integer RandomLen draws
ceil(l/8). The word RandomBnd overload uses bits(bound-1); the ZZ overload uses
bits(bound). For more than three bytes, ZZ RandomBnd draws its high two bytes
first and rejects excessive prefixes before drawing the low bytes. All of these
choices determine subsequent random state, not merely the returned value.

**Port behavior:** The mirrored ZZ.ts helpers take an explicit RandomStream.
RandomWord(stream), VectorRandomWord(k, stream), RandomBits_long(l, stream),
RandomBits_ulong(l, stream), RandomLen_long(l, stream), RandomBits(l, stream),
RandomLen(l, stream), RandomBits_ZZ(l, stream), and RandomLen_ZZ(l, stream) retain
the corresponding native algorithms. VectorRandomWord returns a bigint array;
other samplers return bigint. Counts are exactly represented integer numbers.
RandomBnd(bound, stream, options?) defaults to the arbitrary-integer overload;
options.word=true selects the signed-word overload. Out-of-range word bounds
throw an explicit RangeError before sampling instead of undergoing a lossy host
conversion. Those six comparison controls are adapter guards, not native errors.

Nonpositive bit lengths and bounds at most one produce zero. Exact-length one
produces one without consuming bytes. Native word and arbitrary-integer overflow
checks and their error messages are retained. Rejected candidates, errors, clones,
aliases and copy assignment preserve the stream transitions. Arbitrary-integer
conversion uses a hexadecimal buffer conversion rather than quadratic repeated
BigInt shifts. Each returned vector owns its storage.

**Rationale and trade-offs:** Explicit streams expose deterministic native
behavior without introducing partially implemented global FFT/RNG interactions.
The RandomBnd option distinguishes native overloads erased by bigint type mapping.
No implicit global seeding or native FFT-cache side effect is promised by these
helpers. The existing ZZ wrapper's static stubs remain separate. JavaScript and
native allocation limits differ; huge allocations are not a shared exception
boundary, although the native pre-allocation bit-length overflow guard is retained.

**Behavioral impact:** All 2,626 defined native operation traces match, plus six
explicit signed-word guard traces. Exact original WordFromBytes-through-RandomBnd
bodies are extracted from bundled NTL 11.6.0 ZZ.cpp alongside its generic byte
stream implementation and public stream class. The oracle supplies the explicit
stream to LocalGetCurrentRandomStream without modifying the sampler bodies;
installed NTL 11.4.3 supplies primitive types. Comparisons include every returned
integer and all bytes consumed by explicit reads, plus 64 subsequent bytes from
each of three streams. They cover 8,193-bit requests, 4,096-bit bounds, rejection
sampling, all byte/word boundaries, errors and interleaved copy/alias/assignment.
No unexpected output or state mismatch was observed.


## NTL word polynomial element minimum polynomials

**NTL behavior:** lzz_pX1.cpp DoMinPolyMod projects successive powers through a
linear functional and calls MinPolySeq. ProbMinPolyMod draws n word-field
coefficients and returns a probable divisor of the element's minimum polynomial.
MinPolyMod verifies that candidate by composition and, when necessary, repeats
random projection after transposed multiplication to recover missing factors.
IrredPolyMod uses the fixed functional [1], assuming the element's minimum
polynomial is irreducible. The degree bound defaults to F.n and must bound the
true minimum-polynomial degree for the documented guarantees. Native randomized
routines use the current thread-local stream.

**Port behavior:** The mirrored lzz_pX1.ts exports DoMinPolyMod(g, F, m, R),
IrredPolyMod(g, F, m=F.n), and ProbMinPolyMod/MinPolyMod with overloads
(g, F, stream) and (g, F, m, stream). Polynomial/vector values are dense bigint
arrays with an explicit zz_pXModulus. Randomized routines take a RandomStream,
following the explicit integer-sampling context adapters. All native projection,
sequence reconstruction, composition, multiplier and retry steps are retained.
Invalid degree bounds precede random draws; invalid g can fail after the initial
projection vector has already consumed n field draws. Output arrays own storage.

**Rationale and trade-offs:** Explicit contexts preserve deterministic stream
comparisons without introducing partial global FFT/RNG behavior. Existing native
caller preconditions remain: the coefficient modulus is a field, g is reduced,
and m bounds the true degree. IrredPolyMod's irreducibility assumption is not
silently replaced with a different algorithm. Some comparative controls probe
native results outside degree/irreducibility assumptions; those outputs are not
claimed to satisfy the mathematical guarantee. The original randomized retry
loop has no arbitrary retry cap in the port.

**Behavioral impact:** All 11,713 native records match, including output polynomials,
exception types/messages and the following 64 bytes of the explicit stream.
The oracle compiles complete bundled NTL 11.6.0 word-polynomial/matrix modules
and extracts the original byte-stream and integer-sampling bodies. Its context
accessor points at an explicit stream established after modulus construction,
so ambient FFT initialization is outside this comparison. The installed NTL
11.4.3 public RandomStream class must exactly match the bundled class body;
the builder verifies this before compilation. Source/compiler/configuration/
machine-header hashes and generic-profile checks protect the oracle cache.

Coverage includes all elements of small monic quotient algebras, nonreduced and
nonmonic quotients, all three word-modulus cutoffs, dense degree-192 elements,
constant-element controls through degree 1601, errors after random draws, and
all retry exits. With a zero-key stream, ProbMinPolyMod(X, X^2+X+1 over F_2)
returns [1] after two bytes; MinPolyMod retries and returns [1,1,1] after four.
A degree-four quotient control recovers a degree-two minimum polynomial after
retrying, exercising verification when the result degree is below the bound.
No output or stream-state mismatch was observed.


## NTL cached word samplers and root products

**NTL behavior:** ZZ.h RandomBndGenerator caches p, a byte count, a bit mask and
a pointer to the current stream. Failed build calls leave the previous cache
untouched. Copy construction/assignment copy that stream pointer. VectorRandomBnd
returns zeros for bounds at most one, and otherwise uses a cached generator.
lzz_pX.cpp random draws n field coefficients and strips trailing zeros.
lzz_pX1.cpp BuildFromRoots constructs a monic polynomial, preserving zero and
repeated roots. Its 60-bit context uses multiplication crossovers 150/300/500:
small inputs use IterBuild, while larger inputs use a balanced product tree,
padded with zero roots to the next power of two.

**Port behavior:** ZZ.ts exports RandomBndGenerator(bound: bigint|null, stream),
a copy constructor, build(bound), next(): bigint, and assign(other): this.
The stream is explicit and copies/assignment retain its reference. The cache's
p, nb and mask fields are available; pointer storage is private. A null-bound
object is uninitialized. Calling next before successful build throws an explicit
adapter error rather than reading a native uninitialized pointer. Unrepresentable
signed-word bounds use the existing RangeError guard before sampling. Native
bound-validation errors and failed-build state are retained.

VectorRandomBnd(k, n, stream) returns a bigint array. The mirrored lzz_pX.ts random
function takes (n, p, stream), uses the native cached word sampler and normalizes
the output polynomial. Negative lengths and the 64-bit profile's n*8>=2^60 vector
resource check occur before random draws. lzz_pX1.ts BuildFromRoots(a, p) takes a
dense root vector and returns coefficients from constant upward. Empty input
returns [1]; trailing zero roots are retained. The native small-input iteration
and balanced padded product tree are preserved. Full monic coefficient arrays
replace the native tree's implicit monic coefficients and cyclic FFT correction,
using the existing exact portable multiplication boundary.

**Rationale and trade-offs:** Explicit stream contexts preserve deterministic
ownership/state behavior without claiming global FFT/RNG parity. Returned vectors
own storage instead of requiring output pointers. The full-coefficient product
tree preserves the native algebra and asymptotic structure while using portable
integer buffers. Native allocation limits and pointer preconditions do not map
exactly to JavaScript: short/uninitialized storage guards are identified separately,
and large allocation failures are not a shared exception boundary. Copying native
uninitialized generator fields is excluded from defined-input comparisons.

**Behavioral impact:** All 2,933 defined native records match, plus 406 traces with
explicit uninitialized-generator or signed-word guards. The native oracle compiles
complete bundled NTL 11.6.0 word-polynomial/matrix modules and original byte-stream
and sampler bodies. The installed NTL 11.4.3 RandomStream and cached sampler header
bodies must exactly match the bundled versions; the builder checks both. Explicit
test-context accessors supply the streams after native modulus setup. Comparisons
include roots through length 1025, every algorithm cutoff, repeated/zero roots,
composite moduli, vector resource checks, failed rebuilds, aliases and assignment
between samplers attached to different streams. Outputs, cache state, errors and
subsequent stream bytes all match. No unexpected mismatch was observed.


## NTL word factor recovery with explicit streams

**NTL behavior:** lzz_pXFactoring.cpp finds prime-field roots by randomized
power/GCD splitting. FindRoot keeps the smaller successful split; FindRoots
recurses in native split order. RootEDF converts roots into linear factors.
FindFactors uses a balanced tree of root products and modular composition to
recover factor groups. EDFSplit generates a random polynomial, computes its trace,
recovers a certified minimum polynomial, finds its roots and partitions factors.
EDF recursively refines those groups. SFCanZass2 completes distinct-degree groups;
SFCanZass performs the initial Frobenius/DDF stages as well. These algorithms
assume the documented monic, square-free, split/equal-degree and prime-field
conditions, and use the current native stream.

**Port behavior:** The mirrored lzz_pXFactoring.ts exposes those names with dense
bigint polynomial arrays, an explicit bigint coefficient modulus and a RandomStream.
Ordered roots/factors, native smaller-split selection, trace/minimum-polynomial
recovery, recursive refinement, early returns and random-draw order are retained.
FindFactors is deterministic and does not take a stream. Output arrays own their
storage. Native signed degree division uses truncation toward zero. Negative or
otherwise invalid degrees retain defined native early returns/errors rather than
being replaced with a blanket validation rule.

A zero factor degree is explicitly rejected at the native division boundary.
EDF preserves its earlier monic check; EDFSplit preserves its earlier modulus build;
SFCanZass2 preserves work and errors in earlier groups before rejecting a later
zero degree. These guards cover undefined native integer division, not matching
native exceptions. They are individually labeled in the shared case list. Native
loops retain their original probabilistic retries without an arbitrary retry cap.

**Rationale and trade-offs:** Explicit streams match the dependency's established
context adapters while preserving every observable byte consumed by these routines.
Native ambient FFT-table initialization happens before the oracle's explicit
stream is established and remains outside this contract. Private native helpers
FindFactors and EDFSplit are available through the deep module for composition
and comparison. Preconditions remain necessary for termination and mathematical
factorization guarantees. Some tests exercise defined early behavior outside
those preconditions without extending those guarantees.

**Behavioral impact:** All 4,433 defined native records match, plus 120 explicit
zero-degree guard records. The oracle compiles the original word-polynomial and
matrix modules, including the complete original factoring source in its translation
unit so private helpers are callable without rewriting their algorithms. Original
byte-stream and integer-sampler bodies and verified native header interfaces
supply the explicit stream. Tests compare factor/root order, errors and the next
64 stream bytes, including equal-degree products through degree 128, mixed groups,
small-field recursive retries and large prime-field root searches.

The integration also exposed an existing word-adapter error mismatch: shared
capital-polynomial division reported `ZZ_pX: division by zero`, whereas native
word division reports `zz_pX: division by zero`. The word boundary now translates
that error. Seventy observed mismatches and twenty passing controls are archived
and permanently compared; all ninety match after the correction. The difference
was in the error text, not the returned polynomial or consumed random bytes.


## NTL degree-dependent word contexts

**NTL behavior:** lzz_p.cpp constructs a word-modulus context with a requested
maxroot (default 25). It checks a negative request before checking p, then chooses
FFT primes until their product exceeds p squared times 2^(maxroot+4). More than
four primes is an error. The chosen prime count determines polynomial algorithm
crossovers; the effective FFT root bound is capped by the last prime's support.
ZZXFactoring.cpp requests ceil(log2(deg(f)))+1 during small-prime selection.
Polynomial modulus construction, root products, half-GCD matrix transforms and
FFT division enforce their native transform sizes, with scalar inverse errors
occurring at their original source positions.

**Port behavior:** zz_pXModulus accepts a final optional zz_pXOptions object
with maxroot?: number and optional PolynomialProductState. The same final options
argument is available on ordinary mul/sqr, InvMod/InvModStatus, GCD/XGCD,
MinPolySeq, BuildFromRoots, NewDDF, SFCanZass1 and the root/EDF/SFCanZass helpers.
Algorithms that receive F inherit its captured settings, including sequence
reconstruction inside element minimum polynomials. Omitting options preserves
the prior default. Native context data is represented by immutable captured
settings and internal PrimeCnt/MaxRoot/crossover fields; it is not ambient state.
Counts must be safe integers. Six explicitly labeled comparison records exercise
that numeric adapter guard before native context validation.

**Rationale and trade-offs:** The portable polynomial kernels do not require
native FFT tables or CRT storage. The first four deterministic FFT primes suffice
to compute the profile and determine whether more than four would be required.
Without explicit state the port reports the native too-many-primes error without
allocating additional tables. With state, the cold initialization entry below
supersedes that restriction and preserves all preceding cache/random effects.
Allocator failures for enormous requests remain outside the adapter. The original
word-context oracle starts explicit streams after context construction; the newer
cold-context oracle compares initialization too. Portable polynomial products
retain their existing representation and arithmetic implementation. Explicit
state does not establish ambient NTL context save/restore or process-global
FFT-cache/random-state equivalence.

**Behavioral impact:** 6,500 defined native comparisons and six labeled numeric
adapter guards match. Cases cover all four prime counts, exact threshold neighbors,
negative-argument priority, small FFT limits, nonunit inverse timing, long/sparse
GCD and sequence inputs, root products, nested minimum polynomials and ordered
factorization with subsequent stream bytes. An initial uncommitted prototype
omitted half-GCD matrix and FFT-remainder size checks, causing 66 discrepancies.
All are archived as prototype regressions and repaired. These are not described
as failures of the older committed API, which did not accept custom settings.
The oracle compiles the original context-constructor body and full original
polynomial/matrix modules, using verified installed primitive/context layouts
and the shared native TLS pointer. A discarded wrapper duplicated TLS storage
and crashed; those harness failures are excluded from the comparative evidence.


## NTL stateful small-prime factor selection

**NTL behavior:** ZZXFactoring.cpp retains a LocalInfoT containing a segmented
prime sequence, degree-pattern vectors and an admissibility bitset. Its static
SmallPrimeFactorization skips primes dividing the leading coefficient or yielding
a repeated modular factor, then compares distinct-degree patterns at good primes.
A single modular factor or a two-bit degree intersection proves irreducibility.
Otherwise it retains the first prime with the fewest factors, removes that prime
from the active information prefix, restores its context and completes ordered
modular factorization. The native defaults are seven initial good primes and a
maximum of fifty; the latter does not cap unsuitable-prime attempts in this stage.
Native vector shrinking retains initialized storage for subsequent growth.

**Port behavior:** The mirrored ZZXFactoring.ts exports LocalInfoT and
SmallPrimeFactorization(info,f,stream,options?). The result is null for the native
irreducibility exit, or ordered bigint coefficient arrays over info.context.p.
The context's requested maxroot is stored explicitly with p. Optional integer
InitNumPrimes/MaxNumPrimes fields replace the native tuning globals. Scalar reset,
validation order, persistent prime sequence, bitset updates, vector contents and
best-prime removal follow the source. Native SetBit retains the sign while editing
the magnitude; bit_and intersects magnitudes and produces a nonnegative result.

**Rationale and trade-offs:** This static native helper is exposed through the
dependency's deep module so the integer factorization pipeline can compose it.
Explicit streams and contexts replace ambient state. Scalar constructor defaults
(n=-1, counts=0) are explicit initialized adapter values; native uninitialized
scalar fields are not read as expected results. p and pattern accessors return
snapshots, and their setters copy the visible prefix while preserving initialized
storage beyond a shortened length. Previously unwritten prime slots are undefined
in TypeScript; native uninitialized long storage has no portable expected value.
Comparisons start with explicitly initialized hidden prime-vector storage on both
sides so every inspected scalar slot is defined. The active prefix is determined
by NumPrimes, not the vector's retained logical length. The public stream and
PrimeSeq objects preserve their existing explicit ownership rules.

The oracle compiles the complete original integer/word factoring modules, original
polynomial/matrix modules, PrimeSeq, stream and sampler bodies. It compiles the
original word-info constructor, context constructor and init entry while retaining
the installed shared TLS save/restore storage; duplicating that storage would
invalidate library calls. Context, sampler and prime-sequence interfaces are
checked against bundled headers. The first four FFT tables are warmed before
starting the explicit stream. These historical comparisons cover a warmed-cache
profile. The optional context and cold initialization comparisons below now cover
FFT initialization and its random effects using the full bundled native library.

**Behavioral impact:** All 2,175 defined native stateful transcripts match. These
include 724 exhaustively enumerated small squarefree monic polynomials, degree-512
inputs, large coefficients, more than fifty unsuitable leading-coefficient or
repeated-factor primes, tuning boundaries, vector shrink/grow sequences, saved
contexts, exhaustion/restart and subsequent stream bytes. Defined controls outside
the primitive/squarefree or consistent-state premises do not extend mathematical
factorization guarantees. Reusing arbitrary information for an unrelated polynomial
retains the original behavior rather than resetting its admissibility bitset.

An uncommitted prototype used JavaScript signed bit operations for explicitly
seeded negative degree state. It differed in 102 of 168 controls, including retained
bitsets and early exits. The native magnitude/sign behavior repairs all 102; the
full control set and observed prototype outputs are permanent comparative tests.
This is a prototype regression, not a claim that the previously committed API
already offered this internal stage.


## NTL local factor-information updates

**NTL behavior:** ZZXFactoring.cpp UpdateLocalInfo subtracts newly recovered integer
factors from the retained modular degree patterns, intersects the admissibility
bitset, advances NumFactors and recomputes suffix subset degrees. With van Hoeij
disabled it can add one more good prime, subject to a strict NumPrimes+1 comparison
against MaxNumPrimes. Unsuitable-prime attempts are not bounded by that setting.
Native backups restore the ambient word context on normal and exceptional exits.
Pattern subtractions preceding an error, including the failing negative slot,
remain in LocalInfo; NumFactors advances before suffix-count validation.

**Port behavior:** UpdateLocalInfo(info,W,factors,f,count,modulus,options?) updates
LocalInfoT in place and returns a new bigint suffix vector, or null if the native
suffix cache was unchanged. A returned empty vector is a completed recomputation.
The options van_Hoeij and MaxNumPrimes replace operation-local values of native
tuning globals, defaulting to 1 and 50. Exactly represented integer tuning values
follow native truthiness and comparisons; this routine adds no tuning validation.
The lifted modulus must be at least two, and coefficient arrays use the existing
integer/modular normalization. Temporary word contexts preserve info.context,
including its object identity, on both success and failure.

**Rationale and trade-offs:** The native suffix cache is an output parameter.
Returning its successful replacement or null lets the caller retain the previous
cache when no recomputation occurred, without introducing an ambient mutable
vector. Output-parameter contents from a failed call are not returned; retained
LocalInfo mutations remain observable. Input vectors are copied. Verbose progress
printing and its timing counters are not exposed. Optional PolynomialProductState
now preserves FFT-table creation and random draws, as described in the recombination
state entry below; omitting state retains the earlier warmed adapter.
Native mathematical inputs have consistent original degree, prime patterns and
recovered factors; defined error controls outside these premises do not extend
factorization guarantees. Native out-of-bounds accesses and unwritten scalar
storage are not treated as expected comparison values.

**Behavioral impact:** All 1,216 native stateful comparisons match, including
incremental recovery, reordered and nonmonic factors, negative admissibility
bitsets, partial pattern failures, counter/error order, context restoration,
prime limits, vector growth and finite exhaustion/restart. The traces compare
successful suffix vectors, retained state, errors, subsequent prime values and
64 following random-stream bytes. No prototype discrepancies were observed in
this stage. The native oracle compiles the complete bundled factorization modules
and the same verified primitive/context interfaces used for small-prime selection.


## NTL integer-factor cardinality search

**NTL behavior:** ZZXFactoring.cpp searches subsets of lifted modular factors in
cardinality order. CardinalitySearch applies admissible-degree, constant-term,
coefficient-bound and exact-division tests. CardinalitySearch1 additionally uses
unsigned-word first/second coefficient ratios, cached prefix products and sums,
f(1) divisibility and iterative subset-pruning tables. The filters and counters
follow a specific order; after more than two million counted operations they
refresh local prime information. Recovered factors remove selected modular
entries, while the remaining integer polynomial and selected caches persist.
FindTrueFactors starts with cardinality one, uses the optimized search above one,
and appends the final remaining polynomial without sorting the recovered factors.

**Port behavior:** The mirrored module exports CardinalitySearch and
CardinalitySearch1 returning copied [recoveredFactors, remainingPolynomial,
remainingModularFactors] tuples. FindTrueFactors returns the recovered factor
vector. LocalInfoT partial mutations remain visible, including on exceptions.
The explicit modulus replaces native ambient ZZ_p state. Operation-local options
van_Hoeij/MaxNumPrimes/MaxPrune replace tuning globals with defaults 1/50/10.
The original unsigned 64-bit wraparound, floating log/ceil pruning heuristics,
filter order, counter costs and cache invalidation schedule are preserved.
In particular, the native f(1) target and coefficient-ratio arrays are not
recomputed after a factor is removed. Product trees, centering, primitive-part
extraction, exact division and local-information updates use the NTL port's
existing kernels. Temporary contexts retain the saved explicit word context.

**Rationale and trade-offs:** These static native helpers are exposed through the
dependency's deep module for integer-factorization composition. Successful output
parameters become independent returned arrays; failed output-parameter contents
are discarded, while retained LocalInfo partial writes are still observable.
Verbose progress and timing output is not exposed. Inputs must use an integer
modulus at least two, exactly represented native integer parameters, positive-degree
lifted monic factors and consistent original degree/information for mathematical
factorization guarantees. Defined raw controls outside the factor-product premise
exercise native error and counter behavior without extending those guarantees.
Extreme allocator failures remain outside the explicit-context adapter. Earlier
comparisons preload four FFT tables before the explicit stream. The optional shared
state and cold comparisons below now cover initialization and its random draws.

The port explicitly rejects two undefined native indexing paths: zero cardinality
with an empty factor vector in CardinalitySearch, and an empty vector after the
optimized search's cardinality checks. Native defined errors retain precedence,
including basic vector-length overflow and optimized cardinality overflow. Five
shared cases label these adapter guards; their errors are not claimed as native
outputs. Other undefined indexing or uninitialized scalar storage is excluded.

**Behavioral impact:** All 2,879 defined native comparisons and five labeled guards
match. Cases cover factor order, complementary products, failed filters, coefficient
bounds, exact-division rejection, pruning tables, signed degree bitsets, preserved
contexts and parameter boundaries. Six long-search controls exercise continued
search, admissibility exit and partial-pattern failure after the native two-million
counter threshold. They compare updated factor counts and bitsets in addition to
results, errors, subsequent primes and 64 following stream bytes. No prototype
mismatches were observed. The oracle additionally compiles the complete original
ZZX1.cpp integer-polynomial dependency alongside the bundled factorization modules.


## NTL probable-prime stream contexts

**NTL behavior:** ZZ.cpp supplies signed-word and arbitrary-integer probable-prime
checks, trial-division and probability bounds, and random/error-bound-controlled
prime generators. Primality checks consume nonzero Miller–Rabin witnesses after
base two and native divisibility filters. RandomPrime uses sequential forced-odd
candidates below 256 bits and seeded candidate blocks from 256 bits. The historical
generator remains sequential; the word generator does not force candidates odd.
The block schedule consumes all final witnesses before testing them, restoring
the calling thread's random stream after worker candidate generation.

**Port behavior:** The eight deep-module helpers ComputePrimeBound, ErrBoundTest,
ProbPrime, RandomPrime, OldRandomPrime, RandomPrime_long, GenPrime and GenPrime_long
use bigint values and explicit RandomStream objects. Options replace NumTrials
(default 10, negative values clamped to zero) and err (default 80, clamped to
1–512). ProbPrime's word option selects the signed-word overload. Bound helpers
retain native floating expressions and the 64-bit GMP-limb, 60-bit small-prime,
64-bit nonce, 2^52 floating-precision profile. Modular arithmetic uses exact
BigInt multiplication and logarithmic binary powering. Candidate blocks preserve
the original order, nonce resets, seed derivation and caller-stream consumption.

**Rationale and trade-offs:** Explicit streams make dependency calls reproducible
without a process-global RNG. The native oracle runs the bundled ZZ.cpp bodies
against the verified installed primitive interfaces with one thread and generic
ChaCha20. The port runs the same candidate schedule sequentially; native thread
scheduling performance and ambient RNG state are not reproduced. Saved worker
contexts become private derived streams, retaining caller aliases and assignment
identity. Number-valued parameters must be exactly representable native signed
integers; native signed-overflow and allocator-failure inputs are excluded.
The word overload explicitly rejects bigint values outside signed 64-bit range,
which cannot be passed to the native overload. Its three shared guard cases are
labeled and are not attributed to native output.

**Behavioral impact:** All 4,621 defined native transcripts and three labeled
adapter guards match. Comparisons include primes, composites and pseudoprimes,
small-prime filters, word/bigint and 255/256-bit dispatch, byte boundaries, all
probability-bound return codes, exact argument errors, default/negative trials,
error clamps, copied/aliased/assigned streams, nonce resets and 64 subsequent
bytes from each stream. The enormous-length RandomPrime error consumes its seed
before failing; OldRandomPrime does not. No prototype mismatch was observed.
The candidate-block pseudoprime retry and 63-bit nonce-counter exhaustion are
preserved from source but are not reached by these ordinary seeded transcripts;
this validation does not claim exhaustive branch or probabilistic guarantees.

The two pre-existing ProbPrime/RandomPrime class-wrapper signatures are retained
as overloads and still delegate to the unimplemented ZZ class methods. Their
package-root re-exports select the same implementations. This adds usable
bigint/stream overloads without changing old wrapper failures.


## NTL word matrix linear algebra

**NTL behavior:** mat_lzz_p.cpp computes Gaussian elimination, image/left-kernel
bases, inverses, determinants and both left/right solves over the ambient zz_p
context. Elimination chooses the first nonzero pivot and leaves pivot rows
unnormalized. Its full output preserves the remaining columns beyond a requested
pivot width. Kernel rows retain the source multiplier/back-substitution and reverse
swap order. Strict inverse/determinant/solve routines throw at the first nonunit
pivot; relaxed routines skip such pivots and return zero status if none is usable.
Singular status calls retain their previous output matrix/vector. All these
routines use cubic Gaussian arithmetic, with packed/32-column kernels selected
at native size and modulus thresholds.

**Port behavior:** The mirrored dependency module exposes gauss, image, kernel,
relaxed_inv, inv, relaxed_determinant, determinant, relaxed_solve and solve, using
bigint rows and explicit modulus p. Options carry columns for empty matrices,
w for gauss, relax (default true for relaxed entries), left (default false,
meaning A*x=b), and previous inverse/solve output. gauss returns [rank, matrix],
relaxed_inv returns [determinantStatus, inverseOrPrevious], and solve entries
return [determinantStatus, solutionOrPrevious]. inv is the strict matrix-returning
form and throws on singularity. Native arithmetic, pivot order, kernel signs,
reverse permutations, error ordering and retained output values are preserved.

**Rationale and trade-offs:** The adapter uses the exact scalar Gaussian kernels
for all sizes, retaining their cubic arithmetic complexity while replacing packed
word accumulation and SIMD/thread scheduling with exact BigInt modular reductions.
Native block thresholds affect performance rather than the tested outputs.
Explicit modulus replaces ambient context/FFT initialization and its RNG effects;
those ambient effects are not included in these deterministic comparisons.
Output parameters become independent returned arrays, and failed output vectors
are not returned. Arrays carry rows; an empty result does not itself retain a
column-count field, so callers keep that dimension separately. Inputs must be
rectangular, with exactly represented native dimensions and 1 < p < 2^60.
Native allocator failures and unsafe dimension products are excluded. Eight
labeled guards cover moduli outside the signed-long domain, which cannot be passed
to the original constructor. All other tested invalid moduli execute its original
constructor, preserving modulus validation before dimension checks.

**Behavioral impact:** All 5,086 defined native comparisons and eight labeled
guards match. Controls include exhaustive small matrices, prime and composite
moduli, rank deficiency, missing/nonunit pivots, normalized signed coefficients,
left/right solves, retained singular outputs, and 15/16/17 and 127/128/129 native
dispatch boundaries. Zero-row matrices with enormous declared column counts return
immediately. Fifty-four native matrix-construction error controls compare ragged
rows and negative columns through each public entry. The public tri dispatcher
preserves its actual nonsquare message, "inv: nonsquare matrix". No prototype
mismatch was observed. The native oracle compiles complete bundled mat_lzz_p.cpp
plus the original context constructors against verified installed primitives.
General composite-modulus outputs follow NTL's pivot algorithm; zero status need
not be the mathematical determinant over every composite ring. Native transcripts
are stored with deterministic lossless gzip compression for unit replay; every
coefficient, error and retained input is compared after decompression.


## NTL FFT-prime cache contexts

**NTL behavior:** FFT.cpp uses a global lazy FFT-prime table and static descending
candidate-search state. IsFFTPrime searches for a primitive two-power root, applies
ProbPrime with five random witnesses when needed, then performs its special-factor
trial divisions. It reduces the returned root to the native maximum order. Failed
checks retain the root output parameter. NextFFTPrime rolls back candidate state
when an index is retried. UseFFTPrime extends the table through the requested
index; cached accesses consume no random bytes. Initialized entries contain root,
inverse-root and inverse-power-of-two tables and a floating modulus reciprocal.

**Port behavior:** FFT.ts exposes the corresponding functions with an explicit
FFTPrimeContext and RandomStream. IsFFTPrime returns [status, rootOrPrevious].
NextFFTPrime returns [modulus, root] without populating the table. UseFFTPrime
populates immutable FFTPrimeInfo entries; context.length() and context.get(index)
map the native lazy-table size and const indexed access. GetFFTPrime and
GetFFTPrimeRecip return cached values. InitFFTPrimeInfo builds the mathematical
root tables, and CalcMaxRoot preserves the 25-level cap. The fixed profile uses
64-bit signed words, 60-bit primes and the original 20,000-entry resource limit.
The source's misspelled negative-index error, "invalud FFT prime index", is retained.

**Rationale and trade-offs:** Explicit cache objects isolate the native global
state and make its stream consumption reproducible. Separate contexts can coexist;
callers share one when they need native shared-cache behavior. The port omits
word-preconditioning, SIMD and large multiplication-table metadata, which are
performance mechanisms replaced by exact BigInt arithmetic. Root tables and IEEE
reciprocal bits remain equivalent. Returned information is deeply frozen, matching
const native entries and stable pointer identity across table growth. Canonical
initialization inputs require 1 < q < 2^60 and 0 <= root < q. Number-valued indices
must be exactly represented native signed integers; previous root outputs must
fit a signed word. BigInt square roots replace the native GMP sqrt primitive.

Twenty-nine shared cases label explicit adapter guards: three unrepresentable
IsFFTPrime inputs; four CalcMaxRoot nontermination/signed-subtraction cases; fifteen
uninitialized table reads; and seven invalid initialization inputs. Guard errors
are not claimed as native outputs. Native partial allocation failures and the
uppermost successful cache index are not exercised. Explicit cache behavior does
not establish ambient Sage/NTL context fidelity across other ported modules.

**Behavioral impact:** All 3,230 defined native comparisons and 29 labeled guards
match. Traces compare results, errors, cache length, every root-table entry,
reciprocal IEEE bits and 64 subsequent bytes after each operation. They cover
cached reuse, nonce resets, extension through index 256, repeated-index rollback,
small and large two-adic candidates, raw exhaustion and nonunit initialization.
The raw NextFFTPrime(-1) control reaches native rollback/exhaustion outside normal
cache use; it is not a normal public UseFFTPrime request. A temporary native probe
found raw key counter 1,699,152, which reaches special-factor rejection of 2047
after the probabilistic checks. That key and its neighbors were then compared
against the unmodified original. No prototype discrepancy was observed.

The oracle compiles bundled FFT-prime and ZZ stream/primality bodies against
verified installed interfaces, including the unchanged FFT header and native
table-initialization primitive. A fresh native child process per trace preserves
original global/static initialization while retaining shared cache state within
the trace. This avoids replacing the native cache algorithm with a test model.


## NTL integer CRT and modular determinant adapters

**NTL behavior:** ZZ.cpp and mat_ZZ.cpp implement incremental CRT with the
interval -a/2 < g <= a/2, normalizing initially unbalanced residues and retaining
the positive half on ties. Even new moduli require a sign adjustment when the old
residue is positive and the correction equals p/2. Matrix CRT computes its shared
modular inverse once. DetBound multiplies bounds from every row norm, using
floor(sqrt(norm)) + 1 for norms above one, then returns the product's bit length.
The arbitrary-modulus matrix determinant uses the first nonzero pivot, delayed
integer reductions and native modular inversion, which may fail over composites.

**Port behavior:** ZZ.ts exposes CRTInRange(g, a) and CRT(g, a, G, p, {word?}).
The latter returns [modified, residue, modulus]. mat_ZZ.ts exposes the matrix CRT
with the same status tuple and DetBound; mat_ZZ_p.ts exposes determinant. Matrices
are arrays of rows with BigInt entries, copied on input and output. Explicit
columns and residueColumns options represent zero-row shapes. A modulus argument
replaces the ambient native modular context, preserving constructor validation.
The source's word and big-integer CRT schedules coincide on their defined inputs;
word inputs at or above 2^60 follow the big-integer result. GMP modular inverse,
integer square root and arithmetic primitives use exact BigInt Euclidean/Newton
operations. Modular determinant retains cubic elimination with delayed reductions;
thread scheduling and native capacity preallocation have no port counterpart.

**Rationale and trade-offs:** Return tuples replace C++ output references, and
plain arrays replace native matrix containers. Input arrays stay unchanged, and
output rows are independently owned. Context and dimension errors follow native
construction and call order. The scalar CRT domain is a > 0, p > 1, 0 <= G < p;
noncoprime moduli retain the native inverse error. Matrix residues are normalized
by modular construction, and matrix CRT requires a positive accumulated modulus.
The twenty-two labeled guards comprise six invalid accumulated moduli, six invalid
scalar new moduli, four noncanonical scalar residues, four unrepresentable scalar
word inputs and two unrepresentable matrix word moduli. These guards are explicit
adapter behavior, not assertions about undefined native arithmetic. Other tested
modulus, dimension and inverse failures execute the original native code.

**Behavioral impact:** All 8,881 defined native comparisons and 22 labeled guards
match. Tests include exhaustive small matrices/residues, half ties, native limb
boundaries, signed coefficients, empty shapes, nonunit pivots and delayed-reduction
scheduling thresholds. Determinant over a composite modulus can throw even when a
later pivot would be invertible; the port preserves that behavior. DetBound([])
is one and a zero row makes the bound zero. No new native behavioral discrepancy
was observed. These primitives do not yet complete integer determinant/inverse
reconstruction or integer polynomial factorization routing.


## NTL integer matrix reconstruction contexts

**NTL behavior:** mat_ZZ.cpp computes integer determinants and inverses with
successive FFT-prime reductions, incremental CRT and a Hadamard stopping bound.
The default schedule can confirm stabilized residues with generated probable
primes. The deterministic flag selects the original full-bound schedule. Inverse
reconstruction stabilizes both the determinant and adjugate, then checks the exact
matrix product x*A = d*I. A later determinant change rescales the certified output.
Status inverse retains its previous matrix output on singularity; empty input
clears it. Strict inverse accepts determinants 1 and -1, negating the adjugate for
-1, and otherwise throws. Native word and big modular contexts are saved/restored.

**Port behavior:** mat_ZZ.ts exposes determinant(A, FFTPrimeContext, RandomStream,
options?) and inv(A, FFTPrimeContext, RandomStream, options?). Determinant options
include columns and deterministic. Strict inverse returns a matrix and accepts
columns; the status:true overload returns [determinant, adjugateOrPrevious] and
also accepts previous and deterministic. The native probabilistic default remains
unchanged. FFT-prime cache and random stream objects are explicit and shared by
calls. Modular arithmetic contexts are explicit modulus parameters internally,
so there is no ambient port modulus to modify. Matrices are arrays of BigInt rows.
Native cubic exact certification and the original CRT/prime-generation schedules
are preserved. Returned rows are independently owned.

**Rationale and trade-offs:** Explicit state objects replace native globals, and
return tuples replace output references. Prior status output is copied; its
zero-row column metadata is not embedded in returned arrays. Column options
represent input shapes. Context restoration is asserted in the native oracle;
the port's explicit modulus arguments avoid ambient state mutation. This does
not establish ambient-context equivalence for other modules. The strict inverse
has no deterministic option, matching the original overload; the status overload
provides the native deterministic path. Native storage allocation and threading
are replaced by ordinary arrays and exact BigInt arithmetic.

**Behavioral impact:** All 3,595 defined native comparisons match full outputs,
errors, FFT cache root/reciprocal tables and following stream bytes. They cover
signed/singular/unimodular matrices, modulus divisibility, matrix scheduling
thresholds, nonzero singular rows with large bounds, rank-one shears, corrected
probabilistic checks, explicit nonce/cache state and native context restoration.
Three valid seed-adapted inputs also reproduce the original Monte Carlo failure
outcomes: a default determinant returns a stabilized residue, and default inverse
calls can classify an invertible matrix as singular. Deterministic calls recover
the full outputs. These are original seeded behaviors, not new port discrepancies;
expected values come from unmodified bundled native routines. The conditional
post-certificate rescaling path is source-reviewed but remains unexercised, as does
its defensive inexact-division error. No behavioral prototype mismatch was observed.


## NTL certified integer elimination

**NTL behavior:** The internal gauss helper in ZZXFactoring.cpp accepts a nonempty
integer matrix whose rows are linearly independent. It samples a 60-bit prime,
initializes its word-modulus context, eliminates modulo that prime and selects
pivot columns. An integer inverse of those original columns scales the input to
R. It retries on modular rank loss, a zero inverse result or failure of its exact
row-form checks. Returned d is the signed determinant of the selected columns,
and R/d is reduced row echelon form. The original restores the prior word context.

**Port behavior:** ZZXFactoring.ts exposes gauss(M, FFTPrimeContext, RandomStream,
{columns?}), returning [d, R] with BigInt row arrays. It delegates sampling, modular
elimination and integer inversion to the corresponding NTL ports and preserves
all retries and their random-state order. The observable word-context setup calls
UseFFTPrime until the product exceeds p^2 * 2^29, as the original 25-level profile
requires. Word preconditioner metadata is unused by exact BigInt arithmetic and
is omitted. The explicit modulus avoids ambient context mutation. Input arrays
remain unchanged and returned rows are independently owned.

**Rationale and trade-offs:** Explicit cache/stream objects and tuple outputs
follow the existing dependency adapters. Column options preserve zero-row shapes.
The native independent-row precondition is retained: dependent nonempty inputs
are outside the supported domain and can retry indefinitely, as in NTL. Empty
inputs retain the original internal-error message. No rank precheck, arbitrary
retry limit or substituted rational elimination algorithm is introduced. Native
context restoration is asserted in the oracle; this does not establish ambient
state fidelity across other ported modules.

**Behavioral impact:** All 3,250 defined native comparisons match full results,
errors, FFT root/reciprocal tables and following stream bytes. Tests include
exhaustive full-rank small matrices, signed pivot scaling, rectangular shapes,
large coefficients, cache/nonce/context controls and native construction errors.
Three seed-adapted controls exercise modular-rank, exact-form and inverse-zero
retries. A temporary native diagnostic confirms counts [1,0,0], [0,0,1] and
[0,1,0] respectively; its outputs equal the unmodified oracle. Permanent expected
results come only from unmodified bundled routines. No prototype discrepancy was
observed. Defensive checks guaranteed by the exact inverse certificate are
source-reviewed; line execution does not imply every branch was taken.


## NTL stateful modular polynomial products

**NTL behavior:** ZZ_pX.cpp selects multiplication and alias-square algorithms using
coefficient limb counts, degree and ZZX1.cpp's SSRatio. FFT products lazily install
ZZ_p.cpp's modulus-specific FFT information, extending the global prime cache until
the prime product exceeds p squared times 2^29. Searching for those primes consumes
random bytes. Classical, Karatsuba and Schoenhage–Strassen products do not initialize
this cache. Squaring and distinct equal operands have different crossover thresholds.

**Port behavior:** ZZ_pX.mul and ZZXFactoring.mul accept an optional explicit
PolynomialProductState containing FFTPrimeContext and RandomStream. Stateful calls
follow the original dispatch and cache initialization before computing the exact
coefficient product with the existing portable BigInt kernel. Product-tree operands
remain distinct copies, even when selected indices repeat. SSRatio uses exact integer
intermediates followed by the original floating division. Its Number arguments must
be exactly represented signed native integers. Undefined dimension, bound and padded
bound overflows are guarded; the native NextPowerOfTwo overflow error is retained.

**Rationale and trade-offs:** Explicit state enables reproducible shared-cache behavior
without process-wide mutable globals. Cached prime reuse needs no additional randomness.
The existing three-argument coefficient-only multiplication API remains available.
Native preconditioners, FFT buffers, allocation behavior, threading and asymptotic
performance are still replaced by the documented portable product implementation.
Other polynomial operations and ambient Sage callers do not yet share this state.
SSRatio guards define errors outside the original's defined native arithmetic domain;
seven shared cases identify these adapter guards separately from native results.

**Behavioral impact:** The first 663 traces exposed 99 state discrepancies in the old
coefficient-only calls. The stateful implementation matches those original outputs,
including every prime/root/reciprocal table and the next 64 random bytes after every
operation. Additional boundary cases cover native overflow, seven labeled guards and
executed documentation examples. Large-modulus tests exercise both FFT and SS routes;
rare high-degree crossover combinations and native resource/allocation exhaustion
are not claimed as covered. This is not a claim of ambient state fidelity for all
polynomial operations.

The oracle builds the complete unmodified bundled NTL 11.6.0 static library in a
content-addressed temporary directory and passes its own make check suite. It uses
the fixed 64-bit, non-threaded, GMP/gf2x, generic ChaCha profile used by the earlier
oracles. Including the original ZZXFactoring.cpp exposes its static product-tree
helper while linking all dependencies from that same build. No installed NTL globals
or TLS contexts are mixed into it. A pilot matches 3,247 earlier certified-elimination
traces, plus four smoke traces; three raw static NextFFTPrime controls are excluded
from that pilot because the full library does not export that internal symbol.


## NTL stateful integer polynomial division

**NTL behavior:** ZZX1.cpp dispatches exact division to PlainDivide when the divisor
degree or quotient-degree bound is at most eight, otherwise to HomDivide. PlainDivide
checks polynomial degrees before content, so a zero numerator and a nonconstant
divisor report failure. HomDivide handles zero numerators first and succeeds with
zero quotient. Both retain the previous output on failure. HomDivide initializes
each FFT-prime context before rejecting primes dividing the divisor's leading
coefficient, reconstructs the quotient with balanced CRT, and certifies stabilization
using the original coefficient bound. A later residue can undo early stabilization.

**Port behavior:** The three same-named array adapters return [status, quotientOrPrevious].
The ordinary divide also accepts a bigint scalar divisor. Optional previous arrays
are copied, and optional PolynomialProductState supplies the shared FFTPrimeContext
and RandomStream. The existing internal null/quotient adapter delegates to this native
dispatch and accepts the same optional state. It no longer treats zero divided by a
nonconstant polynomial as success. Prime initialization precedes the leading-coefficient
check, and native reconstruction retry/rejection order is preserved.

**Rationale and trade-offs:** Explicit owned results replace mutable output parameters;
explicit state isolates process-global NTL contexts. Without state, the historical
coefficient-only deterministic FFT-prime list remains in use, including in existing
GCD and squarefree callers. Modular polynomial arithmetic uses the existing exact
portable kernels instead of native word/FFT buffers. Native allocation limits,
threading and performance are not reproduced. This does not establish ambient Sage
state equivalence or finish factor recombination.

**Behavioral impact:** Comparing the original helper's first call found 181 value/status
discrepancies and 172 prime/cache/stream discrepancies across the initial 2,865 traces.
The final 2,870 defined native comparisons cover zero/scalar inputs, native degree
crossovers, signed content, failed divisibility, retained outputs, large coefficients,
cache warming, nonce changes and prior modular contexts. Complete root/reciprocal
and inverse-power tables and the next 64 bytes are compared after every operation.
The oracle asserts that division restores both prior modular context pointers.

A temporary diagnostic native build confirms a skipped leading-coefficient prime,
reconstruction destabilization after apparent CRT stability, a nonzero low-degree
verification remainder, a modular division failure, a nonunit plain quotient step
and coefficient-bound rejection. Its full outputs equal the unmodified original;
only the unmodified outputs are permanent expected results. The complete bundled
NTL static library supplies all native dependencies. Resource allocation failures
are not claimed as covered.


## NTL van Hoeij factor acceptance

**NTL behavior:** ZZXFactoring.cpp's internal GotThem certifies row reduction, checks
that each column belongs to exactly one factor group, rejects groups of at most three
lifted factors, then checks constant terms and coefficient bounds. It multiplies and
makes primitive all but the final group, divides the input polynomial exactly by each
candidate and appends accepted factors plus the remaining quotient to the output prefix.
Failure retains the prefix. Trial-division failure unconditionally prints X and a newline,
even with verbose disabled. The degree accumulator includes every lifted polynomial
for every row, so all degree keys are equal and its sort leaves row order unchanged.

**Port behavior:** GotThem returns [status, factors], owning all output arrays. The
explicit FFTPrimeContext and RandomStream are shared by certified elimination,
modular products and integer polynomial division. It retains native factor order,
constant-term balancing, bounds before primitive-part removal, the last-factor special
case and console.error('X') on trial-division failure. The equal-key sort is omitted
because it is a no-op for all defined input sizes. The bound must be an exactly
represented signed native integer; two shared guard cases cover out-of-range bounds.

**Rationale and trade-offs:** Rows and owned coefficient arrays replace native matrix
and mutable vector objects. Nonempty independent matrix rows and enough lifted factors
for the columns are required by the original algorithm. Dependent matrices can cause
native certified elimination to retry indefinitely and are not supplied. Empty matrix
errors are retained. Verbose timing/logging mode is not exposed. Existing portable
polynomial kernels replace machine-specific FFT buffers, with explicit observable
cache/stream state. Allocation limits and platform-specific performance can differ.

**Behavioral impact:** The complete bundled NTL library supplies the unmodified
acceptance helper and dependencies. The initial 1,657 traces match all outputs/errors,
cache root/reciprocal/inverse-power tables, following 64 bytes and diagnostic bytes.
They include 805 accepted and 2,452 rejected calls, plus 18 empty-matrix errors, with
280 X-newline diagnostics. Further boundary tests cover invalid moduli, word bounds,
ignored extra lifted factors, retained prefixes and executed API examples. Unequal-
degree groups explicitly preserve native order. Degree-99/100/101 lifted polynomials
exercise the modular FFT crossover while sharing state through integer division.
This helper does not complete van Hoeij iteration or the overall factorization audit.


## NTL cold word-context initialization

**NTL behavior:** lzz_p.cpp calls UseFFTPrime for each prime required to exceed
p squared times 2^(maxroot+4), then rejects contexts needing more than four primes.
Generating missing prime tables consumes the current RandomStream. In integer
small-prime selection this occurs before randomized factor recovery and can change
both factor order and subsequent bytes. Reused tables consume no new random draws.

**Port behavior:** zz_pXOptions accepts optional state: PolynomialProductState
(the existing shared FFTPrimeContext/RandomStream pair). The constructor and nested
word-polynomial algorithms retain that state. SmallPrimeFactorizationOptions accepts
an optional context and combines it with its existing stream. Failed initialization
retains the native cache prefix and stream advancement. Captured settings are frozen;
the explicitly shared cache and stream remain mutable. When a word-factorization
entry also takes a stream, callers must pass that same stream in state.

**Rationale and trade-offs:** Explicit ownership replaces native global state and
allows independent jobs. Omitting state preserves the earlier mathematical/warmed
adapter, whose fixed first four primes need no random draws. These are opt-in native
state effects; this change does not establish global Sage random-state equivalence
or wire the complete integer-factorization pipeline. Huge allocator/exhaustion
behavior is not emulated. maxroot retains its documented safe-integer guard.

**Behavioral impact:** 3,331 comparisons against the complete unmodified bundled
NTL match outputs, errors, retained LocalInfo data, all FFT-prime tables and following
64 stream bytes. The previous implementation differs in 2,786 trace cache/stream
results, including 864 traces with factor-output or retained-information differences.
The test set replays all 2,175 earlier small-prime scenarios from cold state, adds
warm/cold seeded factor-order controls, all context prime counts, error precedence,
requests up to maxroot=1024, FFT construction crossovers and cache reuse after
failure. Native uninitialized long-vector slots are initialized to the documented
sentinel before comparison. Full-library fresh child processes replace the older
mixed-library warmed oracle for these new tests. The old warm tests remain valid.


## NTL factor recombination cache state

**NTL behavior:** Local factor-information updates initialize word contexts both
when subtracting recovered factor patterns and when adding a new prime. Cardinality
search multiplies selected or complementary modular factors and certifies integer
divisibility. These nested operations share the original global FFT cache and
RandomStream, so a later randomized operation observes all preceding draws.

**Port behavior:** UpdateLocalInfoOptions accepts optional state:
PolynomialProductState; FactorRecombinationOptions inherits it. Both cardinality
searches pass this state to modular products, complementary products, exact integer
division and local information updates. InvMul accepts the same state as its fourth
argument. FindTrueFactors forwards its options throughout the search. Existing
return shapes, partial LocalInfo updates and owned input/output arrays are retained.

**Rationale and trade-offs:** The optional explicit cache/stream pair allows this
native helper pipeline to compose without process-global mutable state. Omitting
state preserves the earlier warmed mathematical adapter. This does not model
verbose progress timing, allocation failures, or ambient Sage state. Complete
Hensel lifting and the van Hoeij driver use this interface as described below.

**Behavioral impact:** All 4,344 shared comparisons match: 4,339 defined native
traces and five labeled guards for undefined native empty-vector indexing. The
complete bundled library runs in a fresh child for each trace. Every comparison
checks outputs/errors, retained LocalInfo data, full FFT-prime information and
following 64 stream bytes. The initial 4,252 traces exposed 112 discrepancies in
the earlier implementation, including 80 with different later outputs or retained
information. All are permanently retained. Coverage includes cold and warmed
updates, later factorization order, cache reuse, degree-99/100/101 modular-product
crossovers, selected/complementary products, integer-division initialization and
coefficient-bound rejection before division. Existing mathematical domain and
explicit native-indexing guard limitations still apply.


## NTL stateful integer polynomial products

**NTL behavior:** ZZX1.cpp chooses plain or Karatsuba arithmetic for smaller
products and then calls ChooseSS to select Schönhage–Strassen or homomorphic
multiplication. The latter's NewFastCRTHelper initializes FFT primes until their
product's bit length exceeds the coefficient bound plus two certification bits.
Integer SS multiplication uses a separate Fermat-modulus transform and consumes
no FFT-prime random draws. Shared-input multiplication dispatches to squaring.

**Port behavior:** ZZX1.mul and sqr accept optional PolynomialProductState. They
preserve the native initialization dispatch, shared alias handling and cache reuse.
The existing exact coefficient kernel remains in use. ChooseSS is exported with
the native name and four integer-number parameters so its complete decision chain
can also be compared without allocating polynomials of extreme dimensions. The
predicate delegates padding to SSRatio and checks native coefficient-size overflow.

**Rationale and trade-offs:** Explicit cache/stream ownership replaces global state;
omitting it preserves the earlier coefficient-only behavior. The GMP 64-bit limb,
60-bit word and single-thread native profile determines dispatch thresholds. The
existing portable coefficient arithmetic still differs from NTL's machine-specific
FFT buffers and performance. Verbose timing, allocation failures, Hensel state and
ambient Sage random-state composition remain outside this change. Seven shared
cases explicitly guard undefined native parameter arithmetic; non-integer or
nonfinite Number inputs are rejected before native-word interpretation.

**Behavioral impact:** 1,497 shared traces comprise 1,490 defined native cases and
seven labeled arithmetic guards. Complete coefficient/predicate results, errors,
FFT primes/roots/reciprocals/inverse-power tables and following 64 bytes match the
unmodified complete bundled library. The initial 540 product traces expose 69
previous cache/random-byte discrepancies and no coefficient or error discrepancies.
Both native context pointers are asserted restored after each product or square.
The additional 957 predicate controls target all five SS criteria, nonselection,
signed parameters and boundary arithmetic, including both intermediate additions
in the native limb-count expression. All multiplication/squaring dispatch
branches and all SS criteria are checked with a separate diagnostic build whose
outputs are compared back to the unmodified original; diagnostic counters are not
used as permanent expected values. The test set includes scalar/zero products,
trailing zeros, signed coefficients, limb/length thresholds, sparse inputs through
degree 4,096, cold alias squaring, copied equal operands and cache/nonce reuse.


## NTL word modulus rebuild state

**Sage/NTL behavior:** `build(zz_pXModulus, f)` writes public f/n before validation,
then updates FRep before computing the inverse and HRep. A failed inverse retains
HRep. A failed multiplier build can retain B1 after updating B2, or retain both
cached representations after updating its public polynomial. Later operations use
these exact cached representations and native plain/FFT thresholds. The native
truncated inverse can produce root-dependent coefficients when stale caches exceed
the expected product degree; this is deterministic and does not access uninitialized
transform entries.

**Port behavior:** private caches retain the independently updated representations,
their transform sizes and their coefficient equivalents. Failed operations preserve
native overwrite order, errors and later cyclic reductions. Raw multiplication and
squaring retain native transform checks even if coefficients cancel modulo p. The
cached-multiplier plain fallback also validates ordinary FFT multiplication before
reduction when a failed build leaves an oversized public b with a plain flag.
Products exceeding the initialized inverse-transform prefix delegate to the portable
NTL forward/inverse TFT and the original nearest-integer CRT reconstruction. Optional
PolynomialProductState supplies the actual native FFT roots. Without explicit state,
this exceptional interpolation uses a private FFT context with a 32-byte zero key.

**Rationale:** clearing a reciprocal after failure changes subsequent results and
errors. Coefficient arithmetic alone cannot reproduce an incomplete inverse transform.
The private zero-key context gives the existing stateless adapter a reproducible
choice where native behavior depends on ambient random state.

**Trade-offs:** coefficient caches replace packed FFT buffers and machine-word
optimizations. Explicit state is required to match a caller's particular ambient
root choices after invalid rebuilds; the stateless zero-key convention cannot infer
another native process's random state. Native transform-size and stale-cache errors
are preserved. A failed constant-modulus rebuild retains n=0 and earlier FFT
selection; large reduction then preserves the native negative-buffer-length error
instead of entering a loop with negative capacity. Mutating internal fields directly is not a supported cache operation.

**Behavioral impact:** 5,210 fixture traces, represented by 5,186 distinct shared
rows, cover rebuild failures, later reductions,
products and recovery, independent cache sizes, degree cancellation, zero-key fallback,
and direct transform primitives. Defined native outputs/errors and supplied random/
FFT state match. No differences are introduced for valid quotient arithmetic.

The word multiplier's allocated FFT-prime count also survives plain and failed
rebuilds. A different count raises `fftRep: inconsistent use` while resizing B2,
before overwriting either retained representation. Same-count rebuilds may retain
coefficients larger than the current coefficient modulus when an earlier quotient
is used again. Complete products normalize those coefficients before entering the
packed BigInt kernel; an incomplete inverse transform uses their original FFT-prime
residues instead. Normalizing a retained representation before that interpolation
changes the result. The additional 24.45.1 comparisons cover 740 word traces,
including native error timing, recovery, projection controls and one-to-four-prime
numeric boundaries. They fix 200 missing-error traces and 28 coefficient discrepancies.

## NTL portable truncated transforms

**Sage/NTL behavior:** FFTFwd_trunc and FFTRev1_trunc use the recursive truncated
transforms in FFT.cpp, with bit-reversed frequency order, admissible prefix lengths,
and normalized inverse results. FFT_impl.h rounds prefixes for the generic
NTL_FFT_RDUP=4 profile. Lazy modular butterflies and precomputed word tables optimize
the implementation; malformed buffer dimensions are native precondition violations.

**Port behavior:** FFT.ts exports array-returning FFTFwd_trunc(a,k,info,yn,xn) and
FFTRev1_trunc(a,k,info,yn). FFT_impl.ts exports FFTRoundUp(xn,k). The forward and two
inverse recurrences mirror new_fft_short/new_ifft_short1/new_ifft_short2, with exact
BigInt modular butterflies and O(N log N) work. Inputs are normalized modulo info.q,
outputs are independent arrays of length yn. Invalid root exponents, inadmissible
lengths, short buffers and unsafe numeric dimensions raise labeled RangeErrors.

**Rationale:** these primitives are necessary to preserve deterministic cache-failure
behavior, including incomplete inverse transforms. BigInt replaces lazy machine-word
arithmetic while retaining the original recursion and mathematical results.

**Trade-offs:** no SIMD, lazy-word storage or process-global multiplier tables.
The adapter owns arrays and explicitly supplies prime information. Noninteger or
unsafe dimensions receive defined guards rather than invoking native undefined
behavior. FFTRoundUp uses the bundled generic rounding profile, not an AVX profile.

**Behavioral impact:** direct comparisons cover four native FFT primes, independent
forward/inverse prefixes, arbitrary inverse inputs, signed coefficients and transforms
through 32,768 points. Prefix rounding covers both native rounding regimes and signed
word boundaries. CRT's binary64 estimate intentionally mirrors the native algorithm;
polynomial coefficients and transform butterflies use BigInt.


## NTL word CRT reconstruction

**Sage/NTL behavior:** the vector FromModularRep in lzz_pX.cpp reconstructs
coefficients from one through four FFT-prime residue rows. One prime uses exact
integer centering at floor(q/2). Multiple primes use a nearest-integer CRT estimate;
the bundled AArch64 clang build contracts each update into binary64 FMADD before
adding 0.5 and truncating. This profile is confirmed in the native object code.

**Port behavior:** FromModularRep(residues,F) returns independent bigint coefficient
columns, with one row per F.NumPrimes. It normalizes signed residue inputs modulo
their FFT primes, preserves zero columns, and validates the rectangular shape.
One-prime centering uses BigInt. The multi-prime estimate emulates each fused update
by aligning exact dyadic integers and rounding once to binary64. The internal helper
uses only the nonnegative normal/zero operands and bounded sums produced by CRT;
it is not exposed as a general-purpose floating-point function. Truncated polynomial
products delegate to this reconstruction helper. In FFT-prime coefficient contexts,
it implements native FromfftRep's direct residue copy instead: the native static
FromModularRep CRT helper is not valid for that context.

**Rationale:** ordinary JavaScript multiply-plus-add rounds twice. The approximate
path also differs from native exact one-prime centering near half a prime. These
errors become observable after truncated interpolation with stale transform caches.

**Trade-offs:** the explicit matrix adapter replaces native packed FFT buffers.
Malformed shapes receive defined RangeErrors rather than native buffer misuse.
The fused estimate follows the bundled AArch64 artifact; an unfused native compiler
profile may differ for residues outside the coefficient certification bound. Native
fixture generation rejects a non-AArch64 host for this oracle. Arrays and coefficient
arithmetic remain BigInt; only the native estimate has binary64 semantics.

**Behavioral impact:** 232 additional shared traces cover one through four primes,
empty columns, signed/random residues, exact half-prime and CRT half-product
neighborhoods, zero terms within a fused sum and shape guards. The first 224 expose
84 differing traces and 626 differing coefficients in the earlier formula. All new
native results and the earlier 4,978 results match after the fix.

## NTL portable transposed transforms

**Sage/NTL behavior:** FFTFwd_trans and FFTRev1_trans apply the transposes of
NTL's forward and normalized inverse transforms. The flipped recurrences invert
evaluation points while retaining bit-reversed frequency order. The forward
transpose has no scaling; the inverse transpose scales by 2^-k.

**Port behavior:** FFT.ts exposes FFTFwd_trans(a,k,info) and FFTRev1_trans(a,k,info),
returning independent full-length bigint arrays. The existing divide-and-conquer
kernels run with exchanged root tables and native scaling. Root powers are cached
by prime-info identity and extended lazily; separate seeded roots never share tables.
Invalid exponents and short input buffers receive the same labeled adapter guards
as the truncated primitives.

**Rationale:** transposed modular multiplication requires these precise frequency
operations, including after a quotient rebuild. Reusing the native FFT recurrences
retains O(N log N) work instead of introducing a slower coefficient convolution.

**Trade-offs:** BigInt butterflies replace SIMD and lazy word arithmetic. Array
ownership and explicit FFTPrimeInfo replace pointer aliasing and global tables.
Malformed buffers are guarded rather than passed to unchecked native storage.

**Behavioral impact:** direct native comparisons cover all four FFT primes, arbitrary
signed inputs, zero/single-coefficient vectors and powers of two through 32,768.
Outputs and supplied random/FFT state match the bundled generic NTL profile.

## NTL word projection rebuild state

**Sage/NTL behavior:** UpdateMap chooses its path from B.UseFFT, independently of
the modulus's current plain/FFT flag. It uses retained FRep, B1 and B2, along with
the independently retained transform sizes. A plain modulus rebuild can change n
while leaving these transforms allocated. Failed builds can leave mismatched sizes.
The plain path's LeftShift accepts -1 when the public multiplier degree equals n;
a greater degree writes before the output vector in the bundled unchecked build.

**Port behavior:** UpdateMap follows the original RevTofftRep, pointwise products,
RevFromfftRep and AddExpand sequence, including intermediate CRT reconstruction,
coefficient negation, cyclic shifts and frequency-prefix expansion. Private retained
coefficient arrays lazily reconstruct full forward representations using their exact
seeded prime roots. The internal _zz_pX_transform_kernels(F,B) adapter shares these
operations with lzz_pX1.ts. The plain path preserves the defined right-shift boundary.
After native input-length validation, an oversized plain multiplier receives
RangeError("UpdateMap: multiplier degree exceeds native buffer bounds").

**Rationale:** current public polynomial coefficients cannot replace an older cached
transform. Native error order and intermediate coefficients remain observable even
with empty projection vectors. The guard avoids executing native buffer corruption.

**Trade-offs:** explicit arrays and private WeakMap caches replace packed native
FFT buffers. Direct mutation of modulus/multiplier internals or prime tables remains
outside the adapter contract. Without explicit PolynomialProductState, the existing
private zero-key convention supplies roots for exceptional stale-cache behavior.

**Behavioral impact:** permanent native comparisons cover successful plain rebuilds,
failed inverse/size/degree builds, retained multiplier data, recovery, exact-degree
plain shifts, empty/trailing-zero vectors and one through four FFT-prime contexts.
Defined results, errors and supplied state match. Buffer-overrun cases are labeled
adapter guards and are excluded from claims of native behavioral equivalence.

## NTL arbitrary-modulus polynomial quotient adapters

**NTL behavior:** `ZZ_pXModulus`, `ZZ_pXMultiplier`, `build`, `rem`, `MulMod` and
`SqrMod` use the ambient `ZZ_pContext` and NTL coefficient/FFT buffers. Quotient
construction starts FFT work at degree 22, installs modulus-specific prime data,
and builds the reversed inverse with the native degree-45 Newton crossover.
Failed builds retain the fields already written. Multiplier FFT buffers remember
their allocated prime count even after a plain rebuild; changing that count raises
`FFTRep: inconsistent use` when B2 is resized, before B1 is overwritten.

**Port behavior:** Array overloads in `ntl-ts/src/ZZ_pX.ts` accept explicit quotient
objects carrying `p` and optional `PolynomialProductState`. Coefficients are copied
and normalized modulo `p`. `null` creates a default, uninitialized quotient;
`val()` returns an independent array. The legacy object-valued `ZZ_pX` overloads
remain separate. Private caches preserve the native public-field update order,
retained representations, inverse errors, mismatch errors and block-remainder
boundaries. Failed constant builds preserve the negative-buffer error before the
remainder loop. Shared state uses the existing FFTPrimeContext/RandomStream;
otherwise a private context starts with a 32-byte zero key.

**Rationale:** Explicit contexts make deterministic state observable without
process-wide globals. Array ownership maps native value copies into JavaScript.
The implementation retains NTL's reciprocal/Newton/block algorithms and the
existing portable exact product boundary. Products with stale truncated inverse
representations execute the port's actual truncated forward and inverse transforms,
where a full coefficient convolution would change native results.

**Trade-offs:** Native allocation, threading, word preconditioners and ordinary FFT
product buffers are replaced by existing BigInt product kernels. Consequently this
is not a claim of native multiplication performance or allocation-failure parity.
The stateful `ZZ_pX` coefficient class and downstream Hensel callers remain separate
audit work. Objects retain their own coefficient modulus; changing an object's
ambient modulus implicitly is not an exposed TypeScript operation.

**Behavioral impact:** Shared comparisons cover degree/crossover boundaries,
composite moduli, failed quotient/multiplier rebuilds, context restoration,
constructor failures, mixed and aliased operands, and full prime/root/reciprocal
state plus subsequent random bytes. The first 931 traces identified 28 changed-
prime-count multiplier discrepancies; a permanent comparative test fails on all
28 before the guard and passes after it. The constant-modulus cohort runs under a
fixed subprocess watchdog to detect nontermination without hanging the test suite.
This coverage does not establish equivalence over every degree or resource limit.

Additional context-reuse comparisons in 24.45.1 cover 88 arbitrary-modulus traces.
The first 48 expose 28 coefficient discrepancies: complete products now normalize
retained B1/B2 coefficients for the current modulus before packed multiplication,
while partial inverse FFTs read the retained representation under the FFT primes
without changing its original coefficient modulus. Forty additional native controls
cover larger coefficients and two-to-five-prime contexts. Same-count cache reuse
preserves native transform data; it does not require changing the quotient's `p`.

## NTL arbitrary-modulus CRT reconstruction

**NTL behavior:** `ZZ_p.cpp::DoInstall` chooses FFT primes until their product
exceeds `p^2 * 2^29`. `ZZ_pX.cpp::FromModularRep` uses a binary64 quotient estimate
through 800 primes, including the one-prime case. The bundled AArch64 build fuses
multiply-adds. Above 800 primes, `lip.cpp` selects an exact centered CRT tree,
with the native leaf partition determined by the prime count.

**Port behavior:** `FromModularRep(residues,F)` uses one row per modulus-specific
prime and returns one bigint per column. Signed residues are normalized before
reconstruction. The small-CRT path reproduces fused binary64 rounding with exact
dyadic intermediates; the large path uses the original product/remainder tree and
centered reconstruction. Lazy FFT initialization precedes row-count validation.
Ragged rows and wrong prime counts raise explicit `RangeError` adapter guards.

**Rationale and trade-offs:** JavaScript lacks native fused multiply-add and the
original Montgomery/preconditioned integer storage. Exact integer intermediates
preserve the selected native platform's rounded output. CRT trees use BigInt
rather than GMP buffers. Guards provide defined behavior for malformed matrix
shapes that do not satisfy the original internal function's buffer preconditions.
They are explicitly reproduced by the native test adapter.

**Behavioral impact:** Midpoint rounding agrees with bundled AArch64 NTL even when
it differs from exact centered CRT: with `p=6`, residue `441352763482308608` produces
`3`. This differs from the separate word-modulus one-prime adapter, which centers
exactly. Native comparisons cover signed/midpoint inputs, empty columns, shape
guards, and prime counts 799, 800, 801, 1023 and 1024. Other machine/compiler
floating-contraction profiles and resource exhaustion are outside this profile.


## NTL stateful Hensel lifting

**NTL behavior:** ZZXFactoring.cpp builds a minimum-degree word-polynomial factor
tree and obtains Bezout coefficients through word XGCD. TreeLift installs the
arbitrary-modulus context for each precision stage. HenselLift computes an integer
product and exact scalar division, builds both polynomial quotient moduli, then
both multiplier caches, applies factor corrections and lifts inverses. Its final
HenselLift1 stage uses raw products and skips the inverse caches/corrections.
AdditionalLifting normalizes the target before installing a word context with
NextPowerOfTwo(deg(f))+1 as maxroot. Completed operations retain their global FFT-table
and random-stream effects, including when a later operation fails.

**Port behavior:** MultiLift accepts a final optional zz_pXOptions object. Its
factor tree delegates to word mul/XGCD, its integer products to stateful ZZX1.mul,
and its corrections to ZZ_pXModulus/ZZ_pXMultiplier operations in native order.
AdditionalLifting accepts an optional PolynomialProductState as its final argument,
computes the native degree-dependent word setting and passes both into MultiLift.
The existing native context, constant-quotient, scalar-inverse and inexact-division
errors are preserved. The exponent-one shortcut still copies normalized factors.

**Rationale and trade-offs:** Optional explicit state follows the existing portable
NTL adapters. Coefficient arrays and independent returned tuples replace mutable
output references; AdditionalLifting does not destroy the caller's factor array
on failure. Native ambient context save/restore is represented by locally captured
coefficient settings, without changing process-global JavaScript state. Allocator
limits, verbosity, SIMD, native FFT buffer layouts and thread scheduling are not
reproduced. Complete products use the existing exact portable multiplication
boundary, with its documented complexity limits.

**Behavioral impact:** 362 native traces compare full output/error sequences,
all generated FFT-prime information and the next 64 random bytes, including cold
AdditionalLifting calls and prewarmed integer products. The oracle extracts the
seven original bundled NTL 11.6.0 bodies verbatim, renaming MultiLift for linkage,
and links the complete matching bundled library. Its C++ error callback follows
the existing JavaScript Error adapter. Cases include staged exponents through 128,
word/quotient/product crossover neighbors, mixed-degree trees, constants,
non-coprime and composite inputs, low word FFT caps, native word-modulus limits,
leading-coefficient normalization and the exponent-one shortcut. Results, errors
and supplied state match on these inputs. This coverage is bounded by those inputs
and the fixed 64-bit native profile. The complete integer factor driver uses
this lifting interface.

## NTL FFT-prime word coefficient contexts

**NTL behavior:** zz_p::FFTInit(index) uses that cached FFT prime as the coefficient
modulus, with PrimeCnt=0 (algorithm selector), NumPrimes=1 (transform rows), and
MaxRoot derived from the prime. Transforms use its root table directly and bypass
CRT reconstruction. Sequence minimal polynomials switch algorithms at 400.

**Port behavior:** zz_pXOptions adds `fftPrime?: number`. Supply p equal to that
prime and use the same coefficient context for all operands and retained caches,
as required by native NTL. This option ignores maxroot, which FFTInit does not
accept. Nested operations retain the selected prime and shared state; omitted
state creates a private zero-key stream/cache. Negative/out-of-range native
indices preserve native errors. Non-safe-integer indices and mismatched explicit
p raise RangeError, because those inconsistencies cannot be represented by the
native constructor.

**Rationale:** An explicit option represents native context selection without
adding mutable ambient state. Separating selector and row count preserves native
FFT thresholds, allocation checks and direct reconstruction.

**Trade-offs:** Explicit p and index require a consistency check. Arbitrary
UserFFTInit contexts and native out-of-context arithmetic are unsupported.

**Behavioral impact:** Supported native contexts reproduce coefficients, errors,
FFT tables and random-stream continuation. Adapter guards define errors only for
inputs with no corresponding native context construction.


## NTL complete integer factorization driver

**NTL behavior:** ZZXFactoring.factor removes signed content and squarefree
multiplicities, then applies SFFactor's deflation/inflation schedule. ll_SFFactor
extracts X and X-1, optionally reverses the polynomial, selects small primes,
computes native lifting bounds and performs multifactor Hensel lifting. Bounded
cardinality search precedes the adaptive van Hoeij trace/lattice loop. The loop
retains native sparse/dense decisions, signed random trace matrices, additional
lifting, exact LLL, CutAway and GotThem certification, including the power-hack
abandonment schedule. Returned factor order follows these stages.

**Port behavior:** The same functions delegate to the existing NTL stage ports.
Typed results replace output references: factor returns signed content and
[coefficient-array, numeric-multiplicity] pairs; squarefree drivers return copied
coefficient arrays. IntegerFactorizationOptions represents native tuning globals
and optional PolynomialProductState. The three trace-precision helpers accept
bigint primes as well as existing safe integer numbers.

**Rationale:** Preserve the native algorithms, dependency structure, factor order,
and state effects without ambient mutable C++ globals.

**Trade-offs:** Omitting state creates a private zero-key stream/cache; supply shared
state to reproduce native continuation. Verbose timing, native allocation failures,
and threaded/alternate FFT profiles are not modeled. Direct van Hoeij callers must
provide the original driver's certified lifting precision and local information;
under-lifted factors can make native certification stall. Returned arrays do not
model native destructive output-vector updates on exceptional exits.

**Behavioral impact:** The supported native profile matches complete factors,
multiplicities, errors, FFT-prime tables and stream continuation, including dense
lattice steps, additional lifting, failed certification and early abandonment.
The number-field constructor now uses this backend in Sage's NTL degree ranges.


## Number-field automorphism ordering and conjugate backend

**Sage:** For unembedded fields of degree greater than two, element comparison
preserves the PARI conjugate order. The integral generator uses the signed leading
coefficient of the original denominator-cleared defining equation. PARI compares
conjugate polynomials by degree and then coefficients from highest to lowest.
put_natural_embedding_first swaps the natural map with entry zero.

**Port:** NumberField.automorphisms reconstructs that polynomial order from the
local positive-scale integral model, reversing coefficient comparisons at even
indices when Sage's scale is negative. It then performs the same identity swap.
The one- and two-dimensional cases retain their identity-first shortcuts.

**Rationale:** Preserve the ordered public result while retaining the current
integral-model and automorphism object representation.

**Trade-offs:** The conjugates still come from the local certified lattice kernel,
whose delegation and resource limits remain open under Number-Field Kernel Not
Delegated to parigp-ts. Distinguished real embeddings, Sage's cached immutable
Sequence and backend random-state fidelity are not supplied by this repair.
PARI flag four's identity fallback cannot replace flag zero's nfroots fallback.

**Behavioral impact:** Ordered generator images agree with the bundled Sage
scaling/conversion rules for the registered unembedded fields, including rational
and negative defining-equation scales. The comparison oracle executes the bundled
_pari_absolute_structure method and explicitly applies its beta map because
installed Sage 10.3 still uses polredbest and a different conversion map.


## PARI conjugate-count bound and word derivative adapters

**PARI:** numberofconjugates tests unramified primes using Flx_is_squarefree,
counts distinct-degree factors, combines their orbit bounds and applies the
native stopping rules. Its signed starting-prime argument converts to an unsigned
word for forprime. Flx_deriv computes modular coefficient derivatives and
Flx_is_squarefree tests whether the derivative GCD has degree zero.

**Port:** The galconj function delegates to the word polynomial and prime-iterator
ports. The legacy pari_nf helper delegates to it, preserving its bigint result.
Word derivatives use exact BigInt multiplication in both native word-width paths.
All routines retain PARI's shared random state and leave input arrays unchanged.

**Rationale:** Preserve the dependency route, degree-count algorithm, prime
boundaries and early stopping while representing native arrays as coefficient lists.

**Trade-offs:** The public galconj bound returns a number (bounded by array degree),
while the legacy wrapper returns bigint. pinit represents a signed 64-bit value;
negative values undergo the same unsigned conversion as in C. The bound requires
a monic irreducible polynomial of positive degree, as in the number-field caller.
Word coefficients must be reduced modulo a word prime; their existing adapter
validation remains in force. Native stack lifetimes and debug timing are omitted.

**Behavioral impact:** Registered bounds, derivatives, squarefree predicates and
random states match bundled PARI. Negative starting bounds do not initiate a
signed scan from a negative integer. The bound remains an upper bound, which may
be the original degree when the unsigned prime range is already exhausted.


## PARI finite-field universal comparison

**PARI:** cmp_universal compares finite-field backend tags, raw values, defining
polynomials and characteristics in that order. Raw Flx coefficients and packed
F2x words are signed native long slots, compared low first after their lengths.
Large-prime FpX values compare polynomial length and coefficients low first.
Sage's PARI-backed extension elements delegate comparison to this function.

**Port:** gen2.c's finite-field branch is exposed as
`cmp_universal(x: PariFfelt, y: PariFfelt): number`. The polynomial factor comparator
delegates to it for the documented PARI extension-element profile. Reduced
coefficient arrays replace raw GEN buffers, and binary coefficients are packed
into 64-bit words before comparison. Scalar prime-field values use the defining
polynomial X when none is supplied. Inputs are not mutated.

**Rationale:** Preserve native factor/root ordering and signed-word comparisons
without introducing process-dependent native memory representations.

**Trade-offs:** The adapter covers finite-field GENs at one shared variable index,
not arbitrary GEN types or different variable priorities. Inputs require reduced
coefficients and their matching defining polynomial. The selected profile uses
64-bit words. The polynomial extension factorization algorithm itself remains a
local pipeline; this comparator does not establish full backend/state equivalence.

**Behavioral impact:** Roots now preserve factor order, including zero and
extension coordinates. Prime-field roots use factorization even for small fields,
so the PARI random-state effects of binary factorization are retained. The native
comparisons cover same/cross-field values, zero padding, binary word boundaries,
signed-word coefficients and large-prime representations. Different native word
sizes, Givaro element order and native memory allocation failures are not modeled.

## Modular Polynomial Roots and Hensel Lifting

**Sage and the port:** IntegerModRing's root hook delegates nonzero prime-modulus
polynomials with multiplicities to its cached finite field. Returned coefficients
therefore belong to GF(p), not Zmod(p). Distinct roots remain in Zmod(n): zero
lists every residue, linear cases use the bundled congruence recursion, prime
moduli use distinct finite-field roots, and composite moduli use factorization,
precision-doubling Newton lifts or singular digit lifts, followed by CRT in
Cartesian-product order. The port now follows those branches and exact errors.
`field()` and `factored_order()` cache their results; the latter delegates to the
existing PARI-backed integer factorization port.

**Known upstream bug, retained for compatibility:** bundled
`integer_mod_ring.py:1983-1997` computes a/g but recurses on f without dividing
its coefficients by g. For 2*x+2 over Zmod(4), the bundled method returns [0,2],
although the mathematical roots are [1,3]. Fifteen of the initial 3,990 shared
root cases expose this problem; the other successful small-modulus distinct
results agree with exhaustive mathematical root sets. The compatibility tests
preserve the bundled outputs, not a claim that those fifteen outputs are roots.
This is upstream behavior, not a TypeScript arithmetic discrepancy.

**Rationale:** the audit targets the bundled implementation's deterministic
behavior, including root parents, ordering and observable upstream bugs. General
polynomial ring overrides remain outside the exposed roots API. The ring hook
accepts only its own ring or null and ignores its algorithm option, like Sage.
The helper `_lift_residue_field_root` is a static TypeScript method; Python also
allows accessing its static method through an instance. Factorization results
retain the project's array-of-bigint-pairs adaptation.

**Trade-offs and behavioral impact:** the nonunit-linear upstream bug is not
corrected mathematically. Native process-global NTL state and cross-parent global
GF factory identity are not established by these cases. Per-ring cached field
identity, cached factorization identity, input ownership, exact errors, Newton
and singular lifting, and ordered CRT results are compared. Large composite
regressions have external ten-second watches; timeouts fail their tests.

Fresh TypeScript direct constructors are compared with native generic parents
initialized with the factory's refined category. A raw Sage generic parent lacks
factory metadata and its first prime-field category refinement raises
TypeError("'NoneType' object is not subscriptable"). The port has no Sage category
framework and does not reproduce that metadata failure. This explicit adapter
matches the existing category-free construction profile; category/proof options
remain outside it. See [Modular Integer Coercion and Factories](#modular-integer-coercion-and-factories).


## Weierstrass Isomorphism Root Dispatch

- **SageMath and port:** all `_isomorphisms` characteristic branches construct
  polynomials and request distinct roots. Characteristic two with j=0 uses a
  quartic in s. Root order is inherited directly; arguments must be generic
  elliptic curves. The former local square/cube helpers incorrectly used the
  characteristic as the field order, compared extension elements as strings,
  and failed over QQ. Those discrepancies are repaired.
- **Rationale:** use the shared polynomial algorithms and backend delegation,
  preserving the source's equations and first-isomorphism selection.
- **Trade-offs:** this does not add a common-parent coercion framework or new
  number-field root backends. Existing finite-field factorization/backend-state
  limitations still apply. Comparative extension parents explicitly use PARI
  elements and fixed moduli; alternative Sage element backends may order roots
  differently.
- **Behavioral impact:** native ordered tuples, parent identity, transformation
  validity, constructor selection/errors and rejected argument types match the
  supported QQ and finite-field comparisons. Broader field/coercion profiles
  remain open; this is not a claim of complete elliptic-curve equivalence.

## Generic Curve Isomorphism Ordering

- **SageMath and port:** first-isomorphism and boolean queries stop at the first
  generator result. Isomorphism and automorphism lists use the source's
  identity/negation-pair sorting key. Scalar ordering delegates to rational
  comparison or PARI's finite-field comparator; prime residues compare exactly.
  Curve arguments, base parents and morphism domains/codomains are checked.
- **Representation:** the generic-curve APIs retain their documented tuple
  results instead of returning Sage morphism objects. Morphism comparisons accept
  string operator names/punctuation instead of native integer opcodes; unknown
  operands return `null` for Python's `NotImplemented`. Ordering unequal QQ or
  finite-field curves raises the native TypeError, including category class names.
- **Rationale:** share the generator and comparison algorithm while preserving
  the existing TypeScript API and supporting its separate prime-field classes.
  Parent equality uses definitions because the port lacks Sage's global factory
  and category interning. Equivalent-parent coefficients are coerced into the
  receiving implementation before arithmetic and comparison.
  Prime parents also compare their chosen generators: GF(p) with distinct
  degree-one moduli represents distinct native factory parents, even though
  their underlying fields are isomorphic. Explicit x-1 agrees with the default.
- **Trade-offs:** comparing distinct extension-parent objects inspects the
  defining polynomial rather than native identity alone. Manually initialized
  native generic parents whose identities differ from factory parents are outside
  this adapter. Comparisons use default native proof/check metadata and the PARI
  extension backend; other factory metadata are not represented by the legacy
  parent-definition adapter. Optional field-extension arguments, general cross-field injections,
  unknown-field ordering and global interning remain outside this repaired profile.
- **Behavioral impact:** supported factory QQ/finite-field tuples, sort order,
  parent distinctions, error messages, invalid operands and early-exit root-call
  traces match native comparisons. Scalar-root callers have their separate
  audited profile below; general cross-field lifting remains open.


### PARI finite-field scalar square-root adapter

- **Source and port:** `FF.c:913` uses an optional output pointer. The TypeScript
  `FF_issquareall(PariFfelt)` represents only the requested-root form, returning
  a field-element record or null. Inputs are reduced coefficients over an
  actual finite field, using the existing variable-zero, 64-bit PARI profile.
  `finite_field_extension.FiniteFieldElement.sqrt` delegates to this root call
  and preserves the native default extend=false, all-root ordering and errors.
- **Rationale:** JavaScript has no C output pointers or native GEN memory.
  The existing record adapter carries the characteristic, defining polynomial
  and coefficients explicitly. It returns coefficient arrays even for constants.
- **Trade-offs:** This does not provide general FF_ispower exponents or native memory aliasing.
  The omitted-pointer predicate is exposed separately as FF_issquare. Invalid raw GEN
  layouts and nonfield quotient rings are outside this field-element contract.
- **Behavioral impact:** Valid root values, root ordering and retained PARI
  random state match the bundled dependency. Word-prime extensions deliberately
  retain the source gequal1(t_VECSMALL) behavior: even sqrt(1) samples randomness
  and may return -1. Binary and large-prime polynomial-one shortcuts differ.
  Generic curve scalar-root callers now delegate through these methods; their
  separate input-domain limits are recorded below.


### PARI finite-field norm and square-predicate adapters

- **Source and port:** FpX.c/Flx.c resultants use Euclidean basecases and native
  half-GCD thresholds with deferred leading-coefficient/sign accumulation.
  FpXQ/Flxq norms apply the source leading-modulus correction. FF.c dispatches
  the corresponding norms and quadratic-symbol predicates, with the binary
  and large-prime constant shortcuts. Extension is_square calls FF_issquare.
- **Rationale:** Existing polynomial arrays and PariFfelt records replace GEN
  layouts and temporary stack storage. Returned residues are bigint values and
  C truth values are booleans. The existing 64-bit backend profile is retained.
- **Trade-offs:** Raw malformed GENs and invalid finite-field definitions are
  outside the adapters. Word polynomial inputs must have reduced coefficients
  and positive word moduli, following the existing word-kernel input contract.
  Low-level resultants/norms assume the native valid polynomial/modulus domain.
- **Behavioral impact:** Compared resultants, norms and predicate values agree
  with bundled PARI and consume no random state. Nonmonic and reducible
  positive-degree modulus norms follow the native resultant formula; square
  predicates require an actual finite field. The Sage extension norm/trace
  methods, which have separate source delegation paths, remain under audit.


### Generic curve scalar-root callers

- **Source and port:** is_x_coord first invokes the base-field constructor on
  every input and replaces only conversion TypeError with the caller's native
  message. It then uses the scalar square predicate or default characteristic-two
  polynomial roots. lift_x checks canonical parent maps, embeds the coordinate
  in the curve's field when possible, or changes the curve to the coordinate's
  field when the canonical map goes the other way. Its scalar-input stage preserves
  missing-parent AttributeError. It uses distinct characteristic-two
  roots, or a guarded scalar sqrt(all=True), and sorts the y coordinates using
  the field's native order. Montgomery model construction preserves both
  distinct-root lists and evaluates every candidate's predicate before selecting
  the first maximum. QQ uses rational comparisons; PARI extension order delegates
  to cmp_universal. Coefficient representation follows _coeff_repr's parentheses
  and whitespace rules for extension elements.
- **Rationale:** The existing field methods and dependency ports implement Sage's
  scalar algorithms. Removing the curve-local nonresidue search repairs its
  even-degree extension failure and its invalid QQ cardinality assumption.
- **Trade-offs:** This profile covers coordinates already in the curve's base
  field, including QQ, prime implementations and explicit PARI extension parents.
  The existing positional all flag and point representation are retained.
  is_x_coord conversion is checked across QQ, prime/extension implementations,
  scalar wrappers, strings, null and unsupported objects. lift_x parent maps are
  checked for ZZ/QQ, the prime implementations and explicit finite extensions;
  coefficient conversion and unchanged-field curve identity match base_extend.
  Integral JavaScript numbers retain the port's integer interpretation; Python
  floats instead belong to RDF. extend=True now constructs the native-shaped
  quadratic number field (QQ), finite field (prime base), or polynomial quotient
  (explicit finite extension). Root ordering follows the default finite backend
  or the standard quadratic embedding, using exact rational comparisons.
  Real-coordinate scalar roots, global field/curve factories and other starting
  parent families remain under audit. Explicit named extensions are not automatically embedded merely because
  their degrees divide. Twisted
  Montgomery plane curves and their rational maps
  remain explicitly unsupported. Unknown coefficient fields need their own
  scalar hooks and native comparison adapter.
- **Behavioral impact:** Compared points, deterministic root order, error messages
  and dependency-call traces agree with the bundled caller bodies. Odd extension
  scalar state also agrees with bundled PARI. Polynomial-root random streams are
  covered by their separate dependency oracles rather than reconstructed by this
  caller oracle. The port-only finite torsion_points alias still has a separate
  enumeration/domain gap; sharing the scalar helper does not repair that method.


### Generic curve y-coordinate extensions

- **Source:** lift_x forms y^2 + b*y - f over the coordinate field, calls its
  fraction-field extension constructor, changes the curve and sorts the two
  conjugate roots. QQ uses NumberField; prime fields use GF; non-prime finite
  bases use the generic polynomial quotient constructor.
- **Port:** the curve routes those branches directly to the existing NumberField,
  FiniteFieldExtension and QuotientRing constructors. Supported bases are already
  fields, so fraction_field is the identity. The third positional extend flag
  preserves the existing all flag. NumberFieldElement additionally exposes the
  generic isZero spelling of is_zero; RationalPolynomial retains its variable.
- **Rationale:** these are the available field implementations behind the port's
  structural curve interface; no new inheritance hierarchy or root enumeration
  is needed. Finite extensions retain the port's existing PARI-backed arithmetic.
  For newly adjoined prime-base fields, curve-root comparison selects Sage's
  default integer encoding order for small/Givaro or binary/NTL fields, otherwise
  PARI order. This does not implement a Givaro arithmetic kernel.
- **Trade-offs:** the tests establish exact defining equations, roots, single/all
  results, promotion, point equations and coordinate-parent identity within a
  call. Cross-call global factory identity, arbitrary starting number fields,
  further lifting on new algebraic parents, and backend random-state equivalence
  for extension construction are not established. Structural point types remain
  broader than the individual algebraic element APIs (e.g. number-field parent()).
- **Behavioral impact:** the tested QQ and finite-field calls now return points
  where extend=True was previously ignored. Existing roots do not create a new
  extension; errors and extend=False behavior are preserved. The gaps above remain
  explicit audit work and are not a full-equivalence claim.

### Optimized finite-curve coordinate adapters

- **Source:** EllipticCurve_finite_field inherits is_x_coord and lift_x from the
  generic curve implementation, including canonical promotion and extend=True.
- **Port:** the optimized prime-field short-model class delegates these calls to
  its cached generic curve. Results over the original curve are adapted back to
  optimized points; promoted/extended results retain the generic point class.
- **Rationale:** reuse the audited scalar and polynomial dependencies while
  preserving the existing optimized point API for its supported prime field.
- **Trade-offs:** the two existing point classes have different coordinate access
  APIs (properties versus methods), so newly supported promotion/extension results
  have a union type. Full constructor/model conversion, arbitrary starting fields,
  global factory identity and the generic extension backend gaps remain open.
- **Behavioral impact:** tested inputs now match Sage's conversion, errors, roots,
  parent promotion and original-curve identity. Previously ignored extend flags
  now take effect. Returned coordinates follow the established port point APIs.


### Hyperelliptic root callers

- **Source:** hyperelliptic_generic.py:499-515 calls the scalar sqrt with
  all=True/extend=False, requests distinct polynomial roots in characteristic two,
  and sorts final coordinates using the element backend. odd_degree_model:605
  takes the first distinct root; jacobian_morphism.py:203 takes roots()[0][0].
- **Port:** supported scalar classes supply their own sqrt operation; explicit
  finite extensions use PARI universal comparison. Caller root order is preserved
  without an extra factor-order sort. The structural field helper still retains
  its earlier Tonelli-Shanks/iteration fallback for custom fields without sqrt.
- **Rationale:** shared helpers bridge separate TypeScript element classes while
  native scalar and polynomial dependencies determine the supported-field results.
- **Trade-offs:** custom-field fallback routing has not been audited; generic
  hyperelliptic coordinate common-parent promotion still needs work. Cross-parent
  and real/complex inputs are not covered by the same-parent root profile.
- **Behavioral impact:** exact sorted roots, scalar random state, dependency-call
  options, odd-degree model coefficients and valid Cantor reductions match the
  live profiles over QQ and the tested finite fields. The previously documented
  ValueError adaptation for invalid Cantor inputs remains intentional.


### PARI finite-field trace adapters

- **Source:** element_base.pyx:660-687 delegates trace to PARI FF_trace. FF.c:983
  dispatches to F2xq_trace, Flxq_trace or FpXQ_trace. Each kernel multiplies by the
  modulus derivative, reduces modulo that modulus and extracts one coefficient.
- **Port:** the same dispatch and formula use existing polynomial quotient
  multiplication/reduction kernels. F2x polynomials are packed bigints; word and
  arbitrary-prime polynomials are ascending bigint arrays rather than native GENs.
- **Rationale:** preserve native arithmetic and dispatch through the established
  TypeScript polynomial representation, then restore the Sage prime-field parent.
- **Trade-offs:** valid reduced inputs and separable positive-degree moduli are
  assumed by low-level kernels; malformed native-memory representations are not
  modeled. Native cached GEN modulus wrappers are represented by existing internal
  quotient contexts, not accepted as a new public input type.
- **Behavioral impact:** values, unchanged random state, return parents and public
  FF_trace delegation match the live tests. Direct FpXQ_trace retains native
  derivative-degree indexing, including small-characteristic calls; FF_trace uses
  the word/binary route there. Norm/charpoly routing is covered separately below.

### PARI bivariate polynomial storage

- **Source:** element_base.pyx:632 obtains norm from charpoly's signed constant
  term; element_pari_ffelt.pyx:982 calls FF_charpoly. Native FpXQ/Flxq charpoly
  calls polarit3.c bivariate resultants, choosing interpolation below the word
  characteristic bound and polynomial subresultants otherwise. Arbitrary-prime
  resultants interpolate. Product/remainder trees implement interpolation.
- **Port:** those algorithms and public caller routing are now ported. Existing
  Q storage is preserved: Q[eliminated degree][retained degree], with ascending
  bigint coefficients. Native variable tags and transposed GEN storage are absent.
- **Rationale:** retain the existing TypeScript resultant API and callers while
  matching native arithmetic and dispatch through the dependency package.
- **Trade-offs:** low-level inputs must follow native reduced-coefficient and
  positive-prime contracts. Charpoly uses a monic positive-degree modulus;
  resultants also accept nonmonic T. Arbitrary-prime resultants require nonzero
  T and Q. Native malformed-memory inputs/cached GEN wrappers are not modeled.
  Interpolation requires equal nonempty vectors (a RangeError guards bad lengths).
- **Behavioral impact:** live results, variable names, prime-field parents, norm
  dependency calls and RNG state agree over tested domains. Native degree-one
  binary charpoly reaches repeated interpolation points and raises Fl_inv; the
  port preserves that native error, including for scalar FF records. Sage-facing
  extension fields of degree >= 2 do not take that failing branch. Interpolation
  otherwise assumes distinct points. No new broad finite-field parity claim.

### PARI composed-sum algorithm

- **Source:** FpX_composedsum/Flx_composedsum use Newton sums and truncated series,
  with a p-adic lift when the characteristic is too small for factorial inverses.
- **Port:** FpX_composedsum still computes Res_y(P(y), Q(x-y)), now using the native
  bivariate interpolation/subresultant implementation.
- **Rationale:** preserve the existing exact composed-sum result while its separate
  native Newton/Laplace dependencies remain unported.
- **Trade-offs:** composed-sum algorithm and complexity are not native-equivalent.
- **Behavioral impact:** supported values agree in the existing composed-sum and
  ffinit tests. Repairing resultants and norm does not close this algorithm gap.


### Elliptic torsion caller and group-relation adapters

- **Source:** ell_generic.py:3465 uses generic.linear_relation to choose a p-torsion
  basis and repeatedly divides points. ell_point.py:1349 removes repeated factors
  for nonzero 2-torsion targets, propagates known orders, and sorts by (Z,X,Y).
  Finite point scalar actions preserve order/gcd(order,multiplier) (line 4402).
- **Port:** those caller steps now match, using the shared curve-coordinate
  comparator (PARI universal ordering for explicit finite extensions). The added
  generic linear_relation uses the bundled divisor/BSGS algorithm. IntegerMod's
  additive_order follows integer_mod.pyx:1769; point additive_order aliases order.
- **Rationale:** keep dependency calls, exact selected points/polynomials and
  cached orders aligned, instead of comparing only subgroup size or root sets.
- **Trade-offs:** generic finite-field point-order backend routing is repaired for prime
  fields p > 3; the other finite-field backends remain open. Generic multiplication
  now delegates over supported prime fields of characteristic > 3, while small
  characteristics/extensions still use the local path. Infinite-order
  cache values and broad number-field torsion computation are not covered here.
  The finite torsion_points convenience alias still has its enumeration gap.
- **Behavioral impact:** tested p-primary bases, division-point lists, reduced
  polynomials, validation and known-order state agree. Bundled Sage's coprime-order
  linear_relation shortcut returns (order(Q),0); the port preserves it. The oracle
  bridges Sage 10.3's missing homset zero() method and binds leading-coefficient
  sqrt to bundled PARI to avoid comparing different PARI root representatives.
  Scalar-body selection follows bundled field classification; installed Sage 10.3
  also selects finite-field points over composite residue rings.
  It does not rewrite caller branches or normalize distinct output polynomials.


### Finite elliptic point-list and group containers

- **Source:** ell_finite_field.py:116/192 enumerates points from a group basis,
  sorts them and caches an immutable Sequence. abelian_group (line 851) caches
  an AdditiveAbelianGroupWrapper and updates the gens cache after basis correction.
- **Port:** prime short-model curves now follow that enumeration and cache flow.
  The existing AbelianGroupStructure record remains the group representation;
  points returns a frozen array proxy that raises the native ValueError on edits.
- **Rationale:** preserve the established TypeScript record/array API and native
  ordering, cache identity, mutation protection and generator selection.
- **Trade-offs:** the record is not Sage's group parent and does not expose its
  element wrappers or category machinery. The array has JavaScript methods;
  use slice() for a mutable copy. Arbitrary-property writes also raise ValueError.
- **Behavioral impact:** exact prime-field point lists, repeat-call identity and
  attempted entry-assignment errors match live Sage comparisons. Group generator
  choices remain randomized; tests check basis/cache consistency, not exact RNG
  parity. Generic torsion_points and extension-field group backends remain open.


### Constructor model routing

- **Source:** constructor.py:449–462 validates list length and coerces every
  coefficient; ell_generic.py:145 owns a tuple of coerced coefficients and checks
  its discriminant. Native finite curves inherit that general equation and add
  their PARI point/group backend.
- **Port:** the default factory's newly declared five-coefficient overload returns
  EllipticCurveGeneric, preserving the input equation instead of silently changing
  it to a short model. Two-coefficient prime inputs retain the optimized class.
  Generic construction copies/coerces its coefficient tuple. Both factories match
  native length validation and singular-equation messages for tested inputs.
- **Rationale:** use the existing general equation implementation to correct model
  preservation and coefficient parents without duplicating its arithmetic.
- **Trade-offs:** the two existing curve/point APIs remain distinct: general curves
  use base_ring and x()/y(); optimized curves use field and x/y properties. The old
  five-coefficient runtime path returned an optimized curve on a different model;
  callers of that undeclared overload must now use the generic API. Specialized
  finite-field cardinality/group operations still require backend integration.
- **Behavioral impact:** exact equations, coefficients, discriminants, j-invariants,
  base-ring parents, validation, input-array independence and tested point lifts,
  doubling and negation agree over QQ, prime fields and explicit finite extensions.
  This does not claim general/extension finite-field backend parity or factory
  identity caching. Default two-coefficient extension dispatch also remains open.


### General-model PARI scalar multiplication

- **Source:** ell_point.py:4402 calls pari.ellmul and propagates cached order/gcd.
  elliptic.c:2289 maps a general Fp curve into [-27*c4,-54*c6] using
  [6,3*b2,3*a1,108*a3], calls FpE_mul, and maps the result back.
- **Port:** generic points over prime fields of characteristic > 3 now follow
  that delegation. The new cached pari_curve/__pari__ adapter keeps a general
  ellinit record; FpE_changepoint and FpE_changepointinv port FpE.c:190–233.
  mul/rmul coerce with ZZ, preserving wrapped integers and fractional-input errors.
- **Rationale:** preserve both the supplied model and native dependency routing,
  using the existing PARI Jacobian scalar kernel rather than duplicating arithmetic.
- **Trade-offs:** records and point unions replace native GEN/modular wrappers.
  The low-level change functions require reduced coordinates and invertible u.
  General ellmul records support Fp with p > 3 and binary/odd FFEllipticCurve
  records; the existing short record API stays available. Generic point order now delegates on the same supported prime parents,
  while number-field scalar paths,
  and extension/small-characteristic order/group backends remain open. Native PARI error/fallback behavior outside valid prime-field
  point inputs is not established by this batch.
- **Behavioral impact:** native kernel values, caller dependency arguments, cache
  identity, known-order propagation and scalar coercion errors match live checks.
  Fresh transformed-model inputs vary all five coefficients and preserve known
  points; word and 127-bit primes exercise the native coordinate kernels.


### General-model PARI point orders

- **Source:** ell_point.py:4873 caches ellcard on the curve, then calls ellorder
  with that multiple. elltors.c:745 changes coordinates before FpE_order; when
  no multiple is supplied, ellff_get_o caches the factored group exponent.
- **Port:** generic prime-field points with p > 3 now follow the default/explicit
  PARI caller, including zero-point shortcuts, point/curve caches and unknown
  algorithm errors. General ellinit records are accepted by ellcard and ellorder;
  a WeakMap retains the transformed model, cardinality and factored exponent.
- **Rationale:** reuse the existing PARI group kernels while preserving the
  original curve equation and the caller's dependency/cache behavior.
- **Trade-offs:** this repairs the general-model adapter, not every group kernel.
  The documented large-prime SEA dispatch difference remains. Hybrid now follows Sage's bounded-search/partial-factor schedule (25.3.0).
  Other base fields retain their existing incomplete order paths. Default generic
  prime cardinality is integrated; additional algorithms and group APIs remain open. Legacy
  short records still use cardinality when no point-order multiple is supplied;
  general records use the native group exponent. Native group RNG parity and
  arbitrary invalid point/multiple inputs have not been established.
- **Behavioral impact:** tested general-model cardinalities and point orders,
  caller arguments, cache reuse and errors match the bundled implementations.
  General models no longer produce incorrect orders by treating their original
  a4/a6 and coordinates as a short model. Compact regressions and fresh generators
  compare bundled PARI directly and execute bundled Sage caller bodies, with a
  thin curve/point adapter exposing PARI and cached state.


### Bounded factorization for elliptic hybrid orders

- **Source:** ell_point.py:794 expands bounded searches, retrieves curve order
  when sqrt_ub reaches 5000, and factors with limit=sqrt_ub. It switches to PARI
  only after all factor bases are irreducible or units. integer.pyx:4127 delegates
  bounded factorization to factorint.pyx:251, using trial_division and mpz_remove.
- **Port:** the generic point hybrid caller now follows that schedule, including
  retrying only ValueError and preserving point-cache behavior. Integer.factor
  accepts a limit and delegates to the new mirrored factorint module. Existing
  GMP removal and trial-division dependencies are reused. Generic prime curves
  expose default cardinality/order and cardinality_pari through their PARI model.
- **Rationale:** preserve the actual bounded searches and factor-completeness
  decisions, instead of substituting unbounded BSGS or full factorization.
- **Trade-offs:** factorizations remain arrays, with a negative unit represented
  by [-1,1]; the caller evaluates native completeness over those bases rather
  than constructing Sage's Factorization parent. Limits use IntegerLike and a
  signed 64-bit C-long range. The new generic cardinality methods cover the
  default prime-field path only: explicit algorithms, extension-degree options,
  other parent backends remain open; optimized-point algorithm options were
  integrated in 25.4.0.
- **Behavioral impact:** exact values, errors, bounded-search intervals, factor
  limits, PARI transition, and curve/point caches agree in live comparisons.
  A real curve with cardinality 999963794 retains composite cofactor 499981897
  at bound 8192, then factors completely at 32768; the next action is PARI.
  The oracle executes bundled point-order and factor-completeness bodies with
  observed dependencies. Its live integer factorization uses installed Sage;
  the bundled factorint.pyx/integer.pyx algorithms were reviewed directly.


### Optimized finite-point order dispatch

- **Source:** ell_point.py:745/4873 checks the point cache and infinity before
  dispatch, shares ellcard with the curve cache, delegates generic_small/hybrid
  to the base point implementation, and rejects unknown algorithms. Its
  additive_order alias accepts the same argument.
- **Port:** optimized short prime-field points now follow those decisions and
  expose order(options), additive_order(options), and _compute_order(algorithm).
  Their existing null-based curve cache is shared with the point-order caller.
- **Rationale:** preserve the established optimized point representation while
  using the already compared base algorithms. The optimized implementation
  invokes the generic point prototype only for generic_small/hybrid; those
  branches depend solely on the common group protocol and dynamically dispatch
  the later PARI transition back to the optimized implementation.
- **Trade-offs:** the optimized class still uses its existing short PARI record
  and arithmetic kernels. This does not establish general-model conversion, RNG
  parity for that separate class, nor close the
  documented large-prime SEA dispatch difference. Its scalar delegation/coercion
  and order-cache propagation are repaired in 25.5.0.
- **Behavioral impact:** exact orders, search intervals, factor limits, PARI
  calls, unknown-option errors, additive alias and point/curve caches agree with
  the bundled Sage caller. Cached points and infinity bypass option validation.


### Optimized finite-point scalar delegation

- **Source:** ell_point.py:4402 coerces the scalar through ZZ before invoking
  pari.ellmul, creates a fresh point and propagates order/gcd(order,k). PARI's
  ellffmul transforms a prime curve using its ellinit invariants, invokes FpE_mul
  and changes back to the original model.
- **Port:** optimized mul now follows that caller, including zero/infinity,
  wrapped/integral numeric scalars, fractional errors and known-order propagation.
  Its new pari_curve alias returns the existing cached short record. PARI ellmul
  adapts short records to native ellinit and selects the FpE or finite-field
  scalar branch. A WeakMap caches those invariants and detects changed short
  coefficients/characteristic to avoid stale results on mutable records.
- **Rationale:** share the native PARI scalar route while preserving the optimized
  class's established point and curve-record APIs.
- **Trade-offs:** short records encode the zero a1/a2/a3 coefficients implicitly;
  toPari, __pari__ and pari_curve keep their original short representation rather
  than changing consumers to the general record type. Native extension/small-prime
  cardinality/order/group backends, error/fallback behavior on invalid domains,
  and exact group RNG behavior remain open.
- **Behavioral impact:** live scalar values, dependency arguments, known-order
  state, coercion errors and cache identity agree on tested prime parents, including p=3.
  Direct bundled PARI comparisons include 127-bit primes, 80-bit multipliers and
  changing the coefficients of an already-used short record. No large fixture
  files are retained.


### PARI binary elliptic point kernels

- **Source:** F2xqE.c:43–247 implements binary coordinate changes and ordinary/
  supersingular point arithmetic; multiplication delegates to gen_pow_i.
  F2x.c:1058–1078 uses extended-GCD coefficient inversion and inverse-then-multiply
  division. Those inverse algorithms were already in the extension-field port.
- **Port:** F2xqE.ts mirrors the native point kernels, reuses gen_pow_i's exact
  window schedule and adds F2xq_invsafe/inv/div wrappers over the existing native
  coefficient implementation. Ordinary coefficients are packed a2; supersingular
  coefficients are [a3,a4,a3^-1]. Coordinates use the existing EllipticPoint union.
- **Rationale:** establish the actual binary dependencies before wiring Sage's
  small-characteristic/extension scalar calls. No local general-curve multiplication
  is substituted for this native backend.
- **Trade-offs:** packed bigint polynomials omit native GEN variable tags and
  ownership. Points/coefficient records are treated as mathematical values; native
  GC and object aliasing are not modeled. Inverse errors use the existing polynomial
  display convention. Curve kernels expect valid reduced field elements and the
  supplied supersingular inverse; they do not validate the curve or its a6.
  Binary FF model conversion and scalar caller integration are now available;
  binary cardinality/order/group backends remain open.
- **Behavioral impact:** tested kernel outputs and inverse errors agree with
  bundled PARI, including infinity, ordinary x=0 doubling, supersingular arithmetic,
  signed scalars, coordinate changes and nonunit inverses. Fresh valid-curve
  constructions exercise exponents across native word/arbitrary window thresholds.
  No claim of complete binary-curve caller support is made by these kernels alone.


### PARI binary elliptic model adapters

- **Source:** ff.c:1200/1363/1466 converts ordinary and supersingular binary
  models, initializes finite-field records (including singular records with j=0),
  and changes point coordinates before and after native scalar multiplication.
  elliptic.c:798 rejects singular models in ellinit_Fq. ell_point.py:4402 delegates
  finite-point multiplication to PARI and propagates cached order/gcd.
- **Port:** these binary paths are implemented over valid supplied fields;
  FF_ellinit accepts twelve precomputed invariants, ellinit_Fq computes them for
  `[j]`, `[a4,a6]` or five coefficients, and ellmul dispatches FFEllipticCurve to FF_ellmul. Generic
  Sage points use the cached binary model for GF2Field, prime parents modulo two
  and explicit binary extensions, reconstructing coordinates in their parent.
- **Rationale:** reuse native binary kernels and the source conversion formulas
  while keeping field elements distinguishable from packed polynomial integers.
- **Trade-offs:** separate typed records replace GEN arrays and native ownership.
  The low-level adapters assume valid same-field inputs; native mixed-field and
  malformed-input coercion/error behavior remains unaudited. Native GEN record
  vectors are not accepted by the typed ellinit_Fq coefficient-tuple adapter.
  Finite-extension cardinality, order and group dependencies remain unfinished. Generic binary point orders still
  use the previously documented generic fallback. Native memory aliasing is not
  reproduced; cached Sage model identity is preserved.
- **Behavioral impact:** supported scalar results, model/invariant values, singular
  handling, integer-vs-polynomial coercion, caller dependency arguments, parent
  identity, scalar errors and cached-order propagation match live comparisons.
  The generic pari_curve/__pari__ return type widens to a union; TypeScript
  consumers of prime-only APIs must narrow the record first.


### PARI odd-extension elliptic kernels

- **Source:** FlxqE.c:40–246 and FpE.c:1519–1686 implement coordinate changes,
  affine arithmetic and signed scalar multiplication over polynomial fields.
  FlxqE has a separate ordinary characteristic-three `[a2]` model; FpXQE takes
  a short-model polynomial a4. Both delegate multiplication to gen_pow_i.
- **Port:** mirrored FlxqE/FpE modules share the identical affine formulas, select
  the existing word/arbitrary-prime coefficient backends, retain the ternary
  branch and use the same powering schedule. No naive scalar loop is substituted.
- **Rationale:** reuse audited polynomial dependencies and identical formulas
  while preserving each native backend's operations and inverse diagnostics.
- **Trade-offs:** ascending arrays replace tagged GEN polynomials. Inputs require
  canonical reduced coordinates and a valid supplied modulus; variable tags,
  native GC and memory ownership/aliasing are not reproduced. Native pi reduction
  metadata is handled by the existing word polynomial port, not exposed here.
- **Behavioral impact:** native point outputs and inverse errors agree across
  word-size boundaries, characteristic-three models, coordinate changes, infinity
  and signed scalar window thresholds. Polynomial error display retains the
  existing y-variable convention.

### PARI odd-extension elliptic model adapters

- **Source:** ff.c:1145–1197/1363/1466 converts characteristic-three and larger
  odd models, wraps the finite-field invariants and restores scalar coordinates.
  elliptic.c:463/798 computes invariants and rejects singular models.
- **Port:** FF_ellinit, ellinit_Fq and FF_ellmul now dispatch to the native word
  and arbitrary-prime polynomial kernels; general Sage characteristic-three and
  explicit extension scalar callers use their cached PARI model. Native field
  parents and known point orders are preserved on the returned Sage points.
- **Rationale:** keep the source's ordinary/supersingular ternary distinction and
  select FlxqE versus FpXQE at the native 64-bit word boundary.
- **Trade-offs:** FFEllipticCurve now unites BinaryFFEllipticCurve and
  OddFFEllipticCurve. Consumers of model-specific properties must narrow the
  union. These adapters still assume valid same-field inputs; mixed/malformed
  field coercion and native ownership behavior are not established. Native GEN
  record-vector inputs and native extension/small-prime cardinality/order/group
  backends remain open.
- **Behavioral impact:** native model/invariant values, original-coordinate
  scalar results, singular handling, scalar coercion, caller dependency arguments,
  cache identity, parent identity and known-order propagation agree in live tests.
  The broader record union is a breaking TypeScript model-access change.


### PARI small-prime elliptic initialization

- **Source:** elliptic.c:773–810 initializes degree-one FF records for integer
  domains two and three, returning the empty vector for singular curves.
- **Port:** these domains return FFEllipticCurve or null, using ellinit_Fq and
  the same native invariants/model conversion. Accessors preserve field types;
  short-record scalar callers unwrap degree-one FF coordinates for their API.
- **Rationale:** typed records and null replace GEN vectors while retaining
  finite-field arithmetic and native scalar dispatch.
- **Trade-offs:** callers with a general bigint domain must handle the record
  union and null. Domain-free and p>3 singular initialization still throws the
  legacy EllipticCurveError; closing that older exception/empty-vector mismatch
  remains a separate API migration. Malformed domain/input boundaries and native
  cardinality/order/group backends are not established by this repair.
- **Behavioral impact:** small-prime model values, singular results, scalar
  outputs and tested caller state match live native/Sage comparisons. The public
  ellinit result type changes for domain-bearing calls; this is a major version.


### PARI generic and extension-curve order adapters

- **Source:** bb_group.c:560–578/668–713 parses the order/factorization and computes
  orders recursively. FlxqE.c:313, FpE.c:1727 and F2xqE.c:298 wrap that routine;
  ff.c:1495 converts the original point before dispatch.
- **Port:** gen_order uses callback power/identity operations and typed order
  forms. The three extension kernels and FF_ellorder delegate to it, preserving
  prime-factor recursion, supplied factorizations and native integer-bound errors.
  Existing FpE_order now shares that recursion, including native errors for
  nonpositive supplied bounds.
- **Rationale:** reuse the native schedule without introducing linear searches
  or recomputing orders outside the dependency library.
- **Trade-offs:** factor rows and callback values replace GEN matrices/objects.
  Missing-order error tags cover integer and point/vector representations; other
  abstract callback value types have no native GEN tag. Factors are assumed to
  describe a positive annihilating bound. Native behavior outside that precondition
  (including malformed/inconsistent factor matrices) is not generally established.
  Automatic FF cardinality/group-exponent computation and default ellorder/Sage
  caller integration remain open; these entry points require the supplied bound.
- **Behavioral impact:** live tests compare callback order and results for the
  generic routine, plus actual native binary, ternary, other word-prime and
  127-bit-prime FF point orders. No fallback enumeration is used in production.


### PARI base-field extension cardinality adapters

- **Source:** FpE.c:2046–2062 powers X modulo X²-tX+q to extend a Frobenius
  trace, then counts base-field models over extensions with that trace.
- **Port:** elltrace_extension uses exact coefficient pairs and the native word
  powering schedule. Fp_ffellcard delegates the base count to the existing
  prime-curve ellcard backend before extending its trace.
- **Rationale:** preserve the source's logarithmic degree dependence without
  introducing a general symbolic polynomial object for the quadratic quotient.
- **Trade-offs:** degrees use the port's native number convention for dimensions
  and must be nonnegative safe integers. The base counter still has the documented
  SEA-dispatch/word-Shanks fidelity gaps; this wrapper does not repair those.
  General extension coefficients and FF cardinality dispatch are not implemented
  by these two helpers.
- **Behavioral impact:** trace values and base-field extension counts match live
  native calls, including degree zero, negative traces and large integer inputs.
  Existing base-counter algorithm differences affect performance and native
  operation scheduling, and remain open work before complete cardinality fidelity.
