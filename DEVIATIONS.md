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

56. [Template for New Deviations](#template-for-new-deviations)

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
| `QQ(...)` (`RationalField.__call__`) | `rational.pyx:591-704` accepts `Rational`, `Integer`, `int`, `float`, `str`, `[n, d]`, a length-1 list, `fractions.Fraction`, anything with `_rational_()`, and `(n, d)` | A port of that method with Sage's exact `TypeError("unable to convert {!r} to a rational")` and `ValueError('denominator must not be 0')`. Two additions: the two-argument spelling is kept (Sage passes `d` as a *base*, and string bases are unsupported here), and `{numer, denom}` objects are accepted alongside Sage's `{numerator, denominator}` |
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
| `sampleExact()`, `samplesExact()`, `isIntegral`, `basisExact`, `cExact` | `stats/distributions/discrete_gaussian_lattice.ts` | See [Discrete Gaussian Samplers](#discrete-gaussian-samplers) |
| `Integer.nth_root_mod` | `rings/integer_ring.ts` | Sage exposes this only as `Mod(a, p).nth_root(n)`; there is no `Integer` method, so there is no doctest to match. When `a == 1` and `gcd(n, q-1) > 1` we return 1 where Sage's `_nth_root_common` returns a *primitive* gcd-th root of unity (`Mod(1,11).nth_root(5)` is **9** in SageMath 10.3, deterministically across runs). Every returned value satisfies `r^n == a (mod p)` |
| `sqrt_mod(a, m)` | `arith/misc.ts` | `sage/arith/misc.py:2274` is a **commented-out stub** (`# def sqrt_mod(a, m):`), so there is no upstream function. Returns `null` for a non-residue |
| `Ti`, `GF2n(n)` | `rings/finite_rings/tower_field.ts` | Binary tower fields as used in Binius. There is no `sage/rings/finite_rings/tower_field.py`, so this module has no SageMath counterpart. `GF2n` accepts only powers of two — a property of the Binius construction, not a restriction on `GF(2^n)`, which `FiniteField(2^n)` builds for any `n` |
| `RandState.random()`, `python_random(seed)` | `misc/randstate.ts` | `random()` is a `@deprecated` alias for `c_rand_double()` retained for three call sites; `python_random`'s `seed` parameter lets a caller get a self-contained stream instead of inheriting the global state (see [Random State](#random-state-and-seeding)) |
| `PythonRandom` | `misc/randstate.ts` | A port of CPython's `random.Random`, which Sage reaches through `randstate.python_random()` rather than exposing as a class |
| `NumberFieldElement.is_integral_unit()` | `rings/number_field/` | Sage overloads `is_unit()` on the parent's type; without a distinct element class for orders the two behaviours cannot share one name |
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
| `Rational.valuation()`, `val_unit()`, `RationalField.quadratic_defect()` | `infinity` | `'Infinity'` |
| `pAdicGenericElement.valuation()`, `multiplicative_order()`, `additive_order()` | `infinity` | `Number.POSITIVE_INFINITY`, typed `InfiniteOr<bigint> = bigint \| number` |
| `RealNumber`/`ComplexNumber` `multiplicative_order()`, `additive_order()` | `infinity` | `Number.POSITIVE_INFINITY` |
| `FreeModule.cardinality()`, `indexIn()` | `infinity` | `Number.POSITIVE_INFINITY` |
| `FunctionFieldElement.valuation()`, `FunctionFieldIdeal.valuation()` | `+Infinity` | `Number.POSITIVE_INFINITY` (typed `bigint \| number`) |
| `QuotientRing.cardinality()` | `Infinity` | `Number.POSITIVE_INFINITY` |
| `arith.valuation(0, p)` | `+Infinity` | **Raises** `ValueError('valuation of 0 is infinite')` |

### Rationale

1. **`bigint` has no infinite value**, and widening every arithmetic return type to admit one would
   ripple through every caller.
2. **Numeric comparisons work** — JavaScript relational comparisons between `bigint` and `number`
   behave as expected (`Infinity > 0n`), so guards read naturally.
3. **Module-local consistency** — the string sentinel predates this register in `rings/rational*.ts`;
   changing it there is a larger breaking change than the inconsistency costs.
4. **`arith.valuation` returns `bigint`**, so introducing a sentinel would change the type for
   every caller; it throws instead.

### Trade-offs

- Three sentinels for one concept; callers must know which module they are in.
- `valuation(0, p)` fails where SageMath answers.

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
| `Integer.real_log` | MPFR | Exact above 2^53 (`ln(10^k) = k·ln10`), `Math.log` below |
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
| PARI 2.18.1's Schoenhage fast reduction for binary quadratic forms | `qfi_red` (`Qfb.c:975-987`), `qfi_redsl2` (`:857-880`), `qfr_red_i` (`:914-940`) and `qfr_redsl2` (`:825-855`) negate `b` before `pqfbred_rec` (which requires non-negative coefficients) and never negate it back. For `b < 0` that returns the reduced form of the **conjugate** class, and `qfi_redsl2` negates only the second *row* of `U`, giving `det U = -1` | The intermediate form is conjugated back and `D·U·D` (`D = diag(1,-1)`) is used instead of a row negation, so `det U == 1` and `Q ∘ U == result`. `qfr_redsl2` additionally falls back to the base case when `a < 0` | The path is gated at 9000 bits of excess size, so it is unreachable for anything Sage currently does — but "never return a plausible wrong answer" outranks bug-compatibility. Verified by forcing the threshold to −∞ on 80 sample forms |
| PARI's `matdetmod` documentation | `reference/pari/src/functions/linear_algebra/matdetmod` claims `matdetmod([4,2,3;4,5,6;7,8,9],27) == 9` | Returns `18` | The determinant is −9, i.e. 18 mod 27, and PARI's **own regression output** (`test/32/bbhnf` line 240) records 18. A stale doc block. Sage's `qfbsolve` doc block is stale in the same way — it predates `allsols`' `v >= 0` normalisation and the lexsort at `Qfb.c:1930` — so we follow the source and assert solution *sets* |
| `arith.CRT_basis(moduli, false)` | `arith/misc.py:3694-3726` appends to `cs` inside the coprime loop and, after catching the `ValueError`, does `cs.extend(...)` onto that same non-empty list. Running the vendored body verbatim under Sage for `[7,6,10]` gives `cs = [120, 120, -140, 21]` — **four** entries for three moduli | The partial entries are discarded, so the result always has exactly `moduli.length` coefficients (`[120n, -140n, 21n]`) | A basis longer than the modulus list makes `CRT_vectors` index out of range or silently mis-combine. The documented `[60,90,150]` case is unaffected (Sage bails on the first modulus, so `cs` is empty) |
| `groups.generic.discrete_log_lambda` with `N = 1` | `k = 0`; the loop body never runs and `hash % k` raises `ZeroDivisionError` | `k` is forced to at least 1 (`groups/generic.ts:1473`) | Avoids a crash where Sage crashes; results for all `N >= 2` are unaffected |
| `ell_finite_field.twists()` `break` placement | `ell_finite_field.py:1940-1944` puts the `break` at the for-loop level, so only `twists[0]` is ever tested for isomorphism with `self` | **Replicated verbatim**, with a comment citing the line numbers | CLAUDE.md requires behavioural equivalence with the vendored Sage. Implementing the docstring's stated intent would change the returned ordering for `j = 0`/`1728` curves. The returned *set* is complete and pairwise non-isomorphic either way |
| `FunctionFieldIdeal.is_zero()` | Inherits from `Element`, whose `__bool__` falls back to `True` when the parent has no zero — and the parent is the *multiplicative* `IdealMonoid`. So `O.ideal(K(0)).is_zero()` is `False` and the zero ideal prints as `Ideal (0) of …`, never `Zero ideal of …` | Returns `false` unconditionally, with the upstream reason cited in the docstring | This **is** upstream's observable behaviour (executed on 10.3) and is visible through `_repr_` and through `divisor()` never taking its "not defined for zero ideal" branch. Recorded so nobody "fixes" it later |
| Factoring the zero element in a function field | Backend-dependent: over GF(2) (NTL `GF2X`) `K(0).factor()`, `O.ideal(0).factor()` and `.divisor()` all raise `ArithmeticError: factorization of 0 is not defined`; over GF(3), GF(5), GF(7), GF(11), GF(13) (the `FpT` backend) they return the empty factorization with unit 0 and `divisor()` returns 0 | Always the odd-characteristic behaviour: `{unit: 0, factors: []}`, zero divisor `0` (finite order) / `Place (1/x)` (infinite order, which is what Sage's `den.degree() - num.degree() = 0 - (-1) = 1` gives) | Upstream is internally inconsistent — the difference comes purely from the polynomial representation — and the majority behaviour is the non-throwing one. Over a GF(2) constant field we return an empty factorization where SageMath raises |
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

- Code written against SageMath's buggy `CRT_basis` output length would behave differently here.
- Divergence from the vendored 2.18-dev PARI source in `gen_ellgroup`.
- Code catching Python's incidental `IndexError`/`AttributeError` will not catch our named errors.

### Behavioral Impact

`CRT_basis` returns `moduli.length` coefficients always (a regression test asserts the invariant
over four non-coprime modulus lists). `discrete_log_lambda` answers where Sage raises.
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
| Zassenhaus factorisation over Z in `pari_nf.ts` | Hensel lift | One big prime, throwing above 2^200 | `NotImplementedError` rather than a wrong factorisation |
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

> **One consumer is off by one draw.** `ZZ.random_element` is **not** currently bit-exact:
> `integer_ring.pyx:800-801` draws `den = rstate.c_random() - SAGE_RAND_MAX/2` **unconditionally**,
> before branching on the distribution, whereas `rings/integer_ring.ts:220-243` draws it only inside
> the `'1/n'` branch. Every seeded `ZZ.random_element(x)`, `ZZ.random_element(x, y)` and
> `random_element(distribution='uniform')` stream is therefore shifted by one word per call:
> `set_random_seed(0)` then four `ZZ.random_element(10^30)` gives Sage
> `670431516147804558529383265611, 772308321268490156498894882619,
> 551349305655019862415052218319, 369074466760966749087383998069` and the port
> `278091837517481385900793178228, 369074466760966749087383998069, …`. Inserting one `c_random()`
> before each `random_below` reproduces Sage's four values exactly. **This is a one-line fix, not an
> accepted deviation**; it is stated here rather than in Part II because it is a defect in an
> otherwise-verified entry.

---

## Vendored SageMath 10.9.beta4 vs Installed 10.3

CLAUDE.md directs the port at the vendored tree under `reference/sage`, which is 10.9.beta4. The
locally installed oracle is SageMath 10.3. Where the two disagree, the port follows the **vendored**
source, so a reader comparing against a stock 10.3 will see differences that are version drift, not
port defects.

| Site | Installed 10.3 | Vendored 10.9.beta4 / port |
|------|----------------|----------------------------|
| `Rational.round()` default mode | `'away'`, with a `DeprecationWarning` | `rational.pyx:3403` is `def round(Rational self, mode="even")`. The port defaults to `'even'` (`rational.ts:587`). All six modes (`toward`, `away`, `up`, `down`, `even`, `odd`) exist in both and agree value for value |
| `BinaryQF.is_reduced(D=0)` | Returns `False` | Raises `ValueError('the quadratic form must be non-singular')` (Sage issue #37635 rewrote `is_reduced` after 10.3). The port implements the 10.9 behaviour; the property-test cases route `D = 0` through other methods and short-circuit `canonical` to the literal string `'singular'` |
| `discrete_gaussian_lattice.py` | Brute-force `_normalisation_factor_zz(tau=3)`; `c` a property; repr `Discrete Gaussian sampler with σ = %f`; `_c_in_lattice` ignores whether `_G == 1` | The vendored module, loaded by path in the property-test oracle. Measured differences: `nf(ZZ^8, 0.5)` is `6.81052960784091` in 10.3 vs `6.82492448921763` in 10.9; `nf(ZZ^3, 1.0)` `15.5284660320764` vs `15.7496101985309`; sampling on basis `[[1,3,0],[-2,5,1],[3,-4,2]]` takes `_call_in_lattice` in 10.3 and `_call` in 10.9. Comparing against 10.3 would have reported four false port defects |
| `groups.generic.order_from_multiple` | A bare `assert` (message-less `AssertionError`) | `ValueError(f"The order of P(={P}) does not divide {M}")` (`generic.py:1361`). Also `order_from_bounds(P, None, …)` and `IntegerMod.log(b, order, check=)` (`integer_mod.pyx:795-798`) exist only in 10.9 |
| `basis_for_quaternion_lattice` | `reverse=False` (deprecated); `maximal_order` fails with `ValueError('basis must have rank 4')` for invariants such as `(-4,-28)`, `(-292,-732)`, `(-48,-564)`, `(-436,-768)` (Sage issue 37417); no `is_definite`/`ramified_places`/`order_with_level`/`pushforward`/`pullback`/`reduced_basis`/`is_principal(certificate)` | `reverse=True` default, issue 37417 fixed, all newer methods present. Printed ideal/order bases can differ from 10.3's; the **lattices** are identical, verified by comparing echelon basis matrices |
| `discrete_gaussian_integer` `precision='qq'` message | Adds a trailing period | `discrete_gaussian_integer.pyx:400` has none; the port matches the vendored source |
| `Frac(GF(p)[x])` element comparison for small odd `p` | `FpTElement._richcmp_` compares `(numerator, denominator)` lexicographically; its own docstring says "the ordering is arbitrary" (`fraction_field_FpT.pyx:376`) | The generic `richcmp(a.num*b.den, a.den*b.num)` (`fraction_field_element.pyx:994`), which is what GF(2), GF(4), GF(65537), GF(131101) and QQ use upstream. See [Function Fields](#function-fields) |
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

Covers `roots()`/`factor()` over the supported rings and the residual differences between our
integer-polynomial factorization and FLINT's `fmpz_poly_factor`.

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `roots()` ring support | Many exact/approximate rings (ZZ, QQ, finite fields, RR/CC, p-adics, number fields, symbolic) | Base ring finite fields, ZZ, QQ only |
| `roots()` options | `ring`, `multiplicities`, `algorithm` | No ring override; always returns multiplicities |
| `factor()` ring support | Broad, via FLINT/NTL/Singular/PARI | Finite fields, ZZ, QQ only |
| `factor()` return type | A `Factorization` holding the unit separately (`F.unit()`); iterating yields only the non-unit factors | `Array<[Polynomial, number]>`; when the unit is not 1 it is included as an extra degree-0 factor, sorted first (`6x^2+x-2` gives `(6)^1 (x + -1/2)^1 (x + 2/3)^1`) |
| Integer/rational factorization pipeline | FLINT `fmpz_poly_factor`: Zassenhaus, with van Hoeij/LLL recombination for hard inputs | The **same pipeline**, transcribed routine by routine from `factor_zassenhaus.c`, `factor_van_hoeij.c`, `CLD_mat.c`, `CLD_bound.c`, `van_hoeij_check_if_solved.c`, `next_col_van_hoeij.c`, `col_partition.c`, `zassenhaus_subset.c`, with the same dispatch (`r == 1` irreducible, `r <= 8` Zassenhaus recombination, `r > 8` van Hoeij) |
| Mignotte coefficient bound | `_fmpz_poly_factor_mignotte` | Transcribed verbatim **including** upstream's `fmpz_set_ui(b, m-1)` initialisation, which makes the returned bound `m−1` times the textbook Mignotte bound — still valid, so it is reproduced rather than "fixed" |
| LLL inside van Hoeij | `fmpz_lll_wrapper_with_removal_knapsack`: doubles, then a heuristic, then MPF, each verified against an exact predicate | Cohen Algorithm 2.6.7, **exact integral** Gram-Schmidt with `delta = 0.99`, so FLINT's own `is_reduced_with_removal` predicate holds by construction |
| Hensel lifting | Binary product tree with quadratic steps (`hensel_build_tree.c`) | The linear multifactor lift (von zur Gathen & Gerhard, Algorithm 15.17). Output is identical — monic `H_i` with `prod H_i = f/lc(f) mod p^a` in the symmetric range; only the cost differs. When van Hoeij doubles the precision we re-lift from `p` rather than continuing the tree |
| `fmpz_mat_col_partition` | Numbers its classes in `qsort` order of a hash | Numbers them by first occurrence, comparing columns exactly. Same partition, different labels; the trial factors are sorted by degree immediately afterwards |
| `fmpz_poly_CLD_bound` | Explicitly inexact (doubles) | **Kept in doubles.** `N·max(B_1, B_2)` is a valid bound for every `r > 0`; a looser bound costs an extra Hensel doubling and cannot produce a wrong factorization |
| Degree 2 / 3 fast paths | `_fmpz_poly_factor_quadratic` / `_cubic` closed forms | Not ported; the general route gives the same factors |
| Final check | — | The product of the returned factors is verified against the input; a mismatch raises `ArithmeticError` rather than returning a wrong factorization |
| Affected modules | `sage/rings/polynomial/polynomial_element.pyx`, `flint/fmpz_poly_factor/*` | `packages/sagemath-ts/src/rings/polynomial/polynomial_element.ts` |

### Rationale

1. **No `Factorization` type** — returning the unit as a degree-0 factor restores
   `prod(factors) === f`, which several call sites and tests rely on.
2. **Exact LLL rather than FLINT's floating-point chain.** FLINT's own exact predicate
   (`gr_mat_is_row_lll_reduced_with_removal_naive`) is the specification its three approximate
   implementations are measured against; satisfying it directly is simpler and inside CLAUDE.md's
   no-floating-point rule.
3. **Reproduce upstream's arithmetic where the value is observable** (the Mignotte initialisation,
   the CLD bound), and be exact where it is an implementation artefact (the LLL).

### Trade-offs

- Roots and factorization over rings outside the supported set are unavailable, and algorithm
  selection / ring overrides are not exposed.
- Callers must skip the degree-0 entry to iterate "the factors" the way Sage does.
- The reduced basis our LLL returns can differ from FLINT's on inputs with several near-optimal
  reductions. van Hoeij only needs *some* reduced basis, and every candidate factor is certified by
  exact trial division, so this cannot affect the answer.

### Behavioral Impact

Factorization over ZZ and QQ is **correct and complete**, verified factor-for-factor against
SageMath 10.3's own `factor()` (which *is* FLINT's `fmpz_poly_factor`) on 729 polynomials with
0 mismatches, and against an independently written Kronecker-method oracle on 200 random QQ
polynomials of degree <= 6, 400 composite products and 150 products of oracle-certified
irreducibles up to degree 14. `x^105 − 1` factors into its eight cyclotomic pieces in 109 ms;
Swinnerton-Dyer polynomials of degree 16 and 32 are proved irreducible in 12 ms and 23 ms;
`x^2 − primorial(10007)` in 105 ms. The LLL is checked on 300 random lattices against a from-scratch
port of FLINT's exact `is_reduced_with_removal` predicate, and the CLD bound contract
`|[x^n] f g'/g| <= CLD_bound(f, n)` on 4065 checks.

For the printed form of a polynomial and the handling of the integer `factor()` unit, see
[Polynomials](#polynomials--printing-factor-shape-term-orders-and-base-rings).

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
`modulus=` keyword is **not** accepted by either — see
[Finite Fields](#finite-fields--conway-table-and-minimal-polynomials).

---

## Generic Group API and DLP

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `discrete_log` options | Supports `bounds`, `algorithm` (`bsgs`, `rho`, `lambda`), `verify`, and `ord=oo` | Implements exactly Sage's `bounds=None, algorithm='bsgs', verify=True` path — including the repair of an `ord` that is a proper multiple of the base's order, and the `<30` linear branch of `bsgs`. `ord` is optional (`generic.ts:703`) and falls back to the base element's order |
| `discrete_log_rho` options | `ord` optional; configurable `hash_function` (`generic.py:659`) | `ord` **required and must be prime** (`generic.ts:1573-1592`) — stricter than Sage; no custom hash function |
| Hashing in DLP algorithms | Element hashing/equality | `bsgs` and `discrete_log_rho` hash via string representations |
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
- `discrete_log_rho` refuses composite orders that Sage accepts.
- Potential hash collisions for group elements with non-unique `toString()` outputs.
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
| `change_ring(matrix, ring)` | `matrix0.pyx:1666-1715`, relying on the coercion framework | No coercion framework, so `_coerce_entry` asks the target ring to convert and, on failure, builds the canonical morphism only in the two cases where one **provably** exists: `QQ -> R` (`n/d ↦ R(n)·R(d)^-1`, raising if `R(d)` is not a unit) and `Z/mZ -> R` when `char(R) \| m` (or `m = 0`). Anything else raises `TypeError`, so `Z/8 -> GF(7)` is correctly refused |
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
| Weierstrass models | General form in all characteristics | `ell_generic.ts` stores all five a-invariants and `ell_point.ts:235-331` implements the full general-Weierstrass chord/tangent formulas with `a1,a2,a3`, in **every** characteristic: `EllipticCurve(GF(2),[1,0,0,0,1])` builds, gives `disc = 1`, `j = 1`, and `P = (0:1:1)` doubles to `(0:1:0)`, `3P = (0:1:1)`, `4P = (0:1:0)` — identical to Sage. The **other** curve class, `EllipticCurveFiniteField` (`ell_finite_field.ts`, what `index.ts:77` exports as the default `EllipticCurve`), stores only `a`, `b` and converts via `c4`/`c6`, so it raises `ValueError('General Weierstrass form in characteristic 2 not yet supported')` (and the same for 3) at `:1054-1062` |
| `is_j_supersingular` | Checks `supersingular_j_polynomial(p)(j) == 0` when `p` is in the precomputed table, giving an exact answer even with `proof=False` | Skips the table (`supersingular_j_polynomial` is not ported) and always falls through to the 10 random-point tests (`ell_finite_field.ts:1463`), plus the trace-of-Frobenius check when `proof` is set (the default). With `proof=True` — Sage's and our default — the answer is identical and proved |
| `division_points` for 2-torsion `P` | `ell_point.py:1531-1557` replaces `g` by `gcd(g, g')·sqrt(lc(g))` (times `(x − x(P))` for odd `m`) | Uses `g` unreduced and iterates its distinct roots. The reduction only strips repeated factors, so the *set* of roots — the only thing the algorithm consumes — is unchanged; verified equivalent on 10 512 brute-force cases |
| `montgomery_model` representative | `EllipticCurveIsogeny(GF(7) j=1728 curve, (0,0), model='montgomery')` reports `A = 1` | Returns `A = 6`, the other root of the defining cubic. Both are valid Montgomery forms; see the root-ordering row under [Polynomial Roots](#polynomial-roots-and-factorization) |
| `possible_isogeny_degrees(E)` over Q | Billerey/Larson bounds | Mazur's list `[2,3,5,7,11,13,17,19,37,43,67,163]`, optionally intersected with the degrees for which `isogenies_prime_degree` finds an isogeny. Correct as a **superset** over Q; **not valid over larger number fields** |
| `isogeny_degrees_cm(E)` | Exact | Ported including the horizontal-primes step (`isogeny_class.py:1309-1317`) and the `n/(2h)` downward-ramified test. The function's contract ("this list is not necessarily minimal") holds |
| `Frobenius_filter` good-reduction test | `E.has_good_reduction(p)` on the **minimal** model (Tate's algorithm via `local_data`) | `v_P(disc) == 0` of the *global integral* model built in place (Laska-Kraus-Connell minimisation is not ported). A non-minimal model makes a few extra primes look bad; those are skipped, which can only make the filter **weaker** (a superset), never unsound |
| `isogenies(fill=true)` | `isogeny_class.py:369-370` raises `NotImplementedError` | The same — this is **not** a deviation, it is upstream's behaviour, recorded here so it is not re-reported |
| `qf_matrix()` | `ValueError('qf_matrix only defined for isogeny classes with rational CM')` (`isogeny_class.py:329-330`) | The same `ValueError` (`isogeny_class.ts:340-345`) |
| Affected modules | `sage/schemes/elliptic_curves/` | `packages/sagemath-ts/src/schemes/elliptic_curves/` |

### Rationale

1. **Unported dependencies** — `supersingular_j_polynomial`, Laska-Kraus-Connell minimisation.
2. **Same answer, cheaper route** — the `division_points` root set is unchanged; only the algorithm
   differs.
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
| `minimal_element()` | PARI `qfminim(q, NULL, NULL, 1)`, i.e. the column of PARI's `qflllgram` transform with the smallest reduced-Gram diagonal entry | The same algorithm (`bibli1.c:1355-1365`) but on `parigp-ts`'s `lllgramint`, which returns a different valid LLL transform on some inputs. In 4 of 254 sweep cases the returned element differs from SageMath's; the reduced norm is always the true minimum and the element always lies in the ideal (both checked against Sage) |
| `P1List` for prime level | `cyclic_right_subideals` uses `sage.modular.modsym.p1list.P1List(p)` | A private `p1list(p)` / `p1_normalize(p, u, v)` valid for prime `p` (the only case `cyclic_right_subideals` supports). Reproduces `P1List(3).list()` and `P1List(5).list()` and, end to end, the exact order of the returned subideals for `p = 3, 5, 7, 13` |
| `intersection_of_row_modules_over_ZZ` | `s.right_kernel_matrix(algorithm='pari', basis='computed')`, whose particular ZZ-basis appears in the doctest matrix | The repo's Smith-form kernel, so the returned 4×4 matrix can differ by a unimodular factor. The **row module** — the only thing every caller uses, via `ZZLattice.span` — is identical |
| Pickling, Magma, hashing | `unpickle_QuaternionAlgebra_v0`, `_magma_init_`, `__hash__` | Omitted; algebra identity comes from the factory cache and equality from `eq()` on canonical lattices |
| Affected modules | `sage/algebras/quatalg/` | `packages/sagemath-ts/src/algebras/quatalg/` |

### Rationale

1. **The dependencies are not ported** — `sage.quadratic_forms.QuadraticForm`,
   `sage.modules.free_module` in a form usable for ZZ-spans of rational vectors, ring morphisms and
   `sage.modular`.
2. **The lattice is what every consumer uses**, and it is identical; only the presentation differs.
3. The `minimal_element` divergence has its correct fix in `parigp-ts`'s `lllgramint`, not here.

### Trade-offs

- No `Factorization`, `QuadraticForm` or `FreeModule` objects to hand back.
- `minimal_element` affects the *certificates* of `is_principal` / `is_right_equivalent` /
  `is_left_equivalent`, the conjugator from `isomorphism_to` (sign/choice), and the particular basis
  returned by `reduced_basis`.
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
| Element and ideal comparison over small odd prime constant fields | `Frac(GF(p)[x])` for small odd `p` uses `FpTElement._richcmp_`, which compares `(numerator, denominator)` lexicographically and whose own docstring says "the ordering is arbitrary" (`fraction_field_FpT.pyx:376`). GF(2), GF(4), GF(65537), GF(131101) and QQ use the generic `richcmp(a.num·b.den, a.den·b.num)` (`fraction_field_element.pyx:994`) | Always the generic cross-multiplication rule |
| `element.factor()` / `ideal.factor()` | A `Factorization` with a `.unit()`; `Factorization(self._factor(), cr=True)` | `{ unit, factors }` and a sorted array of `[ideal, exponent]`. Both reproduce `Factorization.sort`'s ordering: `(degree, exponent, prime)` for elements and the prime itself for ideals (`factorization.py:671`) |
| `place.residue_field()`, `valuation_ring.residue_field()`, `divisor.function_space()` | `(k, mor_from_k, mor_to_k)` as `FunctionFieldRingMorphism` objects; `(V, mor_from_V, mor_to_V)` with a `VectorSpace` | `[k, from_k, to_k]` and `[dimension, from_V, to_V]` with plain functions. `to_V(from_V(v)) == v` is tested |
| `RationalFunctionField.field()`, `element()` | `Frac(k[x])`; the underlying `FractionFieldElement` | `k[x]` (the ring whose fraction field it is) and `f` itself. `FunctionFieldElement_rational` stores the normalised numerator/denominator pair directly, and `numerator()`/`denominator()` give the same data |
| `FunctionFieldMaximalOrderInfinite` | Multiple inheritance from `FunctionFieldMaximalOrder` and `FunctionFieldOrderInfinite` | Extends `FunctionFieldOrderInfinite` only and overrides `_repr_`. That is the branch `FunctionFieldPlace._richcmp_` depends on (`isinstance(prime.ring(), FunctionFieldOrderInfinite)`, `place.py:166`), so place ordering is preserved |
| `_place_class` | A class attribute set in `RationalFunctionField.__init__`, consumed by `PlaceSet` as `self.Element` | An abstract factory method `_place_class(parent, prime)`. TypeScript has no `element_class` mechanism; the constructed objects are identical |
| `ConstantField` / `ConstantFieldElement` | Any object in `Fields()`; the coercion framework does the rest | A new `constant_field.ts` (no upstream counterpart) declaring the structural interface a constant field must satisfy, plus `constant_field_characteristic/_is_finite/_cardinality/_element_list`, `compare_constants` and `divide_constants`, tolerating `characteristic` as a property or a method |
| `IndexError` | `K.gen(1)` raises `IndexError("Only one generator.")`; `O.gen(1)` raises `IndexError("there is only one generator")` | JavaScript's `RangeError` with SageMath's exact message (`errors.ts` has no `IndexError`) |
| Valuations | `+Infinity` | `bigint \| number` with `Number.POSITIVE_INFINITY`; see [Infinity Representation](#infinity-representation) |
| Affected modules | `sage/rings/function_field/` | `packages/sagemath-ts/src/rings/function_field/` |

### Rationale

1. **One representation means one comparison rule**, and the generic cross-multiplication rule is
   the mathematically meaningful one, is the default in `sage.rings.fraction_field_element`, and is
   what SageMath itself uses for every constant field except small odd primes.
2. **No morphism category, no `VectorSpace`, no `Factorization`, no coercion framework** — the maps
   and the ordered content are the usable payload.
3. **No multiple inheritance in TypeScript**; the branch that carries place ordering was kept.

### Trade-offs

- `a.cmp(b)` and `I.cmp(J)` can differ in **sign** from SageMath over `GF(p)` for small odd `p`. It
  does **not** affect place ordering (places are split infinite-first, and finite places compare
  monic polynomials where both rules agree), divisor ordering, supports, or any Riemann-Roch output —
  all byte-identical. In a 34.8k-line transcript run the only divergences were 164 line pairs, all
  `cmp`/`Ilt`, all over GF(3)/GF(5)/GF(7)/GF(11).
- No `Factorization` repr and no `.prod()`; no morphism domain/codomain introspection.
- `field()`'s repr differs; `element()` is the identity.
- `compare_constants` falls back to string comparison for constant fields whose elements expose no
  integer lift, which would order such constants differently from SageMath. Prime fields (the tested
  cases) all expose `toBigInt()` and match.
- The `_place_class` factory and `constant_field.ts` are architectural and arguably belong in
  `DESIGN.md`; they are listed here because `constant_field.ts` has no mirrored upstream path.

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
| `DiscreteGaussianDistributionLatticeSampler.__call__` | Returns a vector over the base ring of the basis matrix — integers for `ZZ^n`, rationals for a `QQ` basis | `sample()`/`samples()` keep the `bigint[]` signature and throw `ValueError('lattice basis is not integral; use sampleExact() for exact rational samples')`; `sampleExact()`/`samplesExact()` return exact `Rational[]` for any basis. Public `isIntegral`, `basisExact`, `cExact` were added |
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
- **Non-positive-definite `Σ` raises the built-in `Error`**, not `RuntimeError`, with Sage's message
  text byte-identical including its column padding.
- **`sigma` and `c` are methods** (`D.sigma()`, `D.c()`), matching the vendored source, rather than
  public properties. `cNumeric()` is kept as a float convenience view.
- **`BernExpMp` stores `p = 0` where dgs leaves memory uninitialised.** `dgs_bern.c:121-124` breaks
  out of the table build when `exp()` underflows but still sets `l = i+1`, so index `i` may be read
  while `p[i]`/`B[i]` were never initialised — undefined behaviour in C. We store an explicit zero,
  which deterministically rejects.
- **`sigma2+logtable` rejects `sigma < sigma_2`.** dgs computes `k = round(sigma/sigma_2)` with no
  lower guard, so `k = 0` makes `mpz_urandomm` divide by zero. We throw `ValueError` naming
  `sigma_2 = 0.849322`. An added error, not a changed result.
- **The lattice sampler's constructor requires `sigma`** where upstream defaults it to 1
  (`(B, sigma=1, c=None, r=None, precision=None, sigma_basis=False)`), and invents three messages
  Sage does not have: `TypeError: sigma is required` (`discrete_gaussian_lattice.ts:1164`),
  `'sigma must be a finite number, got X'` (`:1174`) and `'sigma must be > 0, got X'` (`:1177`).
  **This should be fixed** by defaulting `options.sigma` to 1; it is recorded here rather than in
  Part II because it is a one-line signature change inside an otherwise-faithful entry.

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
| `from_png(filename)` / `to_png(A, filename)` | File I/O | Retained but throw `NotImplementedError` with guidance to the data-based alternatives (`matrix_mod2.ts:1043-1049`) |
| Substitutes | — | `from_png_data(width, height, pixels)` and `to_png_data(A)` (`:1007-1028`, `:1064-1088`), with the black = 1 / white = 0 convention SageMath's libgd path uses (`matrix_mod2_dense.pyx:2613-2700`) |
| Affected modules | `matrix/matrix_mod2_dense.pyx` | `matrix/matrix_mod2.ts` |

### Rationale

1. **No native image libraries** — JavaScript environments have no universal PNG library like libgd.
2. **Environment portability** — browser, Node.js and Deno handle PNG differently.
3. **Separation of concerns** — raw pixel data lets users choose their PNG encoder.

### Trade-offs

- No direct file I/O; users must handle PNG encoding/decoding themselves (e.g. `pngjs` on Node.js,
  the Canvas API in browsers).
- `to_png_data` throws `TypeError('cannot create image with dimensions 0 x 0')` where Sage's message
  is `cannot write image with dimensions {c} x {r}` (`matrix_mod2_dense.pyx:2684`). The wording
  should be corrected to Sage's; the test at `matrix_mod2.test.ts:488` pins the current text.

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
| Pollard-Brent's `n < 2^96` size gate | `ifactor1.c:1361-1367` declines there because MPQS is faster in that range | **Removed.** Restoring it would only change which stage of `ifac_crack` splits a composite, never the factorization returned, and it would invalidate the existing `pollardbrent` range tests. Rho is a complete algorithm, so only the time changes |
| Pollard-Brent's round budget | PARI's `c0` formula (`:1369-1372`) goes negative below ~60 bits, which PARI never evaluates because it declines first | Clamped to at least `tune` (14) |
| SQUFOF acceptance bound | 2^46 on 64-bit builds (`:1487`) | 2^59. The 2^46 cut-off is a tuning choice that hands 2^46..2^59 to MPQS; the algorithm and PARI's own comment (`:1492`) document validity for `5 < n < 2^59`. Measured 0.5 ms per 56-bit semiprime, 200/200 split |
| ECM | Up to 64 curves in parallel with Montgomery batched inversion, Montgomery's PRAC chain, and a stage 2 over a 48-entry helix of residue classes mod 210 with a 1024-entry baby-step table | Same curve family, `nbc`/`dsn`/seed derivation, `TB1`/`TB1_for_stage` schedules and B1 phase, but serial with a binary ladder and an additive stage 2. Detection is equivalent (`[p]Q` vanishes mod a prime divisor exactly when the denominator does). PARI's *insisting* mode loops forever; ours is bounded by `FactorOptions.ecmRounds` (default 4) |
| `is_357_power`'s residue sieve | mod 211/209/61/203, then 117/31/43/71 through a 106-entry mask table | Omitted — a pure speed filter; the exact k-th root is taken directly, keeping PARI's mask semantics. Verified against brute force exhaustively on 2..20000 |
| Perfect powers / `isprimepower` | `ispower.c`: `Z_issquareall`, `is_357_power`, `is_kth_power`, `is_pth_power`, `Z_isanypower_101`, `Z_isanypower`, `isprimepower` | All ported, using **exact integer k-th roots** (Newton) instead of PARI's `sqrtnr`/`mpexp` float guesses. `isprimepower` never factors `n` |
| `pollardbrent_i` on 4 unlucky restarts | `pari_err_BUG` (`:1330`) | Returns `null`. A can't-happen path; degrading to "this engine failed" cannot produce a wrong factorization |
| File placement | `ispower.c` | `isprimepower`, `Z_isanypower`, `Z_isanypower_101`, `Z_issquareall`, `is_kth_power`, `is_pth_power`, `is_357_power` and `Z_iroot` live in `ifactor.ts`; they should move to `parigp-ts/src/ispower.ts`. Each carries a JSDoc naming its upstream file and line |
| Affected modules | `pari/src/basemath/ifactor1.c`, `ispower.c`, `mpqs.c` | `packages/parigp-ts/src/ifactor.ts`, `mpqs.ts` |

### Rationale

1. **No floating point** — PARI's root guesses are heuristics followed by an exact check, so
   replacing the guess with an exact integer root cannot change the verified answer.
2. **Tuning thresholds are not values.** Where a PARI gate exists only to hand work to a different
   stage, moving it changes the running time and nothing else, and each such move is justified by a
   measurement rather than by preference.

### Trade-offs

- Memory ceiling on very large MPQS inputs, where PARI spools to disk.
- ECM's constant factor is worse than PARI's (serial curves, additive stage 2).
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
| GF(2) kernel | `F2Ms_ker` (`F2v.c:1063`): dense `F2m_ker_sp` for `nbrow <= 640`, randomized block Lanczos above | Always the dense transcription. Same kernel *space*, so no output can differ. Relation collection is the binding constraint at every reachable size — the 79-digit case spends 192 s collecting ~600 relations whose elimination is milliseconds |
| `Fl_sqrt` | Random search for a generator of the 2-Sylow, so *which* root comes back is not deterministic upstream either | Deterministic Tonelli-Shanks. Both roots give a correct polynomial family; verified exhaustively for every odd prime `p < 500` |
| Sieve inner loops | `mpqs_sieve_p`/`_p1`/`_p2` are 4x/8x unrolled and interleave the two progressions | Two plain loops per factor-base entry. The multiset of byte additions is identical, so the sieve array is bit-for-bit the same |
| `mpqs_eval_sieve` bit array | `__v2di` (16 bytes) with SSE2, else `ulong` | The 8-byte scalar layout. The threshold is always `>= 128`, so both collect exactly the bytes `>= threshold`, in increasing order |
| `relaprimes` / `relp` buffers | Fixed `MAX_PE_PAIR = 60`; a candidate with more distinct factor-base divisors overruns them | Sized to the factor base / growable. Strictly a bounds fix |
| Relation exponent packing | `pi \| (ei << 20)` in a 64-bit long, so `\|ei\| < 2^43` | The same packing in 32 bits, so `\|ei\| < 2^11`. Factor-mode exponents are bounded by `log2(4·A·Q(x))` (tens); negative exponents occur only in the unreachable class-group mode |
| Factor base layout | 32-byte union on 64-byte boundaries | Parallel typed-array columns; every field keeps its type, **including** the C `float` `fbe_flogp` (`Float32Array`), whose rounding participates in `mpqs_locate_A_range` and `mpqs_si_choose_primes` |
| `MpqsOptions.maxPolys` / `.debug` | No such knobs (`MPQS_DEBUG` is a compile-time `-D`) | Added. Defaults (`0`, `false`) reproduce upstream exactly; `debug` turns upstream's own `mpqs_check_rel` (`:1069`) and post-Gauss check (`:1525`) into a permanent test oracle |
| `mpqs_class_init` / `mpqs_class_rels` (`:1775`, `:1815`) | Present | **Absent.** Their only caller is `buch2.c`, which is not ported, so there is no oracle and no reachable call site. Every `MPQS_MODE_CLASSGROUP` branch of the *shared* routines is transcribed and is therefore **untested code** |

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
| `t_REAL` | `mp.c`'s kernel, with AGM/Newton for `logr_abs` | A **semantics emulation** — `{sign, normalized mantissa, exponent, bit precision}` with round-to-nearest, `sqrtr` by integer square root, `logr_abs` by argument reduction + the `atanh` series. The representation, `expo()`, the precision discipline and `truncr`/`gcvtoi` are PARI's. Accuracy is asserted by the algorithm itself, exactly as upstream: `get_R` (`buch1.c:996`) only accepts a regulator when `h·R·invhr` lies in `(0.8, 1.3)` |
| `ZM_pivots` | Modular rank profile (`Flm_pivots`) certified by exact linear algebra | One-step fraction-free (Bareiss) elimination with the "first unused row" rule — the same canonical row rank profile PARI certifies |
| Randomness | `pari_rand` | Seeded xorshift32 (`setBuchRandomSeed`). `no` and `cyc` are canonical (which is why PARI's golden values match); the **generators** can differ from PARI's for a given `setrand`, as PARI's own do |
| MPQS relation collection for `\|D\| >= 2^60` | `mpqs_class_init`/`mpqs_class_rels` | `use_mpqs` is permanently `false`, i.e. always upstream's own fallback `imag_relations` (`buch1.c:746`). Same relations, same class group; slower (25 digits 0.4–2.2 s, 34 digits 12 s) |
| `gcvtoi`'s error exponent for an exact integer | `expo(0) = -HIGHEXPOBIT` | `-2^30`. Only ever consumed as the predicate `e > 0` |
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
| `galoisborne` | Computes the complex roots with `QX_complex_roots` at `t_REAL` precision, forms the inverse Vandermonde numerically and takes its operator sup-norm | An **exact Hadamard/Cramer bound**: `borne = ceil(den·n^(n/2)·B^(n(n−1)/2) / floor(sqrt\|disc T\|))` with `B` Cauchy's root bound. Any *upper* bound is correct here — it only sizes the `l`-adic accuracy and a rejection threshold |
| `indexpartial` | Refines `p^e \|\| disc` with `ZpX_reduced_resultant_fast` | Keeps `p^(e/2)`. Still a multiple of the denominator of an algebraic integer in the power basis |
| "Combinatorics too hard" | `frobeniusliftall` and `testpermutation` print a warner and then give up or `return identity_perm(n)` — PARI can return a **wrong/partial group with a warning** | Both throw. There is no warning channel here, and the thresholds (10^15 and 10^14 tests) are never approached |
| Output certification | Trusts its p-adic bounds | An extra `certify` pass: every element must be a distinct permutation of the `l`-adic roots induced by its own polynomial, and the set must be closed under composition, else `PariBugError`. It never fired |
| `FpX_factor_squarefree` | Handles arbitrary input | Only the squarefree case (the only one `galconj.c` needs); non-squarefree input raises instead of looping |
| Root ordering / choice of `l` | Whatever the splitting algorithm produces; cyclotomic `T` short-circuits through `galoiscyclo` | Roots sorted increasingly; no cyclotomic short-circuit. The permutation **labels** and `gal.roots`/`gal.p` can differ from PARI's by a relabelling (e.g. `l = 17` vs PARI's 41 for `polcyclo(8)`); everything label-independent is identical. `galoisinit(x^4-x-1)` returns `null` where live PARI returns `0` |
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
| `lllgramint` | `ZM_lll(x, 0.99, LLL_IM\|LLL_GRAM)`, ~2600 lines of floating-point/`flatter` hybrid | Cohen Algorithm 2.6.3, exact rational Gram-Schmidt, same `delta = 0.99`. Only the *speed* of `qfrep0` depends on the reduction — representation numbers are invariant under any unimodular change of basis — and `qfrep0` recomputes `det(u)` and falls back to the unreduced form if it is not `±1` |
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
| `t_REAL` rounding | One guard word per primitive, rounded up when its top bit is set; `divrr` is a truncating long division with a partial correction — i.e. **not** correctly rounded | The exact result rounded to nearest, ties away from zero (the rule PARI's guard word implements), with PARI's per-function *output precision* rules preserved. On the first value of PARI's own regression file PARI's printed 38 digits are 2.4 ulp above the truth and ours are 0.47 ulp above, so the 38th printed digit can differ by one |
| `qfrpow` for `n <= -2` | Inverts the form and then passes the **signed** `n` to `qfr5_pow`/`qfr3_pow`, which invert again — so PARI returns `x^\|n\|`. Verified on a live PARI 2.15.4: `qfbpow(f,-6) == qfbpow(f,6)` | Inverts once and raises to `\|n\|`, so `x^-n` is the inverse of `x^n` |
| `qfr5_pow` exponent loop | Loops over the machine **words** of `n` with `if (m == 1 && i == 2) break;` and an arithmetic `m >>= 1`; a word with leading zero bits skips squarings it owed, and a word with its top bit set never terminates. `qfbpow(f, 2^64+1)` returns `qfbpow(f,3)`; `qfbpow(f, 2^63)` overflows the PARI stack | A plain right-to-left binary chain over the whole `bigint` exponent — identical for every single-word `n`. Our `f^(10^20)` distance is congruent mod the regulator to its form's cycle distance (residual `< 1e-20` in exact rationals); PARI 2.15.4's answer has residual 1057.8 with `R = 2641.55` |
| `qfrpowraw` distance sign | Forms the distance *after* negating `n`, so `x^-k` reports `+k·d`, contradicting its own form | Uses the original signed exponent. Ten of the 360 oracle values differ from PARI by exactly this sign flip |
| `qfr5_to_qfr`'s `mplog2(lg(d0))` | Passes a **word length** where a bit precision is expected — a call site missed when PARI 2.16 changed `prec` from words to bits (`Qfb.c:488` vs `:428`) | `mplog2(precision(d0))`. The branch fires whenever `fix_expo` has (reachable: `qfbpow([f,0.], 10^8)` fires it three times); with upstream's `lg(d0)` the distance would be wrong by ~1e5 |
| `qfr_1_fill`'s `subiu(y,1)` | Reads the container `y` as an integer — a typo | `y2 − 1`, which is what `qfr_1_by_disc` computes for the same discriminant. Unreachable from every public entry point, hence untested |
| `logr_abs` AGM branch | Taken above `LOGAGM_LIMIT` | Series branch always. `LOGAGM_LIMIT` is far above the 128–512 bits this module uses and both branches compute the same logarithm |

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

Two return types are unions, `SqrtInteger` and `UnsignedInfinityElement` are port-only value types,
and the special functions remain bounded by IEEE-754 output precision.

### Behavioral Impact

Values agree with the executable oracle on the covered domain. TypeScript callers must handle the
documented unions and port-only exact wrappers.

---
# Part II — Open Fidelity Gaps

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

## Number-Field Kernel Not Delegated to parigp-ts

| Aspect | CLAUDE.md requires | sagemath-ts |
|--------|--------------------|-------------|
| Where the PARI `nf` routines live | Delegation to `parigp-ts` | `packages/sagemath-ts/src/rings/number_field/pari_nf.ts` holds `nfbasis` (`:1185`), `nfdisc` (`:1224`), `nfgaloisconj` (`:1385`), `primedec`, `quadunit` and `polisirreducible`. `packages/parigp-ts/src/` has **no `nf` module at all** (its directory holds only `buch`, `elliptic`, `ff`, `ffinit`, `ifactor`, `matkermod`, `mpqs`, `polmodular`, `qfb`, `qfrep`) |
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

Three call sites compute quadratic class numbers by enumeration or by table, when the real
subexponential algorithm is already ported, exported and oracle-tested in this repo.

| Site | SageMath | sagemath-ts |
|------|----------|-------------|
| `class_group.ts` | PARI `bnfinit`/`quadclassunit` | Enumerates reduced primitive forms of the field discriminant (definite for `D < 0`, rho-cycles of indefinite forms for `D > 0`), composes with Dirichlet composition (Cohen 5.4.7) and reads the elementary divisors off the `\|G[p^k]\|` counts. Guarded by `CLASS_GROUP_DISC_BOUND = 2_000_000n` (`:501`), above which `:704` throws `SAGE_NOT_IMPLEMENTED: class group of discriminant D requires PARI bnfinit/quadclassunit`. `QuadraticField(10000003).class_number()` is 2 in Sage and throws here. `class_group.ts:489` still carries a stale comment claiming `parigp-ts` has neither routine |
| `cm.ts` | `D.class_number(proof)`, i.e. PARI `qfbclassno` | Counts reduced primitive positive-definite forms (`-a < b <= a <= c`, `b >= 0` when `a == c`, `gcd(a,b,c) = 1`), memoized in a module-level `h_dict` (`:359`, `:376`). `OrderClassNumber` is used for non-fundamental discriminants exactly as in Sage. `O(\|D\|)` rather than PARI's `O(\|D\|^(1/4))` |
| `integer_ring.ts` `Integer.class_number()` | `integer.pyx:5857-5862` | A **hardcoded lookup table** of 80 imaginary and 30 real discriminants (`knownImaginary` at `:818`, `knownReal` at `:903`), returning table values for those and throwing `NotImplementedError('class_number: computation for discriminant D requires PARI integration')` at `:949` otherwise. All 110 tabulated values are correct (0 mismatches against Sage), but a magic table is exactly what CLAUDE.md's Algorithm Fidelity section forbids. The docstring also states the **wrong convention** for `D > 0` ("this is the narrow class number"); `integer.pyx:5857-5862` says the opposite |

### Why this is still open

`parigp-ts/src/index.ts:529-531` **already exports** `Buchquad`, `quadclassunit0` and `quadclassno`
— the real McCurley-Buchmann index calculus, verified against all 608 discriminants of PARI's own
`test/32/quadclassunit` regression output and against `quadclassunit(±2^81+c)`. Executed:
`quadclassno(-10000003n)` returns `706n`, matching `QuadraticField(-10000003).class_number()`, and
`Buchquad(40000012n)` returns `no = 2`, `cyc = [2]` instantly. `class_group.ts` imports nothing from
`parigp-ts` (`:13-15`). This is a **wiring gap, not a capability gap**; the delegation was
deliberately not made in the same pass as the port.

### How to close

Import `quadclassno`/`quadclassunit0` in `class_group.ts`, `cm.ts` and `integer_ring.ts`; delete
`formClassGroup`, `reducedDefiniteForms`, `h_dict`, `CLASS_GROUP_DISC_BOUND` and the 110-entry table;
fix the `class_number` docstring's convention; re-run the `discriminants_with_bounded_class_number`
and `cm_j_invariants` doctests.
**Effort: about a day** — two or three call sites to rewire and ~300 lines to delete. The
`Integer.class_number()` half alone is **trivial**: one import and one call.

### Trade-offs of leaving it open

Discriminants beyond 2 000 000 raise; `Integer.class_number()` answers only for 110 hardcoded inputs;
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
| `isprimepower` / `Z_isanypower` / `Z_isanypower_101` | PARI `basemath/ispower.c`, via `Integer.is_prime_power` -> `__pari__().isprimepower()` | **Two copies.** The real port is `parigp-ts/src/ifactor.ts:659` (`Z_isanypower_101` at `:496`, `Z_isanypower` at `:534`) but is **not re-exported** from that package's `index.ts`, so `sagemath-ts/src/arith/misc.ts:819` keeps a local `isprimepower` with `anyPower101` at `:745` | Add the three names to `parigp-ts/src/index.ts`, import them in `arith/misc.ts`, delete `misc.ts:745-845`. **Effort: under an hour** — one export line plus ~100 deleted lines |
| `matfrobenius` (`alglin2.c:428-720`) | PARI | `packages/sagemath-ts/src/matrix/matrix_integer.ts:1440` `frobenius_form_integer`, with a docstring citing `alglin2.c:617`/`:688` | Add a matrix module to `parigp-ts` and move it. **Effort: medium** — the code is written and verified; it needs a home |
| `dilog`, `incgam` | PARI | `packages/sagemath-ts/src/rings/complex_mpfr.ts` | Add a transcendental-functions module to `parigp-ts`. **Effort: medium** |
| `qfgaussred` | `self.__pari__().qfgaussred()` | `quadratic_forms/quadratic_form__local_field_invariants.ts`, a line-for-line port of `alglin2.c:1650-1749`; `parigp-ts` exports only `qfgaussred_positive` | Move it into `parigp-ts` beside `qfgaussred_positive`. **Effort: small** — the output is already pinned against real PARI in a test, so the move is safe |
| NTL `SFCanZass` / `CanZass` factoring over GF(2) | NTL, via `polynomial_template.pxi` | `polynomial_gf2x.ts` keeps local `squareFreeDecomp` (`:552`), `equalDegreeFactorization` (`:642`) and `factor` (`:680`), because `ntl-ts`'s versions all throw | See [ntl-ts GF2X Factoring Stubs](#ntl-ts-gf2x-factoring-stubs) |
| `parigp-ts` error classes | — | `PariTypeError`, `PariDomainError`, `PariDimError`, `PariInvError`, `PariPrimeError`, `PariSqrtnError` and `PariFlagError` are defined once in `matkermod.ts:78-120` and imported by `ffinit.ts` and `qfb.ts`; `packages/parigp-ts/src/errors.ts` does not exist | Create `parigp-ts/src/errors.ts` re-exporting the same names. **Effort: trivial** |

### Why this is still open

File ownership: each fix landed in a work unit that did not own the dependency package.

### Trade-offs of leaving it open

Two live copies of `isprimepower` can drift. The dependency graph does not mirror SageMath's, so the
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
| `hilbert_symbol(a, b, p, algorithm)` | PARI (`algorithm='pari'`, the default) or a direct algorithm; `'all'` cross-checks and raises `RuntimeError` on disagreement; accepts rationals via `a = QQ(a).numerator()*QQ(a).denominator()` (`misc.py:4985-4986`) | `misc.ts:3546-3634` is a faithful transcription of `misc.py:4994-5026` (including Sage's `ValueError('p must be prime or -1')` at `:3556` and the archimedean case at `:3566-3568`), but `'pari'` and `'all'` both **fall through to the direct code** (`:3571-3577`), so there is no cross-check; rational input dies with a raw JS `TypeError: Invalid mix of BigInt and other type in remainder`. 200 pseudo-random `(a,b,p)` with `p ∈ {2,3,5,7,11,13}` agree with Sage's default: 0 mismatches | Port `hilbertii` (`reference/pari/src/basemath/arith1.c:587-612`, **25 lines**, `Z_pvalrem` + `kronecker`; `mphilbertoo` at `:579-585` is already reproduced) into `parigp-ts`, route `'pari'`/`'all'` through it with Sage's `RuntimeError`, and accept `Rational` via numerator·denominator. **Effort: under an hour** |
| `bernoulli(n, algorithm)` | Multiple backends (FLINT/Arb/PARI/bernmm) with heuristics; `B_1000` is instant | `misc.ts:2149-2187` is the **classical binomial recurrence** `B_m = -1/(m+1)·sum C(m+1,k)·B_k` over rationals with a gcd after every term — *not* the Akiyama-Tanigawa the comment at `:2135` claims, and the Von Staudt-Clausen comment at `:2113-2114` is dead text. `algorithm` is accepted and ignored ('pari'/'flint'/'arb'/'bernmm'/'gap' all return `-691/2730` for `n = 12`). The operative limit is not the safe-integer range but the `O(n^2)`-rational cost: measured `B_200` 231 ms, `B_400` 2.66 s, `B_1000` 76.1 s, so the practical ceiling is ~`n = 500`. Returns a plain `{numerator, denominator}` object, not the port's `Rational` class | Port FLINT's exact multi-modular path: `reference/flint/src/arith/bernoulli_number_denom.c` (68 lines, von Staudt-Clausen) + `bernoulli/mod_p_harvey.c` + `fmpq_ui_multi_mod.c` + CRT. Avoids the Arb/zeta path entirely and is exact integer arithmetic, so it satisfies CLAUDE.md's no-floating-point rule. **Effort: moderate, ~400 lines.** Also fix the algorithm comment and the return type |
| `dedekind_sum(p, q, algorithm)` | FLINT by default; `reference/flint/src/fmpq/dedekind_sum.c:80` returns **0** whenever `fmpz_cmp_ui(k,2) <= 0`, i.e. for every `k <= 0` | `misc.ts:4669-4680` takes `k = \|q\|` and returns `s(p, \|q\|)`. A 2400-pair sweep (`p ∈ [-8,40] × q ∈ [-6,40]`) found **109 mismatches, all with `q < 0`**; for `q >= 0` the port agrees with Sage's FLINT default on every pair, **including** every non-coprime pair and the large doctests `3^54−1`, `2^93+1`. (The previously documented "non-coprime inputs may differ from PARI" divergence does **not** occur: `pari` vs `flint` in Sage over `q ∈ [1,40)`, `p ∈ [-20,40)` gives 0 differences.) `algorithm` is ignored; returns `{numerator, denominator}` | Add `if (q <= 2n) return 0` before the gcd reduction — one line, matching `fmpq_dedekind_sum.c:80`. **Effort: trivial.** Full parity (the machine-word branch + the continued-fraction branch) is ~120 lines |
| `gauss_sum(char_value, finite_field)` | Generic over `char_value.parent()`; the doctests feed it `UniversalCyclotomicField().zeta(q-1)` | `misc.ts:4793-4853` is a near-verbatim transcription of `misc.py:6428-6446` (same `resu += zq_power * zeta_p_powers[gen_power.trace().lift()]`, same `gen_power *= gen`, `zq_power *= zeta_q`) — **not** a "numeric-only fallback". The real blocker is that **no ring in this repo satisfies the `CharacterValue` interface** at `:4760-4764`: there is no `UniversalCyclotomicField` anywhere in `packages/sagemath-ts/src`, so the function has zero call sites and zero tests, violating CLAUDE.md's mandatory property-test rule. Two undocumented hazards: `:4838` `if (resu.add && zq_power.mul)` silently **skips terms** and returns `ring.zero()` when a ring lacks the optional methods, and `:4817` seeds `zq_power` from `ring.zeta(1n).powers(1)[0]` where Sage uses `ring.one()` | Make the missing-method branch **throw** (immediately — it violates the [Honest Failure](#honest-failure-instead-of-silent-approximation) policy), fix the `zq_power` seed, and either port a cyclotomic ring or delete `gauss_sum` until a usable ring exists. `reference/sage/src/sage/rings/universal_cyclotomic_field.py` is vendored (1740 lines) but is libgap-backed, so it is not a straight port. **Effort: minutes for the throw; large for the ring** |
| `algebraic_dependency(z, degree, height_bound)` | `misc.py:210-263` handles `RealField`/`ComplexField` **entirely in Sage** with `M.LLL(delta=.75)`, falling through to `pari(z).algdep(degree)` only at `:269-270` for other types (p-adics). The irreducible-factor selection at `:272-274` applies to **every** branch. Returns `None` when the height bound is violated (`:258`) | `arith/misc.ts` matches Sage's real path including `delta = 0.75` (`:2031`) — verified: `algdep(1.888888888888888, 1)` -> `9x−17`, `algdep(sqrt2, 2)` -> `x^2−2`, `algdep(phi, 2)` -> `x^2−x−1`. Two divergences: the irreducible-factor selection is performed only in the **duplicated** complex branch inside `complex_mpfr.ts:1201-1215` (with `prec = min(prec,53)−6`), not in `arith/misc.ts`; and `:2054` **throws** `ValueError('no polynomial found within height bound')` where Sage returns `None`. (The previously documented "no PARI binding for algdep" rationale was wrong: Sage does not use PARI on this path.) No p-adic support | Move the complex branch's irreducible-factor selection into `arith/misc.ts` so both callers share it, and return `null` rather than throwing on a height-bound miss. **Effort: small** |
| `Integer.class_number()` | See [Quadratic Class Numbers](#quadratic-class-numbers-not-delegated-to-buchquad) | — | — |
| `Integer.__invert__()` | `~Integer(5)` returns `1/5` with parent `Rational Field` | `integer_ring.ts:1885-1894` returns `this` for `±1` and throws `ArithmeticError` otherwise, with an in-code comment "Since we don't have Rational type here" that is **false**: `rings/rational.ts` exists and is used throughout, and does not import `integer_ring`, so wiring it is safe. (This was previously recorded as "mathematical — integers are not invertible in ZZ", which is not what SageMath does.) | Return `new Rational(1n, this.value)` per `integer.pyx.__invert__`. **Effort: trivial** |
| `number_of_partitions(n)` | FLINT, handles large `n` | `integer_ring.ts:2001` throws above `n = 10 000`; `p(10000)` matches `Partitions(10000).cardinality()` | Replace the `O(n^2)` DP with FLINT's exact pentagonal-number recurrence, `reference/flint/src/arith/number_of_partitions_vec.c` (~60 lines, no Arb needed), and delete the cap. **Effort: small.** Matching Sage for astronomically large `n` (HRR via Arb) stays out of reach and is the only part worth documenting |
| `prime_pi(n)` | primecount / FLINT `n_prime_pi` | `integer_ring.ts:2073` throws above `10^7`; `:2077-2082` is a per-integer `_is_prime` loop — exactly the naive pattern CLAUDE.md's Algorithm Fidelity section forbids (433 ms at 10^6, so ~6 s at the cap). `prime_pi(100000) = 9592` matches Sage | Replace with FLINT's sieve: `reference/flint/src/ulong_extras/prime_pi.c` + `prime_pi_bounds.c` (~200 lines). **Effort: small-moderate.** Note the module attribution: Sage has these in `combinat/partition.py:9762` (via `libs/flint/arith_sage.pyx:175`) and `functions/prime_pi.pyx:39` (primecountpy); `prime_pi` is **not** an `Integer` method in Sage at all |

### Trade-offs of leaving these open

Rational Hilbert symbols, large Bernoulli numbers, negative-modulus Dedekind sums, any Gauss sum,
integer inversion, large partition counts and large prime counts are unavailable or wrong.
`gauss_sum` additionally violates the port's own honest-failure policy.

### Behavioral Impact

`dedekind_sum(5, -7)` gives `-1/14` where Sage's FLINT default gives `0` (and its PARI option gives
`1/14`). `algebraic_dependency` throws where Sage returns `None`. `Integer(5).__invert__()` throws
where Sage returns `1/5`. `gauss_sum` can silently return `ring.zero()`. Everything else raises.

---

## Finite Fields — Conway Table and Constructor Algorithms

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Conway polynomial availability | The `conway_polynomials` package: every `(p, n)` in Frank Luebeck's `CPimport.txt` | A hand-pasted subset in `conway_polynomials.ts:36` `CONWAY_LOW_COEFFICIENTS`: p = 2 to n = 64, p = 3 to 24, p = 5 to 18, p = 7 to 14, p = 11/13 to 12, p = 17/19/23/29/31 to 10. `has_conway_polynomial` returns false where SageMath's does not — Sage's own `exists_conway_polynomial` is `True` for `37^2`, `97^2`, `787^5`, `2^64`, `2^65`, `2^100`, `19^21`, `3^24`, `3^25`, none of which we have |
| Default modulus for `GF(p^n)` | `irreducible_element` (`polynomial_ring.py:3560-3626` -> `:2628-2681`): `n == 1` -> `x − 1`; Conway if available; NTL `GF2X_BuildSparseIrred` for `p = 2`; else PARI `ffinit` | **Identical branch chain**, delegating to `ntl-ts` and `parigp-ts`. Outside the tabulated range the modulus is `ffinit`'s / `BuildSparseIrred`'s rather than the Conway polynomial Sage would use, so element representations are not interoperable with SageMath's and the generator need not be primitive |
| `modulus=` / algorithm keyword on the exported `GF`/`FiniteField` | `GF(p^n, 'a', modulus='minimal_weight')` | **Not accepted.** `FiniteField = GFExtended` has signature `(q, variableName)` (`finite_field_extension.ts:1227`), so passing Sage's keyword silently becomes a *variable name* — the generator then prints as `[object Object]^8 + …`. A live bug. The algorithm strings are reachable only via `GFpn` / `new FiniteFieldExtension` |
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
- **`modulus=` keyword:** widen `GFExtended`'s signature to accept an options object carrying
  `modulus`/`name`, and reject a non-string where a name is expected instead of stringifying it.
  **Effort: small.**

### Trade-offs of leaving these open

`GF(p^n)` element representations diverge from SageMath's outside the tabulated range; a caller
passing Sage's `modulus=` keyword gets a nonsense generator name with no error.

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
| `factor(0)` | `ArithmeticError` | `ValueError` with the identical message — contradicting the general claim that error names match SageMath |
| Term orders | `TermOrder` supports thirteen orders (lex, invlex, deglex, degrevlex, neglex, negdegrevlex, negdeglex, degneglex, wdeglex, wdegrevlex, negwdeglex, negwdegrevlex) plus block orders and weighted gradings | `multi_polynomial_element.ts:29` declares `export type TermOrder = 'lex' \| 'deglex' \| 'degrevlex'`; `:52` throws `unknown term order '<name>'`, matching `term_order.py:796` exactly. Weighted grading and `degree(std_grading=…)` have no meaning here. Roughly 19 of Sage's `MPolynomial` methods are honest `NotImplementedError` stubs |
| `MPolynomialRing` over ZZ | Works | **Cannot be constructed at all.** `MPolynomialRingConstructor(ZZ, ['x','y'])` dies inside `gens()` with a raw JS `TypeError: coeff.isZero is not a function` (`multi_polynomial_ring.ts:154`) before any deviation check runs. The only "ZZ" exercised in the test suite is a hand-rolled mock (`multi_polynomial_ideal.test.ts:130` `makeZZ`) |
| `groebner_basis()` over Zmod(n) | `MPolynomialIdeal.groebner_basis` works over ZZ and Zmod(n) via Singular's `std` (doctests at `multi_polynomial_ideal.py:4593-4610`) | Raises Sage's `TypeError('Can only reduce polynomials over fields.')` (`multi_polynomial_ideal.ts:124`, matching `multi_polynomial_element.py:2487-2488`). Without the field guard the quotient coefficient rounds to zero, the subtrahend is zero, and the reduction loop never terminates (reproduced as a timeout) — so the guard is correct *given* the missing Singular backend, but Sage does answer here. `dimension()` over a non-field raises `NotImplementedError('implemented only over fields')` |
| Gröbner engine | Singular/FGb with F4/F5-style algorithms | Naive Buchberger with an iteration budget (see [Bounded Search Budgets](#bounded-search-budgets-and-measured-thresholds)). `dimension()` is the ported Cox-Little-O'Shea algorithm (`multi_polynomial_ideal.py:1128-1192`); all five of Sage's `dimension()` doctests reproduce (`1, -1, 1, 1, 2`) |

### How to close

- **Printing:** fold negative coefficients into the sign in `Polynomial.toString`, emitting Sage's
  `' - c*x^k'` form. **Effort: small**, and it removes a false claim elsewhere in this document.
- **Integer `factor()` unit and `factor(0)` error class:** move the unit out of the list (or document
  it in the same row as the polynomial case) and throw `ArithmeticError`. **Effort: trivial.**
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
| Parent identity | `UniqueRepresentation`: `GF(7) is GF(7)` is `True`, and Sage's coercion framework is built on `parent is parent` checks | `GF(7) === GF(7)` is `false`. Repeated constructor calls return distinct, non-identical parents, so `is`-style comparisons that Sage relies on cannot be reproduced. `PowerSeriesRing` works around it with a structural `is_identical_to` (see [Power Series](#power-series-laurent-series-and-multivariate-series)) |
| Caching generally | Extensive (`@cached_method`, `UniqueRepresentation`) | Minimal (zero/one caching, selective caches) |

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
| `minimal_element()` representative | PARI `qfminim(q, NULL, NULL, 1)` | The same algorithm on `parigp-ts`'s `lllgramint`, which returns a different valid LLL transform on some inputs — see [Quaternion Algebras](#quaternion-algebras). The correct fix belongs in `parigp-ts` |

### Why this is still open

The order/ideal/maximal-order machinery this module targets (Deuring correspondence, isogeny-based
crypto) is QQ-only in SageMath too. Supporting other base rings needs number-field ideals, real
embeddings and PARI `alginit`, none of which is ported.

### How to close

Port number-field ideals and real embeddings (largely the same `nf` layer as
[Number Field Class Groups](#number-field-class-groups-units-and-galois-closure)) and PARI's
`alginit` (`reference/pari/src/basemath/alglin*.c`). For `minimal_element`, make `lllgramint` return
PARI's transform. **Effort: large for the base rings; small for `lllgramint` once `parigp-ts` is
touched.**

### Trade-offs of leaving it open

No GF(p) or number-field quaternion algebras, hence no `ramified_places` over number fields.

### Behavioral Impact

None for QQ, hard error (never a wrong answer) elsewhere. `minimal_element` returns a different
element of minimal norm in 4 of 254 sweep cases; the norm is always the true minimum.

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
