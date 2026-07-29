# Deviations from SageMath

> **📋 This is a living document.** Anyone porting SageMath functionality MUST update this file when their implementation differs from SageMath behavior.

This document tracks **intentional differences** between `sagemath-ts` and the original SageMath implementation. Every deviation must include a rationale explaining why the difference exists.

## How to Use This Document

**For library users:**
- Review relevant sections before using a module
- Check if your use case is affected by documented deviations
- Deviations are linked from function docstrings via `@see Deviation:`

**For contributors:**
- Update this document when your implementation differs from SageMath
- Use the [template](#template-for-new-deviations) at the bottom
- Link to this document from affected function docstrings

**Finding functions affected by deviations:**
```bash
# Find functions with deviation notes
grep -r "@see Deviation" packages/
```

## Table of Contents

1. [No JavaScript Number Coercion](#no-javascript-number-coercion)
2. [Arbitrary Precision Integers](#arbitrary-precision-integers)
3. [Error Handling - Null Returns vs Exceptions](#error-handling---null-returns-vs-exceptions)
4. [Algorithm Implementation Choices](#algorithm-implementation-choices)
5. [Tower Field Extension Degree Limitation](#tower-field-extension-degree-limitation)
6. [Extended Rational Rounding Modes](#extended-rational-rounding-modes)
7. [Empty Collection Identity Values](#empty-collection-identity-values)
8. [Type System Adaptations](#type-system-adaptations)
9. [Random State and Seeding](#random-state-and-seeding)
10. [Number Field Implementation Without PARI](#number-field-implementation-without-pari)
11. [p-adic Number Implementation](#p-adic-number-implementation)
12. [Unimplemented Number-Theoretic Functions](#unimplemented-number-theoretic-functions)
13. [Gauss Sum Simplified (numeric-only)](#gauss-sum-simplified-numeric-only)
14. [Hilbert Symbol Direct Algorithm Only (integer-only)](#hilbert-symbol-direct-algorithm-only-integer-only)
15. [Bernoulli Numbers (single algorithm, size limits)](#bernoulli-numbers-single-algorithm-size-limits)
16. [Dedekind Sum Algorithm Differences](#dedekind-sum-algorithm-differences)
17. [Elliptic Curve p-adic L-series and Isogeny Class Partial Implementation](#elliptic-curve-p-adic-l-series-and-isogeny-class-partial-implementation)
18. [Elliptic Curve Torsion Over Number Fields Not Implemented](#elliptic-curve-torsion-over-number-fields-not-implemented)
19. [Elliptic Curve Isogeny Algorithms Limited](#elliptic-curve-isogeny-algorithms-limited)
20. [GF(2) Matrix PNG Functions](#gf2-matrix-png-functions)
21. [Language & API Adaptations](#language--api-adaptations)
22. [Return Type Differences in Arithmetic](#return-type-differences-in-arithmetic)
23. [Ring/Field Iteration and Coercion](#ringfield-iteration-and-coercion)
24. [Polynomial Representation Differences](#polynomial-representation-differences)
25. [Polynomial Roots and Factorization Limited](#polynomial-roots-and-factorization-limited)
26. [Integer Polynomial Factorization Simplified](#integer-polynomial-factorization-simplified)
27. [Groebner Basis Algorithms Simplified](#groebner-basis-algorithms-simplified)
28. [Generic Group API and DLP Limitations](#generic-group-api-and-dlp-limitations)
29. [Finite Field Constructors and Display](#finite-field-constructors-and-display)
30. [Conway Polynomial Database Limited](#conway-polynomial-database-limited)
31. [Finite Field Extension Minimal Polynomial Simplified](#finite-field-extension-minimal-polynomial-simplified)
32. [Error Class Parity](#error-class-parity)
33. [Caching and Import Paths](#caching-and-import-paths)
34. [Elliptic Curve Short Weierstrass Form Only](#elliptic-curve-short-weierstrass-form-only)
35. [Real/Complex Numerical Approximations](#realcomplex-numerical-approximations)
36. [Matrix and Lattice Algorithm Simplifications](#matrix-and-lattice-algorithm-simplifications)
37. [Shortest Vector Problem (SVP) Implementation](#shortest-vector-problem-svp-implementation)
38. [Closest Vector and Voronoi Cell](#closest-vector-and-voronoi-cell)
39. [Algebraic Dependency Approximation](#algebraic-dependency-approximation)
40. [Combinatorial Function Limits](#combinatorial-function-limits)
41. [Real Matrix Decompositions (SVD, QR, LU) Using IEEE 754 Doubles](#real-matrix-decompositions-svd-qr-lu-using-ieee-754-doubles)
42. [PARI Factorization Algorithms Limited (parigp-ts)](#pari-factorization-algorithms-limited-parigp-ts)
43. [PARI Elliptic Curve Advanced Algorithms Missing (parigp-ts)](#pari-elliptic-curve-advanced-algorithms-missing-parigp-ts)
44. [Discrete Gaussian Samplers](#discrete-gaussian-samplers)
45. [PARI/NTL Routines Ported In Place Instead of Delegated](#parintl-routines-ported-in-place-instead-of-delegated)
46. [Infinity Representation](#infinity-representation)
47. [Exact Arithmetic Where SageMath Uses Floating Point](#exact-arithmetic-where-sagemath-uses-floating-point)
48. [Upstream Behaviour Deliberately Not Reproduced](#upstream-behaviour-deliberately-not-reproduced)
49. [Honest Failure Instead of Silent Approximation](#honest-failure-instead-of-silent-approximation)
50. [Port-Only APIs With No SageMath Counterpart](#port-only-apis-with-no-sagemath-counterpart)
51. [Keyword Arguments, Return Shapes and Signature Adaptations](#keyword-arguments-return-shapes-and-signature-adaptations)
52. [Multivariate Polynomial Term Orders Restricted](#multivariate-polynomial-term-orders-restricted)
53. [Free Module Exactness and Coordinate Types](#free-module-exactness-and-coordinate-types)
54. [Matrix Module Algorithm Substitutions](#matrix-module-algorithm-substitutions)
55. [Binary Quadratic Forms](#binary-quadratic-forms)
56. [Elliptic Curve and Isogeny Deviations](#elliptic-curve-and-isogeny-deviations)
57. [CM and Class Number Computation](#cm-and-class-number-computation)
58. [Coding Theory Deviations](#coding-theory-deviations)
59. [Crypto Module Deviations](#crypto-module-deviations)
60. [ZK Sumcheck and Multilinear Extensions](#zk-sumcheck-and-multilinear-extensions)
61. [ntl-ts GF2X Representation](#ntl-ts-gf2x-representation)
62. [Matrix Special Constructors](#matrix-special-constructors)
63. [Bounded Search Budgets and Measured Thresholds](#bounded-search-budgets-and-measured-thresholds)
64. [Template for New Deviations](#template-for-new-deviations)

> **Note on scope:** entries 44-61 were added in 0.0.11 as part of the July 2026 audit
> (`AUDIT-2026-07.md`) and the fix pass that followed; entries 9, 10, 11, 17, 19, 25, 26, 27, 28,
> 29, 30, 36, 37, 38, 39 and 43 were **corrected** in the same pass because they described
> behaviour the port no longer has — in several cases claiming a function throws when it did not,
> or claiming parity that did not hold.
>
> **0.0.12 (deferred-work pass):** entries 62 and 63 are new. Entries 9, 10, 17, 26, 30, 36, 42,
> 43, 44, 45, 49, 50, 51, 53, 54, 55, 57 and 59 were **rewritten**, in most cases because the
> deviation they described no longer exists: a dependency that was stubbed (PARI `ffinit`,
> `matkermod`, the `Qfb` family, real integer factorisation, the elliptic formal group) has since
> been implemented, so the workaround it forced is gone. Where a gap remains it is now stated in
> terms of the *specific* upstream routine that is still unported, not as a blanket limitation.

---

## No JavaScript Number Coercion

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Integer input types | Python `int` (arbitrary precision) | `bigint` or `Integer` only |
| JavaScript `number` acceptance | N/A | **Rejected with TypeError** |
| Affected modules | All functions accepting `IntegerLike` | `types/coercion.ts`, all integer APIs |

### Rationale

1. **Precision safety** - JavaScript `number` is IEEE 754 double-precision float, which silently loses precision for integers > 2^53-1. In Python, `int` is already arbitrary precision, so coercion is safe.
2. **Cryptographic correctness** - This library targets cryptographic applications where silent precision loss could cause security vulnerabilities or incorrect proofs.
3. **Source code hazard** - Even if we validated at runtime, a literal like `9007199254740993` in source code already loses precision before reaching our code (becomes `9007199254740992`).

### Example of the Problem

```typescript
// This looks correct but silently corrupts the value!
const p = 9007199254740993; // Actually becomes 9007199254740992
gcd(p, 2n); // Wrong input, wrong result

// Safe: use bigint literals
const p = 9007199254740993n; // Correct value preserved
gcd(p, 2n); // Correct
```

### Trade-offs

- Less convenient than SageMath's flexible coercion
- Users must write `123n` instead of `123` for all integer inputs
- Existing JavaScript code using `number` must be updated

### Mitigation

- Clear error message explains the issue and suggests using bigint literals
- TypeScript compiler catches type mismatches at compile time
- The `n` suffix is a minor inconvenience compared to the risk of silent data corruption

### Behavioral Impact

- Functions throw `TypeError` when passed JavaScript `number` instead of silently accepting
- All integer literals in user code must use the `n` suffix (e.g., `123n`)
- This is a **breaking change** from typical JavaScript conventions but essential for correctness

### Bigint to Number Conversion

When internal code must convert `bigint` to `number` (e.g., for floating-point math or array indexing), use `toSafeNumber()` which throws `RangeError` if the value exceeds safe integer range:

```typescript
import { toSafeNumber } from './types/coercion.js';

// Safe: throws if value too large
const idx = toSafeNumber(bigintValue);

// Unsafe: silently loses precision
const idx = Number(bigintValue); // Don't do this!
```

---

## Arbitrary Precision Integers

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Library | GMP (via MPIR, FLINT, NTL, Cython `mpz`) | Native TypeScript `bigint` |
| Affected modules | All integer arithmetic | All integer arithmetic |

### Rationale

1. **Zero dependencies** - No WebAssembly or external libraries required; works in any modern JS environment (browsers, Node.js, Deno, Bun)
2. **Engine-optimized** - V8/SpiderMonkey implement arbitrary precision arithmetic with highly tuned C++ internally
3. **Simpler codebase** - More readable and maintainable TypeScript code
4. **Sufficient for target use cases** - Cryptographic operations typically use 256-512 bit numbers, not thousands of digits

### Trade-offs

GMP would be faster for:
- Very large numbers (thousands of digits) due to FFT-based multiplication (Karatsuba, Toom-Cook, Schönhage–Strassen)
- Specialized operations with hand-optimized routines (`mpz_powm`, `mpz_probab_prime_p`, etc.)
- Architecture-specific assembly optimizations

### Mitigation

If performance becomes critical for specific operations, WebAssembly bindings like `gmp-wasm` could be introduced for hot paths without changing the public API.

### Behavioral Impact

None - arithmetic results are identical, only performance characteristics differ.

---

## Error Handling - Null Returns vs Exceptions

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Non-residue in `sqrt_mod()` | Raises `ValueError` | Returns `null` |
| Affected modules | `arith/misc.ts` | `arith/misc.ts` |

### Rationale

1. **TypeScript idioms** - Returning `null` for "not found" cases is idiomatic in TypeScript/JavaScript
2. **Type safety** - Union types (`bigint | null`) make callers handle the failure case explicitly
3. **No exception overhead** - Avoids try/catch blocks for expected failure cases

### Trade-offs

- Code that catches `ValueError` in SageMath needs to check for `null` instead
- Less informative than exception messages

### Mitigation

None needed - this is an intentional design choice. Callers should use nullish coalescing or explicit null checks.

### Behavioral Impact

Functions return `null` instead of throwing exceptions for mathematically undefined cases.

---

## Algorithm Implementation Choices

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| GCD algorithm | Euclidean (typically) | Stein's binary GCD |
| Integer square root | Various | Newton's method |
| Matrix determinant | Various (FLINT, etc.) | Bareiss algorithm (n > 3) |
| Modular square root | Tonelli-Shanks | Tonelli-Shanks |
| Affected modules | `arith/misc.ts`, `matrix/matrix_integer.ts` | Same |

### Rationale

1. **Performance** - Stein's algorithm avoids division, faster for bigints
2. **Numerical stability** - Bareiss avoids fractions in integer matrix operations
3. **Correctness** - All algorithms produce identical mathematical results

### Trade-offs

- Performance characteristics may differ from SageMath for specific inputs
- Users expecting particular algorithm behavior may be surprised

### Mitigation

None needed - these are implementation details. Results are mathematically identical.

### Behavioral Impact

None - outputs are identical, only internal computation differs.

---

## Tower Field Extension Degree Limitation

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `GF(2^n)` extension degree | Any positive integer | Powers of 2 only |
| Affected modules | Finite fields | `rings/finite_rings/tower_field.ts` |

### Rationale

1. **Tower construction** - Current implementation builds `GF(2^n)` as tower of quadratic extensions
2. **Simplicity** - Avoids complex irreducible polynomial selection for arbitrary degrees

### Trade-offs

- Cannot construct `GF(2^3)`, `GF(2^5)`, `GF(2^7)`, etc.
- Limits compatibility with cryptographic standards using non-power-of-2 extensions

### Mitigation

**TODO: Consider fixing.** Implement general extension field construction:
- Store precomputed Conway polynomials for common degrees
- Or implement irreducible polynomial generation

### Behavioral Impact

Throws error for non-power-of-2 extension degrees.

---

## Extended Rational Rounding Modes

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `Rational.round()` modes | Limited | 6 modes including `'odd'` |
| Affected modules | Rational numbers | `rings/rational.ts` |

### Rationale

1. **Completeness** - Support all IEEE 754 rounding modes plus extras
2. **Flexibility** - Different applications need different rounding behaviors

### Trade-offs

- `'odd'` mode is non-standard extension (rounds ties to nearest odd)
- Code using SageMath-specific rounding may behave differently

### Mitigation

None needed - this is a superset of SageMath functionality. Default mode (`'even'`) matches standard banker's rounding.

### Behavioral Impact

Additional functionality, not a breaking change. Default behavior matches expectations.

---

## Empty Collection Identity Values

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `gcd([])` | Error or 0 (context-dependent) | `0n` |
| `lcm([])` | Error or 1 (context-dependent) | `1n` |
| Affected modules | `arith/misc.ts` | `arith/misc.ts` |

### Rationale

1. **Mathematical identity** - `gcd()` identity is 0 (gcd(0, x) = x for all x)
2. **Mathematical identity** - `lcm()` identity is 1 (lcm(1, x) = x for all x)
3. **Consistency** - Follows reduce/fold conventions with identity elements

### Trade-offs

- SageMath behavior may differ in some contexts
- Users expecting errors for empty input will get silent "success"

### Mitigation

None needed - these are mathematically correct identity values.

### Behavioral Impact

Empty arrays return identity elements instead of raising errors.

---

## Type System Adaptations

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Type system | Duck typing | Static TypeScript types |
| Affected modules | All | All |

### Rationale

1. **Type safety** - TypeScript's static types catch errors at compile time
2. **IDE support** - Better autocomplete and documentation
3. **API clarity** - Explicit types document expected inputs/outputs

### Type Mapping Reference

| SageMath | TypeScript | Notes |
|----------|------------|-------|
| `Integer` | `bigint` | Native arbitrary precision |
| `Rational` | `Rational` class | Custom implementation |
| `bool` | `boolean` | Standard |
| `None` | `null` | **Not** `undefined` |
| `list` | `T[]` | Generic arrays |
| `tuple` | `[T, U, ...]` | Fixed-length tuples |
| `dict` | `Map<K,V>` or `Record` | Prefer `Map` for non-string keys |

### Trade-offs

- Some SageMath flexibility is lost (e.g., automatic type coercion)
- Function overloads needed for multiple signatures

### Behavioral Impact

Generally none for well-typed inputs. May throw `TypeError` for inputs that SageMath would coerce.

---

## Random State and Seeding

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| RNG core | `sage.misc.randstate` wraps `gmp_randstate_t` from `gmp_randinit_default` (GMP's MT19937) | `RandState`, a verbatim port of GMP 6.3.0 `rand/randmt.c` (`__gmp_mt_recalc_buffer`, `__gmp_randget_mt`, the 624-word `default_state` table) |
| Seeding | GMP `randseed_mt` (`rand/randmts.c`): seed mod 2^19937-20027, `+2`, `mangle_seed` (`r^1074888996 mod 2^19937-20023`), `mpz_export` into `mt[1..623]`, bit 19936 into bit 31 of `mt[0]`, 3 warm-up recalcs, `mti = 2000 % 624` | **Identical** — `mangle_seed` and `randseed_mt` are ported line for line, including GMP's non-canonical `reduce:` loop |
| Seed `0` | `randstate.pyx:556` skips `gmp_randseed` entirely when `seed` is falsy, so GMP's *unseeded* `default_state` buffer is used | Identical: seed 0 lands on the vendored `default_state` table, not on `gmp_randseed(0)` |
| `mpz_urandomm` | `mpz/urandomm.c`: bit length minus the power-of-two adjustment, 80-iteration rejection, then subtract `n`; `n = 1` returns 0 without drawing | Identical port |
| Second generator | `randstate.python_random()` returns a CPython `random.Random` seeded from `ZZ.random_element(1<<128)` | `PythonRandom`, a port of CPython `Modules/_randommodule.c` + `Lib/random.py` (`init_by_array`, `genrand_uint32`, `random_random`, `getrandbits`, `_randbelow_with_getrandbits`, `randrange`, `randint`, `normalvariate`) |
| `random()` | Module-level 31-bit integer; `randstate.c_rand_double()` for `[0,1)` | Both present; `RandState.random(): number` retained as a `@deprecated` alias for `c_rand_double()` |
| Affected modules | Any `.random_element()` / random sampling | Same |

### Rationale

1. **Bit-exact parity is now achievable** - GMP is not vendored under `reference/`, which is why
   0.0.11 used the canonical MT19937 seeding and accepted differing streams. For 0.0.12 the
   GMP 6.3.0 sources were obtained and `randmts.c` / `randmt.c` / `urandomm.c` were ported
   directly, so seeded streams are now **identical to SageMath's**, not merely
   same-distribution.
2. **Centralization** - Matches Sage's single randstate model
3. **Alias retained for callers** - `crypto/boolean_function.ts`,
   `rings/polynomial/polynomial_gf2x.ts` and `rings/polynomial/polynomial_element.ts` call
   `RandState.random()`; the alias delegates to `c_rand_double()`. A follow-up should rename the
   call sites and drop the alias.

### Trade-offs

- MT19937 is not cryptographically secure (neither is Sage's)
- A method (`RandState.random()`) exists that Sage's `randstate` does not have
- `python_random()` reproduces an upstream quirk that is easy to mistake for a bug: it derives
  its seed from `current_randstate()`, **not** from the receiver (`randstate.pyx:623` calls
  `ZZ.random_element`, which reads the global state). So
  `set_random_seed(0); randstate(314159).python_random().random()` is a *seed-0* value. A
  `seed` parameter is offered for callers who want a self-contained stream
- The GMP sources are not vendored under `reference/`, so the port cannot be re-diffed against
  them in-tree; the ported functions carry `gmp-6.3.0/<file>:<line>` citations instead

### Behavioral Impact

Seeded streams now match SageMath **exactly**, verified against a C oracle linked to the
installed libgmp 6.3.0 (7 seeds × 8 measurement families, all matching) and against SageMath
10.3 itself (`c_random()` after `set_random_seed(1207)` = 2008037228; `c_rand_double()` after
`set_random_seed(2718281828)` = 0.22437207488974298; the four `ZZ.random_element(10^30)`
doctest values; `python_random().random()` values; a 24-value `IntegerModRing(11)` stream).

Two consequences for callers that relied on 0.0.11 behaviour: `random_below` on a power of two
now consumes one bit fewer, `random_below(1)` consumes nothing, and negative seeds are reduced
mod 2^19937-20027 rather than absolute-valued.

Prior to 0.0.11 the generator took the *low* bits of a 64-bit MMIX LCG, so bit *k* had period
2^(k+1) and every low-bit consumer emitted a short deterministic cycle (`random_below(2n)`
alternated `1,0,1,0,…`). That is long fixed; as of 0.0.12 the remaining gap — a different
seeding step — is closed too.

---

## Number Field Implementation Without PARI

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Backend library | PARI/GP (nfinit, bnfinit) | Pure TypeScript |
| Affected modules | `rings/number_field/` | `rings/number_field/` |

### Rationale

1. **Zero dependencies** - No external PARI/GP bindings required
2. **Portability** - Works in any JS runtime without native dependencies
3. **Incremental approach** - Core functionality first, advanced features later

### Implemented Features

The following operations are implemented in pure TypeScript:

- **NumberField construction** - From defining polynomial over Q
- **Element arithmetic** - Addition, subtraction, multiplication, division, powers
- **Norm and trace** - Via characteristic polynomial computation
- **Characteristic polynomial** - Using direct formulas (n <= 3) or Faddeev-LeVerrier
- **Minimal polynomial** - Via GCD with derivative factoring
- **Signature** - Using Sturm's theorem for root counting
- **QuadraticField** - With discriminant computation
- **CyclotomicField** - With cyclotomic polynomial computation

Since 0.0.11 the module additionally carries a **number-field kernel**,
`rings/number_field/pari_nf.ts`, which ports the PARI routines SageMath delegates to:

- **`nfbasis` / `nfdisc`** - Pohst-Zassenhaus Round 2 (Cohen, Algorithm 6.1.8 — the algorithm
  behind PARI's `maxord`), giving real maximal orders, integral bases and field discriminants
- **`idealprimedec`** - Dedekind-Kummer prime decomposition **and**, as of 0.0.12, PARI's
  Buchmann-Lenstra round 4 (`base2.c:2248` `primedec_aux`, `:2150` `pradical`, `:2185`
  `pol_min`), so `decomposition(p)` branches exactly as PARI does and handles **inessential
  discriminant divisors**. Dedekind's classical example `x^3 - x^2 - 2x - 8` at `p = 2` (which
  previously threw) now splits completely into three primes with `e = f = 1`
- **`nfgaloisconj`** - as of 0.0.12, LLL-based reconstruction with **no degree cap**. Each
  conjugate is recovered from a single p-adic root by an exact integral lattice reduction, with
  two independent guarantees: `g(beta) = 0` is verified exactly (no false positives) and a
  Gram-Schmidt lower bound on `lambda_1` certifies non-existence (no false negatives), escalating
  precision when inconclusive. `numberofconjugates` (`galconj.c:3113`) is ported as an early exit
- **`quadunit` / `quadunitnorm`** - as of 0.0.12, PARI's `quadunit_uv_basecase` (`quad.c:281`),
  i.e. the continued fraction of `sqrt(D)` / the fundamental solution of Pell's equation. So
  `fundamental_units()`, `units()` and `regulator()` now answer for real quadratic fields
- **`polisirreducible`** - Cantor-Zassenhaus over F_p and Zassenhaus over Z

On top of that, quadratic class groups are computed exactly by enumerating binary quadratic
forms, and Galois groups are built from genuine automorphisms.

### Not Yet Implemented (still throw `NotImplementedError`)

- **Class group / class number for degree > 2** - Requires PARI `bnfinit`. One rigorous
  sub-case *is* implemented (0.0.12): every ideal class contains an integral ideal of norm at
  most the Minkowski bound `M_K`, so if **no** prime ideal has norm `<= M_K` the class group is
  provably trivial. `_classGroupIsProvablyTrivial()` computes a rational upper bound for `M_K`
  and decomposes every prime below it. This gives Sage's `h = 1` for `x^3-x-1`, `x^4-x-1`,
  `Q(zeta_5)`, `Q(zeta_7)`, `Q(zeta_12)`. Everything else — **including fields whose true class
  number is 1**, such as `Q(2^(1/3))`, `Q(zeta_8)` and `Q(zeta_23)` — throws rather than
  guessing, because we cannot *prove* the answer without deciding principality of the small
  primes. See [Bounded Search Budgets](#bounded-search-budgets-and-measured-thresholds)
- **Galois group of a non-Galois field** - SageMath returns the group of the Galois closure via
  `galoisinit`/`polgalois`; we throw when `|Aut(K/Q)| < [K:Q]` rather than fabricating a group.
  `GaloisGroup` methods needing `galoisinit` proper (`fixed_field`, `frobenius`, `inertia_group`
  above degree 2) likewise throw
- **Archimedean embeddings** - `embeddings()`, `real_embeddings()`, `complex_embeddings()`.
  These are what `bnfinit` would need (LLL reduction of ideal lattices under the `T2` form), so
  the two gaps are the same gap

### Trade-offs

- The kernel lives in `rings/number_field/pari_nf.ts` rather than in `parigp-ts`, which has no
  `nf` module at all; this is an **architectural** deviation from CLAUDE.md's delegation rule.
  If an `nf` module is added to `parigp-ts`, `number_field.ts` should be re-pointed at it and
  `pari_nf.ts` deleted
- Quadratic class groups are enumerated (exponential in `log|D|`) where PARI is subexponential;
  guarded by `CLASS_GROUP_DISC_BOUND = 2,000,000`, above which we throw rather than run slowly
- Zassenhaus factorisation over Z uses one big prime instead of a Hensel lift, and throws
  above 2^200 rather than producing a wrong factorisation
- Real embeddings used for the regulator are double precision (the regulator is transcendental,
  so no exact representation is available with the primitives here). `regulator()` uses a new
  overflow- and cancellation-safe `quadraticLogAbs` (log-sum-exp over bigint bit lengths), so the
  ~250-digit unit of `Q(sqrt(1000003))` gives a finite `R = 576.646` where a naive double
  evaluation returns `Infinity`
- Class-group generators come from enumerating reduced forms and greedily taking an element of
  maximal order, **not** from PARI's `quadclassunit`. The subgroup generated is always the whole
  class group (verified on 38 discriminants against the known structures), but the particular
  generators — and hence the particular primes derived from them — can differ from PARI's
- `nfgaloisconj` reaches the same proved answer as PARI but by `n` independent lattice
  reductions rather than from group generators, because `galoisgen`/`galoisinit` (~2000 lines of
  permutation-group machinery) is not ported. Slower on large degrees (`Q(zeta_25)`, degree 20:
  3.1 s), identical output
- `quadunit` uses only PARI's basecase, not the product-tree variant `quadunit_uv`
  (`quad.c:429`) that PARI switches to at `D >= 2000000`. Identical `[u,v]`; the product tree is
  a big-integer-multiplication speed optimisation (measured 0.5 ms for the ~800-digit unit of
  `D = 511681`)
- `primedec` returns `{gens, e, f}` rather than PARI's 5-component prime structure with a
  uniformizer and anti-uniformizer; `number_field.ts` then searches for a two-element
  representation `(p, alpha)` and **certifies** it with `N((p,alpha)) == p^f`, falling back to
  the full generating set. No wrong ideal can be returned; the printed form can have more
  generators than Sage's
- `pari_nf.ts` now depends on `matrix/matrix_integer.ts` for the exact integral LLL (no import
  cycle: `matrix_integer -> matrix_generic -> matrix_operations` imports nothing from
  `rings/number_field`)
- No proof flags for GRH-conditional results

### Mitigation

When parigp-ts adds nfinit/bnfinit support, we can:
1. Update NumberField to call parigp-ts for heavy computations
2. Delete `pari_nf.ts` and keep the TypeScript fallback for basic operations
3. Match SageMath's architecture exactly

### Behavioral Impact

- Basic element operations (norm, trace, arithmetic) match SageMath exactly
- Maximal orders, integral bases, field discriminants, prime decomposition, ideal arithmetic,
  automorphisms and quadratic class groups now **compute** and agree with PARI on every value
  tested; prior to 0.0.11 this section claimed all of them throw, which was wrong in the other
  direction — they were implemented but returning `Z[alpha]`-based approximations
- **Changed in 0.0.12:** `is_galois()` for degree > 2 used to return `false` unconditionally
  (because `automorphisms()` only ever returned the identity above degree 8); it now returns the
  true answer, and `galois_group()` consequently succeeds for Galois fields of any degree.
  `quadraticUnitGroup(K)` for a real quadratic `K` now carries a fundamental unit, so `ngens()`
  goes 1 -> 2 and `gens()` is `[-1, eps]` instead of `[-1]`
- A new exported helper `realQuadraticFundamentalUnit(K)` has no SageMath counterpart (Sage
  exposes the unit only through `K.units()`); it exists because the quadratic case bypasses
  `bnfinit` entirely, exactly as PARI's own `bnfinit` reduces to `quadunit` there
- The zero ideal is no longer reported prime (Sage's `idealismaximal` does not accept it at all)
- `NumberFieldElement.is_unit()` now implements Sage's *field* branch (true for any nonzero
  element); the ring-of-integers test is exposed separately as `is_integral_unit()`
- `QuadraticField(D)` uses `x^2 - D` verbatim as Sage does; `.D` is Sage's `D` and the
  port-only `.d` (squarefree part) is retained for legacy callers and should be treated as
  deprecated

---

## p-adic Number Implementation

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Backend | PARI/GP for some operations | Pure TypeScript |
| Precision model | Multiple (capped-rel, capped-abs, fixed-mod, floating-point) | Capped-relative only |
| Affected modules | `rings/padics/` | `rings/padics/` |

### Rationale

1. **Zero dependencies** - No PARI/GP bindings required
2. **Simplicity** - Capped-relative precision is most commonly used
3. **Performance** - Pure bigint arithmetic is sufficient for typical use cases

### Implemented Features

The following operations are implemented in pure TypeScript:

- **Zp, Qp construction** - p-adic rings and fields
- **Element arithmetic** - add, sub, mul, div, neg, pow, inv
- **Valuation** - valuation(), ordp(), unit_part(), normalized_valuation()
- **Precision** - precision_absolute(), precision_relative(), add_bigoh(), lift_to_precision()
- **Expansion** - expansion(), list(), residue(), __getitem__(), slice()
- **Predicates** - is_zero(), is_one(), is_unit(), is_integral(), is_square()
- **Square root** - `square_root()`/`sqrt()` ported from `padic_generic_element.pyx`, including
  SageMath's deterministic sign selection and the `extend`/`all` keywords
- **n-th roots** - `nth_root(n)` ported in full (`n = p^v * m` split, Newton on the inverse
  root, per-p-th-root precision loss); the residue-field root delegates to
  `Integer.nth_root_mod` (our port of `_nth_root_common`) and every factorisation goes through
  PARI's `Z_factor`
- **Teichmuller lift** - teichmuller() using Newton iteration
- **Log/Exp** - log() for units, exp() for convergent inputs
- **Orders** - `multiplicative_order()`, `additive_order()`
- **Expansion / slice / repr** - ported from `padic_template_element.pxi`,
  `local_generic_element.pyx` and `padic_printing.pyx`
- **Norm/Trace** - Trivial for base field (identity maps)

### Not Yet Implemented

- **Extension fields** - Zq, ramified/unramified extensions. `pAdicExtension` carries a degree
  but no element arithmetic, so `nth_root`'s p-th-root extraction is written only for
  absolute degree 1 (e = f = 1), where it reproduces SageMath's result exactly
- **minimal_polynomial(), charpoly()** - Meaningful only with extension field support

### Trade-offs

- Only capped-relative precision model is available
- SageMath uses PARI for gamma(), sqrt() of extension fields, etc.
- `square_root`/`nth_root` take an options object (`{extend, all}`) per DESIGN.md's keyword
  mapping, and additionally expose `square_root_all()` / `nth_root_all()` for well-typed
  access to the list form. The `algorithm` keyword is not offered: this port has a single
  implementation (SageMath's `'sage'` path), and the PARI path returns the same value
- `+Infinity` is the IEEE double `Number.POSITIVE_INFINITY` (typed `InfiniteOr<bigint>`)
  rather than `sage.rings.infinity.PlusInfinity`; see [Infinity Representation](#infinity-representation)

### Mitigation

When parigp-ts adds padic support, we can delegate operations like:
- gamma() - p-adic gamma function
- sqrt() for extension fields
- Extension field arithmetic

### Behavioral Impact

- Basic arithmetic matches SageMath for capped-relative elements. **Prior to 0.0.11 this claim
  was false**: `add()` multiplied by `p^v` twice for operands of unequal valuation, so most
  p-adic sums were wrong. Fixed by porting `CR_template.pxi:cadd`
- This section previously listed `nth_root` (n != 2), branch `log`, `artin_hasse_exp` and
  `minimal_polynomial`/`charpoly` as unimplemented; all four had real bodies at the time
- Precision semantics should match, but edge cases may differ
- Extension field operations throw NotImplementedError

---

## Unimplemented Number-Theoretic Functions

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Backend library | PARI/GP, FLINT, specialized algorithms | Pure TypeScript with parigp-ts |
| Affected modules | `rings/integer_ring.ts` | Same |

### Rationale

The following functions remain as stubs (throw `NotImplementedError`) because they require either:
1. PARI/GP functions not yet ported to parigp-ts
2. Complex algorithms beyond the current scope
3. Dependencies on other unimplemented modules (e.g., number fields, complex analysis)

### Unimplemented in `integer_ring.ts`

| Function | Reason | PARI Equivalent |
|----------|--------|-----------------|
| `class_number()` | Requires PARI's `qfbclassno` for binary quadratic form class number | `qfbclassno` |
| `__invert__()` | Integers are not invertible in ZZ; would need rational field | N/A (mathematical) |

### Trade-offs

- Some number-theoretic computations are not available
- Users relying on these specific functions will encounter errors
- Some cryptographic or research applications may be blocked

### Mitigation

These functions can be implemented when:
1. parigp-ts adds the corresponding PARI functions (e.g., `qfbclassno`, `bernfrac`)
2. The necessary supporting modules are complete (e.g., complex numbers, characters)
3. Community contribution provides specialized implementations

### Behavioral Impact

Functions throw `NotImplementedError` with descriptive messages like:
```
'SAGE_NOT_IMPLEMENTED: class_number'
```

---

## Gauss Sum Simplified (numeric-only)

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Gauss sum over finite fields | Full character/cyclotomic support via Sage libraries | Numeric-only fallback requiring minimal `add`/`mul`/`zeta()` interfaces |
| Affected modules | `sage/arith/misc.py` | `packages/sagemath-ts/src/arith/misc.ts` |

### Rationale

1. **Missing character/cyclotomic infrastructure** - Full Gauss sums require cyclotomic fields and character theory
2. **Pragmatic fallback** - Provide a minimal numeric implementation for common cases

### Trade-offs

- Not correct for general character values or non-numeric rings
- Limited to environments where `add`, `mul`, and `zeta()` are meaningfully defined

### Mitigation

Implement cyclotomic fields and multiplicative character support, then replace the fallback with a faithful implementation.

### Behavioral Impact

Results may be incomplete or incorrect for general characters/fields; numeric-only cases can work but do not guarantee SageMath parity.

---

## Hilbert Symbol Direct Algorithm Only (integer-only)

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `hilbert_symbol` backend | PARI (`algorithm='pari'`) or direct algorithm; accepts rationals | Direct algorithm only; integer inputs; `algorithm='pari'` maps to direct |
| Affected modules | `sage/arith/misc.py` | `packages/sagemath-ts/src/arith/misc.ts` |

### Rationale

1. **No PARI binding** - PARI's `hilbert` is not exposed in parigp-ts
2. **Scope control** - Integer-only implementation covers common cases

### Trade-offs

- Rational inputs are not supported
- `algorithm='pari'`/`'all'` does not cross-check against PARI

### Mitigation

Add PARI `hilbert` support to parigp-ts and extend inputs to rationals.

### Behavioral Impact

Integer inputs follow SageMath's direct algorithm; rational inputs and PARI-specific behavior are unavailable.

---

## Bernoulli Numbers (single algorithm, size limits)

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `bernoulli(n)` algorithms | Multiple backends (FLINT/ARB/PARI/bernmm) with heuristics | Single recurrence (Akiyama-Tanigawa); `algorithm` parameter ignored |
| Input range | Large `n` supported via FLINT/bernmm | Uses `toSafeNumber(n)` (limited to JS safe integer range) |
| Affected modules | `sage/arith/misc.py` | `packages/sagemath-ts/src/arith/misc.ts` |

### Rationale

1. **Dependency gap** - FLINT/ARB/NTL backends are not yet available in TypeScript
2. **Simplicity** - Recurrence is straightforward and correct for small/medium `n`

### Trade-offs

- Performance degrades rapidly for large `n`
- `algorithm` options are accepted but not honored
- `n` must fit into JavaScript safe integer range

### Mitigation

Implement FLINT/ARB/bernmm backends (or WASM bindings) and dispatch on `algorithm` as SageMath does.

### Behavioral Impact

Small `n` returns correct Bernoulli numbers; large `n` may be rejected or slow compared to SageMath.

---

## Dedekind Sum Algorithm Differences

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `dedekind_sum(p, q)` backend | FLINT by default; PARI option (differs for non-coprime) | Knuth’s algorithm for coprime inputs after reducing by `gcd(p, q)` |
| `algorithm` parameter | Selects FLINT/PARI | Accepted but ignored |
| Affected modules | `sage/arith/misc.py` | `packages/sagemath-ts/src/arith/misc.ts` |

### Rationale

1. **Dependency gap** - No FLINT/PARI backends yet
2. **Determinism** - Coprime reduction plus Knuth algorithm gives consistent results

### Trade-offs

- Results for non-coprime inputs may differ from SageMath’s PARI behavior
- No backend selection via `algorithm`

### Mitigation

Wire `dedekind_sum` to FLINT/PARI backends when available and preserve SageMath’s algorithm semantics.

### Behavioral Impact

For coprime inputs, results should match; for non-coprime inputs, behavior may differ from SageMath’s PARI path.

---

## Elliptic Curve p-adic L-series and Isogeny Class Partial Implementation

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Backend | Full modular symbols, PARI, p-adic fields | Partial: structure only |
| Affected modules | `schemes/elliptic_curves/padic_lseries.ts`, `schemes/elliptic_curves/isogeny_class.ts` | Same |

### Rationale

1. **Missing dependencies** - Full p-adic L-series computation requires modular symbols (sage.modular.modsym) and p-adic arithmetic (sage.rings.padics.Qp)
2. **Incremental approach** - Provide class structure and basic operations now, full implementation later
3. **API completeness** - TypeScript types and interfaces are fully defined

### Implemented Features

**padic_lseries.ts:**
- `pAdicLseries`, `pAdicLseriesOrdinary`, `pAdicLseriesSupersingular` classes
- `elliptic_curve()`, `prime()` - basic getters
- `is_ordinary()`, `is_supersingular()` - reduction type (in subclasses)
- `toString()` - string representation
- `rational()` helper - rational number creation with reduction
- Input validation for prime, implementation, normalize parameters
- Series parameter validation (n, prec, eta, quadratic_twist)

**isogeny_class.ts:**
- `IsogenyClass`, `IsogenyClassNumberField`, `IsogenyClassRational` classes
- `length()`, `get()`, `index()` - curve access
- `[Symbol.iterator]()` - iteration over curves
- `matrix()` - isogeny degree matrix with fill/unfill
- `reorder()` - reorder by 'lmfdb' or custom order
- `copy()`, `contains()` - utility methods
- `isogenies()` - partial (returns structure, not computed isogenies)

**padic_lseries.ts (additionally implemented and correct since 0.0.11):**
- `teichmuller(prec)` - implemented (this was previously listed as a stub; it never was one)
- `alpha(prec)` - correct for good-ordinary and multiplicative reduction
- `_e_bounds` - reproduces SageMath's doctest vectors exactly

**padic_lseries.ts (new in 0.0.12):**
- `bernardi_sigma_function(prec)` - implemented, a step-for-step port of
  `padic_lseries.py:1613-1641`. Reproduces SageMath's doctest for curve 14a verbatim:
  `z + 1/24*z^3 + 29/384*z^5 - 8399/322560*z^7 - 291743/92897280*z^9 + O(z^10)`. Independently
  cross-checked against the Weierstrass p-function recursion (`h'' = -sum c_n z^(2n)` with
  `c_1 = g2/20`, `c_2 = g3/28`) on four short-Weierstrass curves, coefficient for coefficient.
  The 0.0.11 entry describing this as throwing, and the source comment claiming
  `formal_group.log` is a placeholder, were both **stale**
- `alpha(prec)` for **supersingular** reduction - implemented, via two new classes
  `pAdicEisensteinQuadraticExtension` / `pAdicEisensteinQuadraticElement` that build the totally
  ramified quadratic extension `A = K[x]/(x^2 - a_p x + p)` SageMath constructs with
  `K.extension(f, names='alpha')` (`padic_lseries.py:513-518`). Reproduces Sage's 37a/`p = 3`
  doctest: `alpha(10)` prints `alpha + O(alpha^21)` and `alpha^2 - ap*alpha + 3` is
  `O(alpha^22)` and `is_zero()`

### Not Yet Implemented (throw NotImplementedError)

**padic_lseries.ts:**
- `modular_symbol(r, sign, quadratic_twist)` - requires sage.modular.modsym
- `measure(a, n, prec, quadratic_twist, sign)` - requires modular symbols + p-adic arithmetic
- `series(n, quadratic_twist, prec, eta)` - requires modular symbols + power series
- `order_of_vanishing()` - requires series computation
- `frobenius(prec, algorithm)` - requires Monsky-Washnitzer cohomology
- `Dp_valued_series` / `Dp_valued_height` / `Dp_valued_regulator` - depend on `series()` above
- `_c_bound()` - requires `E.galois_representation()` and modular-symbol denominators.
  It previously returned a literal `0`, which is **not** conservative: `c` is subtracted from
  the e-bounds, so `_prec_bounds` claimed more precision than justified

**isogeny_class.ts:**
- `isogenies(fill=true)` - requires full isogeny computation
- `qf_matrix()` - only for CM curves, requires CM discriminant
- `Frobenius_filter(E, L)` **over base fields other than QQ** - needs `K.primes_above(p)`,
  `E.reduction(P)` over a number field and a number-field global integral model, none of which
  exist (`sage.rings.number_field` has no `reduction`/`primes_above` here). Over QQ the filter is
  fully ported (see below)

### Implemented but heuristic (previously mis-listed here as throwing)

- `graph()` - implemented
- `isogeny_degrees_cm(E)` - **corrected in 0.0.12**. Two real defects were fixed: the
  horizontal-primes step (`isogeny_class.py:1309-1317`, a prime represented by each class group
  generator) was **absent entirely**, which made the returned list potentially *too small* — a
  soundness gap, not merely non-minimality — and the downward-ramified test used `n` instead of
  `n/(2h)`. SageMath's `d = -23` verbose trace is now reproduced line for line
  (`{2,3,5}` before filtering). `Frobenius_filter` (`gal_reps_number_field.py:492-586`) is now
  ported and applied wherever the base field is QQ, reproducing both of SageMath's own doctests
  (`11a1` -> `[5]`; `EllipticCurve_from_j(2268945/128)` over `[7,11]` -> `[7]`). Over a number
  field the filter is inapplicable, so `isogeny_degrees_cm` detects that and returns the
  **unfiltered** candidate set rather than a wrong one — so Sage's headline `d = -23` example
  (which lives over a degree-6 field) still gives `[2,3,5]` where Sage gives `[2,3]`. The
  function's contract ("this list is not necessarily minimal") holds either way
- `possible_isogeny_degrees(E)` - implemented, but over Q returns Mazur's list
  `[2,3,5,7,11,13,17,19,37,43,67,163]` (optionally intersected with the degrees for which
  `isogenies_prime_degree` finds an isogeny) instead of SageMath's Billerey/Larson bounds.
  Correct as a superset over Q; **not valid over larger number fields**

### Placement note

`Frobenius_filter` lives in `schemes/elliptic_curves/isogeny_class.ts`, not in a
`gal_reps_number_field.ts` (which does not exist in this port). It should move when that module
is ported. Likewise `pAdicEisensteinQuadraticExtension` belongs in `rings/padics/`, whose
`pAdicExtension` class is currently a shell with no element type, no `gen()` and no arithmetic.

### Also fixed in 0.0.12: `Frobenius_filter` good-reduction test

SageMath calls `E.has_good_reduction(p)` on the **minimal** model (Tate's algorithm via
`local_data`). Laska-Kraus-Connell minimisation is not ported, so we test `p` against the
discriminant of the *global integral* model built in place. A non-minimal model makes a few
extra primes look bad; those primes are skipped, which can only make the filter weaker (return a
superset), never unsound. Both of SageMath's doctests still come out exactly right, including
the deliberately non-minimal `EllipticCurve_from_j(2268945/128)` model.

### Trade-offs

- Cannot compute actual p-adic L-function values
- Cannot compute full isogeny classes (only trivial class with initial curve)
- The two heuristic isogeny-degree functions over-report candidates rather than throwing
- Users get clear NotImplementedError with explanation

### Mitigation

When the following dependencies are available, full implementation becomes possible:
1. `sage.modular.modsym` port - enables modular symbol computation
2. `sage.rings.padics.Qp` completion - enables p-adic arithmetic
3. `sage.rings.power_series_ring` - enables power series
4. Isogeny computation from `ell_curve_isogeny.ts` - enables full class computation
5. `gal_reps_number_field.py` / Billerey bounds - makes the isogeny-degree sets exact

### Behavioral Impact

- Basic class structure and getters work identically to SageMath
- Advanced methods throw `NotImplementedError` with descriptive messages
- `_prec_bounds` now propagates `_c_bound`'s `NotImplementedError` instead of over-reporting
- Test coverage validates implemented functionality

---

## Elliptic Curve Torsion Over Number Fields Not Implemented

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Torsion over number fields (incl. Q) | PARI/GP `elltors`, division polynomials, height methods | Throws `NotImplementedError` (finite fields only) |
| Affected modules | `sage/schemes/elliptic_curves/ell_torsion.py` | `packages/sagemath-ts/src/schemes/elliptic_curves/ell_torsion.ts` |

### Rationale

1. **Dependency gap** - PARI `elltors` is not exposed in parigp-ts
2. **Complexity** - Full number-field torsion requires division polynomials and height machinery

### Trade-offs

- Torsion subgroup over Q/number fields is unavailable

### Mitigation

Expose PARI `elltors` in parigp-ts or implement division-polynomial/height-based torsion computation.

### Behavioral Impact

Constructing torsion subgroups for characteristic 0 curves throws `NotImplementedError`.

---

## Elliptic Curve Isogeny Algorithms Limited

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Kernel polynomial algorithms | Kohel/BMSS/Stark; `compute_isogeny_kernel_polynomial` picks Stark for `ell < 10`, BMSS otherwise | Kohel and BMSS ported in full; **Stark throws** (`compute_isogeny_stark` needs `E.weierstrass_p()` from the unported `sage/schemes/elliptic_curves/ell_wp.py`), so BMSS is always used |
| BMSS differential equation | `fastElkies'` Newton doubling via `solve_linear_de` | Same ODE and initial condition `S = x + O(x^2)`, solved by a direct coefficient recurrence |
| Kernel polynomial validation | `is_kernel_polynomial` from `isogeny_small_degree.py` when `check=True` | Even-degree check ported with Sage's message; **odd-degree check skipped** (`isogeny_small_degree` is not ported) |
| Dual / rational maps | Full support | `dual()` and the pre/post-isomorphism plumbing are ported; `dual()` throws where Sage does, i.e. when `char < 4*deg + 4` |
| Affected modules | `sage/schemes/elliptic_curves/ell_curve_isogeny.py` | `packages/sagemath-ts/src/schemes/elliptic_curves/ell_curve_isogeny.ts` |

### Rationale

1. **Missing infrastructure** - `ell_wp.py` and `isogeny_small_degree.py` are not ported
2. **Equivalent results** - Stark and BMSS both compute the kernel polynomial of the unique
   normalized degree-`ell` isogeny, so the returned polynomial is identical; BMSS's
   precondition (`char >= 4*ell + 4`) is exactly the condition
   `compute_isogeny_kernel_polynomial` has already checked
3. **Uniqueness of the ODE solution** - the solution of `S'^2 = G(x, S)` with
   `S = x + O(x^2)` is unique, so the coefficient recurrence gives the same series as Newton
   doubling; it is `O(l^3)` rather than quasi-linear, which is negligible at the degrees
   reachable under `char >= 4l + 4`

### Trade-offs

- An explicit `algorithm='stark'` request throws instead of running
- A polynomial that is not a genuine odd-degree kernel polynomial produces garbage instead of
  SageMath's `ValueError`
- BMSS is asymptotically slower than upstream

### Mitigation

Port `ell_wp.py` (for `weierstrass_p`) and `isogeny_small_degree.py` (for
`is_kernel_polynomial`); replace the coefficient recurrence with `solve_linear_de`.

### Behavioral Impact

Isogeny *values* now match SageMath. Prior to 0.0.11 Vélu's y-coordinate carried a sign error
that put isogeny images off the codomain; see
[Elliptic Curve and Isogeny Deviations](#elliptic-curve-and-isogeny-deviations) for the
remaining intentional differences. `dual()` on small fields (e.g. GF(11), degree 7) raises
`NotImplementedError` — as SageMath does, via `compute_isogeny_kernel_polynomial`.

---

## GF(2) Matrix PNG Functions

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| PNG I/O library | libgd (C library) | No external library; data-only functions |
| Affected modules | `matrix/matrix_mod2_dense.pyx` | `matrix/matrix_mod2.ts` |

### Rationale

1. **No native image libraries** - JavaScript/TypeScript environments don't have a universal PNG library like libgd
2. **Environment portability** - Different environments (browser, Node.js, Deno) have different PNG handling approaches
3. **Separation of concerns** - Providing raw pixel data allows users to choose their preferred PNG encoding library

### Implemented Alternative

Instead of `from_png(filename)` and `to_png(A, filename)`, we provide:

- `from_png_data(width, height, pixels)` - Convert grayscale pixel array to GF(2) matrix
- `to_png_data(A)` - Convert GF(2) matrix to `{width, height, pixels}` object

The original `from_png()` and `to_png()` functions are retained for API compatibility but throw `NotImplementedError` with guidance to use the data-based alternatives.

### Trade-offs

- No direct file I/O for PNG images
- Users must handle PNG encoding/decoding themselves
- Requires additional library (e.g., `pngjs` for Node.js, Canvas API for browsers)

### Mitigation

For users who need file I/O:

```typescript
// Node.js example with pngjs
import { PNG } from 'pngjs';
import fs from 'fs';
import { from_png_data, to_png_data } from 'sagemath-ts/matrix';

// Reading PNG
const png = PNG.sync.read(fs.readFileSync('input.png'));
const grayscale = new Uint8Array(png.width * png.height);
for (let i = 0; i < grayscale.length; i++) {
  grayscale[i] = png.data[i * 4]; // Use red channel as grayscale
}
const matrix = from_png_data(png.width, png.height, grayscale);

// Writing PNG
const { width, height, pixels } = to_png_data(matrix);
const outPng = new PNG({ width, height });
for (let i = 0; i < pixels.length; i++) {
  outPng.data[i * 4] = outPng.data[i * 4 + 1] = outPng.data[i * 4 + 2] = pixels[i];
  outPng.data[i * 4 + 3] = 255; // alpha
}
fs.writeFileSync('output.png', PNG.sync.write(outPng));
```

### Behavioral Impact

- `from_png()` and `to_png()` throw `NotImplementedError` instead of performing file I/O
- `from_png_data()` and `to_png_data()` provide equivalent functionality with raw data
- Pixel convention matches SageMath: black (0) = 1, white (255) = 0

---

## Language & API Adaptations

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Callable parents | `F(3)` | `F.__call__(3n)` (factory helpers like `GF()` preferred) |
| Operator overloading | `a + b`, `a * b`, `a ** n` | `a.add(b)`, `a.mul(b)`, `a.pow(n)` |
| Keyword arguments | Python kwargs | Positional params or options objects |
| Generic ring elements | Duck typing | Explicit TypeScript interfaces |
| Affected modules | All | All |

### Rationale

1. **Language constraints** - TypeScript lacks Python's `__call__` and operator overloading
2. **Type safety** - Explicit interfaces/methods and options objects make types predictable
3. **Interoperability** - Using native `bigint` and explicit calls integrates with JS tooling

### Trade-offs

- Call sites are more verbose than SageMath
- Some dynamic behaviors (implicit coercion, kwargs) are not available

### Behavioral Impact

APIs are structurally equivalent but require explicit method calls and options objects.

---

## Return Type Differences in Arithmetic

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `gcd()` | Returns `Integer` | Returns `bigint` |
| `factor()` | `Factorization` object | `Array<[bigint, bigint]>` |
| `xgcd()` | Python tuple | Typed tuple `[bigint, bigint, bigint]` |
| `is_prime_power(..., get_data=True)` | `(p, k)` or `(n, 0)` | `[boolean, p, k]` |
| Affected modules | `arith/misc.ts`, `rings/integer_ring.ts` | Same |

### Rationale

1. **TypeScript idioms** - Use native `bigint` and tuple types
2. **Serialization** - Simple arrays are easier to serialize and test
3. **Type clarity** - Single return shape avoids overloaded return types

### Trade-offs

- Callers expecting SageMath object methods must adapt
- Some convenience methods (e.g., Factorization formatting) are helpers instead

### Behavioral Impact

Mathematical outputs are identical; return types differ.

---

## Ring/Field Iteration and Coercion

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Iteration over finite fields | Yields raw values | Yields element objects |
| Coercion | Implicit in arithmetic | Explicit in methods (accepts raw values) |
| `ZZ` | Module-level constant | Singleton instance with `__call__` |
| Affected modules | `rings/` | `rings/` |

### Rationale

1. **JS iteration model** - Uses `Symbol.iterator`
2. **Type safety** - Elements remain typed objects
3. **Singleton rings** - Matches SageMath semantics with a single `ZZ`

### Trade-offs

- Callers must access `.value` to get raw `bigint`
- Explicit coercion is sometimes required

### Behavioral Impact

Iteration yields element objects rather than raw values.

---

## Polynomial Representation Differences

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Coefficient storage | Context-dependent | `coeffs[i]` is coefficient of `x^i` |
| Trailing zeros | Preserved in some contexts | Trimmed |
| String conversion | `str(p)` | `p.toString()` |
| Affected modules | `rings/polynomial/*` | Same |

### Rationale

1. **Indexing clarity** - Direct degree-to-index mapping
2. **Memory** - Trim trailing zeros by default
3. **JS conventions** - `toString()` is standard

### Trade-offs

- Some SageMath display/ordering expectations differ
- Exact preservation of zero coefficients may differ

### Behavioral Impact

Coefficient arrays are normalized and stored in ascending degree order.

---

## Polynomial Roots and Factorization Limited

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `roots()` ring support | Many exact/approx rings (ZZ, QQ, finite fields, RR/CC, p-adics, number fields, symbolic) | Only base ring finite fields, ZZ, QQ |
| `roots()` options | `ring`, `multiplicities`, `algorithm`, keyword options | No ring override, always returns multiplicities |
| `factor()` ring support | Broad via FLINT/NTL/Singular/PARI | Finite fields, ZZ, QQ only |
| `factor()` return type | A `Factorization` object holding the unit separately (`F.unit()`); iterating yields only the non-unit factors | `Array<[Polynomial, number]>`; when the unit is not 1 it is included as an extra degree-0 factor, sorted first |
| Root ordering | `Polynomial.roots()` over GF(p) returns roots in factorisation order (PARI's `FpX_roots` is sorted, FLINT's is not) | `_poly_roots` and `rootsOverPrimeField` sort ascending by the integer representative |
| Cantor-Zassenhaus exhaustion | `AssertionError(f'no splitting of degree {degree} found for {self}')` | `Error` with the same message text (`errors.ts` exports no `AssertionError`); previously the code silently returned an unsplit factor |
| Affected modules | `sage/rings/polynomial/polynomial_element.pyx` | `packages/sagemath-ts/src/rings/polynomial/polynomial_element.ts` |

### Rationale

1. **Dependency gaps** - We have not yet ported the full PARI/FLINT/Singular stack used for polynomial root and factor algorithms.
2. **Incremental porting** - The current implementation targets common cryptographic cases (finite fields, ZZ, QQ).
3. **No `Factorization` type** - Returning the unit as a degree-0 factor restores
   `prod(factors) === f`, which several call sites and tests rely on.
4. **Deterministic representatives** - The small-field root finder already enumerated ascending
   while the large-field path used factorisation order; sorting makes the choice of
   representative deterministic and independent of the field size. This is what makes
   `montgomery_model` return Curve25519's canonical `A = 486662`.

### Trade-offs

- Roots and factorization over many rings are unavailable.
- Algorithm selection and ring overrides are not exposed.
- Callers must skip the degree-0 entry to iterate "the factors" the way Sage does.
- Root order can still differ from Sage's choice of representative for `isomorphism_to` /
  `montgomery_model` when several equally valid answers exist.

### Behavioral Impact

- `roots()` and `factor()` throw `NotImplementedError` for unsupported rings.
- Users cannot request alternative algorithms or ring coercions.
- **Known gap (not an intentional deviation):** factoring over QQ currently fails because
  `polynomial_element.ts:_factorOverRationals` calls
  `this.parent.base_ring.__call__({ numer, denom })` and `rational_field.ts:QQ.__call__` has no
  object form. This makes `minpoly`/`is_semisimple`/`is_similar` over QQ work only when the
  characteristic polynomial is squarefree. It is a real bug awaiting a fix, recorded here so it
  is not mistaken for a deviation.

---

## Integer Polynomial Factorization Simplified

> **Rewritten for 0.0.12.** The 0.0.11 text ("simplified modular factoring + limited rational
> root search", "`factor()` may return incorrect results") no longer describes the code and is
> retained below only as a record of what was fixed.

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Integer / rational polynomial factorization | FLINT `fmpz_poly_factor`: Zassenhaus, with van Hoeij/LLL recombination for hard inputs | A real **Zassenhaus**: prime selection (fewest modular factors among 3 usable primes), mod-`p` factorization by distinct-degree + Cantor-Zassenhaus equal-degree splitting (trace map at `p = 2`), multifactor linear Hensel lifting to `p^k`, an **exact-integer** Landau-Mignotte bound (integer `sqrt`, no floating point), and symmetric-range subset recombination with a `g(0) | lc*f(0)` filter |
| Recombination for hard inputs | van Hoeij's algorithm: LLL on the knapsack lattice of the lifted factors | Exhaustive subset search with a 200 000-subset budget; on exhaustion, `NotImplementedError` naming van Hoeij/LLL — **never** a partial factorization |
| Final check | — | The product of the returned factors is verified against the input; a mismatch raises `ArithmeticError` rather than returning a wrong factorization |
| Affected modules | `sage/rings/polynomial/polynomial_element.pyx` | `packages/sagemath-ts/src/rings/polynomial/polynomial_element.ts` |

### Rationale

1. **The old algorithm was returning wrong answers.** It only peeled off *integer* roots for
   degree <= 10 and handed the remainder back as a single "irreducible" factor, so
   `(x^2+1)(x^2+2)` and `6x^2+x-2` both came back unfactored. This was the blocker behind the
   `minpoly`-over-QQ gap (audit C8) and had to be fixed, not documented around.
2. **van Hoeij is a separate piece of machinery.** It needs an LLL over the knapsack lattice of
   the lifted factors. Bounding the classical search and failing honestly is strictly better
   than returning a partial answer.
3. **No floating point** - the Landau-Mignotte bound uses an integer square root, per CLAUDE.md.

### Trade-offs

- Recombination is exponential in the number of modular factors. Nothing constructible reached
  the budget: the classic worst case (Swinnerton-Dyer of degree 32, 16 quadratic modular
  factors) finishes in 324 ms, and `x^60 - 1` in 25 ms
- `chooseFactorizationPrime` scans primes below 10 000. For a polynomial whose leading
  coefficient and discriminant are divisible by every prime in that range it returns null and
  factorization throws rather than guessing. Raising the bound is trivial; the honest failure was
  left in place deliberately
- ~450 lines of superseded code (`factorByRationalRoots`, `berlekampFactor`,
  `splitUsingNullSpace`, `modMatrixNullSpace`, `distinctDegreeFactor`, `equalDegreeFactor`,
  `henselLift`, `extendedGcdPoly`, `modPowX`) were deleted

### Behavioral Impact

Factorization over ZZ and QQ is now **correct**, verified against an independently written
Kronecker-method oracle: the whole factorization (not just the product) matches on 200 random
QQ polynomials of degree <= 6, 400 composite QQ/ZZ products, and 150 products of
oracle-certified irreducibles up to degree 14 — 0 failures. `x^n - 1` has exactly one factor per
divisor of `n` for `n <= 60`; Swinnerton-Dyer polynomials of degree 8, 16 and 32 stay
irreducible.

Three further defects surfaced and were fixed in the same pass:

- **`is_irreducible()` over QQ counted the port's degree-0 unit factor**, so `2*x` and
  `2*x^2+2` were reported *reducible*. Sage's `Factorization` keeps the unit apart
  (`len(F) > 1 or F[0][1] > 1`), so degree-0 factors are now filtered over a field. Over ZZ
  nothing changed (the primitivity test already handled `2*x`)
- **Squarefree decomposition had a hard `maxIter = 20`**, so `(x-1)^25` was decomposed wrongly
  and then reached Zassenhaus with a non-squarefree input. The bound is now the degree, and the
  two "inexact division" escape hatches raise instead of silently truncating
- `is_irreducible()` over finite fields was already corrected in 0.0.11 to Rabin's test (the old
  `x^(p^n) = x mod f` check passed every fully-split polynomial)

---

## Groebner Basis Algorithms Simplified

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Gröbner basis engine | Singular/FGb with optimized F4/F5-style algorithms | Naive Buchberger implementation |
| Termination | Delegates to Singular; always terminates, no iteration budget | Keeps a `maxIterations` budget (default 10 000) and raises `ArithmeticError` when the S-pair queue is not drained, instead of returning the partial set |
| Base rings | `MPolynomialIdeal.groebner_basis` works over ZZ and Zmod(n) via Singular's `std` (doctests at `multi_polynomial_ideal.py:4593-4610`) | `groebner_basis()` raises Sage's `TypeError('Can only reduce polynomials over fields.')` for non-fields, matching `reduce()`; an additional exact-division check inside the loop guarantees termination even for base rings that expose no `is_field()` |
| Dimension | Exact Krull dimension via Singular | Ported Cox-Little-O'Shea algorithm (`multi_polynomial_ideal.py:1128-1192`); Sage's five `dimension()` doctests reproduce |
| Reduction strategy | Full-featured normal form over multiple coefficient rings | Basic multivariate division assuming coefficient fields |
| Affected modules | `sage/rings/polynomial/multi_polynomial_ideal.py` | `packages/sagemath-ts/src/rings/polynomial/multi_polynomial_ideal.ts` |

### Rationale

1. **Backend gap** - Singular/FGb backends are not available in TypeScript.
2. **Scope control** - A lightweight Buchberger implementation covers basic use cases.
3. **A truncated set is not a Gröbner basis** - `contains()`, `reduce()` and `dimension()` would
   return silently wrong answers; with a cap of 3, an ideal member demonstrably failed to reduce
   to 0. Throwing is the only honest option short of removing the budget entirely.
4. **Division over ZZ truncates** - Without the field guard, the quotient coefficient rounds to
   zero, the subtrahend is zero, and the reduction loop never terminates (reproduced as a
   timeout). An immediate error is strictly better than a hang.

### Trade-offs

- Much slower for large systems or high-degree inputs.
- Lacks advanced criteria and optimized reductions.
- `groebner_basis()` over ZZ/Zmod(n), which SageMath supports, now raises instead of hanging —
  it never worked here.
- Very large systems may hit the iteration budget and raise where Singular would finish.

### Behavioral Impact

- `groebner_basis()` and `reduce()` can be significantly slower, and raise rather than mislead
  on non-field coefficients or budget exhaustion.
- `dimension()` is now exact where it previously used a leading-term heuristic that could be
  wrong for nontrivial ideals.

---

## Generic Group API and DLP Limitations

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `discrete_log` options | Supports `bounds`, `algorithm` (`bsgs`, `rho`, `lambda`), `verify`, and `ord=oo` | Requires finite `ord`; implements exactly Sage's `bounds=None, algorithm='bsgs', verify=True` path (including the repair of an `ord` that is a proper multiple of the base's order, and the `<30` linear branch of `bsgs`) |
| `discrete_log_lambda` for `N = 1` | `k = 0`, then `hash % k` raises `ZeroDivisionError` | `k` is forced to at least 1, so a width-0/1 interval still works |
| `discrete_log_rho` options | `ord` optional; configurable `hash_function` | `ord` required; no custom hash function |
| Hashing in DLP algorithms | Uses element hashing/equality | `bsgs` and `discrete_log_rho` hash via string representations |
| `order_from_multiple` | `(P, m, plist=None, factorization=None, check=True, operation='+', …)` | `(P, m, factorization?, operation='+', identity?, inverse?, op?, options?: {plist?, check?})` — `check` now defaults to **true** and `plist` is honoured, but they are passed in a trailing options object |
| `has_order` | Accepts a Factorization | Accepts only integers |
| Generic utilities | `linear_relation`, `merge_points`, `structure_description` implemented | Not implemented |
| Affected modules | `sage/groups/generic.py` | `packages/sagemath-ts/src/groups/generic.ts` |

### Rationale

1. **Incremental porting** - The generic group API surface has been narrowed to core cryptographic use cases.
2. **Missing backends** - GAP-backed utilities (structure description) are not available.
3. **JS runtime limits** - There is no standard hash for custom objects, so string keys are used internally.
4. **Positional-argument compatibility** - `plist`/`check` are passed in a trailing options
   object rather than in Sage's positions, because `rings/finite_rings/integer_mod.ts` and other
   modules call `order_from_multiple` positionally as `(a, m, factorization, operation)`;
   inserting the new parameters positionally would have broken them.

### Trade-offs

- Fewer algorithm choices and bounds handling for discrete logs.
- Potential hash collisions for group elements with non-unique `toString()` outputs.
- Argument *positions* for `order_from_multiple` differ from Sage's even though the semantics
  now match.
- Some SageMath utilities are unavailable.

### Behavioral Impact

- Calls using `bounds`, `algorithm`, or `verify` on `discrete_log()` are unsupported.
- `discrete_log_rho()` requires an explicit prime order and cannot accept a custom hash.
- `order_from_multiple()` now honours Sage's `check=True` default. Prior to 0.0.11 it silently
  dropped the check, so `order_from_multiple(Mod(2,7), 5, '*')` returned 5 for an element of
  order 3; it now raises as Sage does.
- `has_order()` still ignores Factorization-style inputs.
- `linear_relation`, `merge_points`, and `structure_description` are missing.
- `discrete_log` in a 2^30-order subgroup runs in milliseconds rather than building a
  2^15-entry table, because the Pohlig-Hellman loop now follows Sage's verbatim.

---

## Finite Field Constructors and Display

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `GF(q)` extension fields | `GF` and `FiniteField` are the same factory object; `GF(p^n)` works in one call | `FiniteField` is aliased to `GFExtended` and handles both; the module-local `GF` remains prime-field-only, returning the narrowly typed `FiniteFieldPrime` |
| Extension field constructor | `GF(p^n, 'a')` | `FiniteField(p^n)` / `GFExtended()` / `GFpn()` |
| Element display | `repr(a)` and `a.lift()` | `a.toString()`, `a.repr()`, `a.value` |
| Affected modules | `rings/finite_rings/*` | Same |

### Rationale

1. **Package-level agreement** - Before 0.0.11 `GF(9)` worked at the package level while
   `FiniteField(9)` threw, even though SageMath makes them the same object. Aliasing
   `FiniteField = GFExtended` fixes that
2. **Type inference** - Retyping the module-local `GF` to a union would break the many
   elliptic-curve and matrix call sites that depend on the concrete `FiniteFieldPrime`
3. **Interop** - `.value` provides direct `bigint` access

### Trade-offs

- The two names are still not interchangeable inside the module: `GF` narrows to prime fields,
  `FiniteField` does not
- Users who want a statically typed prime field must call `GF`

### Mitigation

Retype the module-local `GF` to the union once the concrete-type call sites are migrated, and
then export a single factory under both names.

### Behavioral Impact

`FiniteField(p^n)` now constructs extension fields (this entry previously claimed "GF() only
supports prime fields", which under-described the surface). `GF()` still returns
`FiniteFieldPrime` and throws for prime powers.

---

## Conway Polynomial Database Limited

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Conway polynomial availability | The `conway_polynomials` package: every `(p, n)` in Frank Luebeck's `CPimport.txt` | A subset decoded from FLINT's vendored copy of the same data: p = 2 to n = 64, p = 3 to 24, p = 5 to 18, p = 7 to 14, p = 11/13 to 12, p = 17/19/23/29/31 to 10 |
| Default modulus for GF(p^n) | `PolynomialRing_dense_mod_p.irreducible_element` (`polynomial_ring.py:3560-3626`), falling through to `PolynomialRing_dense_finite_field.irreducible_element` (`:2628-2681`): `n == 1` -> `x - 1`; Conway if available; NTL `GF2X_BuildSparseIrred` for p = 2; else PARI `ffinit` (Adleman-Lenstra) | **Identical as of 0.0.12** — the same branch chain, delegating to `ntl-ts` `GF2X_BuildSparseIrred` and `parigp-ts` `ffinit` |
| Explicit `modulus=` algorithm strings | `'conway'`, `'adleman-lenstra'`, `'minimal_weight'`, `'first_lexicographic'`, `'primitive'`, `'random'`, `'ffprimroot'` (`finite_field_constructor.py:729-734`) | All except `'ffprimroot'`, with SageMath's exact error messages. `'primitive'` resolves to `'conway'` when a Conway polynomial exists (as in Sage) and otherwise throws |
| Affected modules | `sage/rings/finite_rings/conway_polynomials.py` | `packages/sagemath-ts/src/rings/finite_rings/conway_polynomials.ts`, `finite_field_extension.ts` |

### Rationale

1. **Correctness of the table** - The database is now regenerated by porting FLINT's
   `conway.c::conway_polynomial_lt_260` decoder against the vendored bit-packed table, rather
   than transcribed by hand. Prior to 0.0.11 **seven entries were reducible** (so e.g.
   `GF(29^2)` was not a field) and one — a `GF(2^128)` pentanomial — was fabricated outright
   (no Conway polynomial of that degree exists) and has been deleted. Every stored entry is now
   verified irreducible, primitive, correctly normalised, and subfield-compatible
   (`C_{p,m}(x^((p^n-1)/(p^m-1))) == 0 mod C_{p,n}`)
2. **Delegation targets now exist (0.0.12)** - `parigp-ts` gained a full `ffinit` (PARI's
   Adleman-Lenstra, `polarit3.c`) and `ntl-ts` already had `GF2X_BuildSparseIrred` /
   `GF2X_BuildIrred` / `GF2X_IterIrredTest`. `getDefaultModulus` was replaced by a faithful port
   of `irreducible_element` that delegates to both, so the 0.0.11 `first_lexicographic` fallback
   is gone. `polynomial_gf2x.ts` likewise now delegates its whole arithmetic and irreducibility
   layer to `ntl-ts` instead of reimplementing NTL
3. **The remaining gap is only the table size** - not the algorithm

### Trade-offs

- The table is finite, so `has_conway_polynomial` returns false where SageMath's does not, and
  outside the tabulated range the defining polynomial is `ffinit`'s / `BuildSparseIrred`'s
  rather than the Conway polynomial SageMath would use. SageMath's `conway_polynomials` package
  covers far more than ours — verified with Sage's own `exists_conway_polynomial`, which is
  `True` for `37^2`, `97^2`, `787^5`, `2^64`, `2^65`, `2^100`, `19^21`, `3^24`, `3^25`, none of
  which we have. **Consequence:** element representations of `GF(p^n)` for `p > 31` or large `n`
  are not interoperable with SageMath's, and the generator need not be primitive (Conway
  polynomials are; `ffinit`'s and NTL's are not)
- `algorithm='ffprimroot'` throws: Sage's line is
  `self(pari(p).ffinit(n).ffgen().ffprimroot().charpoly())`, and `parigp-ts` has `ffinit` but no
  finite-field element type (`ffgen`), no `ffprimroot` and no `charpoly` over `F_q`. It throws
  rather than returning some other irreducible polynomial, which would silently not be primitive
- `algorithm='random'` / `buildRandomIrred` over GF(2) uses rejection sampling from
  `current_randstate()` rather than NTL's `BuildRandomIrred` (which needs `IrredPolyMod` and
  `GF2XModulus`, absent from `ntl-ts`). This is exactly the fallback SageMath itself takes when
  its NTL import fails (`polynomial_ring.py:3615-3620`). The distribution differs from NTL's and
  the concrete polynomial for a given seed differs from Sage's
- `polynomial_gf2x.ts` keeps four local factoring routines (`squareFreeDecomp`,
  `distinctDegreeFactorization`, `equalDegreeFactorization`, `factor`) because `ntl-ts`'s
  `SquareFreeDecomp` / `DistinctDegFactor` / `EqualDegFactor` / `factor` / `BerlekampFactor` all
  still throw `NTL_NOT_IMPLEMENTED`. Delegating them would turn working functions into stubs

### Mitigation

Decode more of the FLINT Conway table to widen the range; add `ffgen`/`ffprimroot`/`charpoly`
to `parigp-ts`; add `IrredPolyMod`/`GF2XModulus`/`BuildRandomIrred` and the factoring routines
to `ntl-ts`.

### Behavioral Impact

Within the database the polynomial is the real Conway polynomial; outside it, the default is now
**SageMath's own** `ffinit` / `BuildSparseIrred` choice, so it agrees with Sage wherever our
Conway coverage matches Sage's. Verified against a SageMath 10.3 oracle: 45/45 explicit-algorithm
moduli match exactly, and of 130 sampled default moduli there are **0 genuine mismatches** — the
74 divergences were each confirmed, one by one with Sage's `exists_conway_polynomial`, to be
cases where Sage has a Conway entry our table lacks. `GF(1009^8)`, where Sage *also* has no
Conway entry, matches Sage byte for byte. Independently, an from-scratch Rabin oracle
(self-validated against Gauss's exact irreducible counts) confirms every `GF(2^n)` default for
`n` in [2,64] and all 420 defaults for the first 60 primes × `n` in [2,8] are irreducible.

Two smaller behaviour changes: a degree-1 extension's modulus is now `x - 1` rather than `x`
(Sage's `irreducible_element(1)` returns `self((-1,1))`; `GF(5)['x'].irreducible_element(1)` is
`x + 4`), so `gen()` is 1 as in Sage; and the `modulus` constructor parameter was widened to
accept Sage's algorithm strings alongside the existing `Polynomial | number[]` forms.

> **Note for future maintainers:** `GF(2^8).modulus()` is `x^8 + x^4 + x^3 + x^2 + 1`, the
> **Conway** polynomial (`finite_field_givaro.py:69`), *not* NTL's `BuildSparseIrred(8)` value
> `x^8 + x^4 + x^3 + x + 1` (the Rijndael polynomial). Conway wins Sage's branch order. The NTL
> value is reachable via `algorithm='minimal_weight'` and is verified there. Do not "fix" this.

---

## Finite Field Extension Minimal Polynomial Simplified

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `minimal_polynomial` over GF(p) | Computes full minimal polynomial via conjugates | Simplified: uses only constant coefficient of Frobenius conjugates |
| Affected modules | `sage/rings/finite_rings/*` | `packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.ts` |

### Rationale

1. **Implementation gap** - Full minimal polynomial requires polynomial arithmetic over extensions
2. **Incremental correctness** - Provide a placeholder for simple cases

### Trade-offs

- Incorrect for general elements (not just base-field elements)

### Mitigation

Implement minimal polynomial over extension fields using polynomial arithmetic with conjugates.

### Behavioral Impact

`minimalPolynomial()` can return incorrect results for nontrivial extension elements.

---

## Error Class Parity

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Error classes | Python built-ins | Custom classes (`ValueError`, `TypeError`, etc.) |
| Affected modules | All | All |

### Rationale

1. **Familiarity** - Mirrors SageMath exception names
2. **Consistency** - All extend JS `Error`

### Trade-offs

- Error class identity differs from Python's built-ins

### Behavioral Impact

Error names match SageMath, but they are JS subclasses of `Error`.

---

## Caching and Import Paths

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Caching | Extensive (`@cached_method`, `UniqueRepresentation`) | Minimal (zero/one caching, selective caches) |
| Imports | Deep module paths | Package root re-exports plus deep paths |
| Affected modules | All | All |

### Rationale

1. **Simplicity** - Avoid heavy caching until profiling indicates need
2. **Ergonomics** - Root re-exports simplify usage

### Trade-offs

- Possible performance differences for repeated computations
- Tree-shaking users must use deep imports explicitly

### Behavioral Impact

APIs are available but caching semantics differ from SageMath.

---

## Elliptic Curve Short Weierstrass Form Only

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Weierstrass models | General form supported in all characteristics | Short Weierstrass only; general form converted for char > 3 |
| Affected modules | `schemes/elliptic_curves/*` | Same |

### Rationale

1. **Implementation scope** - Short Weierstrass is the common case
2. **Simplification** - Avoids char 2/3 edge cases initially

### Trade-offs

- General form in char 2/3 not supported
- Some model transformations are implicit

### Behavioral Impact

Curves in characteristics 2 and 3 may throw or require explicit handling.

---

## Real/Complex Numerical Approximations

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Real precision | MPFR arbitrary precision | IEEE 754 double (`number`) |
| Complex results from real ops | Promote to complex field | Return `NaN` |
| Special functions | MPFR/Arb exact rounding | Polynomial/series approximations |
| Affected modules | `rings/real_mpfr.ts`, `rings/complex_mpfr.ts` | Same |

### Rationale

1. **Runtime limits** - JS lacks arbitrary-precision floats
2. **Portability** - Pure TypeScript implementations
3. **Pragmatism** - Accept approximations for non-integer special functions

### Trade-offs

- ~53-bit precision limit for reals
- No guaranteed correct rounding
- Some complex-valued outputs are not represented in real fields

### Behavioral Impact

Numerical results are approximate and may differ from SageMath at high precision.

---

## Matrix and Lattice Algorithm Simplifications

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| LLL | fpLLL (default) or NTL, floating-point Gram-Schmidt with precision escalation | Cohen 2.6.7 **exact integral LLL** in bigint (Gram determinants `d_i`, `lambda_ij = d_j * mu_ij`; every division exact). Default `delta` is now Sage's 0.99 (was 0.75); `eta` is accepted but advisory, since the exact algorithm always achieves `|mu| <= 1/2` |
| LLL on dependent rows | fpLLL/NTL run MLLL and return the zero rows first | Dependence is detected exactly, the generating set is replaced by the nonzero rows of its Hermite normal form (same lattice), and `nrows - rank` zero rows are prepended |
| BKZ | Full enumeration/pruning | Repeated LLL passes (approximate) |
| Frobenius form | PARI `matfrobenius`; `flag=2` returns `[F, B]` over QQ with `M = B^-1 F B` | **All three flags** are now a verbatim port of PARI's `RgM_Frobenius` (`alglin2.c:428-720`: Storjohann's Lemmas 9.14/9.18, Ozello's theorem 4, `_frobTransL/D/S`, `_minpoly_polslice/listpolslice/dvdslice`) over exact rational arithmetic, 1-indexed to match `gcoeff(M,i,j)` line for line. `flag=2` returns `[F, B]` as `Rational[][]` |
| `right_kernel_matrix` over Z/nZ | PARI `matkermod` for composite `n` | Prime modulus is a faithful port of `matrix_modn_dense_template.pxi:2072` (all three `basis` formats); **composite `n` now delegates to `parigp-ts` `matkermod`**, mirroring `matrix_modn_dense_template.pxi:2136`'s fallback to `matrix2.pyx:4233` |
| `matrix_modn` determinant | LinBox for prime p > 2; generic path otherwise | `n <= 3` uses Sage's naive formulas; `n >= 4` always uses centered-lift-to-ZZ + fraction-free Bareiss |
| HNF transformation matrix | `fmpz_mat_hnf_transform` | Classical row operations; `U` differs from FLINT's for rank-deficient input (both satisfy `U*A == H`, and `H` matches Sage exactly) |
| `p_minimal_polynomials`, `null_ideal`, `integer_valued_polynomials_generators` | `compute_J_ideal` algorithms | Simplified outputs (charpoly/minpoly-based) |
| Affected modules | `matrix/matrix_integer.ts`, `matrix/matrix_modn.ts` | Same |

### Rationale

1. **Exactness over reproduction** - CLAUDE.md forbids floating point where Sage is exact. The
   previous double-precision LLL stopped producing a basis of the input lattice above 2^53 and
   its `is_LLL_reduced` checker certified bases that were not reduced
2. **Dependency gap closed in 0.0.12** - `parigp-ts` gained `matkermod` (a full port of
   `bb_hnf.c` specialised to the `Z/dZ` black-box Hermite ring: Howell form, `gen_kernel`,
   `gen_matimage`, `gen_inv`, `gen_detops`), so `right_kernel_matrix` delegates instead of
   raising. `matfrobenius` was ported directly into `matrix_integer.ts` (see
   [PARI/NTL Routines Ported In Place](#parintl-routines-ported-in-place-instead-of-delegated))
3. **`flag=2` returns rationals, not an `IntegerMatrix`** - `B` is genuinely rational (e.g.
   `[[1,-2,1],[0,-23/15,14/15],[0,2/15,-1/15]]` for `range(9)`), so the declared-but-unreachable
   `[IntegerMatrix, IntegerMatrix]` signature was replaced by `[Rational[][], Rational[][]]`.
   The port has no rational-matrix class. Nothing could depend on the old type — the code threw

### Trade-offs

- An exact LLL cannot reproduce fpLLL's rounding-dependent choice of representative, so
  individual rows can differ, typically by sign: `matrix(ZZ,3,range(1,10)).LLL()` row 1 is
  `[2,1,0]` in Sage and `[-2,-1,0]` here; `matrix(ZZ,[[1,2,3],[31,41,51],[101,201,301]]).LLL()`
  row 1 is `[-1,0,1]` in Sage and `[1,0,-1]` here. The result is always a
  `(delta, 1/2)`-reduced basis of the same lattice
- For rank-deficient input the reduced rows may differ from MLLL's, because the starting basis
  differs
- The generic determinant path is slower than LinBox for prime moduli (values are identical)
- **`frobenius_form` flags 0 and 1 changed block ORDER in 0.0.12.** The previous code derived the
  invariant factors from the Smith normal form of `xI - A` and emitted them smallest-first;
  PARI emits the minimal polynomial **first**, with `P_{i+1} | P_i`. For `diag(1,1,2)` Sage gives
  `[x^2-3x+2, x-1]` and `F = [[0,-2,0],[1,3,0],[0,0,1]]`; the port now matches. This changes the
  meaning of an already-exported symbol, and one existing test that pinned the reversed order was
  corrected to SageMath's values
- `frobenius_form` with an out-of-range flag now raises `ValueError('incorrect flag in
  matfrobenius')` — PARI's own message text with the port's error class, since there is no
  `PariError` in `sagemath-ts` — where flag 3 previously fell through to the flag-2 stub
- `matkermod` is called with `wantIm = true` even though the image is discarded; see
  [Upstream Behaviour Deliberately Not Reproduced](#upstream-behaviour-deliberately-not-reproduced)
  for why (it disables an unsound PARI shortcut). Extra Howell work on tall matrices

### Mitigation

Move `matfrobenius` from `matrix_integer.ts` into `parigp-ts` and delegate, as
`right_kernel_matrix` now does for `matkermod`.

### Behavioral Impact

Values are now exact **and complete**: `frobenius_form(2)` and composite-modulus
`right_kernel_matrix` answer where 0.0.11 raised. Both were validated against SageMath 10.3 as an
oracle, not merely against invariants — 300 random integer matrices agree with Sage on `F`, the
elementary divisors and **both halves** of `flag=2`, with `B^-1 F B == A` verified exactly on all
300; 300 random composite-modulus matrices match Sage's kernel basis entry for entry, and 2190
brute-force cases confirm the returned rows generate the *full* enumerated kernel (size and
membership), not merely a subset.

Prior to 0.0.11, `smith_form_integer` could loop forever, `matrix_modn` `charpoly` had a
Faddeev-LeVerrier sign error, and `right_kernel_matrix` echelonized the transpose.

---

## Shortest Vector Problem (SVP) Implementation

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `IntegerLattice.shortest_vector` | fpylll `SVP.shortest_vector` or PARI `qfminim` | Exact Fincke-Pohst enumeration in bigint/rational arithmetic |
| Affected modules | `modules/free_module_integer.ts` | Same |

### Rationale

1. **No external dependencies** - Pure TypeScript implementation without fpylll or PARI bindings
2. **Exactness** - The previous implementation was a floating-point Schnorr-Euchner walk with a
   hard coefficient cap of 15, which was neither exact nor safe for large entries
3. **Portability** - Works in any JavaScript runtime without WebAssembly or native bindings

### Implementation Details

- **Rank <= `EXACT_SVP_MAX_RANK` (30)**: exact Fincke-Pohst enumeration
- **Rank > 30**: returns the first row of the LLL-reduced basis (2^((n-1)/2) approximation)

### Trade-offs

- Performance is slower than fpylll's optimized C++ implementation
- No pruning strategies from advanced BKZ variants
- Large lattices (rank > 30) use approximation instead of exact SVP
- `isLLLReduced` is likewise exact now (no `1e-10` fudge, default `delta` 0.99); it skips the
  leading zero rows LLL emits for rank-deficient input, and raises Sage's `'sage'`-algorithm
  `ValueError('linearly dependent input for module version of Gram-Schmidt')` for genuinely
  dependent nonzero rows

### Mitigation

For performance-critical applications with large lattices, users can:
1. Use LLL approximation directly via `lattice.LLL()` first row
2. Use BKZ for better approximation quality
3. Interface with external SVP solvers if needed

### Behavioral Impact

- For rank <= 30, results are exact shortest vectors (same as SageMath)
- For rank > 30, returns an approximation (first LLL basis vector)
- SageMath with algorithm='pari' uses qfminim which may enumerate differently
- The only remaining floating point in this module is the legacy exported `gramSchmidt()`
  helper (used by `bkz.ts` and `discreteGaussianSample`) and the heuristic estimators
  (`hadamardRatio`, `gaussianHeuristic`, `hermiteFactor`, `estimateBKZBlockSize`)

---

## Closest Vector and Voronoi Cell

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `IntegerLattice.closest_vector` | Projects `t` onto the span, then Micciancio-Voulgaris over the diamond-cut Voronoi cell | Exact Fincke-Pohst enumeration seeded with Babai's nearest plane (the projection is unnecessary: the orthogonal component of `t` is a constant offset of the objective) |
| `approximate_closest_vector` | `nearest_plane` / `rounding_off` / `embedding` | Same, with Sage's `embedding` default and round-half-to-even |
| `voronoi_relevant_vectors` | Reads the defining point of each inequality of the diamond-cut `Polyhedron` | Enumerates the `2^r - 1` nonzero cosets of `L/2L` and keeps those with exactly two minimal-length vectors (Voronoi's theorem); results sorted lexicographically |
| `voronoi_cell` | Returns a `Polyhedron` | Returns an H-representation `{normals, offsets}` (gcd-normalised `normals . x <= offsets`, `offsets` exact `bigint`) |
| Affected modules | `modules/free_module_integer.ts` | Same |

### Rationale

1. **Exactness** - Both CVP and SVP are now exact. The previous `closestVector` enumerated
   around the **origin** with coefficients in `[-3, 3]`, so for a distant target it degraded to
   Babai (a rank-3 example gave `d^2 = 125` against the true 98), while the registry claimed it
   was exact for rank <= 4
2. **No Polyhedron class** - The port has no polyhedral geometry, so the facet description is
   the natural representation of the same object; the relevant-vector set is mathematically
   identical and was verified sound and complete against brute force
3. **Cost** - Enumeration is far cheaper than building the Voronoi cell, which is exponential
   in the rank

### Trade-offs

- When several lattice vectors are equidistant from the target, which one is returned may
  differ from SageMath's
- `voronoiCell` rank > 24 raises `NotImplementedError` rather than running out of memory
- Callers wanting a polyhedron object must build it from the inequalities themselves

### Behavioral Impact

Results are exact closest vectors, matching SageMath's value (though not necessarily its choice
among ties). All five of Sage's `approximate_closest_vector` doctest values and the
`voronoi_relevant_vectors` doctest reproduce.

---

## Algebraic Dependency Approximation

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `algebraic_dependency` | PARI `algdep`, supports real/complex/p-adic with proof/height bounds; `complex_mpfr.pyx:3252` calls `sage.arith.misc.algebraic_dependency` for both real and complex inputs | `arith/misc.ts` takes a `number` (real only); the **complex** branch is ported separately inside `complex_mpfr.ts` |
| Precision | `prec = z.prec() - 6`, uncapped | Complex branch uses `prec = min(z.prec(), 53) - 6` |
| Irreducible-factor selection | `min over irreducible factors of |f(z)|` | Performed in the complex branch; `arith/misc.ts` still omits it |
| Affected modules | `arith/misc.ts` | `arith/misc.ts`, `rings/complex_mpfr.ts` |

### Rationale

1. **No PARI binding** - `algdep` not available yet
2. **Type mismatch** - `arith/misc.ts:algebraic_dependency` takes a `number`, so a
   `ComplexNumber` cannot be passed to it; the complex branch (the `is_complex` column layout,
   the constant-first-row fallback, LLL from `matrix_integer`, irreducible-factor selection via
   `PolynomialRing` over ZZ) is ported line-for-line into the caller instead
3. **Precision cap** - the scaling factor `2^prec` multiplies a double, so beyond 53 bits the
   extra bits fed into the lattice are pure rounding noise and can make LLL return a spurious
   relation. Capping keeps the answer correct for any declared parent precision

### Trade-offs

- The same algorithm exists in two places; a follow-up should move the whole thing into
  `arith/misc.ts` and have both modules delegate
- No p-adic support
- Proof/height-bound semantics are not equivalent to SageMath
- `ComplexField(200)` gains no accuracy over `ComplexField(53)`

### Behavioral Impact

Complex inputs now work and reproduce SageMath's doctest (`x^2 - x + 1` with
`p(z) = 1.11e-16`) at both `ComplexField(200)` and `ComplexField(30)`. Real results can still be
approximate or incorrect for subtle inputs; the API accepts options but does not fully match
SageMath's guarantees.

---

## Combinatorial Function Limits

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `number_of_partitions(n)` | FLINT implementation handles large `n` | DP up to `n <= 10000`, then throws |
| `prime_pi(n)` | primecount/FLINT for large `n` | Naive counting up to `n <= 10,000,000`, then throws |
| Affected modules | `rings/integer_ring.ts` | Same |

### Rationale

1. **Algorithm complexity** - Efficient implementations are not yet ported
2. **Reasonable limits** - Provide correctness for small inputs

### Trade-offs

- Large inputs raise `NotImplementedError` instead of computing

### Behavioral Impact

Calls above thresholds fail where SageMath would compute exact results.

---

## Real Matrix Decompositions (SVD, QR, LU) Using IEEE 754 Doubles

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Numeric precision | `RDF` uses doubles, `QQ`/`ZZ` use exact arithmetic, scipy.linalg.svd wraps LAPACK | Pure TypeScript with IEEE 754 doubles only |
| SVD algorithm | scipy.linalg.svd (LAPACK DGESDD/DGESVD) | Jacobi one-sided SVD |
| QR algorithm | scipy.linalg.qr (LAPACK DGEQRF) | Householder reflections |
| Affected modules | `matrix/matrix_double_dense.pyx`, `matrix/matrix2.pyx` | `matrix/matrix_decompositions_additions.ts` |

### Rationale

1. **No LAPACK bindings** - TypeScript lacks direct scipy/LAPACK integration
2. **Pure implementation** - Self-contained algorithms without external dependencies
3. **Target use cases** - Numerical linear algebra for moderate-sized problems fits IEEE 754

### Implementation Details

- **SVD_double()**: Jacobi one-sided algorithm with configurable tolerance and max iterations
  - Computes U, S, V such that A = U * diag(S) * V^T
  - Singular values sorted in descending order
  - Returns full U (m x m) and V (n x n) matrices

- **QR_double()**: Householder reflections (more stable than Gram-Schmidt)
  - Supports both full and reduced QR
  - Returns orthogonal Q and upper triangular R

- **LU_double()**: Gaussian elimination with partial pivoting
  - Returns permutation vector P, lower triangular L (unit diagonal), upper triangular U
  - PA = LU factorization

### Trade-offs

- ~53-bit precision limit (IEEE 754 double)
- No arbitrary precision or exact rational matrix decompositions
- Jacobi SVD converges slower than divide-and-conquer LAPACK algorithms for large matrices
- No support for complex matrices (real only)

### Mitigation

For users needing higher precision:
1. Use the generic matrix decompositions (e.g., `LU`, `QR` from `matrix_decompositions.ts`) with exact ring elements
2. For very large matrices, consider interfacing with external LAPACK via WebAssembly

### Behavioral Impact

- Results are numerically equivalent to SageMath's `RDF` matrix operations within IEEE 754 precision
- Exact rational/integer matrix operations use separate functions in `matrix_decompositions.ts`
- Default tolerance is 1e-14, matching typical double precision expectations

---

## PARI Factorization Algorithms Limited (parigp-ts)

> **Rewritten for 0.0.12**, **completed in 0.0.14.** `ifactor.ts` is PARI's real factoring
> chain, and as of 0.0.14 MPQS (`mpqs.ts`) fills the last slot, so every stage of `ifac_crack`
> is present. The failure mode also changed from a silently wrong answer to an exception.

| Aspect | SageMath (PARI/GP) | sagemath-ts (parigp-ts) |
|--------|--------------------|-------------------------|
| `Z_factor` chain | `ifac_crack` (`ifactor1.c:2786`): trial division -> pure powers -> SQUFOF -> Pollard-Brent rho -> ECM (non-insisting) -> **MPQS** -> ECM (insisting), driven by the `ifac_decomp` worklist | The **same chain, same order**, ported with file:line citations: `tridiv_bound` + gcd-with-primorial fast trial division (`Z_oddprimedivisors_fast`, `:3306`), SQUFOF (`:1474`, incl. `squfof_ambig`), Pollard-Brent (`:1184`/`:1361`, incl. Brent fast-forward, backtracking and multi-factor returns), Lenstra-Montgomery ECM (`ellfacteur`/`ECM_loop`, `:752`/`:1038`) |
| MPQS (`mpqs.c`, ~2600 lines) | Self-initializing multiple-polynomial quadratic sieve with a relation database, GF(2) Gaussian elimination and disk-backed large-prime relations | **Ported** in `packages/parigp-ts/src/mpqs.ts` and wired into `ifac_crack` in PARI's position (behind non-insisting ECM, ahead of insisting ECM). The relation store is in memory rather than disk-backed, and `MpqsOptions.maxPolys` can bound the sieve (0 = unbounded, as PARI). Like PARI, it declines above 107 decimal digits (`mpqs.h:400`); only then does `Z_factor` throw |
| Perfect powers / `isprimepower` | `ispower.c`: `Z_issquareall`, `is_357_power`, `is_kth_power`, `is_pth_power`, `Z_isanypower_101`, `Z_isanypower`, `isprimepower` | All ported, using **exact integer k-th roots** (Newton) instead of PARI's `sqrtnr`/`mpexp` float guesses, per CLAUDE.md. `isprimepower` never factors `n` |
| Affected modules | `pari/src/basemath/ifactor1.c`, `ispower.c` | `packages/parigp-ts/src/ifactor.ts` |

### Rationale

1. **The old behaviour was a wrong answer, not a slow one.** `Z_factor` used to `console.warn`
   and return an unfactored composite **as if it were prime**. An exception is strictly better.
2. **The chain is now complete.** MPQS landed in 0.0.14, so `Z_factor` splits the hard
   semiprimes that used to fall off the end of the chain.
3. **No floating point** - PARI's root guesses are heuristics followed by an exact check, so
   replacing the guess with an exact integer root cannot change the verified answer.

### Trade-offs and intentional divergences from PARI

- **MPQS's relation store is in memory.** PARI spools large-prime relations to disk
  (`mpqs.c`, `pari_unique_filename`); we keep them in a `Map`. Same relations, same combined
  full relations, no behavioural difference — only a memory ceiling on very large inputs
- **The five call sites that used to raise now succeed.** `parigp-ts/src/elliptic/group.ts:463`,
  `parigp-ts/src/elliptic/advanced.ts:119`, `sagemath-ts/src/arith/misc.ts:949`,
  `rings/padics/padic_generic.ts:57` and `rings/padics/padic_generic_element.ts:135` inherit the
  completed chain and no longer throw on hard semiprimes
- **Pollard-Brent's `n < 2^96` size gate is still removed** (`ifactor1.c:1361-1367`). PARI
  declines there because MPQS is faster in that range. The gate was dropped while MPQS was
  absent and has deliberately **not** been restored: restoring it would only change which stage
  of `ifac_crack` splits a composite, never the factorization returned, and it would invalidate
  the existing `pollardbrent` range tests. Rho is a complete algorithm, so only the time changes
- **Pollard-Brent's round budget clamped at the low end.** PARI's `c0` formula
  (`ifactor1.c:1369-1372`) goes negative below ~60 bits, which PARI never evaluates because it
  declines first. Clamped to at least `tune` (14)
- **SQUFOF accepted up to 2^59 rather than the 64-bit build's 2^46** (`ifactor1.c:1487`). The
  2^46 cut-off is a tuning choice that hands 2^46..2^59 to MPQS; the algorithm and PARI's own
  comment (`:1492`) document validity for `5 < n < 2^59`. Kept after MPQS landed, for the same
  reason as the rho gate above: it changes which stage splits `n`, never the answer. Measured 0.5 ms per 56-bit semiprime,
  200/200 split. All arithmetic stays exact
- **ECM is serial with a binary ladder and an additive stage 2.** PARI runs up to 64 curves in
  parallel with Montgomery batched inversion, multiplies with Montgomery's PRAC chain, and runs
  stage 2 over a 48-entry helix of residue classes mod 210 with a 1024-entry baby-step table.
  We keep the same curve family, `nbc`/`dsn`/seed derivation, `TB1`/`TB1_for_stage` schedules and
  B1 phase, but process curves one at a time and step `[p]Q` additively. Detection is equivalent
  (`[p]Q` vanishes mod a prime divisor exactly when the denominator does); these are
  constant-factor optimizations. PARI's *insisting* mode loops forever — ours is bounded by
  `FactorOptions.ecmRounds` (default 4). MPQS now backs it up as in PARI, so the bound only
  affects how much work the insisting stage does before MPQS is reached again
- **`is_357_power`'s residue sieve is omitted** (mod 211/209/61/203, then 117/31/43/71 through a
  106-entry mask table). A pure speed filter; the exact k-th root is taken directly, keeping
  PARI's mask semantics. Verified against brute force exhaustively on 2..20000
- **`pollardbrent_i` returns null where PARI raises `pari_err_BUG`** (`:1330`) after 4 unlucky
  restarts. A can't-happen path; degrading to "this engine failed" cannot produce a wrong
  factorization
- **File placement:** `isprimepower`, `Z_isanypower`, `Z_isanypower_101`, `Z_issquareall`,
  `is_kth_power`, `is_pth_power`, `is_357_power` and `Z_iroot` live in `ifactor.ts`, but PARI
  keeps them in `basemath/ispower.c`. They should move to `parigp-ts/src/ispower.ts`; each
  carries a JSDoc naming its upstream file and line

### Behavioral Impact

`Z_factor` is now **correct and complete** everywhere it answers. Verified: agrees with
brute-force trial division exhaustively on 1..5000 and on 4000 random `n < 1e9`; a
2000-random-semiprime sweep where the product is restored, every factor is BPSW-prime and the
factorization equals exactly `{p,q}`; and the published factorizations of `F6 = 2^64+1`, `M67`,
`M71`, `M101` and `F7 = 2^128+1` (the last found by the ported ECM in 16 s).

`isprimepower` no longer factors its argument: for 24-digit primes `p`, `q`,
`isprimepower(p*q)` is `null` and `isprimepower(p*p)` is `[p,2]`. Consequently
`IntegerMod.log()` now **does** run Sage's
`(Mod(5, 123337052926643^4)^(10^50-1)).log(5)` doctest, returning
`99999999999999999999999999999999999999999999999999` in 67 s — the 0.0.11 claim that it cannot
is stale.

`Z_factor` and `factoru` gained an optional `options?: FactorOptions` (`{ ecmRounds }`) second
parameter. Existing one-argument call sites are unaffected.

---

## PARI Elliptic Curve Advanced Algorithms (parigp-ts)

> **Updated for 0.0.14.** SEA is now ported (`packages/parigp-ts/src/elliptic/ellsea.ts`,
> a port of `ellsea.c`), together with the modular polynomials it needs
> (`packages/parigp-ts/src/polmodular.ts`, a port of `polmodular.c`/`polclass.c`/`volcano.c`),
> so the `seadata` package is no longer required: `Phi_L` is computed on demand and cached in a
> `polmodular_db_*` database. The rows below that describe SEA as absent are superseded.

| Aspect | SageMath (PARI/GP) | sagemath-ts (parigp-ts) |
|--------|--------------------|-------------------------|
| `ellcard` dispatch | naive trace enumeration for `expi(p) < 11`, `Fp_ellcard_CM`, `Fp_ellcard_Shanks` in the middle range, SEA for `expi(p) >= 56` (`FpE.c:1424-1437`) | Same naive branch below p = 2048, then `Fp_ellcard_CM`, then Shanks, then **base Schoof** from `expi(p) >= 96` (a *measured* threshold — see below) |
| SEA (Schoof-Elkies-Atkin) | `ellsea.c`, needs the `seadata` modular-polynomial package | **Ported in full** as `Fp_ellcard_SEA` (`elliptic/ellsea.ts`): Elkies, Atkin and the match-and-sort final step, plus `Fp_elljissupersingular` and the CM branch. `seadata` is replaced by `polmodular.ts`, which computes `Phi_L` on demand. The older `ellcard_sea` in `advanced.ts` still delegates to base Schoof and is retained for its regression tests |
| `Fp_ellcard_CM` | Full CM table (`Fp_ellj_get_CM` + `ec_ap_cm`: `ap_j0`, `ap_j1728`, `ap_j8000`, `ap_j287496`, …) | **All thirteen** class-number-one discriminants, ported line by line from `FpE.c:624-666` and `:1282-1421`, delegating to `qfb.ts`'s `cornacchia2`. Includes PARI's signed-int `(CM&3)==0 -> CM>>=2` semantics and the `case -28: ap_cm(-7, -114, …)` quirk |
| `gen_ellgroup` `m` output | `bb_group.c:1035-1043` writes `*pm = g1` and then overwrites it with `*pm = m` (the `lcm(s,t)` of the *final* iteration) | Returns `m = g1` |
| `Fp_ellcard_Shanks` visibility | `static` in `FpE.c` | `export`ed, so the test suite can exercise the BSGS branch on primes small enough for an exhaustive point-count oracle |
| `random_FpE` | `FpE.c:369-385` returns `Fp_sqrt(rhs, p)`, the canonical smallest root | Same (the previous `Math.random() < 0.5` sign flip is removed) |
| `j` / `ellj` return type | t_INT or t_FRAC | `bigint` when `c4^3` is divisible by the discriminant, else an exact `Ratio {num, den}` (with an exported `isRatio` guard) |
| Advanced functions (`ellisogeny*`, `ellfrobenius`) | Fully implemented | Stubs that throw `PARI_NOT_IMPLEMENTED` |
| Affected modules | `pari/src/basemath/ellsea.c`, `FpE.c`, `bb_group.c`, `ellisog.c` | `packages/parigp-ts/src/elliptic/{group,points,init,advanced}.ts` |

### Rationale

1. **SEA needed the modular polynomials, so they were computed rather than read.** Elkies and
   Atkin both need `Phi_l(X,Y) mod p`, which PARI reads from `$pari_datadir/seadata/sea<l>`
   (`ellsea.c:47-51` `seadata_filename`, `:51-101` `get_seadata`, `:139` `pari_err_FILE`).
   `reference/pari` ships the **reader** but not the **data** (`reference/pari/data` is empty).
   0.0.14 resolved this by porting `polmodular.c`/`polclass.c`/`volcano.c` instead
   (`src/polmodular.ts`), which is how PARI *generates* `seadata` in the first place, and
   caching the result in a `polmodular_db_*` database. `find_trace`,
   `find_trace_Elkies_powerell`, `find_trace_Atkin`, `match_and_sort` and `champion` are now
   ported (`src/elliptic/ellsea.ts`)
2. **Schoof's base algorithm is still present** (it was what was reachable before seadata), and
   implemented in full: division polynomials `psi~_m`, arithmetic in
   `F_p[x]/(psi_l)[y]/(y^2-f)` exploiting the `(X(x), b(x)·y)` shape, automatic restart modulo a
   proper divisor of `psi_l` whenever a denominator turns out to be a zero divisor (valid because
   the Frobenius relation holds identically on `E[l]`), plus PARI's own trace-mod-2/mod-4
   seeding, its `TR_mod > sqrti(16p)` stopping bound and its `smallfact` early abort, ported
   verbatim from `ellsea.c:2000-2059`
3. **The CM gap closed once `cornacchia2` existed** (`qfb.ts`, new in 0.0.12)
4. **Termination over literal fidelity for `m`** - When the primes of `N0` are not all settled
   in a single iteration of `gen_ellgroup` (which is common), the final iteration's `m` need
   not be a multiple of `d2`, and then `gen_ellgens` can never terminate: the Weil pairing of
   two `m`-torsion points has order dividing `m`. Measured on `E/F_43: y^2 = x^3 + 7x + 8`
   (group `[12,3]`, `N0 = 2^2 * 3^2`), about 0.5 % of runs produce `m = 4` with `d2 = 3` —
   4 hangs in 885 runs. `g1` provably satisfies `d2 | g1 | d1`, so `gen_ellgens` always
   terminates, and PARI 2.15.4 (the version shipped with Sage 10.3) never hangs on that curve
   over 4000 fresh `ellgenerators` calls, so `g1` reproduces the *shipping* PARI behaviour; the
   vendored 2.18-dev line appears to be a regression
5. **Exactness** - Returning `Ratio` is how the j-invariant becomes exact; there is no rational
   type in this package

### Trade-offs

- **Base Schoof is `O(log^5 p)` with schoolbook `FpX` arithmetic** where SEA is `O(log^4 p)`
  with a far smaller constant. Measured on this port, single random curve, Schoof vs Shanks:
  56 bits 12.8 s / 0.10 s; 64 bits 21.1 s / 0.39 s; 72 bits 82.6 s / 2.41 s; 80 bits 101.7 s /
  4.85 s; 88 bits 189.6 s / 26.1 s; 96 bits ~358 s / 296 s (9 GB rss). **So `ellcard` keeps
  Shanks below `expi(p) = 96`, not PARI's 56** — copying PARI's threshold would make `ellcard`
  hundreds of times slower over 2^56..2^88, and 96 is where the two cross on time and where
  Shanks' `p^(1/4)` baby-step table stops being memory-feasible. The value returned is
  unaffected; the two algorithms were cross-checked against each other
- `Fp_ellcard_Schoof` does **not** take PARI's `j = 0` / `j = 1728` shortcut
  (`ellsea.c:1990-1993`), because routing back into `ellcard` would be a recursion hazard and
  would remove those curves from the Schoof test oracle. `ellcard` applies the CM shortcut before
  ever reaching Schoof, so the fast path is preserved where it matters
- **`cornacchia2` failure inside the `ap_*` helpers falls back instead of returning 0.** PARI
  writes `(void)cornacchia2(...)` and ignores the return value, leaving the out-parameter as
  `gen_0` — which would report a wrong trace. Ours return `null`, `Fp_ellcard_CM` returns `null`,
  and `ellcard` falls through to Shanks/Schoof. Unreachable in theory (every `d` is a
  class-number-one discriminant and the Kronecker symbol is pre-checked) and never fired in
  ~20 000 tested CM curves, but a plausible wrong cardinality is the one outcome to avoid
- `m = g1` means the pairing is computed at a possibly larger exponent than PARI's, i.e.
  marginally slower. No output difference: group structure and generator validity were verified
  over 10 068 curves
- `Ratio` / `isRatio` are not yet re-exported from `packages/parigp-ts/src/index.ts`, and neither
  are the new `Fp_ellcard_CM`, `Fp_ellj_get_CM`, `Fp_ellj_nodiv`, `ec_ap_cm`,
  `Fp_ellcard_Schoof` or `Fp_elldivpol`
- `Fp_elldivpol(l, a4, a6, p)` is a **new public function with no direct PARI counterpart over
  `F_p`** (PARI's SEA never builds `psi_l` because it uses modular polynomials). Exported so the
  recursion is testable in isolation
- `random_FpE` returns canonical points; `<P>` and `<-P>` are the same subgroup, so order,
  group-structure and pairing consumers are unaffected

### Remaining gap

`ellcard` in `elliptic/group.ts` still dispatches to `Fp_ellcard_Schoof` above
`SCHOOF_BIT_THRESHOLD = 96` bits (`group.ts:1318`, `:1357`). It should call `Fp_ellcard_SEA`,
which is what PARI's `Fp_ellcard` does and which is orders of magnitude faster there
(0.1 s at 101 bits versus minutes for base Schoof). `Fp_ellcard_SEA` is exported from the
package root and can be called directly in the meantime.

### Behavioral Impact

`ellcard_sea` **no longer throws** — it returns exact cardinalities at every size. Verified
exhaustively against brute force on all 121 104 curves over every prime `5 <= p <= 120`, against
Shanks on 375 random curves from 20 to 88 bits, and — most importantly — against **PARI's own
`ellsea` regression vectors** (`reference/pari/src/test/in/ellsea` entries `v[9]`, `v[10]`,
`v[11]`) at 65, 70 and 101 bits, all three exact. The zero-divisor split path was instrumented
and fired 5984 times over 16 308 curves with zero errors.

`Fp_ellcard_CM` covers all thirteen discriminants, verified against brute-force point counting
(15 392 curves), against Shanks (3744 curves to 48-bit primes), against a counting-independent
`[#E]P = O` oracle (936 curves at 64/80/96 bits), and against the **published SECG group orders**
of secp160k1, secp192k1, secp224k1 and secp256k1.

Prior to 0.0.11 `ellcard` returned wrong cardinalities at primes as small as p ≈ 100–1069 and
`ellgroup` was wrong in 85 of 476 cases; both were fixed then and remain verified. Calls to
`ellisogeny`, `ellisogenyapply`, `ellisogenycompose` and `ellfrobenius` still throw.

---

## Discrete Gaussian Samplers

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `DiscreteGaussianDistributionIntegerSampler` algorithms | `uniform+table`, `uniform+online`, `uniform+logtable`, `sigma2+logtable` | **All four**, as of 0.0.12. `dgs_bern.c` (`dgs_bern_uniform_t` with its 32-bit bit pool, `dgs_bern_mp_t`, `dgs_bern_exp_mp_t`) and `dgs_gauss_mp.c`'s `dgs_disc_gauss_sigma2p_t` are ported, including the `sigma <- k*sigma_2` adjustment with `sigma_2 = sqrt(1/(2 ln 2))` and `MPFR_RNDN` rounding of `k` |
| Acceptance test | `mpfr_urandomb` compared against the tabulated probability | Same. 0.0.11 used `c_rand_double() < rho/rhoMax`, which is *not* dgs's test, so even the two implemented algorithms consumed randomness differently from Sage |
| `repr()` | `Discrete Gaussian sampler over the Integers with sigma = %f and c = %f` (six decimals, reporting the **adjusted** sigma for `sigma2+logtable`) | Identical. 0.0.11 printed an invented `DiscreteGaussianDistributionIntegerSampler(sigma=…, c=…, tau=…)` |
| `DiscreteGaussianDistributionLatticeSampler.__call__` | Returns a vector over the base ring of the basis matrix — integers for `ZZ^n`, rationals for a `QQ` basis | `sample()`/`samples()` keep the `bigint[]` signature and throw `ValueError('lattice basis is not integral; use sampleExact()')` for a non-integral basis; `sampleExact()`/`samplesExact()` return exact `Rational[]` for any basis. A public `isIntegral` flag and exact `basisExact`/`cExact` fields were added |
| `DiscreteGaussianDistributionPolynomialSampler` | Lives in `sage.crypto.lwe` with signature `(P, n, sigma)`, returning a polynomial in the ring `P`; there is **no** such class in `sage.stats.distributions.discrete_gaussian_lattice` | `crypto/lwe.ts` carries the faithful `(P, n, sigma)` class; `discrete_gaussian_lattice.ts` *additionally* exports a convenience `(n, options)` form returning a coefficient array |
| Non-spherical Σ | Matrix sigma, Peikert's `r`, Cholesky, offline samples, `_call_non_spherical` | **Implemented** in 0.0.12: covariance matrices (with the scaled-identity collapse and `sigma_basis` for `Σ = S Sᵀ`), `_maximal_r()` by power iteration on `Q Σ^-1`, `_precompute_data`'s Peikert branch (exact `B_inv` over QQ, `r = 0.9999·_maximal_r`, `B2 = chol(Σ - r²Q)ᵀ`), `add_offline_samples`, `_randomise`, `_call_non_spherical` |
| `set_c` / `c()` / `sigma()` / `f()` / `_normalisation_factor_zz` | Present | **All implemented** in 0.0.12, including the Poisson-summation normalisation factor (both the `sigma > 1` Fourier branch and the direct branch), the `tau` shortcut, all four `NotImplementedError` guards with Sage's messages, and the non-spherical LLL + `BOUND`-box branch |
| `_normalisation_factor_zz` theta series | PARI `Q.__pari__().qfrep(B, 0)` | A **local** `qfrep` exported from `discrete_gaussian_lattice.ts` (Fincke-Pohst enumeration over the Cholesky factor, exact integer counts). `parigp-ts` has no `qfrep` |
| `precision='dp'` constructor keyword | Routes through `dgs_gauss_dp.c` (drand48 / libc `random()`) | **Absent** — the option does not exist, and neither does Sage's `ValueError("Parameter precision '…' not supported")` |
| Affected modules | `sage/stats/distributions/discrete_gaussian_{integer,lattice}.py`, `dgs_gauss_mp.c`, `dgs_bern.c` | `packages/sagemath-ts/src/stats/distributions/` |

### Rationale

1. **The Bernoulli machinery is now ported** - `dgs_bern.c`'s `dgs_bern_exp_mp` and
   `dgs_gauss_mp.c`'s `dgs_disc_gauss_sigma2p` were the whole blocker; with them the four
   algorithms are only ~200 lines apart
2. **`dp` carries no oracle value** - Sage itself documents that "in the latter case results are
   not reproducible", so there is nothing to match against
3. **No polymorphic "element of the base ring"** - TypeScript cannot express it, and every
   downstream consumer (LWE, crypto) wants `bigint`. Splitting the accessor keeps the common
   integral case statically typed while making the general case exact, instead of the previous
   silent `BigInt(Math.round(...))` corruption
4. **Avoiding cross-module churn** - Moving or removing the duplicate polynomial sampler would
   break `crypto/index.ts` and `stats/index.ts`, so it is registered instead, with an explicit
   `@see Deviation:` docstring pointing callers at the Sage-faithful version in `crypto/lwe.ts`

### Trade-offs and remaining divergences

- **`_normalisation_factor_zz` is evaluated in double precision.** Sage evaluates in
  `RealField(prec)`, so `prec=100` yields ~30 correct digits. `prec` is accepted (as
  `compute_precision` does) but does not raise the working precision, because there is no
  MPFR-equivalent here. Sage's `round(DGL(ZZ^8, 1000)._normalisation_factor_zz(prec=100))`
  doctest prints the 28-digit `1558545456544038969634991553`; ours is `1.5585454565440373e27`,
  agreeing to 15 significant digits and diverging at the 16th — all a double can carry
- **`qfrep` belongs in `parigp-ts`**, which has none. Output verified identical to PARI's on six
  forms; this is an architecture-fidelity gap, not a correctness gap
- **LLL inside the non-spherical branch is a local exact-rational LLL** (delta 0.99) rather than
  fpLLL. It is used only to pick the enumeration centre; the two agree on the reduced lattice but
  may return different reduced bases on inputs with several near-optimal reductions
- **Non-positive-definite `Σ` raises the built-in `Error`**, not `RuntimeError`, with Sage's
  message text byte-identical including its column padding. (`RuntimeError` was added to
  `errors.ts` in this pass and the throw switched — see the changelog)
- **`sigma` and `c` changed from public properties to methods** (`D.sigma()`, `D.c()`), matching
  Sage. `cNumeric()` is kept as a float convenience view. Source-breaking, but the only consumers
  were this module's own tests
- **`BernExpMp` stores `p = 0` where dgs leaves memory uninitialised.** `dgs_bern.c:121-124`
  breaks out of the table build when `exp()` underflows but still sets `l = i+1`, so index `i`
  may be read while `p[i]`/`B[i]` were never initialised — undefined behaviour in C. We store an
  explicit zero, which deterministically rejects. Unreachable for any realistic sigma
- **`sigma2+logtable` rejects `sigma < sigma_2`.** dgs computes `k = round(sigma/sigma_2)` with
  no lower guard, so `k = 0` makes `mpz_urandomm` divide by zero. We throw `ValueError` naming
  `sigma_2 = 0.849322`. An added error, not a changed result
- Callers of a non-integral lattice must choose a different method name (`sampleExact`)
- A polynomial-sampler class exists in a module where SageMath has none

### Behavioral Impact

**All four integer-sampler algorithms now reproduce SageMath's *seeded sample streams* exactly**,
not merely its distributions: 14 pinned 16-sample streams across three parameter sets and the two
non-integer-centre paths, plus Sage's own `_flush_cache` doctests. Distribution correctness is
additionally confirmed by 200 000-sample chi-squared tests at the 0.999 level over 11 settings.

For the lattice sampler, the locally installed Sage 10.3 predates non-spherical support, so the
vendored source's doctests were the oracle; every one reproduces, including
`_maximal_r() = 0.584028653716433`, the `Σ` recovered from a `sigma_basis`, the exact
`RuntimeError` text, the two `f()` values, the three `__repr__` forms and the five spherical
`_normalisation_factor_zz` values. The empirical 3×3 covariance of 20 000 samples matches `Σ` to
< 0.6 in every entry.

Rounding of a non-integer centre follows dgs's round-half-to-even (`dgs_gauss_mp.c:161-165`), so
`sigma=3, c=1.5, tau=2` has support `[-4, 8]` (`c_z = 2`), not the `floor(c)` window.

One pre-existing test was statistically unsound and was corrected rather than deleted: "every
support point of a narrow sampler is hit" demanded that points with expected count 2.2 always
appear, which a *correct* sampler fails for ~20 % of seeds. Its threshold was raised to an
expected count of 20; it still covers 19 support points over a [-18,18] window.

---

## PARI/NTL Routines Ported In Place Instead of Delegated

CLAUDE.md requires that where SageMath delegates to an external library, we delegate to our port
of that library. In several places the target package has no such module and is owned elsewhere,
so the upstream algorithm was transcribed into the calling module instead. These are
**architectural** deviations: the values are PARI's, only the file location differs.

### Resolved in 0.0.12 — now genuinely delegated

| Routine | Ported into `parigp-ts` / `ntl-ts` as | Call site now delegating |
|---------|--------------------------------------|--------------------------|
| `ffinit` (Adleman-Lenstra, `polarit3.c`) + the supporting `FpX` layer, `FpX_composedsum`, `polsubcyclo` | `parigp-ts/src/ffinit.ts` | `finite_field_extension.ts` `irreducible_element` |
| `matkermod` / `matimagemod` / `matdetmod` / `matinvmod` (`bb_hnf.c`, Howell form over `Z/dZ`) | `parigp-ts/src/matkermod.ts` | `matrix_modn.ts` `right_kernel_matrix` for composite `n` |
| `qfbred`, `qfbredsl2`, `qfbcompraw`, `qfbcomp`, `qfbsqr`, `qfbpow`, `qfbsolve`, `primeform`, `cornacchia`/`cornacchia2`, `Zn_quad_roots` (`Qfb.c`, `quad.c`) | `parigp-ts/src/qfb.ts` | `quadratic_forms/binary_qf.ts` `reduced_form` / `compose` / `solve_integer`; `parigp-ts/src/elliptic/group.ts` `Fp_ellcard_CM` |
| `GF2X` arithmetic, `IterIrredTest`, `BuildIrred`, `BuildSparseIrred` | `ntl-ts` (already present in 0.0.11) | `polynomial_gf2x.ts` (whole arithmetic + irreducibility layer) and `finite_field_extension.ts` |

`binary_qf.ts` shed ~170 lines of transcribed `Qfb.c` in the process, and `polynomial_gf2x.ts`
lost its ad-hoc trinomial/pentanomial scan in favour of NTL's real precomputed minimal-weight
table.

### Still ported in place

| Routine | SageMath delegates to | We implement it in |
|---------|----------------------|--------------------|
| `isprimepower` / `Z_isanypower` | PARI `basemath/ispower.c`, via `Integer.is_prime_power` -> `__pari__().isprimepower()` | **Two copies now:** the real port lives in `parigp-ts/src/ifactor.ts` (see [PARI Factorization](#pari-factorization-algorithms-limited-parigp-ts)) but is not re-exported from that package's `index.ts`, so `sagemath-ts/src/arith/misc.ts:836` still keeps its local `isprimepower`/`anyPower101`. Deleting the local copy is a one-line follow-up once the barrel exports it |
| `matfrobenius` (`alglin2.c:428-720`) | PARI | `packages/sagemath-ts/src/matrix/matrix_integer.ts` — ported verbatim there rather than into `parigp-ts`, which still has no matrix module for it |
| `nfbasis`, `nfdisc`, `idealprimedec` (incl. round 4), `nfgaloisconj`, `quadunit`, `polisirreducible` | PARI `base2.c`, `galconj.c`, `quad.c`, `polarit2.c` | `packages/sagemath-ts/src/rings/number_field/pari_nf.ts` |
| `qfbclassno` / `quadclassunit` | PARI `Qfb.c` | `packages/sagemath-ts/src/rings/number_field/class_group.ts` (counts reduced primitive forms) |
| `qfrep` | PARI | `packages/sagemath-ts/src/stats/distributions/discrete_gaussian_lattice.ts` (Fincke-Pohst; output verified identical to PARI's) |
| `polinterpolate` | PARI | `lagrange_polynomial(points, 'pari')` computes the same unique interpolant with the same divided-difference scheme, locally |
| `dilog`, `incgam` | PARI | `packages/sagemath-ts/src/rings/complex_mpfr.ts` (`parigp-ts` has no transcendental functions at all) |
| NTL `SFCanZass` / `CanZass` factoring over GF(2) | NTL, via `polynomial_template.pxi` | `polynomial_gf2x.ts` keeps local `squareFreeDecomp` / `distinctDegFactorization` / `equalDegreeFactorization` / `factor`, because `ntl-ts`'s versions all throw `NTL_NOT_IMPLEMENTED` |

### Rationale

1. **The delegation target did not exist** - and for four routine families it now does, so those
   rows moved to the table above. What remains is where `parigp-ts` still has no `nf`, matrix or
   transcendental module, and `ntl-ts` no randomized factoring
2. **File ownership** - These fixes landed in parallel work units that did not own the
   dependency packages
3. **Verified equivalence** - Each is a transcription of the cited upstream source, checked
   against PARI/Sage values, not a reinvention

### Trade-offs

- The dependency graph still does not fully mirror SageMath's
- `class_group.ts` deliberately does **not** reuse `quadratic_forms/binary_qf.ts`, so the two
  share no code
- `parigp-ts` has no shared errors module, so `PariTypeError`, `PariDomainError`, `PariDimError`,
  `PariInvError`, `PariPrimeError`, `PariSqrtnError` and `PariFlagError` are defined once in
  `matkermod.ts` and imported by `ffinit.ts` and `qfb.ts`. They should move to
  `parigp-ts/src/errors.ts` when such a file exists, re-exporting the same names
- `matkermod`'s matrices use PARI's **column-major** `ZM` layout (`A[j][i]` is row `i`, column
  `j`), so "basis vectors are the columns" stays true. `zm_from_rows`/`zm_to_rows` converters and
  a row-major `matkermod_basis` wrapper exist so a Sage-side caller cannot silently transpose the
  meaning

### Mitigation

Add an `nf` module, matrix routines and transcendental functions to `parigp-ts`; re-export
`isprimepower` and friends from `parigp-ts/src/index.ts` and delete the `arith/misc.ts` copy; add
the randomized factoring routines to `ntl-ts`.

### Behavioral Impact

None on outputs. `Integer.is_prime_power` no longer needs to re-verify its answer locally, and
`IntegerMod.log()` **does** now run Sage's
`(Mod(5, 123337052926643^4)^(10^50-1)).log(5)` doctest — the 0.0.11 text saying it cannot was
stale and has been removed.

---

## Infinity Representation

SageMath has `sage.rings.infinity` with genuine `PlusInfinity` elements. This port has no
infinity ring, so three different sentinels are in use depending on the module's return type.

| Method | SageMath returns | sagemath-ts returns |
|--------|-----------------|---------------------|
| `Integer.multiplicative_order()` (for `n` other than ±1) | `infinity` | The string literal `'Infinity'` (typed `bigint \| 'Infinity'`) |
| `Rational.valuation()`, `val_unit()`, `RationalField.quadratic_defect()` | `infinity` | `'Infinity'` (pre-existing convention in that module) |
| `pAdicGenericElement.valuation()`, `multiplicative_order()`, `additive_order()` | `infinity` | `Number.POSITIVE_INFINITY`, typed `InfiniteOr<bigint> = bigint \| number` |
| `RealNumber`/`ComplexNumber` `multiplicative_order()`, `additive_order()` | `infinity` | `Number.POSITIVE_INFINITY` |
| `FreeModule.cardinality()`, `indexIn()` | `infinity` | `Number.POSITIVE_INFINITY` |
| `QuotientRing.cardinality()` | `Infinity` | `Number.POSITIVE_INFINITY` (return type widened to `bigint \| number`) |
| `arith.valuation(0, p)` | `+Infinity` | **Raises** `ValueError('valuation of 0 is infinite')` |

### Rationale

1. **`bigint` has no infinite value** - and widening every arithmetic return type to admit one
   would ripple through every caller
2. **Numeric comparisons work** - JavaScript relational comparisons between `bigint` and
   `number` behave as expected (`Infinity > 0n`), so guard conditions read naturally where the
   `number` sentinel is used
3. **Module-local consistency** - the string sentinel was already established in
   `rings/rational*.ts` before this pass; changing it there would be a larger breaking change
   than the inconsistency costs
4. **`arith.valuation` returns `bigint`** - introducing a sentinel would change the type for
   every caller, so it throws instead (pre-existing behaviour, now registered)

### Trade-offs

- Three sentinels for one concept; callers must know which module they are in
- `valuation(0, p)` fails where SageMath answers

### Mitigation

Introduce a single `Infinity` singleton (or a branded type) and migrate all seven sites to it.

### Behavioral Impact

Only the representation of the infinite case differs; all finite values are exactly Sage's.
`Integer(1).multiplicative_order()` is `1n` and `Integer(-1).multiplicative_order()` is `2n`,
as in Sage.

---

## Exact Arithmetic Where SageMath Uses Floating Point

CLAUDE.md forbids floating point where Sage or PARI is exact. In several places SageMath itself
uses a bounded-precision float as a *heuristic* guarded by an exact verification; we replace the
heuristic with exact integer arithmetic. Results are identical within the range where Sage's
float has enough precision, and stay correct beyond it.

| Site | SageMath | sagemath-ts |
|------|----------|-------------|
| `Z_isanypower_101` perfect-power search | `logr_abs`/`mpexp` double precision to guess `y = round(x^(1/p))`, with a mod-30011 filter, then an exact `powiu(y,p) == x` check | Every prime exponent `p <= log_103(x)` tested with an exact bigint Newton k-th root |
| `BinaryQF._reduce_indef`, `_Rho`, `_RhoTau`, `BinaryQF_reduced_representatives` | `D.sqrt(prec=53)` and a floored real quotient | Exact `isqrt(D)`. `floor((sqrt(D)+b)/(2|c|)) == floor((isqrt(D)+b)/(2|c|))` because the numerator bound is integral, and the `|c| >= sqrt(D)` branch boundary yields the same `s` either way |
| `Matrix_integer_dense.LLL`, `IntegerLattice` reduction | fpLLL/NTL floating-point Gram-Schmidt | Exact integral LLL (see [Matrix and Lattice Algorithm Simplifications](#matrix-and-lattice-algorithm-simplifications)) |
| Free module `coordinates`/`echelonize`/`discriminant`/kernels | Exact over the fraction field | Exact fraction-field layer (see [Free Module Exactness](#free-module-exactness-and-coordinate-types)) |
| `Integer.real_log` | MPFR | Exact above 2^53 (`ln(10^k) = k*ln10`), `Math.log` below |
| `complex algebraic_dependency` precision | `prec = z.prec() - 6` | Capped at 53 (see [Algebraic Dependency Approximation](#algebraic-dependency-approximation)) |

### Rationale

1. **CLAUDE.md rule** - "Don't use floating point - use BigInt and rational arithmetic"
2. **Sage's floats are heuristics** - each is followed by an exact check upstream, so replacing
   the guess with an exact computation cannot change the verified answer
3. **Correctness past 2^53** - the float paths fail silently on large inputs; the exact ones do
   not. This also fixes the class of bugs SageMath tracks in its issue 37635

### Trade-offs

- The exact search can be slower than the float heuristic (negligible at the sizes involved:
  perfect-power exponents are at most `log_103(n)`)
- Exact LLL cannot reproduce fpLLL's choice of representative (documented separately)

### Behavioral Impact

Identical `(base, exponent)` pairs, identical reduced forms, identical logarithms — with the
float paths' silent failures removed.

---

## Upstream Behaviour Deliberately Not Reproduced

Places where the vendored upstream is itself buggy or crashes, and we deviate on purpose.

| Site | SageMath / PARI | sagemath-ts | Why |
|------|-----------------|-------------|-----|
| PARI `matkermod`'s `m > 2n` shortcut | `bb_hnf.c:1049` computes the kernel as `shallowtrans(matimagemod(shallowtrans(A), d))`. For `A == 0` the image has zero columns, and a `t_MAT` with zero columns carries no row count, so the transpose collapses to 0×0 and the reported basis generates only `{0}` while the real kernel is everything. Reproduced on SageMath 10.3: `matrix(Zmod(4),3,1,[0,0,0]).right_kernel_matrix()` is `[]` | `_right_kernel_matrix_over_integer_mod_ring` calls `matkermod(A, n, wantIm=true)`. PARI's own condition is `!im && m > 2*n`, so requesting the image disables the unsound branch. `matrix(Zmod(4),3,1,[0,0,0])` gives `[[1]]` | Returning an empty basis for a non-trivial kernel is silently wrong and would make a caller conclude the kernel is trivial. A Sage sweep of 2400 small cases showed the two agree **everywhere except** the zero matrices with `m > 2n` (62 cases, all all-zero); below that threshold Sage itself returns the identity, so the fix restores consistency with PARI's own answer in the sound regime |
| PARI 2.18.1's Schoenhage fast reduction for binary quadratic forms | `qfi_red` (`Qfb.c:975-987`), `qfi_redsl2` (`:857-880`), `qfr_red_i` (`:914-940`) and `qfr_redsl2` (`:825-855`) negate `b` before `pqfbred_rec` (which requires non-negative coefficients) and never negate it back. For `b < 0` that returns the reduced form of the **conjugate** (inverse) class, and `qfi_redsl2` negates only the second *row* of `U`, giving `det U = -1` | The intermediate form is conjugated back and `D*U*D` (`D = diag(1,-1)`) is used instead of a row negation, so `det U == 1` and `Q ∘ U == result`. `qfr_redsl2` additionally falls back to the base case when `a < 0`, where upstream negates all three coefficients without touching `U` | The path is gated at 9000 bits of excess size, so it is unreachable for anything Sage currently does — but "never return a plausible wrong answer" outranks bug-compatibility. Verified by forcing the threshold to −∞ on 80 sample forms: definite results become byte-identical to the base case, indefinite results are genuinely reduced and SL2-equivalent |
| PARI's `matdetmod` documentation | `reference/pari/src/functions/linear_algebra/matdetmod` claims `matdetmod([4,2,3;4,5,6;7,8,9],27) == 9` | Returns `18` | The determinant is −9, i.e. 18 mod 27, and PARI's **own regression output** (`reference/pari/src/test/32/bbhnf` line 240) records 18. A stale doc block, not a behavioural difference. Sage's `qfbsolve` doc block is stale in the same way — it predates `allsols`' `v >= 0` normalisation and the lexsort at `Qfb.c:1930` — so we follow the source and assert the solution *sets* |
| `arith.CRT_basis(moduli, false)` | `arith/misc.py:3695-3725` appends to `cs` inside the coprime loop and, after catching the `ValueError`, does `cs.extend(...)` onto that same non-empty list, so for e.g. `[7,6,10]` it returns **more entries than there are moduli** | The partial entries are discarded, so the result always has exactly `moduli.length` coefficients | A basis longer than the modulus list makes `CRT_vectors` index out of range or silently mis-combine. The documented `[60,90,150]` case is unaffected (Sage bails on the first modulus there, so `cs` is empty) |
| `groups.generic.discrete_log_lambda` with `N = 1` | `k = 0`; the loop body never runs and `hash % k` raises `ZeroDivisionError` | `k` is forced to at least 1 | Avoids a crash where Sage crashes; results for all `N >= 2` are unaffected |
| PARI `gen_ellgroup`'s `m` output | `bb_group.c:1035-1043` overwrites `*pm = g1` with the final iteration's `lcm(s,t)`, after which `gen_ellgens` can fail to terminate | Returns `m = g1` | See [PARI Elliptic Curve Advanced Algorithms Missing](#pari-elliptic-curve-advanced-algorithms-missing-parigp-ts); PARI 2.15.4 (shipped with Sage 10.3) does not hang, so `g1` reproduces the shipping behaviour |
| `ell_finite_field.twists()` `break` placement | `ell_finite_field.py:1940-1944` puts the `break` at the for-loop level, so only `twists[0]` is ever tested for isomorphism with `self` | **Replicated verbatim**, with a comment citing the line numbers | CLAUDE.md requires behavioural equivalence with the vendored Sage. Implementing the docstring's stated intent would change the returned ordering for `j = 0`/`1728` curves. The returned *set* is complete and pairwise non-isomorphic either way |
| `elementary_matrix(row1 == row2)` with no scale | `special.py:1512-1516`'s four assignments collapse to `elem[r,r] = 1`, i.e. the identity; Sage raises only when a scale is *also* given | Same: identity for a self-swap, `ValueError` for the two cases Sage rejects | The audit stated Sage raises in general; the vendored source shows otherwise, so the port follows the actual code |

### Rationale

Reproducing an upstream crash or an out-of-range return value would propagate the defect into
every consumer, with no offsetting fidelity benefit — nothing can depend on a `ZeroDivisionError`
or on a basis of the wrong length. Where the upstream quirk is merely a *choice* (the `twists`
break placement, the `elementary_matrix` self-swap), we reproduce it exactly.

### Trade-offs

- Code written against SageMath's buggy `CRT_basis` output length would behave differently here
- Divergence from the vendored 2.18-dev PARI source in `gen_ellgroup`

### Behavioral Impact

`CRT_basis` returns `moduli.length` coefficients always (a regression test asserts the invariant
over four non-coprime modulus lists). `discrete_log_lambda` answers where Sage raises.
`ellgenerators` always terminates.

---

## Honest Failure Instead of Silent Approximation

A large share of the 0.0.11 pass replaced code that returned a plausible-looking wrong answer
with an explicit error. These are listed together because they share one rationale. **Rows marked
"resolved in 0.0.12" now compute the answer** — they are kept here as a record of what changed
and in which direction.

| Site | Previously | Now |
|------|-----------|-----|
| `matrix_operations.norm(A, 2)` | Returned the Frobenius norm; then (0.0.11) `NotImplementedError` naming the missing SVD | **Resolved in 0.0.12** for exact-rational entries: the largest eigenvalue of `A^H A` is isolated *exactly*. Still raises for inexact real/complex base rings — see [Matrix Module Algorithm Substitutions](#matrix-module-algorithm-substitutions) |
| `matrix_operations.is_diagonalizable(A, base_field)` | Silently ignored `base_field`; then `NotImplementedError` | **Resolved in 0.0.12**: a generic `change_ring` was implemented (a port of `matrix0.pyx:1666-1715`) with a conservative coercion policy |
| `matrix_operations.is_similar(A, B, transformation=true)` | Returned a fabricated 0×0 zero matrix; then `NotImplementedError` when similar | **Resolved in 0.0.12**: Sage's `jordan_form` route is wired up, backed by an intertwining-equation solver. `[false, null]` when not similar |
| `matrix_decompositions.jordan_form(transformation=true)` | `NotImplementedError` | **Resolved in 0.0.12**: a port of `matrix2.pyx:12259-12312` + `_jordan_form_vector_in_difference` (`:20895`), returning Sage's *exact* `P`, not merely a valid one |
| `matrix_decompositions.krylov_kernel_basis(var=…)` | `NotImplementedError` ("needs Popov/approximant bases") | **Resolved in 0.0.12** — and the premise was wrong: Sage does not use approximant bases here, it derives the kernel directly from the Krylov basis as `relation = D·C^-1` (`matrix2.pyx:20343-20478`) |
| `eigenvalues`/`eigenvectors` | Defaulted to `extend=false`, silently returning a short list | Default `extend=true` as in Sage; raises `NotImplementedError` naming the missing algebraic closure when the charpoly does not split over the base ring |
| `matrix_integer.frobenius_form(2)` | Returned the identity as the change of basis; then `NotImplementedError` | **Resolved in 0.0.12**: PARI's `RgM_Frobenius` ported in full |
| `matrix_modn.right_kernel_matrix` (composite modulus) | Returned a wrong kernel; then `NotImplementedError` naming PARI `matkermod` | **Resolved in 0.0.12**: delegates to `parigp-ts` `matkermod` |
| `polynomial factor()` over ZZ/QQ | Returned an under-factored list (a composite reported as one irreducible factor) | **Resolved in 0.0.12**: a real Zassenhaus. `NotImplementedError` only when subset recombination exhausts its budget (van Hoeij/LLL) |
| `padic_lseries.bernardi_sigma_function` / supersingular `alpha` | Hardcoded coefficients wrong from `z^5` up; then `NotImplementedError` | **Resolved in 0.0.12** |
| `PARI Z_factor` on a hard composite | `console.warn` and returned the composite **as a prime factor** | Factored by the full `ifac_crack` chain including MPQS. Only inputs above MPQS's own 107-digit ceiling (`mpqs.h:400`) still raise |
| `random_echelonizable_matrix` / `random_unimodular_matrix` `upper_bound` | Accepted and ignored | `NotImplementedError` (the port's generic constructors have no notion of absolute value) |
| `MPolynomialRing.__call__` from a univariate polynomial | Fell into the plain-object dictionary branch and produced nonsense | `NotImplementedError` naming `_mpoly_dict_recursive` |
| Unknown multivariate term order | Silently fell back to degrevlex | `ValueError("unknown term order 'name'")`, Sage's own message |
| `groebner_basis` budget exhaustion / non-field base ring | Returned the partial set / hung | `ArithmeticError` / Sage's `TypeError` |
| `BCH minimum_distance()` | Returned the *designed* distance | Exact enumeration, cached; `NotImplementedError` once `q^k > 2^17`, pointing at `designed_distance()` |
| `BooleanFunction.truthTable('hex')` with `n < 2` | Returned an unpadded nibble | `ValueError('negative shift count')` — Sage's own failure |
| `pAdicLseries._c_bound` / `bernardi_sigma_function` | Returned a literal `0` / hardcoded coefficients wrong from `z^5` up | `NotImplementedError` naming the missing dependency |
| `compute_isogeny_stark` | — | `NotImplementedError` naming `E.weierstrass_p()` |
| `booleanHypercube(n)` | `1 << n` wrapped around at n = 31/32 | `ValueError` above `MAX_HYPERCUBE_DIM = 25` |
| `voronoiCell` rank > 24, exact SVP rank > 30, `nfgaloisconj` degree > 8, quadratic class group `|D| > 2e6`, Zassenhaus prime bound 2^200, BCH field embedding `|E| > 2^22`, Goppa/RS decode beyond the radius | Ran out of memory, or returned a wrong answer | Explicit `NotImplementedError` / documented approximation |

### Rationale

1. **A wrong answer is worse than no answer** in a library whose stated goal is exact
   behavioural equivalence, and whose consumers are cryptographic
2. **Grep-ability** - CLAUDE.md rule 7 requires unimplemented paths to carry
   `SAGE_NOT_IMPLEMENTED`, so the gaps are discoverable with
   `grep -r "SAGE_NOT_IMPLEMENTED" packages/`
3. **Bounded search where the algorithm is exponential** - an explicit cap documents the real
   reachable range instead of hanging

### Trade-offs

- Calls that previously "worked" now raise. In every case listed the previous answer was
  wrong, ignored an argument, or was a placeholder — but this is a **source-breaking change**
  for callers that were not checking
- Some inputs SageMath handles (ZZ Gröbner bases, large BCH minimum distances,
  `Frobenius_filter` over number fields, degree > 2 class numbers, hard-to-factor integers) still
  fail here
- Going the *other* way is also source-breaking: eight of the rows above changed from "throws" to
  "returns a value" in 0.0.12. Two existing tests that pinned the stub behaviour were rewritten
  to pin SageMath's real values (the 2-norm test in `matrix_operations.test.ts` and the
  negative-`n` `solve_integer` test in `binary_qf.test.ts`); no test was deleted or weakened

### Behavioral Impact

Errors where there used to be plausible garbage. Every message names the missing dependency or
the algorithm that would be needed. As of 0.0.12 the honest-failure set is materially smaller —
see the changelog for the full accounting of what became real and what is still deferred.

---

## Port-Only APIs With No SageMath Counterpart

Symbols that exist in this port and **not** in SageMath. They are listed so their presence is
never mistaken for a SageMath contract, and so `@see Deviation:` docstrings have a target.

| Symbol | Module | Notes |
|--------|--------|-------|
| `polynomial_commitment.ts` (866 lines: `compute_quotient`, `batch_quotient`, `barycentric_weights`, `fri_fold`, `split_poly`, `generate_powers`, …) | `zk/` — **moved out of `rings/polynomial/` in 0.0.12** | KZG/FRI helpers. There is no `sage/rings/polynomial/polynomial_commitment.py`, so living inside the mirrored Sage tree implied a provenance it does not have. `rings/polynomial/index.ts` retains a re-export block of the identical 38 values + 2 types, marked in-file as backwards compatibility only, so the `./rings` and `./rings/polynomial` package subpaths are unchanged. Object identity was verified across all five import surfaces. Note `package.json` has **no `./zk` subpath export** yet, so a direct `@sagemath-ts/sagemath-ts/zk` import is not possible — the symbols are reachable from the package root (`export * as zk`) and via `./rings` |
| `src/zk/` (`sumcheck.ts`, `multilinear.ts`) | `zk/` | Ports of `reference/sage_blueprints/`, not of SageMath |
| `estimateBKZBlockSize`, `bkzRootHermiteFactor`, `qaryLattice`, `qaryDualLattice` | `modules/free_module_integer.ts` | Port-invented. The BKZ estimator uses interpolated Gama-Nguyen/Chen values below `beta = 40` and the standard asymptotic formula above; it is a **heuristic** and says so in its docstring. Previous implementations were simply wrong (`only-2-or-n`, and `q*Z^n`) |
| `fold(codeword, challenge, domain?)`, `fold_domain(length)` | `coding/reed_solomon.ts` | FRI folding has no SageMath counterpart |
| `error_correction_capability()` | `coding/goppa_code.ts` | SageMath's `GoppaCode` registers no decoder at all. Returns the *decoder radius*: `deg(g)` in characteristic 2 with squarefree `g` (Patterson), else `floor(deg(g)/2)` (key equation). `distance_bound()` is exactly Sage's `1 + deg(g)`; a caller wanting `(d-1)/2` must compute it from that |
| `sampleExact()`, `samplesExact()`, `isIntegral`, `basisExact`, `cExact` | `stats/distributions/discrete_gaussian_lattice.ts` | See [Discrete Gaussian Samplers](#discrete-gaussian-samplers) |
| `Integer.nth_root_mod` | `rings/integer_ring.ts` | Sage exposes this as `Mod(a, p).nth_root(n)`; there is no `Integer` method, so there is no doctest to match. When `a == 1` and `gcd(n, q-1) > 1` we return 1 where Sage's `_nth_root_common` returns a *primitive* gcd-th root of unity (`Mod(1,11).nth_root(5)` is 4 in Sage). Every returned value satisfies `r^n == a (mod p)` |
| `NumberFieldElement.is_integral_unit()` | `rings/number_field/` | Sage overloads `is_unit()` on the parent's type; without a distinct element class for orders the two behaviours cannot share one name |
| `QuadraticField.d` | `rings/number_field/` | The squarefree part. `.D` is Sage's `D`. Retained only so existing callers do not silently get the wrong number; **treat as deprecated** |
| `Ratio`, `isRatio` | `parigp-ts/src/elliptic/init.ts` | Exact rational j-invariants; not yet re-exported from the package root |
| `Fp_ellcard_Shanks` (exported) | `parigp-ts/src/elliptic/group.ts` | `static` in PARI; exported so the BSGS branch is testable against an exhaustive oracle |
| `Fp_elldivpol(l, a4, a6, p)` | `parigp-ts/src/elliptic/advanced.ts` | The division polynomial `psi~_l` over `F_p`. PARI's SEA never builds one (it uses modular polynomials), so there is no upstream counterpart. Exported so the recursion is testable in isolation |
| `qfrep(Q, bound)` | `stats/distributions/discrete_gaussian_lattice.ts` | Fincke-Pohst theta series. Sage delegates to PARI; `parigp-ts` has none |
| `matkermod_basis`, `zm_from_rows`, `zm_to_rows` | `parigp-ts/src/matkermod.ts` | Row-major convenience wrappers around PARI's column-major `ZM` layout, so a caller cannot silently transpose the meaning of "the kernel basis" |
| `realQuadraticFundamentalUnit(K)` | `rings/number_field/unit_group.ts` | Sage exposes the fundamental unit only through `K.units()` / `UnitGroup.fundamental_units()`. The quadratic case bypasses `bnfinit` entirely (as PARI's own `bnfinit` does), so it deserves a directly testable entry point |
| `x_list(prec)`, `y_list(prec)` | `schemes/elliptic_curves/formal_group.ts` | Coefficient accessors for the formal `x`/`y`. `x(prec)`/`y(prec)` still return a `LaurentSeriesElement`, but that class exposes no `list()`/`__getitem__` and its `power_series()` throws for negative valuation, so the values were **unreadable**. Purely additive; the right long-term fix is to give `LaurentSeriesElement` arithmetic and accessors in `rings/power_series_ring.ts` |
| `PythonRandom` | `misc/randstate.ts` | A port of CPython's `random.Random`, which Sage reaches through `randstate.python_random()` rather than exposing as a class |
| `pAdicEisensteinQuadraticExtension`, `pAdicEisensteinQuadraticElement` | `schemes/elliptic_curves/padic_lseries.ts` | Sage builds this with `K.extension(f, names='alpha')` in `sage/rings/padics/`. Ours lives beside its only consumer because `rings/padics/padic_generic.ts`'s `pAdicExtension` is a shell with no element type |
| `Frobenius_filter` | `schemes/elliptic_curves/isogeny_class.ts` | Sage's is in `gal_reps_number_field.py`, which is not ported |
| `FractionFieldElement`, `tensorProductVector` | `modules/free_module.ts` | Elements of `Frac(R)` for a Euclidean base ring (e.g. `QQ(x)`), and the elementary tensor. Needed to represent coordinates over `QQ[x]` |
| `change_ring(matrix, ring)` | `matrix/matrix_operations.ts` | Sage's is a `Matrix` method in `matrix0.pyx`, which maps to `matrix_generic.ts`; it landed here for file-ownership reasons and is **not** re-exported from `matrix/index.ts` in 0.0.12 — import it from `./matrix_operations.js` |
| `pivots` | `matrix/matrix_decompositions.ts` | Companion to `pivot_rows` (which now correctly returns **row** indices). Re-exported from `matrix/index.ts` as of 0.0.12 |
| `PadeApproximant`, `BivariatePowerSeries` | `rings/power_series_ring.ts`, `schemes/elliptic_curves/formal_group.ts` | Stand-ins for `Frac(R[z])` and `PowerSeriesRing(R, 2)`; see [Keyword Arguments, Return Shapes and Signature Adaptations](#keyword-arguments-return-shapes-and-signature-adaptations) |
| `GF2X.rep()` | `ntl-ts` | Exposes the packed bigint so `sagemath-ts` can convert cheaply |
| `monomial_coefficient(exponentTuple)` | `rings/polynomial/multi_polynomial_element.ts` | A pure **superset** of Sage's signature (which takes a monomial with the same parent, and is fully supported); the tuple form is this port's internal representation and is what ~32 call sites use |
| `multilinear.ts: sparseMultilinearExtension([i])` | `zk/` | Returns the selector `eq(i, x)`; the blueprint short-circuits to the constant `R(i)` — an inconsistency with every other input shape, and its branch still contains leftover `print` debugging |

### Rationale

1. **ZK-specific functionality** is the project's stated focus and has no SageMath equivalent
2. **Testability** - a few upstream-`static` symbols are exported so oracle tests can reach them
3. **Migration safety** - a few (`.d`, `RandState.random`, the duplicate polynomial sampler) are
   retained deliberately so that existing callers fail loudly rather than silently

### Trade-offs

- Users reading the mirrored file layout can mistake these for SageMath APIs
- Some are marked deprecated but not yet removed

### Behavioral Impact

None on any SageMath-named function. Each symbol above carries an `@see Deviation:` docstring.

---

## Keyword Arguments, Return Shapes and Signature Adaptations

Python keyword arguments map to options objects (DESIGN.md); Python tuples, dicts and
polymorphic parents map to TypeScript equivalents. The cases below change an observable *shape*
and are therefore registered individually.

| Function | SageMath | sagemath-ts |
|----------|----------|-------------|
| `power_series_ring.pade(m, n)` | Returns `u / v` in `Frac(R[z])` | A `PadeApproximant` object with `numerator()`, `denominator()`, `power_series(prec)`, and a `toString()` rendering exactly like SageMath's fraction-field repr |
| `formal_group.group_law()` | Element of `PowerSeriesRing(R, 2, 't1,t2')` | A `BivariatePowerSeries` class (total-degree truncation, `add`/`sub`/`mul`/`pow`/`inv`/`subs`/`coefficient`/`toString`); `toString` reproduces Sage's textual form |
| `padic square_root`, `nth_root` | `square_root(extend=True, all=False, algorithm=None)`, `nth_root(n, all=False)` | Options object `{extend, all}` with overloads, plus `square_root_all()` / `nth_root_all()` |
| `matrix indefinite_factorization` | `(L, vector(L.base_ring(), d))` | `[Matrix, R[]]` (was `[Matrix, Matrix]` with `D` as a diagonal matrix) |
| `matrix krylov_basis` / `krylov_kernel_basis` | `(M, shifts, degrees, output_rows=True, var=None, basis_algorithm=None)`; `output_rows` means "also return the row coordinates" | Same, and the polynomial form now works — its argument is named `variable` because `var` is a reserved word in JavaScript, and it accepts a string or a `PolynomialRing`. Previously `output_rows` meant "return as rows, else transpose" and the return type was always `Matrix` — a misreading of the docstring; `degrees` also now accepts a single integer |
| `QQ(...)` (`RationalField.__call__`) | `Rational.__set_value` (`rational.pyx:591-704`) accepts `Rational`, `Integer`, `int`, `float`, `str`, a length-2 `[n, d]` pair, a length-1 list, `fractions.Fraction`, anything with `_rational_()`, and the two-argument `(n, d)` form | **Rewritten in 0.0.12** as a port of that method, with Sage's exact `TypeError('unable to convert {!r} to a rational')` and `ValueError('denominator must not be 0')`. Two additions: the pre-existing two-argument spelling is kept (Sage passes `d` as `Rational.__init__`'s *base*, and string bases are unsupported here, so there is no conflict), and `{numer, denom}` objects are accepted alongside Sage's `{numerator, denominator}` shape because internal code had adopted the former. **Source-breaking:** code catching `ZeroDivisionError` from `QQ('1/0')` or `QQ(n, 0n)` now sees `TypeError`/`ValueError` — the old messages (`'cannot convert object to Rational'`, `'tuple must have exactly 2 elements'`) were invented |
| `matrix density()` | Exact rational (`2/3`) | Exact `Rational` (was a JS float) |
| `matrix norm(p)` | `RDF` element | JS `number` — which *is* what RDF is |
| `matrix is_similar(..., transformation=true)` | `(False, None)`; the returned `T` satisfies `A == T.inverse()*B*T` (`matrix2.pyx:12831`) | `[false, null]`; return type widened to `boolean \| [boolean, Matrix<R> \| null]`. Sage's convention is preserved — the returned `T` satisfies `A == T^-1 · B · T`, **not** `P^-1 A P == B` |
| `matrix frobenius_form(2)` | Two elements of `MatrixSpace(QQ, n)` | `[Rational[][], Rational[][]]` — `B` is genuinely rational and the port has no rational-matrix class. The empty matrix yields `[[], []]` |
| `formal_group x(prec)` / `y(prec)` | A Laurent series supporting `.list()`, arithmetic and composition | Unchanged (`LaurentSeriesElement`), which as of 0.0.14 **does** support `.list()`, `__getitem__`, `coefficients()`, `exponents()`, arithmetic and composition. The `x_list(prec)` / `y_list(prec)` accessors added while it did not are kept (they return the coefficients from valuation −2, resp. −3) and are purely additive |
| `DGL sigma` / `c` | Methods `D.sigma()`, `D.c()` | Converted from public properties to methods in 0.0.12, matching Sage; `cNumeric()` added as a float view |
| `binary_qf reduced_form(algorithm=…)` | `reduced_form(transformation=False, algorithm='default')` | `reduced_form({ transformation?, algorithm? })`, identical semantics and identical error messages for all three algorithm values |
| `binary_qf solve_integer(n)` | `n` may be an `Integer` **or** a `Factorization` | `solve_integer(n, { factorization? })` — `n` is always an integer and the factorization is passed alongside, mirroring PARI's own `[n, factor(n)]` input. There is no `Factorization` class in the port |
| `Integer.multiplicative_order(n)` | No such signature — Sage's takes no argument | The two-argument form was **removed**; the modular computation lives on as `Mod(a, n).multiplicative_order()` (which already existed and is factorization-based) and as a private helper for `is_primitive_root`. Source-breaking for any caller of the old form |
| `Integer.log(b)` | Exact `Integer` only when `b^k == self`, else a real/symbolic logarithm | Always `floor(log_b(self))` (an alias for `exact_log`) — the port has no symbolic ring. Callers can test exactness with `b ** result === self` |
| `MPolynomial.monomial_coefficients()` | `dict` keyed by `ETuple` | `Map` keyed by the canonical comma-joined exponent string from `exponentToKey()`, padded to `parent.ngens_value`; `keyToExponent()` recovers the tuple. JS `Map` keys compare by identity, so array keys would never collide |
| `MPolynomial.coefficient({x: 1, y: 1})` | Dict keyed by the generator objects | `Record` keyed by variable **name**. Sage's list-with-nulls and monomial forms are byte-for-byte |
| `crypto.gen_lattice(quotient=…)` | A symbolic expression or univariate polynomial; validates the parent (`TypeError: quotient should be a univariate polynomial`) | `IntegerLike[]`, the ascending coefficient list (`x^4 - 1` is `[-1n,0n,0n,0n,1n]`). Degree checks and Sage's other messages are preserved; additionally throws `TypeError('quotient must be monic')` |
| `Matrix.subs` | Transmits to the entry's `subs` | Falls back to the entry's `evaluate` (which is how univariate `Polynomial` exposes substitution here) before raising |
| `Matrix.denominator` | LCM of entry denominators, an element of the denominator ring (ZZ for a QQ matrix) | Reads `denominator` whether method or getter, takes the LCM in `bigint`, returns `ring(lcm)` — so a QQ matrix gives the *rational* 30. Previously it always threw |
| `block_matrix` | Flat list + nrows/ncols, ragged list, or list of lists | List of lists only; a ragged one raises Sage's own `ValueError('list of rows is not valid (rows are wrong types or lengths)')` |
| `permutation_normal_form(check=true)` | `is_permutation_of(MS_max, True)[1]`, a pair of 1-based `PermutationGroupElement`s | 0-based index arrays with the convention `normal_form[i][j] === matrix[row_perm[i]][col_perm[j]]`. When the matrix has non-trivial automorphisms several permutations realise the same normal form; the matrix returned is identical |
| `zk sumcheckVerify` | Round count comes from `len(poly.args())` in the driver; `degree_checks` defaults to `None` | `(proof, claimedSum, polyEvaluator, field, numVars, options?: {degreeCheck?})`. `numVars` is **required** and a proof with a different number of rounds is rejected — a stateless verifier cannot recover the round count from the proof without trusting the prover. `degreeCheck` now defaults to no check (it was hardcoded to 1, which made legitimate higher-degree GKR rounds unverifiable) |
| `zk binaryToInt` | No counterpart (Python ints are unbounded) | Returns `bigint` and throws `ValueError('bits must contain only 0 and 1 (got X)')`. `parseInt(bits.join(''), 2)` lost precision above 2^53 and silently truncated `[1,2,1]` to 1 |
| `NumberFieldIdeal.norm()` | `QQ` | `Rational`. `prime_above(p, {degree})` returns one ideal, `primes_above(p)` the list, `decomposition(p)` the `[P, e]` pairs with `bigint` exponents |
| `ClassGroup.exponent()` | Sage/PARI list invariants in decreasing order (`d_{i+1} \| d_i`), e.g. `(38, 2)` | `max(invariants)` rather than the last entry — correct under either ordering convention |
| `convolution()` | Any commutative ring in which multiplication by two is injective | `bigint` natively, plus ring elements exposing `.parent` with `__call__` and `div`/`inv`; anything else raises `NotImplementedError`. The port's `RingElement` carries no `parent`/`zero`, so `R(x/M)` cannot be expressed generically. Every ring element class in this repo exposes `.parent`, so practical coverage is the same |
| `QuotientRing` non-invertible element | `ArithmeticError('element is non-invertible')` on the Singular path; `ZeroDivisionError(f'element {self} of quotient polynomial ring not invertible')` on the fallback | Always the `ZeroDivisionError` fallback message (there is no Singular path) |
| `IntegerMod.log(b, order)` | `log(self, b=None, order=None, *, check=False)`; `order` is consulted **only** when `check=True` | Same now; previously `order` was used as the base order for a prime-power shortcut that Sage does not have |
| `matrix.companion_matrix`, `toeplitz`, `hankel` | See `sage/matrix/special.py` | Argument conventions corrected to Sage's (full monic coefficient list with negated border; `r` counted from the second column with `ncols = len(r)+1`). The port's own tests had pinned the wrong conventions |

### Rationale

1. **DESIGN.md mapping** - keyword arguments become options objects; tuples become arrays;
   `None` becomes `null`
2. **No polymorphic parents** - TypeScript cannot express "element of the base ring", so
   accessors split by concrete type
3. **JS `Map` semantics** - object and array keys compare by identity, so canonical string keys
   are the only correct sparse representation
4. **Soundness** - `sumcheckVerify`'s required `numVars` closes a forgery: a short proof
   verified nothing

### Trade-offs

- Source-breaking for callers of `Integer.multiplicative_order(n)`, `sumcheckVerify`,
  `binaryToInt`, `indefinite_factorization`, `krylov_*` and `NumberFieldIdeal.norm()`
- Dictionary-keyed `coefficient()` calls need a string instead of a generator
- `PadeApproximant` and `BivariatePowerSeries` should be replaced when real
  `Frac(R[z])` / `MultiPowerSeriesRing` types exist

### Behavioral Impact

Values are SageMath's; only shapes and argument positions differ. Every `toString()` listed
above reproduces SageMath's printed form.

---

## Multivariate Polynomial Term Orders Restricted

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Term orders | `TermOrder` supports twelve orders (lex, invlex, deglex, degrevlex, neglex, negdegrevlex, negdeglex, degneglex, wdeglex, wdegrevlex, negwdeglex, negwdegrevlex) plus block orders and weighted gradings; `TermOrder('royalorder')` raises `ValueError` | Only the three global orders `lex`, `deglex`, `degrevlex`, as a bare string union rather than a `TermOrder` object |
| Unknown order | `ValueError` | `ValueError("unknown term order 'name'")` — Sage's exact message. Previously it silently fell back to degrevlex |
| Weighted grading / `degree(std_grading=…)` | Supported | No meaning here |
| Affected modules | `sage/rings/polynomial/term_order.py`, `multi_polynomial_element.py` | `packages/sagemath-ts/src/rings/polynomial/multi_polynomial_element.ts` |

### Rationale

Silently degrading an unknown order to degrevlex produces leading terms, S-polynomials and
Gröbner bases for a *different* ordering than the caller asked for — wrong answers with no
error. Rejecting is strictly safer and matches Sage's behaviour for names it does not know.
Implementing the remaining nine orders was out of scope for this pass.

### Trade-offs

- Nine of Sage's term orders, block orders and weighted gradings are unavailable
- `TermOrder` is not an object, so its methods are unavailable

### Behavioral Impact

Code passing any order outside the three supported ones now raises instead of computing against
degrevlex. Roughly 19 of Sage's `MPolynomial` methods are honest `NotImplementedError` stubs
rather than implementations.

---

## Free Module Exactness and Coordinate Types

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Linear algebra | Exact over any PID, with echelon forms over e.g. `QQ[x]` | Exact fraction-field layer for `bigint`/`number`/`Rational` entries and for fields providing `div`/`inv`, **plus, as of 0.0.12, univariate polynomial rings over a field** (`QQ[x]`, `GF(p)[x]`): a new `FractionFieldElement` (exact `QQ(x)`), a port of `Matrix._echelon_form_PID` (`matrix2.pyx:17305`) and `_generic_clear_column` (`:20613`), and a port of Sage's generic Smith normal form (`:16732`, `:20730`, `:20537`) for kernels and intersections. Over any *other* ring, echelon forms, coordinates, kernels and linear dependence still raise `NotImplementedError('exact linear algebra is not implemented over this base ring')` |
| `tensor_product` | **No such method on free modules.** The only concrete embedded definition upstream is `free_quadratic_module_integer_symmetric.py:1343`, built on `Matrix.tensor_product` (`matrix2.pyx:9983`) | A port of that: rank `m·n` submodule of `R^(deg1·deg2)` whose user basis is the Kronecker product of the basis matrices, with the Kronecker product of the inner product matrices, plus Sage's `discard_basis=True` variant. **The 0.0.11 `@see Reference: FreeModule_generic_pid.tensor_product` citation was fabricated** — no such method exists |
| `quotient` over a field | `FreeModule_generic_field.__quotient_matrices` (`free_module.py:5366`), returning a `FreeModule_ambient_field_quotient` of dimension `dim V - dim W` with a quotient map and section | Ported exactly (basis extension by pivot rows of `B.stack(S)`, `Q = D[:, n-m:n]`, `L = D^-1[n-m:n, :]`) |
| `quotient` over ZZ | `FGP_Module` (`fg_pid/fgp_module.py:268`): relative matrix, Smith normal form, invariants, reduction mod the invariants | Ported exactly, including `invariants()` with and without ones, the Smith generators, `cardinality()` |
| `quotient` over any other ring | Various | Sage's exact `NotImplementedError('quotients of modules over rings other than fields or ZZ is not fully implemented')` (`free_module.py:4472`) |
| `coordinate_vector` | Element of `FreeModule(R.fraction_field(), rank)` — always a `Rational` for a ZZ module | `bigint` when integral, `Rational` otherwise (`number` over a JS-number base ring) |
| `coordinates(check=False)` | Skips verification and can return a vector that does not reconstruct `v` | Always raises `ArithmeticError('vector is not in free module')` when `v` is outside the span |
| `indexIn()` | Base-field element or `infinity` | `bigint`/`Rational`, or `Number.POSITIVE_INFINITY` |
| `cardinality()` | Sage `Integer` or `+Infinity` | `bigint` (cardinalities routinely exceed 2^53 — `GF(2)^70`) or `Number.POSITIVE_INFINITY` |
| `norm(p)` for irrational results | Symbolic (`sqrt(14)`, `276^(1/5)`) | Exact `bigint`/`Rational` whenever the p-th root is rational, a double otherwise |
| `normalized()` | `v / v.norm(p)`; the base ring changes (Sage's own docstring says so), usually to the symbolic ring | A vector over QQ (`Rational` entries), or over the double field when the norm is irrational |
| `discriminant()` | `FreeModule(R,n)` uses `det(gram)`; `FreeModule(R,n,inner_product_matrix=A)` is a `FreeQuadraticModule` whose discriminant is `(-1)^(rank//2) * det(gram)` | Same split, keyed on whether an inner product matrix was supplied (the port merges `free_quadratic_module.py` into `free_module.ts`, so the class distinction becomes a runtime condition) |
| `submodule(gens, check=True)` | `ArithmeticError('argument gens (= …) does not generate a submodule of self')` | Implemented (possible now that `isSubmodule` is exact); `span()` remains unchecked as in Sage |
| Affected modules | `sage/modules/free_module.py`, `free_quadratic_module.py`, `free_module_element.pyx` | `packages/sagemath-ts/src/modules/free_module.ts`, `free_module_element.ts` |

### Rationale

1. **The alternative was silent double-precision arithmetic** over arbitrary rings, which is
   what produced the previous defects: `span()`/`subspace()` used the *number of generators* as
   the rank, and the float RREF/kernel/determinant helpers rounded doubles back into `bigint`
2. **Consumer expectations** - every caller of a ZZ module wants `bigint`; carrying `Rational`
   everywhere would change the entry type of the whole module
3. **`check=False` cannot be reproduced usefully** - Sage's unchecked partial answer needs the
   rref-pivot transformation machinery and yields a meaningless vector; raising never differs
   when `v` is in the module
4. **Delegation** - Hermite normal form goes to `matrix_integer.hermite_normal_form` and
   saturation to `matrix_integer.saturation`, exactly where SageMath delegates to
   `Matrix_integer_dense`

### Trade-offs and source-breaking changes in 0.0.12

- **`FreeModuleField.quotient()` now returns the quotient space, not a complement.** It
  previously returned an arbitrary complement *subspace of `self`* — mathematically the wrong
  object, carrying no maps. The declared return type changed from `FreeModuleField` to
  `FreeModuleQuotient`, and the result's `degree` is now the quotient dimension. The only call
  site in the repo (a `dimension() === 2` assertion) still passes
- **`FreeModuleQuotient`'s `degree`, `rank` and `basis` changed.** Sage has two distinct classes
  (`FreeModule_ambient_field_quotient` and the torsion-carrying `FGP_Module`); the port folds
  both into one, so over a field `rank = degree = dim V - dim W` (previously `degree` was
  `cover.degree()`, making elements the wrong length) and over ZZ `degree` is the number of Smith
  invariants with `rank` the free rank
- **`tensorProduct`'s degree changed** from `rank(M)·rank(N)` to `deg(M)·deg(N)` (equal for
  ambient modules) and the result now has a basis rather than being ambient
- **`isPID()` now recognises univariate polynomial rings over a field**, so `FreeModule(QQ[x], n)`
  is a `FreeModuleAmbientPID` where it was a bare `FreeModuleAmbient`. Strictly more capable
- `intersection()` over a non-ZZ PID returns the mathematically correct module (250/250 random
  cases verified **equal in SageMath**, using Sage's own module equality) but in ~25 % of cases
  its **echelon basis** differs from Sage's by a unit of the base ring (a nonzero rational
  constant). Two upstream normalisations are not reproducible from the Euclidean interface:
  `_echelon_form_PID` is itself not canonical up to units, and Sage's `integer_kernel`
  (`matrix2.pyx:5646`) scales by `Matrix.denominator()` — for a `QQ[x]` matrix the lcm of the
  *coefficient* denominators, a notion living in `Frac(ZZ)` rather than in the Euclidean structure
  of `QQ[x]`. The module generated is always identical; only the representative differs
- `_echelon_form_PID` omits the reduction above the pivots (`matrix2.pyx:17419-17426`), for the
  same reason Sage omits it for `K[x]`: polynomial ideals have no `small_residue`, so Sage's own
  `except AttributeError` swallows that step
- The internal `FractionFieldArithmetic.denominator()` was widened from `bigint` to a base-ring
  element, so `FreeModulePID.denominator()` over `QQ[x]` returns a polynomial (e.g. `x`)

### Behavioral Impact

Ranks, echelon bases, coordinates, kernels, intersections, complements, discriminants and
cardinalities are exact and match SageMath's values. New in 0.0.12: 700 random `QQ[x]`/`GF(p)[x]`
spans match the vendored SageMath's echelon basis exactly; Sage's issue-#9053 doctest, the
fractional-entry span, `vector_space_span` over `QQ(x)` and the `QQ[x]` quotient error all
reproduce. Quotients were checked against the locally installed Sage on 200 random QQ and 191
random ZZ cases — **every** projection, lift, invariant and cardinality agrees exactly (so the
port even picks the same section), with `project(lift(x)) == x`, `project(W) == 0` and additivity
executed on all 391. The `IntegralLattice("D3")` tensor-product doctest reproduces number for
number, including the Gram matrix and the `discard_basis` variant.

`FreeModuleSubmodule`, `FreeModuleSubmodulePID`, `FreeModuleSubspace`, `FreeModuleQuotient` and
`FractionFieldElement` are now re-exported from `modules/index.ts` (they were missing).

---

## Matrix Module Algorithm Substitutions

Where a SageMath matrix routine delegates to a backend this port lacks, an equivalent exact
algorithm is used instead. The *results* are identical unless stated.

| Function | SageMath | sagemath-ts |
|----------|----------|-------------|
| `is_positive_definite` / `is_positive_semidefinite` | Reads eigenvalue signs off the 1×1 diagonal blocks of the Bunch-Kaufman `block_ldlt` factorization | After the same ring check and `is_hermitian` test, reads them off the characteristic polynomial: with `charpoly = sum c_i x^i`, the elementary symmetric functions are `e_k = (-1)^k c_{n-k}`, and a Hermitian matrix is positive (semi)definite exactly when every `e_k` is `> 0` (`>= 0`). Provably equivalent, exact, one division-free charpoly |
| `is_similar` | `A.rational_form() == B.rational_form()` | Compares, for every monic irreducible factor `h` of the charpoly, the multiset of elementary-divisor exponents recovered from the growth of `dim ker(h(A)^k)` — precisely the data determining the rational canonical form. Also detects a 6×6 counterexample that charpoly+minpoly misses |
| `is_similar(transformation=true)` | `matrix2.pyx:13052-13070`: Jordan forms over the fraction field, then over the algebraic closure, else `RuntimeError('unable to compute transformation for similar matrices')` — its own doctest at `:12918-12922` shows two *provably similar* matrices over `GF(7^2)` for which Sage raises | Tries Sage's Jordan-form formula first (so it will match Sage exactly), then falls back to solving the intertwining equation `B X = X A` as an `n²×n²` homogeneous system and searching the kernel for an invertible solution. **This succeeds in cases where Sage raises.** Every candidate is verified (`B·T == T·A` and `rank(T) == n`) before being returned, so a wrong transformation is impossible. The particular `T` is generally not Sage's — it is not canonical in Sage either |
| `norm(A, 2)` | `matrix2.pyx:16466-16471`: `change_ring(CDF)`, `A^H·A`, a numerical SVD, `max(S).real().sqrt()`, returning `RDF` | **The same route as of 0.0.14**: `_entryToCDF` mirrors `change_ring(CDF)`, then `A^H·A` and an SVD in double precision. RR and CC entries are accepted. Number-field entries still raise `NotImplementedError` (no distinguished complex embedding is wired into `_entryToCDF` yet); rings of positive characteristic raise `TypeError`, as `change_ring(CDF)` does upstream |
| `jordan_form(transformation=true)` | `matrix2.pyx:12259-12312` + `_jordan_form_vector_in_difference` (`:20895`) | Ported line for line, including the detail that makes it reproduce Sage's **exact** `P` rather than merely a valid one: Sage's `right_kernel().basis()` is echelonized, and the chains depend on which kernel vector is picked first. Eigenvalues now come from `A.charpoly().roots()` (Sage's own route, `:12228`) with the old `ring.elements()` enumeration kept only as a fallback — that enumeration is why `jordan_form` could not run over QQ at all |
| `jordan_decomposition` | `matrix2.pyx:12383-12400`: a Newton iteration on the minimal polynomial (`h = f/gcd(f,f')`, `h.xgcd(h')`, then `A -> A - h(A)·q(A)`), which succeeds **even when the eigenvalues are not in the base field** — its own doctest at `:12332` does exactly that | Reads `D` and `N` off the Jordan form. Correct whenever the eigenvalues lie in the base field; otherwise it propagates `jordan_form`'s `ArithmeticError` instead of returning the decomposition. A `catch {}` that silently swallowed every failure and reported "not fully implemented" was removed |
| `krylov_kernel_basis` | `matrix2.pyx:20343-20478` | Ported. **Correction to the record:** the 0.0.11 note claiming this "needs Popov/approximant bases" was wrong — Sage builds the kernel directly from the Krylov basis as `relation = D·C^-1`. Porting it fixed the constant output's **row order** and the third component of `row_profile` as a side effect (see below) |
| `change_ring(matrix, ring)` | `matrix0.pyx:1666-1715`, relying on the coercion framework | The port has no coercion framework, so `_coerce_entry` asks the target ring to convert and, on failure, builds the canonical morphism itself in exactly the two cases where one **provably** exists: `QQ -> R` (`n/d ↦ R(n)·R(d)^-1`, raising if `R(d)` is not a unit) and `Z/mZ -> R` when `char(R) | m` (or `m = 0`). Anything else raises `TypeError` rather than guessing — so `Z/8 -> GF(7)` is correctly refused |
| `QR` | Scales each column by `1/sqrt(<v,v>)`, so `Q` is unitary and `R` has non-negative diagonal; raises `TypeError` when the fraction field has no square roots | Delegates to `gram_schmidt_noscale`: `Q`'s columns are **orthogonal but unnormalized** and `R` differs by the corresponding diagonal factor. `full` now defaults to `true` as in Sage |
| `block_ldlt` pivot selection | Bunch-Kaufman compares `\|A_kk\|`, `omega_1`, `omega_r` against `alpha = (1+sqrt 17)/8` in C doubles | Exactly that rule when the base ring's elements expose `abs()` (QQ, RR, ZZ); over rings with no absolute value (finite fields) an exact rule: 1×1 pivot if `A_kk != 0`, else swap in the first `r` with `A_rk != 0` and take 1×1 if `A_rr != 0`, else 2×2. Pivot choice affects only numerical stability, never correctness — `P^T A P == L D L^T` holds exactly |
| `principal_square_root` | Returns `False` when `check_positivity` and not positive definite; works over the algebraic closure via `eigenmatrix_left` | Skips the positivity check and diagonalizes over the base ring; raises `ArithmeticError` when not diagonalizable there. Over a finite field "the" square root of an eigenvalue is defined only up to sign, so the result is *a* square root rather than *the* principal one |
| `is_permutation_of` / `permutation_normal_form` | `BipartiteGraph.is_isomorphic(..., edge_labels=True)` (bliss/nauty) | Complete backtracking over the row assignment, pruned by the column-multiset invariant of every prefix plus row/column signature filters; `permutation_normal_form` does Sage's row-by-row maximisation with column-block refinement. Both outputs are *uniquely specified* (a boolean; the lexicographic maximum over all row/column permutations), so any complete algorithm agrees — both were checked against exhaustive brute force. Worst case is exponential, as is Sage's |
| `pluq` / `ple` | `P` and `Q` come from M4RI's `mzp_t`, i.e. transposition lists applied in order | Both `P` and `Q` are now transposition lists (`P[pivotRow] = foundRow`). `P` was previously the *composed* permutation — invisible in the doctest because `P` is the identity there, but both entries must use one convention for `A = P*L*U*Q` to hold |
| `matrix_operations.pivot_rows` | Row indices | Row indices (it previously returned pivot **columns**) |
| Random matrix constructors | `sage.misc.prandom` (randint/shuffle) driven by the global randstate | `current_randstate().randint(...)` and a randstate-driven Fisher-Yates shuffle / density fraction (was `Math.random()`), so `set_random_seed` makes results reproducible; the stream still differs from Sage's |
| Affected modules | `sage/matrix/matrix2.pyx`, `matrix0.pyx`, `special.py`, `matrix_misc.py` | `packages/sagemath-ts/src/matrix/` |

### Rationale

1. **The delegation targets are missing or broken** - `rational_form` is a stub;
   `block_ldlt` was itself producing invalid factorizations; there is no graph package and no
   SVD. (A generic `change_ring` now exists, and `jordan_form(transformation=true)` landed in
   0.0.12, so `is_similar` and `is_diagonalizable` no longer raise.)
2. **Provable equivalence** - each substitution computes the same mathematical object by a
   different route, and each was verified against upstream doctests and randomized sweeps
3. **Exactness** - the charpoly criterion for definiteness introduces no floating point, per
   CLAUDE.md. `norm(2)` is the exception: upstream is explicitly inexact there
   (`change_ring(CDF)` + SVD), so following it faithfully means following it in double
   precision

### Trade-offs

- `QR`'s `Q` is not orthonormal, so callers expecting a unitary matrix must normalize
- `permutation_normal_form(check=true)` may return a different (equally valid) permutation
- Bunch-Kaufman's permutation may differ from Sage's over finite fields
- `principal_square_root` returns a non-principal root over finite fields
- `norm(A, 2)` is **inexact, like Sage's**, and does not accept number-field entries: once a
  `NumberField` carries a distinguished complex embedding, `_entryToCDF` should route through it
  as `change_ring(CDF)` does. Rings of positive characteristic raise `TypeError`, as upstream
- `jordan_form` **honours `subdivide`** (0.0.14), but the subdivisions are set directly rather
  than through `block_diagonal_matrix(..., subdivide=…)` (`:12255`), because
  `matrix_special.block_diagonal_matrix` has no `subdivide` option
- `jordan_form` raises `ArithmeticError` where Sage raises `RuntimeError`, and so does
  `is_similar` on the (unreachable) double failure. `errors.ts` gained a `RuntimeError` in this
  pass but these two sites were not switched
- Sage's `jordan_form` doctests over `PolynomialRing(QQ, 'x11,…')` and over
  `FractionField(PolynomialRing(QQ,'a'))` cannot be run: the port cannot build a `Matrix` over a
  multivariate polynomial ring or a rational function field. The first was reproduced by
  instantiating `x31 = 5`, `x21 = 7` over QQ and matches Sage entry for entry
- The shifts-Popov property of `krylov_kernel_basis`'s polynomial output **is now asserted
  independently** against `matrix_polynomial_dense.ts` (`is_popov`, `popov_form`, `is_hermite`,
  `minimal_approximant_basis`), which landed in 0.0.14
- **`Matrix.toString` is not subdivision-aware and pads per column.** SageMath pads every entry
  to one global width and draws `|` / `---+---` separators (`matrix0.pyx:2180`). The port has a
  faithful `matrix_str` in `matrix_decompositions.ts`, but `matrix_generic.ts:429`'s `toString`
  does not delegate to it; `jordan_form` attaches `matrix_str` as a per-instance `toString` on
  the subdivided `J` it returns, so only those matrices print like Sage's. That per-instance
  override is a stopgap. Related: `_subdivisions` is not preserved by `Matrix.copy()`, whereas
  SageMath preserves subdivisions under `copy()`
- `krylov_kernel_basis`'s argument is named `variable`, not Sage's `var`, because `var` is a
  reserved word in JavaScript
- Random matrices differ from Sage's for the same seed
- `is_hermitian` is exported from `matrix_operations.ts` but not re-exported from
  `matrix/index.ts`; neither is the new `change_ring`

### Mitigation

Once `rational_form` and `block_ldlt` exist, swap `is_similar` and the definiteness predicates
for the upstream paths with no visible change. Rewrite `jordan_decomposition` to Sage's
minimal-polynomial Newton iteration so it works for non-split characteristic polynomials.

### Behavioral Impact — including two meaning changes

Values match SageMath's, verified by execution against Sage's own printed output: the
issue-12693 `jordan_form` doctest reproduces Sage's `P = [2 1 0/0 0 1/-2 0 -1]` character for
character; the 10×10 QQ doctest gives Sage's block structure 3,2,2,2,1 with `rank(T) = 10`; a
pinned randomized test covers 120 conjugated Jordan matrices with repeated eigenvalues and
non-trivial chains over QQ, GF(7) and GF(101). `is_similar`'s transformation doctest reproduces
Sage's `T = [[1,0,0],[-2/3,1/6,-5/6],[2/3,0,-1/3]]` entry for entry, and the similar/not-similar
verdict was checked against **exhaustive brute force over every invertible `P`** for GF(2) 2×2
(136 pairs), GF(3) 2×2 (3321) and GF(2) 3×3 (131 328) — 0 mismatches. Every printed matrix in
`krylov_kernel_basis`'s docstring is reproduced verbatim, including the deliberately-wrong
too-small-degree-bound case and the zero-`M` case.

**Two exported symbols changed meaning in 0.0.12:**

1. `krylov_kernel_basis`'s constant output previously returned an RREF pivot basis of the left
   kernel, whose rows are in echelon (arbitrary) order — for the main GF(97) doctest it produced
   Sage's rows *in reverse* — and whose `row_profile` third component was the index inside the
   delta-Krylov matrix (5 where Sage says 6). Both now match SageMath. The RREF basis destroys
   the correspondence between row `i` of the kernel and row `i` of `E`, which is what makes the
   polynomial output meaningful. No call site exists outside the re-export in `matrix/index.ts`.
2. `jordan_form`'s eigenvalue source changed from `ring.elements()` enumeration to
   `charpoly().roots()`.

Prior to 0.0.11, `minpoly` returned the minimal polynomial of `e_0` rather than of the matrix,
`echelon_form` was not the RREF, `LU` built `P` incorrectly, `gram_schmidt` orthogonalized the
wrong axis, `smith_form` never populated `V`, and `is_positive_definite` returned true for
anything symmetric with nonzero minors over GF(7) (where Sage raises).

---

## Binary Quadratic Forms

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Composition / reduction backend | `BinaryQF.__mul__` calls PARI `qfbcompraw`; `reduced_form` calls PARI `qfbred`/`qfbredsl2` for non-square discriminants | **Delegated as of 0.0.12** to `parigp-ts`'s new `qfb.ts` (a port of `pari/src/basemath/Qfb.c` and `quad.c`). The ~170 lines of in-place `_reduce_positive_definite` / `qfb_comp` / `qfb_sqr` / `dvmdii_round` are deleted. Sage's own `_reduce_indef` is retained for **square** discriminants only |
| `algorithm` selection | `algorithm = 'sage' if self.is_reducible() else 'pari'` (`binary_qf.py:947-948`), with `'default'`/`'pari'`/`'sage'` accepted explicitly | Identical, with all three of Sage's error paths. **Fidelity bug fixed:** the previous port routed *every* `D > 0` form through `_reduce_indef`, not just square discriminants. It happened to be output-compatible — on 12 507 random indefinite non-square forms PARI's `qfbred` and Sage's `_reduce_indef` returned the identical form and identical base change |
| Squaring dispatch | PARI's `qfb_comp` squares only when the two GEN **pointers** are identical (`if (x == y)`); Sage converts both operands separately (`self.__pari__()` and cypari2's `objtogen(right)`), so that path **never fires from Sage**, not even for `Q * Q` | The `this === other -> _square()` shortcut was **removed** to match. Behaviourally a no-op: `qfb_sqr` and the general `qfb_comp` agreed on all 716 self-compositions tested |
| `solve_integer` | `binary_qf.py:1608-1806`: negative-definite recursion, an elementary algorithm for square discriminants, `qfbcornacchia` for prime `n` with `disc < 0`, else `qfbsolve` with `_flag` in {1,2,3}; accepts a `Factorization` for `n` | **Newly implemented in 0.0.12** (the port had none). Ported in full; Sage's `Factorization` argument becomes an optional `{ factorization }` option, mirroring PARI's own `[n, factor(n)]` input |
| Shanks distance forms (`qfr5_*`, `qfr5_dist`) | Present | **Ported in 0.0.14.** `qfb.ts` now carries a transcription of PARI's `t_REAL` kernel (`nbits2prec`, `addrr`, `mulrr`, `divrr`, `sqrtr`, `mplog2`, `logr_abs`, `shiftr`, …) and a `QfbExt` type carrying the logarithmic distance; `qfbred`/`qfbcomp`/`qfbcompraw`/`qfbsqr`/`qfbsqrraw`/`qfbpow`/`qfbpowraw` each gained a second overload accepting a `QfbLike`. This is the branch PARI reaches with `flag \|= qf_NOD`. `buch.ts` carries a **second, independent** copy of that kernel; only `qfb.ts`'s is re-exported from the package root (see the note in `packages/parigp-ts/src/index.ts`) |
| Reduction arithmetic | `D.sqrt(prec=53)` | Exact `isqrt(D)` (see [Exact Arithmetic Where SageMath Uses Floating Point](#exact-arithmetic-where-sagemath-uses-floating-point)) |
| Affected modules | `sage/quadratic_forms/binary_qf.py` | `packages/sagemath-ts/src/quadratic_forms/binary_qf.ts`, `packages/parigp-ts/src/qfb.ts` |

### Rationale

The blocker for delegation — `parigp-ts` having no `Qfb` module — is gone. `qfb.ts` was verified
against the **real PARI 2.15.4** (reached through the local SageMath) on golden data: `qfbredsl2`
80/80 exact (form *and* base-change matrix), `qfbcompraw` 300/300, `qfbpowraw` 120/120,
`primeform` 60/60, `qfbcornacchia` 104/104, `qfbsolve` 1670/1800.

Every one of the remaining divergences is a **documented upstream change** between the oracle
(2.15.4) and the source we ported (2.18.1), not an error: `CHANGES-2.16 #45` ("changed `qfbred`
to use standard normalization, same as `qfbredsl2`") accounts for the 6 `qfbred` + 14 `qfbcomp`
+ 4 `qfbpow` cases, and `CHANGES-2.16 #9` ("`qfbsolve`: call `qfbcornacchia` when the argument is
prime and `disc < 0`") plus `allsols`' new `v >= 0` normalisation and the lexsort at
`Qfb.c:1930` account for the 130 `qfbsolve` cases. The 2.18 semantics were verified independently
rather than assumed: our `qfbred` equals our `qfbredsl2` representative on all 80 forms (exactly
what #45 states), every base change satisfies `q ∘ U == reduced` with `det U == 1`, and **all
1800** `qfbsolve` outputs satisfy `Q(x,y) == n` with the same emptiness pattern as PARI.

### Trade-offs

- **`reduced_form` for indefinite non-square discriminants now follows the vendored 2.18.1
  normalisation** (`CHANGES-2.16 #45`), i.e. the same representative as `qfbredsl2`. Code
  validated against a PARI older than 2.16, or against Sage's own output on an older PARI, may
  see a different (equally reduced) indefinite representative
- Shanks distance forms are unreachable (above)
- `solve_integer` for hard-to-factor `n` inherits `ifactor.ts`'s factoring chain (which now
  includes MPQS, so only inputs above MPQS's 107-digit ceiling still fail). The optional
  `factorization` argument is the documented workaround and is exercised by a test on `2^128+1`
- `algorithm` is an options-object field rather than a positional keyword (DESIGN.md convention)

### Behavioral Impact

Equivalence with the pre-delegation code was **proven by execution against a side-by-side import
of the HEAD file**: 29 944 random forms (coefficient magnitudes to 10^12) with 0 differences in
both the reduced form and the SL2 base change; 59 280 compositions across all 400 valid
discriminants in [-400,400] with 0 differences; 1600 `BinaryQF_reduced_representatives` calls with
0 differences; and 500 class-group Cayley tables with 0 differing tables (element orders
reproducing the known structures for `D = -4, -20, -23, -39, -56, -84, -260`).

`solve_integer` reproduces every Sage doctest exactly, and its square-discriminant branch was
cross-checked against exhaustive brute force over a 241×241 box on 600 random cases with 0 wrong
solutions and 0 false nulls.

Prior to 0.0.11, Gauss squaring was wrong and `reduced_form` returned non-reduced forms.

---

## Elliptic Curve and Isogeny Deviations

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `is_j_supersingular` | Checks `supersingular_j_polynomial(p)(j) == 0` when `p` is in the precomputed table, giving an exact answer even with `proof=False` | Skips the table (`supersingular_j_polynomial` is not ported) and always falls through to the 10 random-point tests, plus the trace-of-Frobenius check when `proof` is set (the default). With `proof=True` — Sage's and our default — the answer is identical and proved |
| Anomalous ECDLP | `ell_point.py:4640-4642`: when the base field is prime and `n == p`, `log()` delegates to `padic_elliptic_logarithm` (Smart's attack, `O(log p)`) | Falls through to the generic Pohlig-Hellman/BSGS path, returning the same value in `O(sqrt p)`. Our `padic_elliptic_logarithm` is a stub, so routing there would turn a working call into a failure |
| `division_points` for 2-torsion `P` | `ell_point.py:1531-1557` replaces `g` by `gcd(g, g') * sqrt(lc(g))` (times `(x - x(P))` for odd `m`) | Uses `g` unreduced and iterates its distinct roots. The reduction only strips repeated factors, so the *set* of roots — the only thing the algorithm consumes — is unchanged; verified equivalent on 10 512 brute-force cases |
| `montgomery_model` representative | `EllipticCurveIsogeny(GF(7) j=1728 curve, (0,0), model='montgomery')` reports `A = 1` | Returns `A = 6`, the other root of the defining cubic. Both are valid Montgomery forms; see the root-ordering note under [Polynomial Roots and Factorization Limited](#polynomial-roots-and-factorization-limited) |
| Affected modules | `sage/schemes/elliptic_curves/` | `packages/sagemath-ts/src/schemes/elliptic_curves/` |

### Rationale

1. **Unported dependencies** - `supersingular_j_polynomial` and `padic_elliptic_logarithm`
2. **Same answer, cheaper route** - the `division_points` root set and the anomalous DLP value
   are unchanged; only the algorithm differs
3. **Ties** - several equally valid representatives exist for the Montgomery model

### Trade-offs

- `is_j_supersingular(proof=False)` is probabilistic where Sage would be exact for small `p`
- Anomalous-curve discrete logs are `O(sqrt p)` instead of `O(log p)` — a real performance cliff
  on exactly the curves an attacker would target
- The Montgomery `A` can differ from Sage's printed value

### Mitigation

Port `ell_wp.py`, `isogeny_small_degree.py`, `supersingular_j_polynomial` and
`padic_elliptic_logarithm`.

### Behavioral Impact

Vélu's formulas, `division_points`, `multiplication_by_m`, `_isomorphisms` (all char 2/3/p
branches), `lift_x`/`is_x_coord`, the bivariate `division_polynomial`, `abelian_group`/`gens`,
`set_order`/`has_order`, `torsion_basis`, `twists`, `frobenius_order` and the whole formal group
now reproduce SageMath's doctests. Prior to 0.0.11, Vélu's y-coordinate carried a sign error
that put isogeny images off the codomain, and `_p_primary_torsion_basis` was not Sage's
algorithm. Note that the correct Vélu accumulation is upstream's `v += vQ`,
`w += uQ + xQ*vQ` (with `vQ = 2*gxQ` for non-2-torsion) — the audit's suggested
`w += 2(uQ + xQ*gQx)` double-counts `uQ` and yields a non-isogenous codomain.

---

## CM and Class Number Computation

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `class_number(D)` | `D.class_number(proof)`, i.e. PARI `qfbclassno` (Shanks BSGS in the form class group) for fundamental discriminants, `OrderClassNumber` otherwise | Counts reduced primitive positive-definite binary quadratic forms of discriminant `D` (`-a < b <= a <= c`, `b >= 0` when `a == c`, `gcd(a,b,c) = 1`), memoized in a module-level `h_dict`. `OrderClassNumber` is still used for non-fundamental discriminants exactly as in Sage |
| Quadratic class **group** | PARI `bnfinit` / `quadclassunit` | Enumerates reduced primitive forms of the field discriminant (definite for `D < 0`, rho-cycles of indefinite forms for `D > 0`), composes them with Dirichlet composition (Cohen 5.4.7) and reads the elementary divisors off the `\|G[p^k]\|` counts. For `D > 0` the ordinary class group is the narrow one quotiented by the class of `-f0`. Guarded by `CLASS_GROUP_DISC_BOUND = 2,000,000` |
| `cm_j_invariants` etc. return shapes | Sorted j-invariants; `(D,f,j)` triples; `(|D|, count)` pairs; `(True,(D,f))` or `(False,None)` | All four now match Sage. Previously they returned table order; `[j,{discriminant,conductor}]` pairs; a signed discriminant; and `D*f^2` or `0` |
| Affected modules | `sage/schemes/elliptic_curves/cm.py`, `sage/rings/number_field/class_group.py` | `packages/sagemath-ts/src/schemes/elliptic_curves/cm.ts`, `rings/number_field/class_group.ts` |

### Rationale

1. **No `qfbclassno` in `parigp-ts`** - counting reduced forms *is* the definition of `h(D)`,
   so the values are identical (verified against `cm.py`'s doctests and against Klaise's counts,
   which the tests now assert)
2. **Deliberate non-reuse of `binary_qf.ts`** - at the time this landed, that module had two
   confirmed critical defects being repaired in parallel; depending on it would have made this
   correctness contingent on another unit's. Those defects are long fixed and `binary_qf.ts` now
   delegates to `parigp-ts`, so **these helpers should now be folded in**; `isogeny_class.ts` in
   fact already depends on `binary_qf.ts` for `BinaryQF_reduced_representatives` and
   `BinaryQF.compose` when it builds class-group generators
3. **Explicit size guard** - enumeration is exponential in `log|D|` where PARI is subexponential,
   which brushes against CLAUDE.md's no-naive-algorithm rule, hence a hard bound rather than a
   silent slow path

### Trade-offs

- `O(|D|)` rather than PARI's `O(|D|^(1/4))`: measured 47 ms for the full `hmax = 8` sweep
  (`B = 7987`) and 647 ms for `hmax = 16` (`B = 44683`)
- Discriminants beyond 2 000 000 raise
- Duplicated binary-quadratic-form code across two modules

### Behavioral Impact

None on values — results are exact and match PARI on every discriminant tested; only the
reachable input range differs. `discriminants_with_bounded_class_number(8)` went from 6.26 s to
47 ms with identical, doctest-verified output.

---

## Coding Theory Deviations

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `GoppaCode` decoding | Registers **no decoder** (only a `GoppaCodeEncoder`); Sage falls back to generic syndrome decoding from `AbstractLinearCode` | Binary Goppa uses Patterson; non-binary uses the Sugiyama key equation via `_partialXGCDBalanced`, a direct port of `GRSKeyEquationSyndromeDecoder._partial_xgcd` (`grs_code.py:2145-2157`) followed by Forney `e_i = omega(L_i)/sigma'(L_i)`. Radius `floor(deg(g)/2)`, less than a generic syndrome decoder's true covering radius, but exact within it and failing loudly outside |
| `GoppaCode.distance_bound()` | `1 + deg(g)` | Same (was `2t + 1`) |
| `BCHCode.minimum_distance()` | Inherits `AbstractLinearCode.minimum_distance`, delegating to GAP/Guava's Brouwer-Zimmermann | Exhaustive codeword-weight enumeration, cached; `NotImplementedError` once `q^k > 2^17`. Exact wherever it answers (Golay `[23,12]` -> 7); previously it returned the *designed* distance, which is silently wrong |
| BCH field embedding | `finite_field_base.extension` uses `alpha = E.gen()^((\|E\|-1)/(\|F\|-1))` when both fields are Conway, else `self.modulus().any_root(E)` (randomized Cantor-Zassenhaus) | Tries the Conway power first and accepts it if it is a root of the base modulus; otherwise iterates over the splitting field looking for a root; `NotImplementedError` when `\|E\| > 2^22`. Deterministic, and exactly Sage's choice whenever Conway polynomials are tabulated (which is the port's default modulus) |
| Affected modules | `sage/coding/{goppa,bch,grs,reed_muller}_code.py` | `packages/sagemath-ts/src/coding/` |

### Rationale

1. **No upstream Goppa decoder to port** - so the algorithm was taken from the upstream decoder
   for the closest relative (the GRS key equation) rather than invented; the previous
   implementation merely toggled symbols between 0 and 1
2. **No Brouwer-Zimmermann port** - the choice was between an exact-but-limited algorithm and an
   honest stub; enumeration reproduces Sage's value where it answers and fails loudly otherwise
3. **No `any_root`/Cantor-Zassenhaus and no `is_conway` flag** - testing the Conway power
   directly is both deterministic and exactly Sage's choice in the tabulated case

### Trade-offs

- Non-binary Goppa corrects fewer errors than a generic syndrome decoder would
- Large BCH codes get `NotImplementedError` where Sage answers
- Very large non-Conway splitting fields make the BCH embedding unavailable rather than randomized

### Behavioral Impact

Sage's three `GRSKeyEquationSyndromeDecoder` doctests now reproduce exactly; the BCH generator
polynomial divides `x^n - 1` over GF(4)/GF(8)/GF(16); Forney carries the `X_i^(l-b)` factor so
`b ∈ {0,1,2,3}` and `l ∈ {1,5,7}` all decode; Reed-Muller's recursive Plotkin decoder decodes
`u` from both halves and keeps the closer candidate (0 failures over 1.3M decodes), and its
monomial order matches Sage's `Subsets` enumeration (`1, x0, x1, x2`). Every case tested took
the Conway branch of the BCH embedding and matched Sage's doctests.

---

## Crypto Module Deviations

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `gen_lattice` seeded output, `type='modular'` / `'random'` | `sage.crypto.gen_lattice(m=10, seed=42)` etc. print specific matrices | **Reproduced exactly.** Sage draws the random block with `MatrixSpace(ZZ_q, m-n, n).random_element()`, i.e. `Matrix_modn_dense_template.randomize` = `rstate.c_random() % p` row-major (`matrix_modn_dense_template.pxi:2843`); we now use `c_random()` rather than `mpz_urandomm`. The three doctests at `sage/crypto/lattice.py:81-105` and `:147-157` (modular, random, dual) are pinned in `lattice.test.ts` |
| `gen_lattice` seeded output, `type='ideal'` / `'cyclotomic'` | Prints a specific matrix | Does **not** reproduce it. Sage draws through `PolynomialRing.random_element(degree=n-1)` -> `IntegerModRing.random_element()`, a *different* consumption pattern from the modular path (leading coefficient first, with a redraw loop on zero). Our `PythonRandom` stream itself matches Sage's `python_random().randrange(q)` exactly, but the polynomial-ring layer above it is not ported, and the expected answer is Sage-version dependent (the vendored `reference/` doctest and SageMath 10.x disagree on this input). Structural invariants (block shape, circulant rows, `\|det\| = q^n` / `q^(m-n)`, minrep range, primal/dual relation) all match |
| `gen_lattice` quotient argument | Symbolic expression or polynomial | Coefficient list; monic required (see [Keyword Arguments…](#keyword-arguments-return-shapes-and-signature-adaptations)) |
| `IntegerLattice.gen_lattice` (`crypto/lattice.ts`) randomness | Sage's global randstate | A local seeded xorshift (pre-existing), so entries differ from Sage's doctest matrices even for the same seed |
| `LWE`/`RingLWE` `repr` | Prints `None` for an unbounded `m` | Prints `None` (was `null`). Every numeric field matches Sage's doctests — `LWE(20, 401, <sampler>, 'uniform', None)` with sigma 1.915069 and c 401 — and as of 0.0.12 the **sampler's own repr** matches too: `Discrete Gaussian sampler over the Integers with sigma = 1.915069 and c = 401.000000`, verified against SageMath 10.3 (the 0.0.11 invented form `DiscreteGaussianDistributionIntegerSampler(sigma=…, c=…, tau=…)` is gone) |
| `SBox` LAT | Per-mask Walsh-Hadamard transform | Same (was a naive `O(2^(2m+n))` scan; AES went 175 ms -> 6.5 ms). Values unchanged |
| Affected modules | `sage/crypto/{lattice,lwe,boolean_function,sbox}.py` | `packages/sagemath-ts/src/crypto/` |

### Rationale

1. **Two different generators upstream** - Sage's `modular`/`random` lattice block goes through
   GMP's MT19937 (`c_random() % q`) while the `ideal`/`cyclotomic` block goes through CPython's
   via the polynomial ring. Both generators are now ported bit-exactly (see
   [Random State and Seeding](#random-state-and-seeding)), so the modular/random doctests
   reproduce; what is still missing for `ideal`/`cyclotomic` is the *draw order* imposed by
   `PolynomialRing.random_element(degree=n-1)` (leading coefficient first, with a redraw loop on
   zero), which lives in a layer that is not ported
2. **The expected answer for `ideal`/`cyclotomic` is Sage-version dependent** — the vendored
   `reference/` doctest and SageMath 10.x disagree on that input — so pinning either would be
   pinning a version rather than a behaviour
3. **Verification by invariant** - those two branches are therefore verified with
   seed-independent structural and algebraic oracles (determinant, block identity, circulant
   rows, minrep range, primal/dual relation) rather than by pinning random values

### Trade-offs

- `gen_lattice(type='ideal'|'cyclotomic')` does not reproduce Sage's published matrix for a given
  seed; the other two types now do
- `crypto/lattice.ts`'s `IntegerLattice.gen_lattice` still uses a local seeded xorshift

### Behavioral Impact

`gen_lattice`'s three `modular` / `random` / `dual` doctests
(`sage/crypto/lattice.py:81-105`, `:147-157`) now reproduce **exactly**, verified both against
the vendored reference and against the live Sage install. Sage's doctest values also reproduce
exactly for `LindnerPeikert(20)` (q = 2053, sigma = 3.600954), `RingLindnerPeikert(16)`
(q = 1031, sigma = 2.803372), `Regev(20)` (q = 401, sigma = 1.915069, c = 401), the 3-round MISTY
construction (differential uniformity 8, linearity 64), the 8×8 LAT of `SBox(7,6,0,4,2,5,1,3)`
with all three scalings, `min_degree = 2`, the `'03'`/`'43'`/`'00ab'` hex round trips, the
algebraic-immunity cases and the dual lattice (six rows of 11, det = 11^6).

---

## ZK Sumcheck and Multilinear Extensions

These modules port `reference/sage_blueprints/`, not SageMath itself; there was previously no
`zk` entry in this document at all.

| Aspect | Blueprint | sagemath-ts |
|--------|-----------|-------------|
| Constant round polynomial | `sumcheck_round_prover` raises `ValueError('prover: Layer polynomial is not univariate')` when `len(res.variables()) != 1`, which includes a constant. Running the blueprint on the zero polynomial does not even reach that error — it crashes with `AttributeError: 'IntegerMod_int' object has no attribute 'variables'` | Only a round polynomial depending on a variable **other** than the free one is rejected (`'prover: Layer polynomial is not built from the correct variable'`); a constant round polynomial is returned as-is |
| `sumcheckVerify` signature | Round count from `len(poly.args())`; `degree_checks` defaults to `None` | `numVars` required, `{degreeCheck?}` optional and defaulting to no check — see [Keyword Arguments…](#keyword-arguments-return-shapes-and-signature-adaptations) |
| `interpolate_sparse([i])` | Short-circuits to the constant `R(i)`; that branch still contains leftover `print` debugging | Returns the selector `eq(i, x)`, consistent with the multi-index branch |
| `booleanHypercube(n)` | `Tuples([0,1], n)`, unbounded and lazy | Materializes the list; `ValueError` above `MAX_HYPERCUBE_DIM = 25`, and for non-integer `n` |
| Affected modules | `reference/sage_blueprints/{sumcheck,mle}.sage` | `packages/sagemath-ts/src/zk/` |

### Rationale

1. **The zero function is a legitimate sumcheck instance** - it has a zero round polynomial in
   every round and its sumcheck is sound; the blueprint's strictness is an artefact of Sage
   returning a base-ring element rather than a polynomial when all variables vanish
2. **Consistency** - the sparse-MLE short-circuit is inconsistent with every other input shape
   and is unreachable as an intentional API
3. **Memory honesty** - replacing `1 << n` with `2 ** n` fixes the 32-bit wraparound, but
   `n >= 26` then merely exhausts memory instead of silently returning a wrong-sized list

### Trade-offs

- A caller relying on the blueprint's constant behaviour for `[i]` gets a different polynomial
- A (spurious) error signal for constant round polynomials is lost
- Legitimate-if-impractical large `n` is refused

### Behavioral Impact

`sumcheckRun` on all-zero values succeeds here and crashes in the blueprint. The round
polynomial is now built **symbolically** (exact, arbitrary degree) rather than as a hardcoded
line through `p(0)`, `p(1)`, and variables are addressed by the ring's real names instead of
assumed `x0, x1, …` — so a ring named `[a,b]` gives `-2*a + 11`, and
`x0^2*x1 + x0 + 1` over GF(101) gives `x^2 + 2*x + 2`, both matching the blueprint run under a
real `sage` binary.

---

## ntl-ts GF2X Representation

| Aspect | NTL | ntl-ts |
|--------|-----|--------|
| `GF2X` storage | A `WordVector` (`xrep`) of machine words; `normalize()` strips zero words, `SetMaxLength(n)` preallocates, and the object can be temporarily unnormalized | A single bigint whose bit `i` is the coefficient of `x^i`, so it is **always** normalized: `normalize()` is a no-op and `SetMaxLength(n)` only performs NTL's negative-length check. `SetLength(n)` still truncates coefficients `>= n` exactly as NTL does |
| `XGCD` | Switches to a half-GCD recursion above `NTL_GF2X_GCD_CROSSOVER` (`GF2X1.cpp:3625`) | Plain extended Euclid |
| `PowerMod` | Sliding-window exponentiation (`GF2X1.cpp:1743`) | Binary square-and-multiply |
| `random`, `factor`, `SquareFreeDecomp`, `DistinctDegFactor`, `EqualDegFactor`, `BerlekampFactor`, `BuildRandomIrred` | Implemented | Honest `NTL_NOT_IMPLEMENTED` stubs — they need NTL's ChaCha-based `RandomStream`, `IrredPolyMod`/`GF2XModulus`, or randomized factoring. **Consequence:** `sagemath-ts`'s `polynomial_gf2x.ts` keeps four local factoring routines, and `irreducible_element(n, algorithm='random')` for p = 2 uses rejection sampling from `current_randstate()` instead of `BuildRandomIrred` — which is exactly the fallback SageMath itself takes when its NTL import fails (`polynomial_ring.py:3615-3620`) |
| Affected modules | `ntl/src/GF2X*.cpp`, `GF2XFactoring.cpp` | `packages/ntl-ts/src/GF2X.ts`, `GF2.ts`, `GF2X_irred_tab.ts` |

### Rationale

1. **bigint already is an arbitrary-precision bit vector** with XOR and shifts; hand-rolling a
   word vector would add no fidelity and duplicate what the runtime does
2. **Identical results** - over GF(2) the gcd is unique, and the Bezout pair with
   `deg(s) < deg(b) - deg(d)`, `deg(t) < deg(a) - deg(d)` that NTL returns is exactly the one
   plain extended Euclid produces; sliding-window and square-and-multiply compute the same power

### Trade-offs

- `XGCD` is `O(n^2)` bit operations rather than `O(n log^2 n)` — irrelevant at the degrees
  Sage's default-modulus path uses (n up to a few thousand)
- Factoring and random generation are unavailable

### Behavioral Impact

None on any implemented function. `IterIrredTest`, `BuildIrred` and `BuildSparseIrred` are
line-for-line ports of `GF2XFactoring.cpp` over a vendored copy of NTL's 2049-row
`GF2X_irred_tab`, verified against Sage's `polynomial_gf2x.pyx` doctests
(`BuildIrred_list(2/3/4/33)`, `BuildSparseIrred == BuildIrred` for `n ∈ [1,32]`,
`BuildSparseIrred(33) = x^33 + x^10 + 1`) and against exhaustive brute-force irreducibility over
all 2046 monic polynomials of degree <= 10.

As of 0.0.12, `sagemath-ts`'s `polynomial_gf2x.ts` **delegates** its `add`/`sub`/`neg`/`mul`/
`sqr`/`leftShift`/`rightShift`/`trunc`/`divRem`/`gcd`/`xgcd`/`powMod`/`derivative`/`reverse`/
`is_irreducible`/`buildIrred`/`buildSparseIrred` to this package instead of reimplementing NTL.
Old-vs-new values are identical for every `n ∈ [2,160]`, so nothing regressed, and
`buildSparseIrred` now uses NTL's real precomputed minimal-weight table rather than an ad-hoc
trinomial/pentanomial scan.

---

## Matrix Special Constructors

`sage/matrix/special.py` maps to `packages/sagemath-ts/src/matrix/matrix_special.ts`
(3662 lines). Individual divergences were previously scattered across four other sections of this
document with no entry naming the module; audit item **L44** asked for a consolidated one. This
is it.

| Function | SageMath | sagemath-ts |
|----------|----------|-------------|
| `companion_matrix`, `toeplitz`, `hankel` | See `special.py` | Argument conventions corrected to Sage's in 0.0.11 (full monic coefficient list with negated border; `r` counted from the second column with `ncols = len(r)+1`). The port's own tests had pinned the *wrong* conventions |
| `elementary_matrix(row1 == row2)` with no scale | `special.py:1512-1516`'s four assignments collapse to `elem[r,r] = 1`, i.e. the identity; Sage raises only when a scale is *also* given | Replicated verbatim — identity for a self-swap, `ValueError` for the two cases Sage rejects. The audit claimed Sage raises in general; the vendored source shows otherwise |
| `block_matrix` | Flat list + `nrows`/`ncols`, ragged list, or list of lists | List of lists only; a ragged one raises Sage's own `ValueError('list of rows is not valid (rows are wrong types or lengths)')` |
| `random_echelonizable_matrix` / `random_unimodular_matrix` `upper_bound` | Size control by repeatedly rejecting row operations that push an entry past the bound — **only over ZZ and QQ** | `NotImplementedError`. The port's generic constructors work over any ring with `random_element()` and have **no notion of absolute value**, so the option is refused rather than silently ignored |
| `random_unitary_matrix`, `vector_on_axis_rotation_matrix`, `ith_to_zero_rotation_matrix` | Implemented over RDF/CDF via QR / Haar measure / trigonometric rotations | `NotImplementedError` naming the requirement (`sqrt` and trigonometric functions over an inexact ring). CLAUDE.md forbids introducing floating point for exact results, and there is no RDF/CDF matrix type here |
| `hadamard_bound` | Uses `sqrt` in the base ring | `NotImplementedError` for rings without `sqrt` (i.e. everything but RDF/`RealField`) |
| `lehmer`, `hilbert` | Return matrices over QQ | Require a ring argument supporting `__call__` and division; `NotImplementedError` otherwise |
| `rook_vector` | `ButeraPernici` (default), `Ryser` and `Godsil` algorithms | Naive placement counting only. Throws `NotImplementedError` naming the two faster algorithms once `positions.length > 50` **and** `k > 5`, rather than enumerating `C(p, k)` |
| `is_permutation_of` / `permutation_normal_form` | `BipartiteGraph.is_isomorphic(..., edge_labels=True)` (bliss/nauty) | Complete backtracking pruned by column-multiset and row/column signature invariants. Both outputs are *uniquely specified*, so any complete algorithm agrees; both were checked against exhaustive brute force. Worst case exponential, as is Sage's. `permutation_normal_form(check=true)` returns 0-based index arrays with the convention `normal_form[i][j] === matrix[row_perm[i]][col_perm[j]]` instead of a pair of 1-based `PermutationGroupElement`s |
| Random matrix constructors | `sage.misc.prandom` (`randint`/`shuffle`) driven by the global randstate | `current_randstate().randint(...)` and a randstate-driven Fisher-Yates shuffle / density fraction (was `Math.random()`). As of 0.0.12 the underlying stream is bit-identical to GMP's, so `set_random_seed` now gives Sage's values wherever the *draw order* also matches |
| `random_diagonalizable_matrix` | — | Raises `NotImplementedError('unexpected eigenvector layout')` on an internal invariant violation rather than emitting a matrix that is not diagonalizable |
| Affected modules | `sage/matrix/special.py` | `packages/sagemath-ts/src/matrix/matrix_special.ts` |

### Rationale

1. **No inexact matrix type** - `random_unitary_matrix` and the two rotation constructors are
   defined over RDF/CDF. Implementing them with JS doubles inside a `Matrix<R>` would introduce
   floating point into the exact matrix hierarchy for functions nothing in this port consumes.
2. **Refusing an argument beats ignoring it** - `upper_bound` is the clearest case: SageMath's
   own implementation is ZZ/QQ-only, and the port's constructors are ring-generic.
3. **Bounded naive algorithms are declared** - `rook_vector`'s cap is explicit so the reachable
   range is documented rather than discovered as a hang. See
   [Bounded Search Budgets](#bounded-search-budgets-and-measured-thresholds).

### Trade-offs

- Five constructors (`random_unitary_matrix`, `vector_on_axis_rotation_matrix`,
  `ith_to_zero_rotation_matrix`, `hadamard_bound` off RDF, `rook_vector` at scale) are
  unavailable where SageMath answers
- `block_matrix`'s flat-list and ragged forms are unavailable
- `permutation_normal_form(check=true)` may return a different (equally valid) permutation when
  the matrix has non-trivial automorphisms; the matrix returned is identical

### Behavioral Impact

98 tests cover the module. Values match SageMath's for every implemented constructor. The
divergences above are either honest refusals or index-base/shape adaptations — none returns a
different mathematical object.

---

## Bounded Search Budgets and Measured Thresholds

Several upstream algorithms are unbounded loops backed by a fallback this port does not have, or
are asymptotically better than what is reachable here. Rather than hang or silently degrade, the
port imposes an explicit budget and fails with a message naming what is missing. Collected here
because they share one rationale and are easy to mistake for arbitrary magic numbers.

| Site | Upstream | Budget here | What happens at the limit |
|------|----------|-------------|---------------------------|
| Polynomial factorization recombination | van Hoeij/LLL on the knapsack lattice | 200 000 subsets | `NotImplementedError` naming van Hoeij/LLL. Nothing constructible reaches it: Swinnerton-Dyer degree 32 (16 quadratic modular factors) finishes in 324 ms |
| Factorization prime search | Unbounded | Primes below 10 000 | `NotImplementedError` rather than guessing. Raising the bound is trivial; the honest failure was left deliberately |
| PARI insisting-ECM | Loops **forever**, because MPQS backs it up | `FactorOptions.ecmRounds`, default 4 (MPQS backs it up here too, as in PARI) | Bounded work per insisting round; without the bound `Z_factor` could hang indefinitely |
| `ellcard` Schoof/Shanks crossover | PARI switches to SEA at `expi(p) >= 56` | Schoof from `expi(p) >= 96` | Not a failure — a **measured** threshold. Copying 56 would make `ellcard` hundreds of times slower over 2^56..2^88 with only base Schoof available. Timings are in the source comment and in [PARI Elliptic Curve](#pari-elliptic-curve-advanced-algorithms-missing-parigp-ts) |
| Class group of a degree > 2 field | `bnfinit` (subexponential) | Minkowski bound `<= 10^6`, and only the provably-trivial case | `NotImplementedError` naming `bnfinit` and its missing components. **Fields whose true class number is 1 still throw** if we cannot prove it |
| Quadratic class group | `quadclassunit` (subexponential) | `CLASS_GROUP_DISC_BOUND = 2 000 000` | `NotImplementedError` |
| `nfgaloisconj` precision escalation | PARI escalates inside `galoisgen` | 12 attempts, squaring the p-adic precision each time | `NotImplementedError`. The Gram-Schmidt certificate means an inconclusive result is *known* to be inconclusive, never a wrong answer |
| `is_similar` intertwining search | Sage raises `RuntimeError` instead | 200 pseudo-random kernel combinations | `ArithmeticError` with Sage's message text. Unreachable in testing: the proportion of units in a centralizer algebra over `F_q` is at least `prod(1 - q^-i) >= 0.288` |
| `rook_vector` naive algorithm | `ButeraPernici` / `Ryser` | 50 positions and `k > 5` | `NotImplementedError` naming the two faster algorithms |
| `voronoiCell` / exact SVP / BCH minimum distance / BCH field embedding | Backend-accelerated | rank 24 / rank 30 / `q^k > 2^17` / `\|E\| > 2^22` | `NotImplementedError` (pre-existing, 0.0.11) |

### Rationale

1. **A hang is the worst failure mode** in a library, worse than an exception: it gives the
   caller nothing to act on and no signal that a dependency is missing.
2. **The budget documents the reachable range.** Where the upstream loop is unbounded only
   because a *later* stage catches the hard cases (MPQS behind ECM; van Hoeij behind Zassenhaus),
   removing that stage without adding a bound converts "slow" into "never returns".
3. **Thresholds must be measured, not copied.** PARI's 56-bit SEA crossover is correct *for
   PARI*; transplanting it into a port that only has base Schoof would be a fidelity gesture that
   makes the function unusable. The number that matters is where our two implementations actually
   cross.

### Trade-offs

- Inputs SageMath handles can fail here, and the failure is a hard error rather than a long wait
- The specific constants are tuning choices, not upstream values, and would need re-measuring on
  different hardware or after asymptotic improvements to the underlying arithmetic
- A caller who *wants* to spend more time has an escape hatch only where one was added
  (`FactorOptions.ecmRounds`, `solve_integer`'s `factorization`, `qfbsolve`'s `fa`)

### Behavioral Impact

Every budget exhaustion raises with a message naming the missing upstream routine, so
`grep -r "SAGE_NOT_IMPLEMENTED" packages/` remains the authoritative list of gaps. No budget
silently truncates a result.

---

## Template for New Deviations

Copy this template when adding a new deviation:

```markdown
## [Deviation Title]

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Description | ... | ... |
| Affected modules | ... | ... |

### Rationale

1. **Reason 1** - Explanation
2. **Reason 2** - Explanation

### Trade-offs

- What we lose by deviating

### Mitigation

How we could address the trade-offs if needed.

### Behavioral Impact

Does this change outputs? Edge cases? Error messages?
```
