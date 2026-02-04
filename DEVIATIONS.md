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
25. [Finite Field Constructors and Display](#finite-field-constructors-and-display)
26. [Conway Polynomial Database Limited](#conway-polynomial-database-limited)
27. [Finite Field Extension Minimal Polynomial Simplified](#finite-field-extension-minimal-polynomial-simplified)
28. [Error Class Parity](#error-class-parity)
29. [Caching and Import Paths](#caching-and-import-paths)
30. [Elliptic Curve Short Weierstrass Form Only](#elliptic-curve-short-weierstrass-form-only)
31. [Real/Complex Numerical Approximations](#realcomplex-numerical-approximations)
32. [Matrix and Lattice Algorithm Simplifications](#matrix-and-lattice-algorithm-simplifications)
33. [Shortest Vector Problem (SVP) Implementation](#shortest-vector-problem-svp-implementation)
34. [Closest Vector Approximation](#closest-vector-approximation)
35. [Algebraic Dependency Approximation](#algebraic-dependency-approximation)
36. [Combinatorial Function Limits](#combinatorial-function-limits)
37. [Real Matrix Decompositions (SVD, QR, LU) Using IEEE 754 Doubles](#real-matrix-decompositions-svd-qr-lu-using-ieee-754-doubles)
38. [PARI Factorization Algorithms Limited (parigp-ts)](#pari-factorization-algorithms-limited-parigp-ts)
39. [PARI Elliptic Curve Advanced Algorithms Missing (parigp-ts)](#pari-elliptic-curve-advanced-algorithms-missing-parigp-ts)
40. [Template for New Deviations](#template-for-new-deviations)

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
| RNG core | `sage.misc.randstate` (GMP randstate + Python `random`) | `RandState` with 64-bit LCG |
| Seeding | `set_random_seed()` affects all random consumers | `set_random_seed()` affects all random consumers |
| Affected modules | Any `.random_element()` / random sampling | Same |

### Rationale

1. **Portability** - Works in any JS runtime without GMP bindings
2. **Centralization** - Matches Sage’s single randstate model
3. **Determinism** - Seeded runs are reproducible within sagemath-ts

### Trade-offs

- RNG sequences differ from Sage for the same seed
- LCG is not cryptographically secure

### Mitigation

If exact Sage random sequences are needed, replace `RandState` with an MT19937
or GMP-compatible generator and keep the same API.

### Behavioral Impact

Random outputs differ from Sage for identical seeds, but distribution and
range semantics match Sage.

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

### Not Yet Implemented (Require PARI)

The following operations throw `NotImplementedError`:

- **Ring of integers / maximal order** - Requires integral basis computation
- **Class group and class number** - Requires PARI bnfinit
- **Unit group** - Requires PARI bnfinit for fundamental units
- **Regulator** - Requires numerical computation of logarithms
- **Prime ideal factorization** - Requires Dedekind-Kummer theorem
- **Galois group** - Requires root finding in splitting field
- **Ideal operations** - HNF, LLL for ideal arithmetic

### Trade-offs

- Class group and unit group computations are not available
- Some operations that SageMath computes exactly may be unavailable
- No proof flags for GRH-conditional results

### Mitigation

When parigp-ts adds nfinit/bnfinit support, we can:
1. Update NumberField to call parigp-ts for heavy computations
2. Keep the TypeScript fallback for basic operations
3. Match SageMath's architecture exactly

### Behavioral Impact

- Basic element operations (norm, trace, arithmetic) match SageMath exactly
- Advanced operations (class group, unit group) throw NotImplementedError instead of computing
- Signature computation may have edge cases for higher degree fields

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
- **Square root** - sqrt() using Tonelli-Shanks + Hensel lifting
- **Teichmuller lift** - teichmuller() using Newton iteration
- **Log/Exp** - log() for units, exp() for convergent inputs
- **Norm/Trace** - Trivial for base field (identity maps)

### Not Yet Implemented

- **nth_root() for n != 2** - Only square root is implemented
- **log() with branch** - Branch parameter for non-units not supported
- **artin_hasse_exp()** - Artin-Hasse exponential
- **Extension fields** - Zq, ramified/unramified extensions
- **minimal_polynomial(), charpoly()** - Require extension field support

### Trade-offs

- Only capped-relative precision model is available
- SageMath uses PARI for gamma(), sqrt() of extension fields, etc.
- Some precision handling may differ in edge cases

### Mitigation

When parigp-ts adds padic support, we can delegate operations like:
- gamma() - p-adic gamma function
- sqrt() for extension fields
- Extension field arithmetic

### Behavioral Impact

- Basic arithmetic matches SageMath for capped-relative elements
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

### Not Yet Implemented (throw NotImplementedError)

**padic_lseries.ts:**
- `modular_symbol(r, sign, quadratic_twist)` - requires sage.modular.modsym
- `measure(a, n, prec, quadratic_twist, sign)` - requires modular symbols + p-adic arithmetic
- `alpha(prec)` - requires p-adic fields (Q_p root finding)
- `teichmuller(prec)` - requires p-adic fields
- `series(n, quadratic_twist, prec, eta)` - requires modular symbols + power series
- `order_of_vanishing()` - requires series computation
- `frobenius(prec, algorithm)` - requires Monsky-Washnitzer cohomology
- `bernardi_sigma_function(prec)` - requires formal group computation

**isogeny_class.ts:**
- `graph()` - requires graph data structure
- `isogenies(fill=true)` - requires full isogeny computation
- `qf_matrix()` - only for CM curves, requires CM discriminant
- `isogeny_degrees_cm(E)` - requires CM and class group computation
- `possible_isogeny_degrees(E)` - requires Galois representation analysis

### Trade-offs

- Cannot compute actual p-adic L-function values
- Cannot compute full isogeny classes (only trivial class with initial curve)
- Users get clear NotImplementedError with explanation

### Mitigation

When the following dependencies are available, full implementation becomes possible:
1. `sage.modular.modsym` port - enables modular symbol computation
2. `sage.rings.padics.Qp` completion - enables p-adic arithmetic
3. `sage.rings.power_series_ring` - enables power series
4. Isogeny computation from `ell_curve_isogeny.ts` - enables full class computation

### Behavioral Impact

- Basic class structure and getters work identically to SageMath
- Advanced methods throw `NotImplementedError` with descriptive messages
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
| Kernel polynomial algorithms | Kohel/BMSS/Stark for broad cases | BMSS uses small-field enumeration only; Kohel kernel polynomial not implemented |
| Rational maps | Multiple algorithms | Rational maps only for Velu isogenies |
| Dual/inseparable handling | Full support | Dual for inseparable or char 2/3 not implemented |
| Isogenous check over finite fields | Uses cardinality via PARI/SEA | Small fields only (enumeration); large fields throw |
| Affected modules | `sage/schemes/elliptic_curves/ell_curve_isogeny.py` | `packages/sagemath-ts/src/schemes/elliptic_curves/ell_curve_isogeny.ts` |

### Rationale

1. **Missing infrastructure** - Kohel kernel polynomials and SEA-based cardinality are not yet available
2. **Algorithmic scope** - Implemented Velu and small-field enumeration first

### Trade-offs

- Large-field isogeny workflows are incomplete
- Some advanced isogeny operations raise `NotImplementedError`

### Mitigation

Implement Kohel/Stark algorithms and integrate PARI-based cardinality/SEA for large fields.

### Behavioral Impact

Isogeny computations are limited to small fields and Velu-based paths; some operations throw for larger fields or special characteristics.

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

## Finite Field Constructors and Display

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `GF(q)` extension fields | `GF(p^n)` works in one constructor | `GF()` only supports prime fields |
| Extension field constructor | `GF(p^n, 'a')` | `GFExtended()` / `GFpn()` |
| Element display | `repr(a)` and `a.lift()` | `a.toString()`, `a.repr()`, `a.value` |
| Affected modules | `rings/finite_rings/*` | Same |

### Rationale

1. **API clarity** - Separate constructors improve type inference
2. **Implementation staging** - Extension fields are handled in dedicated modules
3. **Interop** - `.value` provides direct `bigint` access

### Trade-offs

- One-call `GF(p^n)` convenience is limited
- Users must choose the right constructor for extensions

### Behavioral Impact

Prime fields use `GF()`, extension fields require explicit constructors.

---

## Conway Polynomial Database Limited

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| Conway polynomial availability | Large curated database | Limited to a small subset of primes/degrees |
| Default modulus for GF(p^n) | Conway polynomial when available | Conway when available, otherwise random irreducible |
| Affected modules | `sage/rings/finite_rings/conway_polynomials.py` | `packages/sagemath-ts/src/rings/finite_rings/conway_polynomials.ts`, `finite_field_extension.ts` |

### Rationale

1. **Scope** - Only a small subset of Conway polynomials is bundled
2. **Pragmatism** - Random irreducible polynomials are sufficient for many uses

### Trade-offs

- Field generators may differ from SageMath’s canonical choices
- Embeddings between fields may not follow Conway compatibility

### Mitigation

Expand the Conway polynomial database or provide a fetchable optional dataset.

### Behavioral Impact

For primes/degrees not in the local database, GF(p^n) construction may choose a different defining polynomial than SageMath.

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
| LLL | fplll/NTL exact/MPFR variants | Simplified float-based LLL |
| BKZ | Full enumeration/pruning | Repeated LLL passes (approximate) |
| Frobenius form | PARI `matfrobenius` | Basis matrix is identity placeholder |
| `p_minimal_polynomials`, `null_ideal`, `integer_valued_polynomials_generators` | `compute_J_ideal` algorithms | Simplified outputs (charpoly/minpoly-based) |
| Affected modules | `matrix/matrix_integer.ts` | Same |

### Rationale

1. **Dependency gap** - No fplll/NTL/pari bindings in TypeScript
2. **Incremental implementation** - Provide structure and basic outputs first
3. **Complexity** - Full algorithms are substantial

### Trade-offs

- Reduction quality and correctness can differ
- Some outputs are placeholders, not canonical

### Behavioral Impact

Results may be approximate or structurally incomplete compared to SageMath.

---

## Shortest Vector Problem (SVP) Implementation

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `IntegerLattice.shortest_vector` | fpylll SVP.shortest_vector or PARI qfminim | Schnorr-Euchner enumeration in pure TypeScript |
| Affected modules | `modules/free_module_integer.ts` | Same |

### Rationale

1. **No external dependencies** - Pure TypeScript implementation without fpylll or PARI bindings
2. **Algorithm parity** - Schnorr-Euchner is the same algorithm fpylll uses internally
3. **Portability** - Works in any JavaScript runtime without WebAssembly or native bindings

### Implementation Details

- **Dimension <= 4**: Simple exhaustive enumeration
- **Dimension 5-50**: Schnorr-Euchner enumeration with zig-zag coefficient pattern and branch-and-bound pruning
- **Dimension > 50**: Returns first LLL-reduced basis vector (2^((n-1)/2) approximation)

### Trade-offs

- Performance is slower than fpylll's optimized C++ implementation
- No pruning strategies from advanced BKZ variants
- Large lattices (dim > 50) use approximation instead of exact SVP

### Mitigation

For performance-critical applications with large lattices, users can:
1. Use LLL approximation directly via `lattice.LLL()` first row
2. Use BKZ for better approximation quality
3. Interface with external SVP solvers if needed

### Behavioral Impact

- For small lattices (dim <= 50), results are exact shortest vectors (same as SageMath)
- For large lattices (dim > 50), returns an approximation (first LLL basis vector)
- SageMath with algorithm='pari' uses qfminim which may enumerate differently

---

## Closest Vector Approximation

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `IntegerLattice.closest_vector` | Voronoi-cell (exact CVP) | Babai nearest-plane for rank > 4 |
| Affected modules | `modules/free_module_integer.ts` | Same |

### Rationale

1. **Complexity** - Exact CVP is expensive and involved
2. **Practicality** - Babai is fast and often sufficient

### Trade-offs

- May return non-optimal vectors for higher dimensions

### Behavioral Impact

For rank > 4 lattices, results can differ from SageMath's exact CVP.

---

## Algebraic Dependency Approximation

| Aspect | SageMath | sagemath-ts |
|--------|----------|-------------|
| `algebraic_dependency` | PARI `algdep`, supports real/complex/p-adic with proof/height bounds | Simplified LLL on IEEE doubles, real-only |
| Affected modules | `arith/misc.ts` | Same |

### Rationale

1. **No PARI binding** - `algdep` not available yet
2. **Simple fallback** - LLL on doubles gives practical results for small degrees

### Trade-offs

- No complex or p-adic support
- Proof/height-bound semantics are not equivalent to SageMath

### Behavioral Impact

Results can be approximate or incorrect for subtle inputs; API accepts options but does not fully match SageMath guarantees.

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

| Aspect | SageMath (PARI/GP) | sagemath-ts (parigp-ts) |
|--------|--------------------|-------------------------|
| Integer factorization (`Z_factor`) | Trial division + Pollard rho + ECM + MPQS | Trial division + BPSW only; composites beyond the bound may remain unfactored |
| Affected modules | `pari/src/basemath/ifactor1.c` | `packages/parigp-ts/src/ifactor.ts` |

### Rationale

1. **Dependency gap** - ECM/MPQS implementations are not yet available in TypeScript
2. **Performance constraints** - Pure JS bigint arithmetic makes advanced algorithms significantly slower
3. **Incremental porting** - Prioritize correctness for small/medium inputs first

### Trade-offs

- Large composite inputs can remain partially factored
- Factorization performance diverges from PARI for cryptographic-size integers

### Mitigation

Implement Pollard rho/ECM/MPQS in parigp-ts or delegate to a native/WASM backend once available.

### Behavioral Impact

For large composites, `Z_factor` may return a composite as a single factor and emit a warning where PARI would fully factor.

---

## PARI Elliptic Curve Advanced Algorithms Missing (parigp-ts)

| Aspect | SageMath (PARI/GP) | sagemath-ts (parigp-ts) |
|--------|--------------------|-------------------------|
| Point counting for large primes (`ellcard`) | SEA (Schoof-Elkies-Atkin) for large p | BSGS for all p; large primes can be impractical |
| Advanced functions (`ellcard_sea`, `ellisogeny*`, `ellfrobenius`) | Fully implemented | Stubs that throw `PARI_NOT_IMPLEMENTED` |
| Affected modules | `pari/src/basemath/ellsea.c`, `pari/src/basemath/ellisog.c` | `packages/parigp-ts/src/elliptic/group.ts`, `packages/parigp-ts/src/elliptic/advanced.ts` |

### Rationale

1. **Missing modular polynomial infrastructure** - SEA and isogenies depend on heavy polynomial arithmetic
2. **Complexity** - Full implementations are substantial and not yet ported
3. **Prioritization** - Focus on BSGS-based functionality first

### Trade-offs

- Large prime point counts can be too slow to compute in practice
- Isogeny and Frobenius functionality is unavailable

### Mitigation

Add modular polynomial support and implement SEA/isogeny algorithms (or delegate to a native/WASM backend).

### Behavioral Impact

For large primes, `ellcard` may be infeasible. Calls to `ellcard_sea`, `ellisogeny`, `ellisogenyapply`, `ellisogenycompose`, and `ellfrobenius` throw errors.

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
