# Design Decisions

This document explains the architectural decisions made when porting SageMath to TypeScript. It covers type mappings, function signatures, error handling patterns, and dependency architecture.

> **Related Documents:**
> - **DEVIATIONS.md** - Behavioral differences from SageMath (when outputs differ)
> - **SCOPE.md** - Module implementation status and assignments
> - **AGENTS.md** - Contributor workflow and testing guidelines

## Table of Contents

1. [Type System Mapping](#type-system-mapping)
2. [Function Signatures](#function-signatures)
3. [Error Handling](#error-handling)
4. [Module Structure](#module-structure)
5. [Ring Coercion](#ring-coercion)
6. [Dependency Architecture](#dependency-architecture)
7. [TypeScript Limitations](#typescript-limitations)

---

## Type System Mapping

### Integer Types

SageMath has two integer types that work interchangeably:

| Python/SageMath | TypeScript/sagemath-ts | Purpose |
|-----------------|------------------------|---------|
| `int` (primitive) | `bigint` (primitive) | Native arbitrary precision |
| `Integer` (class) | `Integer` (class) | Rich wrapper with methods |

### The `IntegerLike` Type

To match SageMath's flexibility where functions accept both `int` and `Integer`, we define:

```typescript
// src/types/coercion.ts
export type IntegerLike = bigint | Integer;
```

**Why NOT include JavaScript `number`?**
- JavaScript `number` is IEEE 754 double-precision float, which silently loses precision for integers > 2^53-1
- Source code literals like `9007199254740993` already lose precision before reaching our code (becomes `9007199254740992`)
- This library targets cryptographic applications where silent precision loss could cause security vulnerabilities
- The `n` suffix for bigint literals is a small inconvenience compared to data corruption risk

See **DEVIATIONS.md** section "Language and Type-System Adaptations" for the full rationale,
including the three APIs that still accept a raw `number` and should be widened.

### The `toBigInt()` Coercion Function

All functions normalize inputs immediately:

```typescript
import { IntegerLike, toBigInt } from '../types/coercion.js';

export function gcd(a: IntegerLike, b: IntegerLike): bigint {
  const _a = toBigInt(a);  // Normalize to bigint
  const _b = toBigInt(b);
  // ... implementation using _a and _b
}
```

The `toBigInt()` function:
- Returns `bigint` unchanged
- Extracts `.value` from `Integer` instances
- **Throws `TypeError` if given a JavaScript `number`**

### The `toSafeNumber()` Conversion Function

When internal code must convert bigint to number (e.g., for array indices or floating-point math):

```typescript
import { toSafeNumber } from '../types/coercion.js';

// Safe: throws RangeError if value exceeds ±2^53-1
const idx = toSafeNumber(bigintValue);

// Unsafe: silently loses precision - avoid!
const idx = Number(bigintValue);
```

### Design Rationale

| Decision | Rationale |
|----------|-----------|
| `bigint` as return type | Simpler, more efficient than wrapping in `Integer` |
| Accept `IntegerLike` inputs | Matches SageMath's flexibility |
| **Exclude `number`** | Prevent silent precision loss in crypto applications |
| Coerce immediately | Clear semantics, single code path |

### Examples

```typescript
// All of these work:
gcd(12n, 8n)                              // bigint literals
gcd(new Integer(12n), new Integer(8n))    // Integer objects
gcd(12n, new Integer(8n))                 // Mixed types

// This throws TypeError (numbers not accepted):
gcd(12, 8)  // "JavaScript numbers are not accepted due to precision loss risk"
```

---

## Function Signatures

### Naming Convention

**Match SageMath function names exactly.** Even when TypeScript conventions differ:

```typescript
// SageMath: def is_prime(n):
// We write:
function is_prime(n: IntegerLike): boolean  // NOT isPrime

// SageMath: def nth_prime(n):
// We write:
function nth_prime(n: IntegerLike): bigint  // NOT nthPrime
```

### Parameter Types

| Parameter accepts... | Use type |
|---------------------|----------|
| Any integer | `IntegerLike` |
| Must be bigint (internal) | `bigint` |
| Specific class instance | `Integer`, `Polynomial`, etc. |

### Return Types

| Situation | Return type |
|-----------|-------------|
| Computed integer value | `bigint` (not `Integer`) |
| May not exist | `bigint | null` |
| Class instance | The class type |
| Boolean predicate | `boolean` |

**Why return `bigint` instead of `Integer`?**
- Simpler for consumers (no unwrapping needed)
- More efficient (no object allocation)
- Easy to wrap if needed: `new Integer(result)` (or `ZZ.__call__(result)`, which returns a
  bare `bigint` — see [Ring Coercion](#ring-coercion))

### Scalar Constructors

`Integer` and `ZZ.__call__` share the source-defined coercion path; they accept no value
(zero), IntegerLike, integral numbers, integral rationals, booleans, null, strings with an
optional IntegerLike base, digit arrays with a base, and structural integer conversion hooks.
This does not widen `toBigInt` or ordinary IntegerLike arithmetic to floating-point inputs.
An unrecognized plain JavaScript object (including a null-prototype record) is the
bare Python `object()` stand-in for scalar conversion diagnostics: its Integer
rejection uses `unable to coerce <class 'object'> to an integer`. Recognized hooks
and polynomial/fraction objects retain their existing conversion paths.

Rational's pair constructor defaults to zero. Its `from` factory accepts Integer wrappers,
booleans, null and an omitted value. QQ tuple components delegate to Integer coercion, while
constructor-level None is handled separately from nested singleton lists. See the scalar
constructor entry in DEVIATIONS.md for the required tuple and conversion-map adaptations.

### Rational Coercion and Extended Results

`Rational.from(number)` and `QQ.__call__(number)` follow Sage's real-to-rational coercion,
using the existing real field's `simplest_rational`. They are distinct from strict integer
coercion through `toBigInt`. Rational strings contain an integer or a quotient in GMP's
base-zero syntax; decimal values enter through numeric coercion.

Single-argument `[numerator, denominator]` tuples require integer components as in Sage.
The retained two-argument `QQ.__call__(numerator, denominator)` convenience computes a quotient
of rationals; it is explicitly distinct from Sage's second argument (a string base).
`QQ.random_element` delegates all draws to `ZZ.random_element`, including unbounded defaults
and zero-denominator retries, so seeded streams remain aligned.

`Rational.toNumber()` is the shared nearest-even binary64 conversion, including for rational
matrix entries. `ord` and zero numerator valuations use the existing `'Infinity'` sentinel;
`gamma` at integral poles uses `UnsignedInfinityElement`. Consumers must handle the declared
unions. Approximate height precision limits and unimplemented symbolic branches are recorded
in DEVIATIONS.md.

### Options Objects

Replace Python's keyword arguments with TypeScript options objects:

```python
# SageMath
def some_function(n, algorithm='pari', proof=None, limit=None):
    ...
```

```typescript
// TypeScript
interface SomeOptions {
  algorithm?: 'pari' | 'flint';
  proof?: boolean;
  limit?: number;
}

function some_function(n: IntegerLike, options?: SomeOptions): Result {
  const { algorithm = 'pari', proof, limit } = options ?? {};
  // ...
}
```

> These names are illustrative. In particular our `factor(n)` takes **no** options — it
> delegates unconditionally to `parigp-ts`'s `Z_factor`. Check the real signature in the
> source before writing a call.

### Overloads

When SageMath has different behavior based on argument types, use TypeScript overloads:

```typescript
// Different return type based on input — this is the real `sqrt_mod` shape
function sqrt_mod(a: IntegerLike, p: IntegerLike, all_roots?: false): bigint | null;
function sqrt_mod(a: IntegerLike, p: IntegerLike, all_roots: true): bigint[];
function sqrt_mod(a: IntegerLike, p: IntegerLike, all_roots = false): bigint | null | bigint[] {
  // Implementation
}
```

Keep the flag's SageMath name (`all_roots`, `get_data`, `proof`), not a shortened one.

---

## Error Handling

### Exception Mapping

Map Python exceptions to TypeScript custom error classes:

| Python Exception | TypeScript Class | Location |
|-----------------|------------------|----------|
| `ValueError` | `ValueError` | `src/errors.ts` |
| `TypeError` | `TypeError` (native) | Built-in |
| `NotImplementedError` | `NotImplementedError` | `src/errors.ts` |
| `ZeroDivisionError` | `ZeroDivisionError` | `src/errors.ts` |
| `ArithmeticError` | `ArithmeticError` | `src/errors.ts` |

### Error Messages

**Preserve SageMath's error messages** when possible:

```typescript
// SageMath: raise ValueError("n must be positive")
throw new ValueError("n must be positive");

// SageMath: raise ZeroDivisionError("rational division by zero")
throw new ZeroDivisionError("rational division by zero");
```

### Null vs Exceptions

In some cases, we deviate from SageMath by returning `null` instead of throwing:

| Function | SageMath | sagemath-ts |
|----------|----------|-------------|
| `sqrt_mod()` for non-residue | `ValueError` | `null` |
| `discrete_log()` not found | `ValueError` | `null` |

**Rationale:** TypeScript idioms favor returning `null` for "not found" cases with union types (`bigint | null`) to force callers to handle failure explicitly.

Document these in **DEVIATIONS.md** when they affect behavior.

### Unimplemented Functions

Stub all functions with `NotImplementedError`:

```typescript
export function unimplemented_function(n: IntegerLike): bigint {
  throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: unimplemented_function');
}
```

The `SAGE_NOT_IMPLEMENTED:` prefix enables discovery:
```bash
grep -r "SAGE_NOT_IMPLEMENTED" packages/
```

---

## Module Structure

### Directory Mapping

Mirror SageMath's directory structure:

| SageMath Path | TypeScript Path |
|---------------|-----------------|
| `sage/rings/integer.py` | `src/rings/integer.ts` |
| `sage/arith/misc.py` | `src/arith/misc.ts` |
| `sage/rings/finite_rings/finite_field_constructor.py` | `src/rings/finite_rings/finite_field_constructor.ts` |
| `sage/schemes/elliptic_curves/ell_point.py` | `src/schemes/elliptic_curves/ell_point.ts` |

### Re-exports

Re-export commonly used functions at the package level (like SageMath does):

```typescript
// src/index.ts
export { gcd, lcm, factor, is_prime } from './arith/misc.js';
export { Integer, ZZ } from './rings/integer_ring.js';
export { GF } from './rings/finite_rings/finite_field_constructor.js';
```

This allows:
```typescript
import { gcd, ZZ, GF } from 'sagemath-ts';
```

### File Organization

Each module file should contain:

1. Module docstring with `@module` tag
2. Imports (external, then internal)
3. Type definitions
4. Class definitions
5. Exported functions
6. Internal/private helpers (prefixed with `_`)

---

## Ring Coercion

### Ring `__call__` Pattern

SageMath uses `Ring.__call__()` for element creation and coercion. TypeScript classes
cannot be made callable while still carrying their prototype methods, so we expose the
method under its Python name instead of making the ring object itself callable:

```typescript
// Create ring instance
const F = GF(7n);          // Finite field of order 7

// Ring call creates/coerces elements
const a = F.__call__(3n);  // Create element 3 in F
const b = F.__call__(10n); // Coerces 10 to 3 (mod 7)

F(3n)                      // TypeError: F is not a function
```

Every ring follows this: `ZZ.__call__`, `Zmod(n).__call__`, `GF(p).__call__`,
`MatrixSpace(...).__call__`. `Mod(value, modulus)` is a shorthand free function for
building a single `IntegerMod` without naming the ring.

### Automatic Coercion

Ring *elements* are more permissive than the `IntegerLike` free functions: their
arithmetic methods and their ring's `__call__` accept plain `number` as well as `bigint`,
because the value has already been reduced into a bounded ring where the 2^53 precision
cliff cannot be reached silently.

```typescript
const F = GF(7n);
const a = F.__call__(3n);

// These all work:
a.add(F.__call__(4n))   // Element + Element
a.add(4)                // Element + number (coerced)
a.add(4n)               // Element + bigint (coerced)
```

Contrast with the arbitrary-precision free functions (`gcd`, `factor`, `power_mod`, …),
which take `IntegerLike = bigint | Integer` and reject `number` outright.
@see Deviation: no-number-coercion

### Coercion Hierarchy

```
IntegerLike  -->  Ring Element
    |
    v
 toBigInt()  -->  Ring.element()
```

---

## Dependency Architecture

### Core Principle

**When SageMath delegates to an external library, we delegate to our port of that library.**

### Library Mapping

| SageMath uses | sagemath-ts uses |
|---------------|------------------|
| PARI/GP (via cypari2) | `parigp-ts` |
| FLINT | `flint-ts` |
| NTL | `ntl-ts` |
| GSL real integer powers, log and exp | `gsl-ts` |
| M4RI binary matrix multiplication | `m4ri-ts` |
| GMP | Native `bigint` |

### How to Identify Dependencies

When implementing a SageMath function, check the source for external calls:

```python
# Example: SageMath's elliptic curve cardinality
# In sage/schemes/elliptic_curves/ell_finite_field.py:

def cardinality(self, ...):
    return self.__pari__().ellcard()  # <-- Calls PARI/GP!
```

Look for:
- `__pari__()`, `pari(...)`, `cypari2` --> Use `parigp-ts`
- `flint_...`, `fmpz_...` --> Use `flint-ts`
- `ntl_...`, `ZZ_p`, `GF2X` --> Use `ntl-ts`
- `gsl_pow_int`, `gsl_sf_log`, `gsl_sf_exp` --> Use `gsl-ts`

### Correct Implementation

```typescript
// WRONG: Reimplementing the algorithm ourselves
function cardinality(): bigint {
  return this.cardinalityBSGS();  // Our own implementation
}

// CORRECT: Matching SageMath's architecture
function cardinality(): bigint {
  // SageMath calls self.__pari__().ellcard()
  return this.toPari().ellcard();
}
```

### Why This Matters

1. **Behavioral equivalence** - Same algorithms produce same results
2. **Performance parity** - PARI/GP algorithms are heavily optimized
3. **Bug compatibility** - Even quirks are preserved
4. **Maintainability** - Updates to dependency packages benefit all callers

### Implementation Order

If a dependency package doesn't exist yet:

1. Implement the dependency package FIRST (`parigp-ts`, `flint-ts`, etc.)
2. THEN implement the sagemath-ts module that calls it

---

## TypeScript Limitations

This section documents fundamental TypeScript/JavaScript limitations that prevent us from achieving full parity with SageMath's syntax and behavior.

### No Operator Overloading

TypeScript doesn't support operator overloading like Python. In SageMath, operators like `+`, `-`, `*`, `/`, and `^` work seamlessly with custom types and automatic coercion:

```python
# SageMath - implicit coercion works
F = GF(7)
a = F(3)
b = a + 5  # 5 is automatically coerced to F(5)

R.<x> = ZZ[]
f = x^2 + 3  # 3 is automatically coerced
```

In our TypeScript port, we must use explicit method calls:

```typescript
// TypeScript - must use explicit methods, and rings are called via __call__
const F = GF(7n);
const a = F.__call__(3n);
const b = a.add(F.__call__(5n));  // Must explicitly wrap 5
// OR rely on the element-level coercion:
const b2 = a.add(5);              // add() accepts bigint and number

const R = PolynomialRing(ZZ, 'x');
const x = R.gen();
const f = x.pow(2n).add(R.__call__(3n));  // Must explicitly construct
```

**Implication**: Users must explicitly construct ring/field elements or use methods that accept raw values. We cannot achieve the same syntactic convenience as SageMath.

### Polynomial Coefficient Coercion

SageMath automatically coerces mixed coefficient types when constructing polynomials:

```python
# SageMath
R.<x> = QQ[]
f = x^2 + 1/2  # Integer 2 coerced to Rational, then to QQ element
```

In TypeScript, we require explicit construction through the polynomial ring:

```typescript
// TypeScript
const R = PolynomialRing(QQ, 'x');
const x = R.gen();
const f = x.pow(2n).add(R.__call__(new Rational(1n, 2n)));  // Explicit construction needed
```

`QQ.__call__(x)` coerces a single value; `new Rational(num, den)` builds a fraction
directly. There is no callable `QQ(1n, 2n)`.

### Matrix Entry Coercion

Same limitation as polynomials - matrix entries must be explicitly constructed or passed through matrix constructors that handle coercion internally:

```python
# SageMath
M = Matrix(GF(7), [[1, 2], [3, 4]])  # Integers auto-coerced
```

```typescript
// TypeScript
const F = GF(7n);
const M = matrix(F, [[F.__call__(1n), F.__call__(2n)], [F.__call__(3n), F.__call__(4n)]]);
// OR let the constructor coerce — `matrix` routes through MatrixSpace.__call__:
const M2 = matrix(F, [[1n, 2n], [3n, 4n]]);  // equivalent
```

Note the lowercase `matrix(ring, entries)` factory (alias of `MatrixFromEntries`);
`Matrix` itself is a class and needs `new Matrix(ring, nrows, ncols, entries)`. The
generic `Matrix<R>` requires entries that are ring-element *objects*, which is why it
works over `GF(p)` but not over `ZZ` — `ZZ`'s elements are bare `bigint`s with no `.mul`.
Integer linear algebra goes through `IntegerMatrix` instead.

### None vs undefined

Python uses `None` for missing/null values, TypeScript uses `undefined`. Some SageMath functions distinguish between "parameter not provided" and "explicitly passed None". We generally treat both as `undefined`.

```python
# SageMath
def foo(x=None):
    if x is None:  # Could mean "not provided" or "explicitly None"
        ...
```

```typescript
// TypeScript
function foo(x?: SomeType): void {
  if (x === undefined) {  // "not provided"
    ...
  }
}
```

When the distinction matters, document it in the function's JSDoc comments.

### Tuple vs Array Returns

SageMath returns tuples (immutable, sometimes named). We return arrays:

```python
# SageMath
g, s, t = xgcd(6, 4)  # Returns named tuple, can also access .gcd, .s, .t
result = xgcd(6, 4)
print(result.gcd)     # Named access available
```

```typescript
// TypeScript
const [g, s, t] = xgcd(6n, 4n);  // Returns [bigint, bigint, bigint] array
```

Named tuple access is not available in our implementation. When named access is important for API usability, we return objects instead:

```typescript
// Alternative: return object when names matter
interface XgcdResult {
  gcd: bigint;
  s: bigint;
  t: bigint;
}
```

### JavaScript Number Precision

JavaScript `number` has only 53 bits of integer precision (max safe integer: 9007199254740991). Python `int` is arbitrary precision.

`IntegerLike` is `bigint | Integer`. `toBigInt()` rejects every JavaScript number,
including safe integers, to avoid silently accepting precision already lost by the caller.

```typescript
gcd(12n, 8n)    // 4n
gcd(12, 8)      // TypeError: JavaScript numbers are not accepted ...
```

The prime iterator's optional stop also accepts the numeric `Infinity` sentinel.
Finite stops still use `IntegerLike` coercion; this exception does not admit finite numbers.

---

## Summary

| Aspect | Decision |
|--------|----------|
| Integer parameters | Accept `IntegerLike`, coerce with `toBigInt()` |
| Integer returns | Return `bigint` (not `Integer`) |
| Function names | Match SageMath exactly |
| Keyword args | Use options objects |
| Errors | Custom classes matching Python |
| Module paths | Mirror SageMath structure |
| Ring elements | Callable ring pattern |
| External libs | Delegate to TypeScript ports |


## Finite-field scalar overloads and extended integer results

Finite-field element arithmetic accepts scalar `number`/`bigint` operands as well as compatible
field elements. Arithmetic checks parents before operating; explicit parent construction remains
a conversion. Equality returns false for incompatible parents. Concrete element overloads are
retained alongside scalar unions so the existing recursive `RingElement` constraints remain usable.

`Integer.pow` may return `Rational` for negative exponents; `__invert__` always returns `Rational`,
even for ±1. The Catalan convenience returns `Rational(-1/2)` at -1. Infinite integer results use
the existing string-sentinel convention described in `DEVIATIONS.md`; callers must handle the unions.

### Finite-field scalar and parent conversion

Prime-field constructors use the same Integer/Rational scalar conversion rules as Sage's
IntegerMod constructors. Explicit conversion of a Rational reduces its numerator and divides
by its denominator. Arithmetic does not implicitly convert Rational operands: there is no
canonical QQ-to-finite-field map. Integer wrappers and booleans follow integer arithmetic.

Extension parents are compared structurally (characteristic, degree, generator name and
modulus) because factories do not cache object identity. Polynomial inputs explicitly change
the coefficient ring before reduction. `fromInteger` represents base-p coordinates and checks
its range; ordinary `__call__` treats an integer as an element of the prime subfield.

### Finite-field factory options

The public `GF`/`FiniteField` aliases accept a name string, an options object, or a name
followed by options. `FiniteFieldOptions` exposes `name` and `modulus`; factory coefficient
lists include the leading term, as Sage's polynomial constructor does. The older `GFpn`
coefficient shorthand still omits the monic leading 1. A polynomial modulus is converted
to the target coefficient ring and normalized before its degree is checked. Prime-field
`gen()` stores the root of a supplied linear modulus. Order decomposition delegates to
PARI's `Z_isanypower`; names use Sage's `certify_names`/`normalize_names` validation order.

### Prime-field square-root return types

Prime-field `sqrt` maps Sage keyword arguments to `{extend, all}` and follows Sage's
`extend:true` default. Overloads return the original element type for explicit
`extend:false`, an array for `all:true`, and a union with the extension element otherwise.
Elliptic-curve point enumeration explicitly selects `extend:false` because it needs points
rational over the curve's original field. Both prime-element implementations share the
same root computation, while preserving their respective parent on base-field results.

### Expression parser operator dispatch

`misc/parser.ts` mirrors Sage's tokenizer and recursive-descent arithmetic grammar. Its
Parser retains the constructor ordering `(make_int, make_float, make_var, make_function,
implicit_multiplication)` and adds a required sixth `{binary, unary}` callback object to
perform arithmetic on TypeScript values. Variable dictionaries use LookupNameMaker; callable
name factories are preserved. No JavaScript eval or runtime source compilation is used.

Polynomial strings keep integer/rational intermediate values separate from polynomial values,
so an exponent such as `7+1` is evaluated over the integers before field reduction. Fraction
intermediates reuse FunctionFieldElement_rational's existing normalized numerator/denominator
representation. Only a resulting polynomial is accepted by the polynomial constructor; an
extension field then reduces that polynomial by its defining modulus. These internal objects
are not exposed as a new public fraction-field API. Syntax errors carry source/position fields.


### Modular arithmetic result parents

`IntegerMod` operations can change the parent through Sage's canonical quotient maps or
`QuotientFunctor.merge`: use the gcd of the two moduli, except a newly formed trivial
quotient is rejected. Prime-field operands promote to the prime field when its characteristic
divides the modulus. Arithmetic overloads expose the promoted prime-element return type;
integer/modular operands retain IntegerMod. Explicit constructor conversion is broader than
implicit arithmetic coercion. PARI `znorder` accepts separate bigint residue/modulus arguments
instead of a GEN object and retains the prime-power/valuation/lcm algorithm.


### Zero-modulus factory types

`Zmod()`, `Zmod(0)` and its aliases return the ZZ singleton; negative orders normalize only
after selecting the signed factory key. The factory uses WeakRef/FinalizationRegistry to
mirror Sage UniqueFactory's weak-value cache without retaining unused ring objects.
The direct `new IntegerModRing(n)` still requires a positive order. Conditional return types
preserve a finite-ring type for nonzero literals and expose IntegerRing | IntegerModRing for
dynamic orders. Consumers with a proved positive modulus can use the direct constructor.

`Mod(value, 0, parent?)` returns value itself, preserving its type and identity. For a
nonzero modulus, an explicit parent takes precedence and can be a modular ring or either
prime-field implementation. Its generic return type follows both the zero case and the
explicit parent. This corrects previously accepted calls' return types, so the synchronized
version increment is major. LWE retains its existing nonzero-modulus precondition; supporting
its ZZ case would additionally require integer-valued vector/polynomial samplers.


### Prime-field exponentiation

Both prime-element implementations reuse IntegerMod.pow, reflecting Sage's inheritance
rather than maintaining separate modular exponentiation/error paths. The intermediate uses
a positive-modulus IntegerModRing and the result is returned in the original prime parent.
Dispatch retains the original exponent type before Integer coercion: a Rational or string
that converts to a small negative integer still takes the GMP error path. Bigint/Integer and
integral JavaScript numbers can take the native path; the existing JS number mapping remains
explicit in comparative oracles. Powers use modular exponentiation followed by inversion,
with the same asymptotic complexity as Sage's native/GMP implementation.


### Canonical finite-ring arithmetic and Python sequences

The four finite-ring element classes share a bounded canonical-parent adapter in
`integer_mod.ts`, reflecting Sage's IntegerMod inheritance and coercion maps. Explicit
parent calls and implicit arithmetic remain separate: the latter accepts exact integer
scalars, compatible quotient maps, prime-subfield maps and structurally equal named
extensions. It rejects QQ, nonintegral JavaScript numbers and unrelated parents. Maps
are chosen before division, so promoted zero divisors report the extension backend's error.

Arithmetic overloads expose extension results when a compatible extension operand is
supplied. Native modular and both prime classes also expose `mul(string): string` and
`mul<T>(readonly T[]): T[]`, reproducing Python's reflected sequence multiplication via
IntegerMod.__index__. Array repetition is shallow; extension elements have no index
conversion. Comparative adapters use Sage's global `parent(result)` so valid Python
sequences are not mistaken for missing-parent errors. Signed 64-bit index overflow is
checked even for empty sequences; resource exhaustion remains a host-runtime boundary.


### Finite generator indices and integer conversion hooks

`gen(n)` is not an integer-constructor call. Prime and PARI extension parents test Python
truthiness, while quotient rings delegate to ZZ.gen's equality-to-zero test. A shared
bounded helper handles the port's scalar, array and finite-element representations; it
preserves the distinct IndexError messages and accepts arbitrary inputs for validation.

Modern prime and extension elements expose Sage's `_integer_(ZZ=None)` conversion hook.
The extension hook checks the polynomial degree before returning its constant coefficient;
it is separate from base-p `integer_representation()`. Both prime parent constructors use
ZZ coercion first. The legacy `check:false` option bypasses primality, while the inherited
generic-ring positive-order check remains in force, including Sage's order-one parent.

QQ.gen uses the same equality-to-zero rule as ZZ.gen: numeric zero, false and exact
ring elements equal to zero are accepted. It calls the port's eq method for ring objects
instead of interpreting strings, lists or null through the rational constructor.


### PARI quotient kernels as the extension arithmetic backend

Extension elements convert their polynomial coefficients to the existing parigp-ts dense
FpX arrays and reconstruct results in the same parent. FF_add/sub/neg correspond to FpX
operations; multiplication, division, inversion and power use the FpXQ quotient kernels.
This models the FpXQ branch of PARI's FF dispatcher with the port's existing representation.
Native packed polynomial dispatch remains an explicit performance boundary in DEVIATIONS.md.

FpXQ_pow and its FpXQ_powBig alias share signed inversion and the binary/sliding-window
algorithm from bb_group.c. Width selection changes at 512, 2^25, the unsigned 64-bit bound,
and PARI's larger bit-length thresholds. The three-argument FpXQ_inv entry point reuses the
base inverse of the pre-existing Hensel implementation; a fourth prime argument preserves
its q=p^e lifting API. Direct C-kernel comparisons validate both entry points.

Sage's extension exponent comparisons precede Integer coercion. We reproduce this ordering
using the port's numeric and finite-element representations, then pass an exact bigint to
PARI. In particular, zero to a negative fraction raises ZeroDivisionError before Integer
rejects the fraction; a nonzero base reaches the conversion error instead.


### Polynomial constructor coercion and coefficient parents

PolynomialRing and PolynomialRingBase expose `__call__(x?: unknown)` and `gen(n?: unknown)`.
They validate inputs at runtime: coefficient lists pass through the base constructor before
normalization, while exact parent identity and polynomial-base constant embedding are checked
first. The QQ FLINT singleton-list and GF2X recursive-list branches follow their respective
Sage constructors; nested coefficient lists are not QQ's separate TypeScript tuple API.
Generator indices use equality to zero, preceded by the cached Sage method's list hashability
check. This does not implement general Sage parent caching or categorical coercion maps.

The Hilbert-class-polynomial caller wraps public ZZ's bigint results in Integer for all three
coefficient factories (zero, one, and coercion). Polynomial algorithms require RingElement
methods; claiming that primitive-returning ZZ implements that interface left invalid parents.
Both prime-element classes and IntegerMod implement Sage's `_rational_` residue lift; extension
constructor input dispatch now includes modular and legacy prime elements with compatible
characteristic. Nonconstant extension-to-prime errors distinguish the same-characteristic
section map from the fallback QQ conversion.

### Dense FLINT GCD adapters

`flint-ts` exports `_fmpz_poly_gcd`, `_fmpz_poly_gcd_subresultant`,
`_fmpz_poly_gcd_heuristic`, `_fmpz_poly_gcd_modular` and `_nmod_poly_gcd` from its root.
These low-level ports accept dense, constant-first `readonly bigint[]` buffers and return
new normalized arrays. Lengths and output allocation are implicit; the old mutable
`fmpz_poly` / `nmod_poly` object APIs remain separate stubs. The heuristic returns `null`
for its native failure status. `_nmod_poly_gcd(a,b,p)` requires a prime modulus and returns
an unnormalized final remainder, as needed by integer modular reconstruction.

Sage's integer-polynomial GCD and the primitive GCD used by squarefree factorization
call this dependency. Native BigInt implements coefficient and packed integer arithmetic.
The remaining native dispatch differences are recorded in `DEVIATIONS.md` under
“Dense FLINT polynomial GCD kernels.”

### Monic polynomial result parents

Integer polynomial normalization can change the coefficient ring to QQ. The public `monic`
and `_monic` overload for `Polynomial<Integer & RingElement>` therefore returns
`Polynomial<Integer & RingElement> | Polynomial<Rational & RingElement>`; the generic
coefficient overload remains `Polynomial<C>`. The intersections adapt the existing
polymorphic-this coefficient constraint without claiming that arbitrary user-defined
coefficient classes support every Sage base extension. The runtime parent is authoritative.

The FLINT array API also exports `_nmod_poly_make_monic(a: readonly bigint[], n: bigint)`.
It returns a new coefficient array. Inputs must represent a nonzero polynomial with a unit
leading coefficient modulo `n >= 2`; `n` need not be prime. Sage checks its own error cases
before delegation. Native arithmetic representation boundaries match the GCD array adapters.


### Real double coefficient parent and resultant boundaries

`real_double.ts` mirrors Sage's distinct RDF parent; it does not reuse the approximate
MPFR `RealField`. `RealDoubleElement.value` is a binary64 number and `parent` is its RDF
field. Basic coefficient operations use the existing `RingElement` method convention.
Numeric conversion preserves Integer/Rational semantics before storing a double.

FLINT rational resultants accept two integer numerator arrays with positive denominators
and return a reduced numerator/denominator pair. The Sage wrapper restores a Rational.
The PARI real-resultant array adapter accepts constant-first binary64 arrays, converts to
64-bit native-style real mantissas, and converts its result back at the API boundary.
The Buchmann public adapter now delegates binary64 conversion to this kernel.


PARI's binary64 resultant inputs stay at a single 64-bit mantissa word during elimination.
Internal native addition and division modules preserve that storage precision through
cancellation; QFB and Buchmann elementary arithmetic now share these native kernels. The internal
real determinant stores rows and distinguishes exact integer padding zeros (`null`) from
real zeros with an accuracy exponent. Public coefficient arrays contain numbers only.


### Shared derivative argument protocol

The mirrored `misc/derivative.ts` module expands arrays of variables and counts and calls
an object's `_derivative` hook. JavaScript integral numbers and BigInts represent integer
counts; null represents Sage's default variable. A single nested array is already expanded
and is returned by reference. Polynomial aliases expose the same method function; the
runtime receiver determines the polynomial being differentiated.

### Polynomial fraction classes and native modular products

`PolynomialRing.fraction_field()` caches a parent and selects the mirrored generic,
univariate-field or small-prime `FpT` class. The public `FractionElement<C>` interface
covers their common arithmetic and accessor protocol; `_is_fraction_field_element` is an
internal brand used across the cyclic polynomial/fraction dependency graph. Parent
imports from polynomial modules are type-only where possible. `FpTElement` is a separate
implementation because the source exposes no `reduce()` method and uses structural
stored-pair equality. Generic elements retain polynomial references; native accessors
construct fresh polynomials from immutable coefficient arrays.

The native nmod product boundary accepts constant-first BigInt arrays and preserves
input arrays and repeated-operand alias information. FLINT's portable dispatcher selects
classical, KS, KS2 or KS4; binary string packing and BigInt multiplication replace native
limb buffers. KS4 reconstructs overlapping coefficients using the original reciprocal
recovery algorithm. Explicit kernels remain available for direct native differential tests.


### Fraction conversion sections and dependency initialization

Polynomial/fraction constructors distinguish canonical coefficient embeddings from partial
sections back into a polynomial or scalar ring. Exact generic sections with denominator one
return the stored numerator; native FpT sections normalize their input before testing
integrality/constancy and allocate a polynomial result. `has_coerce_map_from` exposes the
supported polynomial parent relation used by these branches.

The fraction parent module must not import Rational or Polynomial as runtime values:
those modules lead back to FpT, whose parent extends the fraction-field class during module
initialization. Structural predicates plus type-only imports keep each public fraction
subpath independently importable. Comparative tests import each module in a fresh process.

Sage's QQ polynomial numerator/denominator storage is represented here by individual
Rational coefficients. The fraction constructor's ZZ-target resolver clears their common
coefficient denominator using exact BigInt GCD/LCM and creates the integer coefficient
vector in the target polynomial parent. This internal adapter supplies the original
`resolve_fractions` algorithm without claiming new public polynomial numerator APIs.


### Scalar sections and category hooks

Generic fraction elements expose `_conversion`, `_integer_` and `_rational_`; native FpT
continues to use its separate partial section protocol. Scalar constructors choose sections
before fallback hooks, matching Sage's parent conversion ordering. Polynomial direct hooks
therefore need not report the same error as `QQ(poly)` or `ZZ(poly)`.

The fraction element module remains free of runtime scalar-parent imports. Its rational
hook calls the coefficient QQ parent when a canonical constant-polynomial section applies;
otherwise it invokes the polynomial rational hook. Integer conversions use native BigInt
unit inverses. The modular RingElement adapter checks `isUnit()` before `inv()` to supply
Sage's `inverse_of_unit` protocol without changing the existing direct `inv()` behavior.
The scalar comparison frame names the requested parent for Rational and Integer results,
whose TypeScript representations have no runtime `.parent` member.

### Integer polynomial power adapters (3.0.0)

`Polynomial.pow` dispatches by the same native coefficient backend as Sage's polynomial
subclasses. FLINT owns ZZ/QQ and word-modular algorithms; NTL owns GF2X, ZZ_pX and ZZ_pEX
power. PARI finite-extension coefficients still use NTL for their polynomial powers, as
Sage's `Polynomial_ZZ_pEX` class does. Generic polynomial fallback retains scalar, generator,
binary and characteristic-Frobenius branches.

The return type is a polynomial/fraction union for dynamic and negative integer exponents;
a nonnegative literal overload preserves fluent polynomial calls. Known nonnegative internal
callers explicitly assert that fact. Native dense arrays use ascending coefficients; QQ
conversion clears denominators on entry and reconstructs via the original coefficient ring's
integer constructor and division on exit. This supports minimal coefficient-ring adapters.
NTL extension arrays use an inner vector modulo a monic polynomial, with Kronecker stride
`2*extension_degree-1`. Portable product limits and oracle signal guards are recorded in
DEVIATIONS.md under “Polynomial Integer Powers and Portable Native Products”.

### Polynomial root and series adapters

Polynomial `pow` accepts an unknown exponent to preserve the source's backend-specific
coercions; the nonnegative literal overload remains available. Exact rational powers over
ZZ/QQ can return roots or fractions. `nth_root`, `_nth_root_series`, `inverse_series_trunc`,
`power_trunc`, `_power_trunc`, `_mul_trunc_` and `multiplication_trunc` mirror the source
methods. The last method resolves common parents; the internal product assumes matching
polynomial backends. Native FLINT kernels use constant-first bigint arrays and rational
numerator/common-denominator pairs. NTL's array extension inverse owns its Newton and
half-GCD algorithms, without making the older stateful ZZ_pX class appear implemented.
Finite coefficient roots use an internal positive-single-root bridge in the mirrored
`finite_rings/element_base.ts`; this does not advertise unsupported public finite-root options.

Truncated multiplication resolves a common canonical parent before converting precision.
Nonintegral JavaScript numbers map to RDF for this coercion, while integer-like primitives
retain the existing ZZ adapter. Polynomial variable towers can promote their aligned suffix
coefficient ring without introducing unrelated variables. `_mul_trunc_` explicitly includes
`null` in its TypeScript operand type to represent the compared Cython `None` argument.
Ordered generic products are shared by truncated products and inexact full multiplication;
real products retain the original term and aliased-square paths to preserve coefficient bits.

### Full polynomial multiplication dispatch

Full products use the same coefficient-backend classification as powers and series:
integer/rational/word modular products call FLINT array adapters, binary products call
GF2X, large-modulus products call ZZ_pX, and extension products call ZZ_pEX. Adapters
preserve input aliasing when selecting native square behavior. Rational reconstruction
uses the supplied coefficient ring's integer constructor and division, preserving strict
ring adapters. Generic exact products retain SageMath's ordered Karatsuba recursion and
ring thresholds (nested polynomial: zero; fraction field: effectively schoolbook; other:
eight). Term actions preserve left/right coefficient order even over noncommutative rings.

### Modular polynomial power adapters

`Polynomial.pow(n, modulus?)` retains its nonnegative-literal polynomial overload and
its dynamic polynomial/fraction union. The modulus is optional and can trigger canonical
parent promotion; runtime validation follows each Sage backend. FLINT modular-power
adapters take constant-first bigint arrays, while NTL extension adapters take arrays of
constant-first coefficient arrays plus the coefficient-field modulus and characteristic.
Reciprocal-series arrays are explicit native precomputation inputs, not polynomial objects.
The small GSL dependency accepts JavaScript numbers because these APIs operate on binary64
values. Public RDF integer exponents retain IntegerLike and normalize with toBigInt.

### Polynomial evaluation and composition adapters (3.4.0)

`evaluate` is the TypeScript spelling of Sage polynomial `__call__`; `compose` invokes
that same dispatch for a polynomial argument. Native QQ composition and scalar evaluation
use FLINT, word-modular evaluation/composition use FLINT, and extension scalar evaluation
uses NTL. Large-modulus scalar evaluation uses the NTL array adapter. The bundled integer
wrapper tests the argument tuple instead of its first item in its composition branch;
we preserve its reachable generic composition path. Integer scalar evaluation still uses
FLINT. Native rational results are reconstructed through the supplied coefficient parent,
so custom coefficient adapters remain usable.

Generic evaluation unwraps one argument-list level, evaluates nested coefficients when
extra arguments are supplied, and resolves a common canonical parent before constant/zero
shortcuts. JavaScript numbers represent Python floats in this API; use bigint or Integer
for integral inputs. Evaluation at a number can return a number or an RDF element, depending
on the coefficient parent. The same-coefficient overload retains its coefficient result;
typed unknown inputs expose the broader result union. Minimal coefficient adapters without
runtime parents are recognized only when they share the coefficient constructor.

The mirrored `polynomial_compiled.ts` exports `CompiledPolynomialFunction` through
`sagemath-ts/rings/polynomial`. It constructs the original sparse Horner instruction graph
and fills gaps with the binary algorithm. A max heap replaces BinaryTree's ordered gap set;
an iterative evaluation order avoids JavaScript call-stack limits. Reference counts release
intermediate values after their last use. A WeakMap caches plans by polynomial. The optional
arithmetic argument to `eval` is an internal bridge for Sage's dynamic scalar coercions;
it does not convert coefficients before deciding their original zero/nonzero graph structure.

`PolynomialRing.gen()` caches its generator. The constructor's optional third `is_gen`
boolean records the original distinguished-generator flag in a WeakSet; `is_gen()` follows
native equality-based overrides for FLINT/NTL template backends and the flag elsewhere.
Generic representations honor that flag even over the zero ring. Polynomial ring objects
are still not globally interned; evaluation recognizes equivalent same-variable coefficient
parents for the original identity shortcut.

FLINT evaluation returns unreduced numerator/denominator pairs at its internal array
boundary. Composition returns normalized arrays (and a positive common denominator for QQ).
NTL exports the reserved native name `eval` through an export alias in `ZZ_pEX.ts`, with
`ZZ_pEX_eval` at the package root; the prime array adapter is `ZZ_pX_evaluate`. These APIs
leave the older stateful NTL polynomial classes' unrelated stubs explicit.

### Polynomial matrix actions (3.5.0)

Matrix evaluation handles the existing generic, integer and modular matrix classes through
explicit scalar/matrix actions. Coercion resolves the result coefficient ring while sparse
power nodes keep the original matrix ring until a coefficient action promotes them. A WeakMap
tracks original default versus explicitly generic matrix parents across intermediate values;
this preserves conversion and error-parent decisions without changing public matrix storage.
The return overload exposes the union of supported matrix classes because promotion can
change the concrete representation. Specialized binary matrices and additional polynomial coefficient domains remain open.

The mirrored `matrix/change_ring.ts` exports `integer_to_real_double_dense` through the matrix
subpath and delegates entries to the new FLINT `fmpz/get.ts` array-independent bigint adapter.
FLINT truncates integer doubles; scalar RDF construction retains its distinct rounding rule.
The result currently uses `Matrix<RealDoubleElement>` storage for Sage's native real matrix.
Existing matrix arithmetic backends and their documented limits are inherited by evaluation.

### Packed binary matrix multiplication and evaluation (3.6.0)

`Matrix_mod2_dense.mul` now follows Sage's M4RI delegation. `_multiply_classical`,
`_multiply_m4rm` and `_multiply_strassen` share a private storage adapter. Matrix entries
remain the existing bit arrays at the public boundary; 64-bit chunks pack/unpack rows in
linear time. The dependency `mzd_t` represents each packed row by a bigint, with bit j at
column j. Owned immutable matrices replace native output pointers and writable windows.

`mzd_mul` follows the native seven-product Bodrato schedule, including a separate square
schedule and word-aligned rectangular tails. M4RM uses eight Gray-code tables and blocked
row traversal; its classical fallback uses transposed packed parities or selected-row XOR.
A fixed 4 MiB cache model supplies portable tuning defaults (4096 Strassen, 2048 row block).
These values match the oracle's cache model but do not claim hardware-adaptive tuning.

Evaluation recognizes the specialized binary class, retains it when its coefficient parent
is GF(2), and exposes it in the matrix return union. An argument-list overload and a generic
conditional for mixed matrix/scalar arguments preserve the possible matrix result in types.
A compiler test verifies that those results cannot silently narrow to coefficient scalars.

### Binary matrix indexing adapters

`Matrix_mod2_dense.get/set` represent Python's two-index subscription operations. Their
indices, and `row/submatrix` indices, accept IntegerLike, Rational and boolean while retaining
integer-valued number arguments used by the existing array API. Conversion preserves the
original signed C-index bounds and rational truncation/exactness distinction. Dense storage
is still a bit array. Native submatrix delegation packs only the source words intersecting
the requested slice and returns an owned bit array, with work proportional to the output
plus at most two boundary words per row.

### Native real literal conversion and binary matrix density

`RealLiteral` preserves the original decimal text and delegates conversion to the new
`@sagemath-ts/mpfr-ts` dependency. A private weak map attaches native sign, kind, precision,
mantissa and exponent to literals and their precision-converted `RealNumber` results.
Formatting, exact rational extraction, sign/mantissa/exponent, absolute value and negation
observe that state. Ordinary real arithmetic retains its documented binary64 limitations.
The bundled `RealLiteral.__float__` reparses the literal at 53 bits; `toNumber()` follows
that override instead of rounding its initially stored higher-precision value twice.

The dependency mirrors MPFR's init2, set_str/strtofr, get_d and get_str boundaries. Its bounded
decimal conversion uses exact BigInt base powers and quotient/remainder rounding instead
of MPFR's adaptive interval/Ziv loop. Native allocated output buffers become a digits string
and exponent pair. Unsupported bases, modes and working sizes fail explicitly; see
`DEVIATIONS.md` for the supported domain.

Binary matrix `density()` returns `Rational`, or `0n` for an empty matrix. `density(true)`
delegates sampling to M4RI's `mzd_density` and wraps the returned double with the original
`create_RealNumber` precision rule. This restores the distinct exact and approximate return
types and the original native word-sampling behavior.

### Binary elimination delegation and caches

Binary matrix rank and elimination now delegate to the M4RI dependency's Russian PLE,
recursive PLE/PLUQ, six-table M4RI elimination, density-based crossover and unit triangular
solvers. Packed immutable dependency adapters return owned matrices with ranks and, for
factorizations, native transposition lists. BigInt interval operations replace native word
copies; panel/table and recursive algorithm structure is retained with a portable 4 MiB
cache tuning model. The explicitly requested classical Sage algorithm remains at the
Sage layer and always produces reduced form.

The matrix class caches rank, pivots, echelon status and echelon forms, following original
validation order. Computed forms are immutable; copies are mutable. Native unsafe lexical
ordering swaps preserve the original cache-retention behavior. `_clear_cache` is public,
while the mutability check remains private. Empty `echelonize` returns the matrix itself;
its exported return type therefore includes `this`.

Standalone binary `ple` and `pluq` now delegate every algorithm to M4RI, including the
explicitly requested naive algorithms. The shared parameter conversion mirrors Cython's
numeric `__int__` protocol: bigint/Integer, legacy numbers, Rational and boolean inputs;
conversion and signed C-int bounds precede algorithm validation. The low-level adapter
returns factors and transposition lists, and the Sage layer assembles a fresh mutable matrix.
The pinned PLE kernel uses seven tables; elimination uses six and triangular solving eight.


Binary row/column swap parameters use the same `BinaryMatrixIndex` adapter as row
extraction: convert both to signed `Py_ssize_t` before mutability and cache invalidation.
Permutation number arrays encode the images of `SymmetricGroup(length)` elements;
normalization precedes the matrix call, and original cycle order preserves partial failures.
Binary inversion delegates to the M4RI port's aligned augmented-window implementation;
the high-level rank check stays in the Sage layer.


M4RI solve adapters retain both native mutated arguments in `{matrix, rhs, status}`;
kernel computation similarly returns `{matrix, kernel}` with the original null-kernel
convention. The Sage layer owns shape validation, square-system dispatch, empty shortcuts,
mutable result construction and requested kernel-basis transformations. Solver padding
copies extend to `max(nrows, ncols)` exactly as in the bundled Sage implementation.


Binary column caches mirror Sage's two levels of ownership: a shared mutable outer list
and immutable column vectors. Frozen number-array proxies preserve entry-assignment errors
without introducing a separate vector class. Default `columns()` copies only the outer
list. M4RI transpose uses recursive cache-sized regions and the original six word masks;
packed BigInt rows replace native pointer strides and small unrolled/SIMD variants.


`matrix/matrix0.ts` provides the original ASCII string-rendering operation over a matrix
snapshot interface. The binary subclass delegates mappings and layout to it; lazy entry
access preserves the original validation/mapping/snapshot order. Public subdivision lists
use bigints; private storage keeps the visible lines without native sentinel endpoints.
Copy/transpose/augmentation propagate that metadata using the original method rules.
JavaScript records use numeric property keys for the two field elements, and mapping
callbacks receive number bits. Separate subdivision arguments avoid tuple/list ambiguity.


`matrix/args.ts` adapts the original MatrixArgs list/scalar iteration for the binary
constructor. It checks excess input before yielding the final item, preserving conversion
error precedence. Field construction is lazy to avoid initializing the coefficient/matrix
import cycle. Dimensions map through positional MatrixSpace integer conversion and native
signed dimension limits. Binary cached rows use immutable number-array proxies; default
row extraction returns an independent mutable array.


Binary randomization reads `current_randstate()` through its existing GMP-compatible
`random_bits`, `c_random` and `c_rand_double` operations. A full word always consumes low
then high 32-bit draws; visible columns read the raw combined word from most significant
to least significant bit, matching native `m4ri_swap_bits`. Factory functions retain their
GF(2)-specific names while following Sage's omitted/explicit density dispatch.


`types/python_float.ts` supplies the Python `float` boundary for binary density arguments.
Sage numeric wrappers retain their own conversion behavior; native bigint models Python
int at this boundary. Numeric strings use JavaScript's binary64 parser after CPython grammar
validation and pinned Unicode 15 digit/space normalization. A pinned printable-character
index supplies error repr independently of the JavaScript engine's Unicode version.
Uint8Array maps to Python bytes. The native test helper uses CPython 3.12 for this boundary.


Gaussian lattice construction uses omitted options/sigma for the original default of 1.
An omitted center maps to `0`; explicit `null` maps to Sage's deferred `None` center.
Center accessors therefore return nullable arrays. The low-level GPV `_call()` preserves
Sage's literal-zero empty-loop result, while `sampleExact()` models the public wrapper.
Exception comparisons require both class and message equality; narrowly documented
version adapters belong on the original-runtime side of an area dispatcher.


Gaussian parameter strings use `types/real_format.ts`: binary64 input is converted through
`mpfr_set_d`, and `mpfr_get_str` supplies Sage's truncated decimal representation. The
integer sampler's repr then performs Decimal-style quantization of that representation;
its `%f` validation errors independently format the original binary float.


### Gaussian numeric options

`DiscreteGaussianOptions.c` and `.tau` accept `IntegerLike | number`, matching the
constructor's existing real-center and numeric-tail paths. The `withOptions` convenience
passes overrides to that constructor before conversion, preserving validation order.
Numeric tail cutoffs still must be integral at runtime. The exported internal options
interface remains a narrower numeric specialization for GPV code.


### Prime ordinal delegation

`nth_prime` accepts `IntegerLike` and delegates to `parigp-ts/prime`, matching Sage's
PARI call. `primes_first_n` also accepts `IntegerLike` while retaining its earlier numeric
count signature. This numeric-count exception is deliberate: Sage's bound checks run
before PARI truncates a positive fractional count. The implementation then uses
`prime_range(nth_prime(n)+1)`, as the original does.


### Dedekind dependency results

`arith.dedekind_sum` selects the original FLINT or PARI backend before coercing IntegerLike
arguments. `flint-ts/fmpq_dedekind_sum` and `parigp-ts/sumdedekind` return reduced BigInt
numerator/denominator tuples; the Sage wrapper preserves its existing public object shape.
The two native algorithms retain independent recurrences and sign conventions.


### Integer rounding inputs

The three integer-rounding helpers accept `FloatInput`, the existing numeric/string/byte
input union used by Python float conversion. Floor/ceiling dispatch to an object's native
method before attempting float conversion, preserving Integer/Rational/RealNumber precision.
`bigint` represents a Sage Integer and returns unchanged. Booleans use Python numeric values.
Truncation performs the original sign test then selects floor or ceiling; strings, bytes and
null fail that comparison before any float conversion. This wider conversion interface is
specific to the rounding helpers and does not change ordinary IntegerLike coercion.


### Observable Python sorting behavior

`types/python_sort.ts` provides the CPython 3.12.5 stable powersort used by complex
display ordering. A JavaScript stable sort is insufficient when keys contain unordered
NaNs: comparison order itself changes the output. The adapter copies its input and
accepts a less-than callback; key construction and Python tuple comparison remain at the
Sage call site. Its natural-run/minrun, balanced merge schedule and galloping retain the
original asymptotic behavior. The adjacent CPYTHON-LICENSE.txt preserves attribution.
The supported binary64 complex representation and runtime-version boundary are documented
in DEVIATIONS.md; native comparisons exercise both output and comparison order.


### Integer magnitude metadata and factor extraction

Integer captures its magnitude bit length once at construction, storing it in a private
field beside the immutable BigInt value. Copy construction from Integer reuses this field.
This represents GMP's stored limb-size metadata and keeps bit_length/nbits queries constant
time; native BigInt otherwise provides only a linear binary-string conversion for that
metadata. The positive-factor mpz_remove adapter in types/gmp.ts preserves GMP's low-bit
and power-doubling algorithms using BigInt arithmetic. Scalar methods perform Sage's
validation and infinity handling before entering this internal finite-factor helper.

The free valuation overload for method-bearing objects precedes the Integer-coercible
overload, so objects exposing both protocols retain their method result type. Missing
JavaScript properties take the original AttributeError fallback explicitly. Integer
coercion runs outside the method lookup/call exception handler, so failed conversion
hooks are invoked once and their exceptions propagate unchanged.


### PARI factorial and wide real exponents

The PARI factorial selector returns its native real representation, not an integer.
`MpReal<E extends number | bigint = number>` stores sign `s`, exponent `e`, mantissa
`m` and bit precision `p`; the represented value is `s*m*2^(e+1-p)`. Existing qfb
kernels keep the default number exponent. Factorial and exponential entry points
return `MpReal<bigint>` so PARI's signed exponent range is preserved without Number
rounding or allocating the represented integer. The default/GMP free factorial
returns bigint; overloads make the backend-dependent return type explicit.

The positive gamma branch performs bounded-precision arithmetic on normalized
mantissas and carries the final binary scale separately as bigint. New real kernels
mirror `trans1.c`, `trans2.c`, `bern.c` and `kernel/none`; balanced products and native
binary/sliding-window schedules live in `bb_group.ts` and `ZV.ts`. The Bernoulli
cache uses the original zeta-based algorithm rather than a quadratic recurrence.


### Native PARI square-root component types

The dependency's `sqrtr` returns `MpReal | MpComplex`, where `MpComplex` has `re` and
`im` components of type `bigint | MpReal`. In particular, a negative real root retains
PARI's exact integer zero in the real component and an inexact real imaginary component.
Literal nonnegative/negative sign overloads narrow the return type; an arbitrary
`MpReal` requires a component check. Existing positive-only native callers use
`sqrtr_abs`. This representation does not imply a general cypari2 Gen API.

The GMP support under `parigp-ts/src/kernel/gmp/sqrtrem.ts` packs 64-bit limbs into
BigInts. It retains the normalized reciprocal-table base cases and Karatsuba
root/remainder split from GMP 6.3.0. BigInt arithmetic absorbs machine carries;
`sqrtremi` maps its optional output pointer to a `[root, remainder]` tuple. Integer
root-only wrappers reuse the remainder-producing kernel, retaining divide-and-conquer
complexity while computing a remainder they discard. Native buffer aliasing and
machine-specific assembly tuning are outside the immutable BigInt representation.


### Buchmann binary64 and comparison delegation

`buch.ts` delegates `dbltor`/`rtodbl` to `kernel/none/mp_indep.ts`, and `cmprr`
to `kernel/none/cmp.ts`. Its optional `dbltor` working precision composes native
`rtor` with `nbits2prec`; omission retains the raw native allocation. Real comparison
orders signs/exponents before aligning packed mantissas by precision, so its storage
cost depends on precision rather than exponent distance. Native binary64 errors use
the exported `PariError`; Sage resultant fallback catches that dependency exception.


### Shared Buchmann elementary real kernel

Buchmann constructors, arithmetic and integer conversion delegate to qfb's native
real entry points. The requested-bit conveniences (`itor`, `real_0`, `real_1`,
`setprec`) round up with `nbits2prec`; stored nonzero real operands use complete words.
Canonical zero construction uses `real_0_bit(e)` with no allocation argument.
The immutable `setprec` wrapper is explicitly an `rtor` conversion, not C's mutable
precision-header macro. `mulur`/`divru` accept BigInt for all unsigned-word operands.

`kernel/gmp/mp.ts:mantissa2nr` converts the signed packed mantissa using shifts.
`truncr` validates that all integer bits are represented before calling it;
`gcvtoi` permits missing bits and reports their exponent. `kernel/none/mp_indep.ts`
implements `trunc2nr` for the scaled, precision-loss-permitting `gtrunc2n` route used
by algebraic dependencies. These entry points intentionally have different guards.


### Shared logarithms and Buchmann transcendentals

Buchmann `sqrtr`, `mplog2`, `logr_abs` and `expr` delegate to the native dependency
kernels. `expr`, exported as public `mpexp`, retains `MpReal<bigint>`; binary64
conversion accepts that exponent representation without passing it through a
precision-losing arithmetic operation. Complex square-root results use `MpComplex`.

`trans1.ts` owns logarithm dispatch, its progressive-precision series recurrence,
log(2) cache and the generic `abpq_init`/`abpq_sum` binary splitter. `trans2.ts:atanhuu`
builds its rational series coefficients and delegates its final conversion to `rdivii`.
Qfb logarithm wrappers delegate to `trans1.ts`. The splitter preserves the native
one/two/three-term bases and balanced merge; output structs become returned records.
PARI's logarithm changes logical precision without rounding stored mantissas, so those
steps use truncated views; assignments and final conversions use native `rtor`.


`trans1.exp1r_abs` uses the same `MpReal<bigint>` result representation as `mpexp`.
During repeated doubling, once the exponent is large enough that adding two cannot
change any retained mantissa bit, it keeps a separate BigInt scale. The mantissa
still follows native `mulrr` with distinct operands, preserving truncated-product
rounding rather than dispatching to the identity-sensitive square kernel. The
range-reduced `mpexp` caller converts the bounded residual result back to a number
exponent before entering the shared elementary arithmetic kernel.


Rational function-field factory keys use a WeakMap of constant-field parent objects
and a map of variable names, mirroring Sage's UniqueFactory within the port's parent
identity model. Variable renaming returns the native triple of field and two maps;
the maps are plain functions, consistent with existing residue-field map adapters.
Their polynomial coefficients are copied into the destination variable's ring.

Function-field polynomial valuation arguments now pass through the existing fraction
parent and polynomial conversion machinery, retaining backend-specific integrality
checks. Numerator and denominator valuations are both evaluated before combining
finite values or the infinity sentinel; zero numerators do not bypass denominator
validation. Polynomial zero-factorization guards are likewise reached through the
existing polynomial implementation.


### Rational function fields and their fraction representation

`RationalFunctionField<C>` caches its polynomial ring's fraction parent in `_field`.
`field()` exposes that parent, and `FunctionFieldElement_rational<C>` stores the actual
`FractionElement<C>` as `_x`. Numerator and denominator accessors delegate each time,
so lazy fraction normalization remains visible. Arithmetic uses the fraction element's
operations; negation and powers retain Sage's generic element algorithms. The existing
polynomial-pair constructor is retained as a convenience, constructing the selected
fraction class before wrapping it. Fraction string fallback shares the polynomial
arithmetic parser through an internal evaluation entry point.


`IdealMonoid` implements Sage's `UniqueRepresentation` with a weak cache keyed by the
order object, including direct constructor calls. Its call operation retains an ideal
already in that monoid, matching generic Sage parent construction. The rational
infinite-order basis uses a `[bigint]` tuple because its bundled source returns a literal
Python integer; the finite-order basis remains an array of function-field elements.


FLINT `ulong_extras` kernels expose word values as BigInt with explicit 64-bit input
bounds. Normalized preinversion uses exact two-limb division; modular products use
BigInt remainders, retaining the public reciprocal argument for the mirrored native
signature. Word square checking uses the native residue filters plus an integer Newton
candidate. Jacobi symbols remain machine-sized numbers, while square predicates are
Boolean, consistent with the existing TypeScript numerical conventions.


Fraction elements expose `is_square()` and `sqrt(extend, all)` so rational function
fields can delegate to their actual representation. The generic fraction class also
supports `is_square(true)` and the inherited optional root name; FpT retains its own
narrower native signature. Polynomial optional roots use `[boolean, Polynomial | null]`.
The integer squarefree path passes primitive BigInt arrays to NTL and reconstructs
factors in the original coefficient parent, with a constant pair carrying the unit.

FLINT polynomial roots and series use dense BigInt arrays and explicit prime moduli.
GR helper signatures are modular instantiations; truncation lengths/cutoffs are safe
integer numbers. NTL integer GCD and squarefree routines similarly return dense arrays;
its existing modular arithmetic helper is shared between half-GCD and integer CRT
reconstruction. Internal helper objects remain explicitly marked `@internal`.


Squarefree decomposition dispatch now follows the polynomial backend. The characteristic-
zero category algorithm and finite-field Frobenius algorithm live in mirrored modules;
word prime coefficients call the FLINT dependency. Their outputs use the existing pair
array adapter, with a constant pair representing nontrivial units (including ZZ zero).
A generic IntegerModRing's public squarefree-method availability is separate from its
ordinary factorization capability. Dynamic category mutations remain outside this model.

The dedicated GF2 parent owns a pair of canonical coefficient objects, mirroring the
native IntegerMod_int table for Sage-integer inputs. Direct element constructors allocate
fresh objects; parent conversion preserves existing elements. Coefficient roots use the
same `{extend, all}` options convention as the other prime-field implementations.


### PARI squarefree component vectors

`parigp-ts/FpX_factor.ts` mirrors the native word-prime Shoup squarefree decomposition
and large-prime Yun algorithm. Dense little-endian coefficient arrays replace PARI
polynomials, and vector index `i` represents multiplicity `i+1`. Polynomial-one entries
preserve missing multiplicities; this API does not split components into irreducibles.
The portable native-word boundary is 64 bits. Exact BigInt modular arithmetic replaces
machine-word coefficient reciprocals; polynomial Barrett inverses remain cached.
Word squarefree arithmetic delegates to the Flx kernels. Full `FpX_factor` and
`Flx_factor` return irreducible `[polynomial, multiplicity]` pairs, using the
native small-degree cases and Shoup EDF/EDF_simple strategy threshold. DDF,
trace, minimal-polynomial and recursive pullback operations retain their quotient
contexts. Galois support callers delegate to this pipeline after their existing
squarefree guard. Sage large-prime factorization delegates to PARI and restores
the leading unit in its own pair-array representation.


### Place identities and polynomial coefficient denominators

PlaceSet uses a weak cache by function-field parent. Valuation rings use a weak parent
cache and the existing place key, preserving equal-place reuse without scanning all
places. Residue triples are cached by name; callable-map adapters perform Sage's
domain conversion before invoking the direct map. The map object/introspection
deviation is unchanged.

QQ polynomial numerator/denominator access uses canonical FLINT numerator storage.
Numerators return Integer object coefficients through a cached ZZ polynomial parent;
the integer scalar API continues to return bigint. Generic coefficient denominators
use native LCM, with QQ polynomial LCM routed through FLINT's primitive integer
polynomial LCM. Sage's integer polynomial LCM deliberately retains the sign of its
product, whereas the FLINT dense LCM kernel normalizes the leading sign.


### Divisor function-space adapters

The existing public Riemann-Roch adapter retains a dimension plus two callable maps.
It now caches this triple as Sage does and applies native coordinate-array length
rules before lifting. Native vector backend selection matters at prime modulus
2147483647: generic dense vectors expand an empty list to zero, while QQ and smaller
prime fields use specialized dense constructors. The empty-sum result at dimension
zero preserves SetMorphism's native TypeError. The direct echelon helper retains
None pivots for zero/dependent rows and exact list-index coordinate errors.

### Lazy polynomial enumeration and projective constructor validation

Polynomial iterators use the mirrored `misc/mrange.ts` Cartesian iterator with
restartable `[Symbol.iterator]` inputs, preserving native pull order and row-copy
semantics without materializing coefficient rings. Native list lengths map to JS
`length`; ring cardinality remains a method, including QQ's inherited infinity.

The hyperelliptic model specialization validates the same homogenized equation as
`ProjectivePlaneCurve.__init__`. If its degree is below two, it constructs the
multivariate numerator after multiplication by `z^(2-d)` and reports the native
coordinate-ring conversion error. The `y^2` term proves there is no cancellation
with the denominator, so this does not need a general multivariate fraction field.
Prime-field display follows Singular's `npWrite` balanced representatives and Sage's
`can_convert_to_singular` characteristic boundary. The original coefficient printer
was read from [Singular modulop.cc](https://github.com/Singular/Singular/blob/spielwiese/libpolys/coeffs/modulop.cc)
on 2026-09-11 and verified against live Sage at both sides of the boundary.

### Rational function-field multiplication matrices

`FunctionFieldElement_rational.matrix()` retains the existing nested-array API for
its one-dimensional matrix. Both arrays are frozen; assignment proxies preserve
Sage's immutable-matrix ValueError. `M.map(row => row.slice())` is the mutable-copy
adapter. Matrix entries are formed by multiplication by the basis image (one),
trace adds that entry to the field zero, and the 1-by-1 determinant is the entry.
These are the native `element.pyx:matrix` and `matrix2.pyx:trace/determinant` steps;
returning the original function element would incorrectly retain unreduced data.
Fresh method lookups return fresh matrices, matching Cython's cached-method caller
behavior. Rational residue elements carry no parent object in this port, so the
comparative adapter identifies their unique parent QQ from the Rational class.

### Function-field predicate parent and category protocols

The bundled predicates perform native instance checks followed by parent/category
membership. The port accepts callable `parent()` protocols and its existing element
`parent` properties. BigInt/Integer and Rational scalars use their unique ZZ/QQ
parents; generic matrices use a non-function-field MatrixSpace token. Structural
ring/order parents retain the native unbound category-method errors when passed as
elements. Parent lookup occurs twice when the first result is not a FunctionField.
Category membership reads `category` and `is_subcategory` once each, including
proxy-supplied methods, and catches only AttributeError during category lookup/call.
The internal category token prints `Category of function fields`; this specialization
does not implement Sage's general category hierarchy. Internal source-module helpers
are not part of the package export map.


The rational multiplication-matrix specialization also retains the native
`free_module(map=True)` parent-cache initialization state in a WeakSet. Its cache
key uses only the map flag, so a successful matrix/trace/norm call causes later
hashable base arguments to be ignored for the same parent; failed validation does
not warm it. Matrix `cached_method` argument hashing still rejects list/dict/set
adapters before consulting that state. This does not cache matrix arrays or add
the currently absent public free-module/vector-space API.

### Variable-map domains and polynomial dictionaries

Variable-renaming maps apply the source field constructor before the native map
operation; identity maps do this as well and retain existing source elements.
Function-field conversion delegates singleton-list handling to the fraction field.
Map failures preserve the native domain-conversion diagnostic.

Polynomial coefficient dictionaries accept Maps and numeric-key Records. Maps retain
insertion order and tuple-shaped exponent keys; Records follow JavaScript's numeric
property enumeration order, which the comparative adapter reproduces. Integer
coefficients continue to use the existing Integer-object ring adapter. QQ uses
coefficient setters, word-prime/NTL templates accumulate polynomial-valued entries,
and the generic dense branch converts a coefficient list. No backend is replaced by
a different factoring or arithmetic algorithm.

### FLINT word random-state dependency

`flint_rand_t` stores the two native 64-bit words as BigInts. Initialization returns
that object; seeding mutates it with unsigned-word wrapping, and seed retrieval
returns a copied pair. The deprecated GMP random-state field is unused by these
kernels. Native clear is a no-op. The generator, high-word bounded sampling and
edge-biased test values follow the bundled FLINT source, including seed consumption.
The comparative oracle compiles source-extracted C bodies with platform type/macro
definitions, avoiding accidental dependence on the older installed FLINT generator.

### Native FLINT factor dependencies

The dense word-polynomial factor interface returns `[unit, Array<[coefficients,
exponent]>]`, preserving C insertion order. Sage's `Polynomial.factor()` converts
these factors to its coefficient ring and applies Sage's factor ordering, including
its existing constant-factor adapter. DDF uses the same tuple representation with
irreducible degree instead of multiplicity. EDF probability returns a candidate or
`null` and mutates the explicit `flint_rand_t`; recursive EDF initializes a fresh
native state for each recursive split. All coefficients use `bigint`; polynomial
lengths, exponents and matrix dimensions use safe integer `number` indices.

Brent–Kung composition shares dense modular matrix products. The matrix API accepts
an optional output column count so zero-inner-dimension products retain their shape.
The bundled C test oracle compiles a small ABI wrapper linked to the actual bundled
FLINT library. It builds FLINT in a temporary cache, or uses the build directory in
`SAGEMATH_AUDIT_FLINT_BUILD`; no library is installed and the reference tree is not
modified. Its development headers and Autotools/GMP/MPFR build prerequisites must
be available when the cache is empty.

### PARI global random state

The 64-bit XORGEN stream in `parigp-ts/random.ts` mirrors `basemath/random.c`.
Its module state models one PARI thread's global state; `getrand()`/`setrand(bigint)`
use the native 66-word positive-integer encoding. Integer word arithmetic wraps
explicitly, while sampled integers remain arbitrary-precision bigint values.
Polynomial samplers live in their mirrored `Flx.ts`/`FpX.ts` modules. Native
comparisons run an isolated executable linked to the bundled PARI library so
Sage's older in-process PARI cannot interfere. `SAGEMATH_AUDIT_PARI_BUILD` can select
a configured source build; otherwise the oracle builds a temporary cached copy.

### Binary PARI factor representation

The `F2x.ts` and `F2v.ts` dependencies use bigint coefficient bits and packed matrix
columns. Their array return values replace PARI's mutable GEN vectors. Factorization
remains in the mirrored `FpX_factor.ts` module and shares `random.ts`; the global
native random stream is compared after each factorization. Binary Sage polynomial
objects are converted at the `Polynomial.factor()` boundary and sorted using the
same Sage comparator as the word-prime branch. The separate GF2X convenience API
is not used by this route.


PARI matrix multiplication is shared through ZV.ts and FpV.ts. The established
buch.ts ZM_mul entry point delegates to ZV.ts without changing its 1-indexed
column-array API. Internal row arrays share the native seven-product Winograd
schedule, with exact BigInt arithmetic and the source's dimension/word cutoffs.
Large balanced integer products reconstruct from high-bit word primes with a
balanced CRT tree; synchronous execution replaces native worker scheduling.
F2v.ts keeps packed-column XOR multiplication for the binary FpM dispatch.


PARI modular composition lives in FpX.ts/Flx.ts, with the existing galconj.ts
entry points delegating to them. RgX.ts supplies the original optimal table-size
search. Shared polynomial helpers build the source's blocked coefficient matrices
and call FpM_mul or Flm_mul according to the native route. Automorphism powers
reuse tables; additive traces use the existing native generic powering schedule.


PARI FpX products and squares now delegate to Flx.ts or ZX.ts; the established
ffinit.ts and galconj.ts entry points retain their signatures. Packing uses linear
coefficient-string assembly and signed centered recovery. A 2^16-bit temporary
budget selects the native Karatsuba splits before large packed intermediates
become expensive or exceed the VM limit, following the portable FLINT/NTL policy. Quotient power
construction calls the separate square kernel at the source's square steps.

PARI polynomial division shares `_polynomial_division.ts` for the native coefficient
back substitution, reciprocal basecase/Newton schedule and Barrett block assembly.
FpX.ts and Flx.ts own their distinct 64-bit GMP cutoffs and conversion/early-return
branches; ffinit.ts retains its public signatures and delegates to those modules.
Newton uses the existing `quadratic_prec_mask` dependency and the native polynomial
multiplication ports. Flat ascending arrays represent polynomials; quotient/remainder
pairs replace the native output pointer. Reduction-object caching remains separate.

PARI gcd and Bézout functions retain native coefficient scale; monic normalization
belongs to callers that require it. `_polynomial_gcd.ts` shares half-GCD recursion,
seven-product 2x2 polynomial matrices and backward cofactor reconstruction while
FpX.ts/Flx.ts select their native thresholds. The small transformation matrix API
uses zero-indexed rows of polynomial arrays; it is distinct from the one-indexed
column convention of the existing integer-matrix adapters. `[M,a,b]` replaces
native optional transformed-polynomial pointers. An internal cofactor-selection
parameter preserves the original inversion call without computing an unused u.
Quotient powers share bb_group.ts's generic nonzero exponent-magnitude schedule;
the FpXQ word dispatch boundary is 2^63, independently of multiplication's 2^64.


PARI minimal-polynomial kernels share a private polynomial reducer through the
power-table and matrix-evaluation helpers. Each invocation computes the native
Barrett inverse once at the GMP tuning threshold and retains it for transposed
products and all quotient operations. Public signatures continue to accept
ascending coefficient arrays; this does not introduce a mutable global cache.

The shared PARI polynomial quotient context now also serves distinct-degree
factorization. It retains a precomputed inverse when converting an FpX reduction
object to word arithmetic, including FpXQ_pow's signed-word boundary and
FpXQ_powers' unsigned-word boundary. Native DDF factor/degree vectors use the
existing F2x-style pair representation; galconj's public FpX_ddf Map adapter is
retained for compatibility.

### PARI modular square-root kernels

`ff.ts` preserves the original `Fp_sqrt`, `Fp_sqrt_i` and `Fl_sqrt` entry points;
private `_modular_sqrt.ts` shares the Tonelli loop while retaining the different
word/large-prime generator order. Large primes use the original signed-small
shortcuts and the cost test `e*(e-1) > 20 + 8*floor(log2(p))` for Cipolla.
Cipolla delegates its fused left-to-right schedule to `bb_group.gen_pow_fold`.
Prime enumeration reuses the existing segmented sieve and next-prime search;
representation and invalid-input adapters are recorded in DEVIATIONS.md.


### Native polynomial root splitting

FpX_factor.ts owns the prime-field root/count/predicate entry points; existing
Galois exports delegate there. The shared root kernel follows PARI's deterministic
shift splitter and word-field square/nonsquare cutout. It uses the same private
nonsquare search as modular square roots. The Galois helper's old local RNG and
Cantor splitter have been removed. Word exponent-two root counts use direct
squaring/remainder before any reciprocal cache is constructed.

Dense coefficient arrays and plain root arrays replace GEN polynomial/column
storage. Native raw-degree shortcuts and signed-word sorting are retained, even
where they look unusual in TypeScript; the public API reference illustrates them.

### PARI extension-polynomial coefficient storage

`FpXX.ts`, `FlxX.ts`, `F2x.ts` and `polarit3.ts` mirror the source locations of
extension-polynomial product, square, reduction and normalization entry points.
A shared private kernel preserves native coefficient tags: `bigint` is GEN INT,
`bigint[]` is GEN POL. Word polynomials have only array coefficients; each binary
coefficient is a packed bigint. There is no implicit conversion between an integer
and a constant polynomial, because PARI selects different algorithms for them.

Generic and word products pack blocks of width `2*deg(T)-1` and delegate to the
ported ZX or Flx kernels. Generic mixed integer/polynomial products transpose by
inner degree as in PARI. Binary products use blocks of `2*deg(T)+1`; binary squares
square individual coefficients and insert zero odd powers. Normalization uses
native extended GCD with only the required Bézout output, retains scalar fast paths
and forces the original leading coefficient representation. These primitives are
dependencies for the extension-factor backend; they do not reroute Sage's existing
extension multiplication or complete extension-factor delegation by themselves.

### Extension-polynomial reciprocal contexts

Extension division preserves the native output-pointer modes internally, with
public quotient/remainder, remainder-only and quotient-only entry points. Shared
field callbacks retain integer/polynomial tags and delegate nonlinear arithmetic
to FpX/Flx/F2x; linear word operations reuse the existing reduced-coefficient
implementation. Normalization shares the same coefficient inverses.

The outer GEN reduction pair maps to an unambiguous object with polynomial and
inverse fields. Native get_red thresholds decide when to build it. Generic
word-prime division converts to FlxqX and drops the outer reciprocal as the
original does; its early remainder branch precedes that conversion. Plain inner
moduli remain the public input convention. Inner cached contexts for the complete
extension-factor backend remain separate work.

Reciprocal inversion uses native basecase/Newton cutoffs and quadratic precision
masks. Barrett division updates a backing array and logical length in place,
matching the original linear block traversal. It does not repeatedly copy the
remaining prefix. The binary remainder path preserves the bundled source's
distinct stale-length branch, documented in DEVIATIONS.md.

PARI extension GCD uses a dedicated coefficient-generic recursive matrix kernel.
Its basecase tracks V on every Euclidean step and derives U by quotient division,
matching all three extension backends; the existing prime-field kernel has
different schedules and is not reused for that step. Public extgcd returns
[gcd,U,V], and ExtensionMatrix<A> is [[A,A],[A,A]] in row-major order. GCD values
retain native scaling. Private V-only output avoids computing an unrequested U.

Binary extension products use the native XOR shift packing schedule. Coefficients
fitting individual blocks retain the existing linear string-packing path. Raw
overlapping blocks use 64-bit word arrays, avoiding quadratic whole-integer shifts;
the guard uses the bundled 64-bit allocation size before any write.

Extension arithmetic, division and GCD kernels carry an optional private inner
reciprocal array alongside the raw modulus. The caller prepares it at the native
get_red phase; a generic-to-word conversion converts an existing reciprocal
without allocating one when none existed. Cached coefficient reduction delegates
to the shared Barrett kernel with the matching FpX/Flx arithmetic callbacks.
Native error display uses a shared 1,600-character payload cutoff, with the GEN
type supplied separately from the rendered value.


### Extension composition and matrix packing

Generic and binary extension composition retain PARI's Brent–Kung block sizes
and reduction schedule. Word composition builds coefficient-block matrices and
uses the native Kronecker width derived from the prime, extension degree and
inner matrix dimension. Packed matrices delegate to ZM_mul, retaining its
classical, Strassen–Winograd and CRT routes. Private matrix helpers omit sentinel
entries; the public FlxqM_mul boundary adds the existing column/row zero slots.
The inner reciprocal travels with the private evaluation context, avoiding a
second reciprocal construction inside the matrix product.


### Public polynomial reciprocal preparation

The shared prime-field power/evaluation helpers distinguish a supplied prepared
reducer from the default public path. The public path creates the native inner
reciprocal at the original source stage; internal contexts retain their existing
cache. A private word flag selects Flx multiplication/squaring for word tables.
Canonicalizing numeric trailing zeros preserves native polynomial degrees while
leaving the first power's coefficient values unreduced.


Automorphism wrappers carry the prepared generic reciprocal through successive
compositions. Direct composition during automorphism squaring retains whole-word
dispatch, including its zero shortcut after reducing coefficients. Generic table
evaluation retains its separate matrix route. The automorphism-power seed table
is constructed from the original modulus, as in FpX.c; successive automorphism
and trace tables instead receive the existing prepared context.

### Coefficient substitution and binary composition

The six FpXY/FlxY/F2xY substitution entries share a private coefficient helper.
It preserves integer coefficient tags and passes polynomial coefficients to the
appropriate prime-field or binary table evaluator. Direct calls construct the
native Brent–Kung table from the inner degree and outer coefficient count before
substitution. Raw calls retain the original modulus for each coefficient; a
provided private reciprocal instead travels through the existing quotient context.

Packed binary scalar composition shares RgX.c's block schedule: small blocks are
reduced separately, giant steps multiply by the supplied table entries, and the
final sum is reduced. The existing binary power generator supplies the native
square-versus-multiply schedule. Public wrappers remain in FpXX.ts, FlxX.ts and
F2x.ts, matching the original dependency paths.


### Extension automorphism tuples and contexts

FpXX.ts, FlxX.ts and F2x.ts expose the eight original automorphism functions.
A private typed tuple kernel shares the source's window powering and Brent–Kung
composition schedules. Power pairs carry inner and outer images; generic/word
trace pairs carry only outer image and additive aggregate. Their autsum triples
multiply the aggregate, while binary trace triples add it. Public signatures
retain these distinct layouts. A prepared inner reducer and outer cache object
travel through coefficient substitution, quotient tables and composition without
rebuilding reciprocals. Source-order cache preparation precedes the unit count.
Signed native long counts use BigInt.asUintN(64, n) after checking signed bounds;
word auttrace instead checks unsigned bounds. Zero is outside the native kernel
contract and receives an explicit adapter error.


### Extension random projections and truncated products

The projection helper preserves native coefficient order and uses the existing
PARI random samplers. Generic dot products accumulate exact tagged coefficients
with ZX multiplication before one field reduction; word dots use Flx products
and modular additions before their final remainder. The existing extension
polynomial module shares its Kronecker pack/unpack routines with truncated
products. Generic products truncate packed coefficients before unpacking, while
word products retain the source's full multiplication then outer truncation.
This distinction preserves both scalar tags and inverse-error timing. Generic
integer-only squares intentionally retain the original nonmodular ZX route.


### Extension minimal-polynomial projection

The two mirrored wrappers share a typed Shoup kernel. It reuses extension random
sampling, delayed dot products and truncated multiplication for transposed power
projection; the recurrence is recovered with the original extension half-GCD.
Its monomial input preserves polynomial zero/one coefficient tags even in the
generic backend. One initial power table is retained across projection retries.
The outer cache is prepared with the original inner modulus; the power generator
prepares a local inner cache without replacing the caller's context. This source
phase differs from automorphism powering and is intentionally kept separate.


### Prime quotient inverse and power contexts

The existing ffinit/galconj compatibility exports delegate to a private shared
quotient helper. Generic inverse and word inverse preserve their distinct safe
scalar-inverse behavior and native error payloads. Signed powering prepares one
reducer after negative inversion, retaining it through the native exponent window.
The three-argument inverse remains the original FpXQ API; the explicit four-argument
convenience initializes through ZpXQ_inv's unsigned-word backend and retains its
existing Newton lift to the requested prime-power modulus. That convenience uses
successive modulus squaring capped at q; the native precision-mask schedule is
not claimed for this existing adapter. Both use logarithmically many lift stages.


### Extension Frobenius routing

The five FpXQX_factor.c wrappers and Flx.c prime Frobenius export share a private
phase dispatcher. Full Frobenius leaves the initial moduli unprepared while
computing the inner and outer p-th powers, then uses PARI's exact integer
composition threshold. Half Frobenius prepares inner then outer reciprocals and
passes both through powering and automorphism aggregation. Generic half Frobenius
converts cached moduli recursively at the unsigned-word boundary; generic full
Frobenius retains the existing signed-power entry and its adapter constraints.
Word characteristic-two Frobenius squares directly before requesting a cache.


### Extension split parts and outer derivatives

The root-count wrappers share one split-part helper, keeping the native
word-conversion-before-degree shortcut and delegating nonlinear work to Frobenius
and unscaled extension GCD. FqX_nbroots dispatches a null inner modulus to the
existing prime-field count. Outer derivatives share a tag-aware helper with no
inner reduction context, so coefficient representation and division-error order
remain visible. Squarefree wrappers evaluate that derivative before invoking GCD.


PARI FpX multiplication/division/GCD word routing follows integer storage width:
`lgefint(p) == 3` maps to `0 < abs(p) < 2^64`, with positive magnitude passed to
Flx. Apply this at the native dispatch point, after any earlier shortcuts; do not
normalize the public p argument globally, since generic constant inversions retain
its sign in errors. This rule is specific to the upstream limb-width test; APIs
using `is_bigint` or explicit domain guards retain their own dispatch contracts.


FpX quotient tables and composition contexts apply magnitude-based word routing
at each native phase, including conversion of supplied reciprocal coefficients.
FpM uses the same limb-width rule and its binary special case at abs(p)=2. Its
empty-dimension shortcuts run first; the generic route computes the integer
matrix product before coefficient reduction, preserving zero-product returns at
p=0. Representation guards must not replace such defined early-return behavior.


Polynomial observation helpers retain native operation ordering. Derivatives call
ZX_deriv then FpX_red; quotient-only linear division omits the unused remainder
calculation. Sparse FpX_eval skips coefficient gaps with logarithmic gen_powu_i
scalar powering and BigInt residue callbacks; generic exponent two squares the
raw integer first. Input GEN normalization is separate from output construction:
FpX_center must retain coefficients that become zero during its native comparison.


Exported Fp scalar arithmetic shares the polynomial kernels' residue and invmod
helpers. Do not identify all remainder operations: native Fp_sqr uses remii and
must reject a zero divisor even for a zero numerator. Fp_div follows its native
unsigned-divisor branch before generic inverse multiplication, including integer
exact division after an inverse modulo the small divisor. Keep native signed
outputs and shortcut order rather than applying one final field normalization.


Fp_pow delegates to the private _scalar_power helper, mirroring arith1.c's
word-power, unsigned-power and arbitrary-integer dispatch. Reuse bb_group's
window/fused schedules. Montgomery uses 64-bit magnitude limbs and Newton
inversion modulo 2^64; modified Barrett uses its native folded reciprocal and
bounded subtractions. Signed right shifts truncate magnitude toward zero, as
shifti does. Keep native final normalization separate from the reduction kernels.


Fp_order's native word branch normalizes the base and uses factoru; its generic
branch validates/factors the order and preserves the raw base for equal1 checks.
Do not move reduction or identity checks across this dispatch. The separate
znorder residue/modulus arguments model gmodulo before the native order routine;
format each coprimality operand independently through pariErrorPayload.


QPoly's num/den model corresponds to a canonical PARI rational polynomial.
QPoly_normalize returns positive minimal common denominators; coefficient export
trims polynomial zeros. QPoly_to_FpX must cancel each rational coefficient and
preserve Rg_to_Fp's zero-before-inverse phase. The Galois vectopol helper uses
polarit2.c's centermodii semantics, distinct from the raw Fp_center helper.


Galois integer helpers delegate to the source-based ispower and language/forprime
modules and the existing factoru backend. Word arithmetic remains exact until its
number-valued API return; modular word overflow precedes Number conversion. Prime
iteration keeps BigInt progress across rounded public results. Its initialized
500000 table limit and the Number projection are explicit compatibility choices.


Integer-polynomial Galois helpers now delegate squarefreeness to QX_factor and
p-adic denominator refinement to base2/hnf_snf. Those modules retain native
column-major matrix representation and word/generic arithmetic. The ifactor
strict-factor helper represents its output pointer as [factors, cofactor|null]
and specializes the initialized default factor limit; it is not a package-root
API. Modular integer GCD uses sequential CRT images with exact division checks.


The Zp dependency owns scalar Newton lifting, stable degree-based factor trees,
quadratic precision schedules, Bezout propagation and quotient Newton lifting.
Its vectors are zero-indexed; galconj preserves existing one-indexed public vectors.
Scalar FpX evaluation, ZX/FpX derivatives and FpX split-part live in their original
source modules, so Zp can delegate without importing the Galois facade. Explicit
precondition-guard transcripts are distinguished from native returned-value records.


Vandermonde interpolation lives in FpX alongside its product/remainder tree and
batch-inverse dependencies. Its matrices use native column storage; galconj alone
adapts them to one-indexed row storage. bb_group shares producttree_scheme between
this tree and gen_product, preserving the latter's breadth-first callback order.


Permutation exponents preserve the checked native unsigned/signed word conversions.
Order helpers keep BigInt state until their existing Number return. Unsigned LCM
wraps at each native word operation; signed zv_prod instead enforces its upstream
no-overflow precondition. The native transcript marks those precondition checks
separately and projects defined final words through binary64 for Number interfaces.


Symmetric-polynomial evaluation preserves PARI's scalar/column result distinction:
all-zero weights produce bigint zero, while nonzero weights produce a one-indexed
integer column, even when they cancel. Callers therefore narrow bigint[] | bigint.
The Galois fixed-field constructors already select nonzero separating weights and
retain their column invariant. Newton powers delegate to arith1.Fp_powu, sharing
the existing word, Montgomery and Barrett schedules rather than linear powering.


Unit-subgroup enumeration now delegates through the mirrored PARI files: char.ts
constructs flag-zero znstar cyclic factors and Conrey generators; subgroup.ts
runs the exact-index Birkhoff engine and its prime-to-p lifts; hnf_snf.ts computes
modular HNF; subcyclo.ts evaluates HNF generators and enumerates their cosets.
These internal matrices/vectors are zero-indexed. The existing Galois facade
retains one-indexed Number vectors. Modular HNF reduces coefficients eagerly
after elimination instead of testing PARI allocation-size thresholds; the same
column operations and canonical output are retained. Native packed subgroup
cells are round-tripped explicitly, including removal of initial zero limbs.


All three subgroup-bound modes share the Birkhoff engine. Scalar bounds first
restrict cyclic factors through ifactor1.absZ_factor_limit_strict; its unresolved
cofactor is a returned pair rather than a C output parameter. The default-limit
helper delegates to the same implementation. Word-prime table initialization is
lazy, allowing bounded factorization and the prime iterator to share dependencies
without an eager module-initialization cycle. Native configured table/primorial
boundaries remain fixed and are directly compared.


### PARI finite-field kernel dispatch

`alglin1.ts: FpM_ker` mirrors the native binary, ternary, word and generic-field
branches. Matrices retain the existing one-indexed column/row adapter. `F2v.ts`
handles binary columns and `F3v.ts` uses two-bit trits packed into BigInt columns.
`Flv.ts: Flm_ker` supplies word-field kernels. The existing row-oriented Galois
helper now only adapts its explicit row/column dimensions to this dependency.

`_matrix_kernel.ts` shares the native generic/word elimination schedules: Gaussian
base cases, recursive column echelon, and recursive unit-triangular solves. The
native cutoffs remain five for generic fields and eight for word fields. Generic
field additions and negations stay unreduced where the original bb_field does;
word operations reduce modulo p. This preserves signed generic-kernel output.
Matrix products delegate to the existing FpM_mul dispatch, retaining its packed
binary, word Winograd and arbitrary-integer multiplication implementations.

The ternary adapter represents a sequence of native packed limbs as one BigInt;
its two-bit addition/subtraction logic and pivot/basis order follow F3v.c. Matrix
indices are numbers; scalar field and packed arithmetic use BigInt throughout.


### PARI factor-stage result adapters

The existing SQUFOF/Pollard–Brent helpers return lists of factors whose product
is the input. Native internal stages instead return a scalar discovered factor
or triples of factor/exponent/status. Comparative adapters expand exponents and
append the residual cofactor, matching the TypeScript result contract; native
status pointers are not numerical outputs. The explicit Z_pollardbrent oracle
compares its native vector directly. Round/seed arguments retain their existing
number signatures; their integer arithmetic is performed with BigInt so the
represented input integers do not lose bits during counter shifts or addition.

Fixed-field prime-replacement comparisons convert native cached Galois groups
through gg_get_std, the same accessor used by their native caller, before
comparing generators and relative orders. Internal reconstruction caches are
not part of the TypeScript GaloisGens representation.


### Perfect-power native adapters (2026-09-12)

The comparison exposes root output pointers as tuples: native failure leaves a
pointer unwritten, so `Z_isanypower` and `is_pth_power` pair exponent zero with the
original input, and `is_357_power` pairs it with root zero and its updated mask.
The TypeScript prime iterator exposes `next()` and exhaustion as `null`; it has no
native `forprime_t.p` field. The oracle maps native `T.p` beyond the configured end
to zero for the test's last-consumed-prime observation, matching iterator
exhaustion. Cutoff primes within the configured interval remain observable.

The older finite-field factory adapter explicitly compares magnitude decomposition
on both sides (`pari(abs(n)).ispower()` / `Z_isanypower(abs(n))`). Signed public
PARI behavior is compared separately by the power-helper cases. Its established
expected transcripts remain unchanged.


### Native ECM storage and stages (2026-09-12)

`ifactor.ts:ellfacteur` creates the internal `_ecm.ts:ECM` state once, then invokes
its rounds with the existing TB1 schedules. The internal module ports batched
addition/doubling/inversion, the nine PRAC transformations and the native B2
helix and baby/giant steps. Products/inverses use the PARI BigInt field helpers.

A numeric address table points to mutable BigInt cells. Copying a coordinate
updates a cell (native affii); copying a GEN pointer aliases the address. This
retains overlapping scratch/baby-step storage and the native Xh/Yh offsets across
rounds. Treating the storage as independent point arrays would change subsequent
operations. Seed advancement is BigInt. The prime-table bound and PRAC heuristic
match the initialized native reference, as documented in DEVIATIONS.md.

Direct comparison cases include single rounds, reused states, the non-insisting
driver, and bounded prefixes of the insisting driver. The bounded native wrapper
copies the original ellfacteur schedule and only stops after the configured
number of complete rounds. No root package export changes are introduced.


### MPQS relation and word-root comparisons (2026-09-12)

The test-only mpqsInternals object exposes newHandle, mpqs_FB_ctor and
combine_large_primes so comparisons can construct a typed handle using the same
factor-base storage as the sieve. The relation fixtures use primes 2,3,5,7,11,13.
Class-group fixtures satisfy Y^2-D = 4*q*product(p_i^e_i); factor-mode fixtures
satisfy Y^2 = q*product(p_i^e_i) modulo N. Packed exponents/indices stay in the
range representable by both native and TypeScript relation words.

The Number word-square-root facade delegates to ff.Fl_sqrt. Its valid factor-base
inputs are reduced quadratic residues; -1 represents the native failure word
when stored in the signed int32 factor-base field. The native comparisons for
this helper use quadratic residues, through the signed int32 prime boundary.


### Sparse binary matrix backend (2026-09-12)

PARI F2Ms_ker belongs in F2v.ts; MPQS delegates to it. Sparse columns use one-based
row indices in zero-based outer arrays. The backend keeps the native 640-row
dispatch boundary, singleton elimination, unsigned 64-bit block words, eight
256-entry byte tables, and the global PARI RNG. Ordinary binary columns remain
BigInt bitsets with zero-based coefficient positions. The MPQS facade converts
these to its one-based 32-bit word layout.

The comparison dispatcher packs each sparse column as its length followed by
its row indices in a single flat list because shared fixedValue generators
accept flat integer arrays. Both runners decode that list before calling their
native/ported APIs. Native comparisons return the column-elimination map, ordered
kernel vectors and complete final random state. The TypeScript side also checks
the MPQS word adapter against the same backend result and seed.


### MPQS scalar helper delegation (2026-09-12)

MPQS's krouu/kroiu facades use ff.kronecker so the implementation retains native
word bits and signed residues. The debug relation factor-back helper uses Fp_neg,
Fp_mul and the audited signed Fp_pow dispatch; its machine-sized packed exponents
have the same results and inverse errors as native Fp_pows. The relation checker
uses Fp_mulu for the final modular product and PariError for native bug diagnostics.
These helpers remain on the test-only mpqsInternals object; package root exports
are unchanged. Native debug comparisons compile the original MPQS_DEBUG code and
capture its result/error while suppressing its verbose diagnostic stderr.


### MPQS initialization comparison boundary (2026-09-12)

The test-only mpqsInternals object exposes the existing parameter, factor-base,
sieve/poly constructor, threshold, prime-location, self-initialization and sieve
helpers. The comparison adapter follows the original class initializer's call
sequence; it does not add a class-group entry point to the production package.
Snapshots include factor-base primes/roots/flags, coefficients, selected primes,
start positions, scaled byte logs, every sieve byte, and ordered candidates.
Uninitialized native scratch fields are not read. Separate factor-base controls
exercise direct-factor returns and segmented-prime-table expansion. The candidate
scanner receives adequate extra capacity in both implementations.


### MPQS relation table (2026-09-12)

The relation table is now an incremental transcription of PARI's hash table
for the specific [integer, Vecsmall] relation-key shape. It compares colliding
keys by integer and word-array equality. Resizing happens before the insertion
that exceeds ceil(0.65*bucket_count); old buckets are visited in ascending order
and their chains are prepended to new buckets. Matrix columns therefore follow
the same traversal as hash_keys_GEN. Table capacity is retained independently
of later increases in the requested relation count. The ordinary large-prime
lookup table keeps its existing lookup behavior; it is not used for iteration.


### MPQS diagnostic comparison (2026-09-12)

The post-Gauss consistency check follows native pari_warn rather than raising an
exception. console.warn carries its core message in TypeScript. The comparison
adapter captures that call synchronously and restores the console afterward;
the C adapter forwards each captured warning to the original pari_warn and
records its count. Both compare the returned ordered factor pairs/null as well.

The test-only mpqsInternals object also exposes the existing mpqs_eval_cand and
mpqs_solve_linear_system helpers. Their full-relation table is RelationTable;
the large-prime lookup table remains a Map. Candidate controls reconstruct the
original class initialization sequence and enable its native debug checks.

### Generic Pollard lambda random stream (2026-09-12)

`discrete_log_lambda` delegates step-size draws to the existing CPython generator
owned by `current_randstate()`, matching Sage's `prandom` dependency. Custom hash
values are reduced to a nonnegative step-table index with BigInt arithmetic.
The shared group-area comparison records every hash input in order, so a changed
walk is detected even when both algorithms happen to return the same logarithm.


### Generic Pollard rho draw order (2026-09-12)

Rho coefficients and its initial exponent delegate to IntegerModRing.random_element,
matching the source dependency and its CPython random stream. All m coefficients
are sampled before any n coefficients. The shared oracle loads the bundled
original generic.py directly and supplies the documented TypeScript DJB2 hash;
it compares result/error and the following CPython random draw. This distinguishes
stream drift even when the discrete logarithm is unique.

### Shared group validation and exact order bounds (2026-09-12)

Generic group entry points validate standard/custom operation arguments through
one internal helper before algorithm-specific numeric checks. The custom binary
multiple path mirrors the source's operation schedule and idempotence shortcut.
Order bounds use BigInt quotient/remainder to implement exact signed ceil/floor;
the completed BSGS result delegates to order_from_multiple with its check disabled.
Comparative adapters load the bundled original and record callback arguments in
order, alongside results and errors. Existing standard element fast paths remain.

### Order-from-multiple valuation delegation (2026-09-12)

Nonempty prime lists call the existing arith.valuation dispatch, matching Sage's
M.valuation(p) and the underlying GMP factor-removal path. Empty factor/prime
lists request fallback factorization. The zero-multiple +Infinity sentinel follows
the source's single-base helper or its specific infinite-cost diagnostics.
SignError retains ArithmeticError inheritance and the native error name. This
local path does not supply a general symbolic infinity API. Oracle list entries
are Sage Integers, consistent with the existing IntegerLike value adapter.

### Generic group iterators and parent representations (2026-09-12)

A shared parent adapter accepts a stored parent or parent() method. The multiples
factory obtains standard identities through that parent when available, copies
both arguments immediately, and uses an explicit next-state transition. The
transition increments its index before the callback and assigns the next value
only after the callback succeeds, matching the original class after an exception.

The returned object is a generator with protocol methods implementing that state
machine; inherited iterator helpers remain available. return/throw close it under
JavaScript conventions. Copying honors __copy__/copy hooks and otherwise preserves
the prototype and shallow own descriptors. A hook is required for private/internal
state. The comparative fixture compares construction and each next step, including
copy counts, external mutations, callback errors and subsequent recovery.


### Number-field scalar actions and generic identity predicates

`NumberFieldElement.mul` accepts a field element or `RationalLike`, coercing the
latter through `toRational` and the existing exact `scalarMul` coefficient kernel.
This restores the scalar action used by additive generic-group calls without
adding repeated-addition loops. The existing number-field backend deviation
remains open. Number-field scalars are unbounded exact values; the bounded-ring
number coercion convention above does not apply to this signature.

Generic group identity checks accept Sage's `is_zero`/`is_one` methods and the
existing camel-case adapters. When neither exists, equality against the parent's
identity supplies the check. Standard group algorithm and backend choices are
otherwise retained.


### Number-field construction and scalar coercion

Number-field addition, subtraction, division and equality use the same
`RationalLike` scalar domain as multiplication. The equality overload for field
operands preserves structural compatibility with the existing generic GroupElement
interface; the union overload accepts scalar and union-typed callers.

`NumberField.__call__` retains its explicit JavaScript number input. Those values
follow the existing binary64 `RR.simplest_rational` kernel, including integral
values above 2^53; truncating a real or turning it directly into bigint does not
match Sage's QQ conversion. Exact bigint/Integer/Rational inputs remain exact.
Unembedded fields reject nonfinite reals with the native conversion diagnostic.

Quadratic scalar division uses rational inversion before coefficient scaling;
general fields coerce the scalar first. This distinction preserves the original
zero-division messages. The existing number-field arithmetic backend deviation
remains open; this repair does not claim full NTL kernel delegation.

### Number-field coefficient representation metadata

Coefficient indexing depends on Sage's quadratic element class. A module-private
WeakMap retains the original defining polynomial and lazily caches the equivalent
representation choice. This preserves the effect of original leading-coefficient
sign/scaling while the existing arithmetic kernel uses a monic polynomial.

The classifier mirrors the source's discriminant normalization with
`Integer.squarefree_part(bound=10000)`: GMP exact divisibility/division maps to
BigInt and the native wheel bounds are retained. Full factorization would choose
a different class when a square factor exceeds the bound. No public field member
or new exported helper is required. Integral JavaScript number indices adapt to
Python integers, consistently with existing array-index APIs; nonintegral and
nonfinite numbers exercise the original float-index branches.

### Integer/rational coercion and group completion

Integer-left addition, subtraction, multiplication and comparisons route rational
operands through the existing Rational kernels, mirroring Sage's coercion to QQ.
IntegerLike arithmetic retains its Integer return overload; rational operands
return Rational even when the denominator is one. Floor division is a separate
existing Integer API and is not part of this coercion change.

Generic GroupElement specifies equality against `this`, so concrete numeric
classes need not accept an arbitrary unrelated GroupElement. Standard inverse
dispatch recognizes Sage's `__invert__` spelling and the existing `inv` adapter.
When a host numeric class omits its parent, its neutral power/action supplies the
identity, as already supported by the multiples iterator.

Integer multiplicative powers and inverses can enter QQ. Overloads expose
`Integer | Rational` for standard multiplicative multiple calls and
`GroupOps<Integer | Rational>` for standard parseGroupOps with an Integer sample.
Closed-group calls retain their existing generic element type.

### Rational equality with host numbers

The explicit number overload of Rational.eq uses the existing RDF coercion and
comparison implementation, matching Sage's QQ/Python-float common parent. RDF
converts rationals through the audited nearest-binary64 kernel, then applies
native double equality. This preserves signed zero, underflow, overflow and NaN
behavior. Rational and IntegerLike operands continue through exact comparisons.


### Real-double scalar arithmetic and generic multiples

RDF arithmetic accepts the existing exact scalar types and explicit host
number/bool adapters. Scalar conversion delegates to the RDF constructor before
native double arithmetic. Existing RDF elements supply their stored double.
ModuleElement._add_long and _mul_long shortcuts map Python integer/bool zero
addition and one multiplication to object reuse; ZZ/QQ/RDF operands still follow
the ordinary arithmetic path. Concrete-this overloads preserve RingElement
structural compatibility while exposing scalar operands.

Generic multiple always follows Sage's binary algorithm, including standard
addition/multiplication. Native scalar multiplication and pow have different
rounding and exceptional-value behavior. This distinction does not apply to
Sage's _power_func: algorithms that use that helper intentionally use native
scalar multiplication/pow for standard operations.


### Generic group operation schedules and dictionary lookup

BSGS and order_from_bounds obtain the standard parent identity before their
bound processing and use generic multiple for their binary powers. BSGS preserves
inverse-before-power evaluation, reuses the final baby step to construct the giant
factor, and multiplies that factor on the left. has_order computes both recursive
powers eagerly before short-circuiting the recursive tests.

Discrete logarithms and order_from_multiple retain _power_func's native standard
powers. Discrete logarithms evaluate zero powers for identity tests and negative
powers directly, matching the source even for inexact arithmetic. BSGS buckets
retain the existing string hash adapter but verify object identity or element
equality, following Python dict collision semantics. Distinct NaNs are unequal.

Generic f-string errors adapt RDF's empty-format output for nan/inf/-inf; its
ordinary Sage representation remains NaN/+infinity/-infinity for percent-formatted
BSGS errors. This internal adapter does not expose a partial public __format__ API.


### Parsed parent identities through nested group algorithms

Order reduction compares elements to the parsed parent identity, matching Sage's
explicit equality tests; it does not substitute is_zero/is_one predicates.
Discrete logarithms resolve the target's parent identity before powers and pass
normalized operations to their nested BSGS calls. Bounds-based order searches
retain one parsed identity through exponential search and final order reduction.

Private helpers carry the parsed identity and the already selected power function.
This retains _power_func's native standard powers after the operation has been
normalized, while BSGS/scaling continue to use binary multiple. No function tags,
new public options or exported signature changes are needed. Existing adapters
for stored parent objects and parentless Integer/Rational values remain in effect.


### Pollard mutation hooks and stored-point equality

Both Pollard routines parse the target's parent before validation and preserve
native standard powers through normalized nested calls. When the base exposes
set_immutable, walk points receive that method call at the original source
positions. Lambda can freeze its input target, and a reused power can freeze the
base even before a singleton-interval error. The binary-matrix port's existing
mutability flag participates directly in this behavior.

Rho's four remembered points retain element references and exponent pairs inside
string-key buckets. Lookup verifies identity/equality; eviction removes only the
specific stored entry, preserving unrelated colliding keys. Walk partition hashes
remain the documented host adapters. Lambda preserves k=0 for singleton bounds:
it invokes the hash callback, raises the native zero-modulus error, and consumes
no random step-size draws.


### Number-field defining equations, reduction and powers

NumberField stores the original RationalPolynomial separately from its monic
arithmetic associate. Public polynomial getters and display use the original;
degree-one generators use the root of the monic associate. Element construction
packs rational coefficients into a common positive denominator and delegates
reduction to flint-ts/src/fmpq_poly/rem. That kernel delegates to FLINT's integer
pseudo remainder and divide-and-conquer pseudo division, with the native degree
splits, cutoff 16 and large-dividend threshold 128. Scalar GMP operations map to
BigInt and polynomial products use the existing FLINT multiplication port.

NumberFieldElement.pow accepts IntegerLike and calls the mirrored arith/power
Sage-element helper. Its least-significant-set-bit schedule deliberately differs
from groups.generic.multiple: no identity multiplication, exponent one returns
self, and subsequent set bits multiply the current square on the left.

Callers that require a monic arithmetic model normalize the defining polynomial
explicitly: quadratic unit conversion/sign/logarithms, the existing equation-order
different calculation and the number-field Frobenius filter. This preserves
support for scaled presentations while public getters retain their exact input.


### Number-field ideal trace data

`NumberField._pari_ideal_data()` caches the `NfIdealData` produced by PARI's
`nfmaxord_ideal_data` adapter from the existing maximal integral basis. The record
contains the integral multiplication table, primitive trace inverse and denominator,
scaled codifferent HNF, two-element codifferent representation and different HNF.
`base4.ts` consumes it for ideal HNF products and trace-dual inversion. These routines
use zero-indexed columns; `_matrix_inverse.ts` uses rows internally and converts at
matrix-product boundaries. Rational ideals carry a positive common denominator.

The Sage layer maps returned HNF columns back through its integral basis and caches
that HNF, bounding generator count by the degree during repeated squaring. The
shared `arith/power.ts` exports internal `generic_power_pos` for monoid callers that
handle zero and negative exponents themselves. These dependency helpers are exposed
through source subpaths; they are not new root package exports. Native buffer,
Hadamard-bound and higher-degree generator-search limitations are recorded in
DEVIATIONS.md under Number-field ideal backend adapters.


Ideal two-generator tuples contain two NumberFieldElement objects: even the rational
intersection generator belongs to the number field in Sage. The old bigint return
lost fractional denominators and is replaced in version 22.0.0. Integral ideal
bases use the same HNF-to-integral-basis mapping as ideal arithmetic. The zero ideal
maps to an empty basis. The port's existing free-module record caches that basis
and rank; it is not the generic FreeModule class (see DEVIATIONS.md).


Fractional ideal numerator and denominator are cached NumberFieldIdeal objects,
not scalar coefficient denominators. Their definitions follow Sage: D=(I+O_K)^-1
and N=I*D. Internal HNF records retain a separate positive bigint denominator.
`base4.idealdiv` takes rational HNF records, inverts the divisor first, and multiplies
through existing PARI ideal helpers. The Sage wrapper preserves single-generator
field division and the zero ideal's separate monoid dispatch. Field ideal factories
normalize empty/all-zero input to one zero generator without caching the ideal.

### PARI fast LLL stage and binary64 contraction profile

`parigp-ts/src/lll.ts` ports `lll.c`'s `fplll_fast` and iterative `Babai_fast`.
The internal adapter takes nonempty rectangular integer columns, copies the basis,
starts U at the identity (or omits U), and returns `[status, B, U]`. A status of -1
includes the original stage's partially reduced basis and transformation. This is
an uncertified stage; the heuristic, DPE, arbitrary-precision and FLATTER wrapper
stages remain separate audit work. It is not a replacement for `ZM_lll`.

Binary64 arithmetic is intentional where PARI itself uses C doubles. The bundled
native audit build is PARI 2.18.1, AArch64 Apple Clang 21, `-O3 -fno-strict-aliasing`.
That compiler contracts scalar dot-product and Gram-Schmidt tails into FMA, while
blocks of eight products are rounded separately. The dedicated square loop uses
FMA for every term after the first; the explicitly assigned final diagonal product
is rounded separately. The port spells out these operations to make that profile
stable across JavaScript runtimes. Changing contraction changes partial failure
states and sometimes whether this stage succeeds. `src/_binary64.ts` supplies
single-rounding `fma`, exact `ldexp` subnormal rounding and `frexp`; integer-to-double
conversion still delegates to PARI's existing `itor`/`rtodbl` ports.

AArch64's signed-word conversion saturates out-of-range doubles and maps NaN to zero;
its signed shift masks the count and wraps the word. These are explicit choices for
upstream undefined-C cases in partial failure states, not a portable C guarantee.
The native oracle compiles the unmodified original source, retains normal floating
contraction, and rejects other compiler architectures. Cross-profile representative
parity remains limited as described in DEVIATIONS.md.

### PARI DPE reduction state

The internal `fplll_dpe` adapter follows `lll.c:1430-2022`, preserving its double
significands with separate signed-word exponents, iterative Babai updates and exact
integer Gram updates. Exponents use bigint, so the normalized-zero sentinel
`-LONG_MAX` and adjacent values remain distinct. Immutable `{d, e}` records replace
native DPE output pointers. Integer-to-DPE and norm conversions delegate to the
existing PARI `itor`, `rtodbl` and `dbltor` ports.

The adapter copies optional B/G columns, starts U at the identity or null, and returns
`[status, G, B, U, norms]`. A null input Gram requests incremental construction;
B may be null if G is supplied. On incremental failure, returned G is null, matching
PARI; a supplied Gram retains its partial changes. Only the native updated triangle
is authoritative: `G[j][i]` for `i <= j`; the other entries may be stale. Optional
norms retain PARI real sign/exponent/mantissa/precision records. This is one reduction
stage, not the complete adaptive LLL wrapper.

### PARI real LLL stages and scalar dependencies

`fplll_heuristic` and `fplll` follow `lll.c:1129-1425` and `2025-2314` using the
existing PARI real arithmetic kernels. The heuristic stage has separate approximate
Gram and Gram-Schmidt precisions; the proved stage keeps an exact integer Gram.
Assignments copy and round into the destination precision, preserving distinct real
objects where native output pointers are distinct. Word-rounded coefficients retain
the low unsigned word, including wrap when rounding reaches 2^64. Larger coefficients
follow native mantissa truncation and power-of-two shifts. The native leading-word
and zero-tail test for the 3/2 shortcut is kept literally.

The bundled PARI 2.18.1 precision ABI counts bits (`prec2nbits` is the identity).
The real LLL adapters require positive multiples of 64, matching the adaptive native
caller's precision sequence and the port's word-sized real kernels. This differs
from the older installed Sage precision convention; native oracle arguments use the
bundled bit ABI directly. `src/gen3.ts` adds native `roundr_safe` (ties toward positive
infinity, with precision loss allowed); `kernel/none/cmp.ts` adds `abscmprr`, whose
zero comparisons ignore accuracy exponents. They are tested directly as well as
through the stages. The raw LLL boundary retains cypari2's full divrr-zero payload;
the lower-level typed divrr dependency retains its existing interface.


### PARI Householder QR dependencies

`bibli1.ts` ports QR_init, R_from_QR and gaussred_from_QR. Matrices contain
zero-indexed columns with bigint or MpReal cells; exact integer zeros remain distinct
from finite-accuracy real zeros. QR uses native mp arithmetic (including real zero
products), Householder sign selection, beta conversion at the working precision,
and the single-word squared-norm precision check. Real inputs retain their own
precision except where the native routine explicitly converts them.

The QR adapter returns `[success, B, Q, L]`, with null outputs on precision failure.
Q exposes the native n-1 public reflectors. For tall matrices the native final
reflector is scratch; its write lies outside the declared Q vector and overlaps B's
GEN header. The C oracle reads B using the known column count, while the TypeScript
adapter uses separate arrays and discards that extra reflector. This preserves the
numerical B/Q/L entries without reproducing allocator metadata corruption. R_from_QR
returns transposed L, and gaussred_from_QR returns the native reduction of x^T*x.


### PARI real Cholesky and triangular inverse dependencies

`alglin2.ts` ports the real-input branches of qfgaussred_positive and RgM_Cholesky;
`alglin1.ts` adds RgM_inv_upper. These zero-indexed column adapters consume immutable
MpReal records or exact gen_0, represented by 0n. Positive Gaussian reduction reads
the upper input triangle and retains the native row-update order; Cholesky keeps
real input precision when taking square roots. Its precision parameter is retained
for native caller compatibility, even though all nonzero cells are already real.
Triangular inversion performs native backward substitution, including products with
exact zero entries beyond the current inverse column.

`_real_matrix.ts` factors out the relevant gen1 numeric branches for MpReal/0n only.
Generic multiplication by exact zero returns exact zero (unlike QR's mp multiplication);
real zeros retain finite accuracy. Generic 0/0 originates in dvmdii, real/real zero
in divrr, and inverses distinguish unallocated real zeros from allocated zeros.
The Gen-boundary adapters render native real-zero payloads; typed lower-level kernels
keep their existing error conventions. A large allocated zero in native invr's Newton
branch reads unspecified payload words, so the existing deterministic port rejection
remains explicit and is excluded from exact native value comparison counts.


### PARI integer and real matrix products

`RgV.ts` implements the integer/real numeric branches of RgV_dotsquare,
RgV_dotproduct, gram_matrix and RgM_mul. Matrices use zero-indexed columns. Dot
products retain native sequential summation and scalar pointer-identity squaring;
generic matrix multiplication skips exact-zero left entries after the first term.
Entirely integer products delegate to ZV.ZM_mul, converting at the existing adapter's
unused-slot-zero boundary. Real/integer products keep native mulir rounding and
exact-zero behavior instead of first converting every coefficient to real.

`polarit2.RgM_rescale_to_int` scans real mantissas for their lowest nonzero bit and
integer coefficients for their native exponent, then scales and rounds. In mixed
matrices integer entries can become fractions during scaling and are rounded too;
the name does not guarantee preservation of every exact integer ratio. All-integer
inputs are copied unchanged. With only zero coefficients the native minimum exponent
remains HIGHEXPOBIT: shifting a negative-accuracy real zero then raises expo overflow.
That boundary is retained explicitly without materializing a huge shift.


### PARI adaptive Gram-Schmidt precision

`lll.ts` re-exports the FLATTER dependencies implemented in `_lll_gso.ts`: drop,
potential, spread, condition_bound, GS_extraprec, gramschmidt_upper,
gramschmidt_dynprec and RgM_Cholesky_dynprec. Exponent statistics use bigint so the
native integer-zero HIGHEXPOBIT sentinel and neighboring values remain distinct.
Real zeros retain their accuracy exponent for gexpo; absolute diagonal comparison
uses the native abscmprr dependency. The supported coefficient domain is integer/real.

Dynamic QR/Cholesky use native initial bit counts, failure doubling and the maximum
of 4/3 growth versus the condition/spread bound, rounding allocations to 64-bit words.
The upper-triangular shortcut converts directly at its selected precision. Arithmetic
continues to delegate to the separately audited QR and real Cholesky implementations.
These routines require a full-column-rank basis or a positive-definite Gram matrix;
they are dependency stages, not a complete LLL wrapper. A checked allocation precision
must fit the exact JavaScript integer range. Returned records retain native precision.


### PARI modular matrix pivots and solves

`Flv.ts` adds Flm_pivots and Flm_gauss; `alglin1.ts` adds ZM_pivots, ZM_rank and
ZM_gauss. These follow the existing dependency convention of an unused slot zero on
both matrix axes. Nonnull public pivot vectors also have an unused first slot; pivot
values are native one-based row indices, with zero for dependent columns. Word
operations consume prime-field coefficients. ZM_gauss returns integer numerator
columns and a positive common denominator, rather than per-entry rational GENs.

Shared row-oriented kernels in `_matrix_inverse.ts` preserve small Gaussian
elimination and the cutoff-eight CUP decomposition, using the existing recursive
triangular solvers and matrix products. Integer pivots try the native small-word prime
sequence and certify nonmaximal modular ranks with exact rational solves. Integer
solves use the big-word prime sequence, doubling CRT batches, a balanced CRT tree,
column-wise rational reconstruction with denominator reuse and exact verification.
Products delegate to ZM_mul. Native worker batches execute sequentially in TypeScript;
the balanced arithmetic and certified stopping conditions remain intact.

Native integer rank can raise `inconsistent dimensions in gauss` when the first trial
primes divide every entry. This source behavior is preserved, including its recursive
zero-rank certificate path, and has exact comparison cases. It is not replaced by an
uncertified modular rank or a local rational elimination fallback.

### Adaptive PARI LLL and FLATTER

`src/lll.ts` re-exports `_lll_wrapper.ts` for the full integer wrapper and typed
integer/real lllfp adapter. These adapters use zero-indexed columns, native numeric
flag values and an untagged matrix or [kernel,image] result (`LllResult`). Direct
flat returns [basis, transformation|null, drop, potential], with bigint statistics.
Stage adapters retain their existing identity-initialized transformation interface;
the wrapper composes each returned transformation, including failed-stage partial
states. Exact matrix products delegate to the existing PARI integer backend.

FLATTER follows native recursive block sizes, overlap, rank augmentation, fixed
threshold tables and precision schedules. qfb.ts retains its row-oriented matrices,
so the wrapper transposes redimagsl2's transformation once at that boundary. The
new redimagsl2 entry is native-compatible even where the existing public qfbredsl2
intentionally corrects an upstream large-input sign bug. BigInt mantissas are
multiplied before their power-of-two shifts in all four reduction stages, preserving
native arithmetic complexity. Input mutation and numeric/resource restrictions are
tracked in DEVIATIONS.md under `PARI adaptive LLL and FLATTER adapters`.


### Number-field ideal construction and intersection

The field factory first constructs a NumberFieldFractionalIdeal and maps ValueError
to a fresh base-class zero ideal, following Sage. A same-field fractional ideal
passes through by identity. Direct constructors accept an array, including a single
nested generator list; the fractional constructor tests zero before field coercion.
The exported NumberFieldIdealGenerator and NumberFieldIdealInput aliases live in
rings/number_field/number_field_ideal.ts and its number_field barrel. Rational
constants coerce across abstract fields; general embedding discovery remains open.

The base4.idealintersect dependency takes integral-basis HNF columns and positive
common denominators. It removes common content, cross-scales, calls ZM_lll with
LLL_KER, truncates the kernel, multiplies and computes modular HNF before restoring
the least denominator. The Sage layer only coerces inputs and converts the returned
HNF into the appropriate ideal class. Tests compare canonical power-basis lattices
to avoid conflating installed-Sage versus bundled-source integral-basis conventions.


### Ideal operator coercion and centered maximal-order bases

NumberFieldInput extends scalar field construction with exact-length arrays of
rational-coercible coefficients. NumberField.__call__ reduces those coefficients
through QQ semantics before constructing the power-basis element. Ideal membership
catches TypeError from field coercion; vector-length ValueError propagates. Addition,
multiplication, divisibility, coprimality and division preserve their distinct Sage
entry paths rather than sharing a blanket same-field check. Native bigint maps to
Python int, Integer to Sage ZZ and host number to Python float, including their
zero-division fallback differences.

The final maximal-order basis now uses PARI ZM_hnfcenter. The shared gen3.diviiround
helper rounds half ties toward positive infinity; centering processes rows from
bottom to top and changes only the preceding coordinates of later columns. The
base4.idealmul adapter removes primitive rational content before choosing the
smaller ideal for two-generator reduction, then restores the numerator and least
denominator. Its zero-indexed column representation matches idealdiv/intersection.


The existing bounded class-number certificate normalizes ideal generators to
nonnegative and centered power-basis HNFs before searching each small coefficient
box. This makes its candidate sets independent of the ambient integral basis and
preserves the exact norm/membership certificate. Native bnfinit remains unported.


### Native byte streams and copy assignment

NTL's RandomStream is exposed through its mirrored ZZ.ts dependency module.
Uint8Array replaces byte pointer/length pairs; returned byte arrays are independent
outputs. Uint32Array represents the ChaCha kernel's 32-bit state words. Counts are
exactly represented integer numbers, while the nonce uses bigint and wraps through
the native unsigned 64-bit conversion.

A RandomStream copy constructor duplicates the buffer and position. Its
assign(other): this method represents C++ copy assignment: it updates the receiver
in place and returns that same receiver, so JavaScript aliases observe the update.
Ordinary JavaScript assignment still aliases objects. This explicit method is
necessary to preserve native operator= state transitions without operator overloading.


### Explicit NTL integer-sampling contexts

The dependency's integer sampling helpers take a RandomStream argument in place
of NTL's implicit thread-local LocalGetCurrentRandomStream. This follows the
existing explicit-modulus adapters: callers supply the complete context whose
state the operation reads and updates. It does not introduce an ambient RNG or
claim that FFT-table initialization has been incorporated into random state.
C++ long and ZZ bounded-sampling overloads both return bigint in TypeScript;
RandomBnd defaults to the ZZ algorithm and options.word selects the signed-long
algorithm. Keeping that choice explicit is necessary because the two algorithms
consume different bytes, including for equal small bounds.


### Cached NTL sampler ownership

RandomBndGenerator uses the same explicit RandomStream context as the other
integer samplers. Constructing from another generator copies its cached bound
and shares its stream reference; it does not clone the stream. assign(other)
represents native copy assignment, updating the existing receiver and returning
it so JavaScript aliases observe the copied bound and stream. A null bound
constructs an uninitialized generator tied to its supplied explicit context.
A subsequent successful build initializes that bound. Native implicit context
capture is represented by the constructor's required stream argument.


### NTL word-context propagation

Word polynomial APIs retain explicit bigint p and accept a final optional
zz_pXOptions ({maxroot?: number}) where native context settings select algorithms.
zz_pXModulus captures this value independently of the caller's options object.
Its requested maxroot is retained separately from effective MaxRoot: clamping the
request before choosing primes or passing it to nested routines would change
PrimeCnt and algorithm cutovers. Functions taking F use that captured context;
raw-modulus recursive helpers forward the options. This avoids ambient mutable
modulus state while representing ZZXFactoring's degree-dependent initialization.
The interface is exported from the existing lzz_pX deep module. Native contexts'
CRT storage and global FFT/RNG initialization remain outside this adapter.


### NTL integer-factorization local information

LocalInfoT represents the native static helper's retained information. Its p and
pattern accessors return independent snapshots; setters copy the visible prefix.
Private initialized storage and separate logical lengths preserve NTL Vec's
shrink/grow behavior without scanning or copying the whole retained capacity on
each factorization call. Unknown native long slots correspond to undefined array
slots, not invented numeric values. Integer degree/count fields use exactly
represented numbers; primes and PossibleDegrees use bigint. The s member owns a
PrimeSeq, and context explicitly carries the last selected/restored p and requested
maxroot. SmallPrimeFactorization returns the native nullable factor array and
mutates this information; optional tuning values replace thread-local globals.


### NTL local-information suffix-cache replacement

UpdateLocalInfo in the mirrored integer-factorization module mutates LocalInfoT
and returns bigint[] | null. The array is the native suffix-cache replacement;
null means no recomputation occurred, so a caller uses `updated ?? previous`.
An empty array remains a valid replacement. This maps the native output parameter
without exposing an ambient vector; failed output-parameter writes are discarded,
while native partial mutations of retained LocalInfo remain observable. Temporary
word contexts preserve the saved explicit context object. The existing native
tuning globals map to operation-local integer options with their original defaults.


### NTL cardinality-search outputs and arithmetic profiles

CardinalitySearch and CardinalitySearch1 return the native three output parameters
as [bigint[][], bigint[], bigint[][]], with independent coefficient arrays.
LocalInfoT remains the retained mutable state. FindTrueFactors composes these
entries and returns factors in native recovery order. Explicit modulus/options
replace ambient integer-modulus context and tuning globals. Unsigned filters use
64-bit bigint wraparound; the original floating logarithms and ceilings are
retained solely for pruning-table sizing. Successful filter caches and the native
operation counter are updated in source order, including intentionally retained
f(1) and ratio values after factor removal. Failed output vectors are discarded,
while LocalInfo mutations preceding an exception remain observable.


### NTL probable-prime generation

ZZ.ts exposes the native prime-generation helpers as free functions, with bigint
values and explicit RandomStream arguments. Native long dimensions/trial counts
use exactly representable number parameters; optional trial/error values become
options objects. Bound helpers preserve source floating-point estimates, while
all modular integer arithmetic remains exact. The 256-bit candidate-block path
uses private streams derived from the same 256-bit seed and nonce schedule,
which preserves saved/restored caller streams without ambient TLS or worker
threads. All final witnesses are sampled before primality decisions. Number of
workers is a performance detail of this adapter's tested single-thread profile.


### NTL word matrix output parameters

The matrix linear-algebra dependency accepts rectangular bigint rows and an
explicit word modulus. WordMatrixOptions carries columns for zero-row inputs;
WordMatrixEliminationOptions adds pivot width w. WordMatrixInverseOptions carries
relax and previous matrix output. WordMatrixSolveOptions carries relax, left and
previous vector output. Status overloads return determinant/output tuples, while
strict inv returns the inverse or throws. Inputs and retained previous outputs
are copied, including coefficient normalization. Matrix rows do not store empty
column metadata; callers retain it separately. The native exact scalar Gaussian
schedule is used at all sizes with the same cubic arithmetic complexity as packed
and blocked implementations. No field pivot is normalized in gauss/image, and
kernel computes the native left nullspace rather than transposing a right basis.


### NTL FFT-prime cache state

FFTPrimeContext holds the original descending-candidate state and a lazy array of
immutable FFTPrimeInfo objects. Number indices preserve exact native signed values;
moduli and roots use bigint. IsFFTPrime's output parameter becomes a returned
[status, rootOrPrevious] tuple, and NextFFTPrime returns [modulus, root]. Explicit
RandomStream objects account for initialization and retries. Cache entries expose
q, qrecip, RootTable and TwoInvTable; word preconditions and native large-table
acceleration metadata are omitted because arithmetic uses BigInt. Root vectors
are deeply frozen, retaining native const-entry identity across growth. The
native comparison oracle isolates complete global cache lifetimes per trace and
serializes floating reciprocal bits to avoid textual float-format differences.


### NTL incremental CRT and arbitrary-modulus determinants

NTL scalar and matrix CRT output references become [modified, residue, modulus]
tuples. Integer and modular matrices remain arrays of BigInt rows; explicit
column counts retain zero-row native shapes. Matrix CRT precomputes one inverse
for all entries, preserving the native complexity and positive-half tie rule.
Arbitrary-modulus determinant follows mat_ZZ_p.cpp's delayed-reduction Gaussian
schedule. The modulus is explicit, and native GMP arithmetic primitives map to
exact BigInt operations. DetBound preserves the original row-norm bound and bit
length, including empty and zero-row cases. Native storage preallocation and
thread scheduling are omitted; they do not change these deterministic results.


### NTL integer matrix reconstruction state

Integer determinant/inverse reconstruction receives the same FFTPrimeContext and
RandomStream used by FFT-prime initialization and prime generation. The port keeps
the native stabilized-CRT and probabilistic/full-bound stopping schedules, including
the order of all random draws. Modular determinant/inverse operations receive an
explicit modulus. Status inverse returns [determinant, adjugateOrPrevious], while
the strict overload returns an integral inverse and preserves NTL's determinant
1/-1 restriction. Exact matrix-product certification follows the native cubic
algorithm. Prior status outputs are copied, and zero-dimensional results retain
the established array/column-options convention.


### NTL certified integer row reduction

The integer factorization gauss helper returns [signedPivotDeterminant, scaledRREF]
and shares FFTPrimeContext/RandomStream with reconstruction. It delegates to the
word-prime sampler, word Gaussian elimination and integer status inverse. Ordinary
word-context initialization is represented by its actual FFT-prime cache growth;
machine-word preconditioners are unnecessary for exact BigInt arithmetic. The
original retry order and exact row-form certification remain unchanged. Matrices
use the existing BigInt row-array and optional-column conventions.

### Modular polynomial root parents

`Polynomial.roots` has residue-coefficient overloads because Sage's default
Zmod(p) roots change coefficient parent: the result is a list of
`[FiniteFieldElement, number]` pairs from the cached `IntegerModRing.field()`.
The concrete field class is the prime-field constructor's `FiniteFieldPrime`.
The distinct option returns `IntegerMod[]` in the original residue ring. The
`this: Polynomial<IntegerMod & RingElement>` overload expresses that change
without widening roots over other coefficient classes. The intersection retains
the existing polynomial interface's `this`-returning arithmetic constraint;
it does not introduce another runtime element class.

The polynomial method delegates to `IntegerModRing._roots_univariate_polynomial`.
The ring's `field()` and `factored_order()` caches use WeakMaps keyed by the parent,
retaining result identity without adding private structural members to its type.
Factorizations use the existing bigint-pair array representation. Sage's static
residue-root lifting helper is exposed as a TypeScript static class method.
See [Modular Polynomial Roots and Hensel Lifting](DEVIATIONS.md#modular-polynomial-roots-and-hensel-lifting)
for the native category adapter and the explicitly retained upstream linear bug.

### Generic curve isomorphism ordering and parent adapters

The generic curve's tuple-returning `isomorphism_to`, `isomorphisms` and
`automorphisms` APIs delegate to `weierstrass_morphism._isomorphisms` and its
native morphism comparison key. The first-isomorphism and boolean predicates
consume only one generator result. Lists construct morphisms, sort with the
source key and expose their tuples through the existing API.

`WeierstrassIsomorphism._comparison_impl` keeps a boolean overload for valid
morphism operands. An unknown-operand overload returns `boolean | null`, with
`null` representing Python's `NotImplemented`, as in `_composition_impl`.
String names and punctuation map to the six native comparison opcodes.

The internal `_same_base_ring` helper in `schemes/elliptic_curves/types.ts`
emulates factory-parent equality across the supported finite-field classes.
It first checks identity; prime implementations compare characteristics and
chosen generators (including custom degree-one moduli), while
extension definitions compare characteristic, degree, generator name and modulus.
Equivalent parents' coefficients are explicitly converted before curve arithmetic
or tuple comparison, so the dedicated GF2 implementation can interoperate with
other prime-field implementations. This does not introduce global parent interning.


## Finite-field scalar square-root output pointers

`FF_issquareall(x: PariFfelt): PariFfelt | null` represents PARI's call with a
non-null root output pointer. A root object represents success; null represents
failure. It delegates to the n=2 specialization of the bundled FF_ispower route,
including the word/large-prime high-extension algorithm and binary inverse
exponentiation. The omitted-output-pointer predicate is a separate deterministic
FF_issquare dependency, not inferred from this randomized root call.

Extension elements expose `sqrt({extend?: boolean, all?: boolean})`. Literal
all=true returns an array; omitted/false returns an element; a runtime boolean
returns their union. The default extend=false matches element_pari_ffelt.
Roots retain their input field parent. The existing shared PARI random stream
is preserved; no extra Sage random seeding is added by this scalar method.


## Finite-field norm and square-predicate dependencies

`FF_issquare(PariFfelt): boolean` is the deterministic, no-root-output predicate.
`FF_norm(PariFfelt): bigint` returns the native prime-field residue. Polynomial
kernels expose FpX/Flx resultants, FpXQ/Flxq norms and their square predicates.
Resultants share the existing half-GCD matrix recursion with an optional private
accumulator; calls without it keep the existing exported basecase dispatch.

The Sage extension `is_square(): boolean` delegates to FF_issquare. This does
not change the separate Sage `norm()` method: the bundled inherited method
selects a characteristic-polynomial coefficient, so repairing that method
requires the characteristic-polynomial dependency rather than guessing a
shortcut to FF_norm. Existing norm/trace caller audits remain separate.


### Generic curve coordinate coercion

Generic curve `lift_x` distinguishes canonical parent maps from element
construction. The finite-field subset follows the native `_coerce_map_from_`
relations; it does not infer a map by attempting to convert one element. The
existing structural parent-equality helper represents equal explicit finite
parents across uncached TypeScript constructors. Unknown parent families are
outside this subset (see DEVIATIONS.md).

The point overload keeps `F` for base-field/integer inputs and returns `F | G`
for another field's element, accounting for curve promotion. Runtime input
validation is available through an unknown-input overload. Curve base change
passes original coefficient elements to the target constructor and returns the
original curve for an unchanged parent, matching Sage's cached curve factory.
