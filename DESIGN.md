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

See **DEVIATIONS.md** section "No JavaScript Number Coercion" for full rationale.

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
- Easy to wrap if needed: `ZZ(result)` or `new Integer(result)`

### Options Objects

Replace Python's keyword arguments with TypeScript options objects:

```python
# SageMath
def factor(n, algorithm='pari', proof=None, limit=None):
    ...
```

```typescript
// TypeScript
interface FactorOptions {
  algorithm?: 'pari' | 'flint';
  proof?: boolean;
  limit?: number;
}

function factor(n: IntegerLike, options?: FactorOptions): Factorization {
  const { algorithm = 'pari', proof, limit } = options ?? {};
  // ...
}
```

### Overloads

When SageMath has different behavior based on argument types, use TypeScript overloads:

```typescript
// Different return type based on input
function sqrt(n: IntegerLike): bigint | null;
function sqrt(n: IntegerLike, all: true): bigint[];
function sqrt(n: IntegerLike, all?: boolean): bigint | null | bigint[] {
  // Implementation
}
```

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

SageMath uses `Ring.__call__()` for element creation and coercion. We implement this as callable classes:

```typescript
// Create ring instance
const F = GF(7);  // Finite field of order 7

// Ring call creates/coerces elements
const a = F(3);   // Create element 3 in F
const b = F(10);  // Coerces 10 to 3 (mod 7)
```

Implementation uses a callable pattern:

```typescript
class FiniteField {
  // Make instances callable
  static create(p: bigint): FiniteFieldCallable {
    const ring = new FiniteField(p);
    const callable = (x: IntegerLike) => ring.element(x);
    return Object.assign(callable, ring);
  }
}
```

### Automatic Coercion

Elements coerce automatically in arithmetic:

```typescript
const F = GF(7);
const a = F(3);

// These all work:
a.add(F(4))     // Element + Element
a.add(4)        // Element + number (coerced)
a.add(4n)       // Element + bigint (coerced)
```

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
// TypeScript - must use explicit methods
const F = GF(7);
const a = F(3);
const b = a.add(F(5));  // Must explicitly wrap 5
// OR use the coercion in the method:
const b = a.add(5);  // Only works if add() accepts number

const R = PolynomialRing(ZZ, 'x');
const x = R.gen();
const f = x.pow(2n).add(R(3n));  // Must explicitly construct
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
const f = x.pow(2n).add(R(QQ(1n, 2n)));  // Explicit construction needed
```

### Matrix Entry Coercion

Same limitation as polynomials - matrix entries must be explicitly constructed or passed through matrix constructors that handle coercion internally:

```python
# SageMath
M = Matrix(GF(7), [[1, 2], [3, 4]])  # Integers auto-coerced
```

```typescript
// TypeScript
const F = GF(7);
const M = Matrix(F, [[F(1), F(2)], [F(3), F(4)]]);  // Explicit
// OR use constructor that accepts raw values:
const M = Matrix(F, [[1, 2], [3, 4]]);  // If constructor handles coercion
```

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

We accept `number` in `IntegerLike` for convenience with small values, but `toBigInt()` throws `RangeError` if the number exceeds `Number.MAX_SAFE_INTEGER`:

```typescript
gcd(12, 8)      // OK - small numbers
gcd(12n, 8n)    // OK - bigint (recommended for large values)
gcd(9007199254740992, 2)  // RangeError - exceeds safe range
```

**Recommendation**: Use bigint literals (suffix `n`) for any values that might be large. This ensures precision and avoids runtime errors.

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
