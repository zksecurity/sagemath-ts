# TypeScript Style Guide

This guide ensures consistency across the sagemath-ts codebase and alignment with SageMath's structure.

## General Principles

1. **Match SageMath exactly** - Function names, file structure, behavior
2. **Explicit over implicit** - Use explicit types, avoid `any`
3. **No floating point** - Use BigInt for integers, rationals for fractions
4. **Preserve error semantics** - Match SageMath's error messages and types

## File Organization

### Mirroring SageMath Paths

```
SageMath:                           TypeScript:
sage/rings/integer.py        →      packages/sagemath-ts/src/rings/integer.ts
sage/rings/finite_rings/     →      packages/sagemath-ts/src/rings/finite_rings/
sage/arith/misc.py           →      packages/sagemath-ts/src/arith/misc.ts
```

### File Structure Template

```typescript
/**
 * @module sage/rings/integer
 * @description Arbitrary precision integer arithmetic
 *
 * Port of: sage/rings/integer.pyx
 * Reference: reference/sage/src/sage/rings/integer.pyx
 */

// Imports from our packages
import { ZZ } from './integer_ring';

// Imports from dependency packages
import { fmpz } from '@sagemath-ts/flint-ts';

// Types
export type Integer = bigint;

// Classes (match SageMath class names)
export class Integer {
  // ...
}

// Functions (match SageMath function names)
export function gcd(a: Integer, b: Integer): Integer {
  // ...
}
```

## Naming Conventions

### Match SageMath Names

```typescript
// SageMath: def is_prime(n, algorithm='BPSW'):
export function is_prime(n: Integer, algorithm: 'BPSW' | 'miller_rabin' = 'BPSW'): boolean;

// SageMath: def power_mod(a, n, m):
export function power_mod(a: Integer, n: Integer, m: Integer): Integer;

// SageMath: class EllipticCurve_finite_field
export class EllipticCurve_finite_field extends EllipticCurve_generic;
```

### TypeScript Adaptations

When SageMath uses Python-specific features, adapt as follows:

```typescript
// Python: **kwargs
// TypeScript: options object
function factor(n: Integer, options?: {
  algorithm?: 'pari' | 'flint';
  limit?: Integer;
  proof?: boolean;
}): Factorization;

// Python: *args
// TypeScript: rest parameters
function gcd(...args: Integer[]): Integer;

// Python: property
// TypeScript: getter
class Integer {
  get is_unit(): boolean { return this.abs() === 1n; }
}
```

## Types

### Core Numeric Types

```typescript
// Arbitrary precision integer - use native BigInt
type Integer = bigint;

// Rational number
interface Rational {
  numerator: Integer;
  denominator: Integer;
}

// Finite field element
interface FiniteFieldElement {
  value: Integer;
  parent: FiniteField;
}
```

### Avoid These

```typescript
// BAD: Never use number for mathematical values
function gcd(a: number, b: number): number;

// BAD: Never use any
function factor(n: any): any;

// BAD: Never use implicit any
function helper(x) { return x + 1; }
```

## Error Handling

### Match SageMath Errors

```typescript
// SageMath: raise ValueError("n must be positive")
throw new ValueError("n must be positive");

// SageMath: raise TypeError("...")
throw new TypeError("...");

// SageMath: raise ZeroDivisionError
throw new ZeroDivisionError("rational division by zero");
```

### Custom Error Classes

```typescript
// errors.ts
export class ValueError extends Error {
  name = 'ValueError';
}

export class ZeroDivisionError extends Error {
  name = 'ZeroDivisionError';
}

export class NotImplementedError extends Error {
  name = 'NotImplementedError';
}
```

## Documentation

### JSDoc Comments

```typescript
/**
 * Return the greatest common divisor of a and b.
 *
 * The result is always non-negative.
 *
 * @param a - First integer
 * @param b - Second integer
 * @returns The GCD of a and b
 *
 * @example
 * ```typescript
 * gcd(12n, 8n)  // 4n
 * gcd(-4n, 6n)  // 2n
 * ```
 *
 * @see Reference: sage/arith/misc.py:gcd
 */
export function gcd(a: Integer, b: Integer): Integer {
  // ...
}
```

## Testing

### Property Test Correspondence

Every function needs corresponding test files:

```
packages/sagemath-ts/src/arith/misc.ts
    ↓
tests/property/python/arith/test_misc.py
tests/property/typescript/arith/test_misc.ts
```

### Test Output Format

Both Python and TypeScript tests must output in the same format:

```
# Operation: gcd(12, 8)
# Result: 4
# Operation: gcd(-4, 6)
# Result: 2
```

## Algorithms

### Reference the Source

Always document which algorithm is being used:

```typescript
/**
 * Primality test using BPSW (Baillie-PSW).
 *
 * Algorithm:
 * 1. Trial division up to small bound
 * 2. Miller-Rabin with base 2
 * 3. Strong Lucas test
 *
 * Reference: sage/arith/misc.py:is_pseudoprime
 * See also: reference/pari/src/basemath/ifactor1.c
 */
export function is_prime_bpsw(n: Integer): boolean {
  // ...
}
```

### Match SageMath's Algorithm Choices

When SageMath offers multiple algorithms, implement all of them:

```typescript
export function factor(n: Integer, options: {
  algorithm?: 'pari' | 'flint' | 'ecm' | 'qsieve';
} = {}): Factorization {
  const algo = options.algorithm ?? 'pari';

  switch (algo) {
    case 'pari':
      return factor_pari(n);
    case 'flint':
      return factor_flint(n);
    // ...
  }
}
```

## Module Exports

### index.ts Pattern

Each directory should have an index.ts that mirrors SageMath's `__init__.py`:

```typescript
// packages/sagemath-ts/src/rings/index.ts
export * from './integer';
export * from './integer_ring';
export * from './rational';
export * from './rational_field';
export * from './finite_rings';
```

### Main Package Export

```typescript
// packages/sagemath-ts/src/index.ts
// Mirror sage's top-level imports
export { ZZ } from './rings/integer_ring';
export { QQ } from './rings/rational_field';
export { GF } from './rings/finite_rings';
export { PolynomialRing } from './rings/polynomial';
export { EllipticCurve } from './schemes/elliptic_curves';
export * as arith from './arith';
```
