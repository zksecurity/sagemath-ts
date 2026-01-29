# SageMath to TypeScript Porting Guide

This guide establishes best practices for porting SageMath modules to TypeScript with high fidelity.

## Table of Contents

1. [Documentation Standards](#documentation-standards)
2. [Type Compatibility](#type-compatibility)
3. [Property Testing with fast-check](#property-testing-with-fast-check)
4. [Deviation Tracking](#deviation-tracking)
5. [Porting Checklist](#porting-checklist)

---

## Documentation Standards

### JSDoc Template (Matching SageMath Structure)

Every function should use this documentation pattern that mirrors SageMath's docstring structure:

```typescript
/**
 * Brief one-line description.
 *
 * Extended description with mathematical explanation. Use LaTeX-style
 * notation where helpful: computes gcd(a, b) such that gcd(a,b) | a and gcd(a,b) | b.
 *
 * ## Input
 *
 * - **a** - `bigint`: First integer
 * - **b** - `bigint` (optional): Second integer. If omitted, `a` should be an array.
 *
 * ## Output
 *
 * - `bigint`: The greatest common divisor, always non-negative
 *
 * ## Algorithm
 *
 * Uses Stein's binary GCD algorithm which avoids division operations.
 * Complexity: O(log(min(a,b))²)
 *
 * @example
 * ```typescript
 * // Basic usage
 * gcd(12n, 8n)  // => 4n
 *
 * // Negative numbers
 * gcd(-15n, 25n)  // => 5n
 *
 * // List form
 * gcd([12n, 18n, 24n])  // => 6n
 *
 * // Edge cases
 * gcd(0n, 5n)  // => 5n
 * gcd(0n, 0n)  // => 0n
 * ```
 *
 * @throws {ValueError} If inputs cannot be coerced to integers
 *
 * @see Reference: sage/arith/misc.py:gcd (lines 1234-1290)
 * @see SageMath docs: https://doc.sagemath.org/html/en/reference/rings_standard/sage/arith/misc.html#sage.arith.misc.gcd
 * @see Deviation: Uses Stein's algorithm instead of Euclidean (DEVIATIONS.md#algorithm-implementation-choices)
 */
export function gcd(a: bigint, b?: bigint): bigint;
export function gcd(values: bigint[]): bigint;
export function gcd(a: bigint | bigint[], b?: bigint): bigint {
  // Implementation
}
```

### Key Documentation Elements

| Element | Required? | Purpose |
|---------|-----------|---------|
| Brief description | ✅ Yes | First line, explains what function does |
| Extended description | If complex | Mathematical context |
| `## Input` section | ✅ Yes | Document all parameters |
| `## Output` section | ✅ Yes | Document return value |
| `## Algorithm` section | If non-obvious | Explain approach |
| `@example` | ✅ Yes | At least 3 examples including edge cases |
| `@throws` | If can throw | Document exceptions |
| `@see Reference:` | ✅ Yes | Link to SageMath source with line numbers |
| `@see Deviation:` | If applicable | Link to DEVIATIONS.md if behavior differs |

### When to Link to Deviations

Add a `@see Deviation:` link when:
- Algorithm differs from SageMath
- Return type differs (e.g., `null` vs exception)
- Edge case handling differs
- Performance characteristics differ significantly

---

## Type Compatibility

### SageMath to TypeScript Type Mapping

| SageMath Type | TypeScript Type | Notes |
|---------------|-----------------|-------|
| `Integer` / `ZZ` | `bigint` | Native arbitrary precision |
| `Rational` / `QQ` | `Rational` class | Custom implementation |
| `int` (Python) | `bigint` | Always use bigint for consistency |
| `bool` | `boolean` | Standard |
| `list` | `T[]` | Generic arrays |
| `tuple` | `[T, U, ...]` | Fixed-length tuples |
| `dict` | `Map<K, V>` or `Record<string, V>` | Prefer Map for non-string keys |
| `None` | `null` | Not `undefined` |
| `str` | `string` | Standard |
| `float` | `number` | Avoid in cryptographic code |
| `GF(p).element` | `FiniteFieldElement` | Custom class |
| `Polynomial` | `Polynomial<T>` | Generic over coefficient ring |
| `EllipticCurve` | `EllipticCurve` | Custom class |
| `FiniteField` | `FiniteFieldPrime` or `FiniteFieldExtension` | Depends on field |

### Return Type Conventions

```typescript
// SageMath: Returns value or raises ValueError
// TypeScript: Return value or null (document in @throws if can throw for other reasons)
function sqrt_mod(a: bigint, p: bigint): bigint | null;

// SageMath: Returns (g, s, t) tuple
// TypeScript: Return tuple type
function xgcd(a: bigint, b: bigint): [bigint, bigint, bigint];

// SageMath: Returns True/False
// TypeScript: Return boolean
function is_prime(n: bigint): boolean;

// SageMath: Can return value or (value, proof)
// TypeScript: Use overloads
function is_square(n: bigint): boolean;
function is_square(n: bigint, root: true): [boolean, bigint];
function is_square(n: bigint, root?: boolean): boolean | [boolean, bigint];

// SageMath: Returns iterator
// TypeScript: Return Generator
function* primes(start: bigint, stop?: bigint): Generator<bigint>;

// SageMath: Returns Factorization object
// TypeScript: Array of [prime, exponent] pairs
type Factorization = Array<[bigint, bigint]>;
function factor(n: bigint): Factorization;
```

### Handling Optional Parameters

```typescript
// SageMath: def func(a, b=None, algorithm='default'):
// TypeScript: Use options object for multiple optional params
interface FunctionOptions {
  b?: bigint;
  algorithm?: 'default' | 'fast' | 'proven';
}

function func(a: bigint, options?: FunctionOptions): bigint;

// Or for simple cases, use optional parameters directly
function next_prime(n: bigint, proof?: boolean): bigint;
```

---

## Property Testing with fast-check

### Setup

```bash
bun add -d fast-check
```

### Writing Property Tests

Create property tests in `tests/property/typescript/` using fast-check:

```typescript
// tests/property/typescript/arith.property.test.ts
import { describe, it, expect } from 'bun:test';
import * as fc from 'fast-check';
import { gcd, lcm, xgcd, is_prime, factor } from '@sagemath-ts/sagemath-ts/arith';

/**
 * Property tests for arithmetic functions.
 * These verify mathematical invariants and compare against SageMath.
 */
describe('arith property tests', () => {
  // ============================================
  // GCD Properties
  // ============================================

  describe('gcd', () => {
    it('gcd(a, b) === gcd(b, a) [commutativity]', () => {
      fc.assert(
        fc.property(
          fc.bigInt({ min: 0n, max: 10n ** 20n }),
          fc.bigInt({ min: 0n, max: 10n ** 20n }),
          (a, b) => gcd(a, b) === gcd(b, a)
        )
      );
    });

    it('gcd(a, gcd(b, c)) === gcd(gcd(a, b), c) [associativity]', () => {
      fc.assert(
        fc.property(
          fc.bigInt({ min: 1n, max: 10n ** 10n }),
          fc.bigInt({ min: 1n, max: 10n ** 10n }),
          fc.bigInt({ min: 1n, max: 10n ** 10n }),
          (a, b, c) => gcd(a, gcd(b, c)) === gcd(gcd(a, b), c)
        )
      );
    });

    it('gcd(a, b) divides both a and b', () => {
      fc.assert(
        fc.property(
          fc.bigInt({ min: 1n, max: 10n ** 15n }),
          fc.bigInt({ min: 1n, max: 10n ** 15n }),
          (a, b) => {
            const g = gcd(a, b);
            return a % g === 0n && b % g === 0n;
          }
        )
      );
    });

    it('gcd(a, 0) === abs(a) [identity]', () => {
      fc.assert(
        fc.property(
          fc.bigInt({ min: -10n ** 15n, max: 10n ** 15n }),
          (a) => gcd(a, 0n) === (a < 0n ? -a : a)
        )
      );
    });
  });

  // ============================================
  // Extended GCD Properties
  // ============================================

  describe('xgcd', () => {
    it('xgcd returns Bézout coefficients: g = s*a + t*b', () => {
      fc.assert(
        fc.property(
          fc.bigInt({ min: 1n, max: 10n ** 12n }),
          fc.bigInt({ min: 1n, max: 10n ** 12n }),
          (a, b) => {
            const [g, s, t] = xgcd(a, b);
            return g === s * a + t * b;
          }
        )
      );
    });

    it('xgcd(a, b)[0] === gcd(a, b)', () => {
      fc.assert(
        fc.property(
          fc.bigInt({ min: 0n, max: 10n ** 12n }),
          fc.bigInt({ min: 0n, max: 10n ** 12n }),
          (a, b) => xgcd(a, b)[0] === gcd(a, b)
        )
      );
    });
  });

  // ============================================
  // Prime Testing Properties
  // ============================================

  describe('is_prime', () => {
    it('small known primes return true', () => {
      const knownPrimes = [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 97n, 101n];
      for (const p of knownPrimes) {
        expect(is_prime(p)).toBe(true);
      }
    });

    it('small known composites return false', () => {
      const composites = [4n, 6n, 8n, 9n, 10n, 12n, 15n, 21n, 25n, 100n];
      for (const n of composites) {
        expect(is_prime(n)).toBe(false);
      }
    });

    it('0 and 1 are not prime', () => {
      expect(is_prime(0n)).toBe(false);
      expect(is_prime(1n)).toBe(false);
    });

    it('products of two primes are not prime', () => {
      const smallPrimes = [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n];
      fc.assert(
        fc.property(
          fc.integer({ min: 0, max: smallPrimes.length - 1 }),
          fc.integer({ min: 0, max: smallPrimes.length - 1 }),
          (i, j) => {
            const p = smallPrimes[i];
            const q = smallPrimes[j];
            if (p === q) return true; // p^2 needs separate test
            return !is_prime(p * q);
          }
        )
      );
    });
  });

  // ============================================
  // Factorization Properties
  // ============================================

  describe('factor', () => {
    it('product of factors equals original', () => {
      fc.assert(
        fc.property(
          fc.bigInt({ min: 2n, max: 10n ** 10n }),
          (n) => {
            const f = factor(n);
            const product = f.reduce((acc, [p, e]) => acc * p ** e, 1n);
            return product === n;
          }
        )
      );
    });

    it('all factors are prime', () => {
      fc.assert(
        fc.property(
          fc.bigInt({ min: 2n, max: 10n ** 8n }),
          (n) => {
            const f = factor(n);
            return f.every(([p, _e]) => is_prime(p));
          }
        )
      );
    });

    it('exponents are positive', () => {
      fc.assert(
        fc.property(
          fc.bigInt({ min: 2n, max: 10n ** 8n }),
          (n) => {
            const f = factor(n);
            return f.every(([_p, e]) => e > 0n);
          }
        )
      );
    });

    it('factors are in ascending order', () => {
      fc.assert(
        fc.property(
          fc.bigInt({ min: 2n, max: 10n ** 8n }),
          (n) => {
            const f = factor(n);
            for (let i = 1; i < f.length; i++) {
              if (f[i][0] <= f[i - 1][0]) return false;
            }
            return true;
          }
        )
      );
    });
  });
});
```

### Automatic Regression Test Generation

When a property test finds a counterexample, save it as a regression test:

```typescript
// tests/property/regressions/arith.regressions.ts
/**
 * Regression tests auto-generated from property test failures.
 * DO NOT EDIT MANUALLY - these are generated by the test runner.
 */
import { describe, it, expect } from 'bun:test';
import { gcd, factor } from '@sagemath-ts/sagemath-ts/arith';

describe('arith regressions', () => {
  // Generated: 2024-01-15 from property test "gcd commutativity"
  // Shrunk counterexample
  it('regression #1: gcd with large numbers', () => {
    const a = 12345678901234567890n;
    const b = 98765432109876543210n;
    expect(gcd(a, b)).toBe(gcd(b, a));
  });

  // Add more regressions as discovered
});
```

### fast-check Reporter for Auto-Saving Failures

```typescript
// tests/property/utils/regression-reporter.ts
import { writeFileSync, readFileSync, existsSync } from 'node:fs';
import { join } from 'node:path';

interface RegressionEntry {
  id: string;
  date: string;
  property: string;
  counterexample: unknown[];
  module: string;
}

const REGRESSIONS_FILE = join(import.meta.dir, '..', 'regressions', 'registry.json');

export function saveRegression(
  module: string,
  property: string,
  counterexample: unknown[]
): void {
  const existing: RegressionEntry[] = existsSync(REGRESSIONS_FILE)
    ? JSON.parse(readFileSync(REGRESSIONS_FILE, 'utf-8'))
    : [];

  const entry: RegressionEntry = {
    id: `${module}-${Date.now()}`,
    date: new Date().toISOString(),
    property,
    counterexample,
    module,
  };

  existing.push(entry);
  writeFileSync(REGRESSIONS_FILE, JSON.stringify(existing, null, 2));

  console.log(`\n⚠️  Regression saved: ${entry.id}`);
  console.log(`    Property: ${property}`);
  console.log(`    Counterexample: ${JSON.stringify(counterexample)}`);
}

// Use in property tests:
// fc.assert(fc.property(...), {
//   reporter: (result) => {
//     if (result.failed) {
//       saveRegression('arith', 'gcd commutativity', result.counterexample);
//     }
//   }
// });
```

---

## Deviation Tracking

### When to Document Deviations

**ALWAYS** document in `DEVIATIONS.md` when:

1. **Algorithm differs** - Different approach than SageMath
2. **Return type differs** - `null` vs exception, tuple vs object
3. **Error handling differs** - Different exception types or messages
4. **Edge case behavior differs** - Empty inputs, zero, negative numbers
5. **Performance differs significantly** - O(n²) vs O(n log n)
6. **Features omitted** - Options or modes not implemented
7. **Features added** - Extra functionality not in SageMath

### Where to Reference Deviations

1. **In function docstring**: Add `@see Deviation:` link
2. **In README.md**: Mention DEVIATIONS.md exists
3. **In AGENTS.md**: Already documented ✓
4. **In error messages**: Reference deviation when throwing

```typescript
// In docstring
/**
 * @see Deviation: Returns null instead of ValueError (DEVIATIONS.md#error-handling)
 */

// In code
if (!result) {
  // Deviation: SageMath raises ValueError, we return null
  // See DEVIATIONS.md#error-handling---null-returns-vs-exceptions
  return null;
}
```

---

## Porting Checklist

Use this checklist when porting a new function:

### Before Implementation

- [ ] Read SageMath source thoroughly
- [ ] Note what libraries SageMath calls (PARI, FLINT, etc.)
- [ ] Check SCOPE.md - mark as 🟡 in progress
- [ ] Identify edge cases from SageMath tests/doctests

### Implementation

- [ ] Create stub with `NotImplementedError` first
- [ ] Implement with correct TypeScript types
- [ ] Match SageMath function signature exactly
- [ ] Handle all edge cases from SageMath

### Documentation

- [ ] Write full JSDoc with Input/Output/Algorithm sections
- [ ] Include 3+ examples covering normal and edge cases
- [ ] Add `@see Reference:` with file and line numbers
- [ ] Add `@see Deviation:` if behavior differs

### Testing

- [ ] Port SageMath doctests to `sage-doctests.test.ts`
- [ ] Write property tests with fast-check
- [ ] Add regression tests for any found counterexamples
- [ ] Run comparison tests against SageMath transcripts

### Finalization

- [ ] Update SCOPE.md - mark as ✅ complete
- [ ] Document any deviations in DEVIATIONS.md
- [ ] Run full test suite: `bun test`
- [ ] Run property tests: `bun run test:property`

---

## Example: Complete Port of a Function

Here's a complete example showing all standards:

```typescript
// packages/sagemath-ts/src/arith/misc.ts

import { NotImplementedError, ValueError } from '../errors.js';

/**
 * Return the integer square root of n.
 *
 * The integer square root of n is the largest integer r such that r² ≤ n.
 * This is equivalent to floor(sqrt(n)) but computed exactly using integer
 * arithmetic only.
 *
 * ## Input
 *
 * - **n** - `bigint`: A non-negative integer
 *
 * ## Output
 *
 * - `bigint`: The largest integer r with r² ≤ n
 *
 * ## Algorithm
 *
 * Uses Newton's method (Heron's method) with integer arithmetic.
 * Convergence is quadratic, typically requiring O(log log n) iterations.
 *
 * @example
 * ```typescript
 * // Perfect squares
 * isqrt(0n)   // => 0n
 * isqrt(1n)   // => 1n
 * isqrt(4n)   // => 2n
 * isqrt(100n) // => 10n
 *
 * // Non-perfect squares (floor)
 * isqrt(2n)   // => 1n
 * isqrt(99n)  // => 9n
 * isqrt(101n) // => 10n
 *
 * // Large numbers
 * isqrt(10n ** 100n) // => 10n ** 50n
 * ```
 *
 * @throws {ValueError} If n is negative
 *
 * @see Reference: sage/arith/misc.py:isqrt (lines 2890-2920)
 * @see SageMath docs: https://doc.sagemath.org/html/en/reference/rings_standard/sage/arith/misc.html#sage.arith.misc.isqrt
 */
export function isqrt(n: bigint): bigint {
  if (n < 0n) {
    throw new ValueError('isqrt requires a non-negative integer');
  }
  if (n === 0n) return 0n;
  if (n === 1n) return 1n;

  // Newton's method: x_{n+1} = (x_n + n/x_n) / 2
  let x = n;
  let y = (x + 1n) / 2n;

  while (y < x) {
    x = y;
    y = (x + n / x) / 2n;
  }

  return x;
}
```

And the corresponding property test:

```typescript
// tests/property/typescript/arith.property.test.ts

describe('isqrt', () => {
  it('isqrt(n)² ≤ n < (isqrt(n)+1)²', () => {
    fc.assert(
      fc.property(
        fc.bigInt({ min: 0n, max: 10n ** 30n }),
        (n) => {
          const r = isqrt(n);
          return r * r <= n && n < (r + 1n) * (r + 1n);
        }
      )
    );
  });

  it('isqrt(n²) === n for perfect squares', () => {
    fc.assert(
      fc.property(
        fc.bigInt({ min: 0n, max: 10n ** 15n }),
        (n) => isqrt(n * n) === n
      )
    );
  });

  it('isqrt is monotonic', () => {
    fc.assert(
      fc.property(
        fc.bigInt({ min: 0n, max: 10n ** 20n }),
        fc.bigInt({ min: 0n, max: 10n ** 20n }),
        (a, b) => {
          if (a <= b) {
            return isqrt(a) <= isqrt(b);
          }
          return isqrt(b) <= isqrt(a);
        }
      )
    );
  });
});
```
