# Agent Guidelines

This document provides instructions for AI agents working on sagemath-ts.

## Required Reading

Before working on this codebase, read these documents:

- **DESIGN.md** - Design decisions for the SageMath to TypeScript port (type mappings, function signatures, dependency architecture)
- **DEVIATIONS.md** - Documented differences from SageMath behavior
- **SCOPE.md** - Module implementation status and assignments

---

## Project Overview

We are porting SageMath to TypeScript with a focus on number theory and cryptography. The goal is **exact behavioral equivalence** with SageMath for deterministic functions.

## Key Principles

### 1. Mirror SageMath Structure Exactly

- **File paths must match**: `sage/rings/integer.py` -> `src/rings/integer.ts`
- **Function names must match**: Use the exact same names as SageMath
- **Module hierarchy must match**: Preserve the import structure

### 2. Reference the Source

Before implementing anything:
1. Read the SageMath source in `reference/sage/src/sage/`
2. Read relevant dependency source (PARI/GP, FLINT, NTL) in `reference/`
3. Understand the algorithm completely before writing TypeScript
4. Check for tests in the reference code and replicate them

### 3. TypeScript Style

See **DESIGN.md** for complete details on type mappings and conventions. Key points:

```typescript
// Accept IntegerLike, return bigint
import { IntegerLike, toBigInt } from '../types/coercion.js';

function gcd(a: IntegerLike, b: IntegerLike): bigint {
  const _a = toBigInt(a);
  const _b = toBigInt(b);
  // ... implementation
}

// Use options objects for keyword arguments
function factor(n: IntegerLike, options?: { algorithm?: 'pari' | 'flint' }): Factorization { ... }

// Preserve SageMath's error messages
throw new ValueError("n must be positive");
```

### 4. Property Testing Requirements

Every implemented function MUST have a property test:

```
tests/property/
  python/rings/test_integer.py    # Python/SageMath test
  typescript/rings/test_integer.ts    # TypeScript test
```

Both must use identical random seeds and output results in identical format.

### 5. Document Deviations (MANDATORY)

When your implementation differs from SageMath, you **MUST** document it appropriately:

| Type of Difference | Document In | Purpose |
|-------------------|-------------|---------|
| **Behavioral differences** (outputs differ from SageMath) | `DEVIATIONS.md` | Track when results are different |
| **Architectural decisions** (type patterns, conventions) | `DESIGN.md` | Explain how we map concepts |

**DEVIATIONS.md entries require:**
1. What SageMath does vs what we do
2. **Rationale** - Why we made this choice
3. **Trade-offs** - What we lose
4. **Behavioral impact** - Does it affect outputs?

**Rules:**
- `DEVIATIONS.md` at the project root is the single source of truth
- Same-change requirement: Code changes that introduce deviations must update `DEVIATIONS.md` in the same commit
- Add `@see Deviation:` in affected docstrings

### 5b. Keep `LLM.md` True (MANDATORY)

`LLM.md` is the public API quick reference — it is what downstream agents and vendored
consumers read instead of the source. If you change an exported signature, a method name,
an import path, or the package `exports` map, update `LLM.md` **in the same commit**.

Its examples are executed by `tests/llm-doc.test.ts`. Never fix a failure there by
loosening the assertion: either the doc is stale (fix the doc) or the change was
unintentionally breaking (fix the code). Only document behavior you have actually run.

### 6. Scope Tracking (MANDATORY)

**Update `SCOPE.md` after completing any work.**

| When | Action |
|------|--------|
| Starting work | Mark as 🟡 with your identifier |
| Completing work | Mark as ✅ with test coverage |
| Blocked | Mark as 🔴 with reason |

### 7. Stub Unimplemented Functions

Stub ALL functions first with `NotImplementedError`:

```typescript
export function unimplemented(n: IntegerLike): bigint {
  throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: unimplemented');
}
```

Find unimplemented functions: `grep -r "SAGE_NOT_IMPLEMENTED" packages/`

---

## Architecture Fidelity

**When SageMath delegates to an external library, we MUST also delegate to our port of that library.**

See **DESIGN.md** for the complete dependency mapping. Before implementing, check if SageMath calls:
- `__pari__()`, `pari(...)` -> Use `parigp-ts`
- `flint_...`, `fmpz_...` -> Use `flint-ts`
- `ntl_...` -> Use `ntl-ts`

---

## Current Focus: Number Theory for Cryptography

Priority modules (implement in this order):

1. **`sage.rings.integer`** - Arbitrary precision integers
2. **`sage.rings.finite_rings`** - Finite fields (GF(p), GF(p^n))
3. **`sage.arith`** - Basic number theory (gcd, lcm, factor, primality)
4. **`sage.rings.polynomial`** - Polynomial arithmetic
5. **`sage.groups.generic`** - Generic group operations
6. **`sage.schemes.elliptic_curves`** - Elliptic curve operations

---

## Workflow for New Module

1. **Check SCOPE.md** - Ensure the module isn't already assigned/complete
2. **Update SCOPE.md** - Mark as 🟡 in progress with your identifier
3. **Study the SageMath source** - Read `reference/sage/src/sage/<path>`
4. **Create mirrored file structure** in `packages/sagemath-ts/src/`
5. **Implement with tests** - Write property tests alongside implementation
6. **Run transcript comparison** - `bun run test:property`
7. **Update SCOPE.md** - Mark as ✅ complete with test coverage

---

## Algorithm Fidelity (CRITICAL)

**NEVER write naive O(n) implementations when SageMath uses O(√n) or O(log n) algorithms.**

Before implementing ANY function:

1. **Read the SageMath source** to understand the algorithm used
2. **Check the complexity** - if SageMath uses BSGS, Pohlig-Hellman, factorization-based methods, etc., we must too
3. **Check for delegation** - if SageMath calls PARI/FLINT/NTL, we delegate to our ports

### Common Algorithm Patterns to Watch For

| Operation | WRONG (naive) | RIGHT (SageMath's approach) |
|-----------|---------------|----------------------------|
| Element order in group | O(n) repeated multiplication | O(√n) BSGS via `order_from_bounds` |
| Discrete logarithm | O(n) brute force | Pohlig-Hellman + BSGS |
| Verify exact order | Just check `n*P = O` | Also check `(n/p)*P ≠ O` for all prime divisors |
| Point order on E/F_q | Compute in TypeScript | Delegate to PARI's `ellorder` |
| Curve cardinality | Naive point counting | Delegate to PARI's `ellcard` (Schoof-Elkies-Atkin) |
| Factorization | Trial division only | Delegate to PARI's `Z_factor` |

### Red Flags in Code

If you see any of these patterns, STOP and check SageMath:

```typescript
// 🚫 BAD: O(n) loop for order computation
while (!current.is_zero()) {
  current = current.add(this);
  n++;
}

// 🚫 BAD: Simple divisibility check for has_order
has_order(n) { return this.mul(n).is_zero(); }

// 🚫 BAD: Brute force enumeration
for (let i = 0n; i < groupOrder; i++) { ... }
```

### Delegation Architecture

```
SageMath                    Our Port
────────                    ────────
sage.groups.generic    →    src/groups/generic.ts (BSGS, Pohlig-Hellman)
cypari2.ellorder()     →    parigp-ts/ellorder()
cypari2.ellcard()      →    parigp-ts/ellcard()
cypari2.factor()       →    parigp-ts/Z_factor()
```

When implementing a function that SageMath delegates:
1. First check if the dependency function exists in our port (parigp-ts, flint-ts, etc.)
2. If not, implement it there first
3. Then have sagemath-ts delegate to it

---

## Avoiding Common Mistakes

- **Don't guess algorithms** - Always verify against SageMath source
- **Don't write naive implementations** - Check SageMath's algorithm complexity first
- **Don't skip edge cases** - SageMath handles many edge cases; we must too
- **Don't change function signatures** - Even if TypeScript conventions differ
- **Don't use floating point** - Use BigInt and rational arithmetic
- **Don't implement without tests** - Property tests are mandatory
- **Don't restrict input types** - Use `IntegerLike` not just `bigint` (see DESIGN.md)
- **Don't implement what PARI/FLINT provides** - Delegate to our ports instead

---

## Dependency Libraries

| SageMath uses | We implement in |
|---------------|-----------------|
| PARI/GP (via cypari2) | `packages/parigp-ts/` |
| FLINT | `packages/flint-ts/` |
| NTL | `packages/ntl-ts/` |
| GMP | Native BigInt + `packages/gmp-ts/` if needed |

---

## Testing Commands

```bash
# Run all tests
bun test

# Run property tests with transcript comparison
bun run test:property

# Run tests for specific module
bun test --filter "rings/integer"

# Generate coverage report
bun test --coverage
```

---

## Versioning and Changelog

- Keep versions in sync across `package.json` and all `packages/*/package.json`.
- Semantic versioning: patch for fixes/docs/audits, minor for new features, major for breaking changes.
- Update `CHANGELOG.md` for every change (code, tests, or docs) with date and concise bullets.
- Version bump and changelog update must be in the same commit as the changes.

---

## Git Commits

- Always commit your work when done
- Use one-line commit messages
- Do not add co-author attribution

---

## Questions?

If unclear about implementation details:
1. Check SageMath documentation: https://doc.sagemath.org/
2. Check the source in `reference/sage/`
3. Run the operation in SageMath to observe behavior
4. Document any ambiguities in the code comments
