# Cryptography Tutorial

A hands-on course teaching the number theory foundations of modern cryptography using TypeScript and the `sagemath-ts` library.

## Structure

The course is organized into 6 parts:

### Part 1: Foundations of Modular Arithmetic
- **01-integers-and-divisibility** - Divisibility, GCD, Euclidean algorithm
- **02-modular-arithmetic** - Congruences, modular operations
- **03-groups** - Group theory basics, cyclic groups
- **04-rings-and-fields** - Rings, fields, primitive roots

### Part 2: Core Number-Theoretic Results
- **05-fermats-little-theorem** - Fermat's theorem and applications
- **06-eulers-theorem** - Euler's totient and theorem
- **07-chinese-remainder-theorem** - CRT and applications
- **08-quadratic-residues** - Legendre/Jacobi symbols, square roots

### Part 3: Computational Problems
- **09-primality-testing** - Fermat test, Miller-Rabin
- **10-factorization** - Pollard's rho, p-1 algorithms
- **11-discrete-logarithm** - DLP, baby-step giant-step

### Part 4: Public-Key Cryptography
- **12-rsa** - RSA encryption and signatures
- **13-diffie-hellman** - DH key exchange, ElGamal

### Part 5: Elliptic Curve Cryptography
- **14-elliptic-curves-intro** - EC basics, point addition
- **15-elliptic-curves-finite-fields** - EC over finite fields
- **16-ecdh-and-ecdsa** - ECDH, ECDSA

### Part 6: Advanced Topics
- **17-pairings** - Bilinear pairings, BLS signatures
- **18-lattices-intro** - Lattice basics, LLL

## Running Tutorials

### Run a single tutorial file

```bash
bun tutorial/part1-foundations/01-integers-and-divisibility/gcd.ts
```

### Run all tutorials

```bash
bun tutorial/run-tutorials.ts
```

### Run a specific part

```bash
bun tutorial/run-tutorials.ts part1-foundations
```

### Run a specific chapter

```bash
bun tutorial/run-tutorials.ts part1-foundations 01-integers-and-divisibility
```

## Testing

Run all tutorials as tests:

```bash
bun test tutorial/tutorial.test.ts
```

## Prerequisites

- Bun 1.x
- The `sagemath-ts` package (installed in the monorepo)

## File Format

Each tutorial file follows a literate programming style:

```typescript
/**
 * # Title
 *
 * Mathematical explanation...
 *
 * ## Background
 * ...
 */

import { ... } from 'sagemath-ts';

// Example 1: Description
console.log("=== Example 1 ===");
// ... code ...

// Example 2: Description
console.log("=== Example 2 ===");
// ... code ...
```

The tutorials are designed to be both:
1. **Readable** - Clear explanations in comments
2. **Runnable** - Execute to see concepts in action
