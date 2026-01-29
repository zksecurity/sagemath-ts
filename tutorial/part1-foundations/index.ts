/**
 * # Part 1: Foundations of Modular Arithmetic
 *
 * This tutorial series covers the mathematical foundations essential for
 * understanding cryptography. Each chapter builds on the previous ones,
 * progressing from basic number theory to abstract algebra.
 *
 * ## How to Use These Tutorials
 *
 * Each file is a self-contained, runnable TypeScript program that teaches
 * a specific concept. Run them with:
 *
 *   bun tutorial/part1-foundations/01-integers-and-divisibility/divisibility.ts
 *
 * ## Chapter Overview
 *
 * ### Chapter 1: Integers and Divisibility
 *
 * The fundamentals of number theory that underpin all cryptographic algorithms.
 *
 * - divisibility.ts    : Basic divisibility concepts, the divides relation
 * - gcd.ts             : Greatest Common Divisor, Euclidean algorithm
 * - extended-euclidean.ts : Extended Euclidean algorithm, Bezout's identity
 * - lcm.ts             : Least Common Multiple, connection to GCD
 *
 * Key concepts: divisibility, GCD, LCM, Bezout's identity, coprimality
 * Cryptographic relevance: RSA key generation, modular inverses
 *
 * ### Chapter 2: Modular Arithmetic
 *
 * Arithmetic in finite sets - the computational heart of cryptography.
 *
 * - congruences.ts       : Modular congruence, equivalence classes
 * - modular-operations.ts : Addition, subtraction, multiplication mod n
 * - modular-inverse.ts   : Computing modular inverses
 * - linear-congruences.ts : Solving ax = b (mod n), Chinese Remainder Theorem
 *
 * Key concepts: congruence, modular operations, inverses, CRT
 * Cryptographic relevance: RSA encryption/decryption, ECDSA signatures
 *
 * ### Chapter 3: Groups
 *
 * The abstract structure underlying discrete log cryptography.
 *
 * - group-axioms.ts    : What is a group, verifying group axioms
 * - cyclic-groups.ts   : Cyclic groups, generators
 * - order-of-elements.ts : Element orders, Lagrange's theorem
 * - subgroups.ts       : Subgroups, cosets, quotient groups
 *
 * Key concepts: groups, cyclic groups, order, subgroups, Lagrange's theorem
 * Cryptographic relevance: Diffie-Hellman, ECDSA, small-subgroup attacks
 *
 * ### Chapter 4: Rings and Fields
 *
 * Algebraic structures with two operations - essential for finite field crypto.
 *
 * - rings.ts           : Ring Z/nZ, zero divisors
 * - fields.ts          : Prime fields Z/pZ = GF(p)
 * - multiplicative-group.ts : The group (Z/nZ)*, structure theorems
 * - primitive-roots.ts : Primitive roots, generators of (Z/pZ)*
 *
 * Key concepts: rings, fields, zero divisors, units, primitive roots
 * Cryptographic relevance: RSA ring structure, ECC field arithmetic
 *
 * ## Prerequisites
 *
 * - Basic programming knowledge (TypeScript/JavaScript)
 * - High school algebra
 * - Curiosity about cryptography!
 *
 * ## Dependencies
 *
 * These tutorials use the sagemath-ts library:
 *
 *   import { ZZ, gcd, xgcd, Mod, ... } from 'sagemath-ts';
 *
 * ## Learning Path
 *
 * 1. Start with Chapter 1 to understand divisibility and the GCD
 * 2. Move to Chapter 2 to learn modular arithmetic
 * 3. Chapter 3 introduces the abstract group structure
 * 4. Chapter 4 ties everything together with rings and fields
 *
 * After completing Part 1, you'll be ready for:
 * - Part 2: Core Number-Theoretic Results
 * - Part 3: Computational Problems
 * - Part 4: Public-Key Cryptography (RSA, Diffie-Hellman)
 * - Part 5: Elliptic Curve Cryptography
 * - Part 6: Advanced Topics (Pairings, Lattices)
 *
 * ## Key Takeaways
 *
 * By the end of Part 1, you will understand:
 *
 * 1. Why gcd(a, n) = 1 is essential for modular inverses
 * 2. How the Euclidean algorithm efficiently computes GCDs
 * 3. Why RSA works (Euler's theorem)
 * 4. The structure of (Z/nZ)* and why it matters for security
 * 5. What makes a good cryptographic group (cyclic, large prime order)
 *
 * Happy learning!
 */

// This file serves as documentation. To run examples, execute the individual chapter files.

console.log("Part 1: Foundations of Modular Arithmetic");
console.log("=========================================\n");

console.log("Chapters:");
console.log("  1. Integers and Divisibility");
console.log("     - divisibility.ts");
console.log("     - gcd.ts");
console.log("     - extended-euclidean.ts");
console.log("     - lcm.ts\n");

console.log("  2. Modular Arithmetic");
console.log("     - congruences.ts");
console.log("     - modular-operations.ts");
console.log("     - modular-inverse.ts");
console.log("     - linear-congruences.ts\n");

console.log("  3. Groups");
console.log("     - group-axioms.ts");
console.log("     - cyclic-groups.ts");
console.log("     - order-of-elements.ts");
console.log("     - subgroups.ts\n");

console.log("  4. Rings and Fields");
console.log("     - rings.ts");
console.log("     - fields.ts");
console.log("     - multiplicative-group.ts");
console.log("     - primitive-roots.ts\n");

console.log("Run individual files with:");
console.log("  bun <path-to-file>");
