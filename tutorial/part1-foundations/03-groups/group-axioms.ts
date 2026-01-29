/**
 * # Groups: The Mathematical Structure Behind Cryptography
 *
 * Groups are the fundamental algebraic structure underlying almost all
 * modern cryptography. Understanding groups is essential for elliptic
 * curve cryptography, Diffie-Hellman, digital signatures, and more.
 *
 * ## Mathematical Background
 *
 * A group (G, *) is a set G with a binary operation * satisfying:
 *
 * 1. CLOSURE: For all a, b in G, a * b is in G
 * 2. ASSOCIATIVITY: For all a, b, c in G, (a * b) * c = a * (b * c)
 * 3. IDENTITY: There exists e in G such that e * a = a * e = a for all a
 * 4. INVERSE: For each a in G, there exists a^(-1) such that a * a^(-1) = e
 *
 * If also: 5. COMMUTATIVITY: a * b = b * a for all a, b
 * then G is called an ABELIAN (or commutative) group.
 *
 * ## Why This Matters for Cryptography
 *
 *   - Diffie-Hellman: Security based on the discrete log problem in groups
 *   - ECDSA: Digital signatures using elliptic curve groups
 *   - RSA: Uses the multiplicative group (Z/nZ)*
 *   - Pairings: Map between different groups for advanced protocols
 */

import { gcd, Mod, euler_phi, Zmod } from 'sagemath-ts';

// =============================================================================
// Example 1: Integers Under Addition - The Simplest Group
// =============================================================================
console.log("=== Example 1: (Z, +) - Integers Under Addition ===\n");

/**
 * The integers Z form a group under addition.
 * Let's verify all four group axioms.
 */
console.log("Verifying group axioms for (Z, +):\n");

// 1. Closure
console.log("1. CLOSURE: For any a, b in Z, a + b is in Z");
console.log("   3 + 5 = 8 (an integer)");
console.log("   -7 + 12 = 5 (an integer)");
console.log("   Integers are closed under addition. CHECK!\n");

// 2. Associativity
console.log("2. ASSOCIATIVITY: (a + b) + c = a + (b + c)");
const a = 3, b = 5, c = 7;
const left = (a + b) + c;
const right = a + (b + c);
console.log(`   (${a} + ${b}) + ${c} = ${left}`);
console.log(`   ${a} + (${b} + ${c}) = ${right}`);
console.log(`   Equal? ${left === right}. CHECK!\n`);

// 3. Identity
console.log("3. IDENTITY: There exists e such that e + a = a + e = a");
console.log("   The identity is 0:");
console.log(`   0 + 5 = 5`);
console.log(`   5 + 0 = 5`);
console.log("   CHECK!\n");

// 4. Inverse
console.log("4. INVERSE: For each a, there exists -a such that a + (-a) = 0");
console.log(`   The inverse of 5 is -5:`);
console.log(`   5 + (-5) = 0`);
console.log(`   The inverse of -3 is 3:`);
console.log(`   -3 + 3 = 0`);
console.log("   CHECK!\n");

// 5. Commutativity
console.log("5. COMMUTATIVITY (bonus): a + b = b + a");
console.log(`   3 + 5 = ${3 + 5}`);
console.log(`   5 + 3 = ${5 + 3}`);
console.log("   (Z, +) is an ABELIAN group. CHECK!");

console.log();

// =============================================================================
// Example 2: (Z/nZ, +) - Addition Modulo n
// =============================================================================
console.log("=== Example 2: (Z/nZ, +) - Addition Modulo n ===\n");

/**
 * Z/nZ under addition is a finite abelian group of order n.
 * This is the additive group of integers modulo n.
 */
const n = 6;
const Zn = Zmod(n);

console.log(`The group (Z/${n}Z, +) has ${n} elements: {0, 1, 2, 3, 4, 5}\n`);

// Build the addition (Cayley) table
console.log("Addition table:");
let header = "  + |";
for (let j = 0; j < n; j++) header += ` ${j}`;
console.log(header);
console.log("----+------------");

for (let i = 0; i < n; i++) {
  let row = `  ${i} |`;
  for (let j = 0; j < n; j++) {
    const sum = (i + j) % n;
    row += ` ${sum}`;
  }
  console.log(row);
}

console.log("\nKey observations:");
console.log("  - Each row/column contains each element exactly once (Latin square property)");
console.log("  - 0 is the identity element");
console.log(`  - Inverse of 1 is 5 (since 1 + 5 = 0 mod 6)`);
console.log(`  - Inverse of 2 is 4 (since 2 + 4 = 0 mod 6)`);

console.log();

// =============================================================================
// Example 3: (Z/nZ)* - Multiplicative Group
// =============================================================================
console.log("=== Example 3: (Z/nZ)* - Multiplicative Group ===\n");

/**
 * (Z/nZ)* is the set of elements in Z/nZ that are coprime to n.
 * Under multiplication, this forms a group.
 * Its order is phi(n) (Euler's totient function).
 */
const m = 10n;
const Zm = Zmod(m);

// Find all units (elements coprime to m)
const units: bigint[] = [];
for (let i = 1n; i < m; i++) {
  if (gcd(i, m) === 1n) {
    units.push(i);
  }
}

console.log(`(Z/${m}Z)* = {${units.join(', ')}}`);
console.log(`Order = phi(${m}) = ${euler_phi(m)}\n`);

// Verify group axioms
console.log("Verifying group axioms:\n");

// Closure
console.log("1. CLOSURE: Products of units are units");
const u1 = units[0]!, u2 = units[1]!;
const prod = (u1 * u2) % m;
console.log(`   ${u1} * ${u2} = ${prod} mod ${m}`);
console.log(`   Is ${prod} a unit? ${gcd(prod, m) === 1n}`);

// Identity
console.log("\n2. IDENTITY: 1 is the identity");
const testUnit = units[2]!;
console.log(`   1 * ${testUnit} = ${(1n * testUnit) % m}`);
console.log(`   ${testUnit} * 1 = ${(testUnit * 1n) % m}`);

// Inverse
console.log("\n3. INVERSES:");
for (const u of units) {
  const inv = Mod(u, m).inv();
  console.log(`   ${u}^(-1) = ${inv.value}  (${u} * ${inv.value} = ${(u * inv.value) % m} mod ${m})`);
}

console.log();

// =============================================================================
// Example 4: Non-Example - Why Z Under Division is NOT a Group
// =============================================================================
console.log("=== Example 4: Non-Example - (Z, /) is NOT a Group ===\n");

/**
 * Let's see why integers under division don't form a group.
 */
console.log("Attempting to form a group (Z, /):\n");

// Closure fails
console.log("1. CLOSURE: Fails!");
console.log("   5 / 3 = 1.666... (not an integer)");

// Associativity fails
console.log("\n2. ASSOCIATIVITY: Fails!");
console.log("   (12 / 6) / 2 = 2 / 2 = 1");
console.log("   12 / (6 / 2) = 12 / 3 = 4");
console.log("   1 != 4");

// Identity exists
console.log("\n3. IDENTITY: 1 works as right identity (a / 1 = a)");
console.log("   But 1 / a != a in general (fails as left identity)");

// Inverses don't exist
console.log("\n4. INVERSES: Don't exist for most elements");
console.log("   For a * a^(-1) = 1, we'd need a / (1/a) = a*a = 1");
console.log("   This only works for a = 1 or a = -1");

console.log("\n=> Division on integers does NOT form a group!");

console.log();

// =============================================================================
// Example 5: Symmetry Groups
// =============================================================================
console.log("=== Example 5: Symmetry Groups ===\n");

/**
 * Symmetries of geometric objects form groups.
 * The dihedral group D_3 is the symmetries of an equilateral triangle.
 */
console.log("D_3: Symmetry group of an equilateral triangle\n");

// Elements: identity, two rotations, three reflections
const D3_elements = ["e", "r", "r2", "s", "sr", "sr2"];

console.log("Elements:");
console.log("  e   = identity (do nothing)");
console.log("  r   = rotate 120 degrees counterclockwise");
console.log("  r2  = rotate 240 degrees counterclockwise");
console.log("  s   = reflect across vertical axis");
console.log("  sr  = reflect + rotate 120");
console.log("  sr2 = reflect + rotate 240\n");

// Group multiplication table
const D3_mult: { [key: string]: { [key: string]: string } } = {
  "e":   {"e": "e",   "r": "r",   "r2": "r2",  "s": "s",   "sr": "sr",  "sr2": "sr2"},
  "r":   {"e": "r",   "r": "r2",  "r2": "e",   "s": "sr2", "sr": "s",   "sr2": "sr"},
  "r2":  {"e": "r2",  "r": "e",   "r2": "r",   "s": "sr",  "sr": "sr2", "sr2": "s"},
  "s":   {"e": "s",   "r": "sr",  "r2": "sr2", "s": "e",   "sr": "r",   "sr2": "r2"},
  "sr":  {"e": "sr",  "r": "sr2", "r2": "s",   "s": "r2",  "sr": "e",   "sr2": "r"},
  "sr2": {"e": "sr2", "r": "s",   "r2": "sr",  "s": "r",   "sr": "r2",  "sr2": "e"}
};

console.log("Multiplication table (row * column):");
console.log("     |  e    r   r2   s   sr  sr2");
console.log("-----+-----------------------------");
for (const a of D3_elements) {
  let row = `  ${a.padEnd(3)}|`;
  for (const b of D3_elements) {
    row += ` ${D3_mult[a]![b]!.padEnd(3)}`;
  }
  console.log(row);
}

console.log("\nNote: D_3 is NOT abelian! r * s = sr2, but s * r = sr");

console.log();

// =============================================================================
// Example 6: Matrix Groups
// =============================================================================
console.log("=== Example 6: Matrix Groups ===\n");

/**
 * GL_n(F) is the group of invertible n x n matrices over field F.
 * Let's look at GL_2(Z/2Z) - 2x2 invertible matrices over the field of 2 elements.
 */
console.log("GL_2(Z/2Z): Invertible 2x2 matrices over {0, 1}\n");

// Represent 2x2 matrices as [a, b, c, d] meaning [[a,b],[c,d]]
type Mat2 = [number, number, number, number];

function matMult(A: Mat2, B: Mat2, mod: number): Mat2 {
  return [
    (A[0] * B[0] + A[1] * B[2]) % mod,
    (A[0] * B[1] + A[1] * B[3]) % mod,
    (A[2] * B[0] + A[3] * B[2]) % mod,
    (A[2] * B[1] + A[3] * B[3]) % mod
  ];
}

function det(A: Mat2, mod: number): number {
  return ((A[0] * A[3] - A[1] * A[2]) % mod + mod) % mod;
}

// Find all invertible 2x2 matrices over Z/2Z
const GL2_F2: Mat2[] = [];
for (let a = 0; a < 2; a++) {
  for (let b = 0; b < 2; b++) {
    for (let c = 0; c < 2; c++) {
      for (let d = 0; d < 2; d++) {
        const M: Mat2 = [a, b, c, d];
        if (det(M, 2) !== 0) {
          GL2_F2.push(M);
        }
      }
    }
  }
}

console.log(`|GL_2(Z/2Z)| = ${GL2_F2.length}\n`);

console.log("Elements (as [a,b; c,d]):");
for (const M of GL2_F2) {
  console.log(`  [${M[0]},${M[1]}; ${M[2]},${M[3]}]  det = ${det(M, 2)}`);
}

console.log(`\nThis is isomorphic to S_3 (the symmetric group on 3 elements)!`);

console.log();

// =============================================================================
// Example 7: The Group Axiom Checker
// =============================================================================
console.log("=== Example 7: Automated Group Axiom Checker ===\n");

/**
 * Let's build a tool to verify group axioms for finite sets.
 */
interface FiniteGroup<T> {
  elements: T[];
  operation: (a: T, b: T) => T;
  identity: T;
  inverse: (a: T) => T;
  equals: (a: T, b: T) => boolean;
}

function verifyGroupAxioms<T>(G: FiniteGroup<T>): boolean {
  const { elements, operation, identity, inverse, equals } = G;
  let allPassed = true;

  // 1. Closure
  console.log("Checking CLOSURE...");
  for (const a of elements) {
    for (const b of elements) {
      const result = operation(a, b);
      const inSet = elements.some(e => equals(e, result));
      if (!inSet) {
        console.log(`  FAILED: ${a} * ${b} = ${result} not in set`);
        allPassed = false;
      }
    }
  }
  if (allPassed) console.log("  PASSED");

  // 2. Associativity
  console.log("Checking ASSOCIATIVITY...");
  let assocPassed = true;
  for (const a of elements) {
    for (const b of elements) {
      for (const c of elements) {
        const left = operation(operation(a, b), c);
        const right = operation(a, operation(b, c));
        if (!equals(left, right)) {
          console.log(`  FAILED: (${a}*${b})*${c} != ${a}*(${b}*${c})`);
          assocPassed = false;
          allPassed = false;
        }
      }
    }
  }
  if (assocPassed) console.log("  PASSED");

  // 3. Identity
  console.log("Checking IDENTITY...");
  let identPassed = true;
  for (const a of elements) {
    if (!equals(operation(identity, a), a) || !equals(operation(a, identity), a)) {
      console.log(`  FAILED: ${identity} is not identity for ${a}`);
      identPassed = false;
      allPassed = false;
    }
  }
  if (identPassed) console.log(`  PASSED (identity = ${identity})`);

  // 4. Inverses
  console.log("Checking INVERSES...");
  let invPassed = true;
  for (const a of elements) {
    const inv = inverse(a);
    if (!equals(operation(a, inv), identity) || !equals(operation(inv, a), identity)) {
      console.log(`  FAILED: ${inv} is not inverse of ${a}`);
      invPassed = false;
      allPassed = false;
    }
  }
  if (invPassed) console.log("  PASSED");

  return allPassed;
}

// Test with (Z/5Z, +)
console.log("Testing (Z/5Z, +):\n");
const Z5_add: FiniteGroup<number> = {
  elements: [0, 1, 2, 3, 4],
  operation: (a, b) => (a + b) % 5,
  identity: 0,
  inverse: (a) => (5 - a) % 5,
  equals: (a, b) => a === b
};

const result = verifyGroupAxioms(Z5_add);
console.log(`\nAll axioms satisfied: ${result}`);

console.log();

// =============================================================================
// Exercise: Verify (Z/7Z)* is a Group
// =============================================================================
console.log("=== Exercise: Verify (Z/7Z)* is a Group ===\n");

/**
 * EXERCISE: Create a FiniteGroup object for (Z/7Z)* and verify the axioms.
 */
const Z7_mult: FiniteGroup<bigint> = {
  elements: [1n, 2n, 3n, 4n, 5n, 6n],
  operation: (a, b) => (a * b) % 7n,
  identity: 1n,
  inverse: (a) => Mod(a, 7n).inv().value,
  equals: (a, b) => a === b
};

console.log("Testing (Z/7Z)*:\n");
verifyGroupAxioms(Z7_mult);

console.log();

// =============================================================================
// Cryptographic Application: Groups in Cryptography
// =============================================================================
console.log("=== Cryptographic Application: Why Groups Matter ===\n");

console.log("Groups are the foundation of modern cryptography:\n");

console.log("1. DIFFIE-HELLMAN KEY EXCHANGE");
console.log("   - Uses the group (Z/p)* for large prime p");
console.log("   - Security: Discrete Log Problem is hard in this group");
console.log("   - Alice sends g^a, Bob sends g^b, shared secret is g^(ab)\n");

console.log("2. RSA CRYPTOSYSTEM");
console.log("   - Uses the group (Z/n)* where n = p*q");
console.log("   - Euler's theorem: a^phi(n) = 1 for a in (Z/n)*");
console.log("   - Decryption works because d*e = 1 mod phi(n)\n");

console.log("3. ELLIPTIC CURVE CRYPTOGRAPHY");
console.log("   - Points on an elliptic curve form a group");
console.log("   - Operation: Point addition (not multiplication!)");
console.log("   - Smaller keys than RSA for same security level\n");

console.log("4. DIGITAL SIGNATURES (ECDSA)");
console.log("   - Uses elliptic curve group");
console.log("   - Signature proves knowledge of private key");
console.log("   - Anyone can verify using public key\n");

console.log("Understanding group theory is essential for understanding cryptography!");
