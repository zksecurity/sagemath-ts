/**
 * # The Chinese Remainder Theorem - Statement and Examples
 *
 * The Chinese Remainder Theorem (CRT) is one of the most important
 * results in number theory, with applications throughout cryptography.
 *
 * ## Historical Note
 *
 * The theorem dates back to the 3rd century CE in China, where
 * mathematician Sun Tzu posed a problem: "Find a number which gives
 * remainder 2 when divided by 3, remainder 3 when divided by 5,
 * and remainder 2 when divided by 7."
 *
 * ## Statement
 *
 * Let m1, m2, ..., mk be pairwise coprime positive integers.
 * For any integers a1, a2, ..., ak, the system of congruences:
 *
 *   x = a1 (mod m1)
 *   x = a2 (mod m2)
 *   ...
 *   x = ak (mod mk)
 *
 * has a unique solution modulo M = m1 * m2 * ... * mk.
 *
 * @module tutorial/part2/crt-statement
 */

import {
  crt,
  CRT_list,
  gcd,
  Zmod,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("THE CHINESE REMAINDER THEOREM");
console.log("=".repeat(70));

// ============================================================================
// Example 1: Sun Tzu's Original Problem
// ============================================================================
console.log("\n=== Example 1: Sun Tzu's Original Problem ===\n");

console.log("Find x such that:");
console.log("  x = 2 (mod 3)");
console.log("  x = 3 (mod 5)");
console.log("  x = 2 (mod 7)\n");

const residues1 = [2n, 3n, 2n];
const moduli1 = [3n, 5n, 7n];

// Verify moduli are pairwise coprime
console.log("Checking pairwise coprimality:");
for (let i = 0; i < moduli1.length; i++) {
  for (let j = i + 1; j < moduli1.length; j++) {
    const g = gcd(moduli1[i]!, moduli1[j]!);
    console.log(`  gcd(${moduli1[i]}, ${moduli1[j]}) = ${g}`);
  }
}

const x1 = CRT_list(residues1, moduli1);
const M1 = moduli1.reduce((a, b) => a * b);

console.log(`\nSolution: x = ${x1}`);
console.log(`Unique modulo M = ${moduli1.join(" * ")} = ${M1}`);

console.log("\nVerification:");
for (let i = 0; i < moduli1.length; i++) {
  const remainder = x1 % moduli1[i]!;
  console.log(`  ${x1} mod ${moduli1[i]} = ${remainder} [expected: ${residues1[i]}] ${remainder === residues1[i] ? "OK" : "FAIL"}`);
}

// Show all solutions in a range
console.log(`\nAll solutions less than ${2n * M1}:`);
for (let k = 0n; k < 3n; k++) {
  console.log(`  x = ${x1} + ${k} * ${M1} = ${x1 + k * M1}`);
}

// ============================================================================
// Example 2: Two Congruences (Basic CRT)
// ============================================================================
console.log("\n=== Example 2: Two Congruences ===\n");

console.log("Find x such that:");
console.log("  x = 1 (mod 4)");
console.log("  x = 2 (mod 5)\n");

const x2 = crt(1n, 2n, 4n, 5n);
console.log(`Solution: x = ${x2}`);
console.log(`Unique modulo ${4n * 5n} = 20`);

console.log("\nVerification:");
console.log(`  ${x2} mod 4 = ${x2 % 4n} [expected: 1]`);
console.log(`  ${x2} mod 5 = ${x2 % 5n} [expected: 2]`);

// Manual verification by listing
console.log("\nManual search (checking all x from 0 to 19):");
for (let x = 0n; x < 20n; x++) {
  if (x % 4n === 1n && x % 5n === 2n) {
    console.log(`  Found: x = ${x}`);
  }
}

// ============================================================================
// Example 3: Understanding Uniqueness
// ============================================================================
console.log("\n=== Example 3: Understanding Uniqueness ===\n");

console.log("CRT gives a UNIQUE solution modulo M = m1 * m2 * ... * mk\n");

// Different systems with same moduli but different residues
const systems = [
  { residues: [0n, 0n], moduli: [3n, 5n], description: "x = 0 (mod 3), x = 0 (mod 5)" },
  { residues: [1n, 0n], moduli: [3n, 5n], description: "x = 1 (mod 3), x = 0 (mod 5)" },
  { residues: [2n, 0n], moduli: [3n, 5n], description: "x = 2 (mod 3), x = 0 (mod 5)" },
  { residues: [0n, 1n], moduli: [3n, 5n], description: "x = 0 (mod 3), x = 1 (mod 5)" },
  { residues: [1n, 1n], moduli: [3n, 5n], description: "x = 1 (mod 3), x = 1 (mod 5)" },
];

console.log("Solutions for different systems mod 3 and mod 5:");
console.log("-".repeat(60));

for (const { residues, moduli, description } of systems) {
  const x = CRT_list(residues, moduli);
  console.log(`${description}`);
  console.log(`  Solution: x = ${x} (mod 15)\n`);
}

console.log("Each of the 3*5 = 15 residue classes mod 15 corresponds");
console.log("to exactly one pair of residue classes mod 3 and mod 5.");

// ============================================================================
// Example 4: The Ring Isomorphism
// ============================================================================
console.log("\n=== Example 4: CRT as Ring Isomorphism ===\n");

console.log("CRT gives an isomorphism: Z/mnZ -> Z/mZ x Z/nZ (when gcd(m,n) = 1)\n");

const m4 = 3n;
const n4 = 5n;
const mn = m4 * n4;

console.log(`Mapping elements of Z/${mn}Z to Z/${m4}Z x Z/${n4}Z:`);
console.log("-".repeat(50));
console.log("x (mod 15) | (x mod 3, x mod 5)");
console.log("-----------+-------------------");

for (let x = 0n; x < mn; x++) {
  console.log(`    ${x.toString().padStart(2)}     |    (${x % m4}, ${x % n4})`);
}

console.log("\nThis mapping is:");
console.log("  - One-to-one (bijective)");
console.log("  - Preserves addition: f(x+y) = f(x) + f(y)");
console.log("  - Preserves multiplication: f(x*y) = f(x) * f(y)");

// Demonstrate addition preservation
const a = 7n;
const b = 11n;
const sum = (a + b) % mn;
console.log(`\nAddition example: f(${a} + ${b}) = f(${sum})`);
console.log(`  f(${a}) = (${a % m4}, ${a % n4})`);
console.log(`  f(${b}) = (${b % m4}, ${b % n4})`);
console.log(`  f(${a}) + f(${b}) = (${(a % m4 + b % m4) % m4}, ${(a % n4 + b % n4) % n4})`);
console.log(`  f(${sum}) = (${sum % m4}, ${sum % n4})`);

// ============================================================================
// Example 5: When CRT Doesn't Apply
// ============================================================================
console.log("\n=== Example 5: When CRT Doesn't Apply ===\n");

console.log("CRT requires PAIRWISE COPRIME moduli.\n");

// Non-coprime moduli
const badM1 = 6n;
const badM2 = 4n;
console.log(`Consider m1 = ${badM1}, m2 = ${badM2}`);
console.log(`gcd(${badM1}, ${badM2}) = ${gcd(badM1, badM2)} != 1 (not coprime)\n`);

console.log("System: x = 1 (mod 6), x = 3 (mod 4)");
console.log("\nChecking for solutions by exhaustive search (mod 12):");
let found = false;
for (let x = 0n; x < 24n; x++) {
  const cond1 = x % 6n === 1n;
  const cond2 = x % 4n === 3n;
  if (cond1 && cond2) {
    console.log(`  x = ${x} satisfies both conditions`);
    found = true;
  }
}
if (!found) {
  console.log("  NO SOLUTION EXISTS!");
}

console.log("\nWhy? If x = 1 (mod 6), then x = 1 (mod 2) (odd number)");
console.log("If x = 3 (mod 4), then x = 3 = 1 (mod 2) (odd number)");
console.log("These are consistent for mod 2, so there's no immediate contradiction.");

// But try another system
console.log("\n\nSystem: x = 1 (mod 6), x = 2 (mod 4)");
console.log("Checking:");
found = false;
for (let x = 0n; x < 24n; x++) {
  const cond1 = x % 6n === 1n;
  const cond2 = x % 4n === 2n;
  if (cond1 && cond2) {
    console.log(`  x = ${x} satisfies both conditions`);
    found = true;
  }
}
if (!found) {
  console.log("  NO SOLUTION EXISTS!");
}

console.log("\nWhy? x = 1 (mod 6) means x is odd.");
console.log("But x = 2 (mod 4) means x is even. Contradiction!");

// ============================================================================
// Example 6: Three or More Congruences
// ============================================================================
console.log("\n=== Example 6: Multiple Congruences ===\n");

const system6 = [
  { a: 1n, m: 2n },
  { a: 2n, m: 3n },
  { a: 3n, m: 5n },
  { a: 4n, m: 7n },
];

console.log("System of congruences:");
for (const { a, m } of system6) {
  console.log(`  x = ${a} (mod ${m})`);
}

const residues6 = system6.map(s => s.a);
const moduli6 = system6.map(s => s.m);
const M6 = moduli6.reduce((a, b) => a * b);

const x6 = CRT_list(residues6, moduli6);

console.log(`\nProduct of moduli: M = ${moduli6.join(" * ")} = ${M6}`);
console.log(`Solution: x = ${x6}`);

console.log("\nVerification:");
for (const { a, m } of system6) {
  const r = x6 % m;
  console.log(`  ${x6} mod ${m} = ${r} [expected: ${a}] ${r === a ? "OK" : "FAIL"}`);
}

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
Chinese Remainder Theorem (CRT):

Statement:
  For pairwise coprime moduli m1, m2, ..., mk, the system
    x = a1 (mod m1)
    x = a2 (mod m2)
    ...
    x = ak (mod mk)
  has a unique solution modulo M = m1 * m2 * ... * mk.

Key Points:
1. Moduli must be PAIRWISE COPRIME
2. Solution is unique modulo the product of moduli
3. Provides a ring isomorphism: Z/MZ -> Z/m1Z x Z/m2Z x ... x Z/mkZ
4. If moduli are not coprime, solution may not exist or may not be unique

Using sagemath-ts:
  crt(a, b, m, n)         - Solve x = a (mod m), x = b (mod n)
  CRT_list(residues, moduli) - Solve system with multiple congruences
`);
