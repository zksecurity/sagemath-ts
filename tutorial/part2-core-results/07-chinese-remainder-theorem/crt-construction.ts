/**
 * # Constructive Proof of the Chinese Remainder Theorem
 *
 * The CRT is not just an existence theorem - it provides an explicit
 * algorithm for constructing the solution.
 *
 * ## The Construction
 *
 * Given the system:
 *   x = a1 (mod m1)
 *   x = a2 (mod m2)
 *
 * where gcd(m1, m2) = 1:
 *
 * 1. Use extended Euclidean algorithm to find s, t with:
 *    s*m1 + t*m2 = 1
 *
 * 2. The solution is:
 *    x = a1*t*m2 + a2*s*m1 (mod m1*m2)
 *
 * Why this works:
 *   - x mod m1 = a1*t*m2 mod m1 = a1*(1 - s*m1) mod m1 = a1
 *   - x mod m2 = a2*s*m1 mod m2 = a2*(1 - t*m2) mod m2 = a2
 *
 * @module tutorial/part2/crt-construction
 */

import {
  crt,
  CRT_list,
  xgcd,
  gcd,
  inverse_mod,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("CONSTRUCTIVE PROOF OF CRT");
console.log("=".repeat(70));

// ============================================================================
// The Algorithm: Two Congruences
// ============================================================================
console.log("\n=== Algorithm for Two Congruences ===\n");

function crtManual(a1: bigint, a2: bigint, m1: bigint, m2: bigint): {
  solution: bigint;
  steps: string[];
} {
  const steps: string[] = [];

  steps.push(`System: x = ${a1} (mod ${m1}), x = ${a2} (mod ${m2})`);

  // Step 1: Extended GCD
  const [g, s, t] = xgcd(m1, m2);
  steps.push(`\nStep 1: Extended Euclidean Algorithm`);
  steps.push(`  ${s} * ${m1} + ${t} * ${m2} = ${g}`);

  if (g !== 1n) {
    steps.push(`  gcd(${m1}, ${m2}) = ${g} != 1, no solution!`);
    return { solution: -1n, steps };
  }

  // Step 2: Construct solution
  steps.push(`\nStep 2: Construct Solution`);
  const term1 = a1 * t * m2;
  const term2 = a2 * s * m1;
  steps.push(`  x = a1 * t * m2 + a2 * s * m1`);
  steps.push(`    = ${a1} * ${t} * ${m2} + ${a2} * ${s} * ${m1}`);
  steps.push(`    = ${term1} + ${term2}`);
  steps.push(`    = ${term1 + term2}`);

  // Step 3: Reduce modulo M
  const M = m1 * m2;
  let x = ((term1 + term2) % M + M) % M;
  steps.push(`\nStep 3: Reduce modulo M = ${m1} * ${m2} = ${M}`);
  steps.push(`  x = ${x}`);

  // Verify
  steps.push(`\nVerification:`);
  steps.push(`  ${x} mod ${m1} = ${x % m1} [expected: ${a1}]`);
  steps.push(`  ${x} mod ${m2} = ${x % m2} [expected: ${a2}]`);

  return { solution: x, steps };
}

// Example
const { solution, steps } = crtManual(2n, 3n, 5n, 7n);
for (const step of steps) {
  console.log(step);
}

console.log(`\nLibrary result: ${crt(2n, 3n, 5n, 7n)}`);

// ============================================================================
// Understanding the Construction
// ============================================================================
console.log("\n=== Understanding the Construction ===\n");

console.log("Key insight: Construct 'selector' terms that are:");
console.log("  - 1 modulo one modulus");
console.log("  - 0 modulo the other modulus\n");

const m1 = 5n;
const m2 = 7n;
const [_g, s, t] = xgcd(m1, m2);

// e1 = t * m2 satisfies: e1 = 1 (mod m1), e1 = 0 (mod m2)
// e2 = s * m1 satisfies: e2 = 0 (mod m1), e2 = 1 (mod m2)
const e1 = ((t * m2) % (m1 * m2) + m1 * m2) % (m1 * m2);
const e2 = ((s * m1) % (m1 * m2) + m1 * m2) % (m1 * m2);

console.log(`With m1 = ${m1}, m2 = ${m2}:`);
console.log(`  s = ${s}, t = ${t} (from xgcd)`);
console.log(`\nSelector e1 = t * m2 = ${t} * ${m2} = ${e1} (mod 35):`);
console.log(`  e1 mod m1 = ${e1 % m1}  [selector for a1]`);
console.log(`  e1 mod m2 = ${e1 % m2}  [zero contribution to a2]`);

console.log(`\nSelector e2 = s * m1 = ${s} * ${m1} = ${e2} (mod 35):`);
console.log(`  e2 mod m1 = ${e2 % m1}  [zero contribution to a1]`);
console.log(`  e2 mod m2 = ${e2 % m2}  [selector for a2]`);

console.log(`\nSolution formula: x = a1 * e1 + a2 * e2 (mod ${m1 * m2})`);
console.log(`This gives: x mod m1 = a1*1 + a2*0 = a1`);
console.log(`            x mod m2 = a1*0 + a2*1 = a2`);

// ============================================================================
// General Algorithm: Multiple Congruences
// ============================================================================
console.log("\n=== General Algorithm: Multiple Congruences ===\n");

function crtGeneralManual(
  residues: bigint[],
  moduli: bigint[]
): { solution: bigint; selectors: bigint[]; steps: string[] } {
  const k = residues.length;
  const steps: string[] = [];
  const M = moduli.reduce((a, b) => a * b);

  steps.push(`System with ${k} congruences:`);
  for (let i = 0; i < k; i++) {
    steps.push(`  x = ${residues[i]} (mod ${moduli[i]})`);
  }
  steps.push(`\nProduct of moduli: M = ${M}`);

  // Compute selectors (also called CRT basis)
  const selectors: bigint[] = [];

  steps.push(`\nComputing CRT basis (selectors):`);
  for (let i = 0; i < k; i++) {
    const mi = moduli[i]!;
    const Mi = M / mi; // Product of all moduli except mi
    const MiInv = inverse_mod(Mi, mi); // Mi^(-1) mod mi
    const ei = (Mi * MiInv) % M;
    selectors.push(ei);

    steps.push(`\n  For i = ${i}, m_${i} = ${mi}:`);
    steps.push(`    M_${i} = M / m_${i} = ${Mi}`);
    steps.push(`    M_${i}^(-1) mod m_${i} = ${MiInv}`);
    steps.push(`    e_${i} = M_${i} * M_${i}^(-1) = ${ei}`);
    steps.push(`    Check: e_${i} mod m_${i} = ${ei % mi} (should be 1)`);
  }

  // Compute solution
  let x = 0n;
  steps.push(`\nComputing solution: x = sum of a_i * e_i`);
  for (let i = 0; i < k; i++) {
    const term = residues[i]! * selectors[i]!;
    x = (x + term) % M;
    steps.push(`  + ${residues[i]} * ${selectors[i]} = ${term}`);
  }
  steps.push(`  x = ${x} (mod ${M})`);

  // Verify
  steps.push(`\nVerification:`);
  for (let i = 0; i < k; i++) {
    steps.push(`  x mod ${moduli[i]} = ${x % moduli[i]!} [expected: ${residues[i]}]`);
  }

  return { solution: x, selectors, steps };
}

const residuesGen = [2n, 3n, 2n];
const moduliGen = [3n, 5n, 7n];

const resultGen = crtGeneralManual(residuesGen, moduliGen);
for (const step of resultGen.steps) {
  console.log(step);
}

console.log(`\nCRT basis (selectors): [${resultGen.selectors.join(", ")}]`);
console.log(`Library result: ${CRT_list(residuesGen, moduliGen)}`);

// ============================================================================
// The CRT Basis as Idempotents
// ============================================================================
console.log("\n=== CRT Basis as Orthogonal Idempotents ===\n");

console.log("The selectors e_i form a system of orthogonal idempotents:");
console.log("  - e_i^2 = e_i (idempotent)");
console.log("  - e_i * e_j = 0 for i != j (orthogonal)");
console.log("  - sum of e_i = 1\n");

const moduli4 = [3n, 5n, 7n];
const M4 = moduli4.reduce((a, b) => a * b);

// Compute selectors
const selectors4: bigint[] = [];
for (const mi of moduli4) {
  const Mi = M4 / mi;
  const MiInv = inverse_mod(Mi, mi);
  selectors4.push((Mi * MiInv) % M4);
}

console.log(`Moduli: [${moduli4.join(", ")}], M = ${M4}`);
console.log(`Selectors: [${selectors4.join(", ")}]`);

// Verify idempotent property
console.log("\nIdempotent property (e_i^2 = e_i):");
for (let i = 0; i < selectors4.length; i++) {
  const e = selectors4[i]!;
  const e2 = (e * e) % M4;
  console.log(`  e_${i}^2 = ${e}^2 mod ${M4} = ${e2} = e_${i} [${e2 === e ? "OK" : "FAIL"}]`);
}

// Verify orthogonality
console.log("\nOrthogonality (e_i * e_j = 0 for i != j):");
for (let i = 0; i < selectors4.length; i++) {
  for (let j = i + 1; j < selectors4.length; j++) {
    const prod = (selectors4[i]! * selectors4[j]!) % M4;
    console.log(`  e_${i} * e_${j} = ${selectors4[i]} * ${selectors4[j]} mod ${M4} = ${prod} [${prod === 0n ? "OK" : "FAIL"}]`);
  }
}

// Verify sum is 1
const sum = selectors4.reduce((a, b) => (a + b) % M4, 0n);
console.log(`\nSum: e_0 + e_1 + e_2 = ${sum} [${sum === 1n ? "OK" : "FAIL"}]`);

// ============================================================================
// Iterative CRT Algorithm
// ============================================================================
console.log("\n=== Iterative CRT Algorithm ===\n");

console.log("Alternative approach: solve two congruences at a time.\n");

function crtIterative(residues: bigint[], moduli: bigint[]): {
  solution: bigint;
  steps: string[];
} {
  const steps: string[] = [];

  let currentX = residues[0]!;
  let currentM = moduli[0]!;

  steps.push(`Start: x = ${currentX} (mod ${currentM})`);

  for (let i = 1; i < residues.length; i++) {
    const ai = residues[i]!;
    const mi = moduli[i]!;

    steps.push(`\nCombine with: x = ${ai} (mod ${mi})`);

    // Solve: x = currentX (mod currentM) and x = ai (mod mi)
    const [g, s, t] = xgcd(currentM, mi);
    const newX = currentX * t * mi + ai * s * currentM;
    const newM = currentM * mi;
    currentX = ((newX % newM) + newM) % newM;
    currentM = newM;

    steps.push(`  New: x = ${currentX} (mod ${currentM})`);
  }

  return { solution: currentX, steps };
}

const residuesIter = [1n, 2n, 3n, 4n];
const moduliIter = [2n, 3n, 5n, 7n];

const resultIter = crtIterative(residuesIter, moduliIter);
for (const step of resultIter.steps) {
  console.log(step);
}

console.log(`\nFinal solution: ${resultIter.solution}`);
console.log(`Verification: ${CRT_list(residuesIter, moduliIter)}`);

// ============================================================================
// Complexity Analysis
// ============================================================================
console.log("\n=== Complexity Analysis ===\n");

console.log("Time Complexity for k congruences with moduli of size O(m):\n");

console.log("Direct method (compute all selectors):");
console.log("  - Compute M = product of moduli: O(k) multiplications");
console.log("  - For each i: compute M_i, M_i^(-1): O(log m) using xgcd");
console.log("  - Total: O(k log m) operations on O(km)-bit numbers\n");

console.log("Iterative method:");
console.log("  - k-1 pairwise CRT computations");
console.log("  - Each step: one xgcd, O(1) multiplications");
console.log("  - Numbers grow: O(m), O(m^2), ..., O(m^k)");
console.log("  - Total: Similar complexity, but may have numerical advantages\n");

console.log("Both methods are efficient and used in practice.");
console.log("The iterative method is often preferred for its simplicity.");

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
Constructive CRT Algorithm:

Two Congruences:
  Given x = a1 (mod m1), x = a2 (mod m2) with gcd(m1,m2) = 1:
  1. Find s,t with s*m1 + t*m2 = 1 (extended gcd)
  2. x = a1*t*m2 + a2*s*m1 (mod m1*m2)

Multiple Congruences (Selector Method):
  1. Compute M = product of all moduli
  2. For each i: compute selector e_i = (M/m_i) * (M/m_i)^(-1) mod M
  3. x = sum of a_i * e_i (mod M)

Multiple Congruences (Iterative):
  1. Start with first congruence
  2. Combine with next congruence using two-congruence CRT
  3. Repeat until all congruences are processed

Key Properties of Selectors:
  - e_i = 1 (mod m_i), e_i = 0 (mod m_j) for j != i
  - e_i^2 = e_i (idempotent)
  - e_i * e_j = 0 for i != j (orthogonal)
  - sum of e_i = 1
`);
