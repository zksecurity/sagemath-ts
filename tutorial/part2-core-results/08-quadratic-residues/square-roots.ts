/**
 * # Computing Square Roots Modulo p
 *
 * Given a quadratic residue a modulo an odd prime p, how do we find
 * x such that x^2 = a (mod p)?
 *
 * ## Simple Case: p = 3 (mod 4)
 *
 * If p = 3 (mod 4) and a is a QR, then:
 *   x = a^((p+1)/4) (mod p)
 *
 * ## General Case: Tonelli-Shanks Algorithm
 *
 * Works for any odd prime p. Time complexity: O(log^2 p)
 *
 * @module tutorial/part2/square-roots
 */

import {
  sqrt_mod,
  power_mod,
  legendre_symbol,
  is_prime,
  Zmod,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("COMPUTING SQUARE ROOTS MODULO p");
console.log("=".repeat(70));

// ============================================================================
// Example 1: The Easy Case (p = 3 mod 4)
// ============================================================================
console.log("\n=== Example 1: Easy Case (p = 3 mod 4) ===\n");

console.log("When p = 3 (mod 4), computing square roots is simple:");
console.log("  If a is a QR, then sqrt(a) = a^((p+1)/4) mod p\n");

const easyCasePrimes = [7n, 11n, 19n, 23n, 31n];

console.log("Primes p = 3 (mod 4):");
for (const p of easyCasePrimes) {
  console.log(`  p = ${p.toString().padStart(2)}: ${p} mod 4 = ${p % 4n}`);
}

// Demonstrate
const p1 = 11n;
const a1 = 5n; // 5 is QR mod 11 (since 4^2 = 16 = 5 mod 11)

console.log(`\nExample: sqrt(${a1}) mod ${p1}`);
console.log(`  Check if ${a1} is QR: (${a1}/${p1}) = ${legendre_symbol(a1, p1)} = 1 (yes)`);

const exp1 = (p1 + 1n) / 4n;
const x1 = power_mod(a1, exp1, p1);
const verify1 = power_mod(x1, 2n, p1);

console.log(`  x = ${a1}^((${p1}+1)/4) = ${a1}^${exp1} mod ${p1} = ${x1}`);
console.log(`  Verify: ${x1}^2 mod ${p1} = ${verify1} [${verify1 === a1 ? "OK" : "FAIL"}]`);
console.log(`  Other root: ${p1} - ${x1} = ${p1 - x1}`);

// Why it works
console.log(`\nWhy does this work?`);
console.log(`  We want x^2 = a (mod p)`);
console.log(`  If x = a^((p+1)/4), then x^2 = a^((p+1)/2)`);
console.log(`  = a^((p-1)/2) * a = (a/p) * a = 1 * a = a (mod p)`);
console.log(`  The last step uses Euler's criterion: a^((p-1)/2) = (a/p) = 1`);

// ============================================================================
// Example 2: The Harder Case (p = 1 mod 4)
// ============================================================================
console.log("\n=== Example 2: Harder Case (p = 1 mod 4) ===\n");

console.log("When p = 1 (mod 4), we need the Tonelli-Shanks algorithm.\n");

const p2 = 13n; // 13 = 1 (mod 4)
console.log(`p = ${p2}: ${p2} mod 4 = ${p2 % 4n}`);

// List QRs and their square roots
console.log(`\nQuadratic residues mod ${p2} and their square roots:`);
console.log("-".repeat(40));

const Zp = Zmod(p2);
for (let a = 1n; a < p2; a++) {
  const leg = legendre_symbol(a, p2);
  if (leg === 1n) {
    const x = sqrt_mod(a, p2);
    const x2 = p2 - x!;
    console.log(`  sqrt(${a.toString().padStart(2)}) = {${x}, ${x2}} (verify: ${x}^2 = ${power_mod(x!, 2n, p2)})`);
  }
}

// ============================================================================
// Example 3: Tonelli-Shanks Algorithm
// ============================================================================
console.log("\n=== Example 3: Tonelli-Shanks Algorithm ===\n");

console.log("Tonelli-Shanks algorithm finds sqrt(a) mod p for any odd prime p.\n");

function tonelliShanksWithSteps(a: bigint, p: bigint): { x: bigint | null; steps: string[] } {
  const steps: string[] = [];

  // Check if a is QR
  const leg = legendre_symbol(a, p);
  if (leg === 0n) {
    steps.push(`a = 0 mod p, sqrt(0) = 0`);
    return { x: 0n, steps };
  }
  if (leg === -1n) {
    steps.push(`a is not a quadratic residue mod p`);
    return { x: null, steps };
  }

  steps.push(`Step 0: Verify a is QR: (${a}/${p}) = 1 [OK]`);

  // Easy case
  if (p % 4n === 3n) {
    const x = power_mod(a, (p + 1n) / 4n, p);
    steps.push(`Easy case (p = 3 mod 4): x = a^((p+1)/4) = ${x}`);
    return { x, steps };
  }

  // Write p - 1 = Q * 2^S with Q odd
  let Q = p - 1n;
  let S = 0n;
  while (Q % 2n === 0n) {
    Q /= 2n;
    S++;
  }
  steps.push(`Step 1: Write p-1 = Q * 2^S`);
  steps.push(`        ${p}-1 = ${Q} * 2^${S}`);

  // Find a quadratic non-residue z
  let z = 2n;
  while (legendre_symbol(z, p) !== -1n) {
    z++;
  }
  steps.push(`Step 2: Find quadratic non-residue z = ${z}`);

  let M = S;
  let c = power_mod(z, Q, p);
  let t = power_mod(a, Q, p);
  let R = power_mod(a, (Q + 1n) / 2n, p);

  steps.push(`Step 3: Initialize:`);
  steps.push(`        M = S = ${M}`);
  steps.push(`        c = z^Q = ${z}^${Q} mod ${p} = ${c}`);
  steps.push(`        t = a^Q = ${a}^${Q} mod ${p} = ${t}`);
  steps.push(`        R = a^((Q+1)/2) = ${a}^${(Q + 1n) / 2n} mod ${p} = ${R}`);

  let iteration = 0;
  while (t !== 1n) {
    iteration++;
    steps.push(`\nIteration ${iteration}:`);

    // Find the least i such that t^(2^i) = 1
    let i = 1n;
    let temp = (t * t) % p;
    while (temp !== 1n) {
      temp = (temp * temp) % p;
      i++;
    }
    steps.push(`  Find i: smallest i where t^(2^i) = 1, i = ${i}`);

    // Update
    const b = power_mod(c, 1n << (M - i - 1n), p);
    M = i;
    c = (b * b) % p;
    t = (t * c) % p;
    R = (R * b) % p;

    steps.push(`  b = c^(2^(M-i-1)) = ${b}`);
    steps.push(`  M = ${M}, c = ${c}, t = ${t}, R = ${R}`);
  }

  steps.push(`\nResult: x = ${R}`);
  steps.push(`Verify: ${R}^2 mod ${p} = ${power_mod(R, 2n, p)}`);

  return { x: R, steps };
}

// Example
const p3 = 41n;
const a3 = 5n;

console.log(`Computing sqrt(${a3}) mod ${p3}:`);
const { x: result3, steps: steps3 } = tonelliShanksWithSteps(a3, p3);
for (const step of steps3) {
  console.log(step);
}

console.log(`\nLibrary result: ${sqrt_mod(a3, p3)}`);

// ============================================================================
// Example 4: Special Case p = 5 (mod 8)
// ============================================================================
console.log("\n=== Example 4: Special Case (p = 5 mod 8) ===\n");

console.log("When p = 5 (mod 8), there's a simpler formula than full Tonelli-Shanks:");
console.log("  v = (2a)^((p-5)/8) mod p");
console.log("  i = 2av^2 mod p");
console.log("  x = av(i - 1) mod p\n");

const p4 = 13n; // 13 = 5 (mod 8)
const a4 = 10n; // 10 is QR mod 13

console.log(`p = ${p4}, ${p4} mod 8 = ${p4 % 8n}`);
console.log(`a = ${a4}, (${a4}/${p4}) = ${legendre_symbol(a4, p4)}`);

const v = power_mod(2n * a4, (p4 - 5n) / 8n, p4);
const i = (2n * a4 * v * v) % p4;
const x4 = (a4 * v * ((i - 1n + p4) % p4)) % p4;

console.log(`\nComputation:`);
console.log(`  v = (2*${a4})^((${p4}-5)/8) mod ${p4} = ${v}`);
console.log(`  i = 2*${a4}*v^2 mod ${p4} = ${i}`);
console.log(`  x = ${a4}*v*(i-1) mod ${p4} = ${x4}`);
console.log(`  Verify: ${x4}^2 mod ${p4} = ${power_mod(x4, 2n, p4)}`);
console.log(`\nLibrary: ${sqrt_mod(a4, p4)}`);

// ============================================================================
// Example 5: Getting All Square Roots
// ============================================================================
console.log("\n=== Example 5: Getting All Square Roots ===\n");

console.log("For any QR mod p, there are exactly 2 square roots: x and p-x");
console.log("(These are negatives of each other mod p)\n");

const p5 = 17n;
const testValues = [1n, 2n, 4n, 8n, 9n, 13n, 15n, 16n];

console.log(`Square roots mod ${p5}:`);
console.log("-".repeat(50));

for (const a of testValues) {
  const allRoots = sqrt_mod(a, p5, true); // Get all roots

  if (allRoots.length === 0) {
    console.log(`  sqrt(${a.toString().padStart(2)}) = no solution (QNR)`);
  } else if (allRoots.length === 1) {
    console.log(`  sqrt(${a.toString().padStart(2)}) = {${allRoots[0]}} (only 0 has one root)`);
  } else {
    console.log(`  sqrt(${a.toString().padStart(2)}) = {${allRoots.sort((x, y) => x < y ? -1 : 1).join(", ")}}`);
  }
}

// ============================================================================
// Example 6: Square Roots Mod p^k
// ============================================================================
console.log("\n=== Example 6: Hensel Lifting (Square Roots mod p^k) ===\n");

console.log("Hensel's Lemma: A square root mod p can be 'lifted' to mod p^k.\n");

function henselLiftSquareRoot(x: bigint, a: bigint, p: bigint, k: bigint): bigint {
  // Newton iteration: x_{n+1} = x_n - f(x_n) / f'(x_n)
  // For f(x) = x^2 - a, f'(x) = 2x
  // x_{n+1} = x - (x^2 - a) / (2x) = (x + a/x) / 2
  // In modular arithmetic: x_{new} = (x + a * x^(-1)) * 2^(-1) mod p^(k+1)

  let currentMod = p;
  let currentX = x % p;

  for (let i = 1n; i < k; i++) {
    const newMod = currentMod * p;
    // x_new = x - (x^2 - a) * (2x)^(-1) mod newMod
    const twoXInv = power_mod(2n * currentX, newMod - newMod / p - 1n, newMod); // Approximation
    // Actually use: x_new = (x + a * inv(x)) / 2
    const xInv = power_mod(currentX, currentMod - 2n, currentMod);
    const sum = currentX + ((a * xInv) % newMod);
    currentX = (sum * power_mod(2n, newMod - newMod / p - 1n, newMod)) % newMod;

    currentMod = newMod;
  }

  return currentX;
}

// Simpler approach: direct computation
const p6 = 7n;
const a6 = 2n;

console.log(`Finding square roots of ${a6} mod ${p6}^k:`);

// First, find root mod p
const x_mod_p = sqrt_mod(a6, p6);
if (x_mod_p !== null) {
  console.log(`  mod ${p6}^1 = ${p6}: sqrt(${a6}) = ${x_mod_p}, verify: ${x_mod_p}^2 mod ${p6} = ${power_mod(x_mod_p, 2n, p6)}`);

  // Lift to higher powers by checking
  for (let k = 2n; k <= 4n; k++) {
    const pk = p6 ** k;
    // Find x such that x^2 = a mod p^k and x = x_mod_p mod p
    let found = false;
    for (let j = 0n; j < p6 ** (k - 1n); j++) {
      const candidate = x_mod_p + j * p6;
      if (power_mod(candidate, 2n, pk) === a6 % pk) {
        console.log(`  mod ${p6}^${k} = ${pk}: sqrt(${a6}) = ${candidate}, verify: ${candidate}^2 mod ${pk} = ${power_mod(candidate, 2n, pk)}`);
        found = true;
        break;
      }
    }
    if (!found) console.log(`  mod ${p6}^${k}: no solution`);
  }
} else {
  console.log(`  ${a6} is not a QR mod ${p6}`);
}

// ============================================================================
// Example 7: Complexity Analysis
// ============================================================================
console.log("\n=== Example 7: Complexity Analysis ===\n");

console.log("Algorithm Complexity:");
console.log("-".repeat(50));
console.log("");
console.log("p = 3 (mod 4) case:");
console.log("  - Single exponentiation: O(log p) multiplications");
console.log("  - Each multiplication: O(log^2 p) or O(log p log log p) with FFT");
console.log("  - Total: O(log^3 p) or O(log^2 p log log p)");
console.log("");
console.log("Tonelli-Shanks (general case):");
console.log("  - Finding quadratic non-residue: O(log p) expected");
console.log("  - Main loop: O(log p) iterations");
console.log("  - Each iteration: O(log p) multiplications");
console.log("  - Total: O(log^2 p) multiplications = O(log^4 p) bit operations");
console.log("");
console.log("Cipolla's algorithm (alternative):");
console.log("  - Extends field to Fp[x]/(x^2 - n)");
console.log("  - Single exponentiation in extension field");
console.log("  - Same asymptotic complexity, sometimes faster in practice");

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
Computing Square Roots mod p:

First, check if a is a QR: (a/p) = a^((p-1)/2) mod p = 1

Easy cases:
  - p = 3 (mod 4): x = a^((p+1)/4) mod p
  - p = 5 (mod 8): x = a * v * (i-1) where v = (2a)^((p-5)/8), i = 2av^2

General case: Tonelli-Shanks algorithm
  1. Write p-1 = Q * 2^S with Q odd
  2. Find quadratic non-residue z
  3. Initialize M, c, t, R
  4. Iterate until t = 1

Properties:
  - Each QR has exactly 2 square roots: x and p-x
  - sqrt_mod(a, p) returns one root (or null if QNR)
  - sqrt_mod(a, p, true) returns all roots

Using sagemath-ts:
  sqrt_mod(a, p)         - One square root or null
  sqrt_mod(a, p, true)   - Array of all square roots
`);
