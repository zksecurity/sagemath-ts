/**
 * # Applications of the Chinese Remainder Theorem
 *
 * CRT has wide applications in:
 * 1. Solving systems of linear congruences
 * 2. Calendar calculations (determining days)
 * 3. Secret sharing
 * 4. Parallel computation
 * 5. Error-correcting codes
 *
 * @module tutorial/part2/crt-applications
 */

import {
  crt,
  CRT_list,
  gcd,
  xgcd,
  euler_phi,
  power_mod,
  inverse_mod,
  Zmod,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("APPLICATIONS OF THE CHINESE REMAINDER THEOREM");
console.log("=".repeat(70));

// ============================================================================
// Application 1: Solving Systems of Congruences
// ============================================================================
console.log("\n=== Application 1: Solving Systems of Congruences ===\n");

console.log("Problem: Find all integers that leave:");
console.log("  - Remainder 2 when divided by 3");
console.log("  - Remainder 3 when divided by 5");
console.log("  - Remainder 4 when divided by 7\n");

const residues1 = [2n, 3n, 4n];
const moduli1 = [3n, 5n, 7n];
const x1 = CRT_list(residues1, moduli1);
const M1 = moduli1.reduce((a, b) => a * b);

console.log(`Solution: x = ${x1} (mod ${M1})`);
console.log(`All solutions: x = ${x1}, ${x1 + M1}, ${x1 + 2n * M1}, ...`);
console.log(`General form: x = ${x1} + ${M1}k for any integer k\n`);

console.log("Verification:");
for (let i = 0; i < moduli1.length; i++) {
  console.log(`  ${x1} mod ${moduli1[i]} = ${x1 % moduli1[i]!} [expected: ${residues1[i]}]`);
}

// ============================================================================
// Application 2: Calendar Calculations
// ============================================================================
console.log("\n=== Application 2: Calendar Calculations ===\n");

console.log("Problem: What day of the week is day X of the year?");
console.log("Days repeat every 7 days, weeks repeat every 7 days...\n");

// More interesting: Find a day that satisfies multiple conditions
console.log("Find a year number Y such that:");
console.log("  - Y mod 4 = 0 (divisible by 4 - leap year candidate)");
console.log("  - Y mod 100 = 1 (1 more than a century)");
console.log("  - Y mod 7 = 0 (year number divisible by 7)\n");

// Note: 4 and 100 are not coprime, so we need to handle this
// Let's use coprime moduli instead
console.log("Simplified problem with coprime moduli:");
console.log("  - Y mod 4 = 0");
console.log("  - Y mod 25 = 1 (using 25 instead of 100)");
console.log("  - Y mod 7 = 0\n");

const calResidues = [0n, 1n, 0n];
const calModuli = [4n, 25n, 7n];
const yearSolution = CRT_list(calResidues, calModuli);
const calM = calModuli.reduce((a, b) => a * b);

console.log(`Solution: Y = ${yearSolution} (mod ${calM})`);
console.log(`First few solutions: ${yearSolution}, ${yearSolution + calM}, ${yearSolution + 2n * calM}`);

// ============================================================================
// Application 3: Secret Sharing (Shamir-like)
// ============================================================================
console.log("\n=== Application 3: CRT-Based Secret Sharing ===\n");

console.log("Share a secret among n parties such that any k parties can recover it.");
console.log("(This is a simplified version; Shamir's scheme uses polynomials)\n");

// Secret to share
const secret = 42n;

// Choose pairwise coprime moduli (must be > secret)
const shareModuli = [101n, 103n, 107n, 109n]; // All primes > 42
const threshold = 3; // Need 3 shares to recover

console.log(`Secret: ${secret}`);
console.log(`Moduli: [${shareModuli.join(", ")}]`);
console.log(`Threshold: ${threshold} shares needed\n`);

// Create shares: share_i = secret mod m_i
const shares = shareModuli.map(m => secret % m);
console.log("Shares:");
for (let i = 0; i < shares.length; i++) {
  console.log(`  Party ${i + 1}: (modulus=${shareModuli[i]}, share=${shares[i]})`);
}

// Recover with 3 shares
console.log("\nRecovery with parties 1, 2, 3:");
const recoveryResidues = shares.slice(0, 3);
const recoveryModuli = shareModuli.slice(0, 3);
const recovered = CRT_list(recoveryResidues, recoveryModuli);
const recoveryM = recoveryModuli.reduce((a, b) => a * b);

console.log(`  CRT solution: ${recovered} (mod ${recoveryM})`);
console.log(`  Recovered secret: ${recovered % (recoveryM)} = ${recovered}`);
console.log(`  Correct: ${recovered === secret}`);

console.log("\nRecovery with parties 2, 3, 4:");
const recoveryResidues2 = shares.slice(1, 4);
const recoveryModuli2 = shareModuli.slice(1, 4);
const recovered2 = CRT_list(recoveryResidues2, recoveryModuli2);
console.log(`  Recovered secret: ${recovered2}`);
console.log(`  Correct: ${recovered2 === secret}`);

// ============================================================================
// Application 4: Parallel Modular Arithmetic
// ============================================================================
console.log("\n=== Application 4: Parallel Modular Arithmetic ===\n");

console.log("Compute (a * b) mod M by breaking into smaller computations.");
console.log("This allows parallel processing and avoiding large numbers.\n");

// Large modulus factored into coprimes
const M4 = 3n * 5n * 7n * 11n; // 1155
const a = 456n;
const b = 789n;

console.log(`Computing ${a} * ${b} mod ${M4}`);
console.log(`Direct: ${(a * b) % M4}\n`);

console.log("Using CRT decomposition:");
const factors4 = [3n, 5n, 7n, 11n];

// Compute in each residue class
const residuesA = factors4.map(m => a % m);
const residuesB = factors4.map(m => b % m);
const residuesProduct = factors4.map((m, i) => (residuesA[i]! * residuesB[i]!) % m);

console.log("  Modulus | a mod m | b mod m | (a*b) mod m");
console.log("  --------+---------+---------+------------");
for (let i = 0; i < factors4.length; i++) {
  console.log(`     ${factors4[i]?.toString().padStart(2)}   |    ${residuesA[i]?.toString().padStart(2)}   |    ${residuesB[i]?.toString().padStart(2)}   |      ${residuesProduct[i]?.toString().padStart(2)}`);
}

const reconstructed = CRT_list(residuesProduct, factors4);
console.log(`\nReconstructed via CRT: ${reconstructed}`);
console.log(`Matches direct computation: ${reconstructed === (a * b) % M4}`);

// ============================================================================
// Application 5: Range Compression
// ============================================================================
console.log("\n=== Application 5: Representing Large Numbers ===\n");

console.log("A large number can be represented as a tuple of small residues.");
console.log("This is useful for hardware with limited word size.\n");

// Represent a large number using multiple small moduli
const largeNumber = 123456789n;
const repModuli = [251n, 253n, 255n, 257n]; // Pairwise coprime
const M5 = repModuli.reduce((a, b) => a * b); // > 4 billion

console.log(`Large number: ${largeNumber}`);
console.log(`Moduli: [${repModuli.join(", ")}]`);
console.log(`Range: 0 to ${M5 - 1n}\n`);

// Represent as residues
const representation = repModuli.map(m => largeNumber % m);
console.log(`Representation: (${representation.join(", ")})`);

// Recover
const recoveredLarge = CRT_list(representation, repModuli);
console.log(`Recovered: ${recoveredLarge}`);
console.log(`Correct: ${recoveredLarge === largeNumber}`);

// Show that arithmetic works
const addend = 987654321n;
const sumResidues = repModuli.map((m, i) => (representation[i]! + (addend % m)) % m);
const sumRecovered = CRT_list(sumResidues, repModuli);
const directSum = (largeNumber + addend) % M5;

console.log(`\nAddition example:`);
console.log(`  ${largeNumber} + ${addend} mod ${M5}`);
console.log(`  Via residues: ${sumRecovered}`);
console.log(`  Direct: ${directSum}`);
console.log(`  Match: ${sumRecovered === directSum}`);

// ============================================================================
// Application 6: Interpolation and Polynomial Reconstruction
// ============================================================================
console.log("\n=== Application 6: Lagrange Interpolation Connection ===\n");

console.log("CRT is analogous to Lagrange interpolation for polynomials:");
console.log("  - CRT: Reconstruct x from x mod m_i");
console.log("  - Lagrange: Reconstruct P(x) from P(x_i)\n");

console.log("The CRT basis elements are like Lagrange basis polynomials:");
console.log("  - CRT: e_i = 1 mod m_i, e_i = 0 mod m_j (j != i)");
console.log("  - Lagrange: L_i(x_i) = 1, L_i(x_j) = 0 (j != i)\n");

// Demonstrate the connection
const crtModuli = [3n, 5n, 7n];
const crtM = crtModuli.reduce((a, b) => a * b);

console.log(`CRT basis for moduli [${crtModuli.join(", ")}]:`);
for (let i = 0; i < crtModuli.length; i++) {
  const mi = crtModuli[i]!;
  const Mi = crtM / mi;
  const MiInv = inverse_mod(Mi, mi);
  const ei = (Mi * MiInv) % crtM;

  console.log(`\ne_${i}:`);
  console.log(`  Formula: (M/m_${i}) * (M/m_${i})^(-1) mod M`);
  console.log(`  Value: ${Mi} * ${MiInv} mod ${crtM} = ${ei}`);
  console.log(`  Evaluations: e_${i} mod m_j for j = 0, 1, 2: ` +
              `[${crtModuli.map(m => ei % m).join(", ")}]`);
}

// ============================================================================
// Application 7: Solving Linear Diophantine Equations
// ============================================================================
console.log("\n=== Application 7: Linear Diophantine Equations ===\n");

console.log("Find integers x, y such that ax + by = c.\n");

function solveDiophantine(a: bigint, b: bigint, c: bigint): {
  hasSolution: boolean;
  x?: bigint;
  y?: bigint;
  general?: string;
} {
  const [g, x0, y0] = xgcd(a, b);

  if (c % g !== 0n) {
    return { hasSolution: false };
  }

  const x = x0 * (c / g);
  const y = y0 * (c / g);

  return {
    hasSolution: true,
    x,
    y,
    general: `x = ${x} + ${b / g}k, y = ${y} - ${a / g}k for any integer k`,
  };
}

const diophantineExamples = [
  { a: 12n, b: 18n, c: 6n },
  { a: 15n, b: 21n, c: 12n },
  { a: 10n, b: 14n, c: 5n }, // No solution
];

for (const { a, b, c } of diophantineExamples) {
  console.log(`Solving ${a}x + ${b}y = ${c}:`);
  const result = solveDiophantine(a, b, c);

  if (result.hasSolution) {
    console.log(`  Particular solution: x = ${result.x}, y = ${result.y}`);
    console.log(`  Verification: ${a}*${result.x} + ${b}*${result.y} = ${a * result.x! + b * result.y!}`);
    console.log(`  General: ${result.general}`);
  } else {
    console.log(`  NO SOLUTION (gcd(${a}, ${b}) = ${gcd(a, b)} does not divide ${c})`);
  }
  console.log();
}

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
CRT Applications:

1. Solving Congruence Systems
   - Find x satisfying multiple modular conditions
   - Useful in scheduling, resource allocation

2. Calendar Calculations
   - Determine dates satisfying multiple periodic conditions
   - Day-of-week calculations

3. Secret Sharing
   - Distribute secret as residues
   - Threshold schemes for secure sharing

4. Parallel Computation
   - Break large modular arithmetic into parallel small computations
   - Each processor handles one modulus

5. Number Representation
   - Represent large numbers as tuples of small residues
   - Arithmetic operations become parallel

6. Interpolation
   - CRT is analogous to Lagrange interpolation
   - Basis elements play similar roles

7. Diophantine Equations
   - Extended GCD and CRT work together
   - Solve linear equations in integers

Key Insight: CRT transforms one problem over a large ring
into independent problems over smaller rings that can be
solved separately and combined.
`);
