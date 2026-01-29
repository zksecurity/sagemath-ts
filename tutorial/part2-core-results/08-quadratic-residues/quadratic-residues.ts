/**
 * # Quadratic Residues
 *
 * A quadratic residue modulo n is a number that is a perfect square
 * modulo n. Understanding quadratic residues is fundamental to many
 * cryptographic constructions.
 *
 * ## Definition
 *
 * An integer a is a quadratic residue modulo n if there exists an
 * integer x such that:
 *   x^2 = a (mod n)
 *
 * If no such x exists (and gcd(a, n) = 1), then a is a quadratic
 * non-residue modulo n.
 *
 * @module tutorial/part2/quadratic-residues
 */

import {
  power_mod,
  gcd,
  is_prime,
  euler_phi,
  legendre_symbol,
  Zmod,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("QUADRATIC RESIDUES");
console.log("=".repeat(70));

// ============================================================================
// Example 1: Definition and Examples
// ============================================================================
console.log("\n=== Example 1: Definition and Examples ===\n");

console.log("A quadratic residue (QR) mod n is a perfect square mod n.");
console.log("A quadratic non-residue (QNR) mod n has no square root.\n");

function findQuadraticResidues(n: bigint): { QR: Set<bigint>; QNR: bigint[] } {
  const QR = new Set<bigint>();

  // Find all squares mod n
  for (let x = 0n; x < n; x++) {
    const square = (x * x) % n;
    QR.add(square);
  }

  // QNR are units not in QR
  const QNR: bigint[] = [];
  for (let a = 1n; a < n; a++) {
    if (gcd(a, n) === 1n && !QR.has(a)) {
      QNR.push(a);
    }
  }

  return { QR, QNR };
}

// Example with p = 7
const p7 = 7n;
console.log(`Modulo ${p7} (prime):`);
console.log("-".repeat(50));

// Show all squares
console.log("\nSquares x^2 mod 7:");
for (let x = 0n; x < p7; x++) {
  const square = (x * x) % p7;
  console.log(`  ${x}^2 = ${x * x} = ${square} (mod ${p7})`);
}

const { QR: QR7, QNR: QNR7 } = findQuadraticResidues(p7);
console.log(`\nQuadratic Residues (QR): {${Array.from(QR7).sort((a, b) => a < b ? -1 : 1).join(", ")}}`);
console.log(`Quadratic Non-Residues (QNR): {${QNR7.join(", ")}}`);

// ============================================================================
// Example 2: QR for Various Primes
// ============================================================================
console.log("\n=== Example 2: QR for Various Primes ===\n");

const primes = [5n, 7n, 11n, 13n];

for (const p of primes) {
  const { QR, QNR } = findQuadraticResidues(p);
  const qrArray = Array.from(QR).filter(x => x !== 0n).sort((a, b) => a < b ? -1 : 1);

  console.log(`p = ${p}:`);
  console.log(`  QR (excluding 0): {${qrArray.join(", ")}} (${qrArray.length} elements)`);
  console.log(`  QNR: {${QNR.join(", ")}} (${QNR.length} elements)`);
  console.log(`  |QR| = |QNR| = (p-1)/2 = ${(p - 1n) / 2n}\n`);
}

console.log("Theorem: For odd prime p, exactly half of the units are QR.");

// ============================================================================
// Example 3: Using Euler's Criterion
// ============================================================================
console.log("\n=== Example 3: Euler's Criterion ===\n");

console.log("Euler's Criterion: For odd prime p and gcd(a, p) = 1,");
console.log("  a is QR iff a^((p-1)/2) = 1 (mod p)");
console.log("  a is QNR iff a^((p-1)/2) = -1 (mod p)\n");

const p11 = 11n;
const exp = (p11 - 1n) / 2n; // 5

console.log(`Testing mod ${p11} using exponent (p-1)/2 = ${exp}:`);
console.log("-".repeat(50));
console.log("  a | a^5 mod 11 | QR/QNR | Verification");
console.log("----+------------+--------+-------------");

for (let a = 1n; a < p11; a++) {
  const eulerValue = power_mod(a, exp, p11);
  const isQR = eulerValue === 1n;

  // Verify by trying to find a square root
  let sqrtFound: bigint | null = null;
  for (let x = 0n; x < p11; x++) {
    if ((x * x) % p11 === a) {
      sqrtFound = x;
      break;
    }
  }

  const status = isQR === (sqrtFound !== null) ? "OK" : "FAIL";
  const sqrtStr = sqrtFound !== null ? `sqrt = ${sqrtFound}` : "no sqrt";

  console.log(` ${a.toString().padStart(2)} |     ${eulerValue === 1n ? " 1" : "-1 (=10)"}    |   ${isQR ? "QR " : "QNR"}  | ${sqrtStr.padEnd(10)} [${status}]`);
}

// ============================================================================
// Example 4: QR Modulo Composite Numbers
// ============================================================================
console.log("\n=== Example 4: QR Modulo Composite Numbers ===\n");

console.log("For composite n = p*q, the situation is more complex.");
console.log("a is QR mod n iff a is QR mod p AND a is QR mod q.\n");

const pp = 3n;
const qq = 5n;
const nn = pp * qq;

console.log(`Consider n = ${pp} * ${qq} = ${nn}`);

const { QR: QR15 } = findQuadraticResidues(nn);
const { QR: QR3 } = findQuadraticResidues(pp);
const { QR: QR5 } = findQuadraticResidues(qq);

console.log(`\nQR mod ${pp}: {${Array.from(QR3).sort((a, b) => a < b ? -1 : 1).join(", ")}}`);
console.log(`QR mod ${qq}: {${Array.from(QR5).sort((a, b) => a < b ? -1 : 1).join(", ")}}`);
console.log(`QR mod ${nn}: {${Array.from(QR15).sort((a, b) => a < b ? -1 : 1).join(", ")}}`);

console.log("\nDetailed analysis:");
console.log(" a | a mod 3 | a mod 5 | QR(3)? | QR(5)? | QR(15)?");
console.log("---+---------+---------+--------+--------+--------");

for (let a = 1n; a < nn; a++) {
  if (gcd(a, nn) === 1n) {
    const amod3 = a % pp;
    const amod5 = a % qq;
    const qr3 = QR3.has(amod3);
    const qr5 = QR5.has(amod5);
    const qr15 = QR15.has(a);
    const expected = qr3 && qr5;

    console.log(`${a.toString().padStart(2)} |    ${amod3}    |    ${amod5}    |   ${qr3 ? "Y" : "N"}    |   ${qr5 ? "Y" : "N"}    |   ${qr15 ? "Y" : "N"}${qr15 !== expected ? " (!)" : ""}`);
  }
}

// ============================================================================
// Example 5: Number of Square Roots
// ============================================================================
console.log("\n=== Example 5: Number of Square Roots ===\n");

console.log("For prime p, each QR has exactly 2 square roots: x and p-x.");
console.log("For n = p*q, each QR has 4 square roots (by CRT).\n");

// Prime case
const p13 = 13n;
console.log(`Square roots mod ${p13} (prime):`);
for (const a of [1n, 3n, 4n, 9n, 12n]) {
  const roots: bigint[] = [];
  for (let x = 0n; x < p13; x++) {
    if ((x * x) % p13 === a) {
      roots.push(x);
    }
  }
  console.log(`  sqrt(${a.toString().padStart(2)}) mod ${p13}: {${roots.join(", ")}} (${roots.length} roots)`);
}

// Composite case
const n21 = 3n * 7n;
console.log(`\nSquare roots mod ${n21} = 3*7 (composite):`);
for (const a of [1n, 4n, 16n]) {
  const roots: bigint[] = [];
  for (let x = 0n; x < n21; x++) {
    if ((x * x) % n21 === a) {
      roots.push(x);
    }
  }
  console.log(`  sqrt(${a.toString().padStart(2)}) mod ${n21}: {${roots.join(", ")}} (${roots.length} roots)`);
}

// ============================================================================
// Example 6: QR and Cryptographic Security
// ============================================================================
console.log("\n=== Example 6: QR and Cryptographic Security ===\n");

console.log("The Quadratic Residuosity Problem:");
console.log("  Given n = p*q and a coprime to n, determine if a is QR mod n.\n");

console.log("This is believed to be computationally hard when:");
console.log("  - The factorization of n is unknown");
console.log("  - p and q are large primes\n");

console.log("Applications:");
console.log("  - Goldwasser-Micali encryption");
console.log("  - Blum-Blum-Shub PRNG");
console.log("  - Various zero-knowledge proofs\n");

// Demonstrate the difficulty
const p_large = 101n;
const q_large = 103n;
const n_large = p_large * q_large;

console.log(`Example: n = ${p_large} * ${q_large} = ${n_large}`);
console.log("Without knowing p and q, determining if a is QR requires");
console.log("either factoring n or using quantum computers.\n");

const testA = 1234n;
const isQR_p = legendre_symbol(testA, p_large) === 1n;
const isQR_q = legendre_symbol(testA, q_large) === 1n;
const isQR_n = isQR_p && isQR_q;

console.log(`a = ${testA}:`);
console.log(`  Legendre(${testA}, ${p_large}) = ${legendre_symbol(testA, p_large)} (QR: ${isQR_p})`);
console.log(`  Legendre(${testA}, ${q_large}) = ${legendre_symbol(testA, q_large)} (QR: ${isQR_q})`);
console.log(`  Is ${testA} a QR mod ${n_large}? ${isQR_n}`);
console.log("\n  (This computation requires knowing the factors p and q!)");

// ============================================================================
// Example 7: Quadratic Character Table
// ============================================================================
console.log("\n=== Example 7: Quadratic Character Table ===\n");

const pTable = 11n;
const Z11 = Zmod(pTable);

console.log(`Quadratic character table for Z/${pTable}Z:`);
console.log("-".repeat(50));

// Build multiplication table showing QR pattern
console.log("\nMultiplication table showing products that are QR:");
console.log("   |  1  2  3  4  5  6  7  8  9 10");
console.log("---+-------------------------------");

const qrSet = findQuadraticResidues(pTable).QR;

for (let i = 1n; i < pTable; i++) {
  let row = `${i.toString().padStart(2)} |`;
  for (let j = 1n; j < pTable; j++) {
    const prod = (i * j) % pTable;
    row += qrSet.has(prod) ? "  Q" : "  .";
  }
  console.log(row);
}

console.log("\nLegend: Q = quadratic residue, . = non-residue");
console.log("\nPattern: QR * QR = QR, QR * QNR = QNR, QNR * QNR = QR");
console.log("This is because (Z/pZ)* / QR is isomorphic to {+1, -1}.");

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
Quadratic Residues:

Definition:
  a is a quadratic residue (QR) mod n if x^2 = a (mod n) has a solution.
  Otherwise (if gcd(a,n) = 1), a is a quadratic non-residue (QNR).

For odd prime p:
  - Exactly (p-1)/2 units are QR
  - Exactly (p-1)/2 units are QNR
  - 0 is considered a degenerate case

Euler's Criterion:
  a^((p-1)/2) = 1 (mod p)  =>  a is QR
  a^((p-1)/2) = -1 (mod p) =>  a is QNR

Square roots:
  - mod p: Each QR has exactly 2 square roots (x and p-x)
  - mod p*q: Each QR has 4 square roots (by CRT)
  - mod p^k: More complex structure

Cryptographic Importance:
  - Quadratic Residuosity Problem is believed hard
  - Foundation for several cryptosystems
  - Related to discrete log and factoring
`);
