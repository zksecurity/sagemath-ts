/**
 * # The Jacobi Symbol
 *
 * The Jacobi symbol generalizes the Legendre symbol to composite moduli.
 * It is essential for efficient computation without factoring.
 *
 * ## Definition
 *
 * For positive odd integer n with prime factorization n = p1^e1 * ... * pk^ek,
 * and integer a, the Jacobi symbol (a/n) is:
 *
 *   (a/n) = (a/p1)^e1 * (a/p2)^e2 * ... * (a/pk)^ek
 *
 * where (a/pi) are Legendre symbols.
 *
 * ## Important Note
 *
 * (a/n) = 1 does NOT mean a is a quadratic residue mod n!
 * The Jacobi symbol only indicates QR status when n is prime.
 *
 * @module tutorial/part2/jacobi-symbol
 */

import {
  jacobi_symbol,
  legendre_symbol,
  kronecker_symbol,
  is_prime,
  factor,
  power_mod,
  gcd,
  formatFactorization,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("THE JACOBI SYMBOL");
console.log("=".repeat(70));

// ============================================================================
// Example 1: Definition via Prime Factorization
// ============================================================================
console.log("\n=== Example 1: Definition via Prime Factorization ===\n");

console.log("For n = p1^e1 * ... * pk^ek:");
console.log("  (a/n) = (a/p1)^e1 * (a/p2)^e2 * ... * (a/pk)^ek\n");

const n1 = 15n; // 3 * 5
const factors1 = factor(n1);
console.log(`n = ${n1} = ${formatFactorization(factors1)}`);

console.log(`\nComputing (a/${n1}) for various a:`);
console.log("-".repeat(60));

for (let a = 1n; a <= 14n; a++) {
  if (gcd(a, n1) !== 1n) continue;

  // Compute via factorization
  const leg3 = legendre_symbol(a, 3n);
  const leg5 = legendre_symbol(a, 5n);
  const productValue = leg3 * leg5;

  // Compute directly
  const jacobi = jacobi_symbol(a, n1);

  console.log(`  (${a.toString().padStart(2)}/${n1}) = (${a}/3) * (${a}/5) = ${leg3.toString().padStart(2)} * ${leg5.toString().padStart(2)} = ${productValue.toString().padStart(2)}  [jacobi_symbol: ${jacobi}]`);
}

// ============================================================================
// Example 2: Jacobi != QR Indicator for Composite n
// ============================================================================
console.log("\n=== Example 2: Jacobi Symbol vs Quadratic Residuosity ===\n");

console.log("IMPORTANT: (a/n) = 1 does NOT mean a is QR mod n!\n");

const n2 = 15n;
console.log(`Consider n = ${n2} = 3 * 5`);

// Find actual quadratic residues
function isQR(a: bigint, n: bigint): boolean {
  for (let x = 0n; x < n; x++) {
    if ((x * x) % n === a) return true;
  }
  return false;
}

console.log("\n a | (a/15) | is QR mod 15?");
console.log("---+--------+--------------");

let contradictions = 0;
for (let a = 1n; a < n2; a++) {
  if (gcd(a, n2) !== 1n) continue;

  const jacobi = jacobi_symbol(a, n2);
  const qr = isQR(a, n2);

  const contradiction = (jacobi === 1n) && !qr;
  if (contradiction) contradictions++;

  console.log(`${a.toString().padStart(2)} |   ${jacobi.toString().padStart(2)}   | ${qr ? "Yes" : "No "}${contradiction ? " <-- (a/n)=1 but not QR!" : ""}`);
}

console.log(`\nFound ${contradictions} values where (a/n) = 1 but a is NOT a QR.`);
console.log("This is why Jacobi symbol alone cannot determine quadratic residuosity.");

// ============================================================================
// Example 3: Properties of the Jacobi Symbol
// ============================================================================
console.log("\n=== Example 3: Properties of the Jacobi Symbol ===\n");

const n3 = 21n; // 3 * 7

console.log(`Properties for n = ${n3}:\n`);

// Property 1: Multiplicativity in top
console.log("1. Multiplicativity in numerator: (ab/n) = (a/n)(b/n)");
console.log("-".repeat(50));

const pairs = [[2n, 5n], [4n, 11n], [8n, 13n]];
for (const [a, b] of pairs) {
  const ab = (a * b) % (2n * n3);  // Keep in reasonable range
  const jacAB = jacobi_symbol(ab, n3);
  const jacA = jacobi_symbol(a, n3);
  const jacB = jacobi_symbol(b, n3);
  const product = jacA * jacB;

  console.log(`  (${a}*${b}/${n3}) = (${ab}/${n3}) = ${jacAB}`);
  console.log(`  (${a}/${n3}) * (${b}/${n3}) = ${jacA} * ${jacB} = ${product} [${jacAB === product ? "OK" : "FAIL"}]\n`);
}

// Property 2: Multiplicativity in bottom
console.log("\n2. Multiplicativity in denominator: (a/mn) = (a/m)(a/n)");
console.log("-".repeat(50));

const m2 = 9n;
const n2b = 5n;
const testA = 7n;

const jacMN = jacobi_symbol(testA, m2 * n2b);
const jacM = jacobi_symbol(testA, m2);
const jacN = jacobi_symbol(testA, n2b);

console.log(`  (${testA}/${m2 * n2b}) = ${jacMN}`);
console.log(`  (${testA}/${m2}) * (${testA}/${n2b}) = ${jacM} * ${jacN} = ${jacM * jacN}`);
console.log(`  Match: ${jacMN === jacM * jacN}`);

// Property 3: (a/n) depends only on a mod n
console.log("\n\n3. Reduction: (a/n) = (a mod n / n)");
console.log("-".repeat(50));

const testValues = [100n, -4n, 256n];
for (const a of testValues) {
  const reduced = ((a % n3) + n3) % n3;
  const jacA = jacobi_symbol(a, n3);
  const jacReduced = reduced !== 0n && gcd(reduced, n3) === 1n ? jacobi_symbol(reduced, n3) : 0n;

  if (gcd(a, n3) === 1n) {
    console.log(`  (${a}/${n3}) = (${reduced}/${n3}) = ${jacA}`);
  }
}

// ============================================================================
// Example 4: Efficient Computation Algorithm
// ============================================================================
console.log("\n=== Example 4: Efficient Computation (without factoring) ===\n");

console.log("The Jacobi symbol can be computed efficiently using:");
console.log("  1. (a/n) = (a mod n / n)");
console.log("  2. (ab/n) = (a/n)(b/n)");
console.log("  3. (-1/n) = (-1)^((n-1)/2)");
console.log("  4. (2/n) = (-1)^((n^2-1)/8)");
console.log("  5. Quadratic reciprocity: (m/n) = (-1)^((m-1)(n-1)/4) * (n/m)\n");

// Manual computation showing steps
function jacobiWithSteps(a: bigint, n: bigint): { result: bigint; steps: string[] } {
  const steps: string[] = [];
  let result = 1n;

  a = ((a % n) + n) % n;
  steps.push(`Reduce: (${a}/${n})`);

  while (a !== 0n && a !== 1n) {
    // Handle factors of 2
    let twos = 0n;
    while (a % 2n === 0n) {
      a /= 2n;
      twos++;
    }
    if (twos > 0n) {
      // (2/n) = (-1)^((n^2-1)/8)
      const twoSymbol = (n % 8n === 1n || n % 8n === 7n) ? 1n : -1n;
      if (twos % 2n === 1n) {
        result *= twoSymbol;
      }
      steps.push(`  Extract ${twos} factors of 2: (2/${n})^${twos} contributes ${twos % 2n === 1n ? twoSymbol : 1n}`);
      steps.push(`  Now: (${a}/${n})`);
    }

    if (a === 1n) break;

    // Quadratic reciprocity: swap a and n
    const recipSign = (a % 4n === 3n && n % 4n === 3n) ? -1n : 1n;
    result *= recipSign;

    [a, n] = [n % a, a];
    steps.push(`  Reciprocity (sign: ${recipSign}): now (${a}/${n})`);
  }

  steps.push(`Final result: ${result}`);
  return { result, steps };
}

const testA2 = 123n;
const testN2 = 455n; // 5 * 7 * 13

console.log(`Computing (${testA2}/${testN2}):`);
const { result: manualResult, steps } = jacobiWithSteps(testA2, testN2);
for (const step of steps) {
  console.log(`  ${step}`);
}

console.log(`\nLibrary result: ${jacobi_symbol(testA2, testN2)}`);
console.log(`Match: ${manualResult === jacobi_symbol(testA2, testN2)}`);

// ============================================================================
// Example 5: Comparison: Legendre vs Jacobi vs Kronecker
// ============================================================================
console.log("\n=== Example 5: Legendre vs Jacobi vs Kronecker ===\n");

console.log("Symbol hierarchy:");
console.log("  - Legendre (a/p): p must be odd prime");
console.log("  - Jacobi (a/n): n must be odd positive integer");
console.log("  - Kronecker (a/n): n can be any integer (generalization)\n");

const values = [
  { a: 3n, n: 7n, desc: "prime n" },
  { a: 3n, n: 15n, desc: "odd composite n" },
  { a: 3n, n: 8n, desc: "even n (Kronecker only)" },
  { a: -3n, n: 5n, desc: "negative a" },
  { a: 5n, n: -7n, desc: "negative n (Kronecker only)" },
];

console.log("  a   |   n   | Legendre | Jacobi | Kronecker | Notes");
console.log("------+-------+----------+--------+-----------+------");

for (const { a, n, desc } of values) {
  const absN = n < 0n ? -n : n;

  let leg = "-";
  if (is_prime(absN) && absN > 2n && n > 0n) {
    leg = legendre_symbol(a, absN).toString();
  }

  let jac = "-";
  if (n > 0n && n % 2n !== 0n) {
    jac = jacobi_symbol(a, n).toString();
  }

  const kron = kronecker_symbol(a, n);

  console.log(`${a.toString().padStart(4)}  |${n.toString().padStart(4)}   |    ${leg.padStart(2)}    |   ${jac.padStart(2)}   |     ${kron.toString().padStart(2)}    | ${desc}`);
}

// ============================================================================
// Example 6: Application - Solovay-Strassen Primality Test
// ============================================================================
console.log("\n=== Example 6: Solovay-Strassen Primality Test ===\n");

console.log("If n is prime, then (a/n) = a^((n-1)/2) mod n for all a coprime to n.");
console.log("If this equality fails for some a, then n is definitely composite.\n");

function solovayStrassen(n: bigint, rounds: number = 10): { isProbablyPrime: boolean; witness?: bigint } {
  if (n < 2n) return { isProbablyPrime: false };
  if (n === 2n) return { isProbablyPrime: true };
  if (n % 2n === 0n) return { isProbablyPrime: false };

  for (let i = 0; i < rounds; i++) {
    // Random a in [2, n-2]
    const a = 2n + BigInt(Math.floor(Math.random() * Math.min(Number(n) - 4, 10000)));

    if (gcd(a, n) > 1n) {
      return { isProbablyPrime: false, witness: a };
    }

    const jacobi = jacobi_symbol(a, n);
    const euler = power_mod(a, (n - 1n) / 2n, n);

    // Convert jacobi to match euler's range [0, n-1]
    const jacobiMod = jacobi === -1n ? n - 1n : (jacobi === 0n ? 0n : 1n);

    if (euler !== jacobiMod) {
      return { isProbablyPrime: false, witness: a };
    }
  }

  return { isProbablyPrime: true };
}

console.log("Testing numbers:");
const testNums = [7n, 15n, 17n, 91n, 561n, 1009n, 1105n];

for (const n of testNums) {
  const actualPrime = is_prime(n);
  const ssResult = solovayStrassen(n, 20);

  let status = "";
  if (actualPrime && ssResult.isProbablyPrime) status = "OK";
  else if (!actualPrime && !ssResult.isProbablyPrime) status = "OK (detected composite)";
  else if (!actualPrime && ssResult.isProbablyPrime) status = "FOOLED";
  else status = "ERROR";

  console.log(`  n = ${n.toString().padStart(4)}: SS says ${ssResult.isProbablyPrime ? "prime" : "composite"}, ` +
              `actually ${actualPrime ? "prime" : "composite"} [${status}]`);
}

console.log("\nNote: 561 = 3*11*17 is a Carmichael number and may fool the test");
console.log("with some witnesses, but Solovay-Strassen is better than Fermat test.");

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
The Jacobi Symbol (a/n):

Definition (n odd positive):
  (a/n) = product of (a/p_i)^e_i for prime factorization n = prod p_i^e_i

Key Properties:
  1. Extends Legendre symbol to composite moduli
  2. Multiplicativity: (ab/n) = (a/n)(b/n), (a/mn) = (a/m)(a/n)
  3. (a/n) = (a mod n / n)
  4. (-1/n) = (-1)^((n-1)/2)
  5. (2/n) = (-1)^((n^2-1)/8)
  6. Quadratic reciprocity: (m/n)(n/m) = (-1)^((m-1)(n-1)/4) for odd m,n > 0

IMPORTANT:
  - (a/n) = 1 does NOT imply a is QR mod n (for composite n)
  - (a/n) = -1 DOES imply a is NOT QR mod n

Computation:
  - Can be computed in O(log n) without factoring n
  - Uses properties and quadratic reciprocity

Applications:
  - Solovay-Strassen primality test
  - Various cryptographic protocols
  - Efficient computation without factoring
`);
