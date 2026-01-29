/**
 * # Extended Euclidean Algorithm and Bezout's Identity
 *
 * The Extended Euclidean Algorithm is one of the most important algorithms
 * in cryptography. It not only computes gcd(a, b) but also finds integers
 * s and t such that: gcd(a, b) = s*a + t*b (Bezout's Identity)
 *
 * ## Mathematical Background
 *
 * Bezout's Identity: For any integers a and b, there exist integers s and t
 * such that: gcd(a, b) = s*a + t*b
 *
 * The coefficients s and t are called Bezout coefficients.
 * They are not unique, but the extended Euclidean algorithm finds one pair.
 *
 * ## Why This Matters for Cryptography
 *
 * The Extended Euclidean Algorithm is THE fundamental algorithm for:
 *   - Computing modular inverses (essential for RSA decryption)
 *   - Solving linear Diophantine equations
 *   - Chinese Remainder Theorem computations
 *   - Lattice basis reduction
 */

import { gcd, xgcd, inverse_mod, Mod } from 'sagemath-ts';

// =============================================================================
// Example 1: Basic Extended GCD
// =============================================================================
console.log("=== Example 1: Extended GCD Basics ===\n");

/**
 * The xgcd() function returns [g, s, t] where:
 *   g = gcd(a, b)
 *   g = s*a + t*b
 */
const pairs: [bigint, bigint][] = [
  [48n, 18n],
  [240n, 46n],
  [35n, 15n],
  [101n, 103n],  // Coprime numbers
];

console.log("Extended GCD computations:");
for (const [a, b] of pairs) {
  const [g, s, t] = xgcd(a, b);
  console.log(`\nxgcd(${a}, ${b}):`);
  console.log(`  gcd = ${g}`);
  console.log(`  Bezout coefficients: s = ${s}, t = ${t}`);
  console.log(`  Verification: ${s}*${a} + ${t}*${b} = ${s*a + t*b}`);
}

console.log();

// =============================================================================
// Example 2: Tracing the Extended Euclidean Algorithm
// =============================================================================
console.log("=== Example 2: Step-by-Step Extended Euclidean Algorithm ===\n");

/**
 * Let's trace through the extended Euclidean algorithm to see how
 * the Bezout coefficients are computed.
 */
function extendedEuclideanTrace(a: bigint, b: bigint): [bigint, bigint, bigint] {
  console.log(`Computing xgcd(${a}, ${b}) step by step:\n`);

  // Track the equations: r_i = s_i * a + t_i * b
  let oldR = a, r = b;
  let oldS = 1n, s = 0n;
  let oldT = 0n, t = 1n;

  console.log(`  Initial: ${oldR} = ${oldS}*${a} + ${oldT}*${b}`);
  console.log(`           ${r} = ${s}*${a} + ${t}*${b}\n`);

  let step = 1;
  while (r !== 0n) {
    const quotient = oldR / r;

    // Update remainders
    const tempR = r;
    r = oldR - quotient * r;
    oldR = tempR;

    // Update s coefficients
    const tempS = s;
    s = oldS - quotient * s;
    oldS = tempS;

    // Update t coefficients
    const tempT = t;
    t = oldT - quotient * t;
    oldT = tempT;

    console.log(`  Step ${step}: quotient = ${quotient}`);
    console.log(`           ${oldR} = ${oldS}*${a} + ${oldT}*${b}`);
    if (r !== 0n) {
      console.log(`           ${r} = ${s}*${a} + ${t}*${b}`);
    }
    console.log();
    step++;
  }

  // Handle negative GCD
  if (oldR < 0n) {
    oldR = -oldR;
    oldS = -oldS;
    oldT = -oldT;
  }

  console.log(`  Result: gcd = ${oldR}, s = ${oldS}, t = ${oldT}`);
  return [oldR, oldS, oldT];
}

extendedEuclideanTrace(240n, 46n);

console.log();

// =============================================================================
// Example 3: Computing Modular Inverses
// =============================================================================
console.log("=== Example 3: Modular Inverses via Extended GCD ===\n");

/**
 * The modular inverse of a modulo n exists iff gcd(a, n) = 1.
 * When it exists, we can find it using the extended Euclidean algorithm:
 *
 * If gcd(a, n) = 1 = s*a + t*n
 * Then s*a = 1 - t*n = 1 (mod n)
 * So s is the modular inverse of a modulo n.
 */
function findModularInverse(a: bigint, n: bigint): bigint | null {
  const [g, s, _t] = xgcd(a, n);

  if (g !== 1n) {
    return null;  // No inverse exists
  }

  // Ensure result is positive
  return ((s % n) + n) % n;
}

console.log("Computing modular inverses:");

const invTests: [bigint, bigint][] = [
  [3n, 7n],    // 3 * 5 = 15 = 1 (mod 7)
  [5n, 12n],   // 5 * 5 = 25 = 1 (mod 12)
  [17n, 43n],  // Both prime
  [6n, 9n],    // gcd(6, 9) = 3 != 1, no inverse
];

for (const [a, n] of invTests) {
  const inv = findModularInverse(a, n);
  const [g, _s, _t] = xgcd(a, n);

  if (inv !== null) {
    console.log(`  ${a}^(-1) mod ${n} = ${inv}`);
    console.log(`    Verification: ${a} * ${inv} = ${a * inv} = ${(a * inv) % n} (mod ${n})\n`);
  } else {
    console.log(`  ${a}^(-1) mod ${n} does not exist (gcd = ${g})\n`);
  }
}

// Using the library's inverse_mod function
console.log("Using library inverse_mod function:");
console.log(`  inverse_mod(3, 7) = ${inverse_mod(3n, 7n)}`);
console.log(`  inverse_mod(17, 43) = ${inverse_mod(17n, 43n)}`);

console.log();

// =============================================================================
// Example 4: RSA Decryption Key Computation
// =============================================================================
console.log("=== Example 4: RSA Key Generation ===\n");

/**
 * In RSA, the private exponent d is the modular inverse of e mod phi(n).
 * This is computed using the extended Euclidean algorithm.
 */
const p = 61n;
const q = 53n;
const n = p * q;
const phiN = (p - 1n) * (q - 1n);
const e = 17n;  // Public exponent

console.log("RSA Key Generation:");
console.log(`  p = ${p}, q = ${q}`);
console.log(`  n = ${n}`);
console.log(`  phi(n) = ${phiN}`);
console.log(`  e = ${e}`);

// Verify e and phi(n) are coprime
const gcdResult = gcd(e, phiN);
console.log(`  gcd(e, phi(n)) = ${gcdResult}`);

if (gcdResult === 1n) {
  // Compute d = e^(-1) mod phi(n)
  const [_g, s, _t] = xgcd(e, phiN);
  const d = ((s % phiN) + phiN) % phiN;

  console.log(`  d = e^(-1) mod phi(n) = ${d}`);
  console.log(`\n  Verification: e * d mod phi(n) = ${(e * d) % phiN}`);

  // Test encryption/decryption
  const message = 42n;
  const cipher = message ** e % n;
  const decrypted = cipher ** d % n;  // Note: This is slow for demo purposes

  console.log(`\n  Test encryption/decryption:`);
  console.log(`    Message m = ${message}`);
  console.log(`    Encrypted: m^e mod n = ${cipher}`);
  console.log(`    Decrypted: c^d mod n = ${decrypted}`);
}

console.log();

// =============================================================================
// Example 5: Linear Diophantine Equations
// =============================================================================
console.log("=== Example 5: Linear Diophantine Equations ===\n");

/**
 * A linear Diophantine equation has the form: ax + by = c
 * where a, b, c are integers and we seek integer solutions (x, y).
 *
 * Solution exists iff gcd(a, b) divides c.
 * If (x0, y0) is a particular solution, general solution is:
 *   x = x0 + (b/g) * k
 *   y = y0 - (a/g) * k
 * for any integer k.
 */
function solveDiophantine(a: bigint, b: bigint, c: bigint): [bigint, bigint] | null {
  const [g, s, t] = xgcd(a, b);

  if (c % g !== 0n) {
    return null;  // No solution
  }

  const multiplier = c / g;
  const x0 = s * multiplier;
  const y0 = t * multiplier;

  return [x0, y0];
}

console.log("Solving linear Diophantine equations ax + by = c:\n");

const equations: [bigint, bigint, bigint][] = [
  [12n, 8n, 20n],     // gcd(12,8)=4 divides 20
  [15n, 6n, 9n],      // gcd(15,6)=3 divides 9
  [14n, 21n, 10n],    // gcd(14,21)=7 does NOT divide 10
];

for (const [a, b, c] of equations) {
  console.log(`Equation: ${a}x + ${b}y = ${c}`);
  const solution = solveDiophantine(a, b, c);
  const g = gcd(a, b);

  if (solution) {
    const [x, y] = solution;
    console.log(`  gcd(${a}, ${b}) = ${g} divides ${c} - Solution exists!`);
    console.log(`  Particular solution: x = ${x}, y = ${y}`);
    console.log(`  Verification: ${a}*${x} + ${b}*${y} = ${a*x + b*y}`);
    console.log(`  General solution: x = ${x} + ${b/g}k, y = ${y} - ${a/g}k`);
  } else {
    console.log(`  gcd(${a}, ${b}) = ${g} does NOT divide ${c} - No solution!`);
  }
  console.log();
}

// =============================================================================
// Example 6: Bezout Coefficients are Not Unique
// =============================================================================
console.log("=== Example 6: Multiple Bezout Representations ===\n");

/**
 * For gcd(a, b) = g = s*a + t*b, the coefficients (s, t) are not unique.
 * All solutions are: (s + k*b/g, t - k*a/g) for any integer k.
 */
const [aa, bb] = [48n, 18n];
const [gg, ss, tt] = xgcd(aa, bb);

console.log(`gcd(${aa}, ${bb}) = ${gg}`);
console.log(`\nAll Bezout representations (s + k*(${bb}/${gg}), t - k*(${aa}/${gg})):\n`);

for (let k = -3n; k <= 3n; k++) {
  const newS = ss + k * (bb / gg);
  const newT = tt - k * (aa / gg);
  console.log(`  k=${k}: ${newS}*${aa} + ${newT}*${bb} = ${newS*aa + newT*bb}`);
}

console.log();

// =============================================================================
// Example 7: Using IntegerMod for Modular Inverses
// =============================================================================
console.log("=== Example 7: IntegerMod Inverse Method ===\n");

/**
 * The IntegerMod class provides an inv() method that uses
 * the extended Euclidean algorithm internally.
 */
const modulus = 17n;
console.log(`Working in Z/${modulus}Z:\n`);

for (let i = 1n; i < modulus; i++) {
  const element = Mod(i, modulus);
  const inverse = element.inv();
  const product = element.mul(inverse);
  console.log(`  ${i}^(-1) = ${inverse.value}  (${i} * ${inverse.value} = ${product.value} mod ${modulus})`);
}

console.log();

// =============================================================================
// Exercise: Implement Extended Euclidean Algorithm
// =============================================================================
console.log("=== Exercise: Implement Extended GCD ===\n");

/**
 * EXERCISE: Implement the extended Euclidean algorithm yourself.
 * Return [gcd, s, t] such that gcd = s*a + t*b.
 */

// TODO: Implement this function
function myXgcd(a: bigint, b: bigint): [bigint, bigint, bigint] {
  let oldR = a, r = b;
  let oldS = 1n, s = 0n;
  let oldT = 0n, t = 1n;

  while (r !== 0n) {
    const quotient = oldR / r;

    [oldR, r] = [r, oldR - quotient * r];
    [oldS, s] = [s, oldS - quotient * s];
    [oldT, t] = [t, oldT - quotient * t];
  }

  // Ensure GCD is positive
  if (oldR < 0n) {
    return [-oldR, -oldS, -oldT];
  }

  return [oldR, oldS, oldT];
}

console.log("Testing your extended GCD implementation:");
const testXgcd: [bigint, bigint][] = [[99n, 78n], [1234n, 5678n], [0n, 15n]];
for (const [a, b] of testXgcd) {
  const [g1, s1, t1] = myXgcd(a, b);
  const [g2, s2, t2] = xgcd(a, b);

  console.log(`\n  xgcd(${a}, ${b}):`);
  console.log(`    Yours:   [${g1}, ${s1}, ${t1}] -> ${s1}*${a} + ${t1}*${b} = ${s1*a + t1*b}`);
  console.log(`    Library: [${g2}, ${s2}, ${t2}] -> ${s2}*${a} + ${t2}*${b} = ${s2*a + t2*b}`);
  console.log(`    GCD match: ${g1 === g2}`);
}

console.log();

// =============================================================================
// Cryptographic Application: Batch Modular Inverse
// =============================================================================
console.log("=== Cryptographic Application: Montgomery's Trick ===\n");

/**
 * Montgomery's trick allows computing n modular inverses with only
 * ONE extended GCD computation plus O(n) multiplications.
 *
 * This is useful when you need many inverses modulo the same number.
 */
function batchInverse(values: bigint[], modulus: bigint): bigint[] {
  const n = values.length;
  if (n === 0) return [];

  // Compute cumulative products: c[i] = values[0] * values[1] * ... * values[i]
  const cumulative: bigint[] = new Array(n);
  cumulative[0] = values[0]! % modulus;
  for (let i = 1; i < n; i++) {
    cumulative[i] = (cumulative[i-1]! * values[i]!) % modulus;
  }

  // Compute inverse of the total product (ONE extended GCD)
  const [g, s] = xgcd(cumulative[n-1]!, modulus);
  if (g !== 1n) throw new Error("Not all values are invertible");
  let totalInv = ((s % modulus) + modulus) % modulus;

  // Work backwards to extract individual inverses
  const inverses: bigint[] = new Array(n);
  for (let i = n - 1; i > 0; i--) {
    // inverse[i] = totalInv * cumulative[i-1]
    inverses[i] = (totalInv * cumulative[i-1]!) % modulus;
    // totalInv becomes inverse of cumulative[i-1]
    totalInv = (totalInv * values[i]!) % modulus;
  }
  inverses[0] = totalInv;

  return inverses;
}

const mod = 101n;
const vals = [3n, 7n, 13n, 19n, 29n];
console.log(`Computing batch inverses mod ${mod} for: [${vals.join(', ')}]`);

const invs = batchInverse(vals, mod);
console.log(`\nResults (using Montgomery's trick):`);
for (let i = 0; i < vals.length; i++) {
  const product = (vals[i]! * invs[i]!) % mod;
  console.log(`  ${vals[i]}^(-1) = ${invs[i]}  (${vals[i]} * ${invs[i]} = ${product} mod ${mod})`);
}

console.log("\nThis is used in elliptic curve point addition for batch verification!");
