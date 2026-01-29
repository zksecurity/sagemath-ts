/**
 * # Modular Inverses
 *
 * The modular inverse is one of the most important concepts in cryptography.
 * It allows us to "divide" in modular arithmetic, which is essential for
 * RSA decryption, elliptic curve operations, and many other algorithms.
 *
 * ## Mathematical Background
 *
 * The modular inverse of a modulo n is an integer a^(-1) such that:
 *   a * a^(-1) = 1 (mod n)
 *
 * The inverse exists if and only if gcd(a, n) = 1 (a and n are coprime).
 *
 * There are several ways to compute modular inverses:
 *   1. Extended Euclidean Algorithm: Most common, O(log n)
 *   2. Fermat's Little Theorem: a^(-1) = a^(p-2) mod p (for prime p)
 *   3. Euler's Theorem: a^(-1) = a^(phi(n)-1) mod n
 *
 * ## Why This Matters for Cryptography
 *
 *   - RSA: Private exponent d = e^(-1) mod phi(n)
 *   - ECDSA: Signature computation requires modular inverse
 *   - Division in finite fields (essential for ECC)
 *   - Lagrange interpolation in secret sharing
 */

import { gcd, xgcd, inverse_mod, euler_phi, power_mod, Mod, Zmod } from 'sagemath-ts';

// =============================================================================
// Example 1: Basic Modular Inverse
// =============================================================================
console.log("=== Example 1: Basic Modular Inverse ===\n");

/**
 * The modular inverse of a mod n satisfies: a * a^(-1) = 1 (mod n)
 */
const modulus = 11n;
console.log(`Finding inverses modulo ${modulus}:\n`);

for (let a = 1n; a < modulus; a++) {
  const inv = inverse_mod(a, modulus);
  const product = (a * inv) % modulus;
  console.log(`  ${a}^(-1) = ${inv.toString().padStart(2)} mod ${modulus}  (${a} * ${inv} = ${a * inv} = ${product} mod ${modulus})`);
}

console.log();

// =============================================================================
// Example 2: When Inverses Don't Exist
// =============================================================================
console.log("=== Example 2: When Inverses Don't Exist ===\n");

/**
 * The inverse of a mod n exists only when gcd(a, n) = 1.
 * Elements without inverses are called zero divisors.
 */
const n = 12n;
console.log(`In Z/${n}Z, not all elements have inverses:\n`);

for (let a = 1n; a < n; a++) {
  const g = gcd(a, n);
  if (g === 1n) {
    const inv = inverse_mod(a, n);
    console.log(`  ${a.toString().padStart(2)} has inverse ${inv.toString().padStart(2)} (gcd = ${g})`);
  } else {
    console.log(`  ${a.toString().padStart(2)} has NO inverse (gcd = ${g})`);
  }
}

console.log(`\nOnly ${euler_phi(n)} elements are invertible (those coprime to ${n}).`);

console.log();

// =============================================================================
// Example 3: Computing Inverse via Extended Euclidean Algorithm
// =============================================================================
console.log("=== Example 3: Extended Euclidean Algorithm Method ===\n");

/**
 * If gcd(a, n) = 1, then xgcd(a, n) returns (1, s, t) where:
 *   1 = s*a + t*n
 * So: s*a = 1 - t*n = 1 (mod n)
 * Therefore: a^(-1) = s (mod n)
 */
function inverseViaXgcd(a: bigint, n: bigint): bigint | null {
  const [g, s, _t] = xgcd(a, n);

  if (g !== 1n) {
    return null;  // No inverse
  }

  // Ensure positive result in [0, n)
  return ((s % n) + n) % n;
}

console.log("Computing inverses via extended Euclidean algorithm:\n");

const examples: [bigint, bigint][] = [
  [3n, 7n],
  [17n, 31n],
  [12345n, 67891n],
];

for (const [a, mod] of examples) {
  const [g, s, t] = xgcd(a, mod);
  console.log(`  xgcd(${a}, ${mod}) = [${g}, ${s}, ${t}]`);

  if (g === 1n) {
    const inv = ((s % mod) + mod) % mod;
    console.log(`    Since 1 = ${s}*${a} + ${t}*${mod}`);
    console.log(`    We have ${a}^(-1) = ${inv} mod ${mod}`);
    console.log(`    Verify: ${a} * ${inv} = ${a * inv} = ${(a * inv) % mod} mod ${mod}\n`);
  }
}

console.log();

// =============================================================================
// Example 4: Fermat's Little Theorem Method
// =============================================================================
console.log("=== Example 4: Fermat's Little Theorem Method ===\n");

/**
 * For prime p and a not divisible by p:
 *   a^(p-1) = 1 (mod p)          [Fermat's Little Theorem]
 *
 * Therefore:
 *   a^(p-2) * a = a^(p-1) = 1 (mod p)
 *
 * So: a^(-1) = a^(p-2) (mod p)
 */
function inverseViaFermat(a: bigint, p: bigint): bigint {
  // Assumes p is prime!
  return power_mod(a, p - 2n, p);
}

const prime = 13n;
console.log(`Using Fermat's Little Theorem for prime p = ${prime}:\n`);
console.log(`  a^(-1) = a^(p-2) = a^${prime - 2n} (mod ${prime})\n`);

for (let a = 1n; a < prime; a++) {
  const inv_xgcd = inverse_mod(a, prime);
  const inv_fermat = inverseViaFermat(a, prime);
  console.log(`  ${a}^(-1): xgcd=${inv_xgcd}, Fermat=${inv_fermat}, match=${inv_xgcd === inv_fermat}`);
}

console.log();

// =============================================================================
// Example 5: Euler's Theorem Method
// =============================================================================
console.log("=== Example 5: Euler's Theorem Method ===\n");

/**
 * Euler's Theorem (generalizes Fermat's Little Theorem):
 *   a^(phi(n)) = 1 (mod n) when gcd(a, n) = 1
 *
 * Therefore: a^(-1) = a^(phi(n)-1) (mod n)
 */
function inverseViaEuler(a: bigint, n: bigint): bigint | null {
  if (gcd(a, n) !== 1n) return null;
  const phi = euler_phi(n);
  return power_mod(a, phi - 1n, n);
}

const compositeN = 15n;
const phiN = euler_phi(compositeN);
console.log(`Using Euler's Theorem for n = ${compositeN} (phi(${compositeN}) = ${phiN}):\n`);

for (let a = 1n; a < compositeN; a++) {
  if (gcd(a, compositeN) === 1n) {
    const inv_xgcd = inverse_mod(a, compositeN);
    const inv_euler = inverseViaEuler(a, compositeN);
    console.log(`  ${a.toString().padStart(2)}^(-1): xgcd=${inv_xgcd}, Euler=${inv_euler}, match=${inv_xgcd === inv_euler}`);
  }
}

console.log();

// =============================================================================
// Example 6: Using IntegerMod Class
// =============================================================================
console.log("=== Example 6: Using the IntegerMod Class ===\n");

/**
 * The IntegerMod class provides convenient methods for modular arithmetic
 * including .inv() for inverse and .div() for division.
 */
const Z23 = Zmod(23);

console.log("Using IntegerMod methods in Z/23Z:\n");

const elem = Z23.__call__(7);
console.log(`  Element: ${elem.value}`);
console.log(`  Inverse: ${elem.inv().value}`);
console.log(`  Verify: ${elem.value} * ${elem.inv().value} = ${elem.mul(elem.inv()).value} (mod 23)`);

// Division using inverse
const num = Z23.__call__(15);
const denom = Z23.__call__(4);
const quotient = num.div(denom);

console.log(`\n  Division: 15 / 4 mod 23`);
console.log(`    4^(-1) = ${denom.inv().value}`);
console.log(`    15 * ${denom.inv().value} = ${quotient.value} mod 23`);
console.log(`    Verify: 4 * ${quotient.value} = ${denom.mul(quotient).value} mod 23`);

console.log();

// =============================================================================
// Example 7: Inverse in Cryptographic Key Generation
// =============================================================================
console.log("=== Example 7: RSA Key Generation ===\n");

/**
 * In RSA, the private exponent d is computed as:
 *   d = e^(-1) mod phi(n)
 *
 * This requires computing a modular inverse.
 */
const p = 61n;
const q = 53n;
const rsaN = p * q;
const phi = (p - 1n) * (q - 1n);
const e = 17n;

console.log("RSA Key Generation:");
console.log(`  p = ${p}, q = ${q}`);
console.log(`  n = ${rsaN}`);
console.log(`  phi(n) = ${phi}`);
console.log(`  e = ${e} (public exponent)`);

// Verify e and phi(n) are coprime
const g = gcd(e, phi);
console.log(`  gcd(e, phi(n)) = ${g}`);

if (g === 1n) {
  const d = inverse_mod(e, phi);
  console.log(`  d = e^(-1) mod phi(n) = ${d} (private exponent)`);
  console.log(`\n  Verification: e * d mod phi(n) = ${(e * d) % phi}`);
}

console.log();

// =============================================================================
// Example 8: Division in Modular Arithmetic
// =============================================================================
console.log("=== Example 8: Division via Multiplication by Inverse ===\n");

/**
 * In modular arithmetic, a/b = a * b^(-1) (when b is invertible).
 * This is how "division" works in finite fields.
 */
function modDiv(a: bigint, b: bigint, n: bigint): bigint | null {
  if (gcd(b, n) !== 1n) return null;
  return (a * inverse_mod(b, n)) % n;
}

const mod = 17n;
console.log(`Division in Z/${mod}Z:\n`);

const divExamples: [bigint, bigint][] = [
  [10n, 3n],
  [5n, 7n],
  [1n, 11n],
];

for (const [a, b] of divExamples) {
  const result = modDiv(a, b, mod);
  console.log(`  ${a} / ${b} mod ${mod}:`);
  console.log(`    ${b}^(-1) = ${inverse_mod(b, mod)}`);
  console.log(`    ${a} * ${inverse_mod(b, mod)} = ${result} mod ${mod}`);
  console.log(`    Verify: ${result} * ${b} = ${(result! * b) % mod} mod ${mod}\n`);
}

console.log();

// =============================================================================
// Example 9: Batch Inverse Computation
// =============================================================================
console.log("=== Example 9: Montgomery's Trick for Batch Inverses ===\n");

/**
 * Computing n inverses normally requires n extended GCD operations.
 * Montgomery's trick computes n inverses with just 1 extended GCD + 3(n-1) multiplications!
 *
 * Key insight: If we know (a1*a2*...*an)^(-1), we can extract individual inverses.
 */
function batchInverse(values: bigint[], n: bigint): bigint[] {
  const len = values.length;
  if (len === 0) return [];

  // Compute cumulative products: c[i] = a1 * a2 * ... * a[i]
  const cumulative: bigint[] = new Array(len);
  cumulative[0] = values[0]! % n;
  for (let i = 1; i < len; i++) {
    cumulative[i] = (cumulative[i-1]! * values[i]!) % n;
  }

  // ONE inverse computation for the total product
  let totalInv = inverse_mod(cumulative[len-1]!, n);

  // Extract individual inverses working backwards
  const inverses: bigint[] = new Array(len);
  for (let i = len - 1; i > 0; i--) {
    inverses[i] = (totalInv * cumulative[i-1]!) % n;
    totalInv = (totalInv * values[i]!) % n;
  }
  inverses[0] = totalInv;

  return inverses;
}

const batchMod = 101n;
const batchValues = [3n, 7n, 13n, 23n, 37n];

console.log(`Batch inverse computation mod ${batchMod}:`);
console.log(`  Values: [${batchValues.join(', ')}]\n`);

const batchResult = batchInverse(batchValues, batchMod);

console.log("Results and verification:");
for (let i = 0; i < batchValues.length; i++) {
  const v = batchValues[i]!;
  const inv = batchResult[i]!;
  const product = (v * inv) % batchMod;
  console.log(`  ${v}^(-1) = ${inv}, verify: ${v} * ${inv} = ${product} mod ${batchMod}`);
}

console.log("\nMontgomery's trick is used in elliptic curve batch verification!");

console.log();

// =============================================================================
// Exercise: Implement Modular Inverse
// =============================================================================
console.log("=== Exercise: Implement Modular Inverse ===\n");

/**
 * EXERCISE: Implement modular inverse using the extended Euclidean algorithm.
 */

// TODO: Implement this function
function myInverseMod(a: bigint, n: bigint): bigint | null {
  // Extended Euclidean Algorithm
  let oldR = a, r = n;
  let oldS = 1n, s = 0n;

  while (r !== 0n) {
    const quotient = oldR / r;
    [oldR, r] = [r, oldR - quotient * r];
    [oldS, s] = [s, oldS - quotient * s];
  }

  // If GCD is not 1, no inverse exists
  if (oldR !== 1n && oldR !== -1n) {
    return null;
  }

  // Handle negative GCD (can happen due to signs)
  if (oldR === -1n) {
    oldS = -oldS;
  }

  // Ensure positive result
  return ((oldS % n) + n) % n;
}

console.log("Testing your modular inverse implementation:");
const testCases: [bigint, bigint][] = [
  [3n, 11n],
  [17n, 23n],
  [6n, 15n],   // No inverse (gcd = 3)
];

for (const [a, mod] of testCases) {
  const yours = myInverseMod(a, mod);
  let library: bigint | null;
  try {
    library = inverse_mod(a, mod);
  } catch {
    library = null;
  }
  console.log(`  ${a}^(-1) mod ${mod}: yours=${yours}, library=${library}, match=${yours === library}`);
}

console.log();

// =============================================================================
// Cryptographic Application: ECDSA Signature
// =============================================================================
console.log("=== Cryptographic Application: ECDSA Signature ===\n");

/**
 * ECDSA signature computation requires a modular inverse:
 *   s = k^(-1) * (hash + private_key * r) mod n
 *
 * Where n is the order of the elliptic curve group.
 */
console.log("ECDSA signature computation (simplified):\n");

// Simplified example with small numbers
const curveOrder = 97n;  // Order of the curve (would be ~256 bits in practice)
const privateKey = 37n;
const k = 42n;  // Random nonce (MUST be unique per signature!)
const hash = 73n;  // Message hash
const r = 55n;  // x-coordinate of k*G

console.log(`  Curve order n = ${curveOrder}`);
console.log(`  Private key d = ${privateKey}`);
console.log(`  Nonce k = ${k}`);
console.log(`  Message hash h = ${hash}`);
console.log(`  r (from k*G) = ${r}`);

// Compute s = k^(-1) * (hash + d * r) mod n
const kInv = inverse_mod(k, curveOrder);
const innerSum = (hash + privateKey * r) % curveOrder;
const s = (kInv * innerSum) % curveOrder;

console.log(`\n  k^(-1) = ${kInv}`);
console.log(`  hash + d*r = ${innerSum} mod ${curveOrder}`);
console.log(`  s = k^(-1) * (hash + d*r) = ${s} mod ${curveOrder}`);
console.log(`\n  Signature: (r=${r}, s=${s})`);

console.log("\nThe modular inverse is essential for ECDSA!");
