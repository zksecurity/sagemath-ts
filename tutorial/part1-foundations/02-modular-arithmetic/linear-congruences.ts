/**
 * # Solving Linear Congruences: ax = b (mod n)
 *
 * Linear congruences are equations of the form ax = b (mod n).
 * Solving them is a fundamental skill in number theory and cryptography.
 *
 * ## Mathematical Background
 *
 * The linear congruence ax = b (mod n) has a solution if and only if
 * gcd(a, n) divides b.
 *
 * If gcd(a, n) = g and g | b, then:
 *   - There are exactly g solutions modulo n
 *   - Solutions differ by n/g
 *   - General solution: x = x0 + k*(n/g) for k = 0, 1, ..., g-1
 *
 * Special case: If gcd(a, n) = 1, there is exactly one solution:
 *   x = a^(-1) * b (mod n)
 *
 * ## Why This Matters for Cryptography
 *
 *   - RSA key generation: Finding d such that ed = 1 (mod phi(n))
 *   - Solving systems of congruences (Chinese Remainder Theorem)
 *   - Lattice-based cryptography: Finding short vectors
 *   - Discrete log problems can be viewed as linear congruences
 */

import { gcd, xgcd, inverse_mod, crt, CRT_list, Mod, Zmod } from 'sagemath-ts';

// =============================================================================
// Example 1: Simple Linear Congruence with Unique Solution
// =============================================================================
console.log("=== Example 1: Linear Congruence with Unique Solution ===\n");

/**
 * When gcd(a, n) = 1, the equation ax = b (mod n) has exactly one solution.
 * Solution: x = a^(-1) * b (mod n)
 */
function solveLinearCongruence(a: bigint, b: bigint, n: bigint): bigint[] | null {
  const g = gcd(a, n);

  // Check if solution exists
  if (b % g !== 0n) {
    return null;  // No solution
  }

  // Reduce to the case where gcd(a', n') = 1
  const aReduced = a / g;
  const bReduced = b / g;
  const nReduced = n / g;

  // Find one particular solution
  const aInv = inverse_mod(aReduced, nReduced);
  let x0 = (aInv * bReduced) % nReduced;
  if (x0 < 0n) x0 += nReduced;

  // Generate all g solutions
  const solutions: bigint[] = [];
  for (let k = 0n; k < g; k++) {
    solutions.push((x0 + k * nReduced) % n);
  }

  return solutions.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
}

// Example with unique solution
const a1 = 7n, b1 = 3n, n1 = 11n;
console.log(`Solving ${a1}x = ${b1} (mod ${n1}):`);
console.log(`  gcd(${a1}, ${n1}) = ${gcd(a1, n1)}`);

const inv = inverse_mod(a1, n1);
const x1 = (inv * b1) % n1;
console.log(`  ${a1}^(-1) = ${inv} mod ${n1}`);
console.log(`  x = ${inv} * ${b1} = ${x1} mod ${n1}`);
console.log(`  Verify: ${a1} * ${x1} = ${a1 * x1} = ${(a1 * x1) % n1} mod ${n1}`);

console.log();

// =============================================================================
// Example 2: Linear Congruence with Multiple Solutions
// =============================================================================
console.log("=== Example 2: Linear Congruence with Multiple Solutions ===\n");

/**
 * When gcd(a, n) = g > 1 and g divides b, there are g solutions.
 */
const a2 = 6n, b2 = 9n, n2 = 15n;
const g2 = gcd(a2, n2);

console.log(`Solving ${a2}x = ${b2} (mod ${n2}):`);
console.log(`  gcd(${a2}, ${n2}) = ${g2}`);
console.log(`  Does ${g2} divide ${b2}? ${b2 % g2 === 0n}`);

const solutions2 = solveLinearCongruence(a2, b2, n2);
if (solutions2) {
  console.log(`  Solutions: x = ${solutions2.join(', ')} (mod ${n2})`);
  console.log(`  Number of solutions: ${solutions2.length}`);

  console.log("\n  Verification:");
  for (const x of solutions2) {
    const product = (a2 * x) % n2;
    console.log(`    ${a2} * ${x} = ${a2 * x} = ${product} mod ${n2}`);
  }
}

console.log();

// =============================================================================
// Example 3: Linear Congruence with No Solution
// =============================================================================
console.log("=== Example 3: Linear Congruence with No Solution ===\n");

/**
 * When gcd(a, n) does not divide b, there is no solution.
 */
const a3 = 6n, b3 = 7n, n3 = 9n;
const g3 = gcd(a3, n3);

console.log(`Solving ${a3}x = ${b3} (mod ${n3}):`);
console.log(`  gcd(${a3}, ${n3}) = ${g3}`);
console.log(`  Does ${g3} divide ${b3}? ${b3 % g3 === 0n}`);

const solutions3 = solveLinearCongruence(a3, b3, n3);
console.log(`  Solutions: ${solutions3 === null ? 'None' : solutions3.join(', ')}`);

// Verify by trying all values
console.log("\n  Trying all residues:");
for (let x = 0n; x < n3; x++) {
  const product = (a3 * x) % n3;
  console.log(`    ${a3} * ${x} = ${product} mod ${n3}${product === b3 ? ' <-- found!' : ''}`);
}
console.log("  (No solution exists)");

console.log();

// =============================================================================
// Example 4: Systematic Algorithm
// =============================================================================
console.log("=== Example 4: Step-by-Step Algorithm ===\n");

/**
 * Algorithm to solve ax = b (mod n):
 * 1. Compute g = gcd(a, n)
 * 2. If g does not divide b, no solution
 * 3. Divide through by g: (a/g)x = (b/g) (mod n/g)
 * 4. Now gcd(a/g, n/g) = 1, so use inverse
 * 5. x0 = (a/g)^(-1) * (b/g) mod (n/g)
 * 6. General solution: x = x0 + k*(n/g) for k = 0, 1, ..., g-1
 */
function solveLinearCongruenceVerbose(a: bigint, b: bigint, n: bigint): void {
  console.log(`Solving ${a}x = ${b} (mod ${n}):\n`);

  // Step 1
  const g = gcd(a, n);
  console.log(`  Step 1: g = gcd(${a}, ${n}) = ${g}`);

  // Step 2
  if (b % g !== 0n) {
    console.log(`  Step 2: ${g} does not divide ${b} -> NO SOLUTION`);
    return;
  }
  console.log(`  Step 2: ${g} divides ${b} -> solution exists`);

  // Step 3
  const aPrime = a / g;
  const bPrime = b / g;
  const nPrime = n / g;
  console.log(`  Step 3: Divide by g: ${aPrime}x = ${bPrime} (mod ${nPrime})`);
  console.log(`          Now gcd(${aPrime}, ${nPrime}) = ${gcd(aPrime, nPrime)}`);

  // Step 4
  const aInv = inverse_mod(aPrime, nPrime);
  console.log(`  Step 4: ${aPrime}^(-1) = ${aInv} mod ${nPrime}`);

  // Step 5
  const x0 = (aInv * bPrime) % nPrime;
  console.log(`  Step 5: x0 = ${aInv} * ${bPrime} = ${x0} mod ${nPrime}`);

  // Step 6
  const solutions: bigint[] = [];
  for (let k = 0n; k < g; k++) {
    solutions.push((x0 + k * nPrime) % n);
  }
  console.log(`  Step 6: ${g} solutions: x = ${x0} + k*${nPrime} for k = 0, ..., ${g - 1n}`);
  console.log(`          x in {${solutions.sort((a, b) => (a < b ? -1 : 1)).join(', ')}}`);
}

solveLinearCongruenceVerbose(12n, 18n, 30n);

console.log();

// =============================================================================
// Example 5: Chinese Remainder Theorem
// =============================================================================
console.log("=== Example 5: Chinese Remainder Theorem ===\n");

/**
 * The Chinese Remainder Theorem (CRT) solves systems of linear congruences:
 *   x = a1 (mod n1)
 *   x = a2 (mod n2)
 *   ...
 *
 * When n1, n2, ... are pairwise coprime, there is a unique solution mod (n1*n2*...).
 */
console.log("System of congruences:");
console.log("  x = 2 (mod 3)");
console.log("  x = 3 (mod 5)");
console.log("  x = 2 (mod 7)\n");

const residues = [2n, 3n, 2n];
const moduli = [3n, 5n, 7n];

const crtSolution = CRT_list(residues, moduli);
const productMod = moduli.reduce((a, b) => a * b, 1n);

console.log(`Using CRT_list: x = ${crtSolution}`);
console.log(`Product of moduli: ${productMod}`);

console.log("\nVerification:");
for (let i = 0; i < residues.length; i++) {
  const remainder = crtSolution % moduli[i]!;
  console.log(`  ${crtSolution} mod ${moduli[i]} = ${remainder} (expected ${residues[i]})`);
}

console.log();

// =============================================================================
// Example 6: CRT Two-Congruence Case
// =============================================================================
console.log("=== Example 6: Two-Congruence CRT ===\n");

/**
 * For two congruences, the crt() function is convenient:
 *   x = a (mod m)
 *   x = b (mod n)
 */
console.log("System:");
console.log("  x = 5 (mod 7)");
console.log("  x = 3 (mod 11)\n");

const x_crt = crt(5n, 3n, 7n, 11n);
console.log(`CRT solution: x = ${x_crt} (mod ${7n * 11n})`);

console.log(`\nVerification:`);
console.log(`  ${x_crt} mod 7 = ${x_crt % 7n} (expected 5)`);
console.log(`  ${x_crt} mod 11 = ${x_crt % 11n} (expected 3)`);

console.log();

// =============================================================================
// Example 7: RSA with CRT Optimization
// =============================================================================
console.log("=== Example 7: RSA-CRT Optimization ===\n");

/**
 * RSA decryption can be sped up using CRT.
 * Instead of computing c^d mod n directly, we compute:
 *   mp = c^(d mod p-1) mod p
 *   mq = c^(d mod q-1) mod q
 * Then combine using CRT.
 */
const p = 61n;
const q = 53n;
const n = p * q;
const phi = (p - 1n) * (q - 1n);
const e = 17n;
const d = inverse_mod(e, phi);

console.log("RSA parameters:");
console.log(`  p = ${p}, q = ${q}, n = ${n}`);
console.log(`  e = ${e}, d = ${d}\n`);

// Encrypt a message
const message = 42n;
const ciphertext = Mod(message, n).pow(e).value;
console.log(`  Message: ${message}`);
console.log(`  Ciphertext: ${ciphertext}\n`);

// Standard decryption
const start1 = performance.now();
const decrypted1 = Mod(ciphertext, n).pow(d).value;
const time1 = performance.now() - start1;

// CRT decryption
const start2 = performance.now();
const dp = d % (p - 1n);
const dq = d % (q - 1n);
const mp = Mod(ciphertext, p).pow(dp).value;
const mq = Mod(ciphertext, q).pow(dq).value;
const decrypted2 = crt(mp, mq, p, q);
const time2 = performance.now() - start2;

console.log("Standard decryption:");
console.log(`  c^d mod n = ${decrypted1}`);

console.log("\nCRT decryption:");
console.log(`  dp = d mod (p-1) = ${dp}`);
console.log(`  dq = d mod (q-1) = ${dq}`);
console.log(`  mp = c^dp mod p = ${mp}`);
console.log(`  mq = c^dq mod q = ${mq}`);
console.log(`  CRT(mp, mq) = ${decrypted2}`);

console.log(`\nBoth methods give: ${decrypted1} (original: ${message})`);

console.log();

// =============================================================================
// Example 8: Solving Polynomial Congruences
// =============================================================================
console.log("=== Example 8: Polynomial Congruences ===\n");

/**
 * Linear congruences generalize to polynomial congruences.
 * f(x) = 0 (mod n)
 *
 * For prime n, a polynomial of degree d has at most d roots.
 */
function findRoots(coeffs: bigint[], n: bigint): bigint[] {
  // coeffs[i] is coefficient of x^i
  const roots: bigint[] = [];

  for (let x = 0n; x < n; x++) {
    let sum = 0n;
    let xPower = 1n;
    for (const c of coeffs) {
      sum = (sum + c * xPower) % n;
      xPower = (xPower * x) % n;
    }
    if (sum === 0n) {
      roots.push(x);
    }
  }

  return roots;
}

// x^2 - 1 = 0 (mod 7) should have 2 roots: 1 and 6
console.log("Solving x^2 - 1 = 0 (mod 7):");
const roots1 = findRoots([-1n, 0n, 1n], 7n);  // -1 + 0*x + 1*x^2
console.log(`  Roots: ${roots1.join(', ')}`);
console.log(`  Verification: 1^2 = 1, 6^2 = 36 = 1 (mod 7)\n`);

// x^2 - 1 = 0 (mod 8) should have 4 roots: 1, 3, 5, 7
console.log("Solving x^2 - 1 = 0 (mod 8):");
const roots2 = findRoots([-1n, 0n, 1n], 8n);
console.log(`  Roots: ${roots2.join(', ')}`);
console.log(`  (More roots than degree! This happens when n is not prime.)`);

console.log();

// =============================================================================
// Example 9: Discrete Logarithm as Linear Congruence
// =============================================================================
console.log("=== Example 9: Discrete Log Connection ===\n");

/**
 * The discrete logarithm problem: find x such that g^x = h (mod p)
 * This can be viewed as solving: x * log(g) = log(h) in the exponent group.
 *
 * For small groups, we can solve by exhaustive search.
 */
function discreteLog(g: bigint, h: bigint, p: bigint): bigint | null {
  let power = 1n;
  for (let x = 0n; x < p; x++) {
    if (power === h) {
      return x;
    }
    power = (power * g) % p;
  }
  return null;
}

const prime = 17n;
const generator = 3n;  // 3 is a primitive root mod 17
const target = 7n;

console.log(`Solving ${generator}^x = ${target} (mod ${prime}):`);

const dlog = discreteLog(generator, target, prime);
if (dlog !== null) {
  console.log(`  x = ${dlog}`);
  console.log(`  Verify: ${generator}^${dlog} = ${Mod(generator, prime).pow(dlog).value} mod ${prime}`);
}

console.log("\nThe discrete log problem is believed to be hard for large primes.");
console.log("This hardness is the foundation of Diffie-Hellman key exchange.");

console.log();

// =============================================================================
// Exercise: Solve a System of Congruences
// =============================================================================
console.log("=== Exercise: Solve a System Manually ===\n");

/**
 * EXERCISE: Solve this system step by step:
 *   x = 3 (mod 4)
 *   x = 1 (mod 5)
 *   x = 4 (mod 7)
 */
console.log("System to solve:");
console.log("  x = 3 (mod 4)");
console.log("  x = 1 (mod 5)");
console.log("  x = 4 (mod 7)\n");

// Step 1: Combine first two
console.log("Step 1: Combine first two congruences");
const x12 = crt(3n, 1n, 4n, 5n);
const m12 = 4n * 5n;
console.log(`  x = 3 (mod 4) and x = 1 (mod 5)`);
console.log(`  -> x = ${x12} (mod ${m12})`);

// Step 2: Combine with third
console.log("\nStep 2: Combine with third congruence");
const x123 = crt(x12, 4n, m12, 7n);
const m123 = m12 * 7n;
console.log(`  x = ${x12} (mod ${m12}) and x = 4 (mod 7)`);
console.log(`  -> x = ${x123} (mod ${m123})`);

// Verify
console.log("\nVerification:");
console.log(`  ${x123} mod 4 = ${x123 % 4n} (expected 3)`);
console.log(`  ${x123} mod 5 = ${x123 % 5n} (expected 1)`);
console.log(`  ${x123} mod 7 = ${x123 % 7n} (expected 4)`);

console.log();

// =============================================================================
// Cryptographic Application: Secret Sharing
// =============================================================================
console.log("=== Cryptographic Application: Shamir's Secret Sharing ===\n");

/**
 * Shamir's Secret Sharing uses polynomial interpolation over a finite field.
 * To reconstruct, we solve a system of linear equations (or use Lagrange interpolation).
 *
 * Secret S is hidden as f(0) where f is a degree k-1 polynomial.
 * We need k points (x_i, y_i = f(x_i)) to recover S.
 */
const secretPrime = 17n;
const secret = 7n;

// Create polynomial f(x) = secret + coef1*x + coef2*x^2 (for k=3 threshold)
const coef1 = 3n, coef2 = 5n;

function f(x: bigint): bigint {
  return ((secret + coef1 * x + coef2 * x * x) % secretPrime + secretPrime) % secretPrime;
}

console.log("Secret sharing example (k=3 threshold):");
console.log(`  Prime: ${secretPrime}`);
console.log(`  Secret: ${secret}`);
console.log(`  f(x) = ${secret} + ${coef1}x + ${coef2}x^2 (mod ${secretPrime})\n`);

// Generate shares
const shares: [bigint, bigint][] = [];
for (let i = 1n; i <= 5n; i++) {
  shares.push([i, f(i)]);
}

console.log("Shares distributed:");
for (const [x, y] of shares) {
  console.log(`  Share ${x}: f(${x}) = ${y}`);
}

// Lagrange interpolation to recover f(0) = secret
function lagrangeInterpolate(points: [bigint, bigint][], x: bigint, p: bigint): bigint {
  let result = 0n;

  for (let i = 0; i < points.length; i++) {
    const [xi, yi] = points[i]!;

    let numerator = yi;
    let denominator = 1n;

    for (let j = 0; j < points.length; j++) {
      if (i !== j) {
        const [xj] = points[j]!;
        numerator = (numerator * ((x - xj + p) % p)) % p;
        denominator = (denominator * ((xi - xj + p) % p)) % p;
      }
    }

    result = (result + numerator * inverse_mod(denominator, p)) % p;
  }

  return result;
}

// Use any 3 shares to recover
const selectedShares = [shares[0]!, shares[2]!, shares[4]!];
console.log(`\nUsing 3 shares to recover: (${selectedShares.map(s => s[0]).join(', ')})`);

const recovered = lagrangeInterpolate(selectedShares, 0n, secretPrime);
console.log(`  f(0) = ${recovered} (original secret: ${secret})`);
