/**
 * # Pohlig-Hellman Algorithm
 *
 * The Pohlig-Hellman algorithm solves the discrete logarithm problem
 * efficiently when the group order has only small prime factors.
 *
 * ## Key Insight
 *
 * If the group order n = p1^e1 * p2^e2 * ... * pk^ek, then:
 * 1. Solve DLP in subgroups of order pi^ei
 * 2. Combine solutions using Chinese Remainder Theorem
 *
 * ## Complexity
 *
 * Time: O(sum of ei * (log n + sqrt(pi))) for n = prod(pi^ei)
 *
 * If all pi are small (n is "smooth"), this is much faster than sqrt(n).
 *
 * ## Security Implication
 *
 * Groups with smooth order are WEAK for cryptography!
 * This is why we use prime-order subgroups or safe primes.
 */

import { power_mod, gcd, is_prime, isqrt, inverse_mod, factor, formatFactorization, crt } from 'sagemath-ts';

console.log("=== Pohlig-Hellman Algorithm ===\n");

// Baby-step giant-step for prime-order groups
function bsgs(g: bigint, h: bigint, p: bigint, order: bigint): bigint | null {
  const m = isqrt(order) + 1n;

  const table = new Map<bigint, bigint>();
  let gj = 1n;
  for (let j = 0n; j < m; j++) {
    table.set(gj, j);
    gj = (gj * g) % p;
  }

  const gmInv = inverse_mod(power_mod(g, m, p), p);
  let gamma = h;

  for (let i = 0n; i < m; i++) {
    const j = table.get(gamma);
    if (j !== undefined) {
      return i * m + j;
    }
    gamma = (gamma * gmInv) % p;
  }

  return null;
}

// Solve DLP in prime power subgroup using Pohlig-Hellman reduction
function dlpPrimePower(g: bigint, h: bigint, p: bigint, order: bigint, prime: bigint, exp: bigint): bigint {
  // Solve g^x = h where the order is prime^exp
  // Uses the "digit extraction" technique

  const n = order; // Full group order (p - 1)
  const q = prime;
  const e = Number(exp);

  let x = 0n;

  // Compute g_q = g^(n/q) which has order q
  const gq = power_mod(g, n / q, p);

  // Extract digits x = x_0 + x_1*q + x_2*q^2 + ... + x_(e-1)*q^(e-1)
  let gamma = h;

  for (let i = 0; i < e; i++) {
    // Compute gamma^(n/q^(i+1))
    const expFactor = n / (q ** BigInt(i + 1));
    const hi = power_mod(gamma, expFactor, p);

    // Solve gq^xi = hi using BSGS
    let xi = 0n;
    if (hi !== 1n) {
      const result = bsgs(gq, hi, p, q);
      if (result !== null) {
        xi = result;
      }
    }

    x += xi * (q ** BigInt(i));

    // Update gamma = gamma * g^(-xi * q^i)
    if (xi !== 0n) {
      const gInv = inverse_mod(power_mod(g, xi * (q ** BigInt(i)), p), p);
      gamma = (gamma * gInv) % p;
    }
  }

  return x % (q ** exp);
}

// Full Pohlig-Hellman algorithm
function pohligHellman(g: bigint, h: bigint, p: bigint): bigint | null {
  const order = p - 1n;
  const factors = factor(order).filter(([prime]) => prime > 0n);

  console.log(`  Group order: ${order} = ${formatFactorization(factors)}`);

  const residues: bigint[] = [];
  const moduli: bigint[] = [];

  for (const [prime, exp] of factors) {
    const primeExp = prime ** exp;

    console.log(`\n  Solving modulo ${prime}^${exp} = ${primeExp}:`);

    const x_mod = dlpPrimePower(g, h, p, order, prime, exp);

    console.log(`    x = ${x_mod} (mod ${primeExp})`);

    residues.push(x_mod);
    moduli.push(primeExp);
  }

  // Combine using CRT
  console.log(`\n  Combining with CRT: x = ${residues.join(', ')} mod ${moduli.join(', ')}`);

  let result = residues[0]!;
  let modulus = moduli[0]!;

  for (let i = 1; i < residues.length; i++) {
    result = crt(result, residues[i]!, modulus, moduli[i]!);
    modulus = modulus * moduli[i]! / gcd(modulus, moduli[i]!);
  }

  return result;
}

// Example 1: Why smooth orders are weak
console.log("Example 1: Smooth Order Weakness");
console.log("-".repeat(50));

console.log(`
  Consider a group of order n. The DLP difficulty depends on:

  Order n            | Largest prime factor | BSGS time  | Pohlig-Hellman
  -------------------|---------------------|------------|---------------
  1000 = 2^3 * 5^3   | 5                   | sqrt(1000) | sqrt(5)
  1024 = 2^10        | 2                   | sqrt(1024) | sqrt(2)
  1009 (prime)       | 1009                | sqrt(1009) | sqrt(1009)

  If the largest prime factor is small, Pohlig-Hellman wins!
`);

// Example 2: Step-by-step Pohlig-Hellman
console.log("\nExample 2: Pohlig-Hellman Step by Step");
console.log("-".repeat(50));

// Choose p where p-1 has small factors
const p = 433n; // p - 1 = 432 = 2^4 * 3^3 = 16 * 27
const g = 5n;   // A generator
const target_x = 47n;
const h = power_mod(g, target_x, p);

console.log(`  p = ${p}, g = ${g}, h = ${h}`);
console.log(`  Target: find x such that ${g}^x = ${h} mod ${p}`);
console.log(`  (We secretly know x = ${target_x} for verification)\n`);

const result = pohligHellman(g, h, p);

console.log(`\n  Final answer: x = ${result}`);
console.log(`  Verification: ${g}^${result} = ${power_mod(g, result!, p)} (should be ${h})`);

// Example 3: Complexity comparison
console.log("\nExample 3: Complexity Comparison");
console.log("-".repeat(50));

console.log(`
  For order n = p1^e1 * p2^e2 * ... * pk^ek:

  Baby-step Giant-step: O(sqrt(n))
  Pohlig-Hellman: O(sum of ei * sqrt(pi))

  Example: n = 2^10 * 3^5 * 5^3 = 31,104,000

  BSGS: sqrt(31,104,000) ~ 5,577 steps

  Pohlig-Hellman:
  - Prime 2: 10 * sqrt(2) ~ 14 steps
  - Prime 3: 5 * sqrt(3) ~ 9 steps
  - Prime 5: 3 * sqrt(5) ~ 7 steps
  - Total: ~ 30 steps

  Speedup: ~186x for this example!
`);

// Example 4: The digit extraction technique
console.log("Example 4: Digit Extraction Technique");
console.log("-".repeat(50));

console.log(`
  For prime power q^e, we extract x digit by digit in base q.

  x = x_0 + x_1*q + x_2*q^2 + ... + x_(e-1)*q^(e-1)

  where each x_i is in [0, q).

  Key identity: h^(n/q) = g^(x * n/q)

  Since g has order n, g^(n/q) has order q.
  And h^(n/q) = g^(x * n/q) = (g^(n/q))^x

  In the order-q subgroup: h^(n/q) = (g^(n/q))^(x_0)

  So x_0 = DLP in a group of order q (easy!).
  Then adjust and repeat for x_1, x_2, ...
`);

// Visualize digit extraction
console.log("\n  Example: q = 2, e = 4, n = 432:");

const order432 = 432n;
const g432 = 5n;
const x_secret = 11n; // 11 = 1 + 1*2 + 0*4 + 1*8 (binary: 1011)
const h432 = power_mod(g432, x_secret, p);

console.log(`  x = ${x_secret} in binary: ${x_secret.toString(2)}`);
console.log(`  Expected digits (little-endian): 1, 1, 0, 1\n`);

let gamma = h432;
const g2 = power_mod(g432, order432 / 2n, p); // g^(n/2) has order 2

for (let i = 0; i < 4; i++) {
  const expFactor = order432 / (2n ** BigInt(i + 1));
  const hi = power_mod(gamma, expFactor, p);

  const xi = hi === 1n ? 0n : 1n; // In order-2 group, only 0 or 1

  console.log(`  i = ${i}: gamma^(${expFactor}) = ${hi}, x_${i} = ${xi}`);

  if (xi !== 0n) {
    const adj = power_mod(inverse_mod(power_mod(g432, xi * (2n ** BigInt(i)), p), p), 1n, p);
    gamma = (gamma * adj) % p;
  }
}

// Example 5: Chinese Remainder Theorem combination
console.log("\nExample 5: CRT Combination");
console.log("-".repeat(50));

console.log(`
  After solving DLP modulo each prime power:
  - x = r1 (mod q1^e1)
  - x = r2 (mod q2^e2)
  - ...

  We combine using CRT to find x mod (q1^e1 * q2^e2 * ...) = n.

  Example: x = 3 (mod 16), x = 2 (mod 27)

  Need to find x mod 432 = 16 * 27.
`);

const r1 = 3n;
const m1 = 16n;
const r2 = 2n;
const m2 = 27n;

const combined = crt(r1, r2, m1, m2);
console.log(`  CRT(${r1}, ${r2}, ${m1}, ${m2}) = ${combined}`);
console.log(`  Verify: ${combined} mod ${m1} = ${combined % m1} (should be ${r1})`);
console.log(`  Verify: ${combined} mod ${m2} = ${combined % m2} (should be ${r2})`);

// Example 6: Security implications
console.log("\nExample 6: Security Implications");
console.log("-".repeat(50));

console.log(`
  CRITICAL: Groups with smooth order are cryptographically weak!

  For (Z/pZ)*:
  - Order is p - 1
  - If p - 1 has only small factors, DLP is easy

  Defense strategies:

  1. Safe primes: p = 2q + 1 where q is also prime
     - Order = 2q, largest prime factor is q ~ p/2
     - DLP difficulty: sqrt(q) ~ sqrt(p)

  2. Schnorr groups: Work in prime-order subgroup
     - Choose large prime q | (p - 1)
     - Generator has order exactly q
     - DLP difficulty: sqrt(q)

  3. Elliptic curves: Group order is typically almost prime
     - #E(F_p) ~ p with small cofactor
     - No Pohlig-Hellman weakness
`);

// Demonstrate weak vs strong primes
console.log("\n  Comparing weak vs strong primes:");

const weakP = 257n;  // 257 - 1 = 256 = 2^8 (very smooth!)
const strongP = 263n; // 263 - 1 = 262 = 2 * 131 (131 is large prime)

console.log(`\n  Weak prime: p = ${weakP}`);
console.log(`    Order = ${weakP - 1n} = ${formatFactorization(factor(weakP - 1n))}`);
console.log(`    Largest factor: 2, Pohlig-Hellman: O(8 * sqrt(2)) ~ 12 ops`);

console.log(`\n  Strong prime: p = ${strongP}`);
console.log(`    Order = ${strongP - 1n} = ${formatFactorization(factor(strongP - 1n))}`);
console.log(`    Largest factor: 131, Pohlig-Hellman: O(sqrt(131)) ~ 12 ops`);

// Example 7: Real-world group sizes
console.log("\nExample 7: Real-World Group Sizes");
console.log("-".repeat(50));

console.log(`
  Recommended group orders and their structure:

  Protocol    | Order bits | Structure           | Pohlig-Hellman?
  ------------|------------|---------------------|----------------
  DH (weak)   | 1024       | p-1 (might be smooth)| Depends!
  DH (safe)   | 2048       | 2q (q prime)        | O(sqrt(q))
  DSA         | 256        | q (prime)           | O(sqrt(q))
  ECDSA       | 256        | n (near-prime)      | O(sqrt(n))

  For 128-bit security:
  - Need largest prime factor >= 2^256
  - Safe prime with 2048+ bits OR
  - Prime-order subgroup with 256+ bits
`);

// Example 8: Practical Pohlig-Hellman attack scenario
console.log("Example 8: Attack Scenario");
console.log("-".repeat(50));

console.log(`
  Scenario: A poorly chosen DH group

  Suppose an implementer chooses p where p - 1 = 2 * 3 * 5 * 7 * 11 * ...
  (product of first k primes, called a "primorial")

  p - 1 = 2# = 2       (1 bit security!)
  p - 1 = 6# = 6       (still trivial)
  p - 1 = 30# = 30     (trivial)
  p - 1 = 2310# = 2310 (still trivial)

  Even 223092870 = 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23
  gives only sqrt(23) ~ 5 operations!

  This is why we NEVER use random p. Always use:
  1. Well-known groups (RFC 2409, RFC 3526)
  2. Safe primes (p = 2q + 1)
  3. Prime-order subgroups (DSA-style)
`);

// Example 9: Relation to index calculus
console.log("Example 9: Beyond Pohlig-Hellman");
console.log("-".repeat(50));

console.log(`
  Pohlig-Hellman exploits smooth GROUP ORDER.

  For large prime-order groups, we need other methods:

  Index Calculus (for (Z/pZ)*):
  - Time: O(exp(sqrt(log p * log log p)))
  - Uses smooth ELEMENTS (not order)
  - Subexponential!

  Pollard's rho:
  - Time: O(sqrt(n)) where n is group order
  - Works in any group
  - For EC groups, this is the best known

  Number Field Sieve variant:
  - Even faster for (Z/pZ)* with large p
  - Why DH needs 2048+ bit primes
  - Doesn't apply to elliptic curves!
`);

console.log("\n=== Summary ===");
console.log(`
  Pohlig-Hellman Algorithm:

  Key insight: DLP in smooth-order groups is easy.

  Algorithm:
  1. Factor the group order: n = prod(pi^ei)
  2. Solve DLP modulo each pi^ei (using digit extraction)
  3. Combine solutions with CRT

  Complexity: O(sum of ei * sqrt(pi))
  - Much faster than O(sqrt(n)) when n is smooth
  - Reduces to largest prime factor

  Security implications:
  - NEVER use groups with smooth order
  - Safe primes: p = 2q + 1 (order ~ 2q)
  - Prime-order subgroups: order is prime
  - Elliptic curves: near-prime order by design

  Key takeaway:
  The security of DLP-based crypto depends on the
  LARGEST PRIME FACTOR of the group order, not just
  the total order size.
`);
