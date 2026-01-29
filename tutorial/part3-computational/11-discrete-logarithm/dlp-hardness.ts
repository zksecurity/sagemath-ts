/**
 * # Why the Discrete Logarithm Problem is Hard
 *
 * This tutorial explores why DLP is believed to be computationally hard
 * and how this assumption forms the foundation of modern cryptography.
 *
 * ## The DLP Assumption
 *
 * Informally: In a properly chosen group G of order n, given g and h = g^x,
 * computing x requires time proportional to sqrt(n) or worse.
 *
 * ## Related Assumptions
 *
 * - CDH (Computational Diffie-Hellman): Given g^a and g^b, compute g^(ab)
 * - DDH (Decisional Diffie-Hellman): Distinguish g^(ab) from random
 *
 * These assumptions are what allow us to build secure cryptosystems.
 */

import { power_mod, gcd, is_prime, isqrt, inverse_mod, factor, formatFactorization, prime_range } from 'sagemath-ts';

console.log("=== Why DLP is Hard: Cryptographic Assumptions ===\n");

// Example 1: The one-way function perspective
console.log("Example 1: Exponentiation as a One-Way Function");
console.log("-".repeat(60));

console.log(`
  A one-way function is:
  - Easy to compute in the forward direction
  - Hard to invert

  Exponentiation in a group:
  - Forward: x -> g^x (O(log x) multiplications)
  - Inverse: h -> log_g(h) (no known efficient algorithm)

  This asymmetry is the foundation of DL-based cryptography.
`);

const p = 101n;
const g = 2n;

console.log(`  Example in (Z/${p}Z)*:`);
console.log(`\n  Forward direction (easy):`);

const start_forward = performance.now();
const results: Array<{ x: bigint; gx: bigint }> = [];
for (const x of [10n, 50n, 99n]) {
  const gx = power_mod(g, x, p);
  results.push({ x, gx });
  console.log(`    ${g}^${x} mod ${p} = ${gx}`);
}
console.log(`    Time: < 1ms`);

console.log(`\n  Reverse direction (harder):`);
const start_reverse = performance.now();
for (const { x, gx } of results) {
  // Brute force search
  let found = 0n;
  for (let i = 0n; i < p - 1n; i++) {
    if (power_mod(g, i, p) === gx) {
      found = i;
      break;
    }
  }
  console.log(`    log_${g}(${gx}) = ${found}`);
}
const reverse_time = performance.now() - start_reverse;
console.log(`    Time: ${reverse_time.toFixed(3)}ms (even for this tiny group)`);

// Example 2: Known algorithms and their complexities
console.log("\nExample 2: DLP Algorithms and Their Limits");
console.log("-".repeat(60));

console.log(`
  Best known algorithms for the DLP:

  Group Type              | Best Algorithm        | Complexity
  ------------------------|----------------------|---------------------------
  Generic group           | Pollard's rho        | O(sqrt(n))
  (Z/pZ)*                 | Index calculus/NFS   | L_p[1/3, c] = exp(c(log p)^(1/3)(log log p)^(2/3))
  Elliptic curve E(F_p)   | Pollard's rho        | O(sqrt(n))

  Key observation:
  - For (Z/pZ)*, subexponential algorithms exist
  - For elliptic curves, only exponential algorithms are known

  This is why EC cryptography can use smaller key sizes:
  256-bit EC curve ~ 3072-bit DH modulus (both ~128-bit security)
`);

// Example 3: Generic group model
console.log("Example 3: The Generic Group Model");
console.log("-".repeat(60));

console.log(`
  In the "generic group model" (Shoup 1997):

  - The only operations allowed are group operations
  - No exploitation of the group's representation
  - Any algorithm must be "black-box"

  Lower bound: Any generic algorithm requires Omega(sqrt(n)) operations.

  Pollard's rho achieves O(sqrt(n)), so it's optimal for generic groups!

  However:
  - (Z/pZ)* is NOT a generic group (has exploitable structure)
  - Index calculus exploits this structure
  - EC groups are "more generic" (no known structure to exploit)
`);

// Example 4: Why structure matters
console.log("Example 4: Why Group Structure Matters");
console.log("-".repeat(60));

console.log(`
  (Z/pZ)* has exploitable structure:

  1. Pohlig-Hellman: Exploits smooth order
     - If p-1 = small primes, DLP is easy

  2. Index calculus: Exploits smoothness of elements
     - Some elements are "smooth" (factor into small primes)
     - Build relations to solve DLP

  3. Number Field Sieve: Algebraic number theory
     - Most powerful method for large p
     - Subexponential complexity

  Elliptic curves lack this structure:
  - No analog of "smooth elements"
  - Index calculus doesn't apply
  - Only generic attacks work
`);

// Demonstrate structure in Z/pZ*
console.log("\n  Smooth elements in (Z/pZ)* for p = 101:");
const smoothBound = 10n;
const smoothElements: Array<{ x: bigint; factors: string }> = [];

for (let x = 2n; x < 30n; x++) {
  const f = factor(x);
  const largest = f.reduce((max, [prime]) => prime > max ? prime : max, 1n);
  if (largest <= smoothBound) {
    smoothElements.push({ x, factors: formatFactorization(f) });
  }
}

console.log(`    ${smoothBound}-smooth elements < 30: ${smoothElements.map(e => `${e.x}(=${e.factors})`).join(', ')}`);
console.log(`    These smooth elements enable index calculus attacks.`);

// Example 5: The Diffie-Hellman assumptions
console.log("\nExample 5: The Diffie-Hellman Assumptions");
console.log("-".repeat(60));

console.log(`
  Three related assumptions (in increasing strength):

  1. DLP (Discrete Logarithm Problem):
     Given (g, g^x), find x.

  2. CDH (Computational Diffie-Hellman):
     Given (g, g^a, g^b), compute g^(ab).

  3. DDH (Decisional Diffie-Hellman):
     Distinguish (g, g^a, g^b, g^(ab)) from (g, g^a, g^b, g^c).

  Hierarchy:
     DDH <= CDH <= DLP
     (If you can break DLP, you can break CDH and DDH)

  Some groups have DDH easy but CDH hard!
`);

// Demonstrate the relationship
console.log("\n  Demonstrating the hierarchy:");

const dh_g = 5n;
const dh_p = 101n;
const a = 13n;
const b = 17n;

const ga = power_mod(dh_g, a, dh_p);
const gb = power_mod(dh_g, b, dh_p);
const gab = power_mod(dh_g, a * b, dh_p);

console.log(`    g = ${dh_g}, p = ${dh_p}`);
console.log(`    a = ${a} (secret), b = ${b} (secret)`);
console.log(`    g^a = ${ga}, g^b = ${gb}`);
console.log(`    g^(ab) = ${gab}`);

console.log(`\n    DLP: From g^a = ${ga}, find a = ?`);
console.log(`    CDH: From (g^a, g^b) = (${ga}, ${gb}), compute g^(ab) = ?`);
console.log(`    DDH: Is (${ga}, ${gb}, ${gab}) a valid DH triple or random?`);

// Example 6: Security levels
console.log("\nExample 6: Security Levels and Key Sizes");
console.log("-".repeat(60));

console.log(`
  Security level = log2(operations needed for best known attack)

  Target Security | EC Curve Bits | DH/RSA Modulus Bits | Notes
  ----------------|---------------|---------------------|----------------
  80 bits         | 160           | 1024                | Legacy, weak
  112 bits        | 224           | 2048                | Minimum today
  128 bits        | 256           | 3072                | Standard
  192 bits        | 384           | 7680                | High security
  256 bits        | 512           | 15360               | Very high

  The disparity comes from different attack complexities:
  - EC: O(sqrt(n)) = O(2^(bits/2)) -> 256-bit curve = 128-bit security
  - DH: subexponential -> need larger modulus for same security
`);

// Example 7: Why we believe DLP is hard
console.log("Example 7: Why We Believe DLP is Hard");
console.log("-".repeat(60));

console.log(`
  Evidence for DLP hardness:

  1. Decades of cryptanalysis:
     - Studied since 1976 (Diffie-Hellman)
     - No polynomial-time algorithm found
     - Best improvements are subexponential (for finite fields)

  2. Reduction from other problems:
     - Breaking DLP in some groups implies factoring
     - Connected to other believed-hard problems

  3. Records show slow progress:
     - DLP record (finite field): ~800 bits (2019)
     - Compare to RSA record: ~830 bits (2020)
     - Progress requires massive computation

  4. NIST/NSA recommendations:
     - Based on extensive analysis
     - Conservative estimates for future progress

  HOWEVER:
  - No proof that DLP is hard (like most crypto assumptions)
  - Quantum computers break DLP (Shor's algorithm)
  - Unexpected algorithmic breakthroughs are possible
`);

// Example 8: Quantum threat
console.log("Example 8: The Quantum Threat");
console.log("-".repeat(60));

console.log(`
  Shor's Algorithm (1994):

  Quantum computers can solve DLP in polynomial time!
  - Time: O((log n)^3) quantum operations
  - Completely breaks DH, DSA, ECDH, ECDSA

  Current status:
  - Small-scale quantum computers exist
  - Largest DLP solved quantumly: very small numbers
  - Cryptographically relevant: estimated 10-30+ years away

  Post-quantum alternatives (not based on DLP):
  - Lattice-based: NTRU, Kyber, Dilithium (NIST standardized)
  - Code-based: McEliece
  - Hash-based: SPHINCS+, XMSS
  - Isogeny-based: SIKE (broken in 2022!), CSIDH

  Transition:
  - "Harvest now, decrypt later" threat
  - Migration to PQ crypto is urgent for long-term secrets
  - Hybrid approaches during transition
`);

// Example 9: Proper group selection
console.log("Example 9: Choosing Secure Groups");
console.log("-".repeat(60));

console.log(`
  DO use:
  - Safe primes: p = 2q + 1 where q is prime
  - Prime-order subgroups (Schnorr groups)
  - Well-vetted elliptic curves (P-256, Curve25519)
  - Standard groups from RFCs

  DON'T use:
  - Random primes (might have smooth p-1)
  - Small groups (< 256 bits for EC, < 2048 for DH)
  - Custom/untested groups
  - Groups without large prime-order subgroup
`);

// Example of checking a group
console.log("\n  Checking group quality:");

const checkPrimes = [
  { p: 23n, label: "Tiny (toy)" },
  { p: 1019n, label: "Small" },
  { p: 4099n, label: "Smooth p-1 (2049 = 3 * 683)" },
  { p: 10007n, label: "Large factor (5003)" },
];

for (const { p, label } of checkPrimes) {
  const order = p - 1n;
  const f = factor(order);
  const largest = f.reduce((max, [prime]) => prime > 0n && prime > max ? prime : max, 1n);
  const securityBits = largest.toString(2).length / 2; // sqrt for generic attack

  console.log(`    p = ${p} (${label})`);
  console.log(`      Order: ${order} = ${formatFactorization(f)}`);
  console.log(`      Largest prime factor: ${largest}`);
  console.log(`      Effective security: ~${securityBits.toFixed(0)} bits`);
  console.log();
}

// Example 10: The bigger picture
console.log("Example 10: DLP in the Cryptographic Landscape");
console.log("-".repeat(60));

console.log(`
  DLP is one of several "hard problems" used in cryptography:

  Problem               | Used In               | Best Classical   | Quantum
  ----------------------|----------------------|------------------|----------
  Integer Factorization | RSA                  | Subexponential   | Polynomial
  DLP in (Z/pZ)*        | DH, DSA, ElGamal     | Subexponential   | Polynomial
  DLP in EC             | ECDH, ECDSA          | Exponential      | Polynomial
  Shortest Vector (SVP) | Lattice crypto       | Exponential      | Exponential*
  Syndrome Decoding     | Code-based crypto    | Exponential      | Exponential*

  * Quantum speedups exist but remain exponential

  The DLP family (including EC) will be obsolete when large
  quantum computers arrive. Plan your crypto accordingly!
`);

console.log("\n=== Summary ===");
console.log(`
  Why the Discrete Logarithm Problem is Hard:

  1. No efficient algorithm is known:
     - Generic groups: O(sqrt(n)) is optimal
     - Finite fields: Subexponential, but still hard
     - Elliptic curves: Only O(sqrt(n)) known

  2. Decades of cryptanalysis:
     - Intensively studied since 1976
     - Progress is slow and incremental
     - Current records require massive computation

  3. Related to other hard problems:
     - Factoring (for some groups)
     - CDH, DDH assumptions

  Caveats:
  - Not proven to be hard (just believed)
  - Quantum computers break DLP
  - Group selection is critical

  The takeaway:
  DLP hardness is a foundational assumption of modern cryptography.
  When choosing groups, ensure:
  - Large enough (256+ bits for EC, 2048+ for DH)
  - Prime or near-prime order
  - Well-vetted and standardized
  - Planning for post-quantum transition
`);
