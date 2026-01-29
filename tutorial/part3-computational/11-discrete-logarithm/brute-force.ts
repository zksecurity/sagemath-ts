/**
 * # Brute Force Discrete Logarithm
 *
 * The naive approach to solving the discrete logarithm problem:
 * try every possible exponent until we find a match.
 *
 * ## Algorithm
 *
 * To find x such that g^x = h in a group of order n:
 * 1. Compute g^0, g^1, g^2, ... until we find g^x = h
 * 2. Return x
 *
 * ## Complexity
 *
 * Time: O(n) group operations
 * Space: O(1)
 *
 * For cryptographic groups (n ~ 2^256), this is completely infeasible.
 * But understanding brute force helps appreciate better algorithms.
 */

import { power_mod, gcd, is_prime, euler_phi, factor, formatFactorization } from 'sagemath-ts';

console.log("=== Brute Force Discrete Logarithm ===\n");

// Naive brute force implementation
function discreteLogBruteForce(g: bigint, h: bigint, p: bigint, maxIterations?: bigint): bigint | null {
  const order = p - 1n;
  const limit = maxIterations !== undefined ? maxIterations : order;

  let current = 1n;

  for (let x = 0n; x < limit; x++) {
    if (current === h) {
      return x;
    }
    current = (current * g) % p;
  }

  return null;
}

// Example 1: Basic brute force
console.log("Example 1: Basic Brute Force Search");
console.log("-".repeat(50));

const p = 31n;
const g = 3n; // Generator of (Z/31Z)*

console.log(`  Finding discrete logs in (Z/${p}Z)* with generator g = ${g}`);
console.log(`  Group order: ${p - 1n}\n`);

// Trace the search
function discreteLogVerbose(g: bigint, h: bigint, p: bigint): bigint | null {
  console.log(`  Searching for x such that ${g}^x = ${h} mod ${p}:`);

  let current = 1n;
  for (let x = 0n; x < p - 1n; x++) {
    if (x <= 10n || current === h) {
      console.log(`    x = ${x}: ${g}^${x} = ${current} ${current === h ? '=== FOUND!' : ''}`);
    } else if (x === 11n) {
      console.log(`    ...`);
    }

    if (current === h) {
      return x;
    }
    current = (current * g) % p;
  }
  return null;
}

discreteLogVerbose(g, 25n, p);

// Example 2: Timing brute force
console.log("\nExample 2: Brute Force Performance");
console.log("-".repeat(50));

function timeDiscreteLog(g: bigint, h: bigint, p: bigint): { x: bigint | null; time: number; steps: bigint } {
  const start = performance.now();

  let current = 1n;
  let steps = 0n;

  for (let x = 0n; x < p - 1n; x++) {
    steps++;
    if (current === h) {
      return { x, time: performance.now() - start, steps };
    }
    current = (current * g) % p;
  }

  return { x: null, time: performance.now() - start, steps };
}

// Test with increasing group sizes
const testCases: Array<{ p: bigint; g: bigint; h: bigint }> = [
  { p: 101n, g: 2n, h: 50n },
  { p: 1009n, g: 11n, h: 500n },
  { p: 10007n, g: 5n, h: 1234n },
  { p: 100003n, g: 2n, h: 54321n },
];

console.log("  Prime p    | Steps    | Time (ms) | Answer");
console.log("  " + "-".repeat(50));

for (const { p, g, h } of testCases) {
  const result = timeDiscreteLog(g, h, p);
  console.log(`  ${p.toString().padEnd(10)} | ${result.steps.toString().padEnd(8)} | ${result.time.toFixed(3).padStart(9)} | x = ${result.x}`);
}

// Example 3: Why brute force is hopeless for crypto
console.log("\nExample 3: Why Brute Force is Hopeless for Cryptography");
console.log("-".repeat(50));

console.log(`
  Cryptographic group sizes:

  Group Type          | Order (bits) | Brute Force Steps | At 10^9/sec
  --------------------|--------------|-------------------|-------------
  Toy example         | 10           | ~1000             | microseconds
  Small example       | 20           | ~10^6             | milliseconds
  Weak (don't use!)   | 40           | ~10^12            | ~17 minutes
  ECDSA (256-bit)     | 256          | ~10^77            | 10^60 years
  DH (2048-bit)       | 2048         | ~10^616           | 10^599 years

  Age of universe: ~4 x 10^17 seconds
  10^9 ops/sec for 10^60 years = 10^67 operations

  Bottom line: Brute force is absolutely impossible for real crypto.
`);

// Example 4: Optimized brute force
console.log("Example 4: Optimizations to Brute Force");
console.log("-".repeat(50));

console.log(`
  Some small optimizations (still O(n)):

  1. Early termination: Stop when we find h
     - Average case: n/2 steps (still linear)

  2. Skip known impossible values:
     - If h is not a quadratic residue but g^2 is, skip odd x

  3. Use Montgomery multiplication:
     - Faster modular multiplication
     - Constant factor improvement

  4. Parallelization:
     - Split range among processors
     - Linear speedup with P processors

  None of these change the fundamental O(n) complexity.
`);

// Parallel simulation (not actually parallel, but shows the idea)
function bruteForceParallel(g: bigint, h: bigint, p: bigint, numWorkers: number): { x: bigint | null; stepsPerWorker: bigint } {
  const order = p - 1n;
  const chunkSize = order / BigInt(numWorkers);

  for (let worker = 0; worker < numWorkers; worker++) {
    const start = BigInt(worker) * chunkSize;
    const end = worker === numWorkers - 1 ? order : start + chunkSize;

    // Compute g^start
    let current = power_mod(g, start, p);

    for (let x = start; x < end; x++) {
      if (current === h) {
        const stepsPerWorker = (x - start) + 1n;
        return { x, stepsPerWorker };
      }
      current = (current * g) % p;
    }
  }

  return { x: null, stepsPerWorker: chunkSize };
}

console.log("\n  Simulated parallel brute force (p = 10007, h = 5000):");
for (const workers of [1, 2, 4, 8]) {
  const result = bruteForceParallel(5n, 5000n, 10007n, workers);
  console.log(`    ${workers} worker(s): found x = ${result.x} (max steps/worker: ${10007n / BigInt(workers)})`);
}

// Example 5: The search space
console.log("\nExample 5: Visualizing the Search Space");
console.log("-".repeat(50));

const smallP = 17n;
const smallG = 3n;

console.log(`  Group (Z/${smallP}Z)* with generator ${smallG}:`);
console.log(`  Order: ${smallP - 1n}`);
console.log();

// Create a table of all discrete logs
const dlogs = new Map<bigint, bigint>();
let val = 1n;
for (let x = 0n; x < smallP - 1n; x++) {
  dlogs.set(val, x);
  val = (val * smallG) % smallP;
}

console.log("  x  |  g^x  | Elements in order of g^x");
console.log("  " + "-".repeat(40));

val = 1n;
for (let x = 0n; x < smallP - 1n; x++) {
  console.log(`  ${x.toString().padStart(2)} |  ${val.toString().padStart(3)}  | ${'*'.repeat(Number(val))}`);
  val = (val * smallG) % smallP;
}

console.log(`
  The sequence g^0, g^1, ... visits all elements of the group
  in a pseudo-random order. Brute force must follow this path.
`);

// Example 6: When brute force might be acceptable
console.log("Example 6: When Brute Force Might Be Acceptable");
console.log("-".repeat(50));

console.log(`
  Brute force can be practical when:

  1. Group order is very small (educational, testing)
     - Orders up to ~10^6: milliseconds
     - Orders up to ~10^9: seconds

  2. Looking for "small" exponents
     - Some protocols use small exponents
     - Check x = 0, 1, 2, ... up to a bound

  3. Bounded discrete log
     - Know that x < B for some small B
     - Only need O(B) time instead of O(n)

  4. Preprocessing step
     - Try small x before using fancier algorithms
     - Catches trivial cases cheaply
`);

// Example: bounded discrete log
function boundedDiscreteLog(g: bigint, h: bigint, p: bigint, bound: bigint): bigint | null {
  let current = 1n;
  for (let x = 0n; x < bound; x++) {
    if (current === h) return x;
    current = (current * g) % p;
  }
  return null;
}

console.log("\n  Bounded discrete log (checking x < 1000):");
const p_bounded = 1000003n;
const g_bounded = 5n;

// Create a target with known small exponent
const knownX = 42n;
const h_bounded = power_mod(g_bounded, knownX, p_bounded);

const result = boundedDiscreteLog(g_bounded, h_bounded, p_bounded, 1000n);
console.log(`    g = ${g_bounded}, h = ${h_bounded}, p = ${p_bounded}`);
console.log(`    Result: x = ${result} (correct: ${result === knownX})`);

// Example 7: Space-time tradeoff preview
console.log("\nExample 7: Preview of Better Algorithms");
console.log("-".repeat(50));

console.log(`
  Brute force: O(n) time, O(1) space

  Can we trade space for time?

  Baby-step Giant-step:
  - O(sqrt(n)) time, O(sqrt(n)) space
  - Store sqrt(n) values in a table
  - Much faster: 2^128 vs 2^256 for 256-bit groups

  Pollard's rho:
  - O(sqrt(n)) time, O(1) space
  - Uses random walk + cycle detection
  - No large table needed

  Pohlig-Hellman:
  - O(sum of sqrt(p_i)) for order n = prod(p_i^e_i)
  - Exploits smooth group orders
  - Why we need prime-order subgroups

  Index Calculus:
  - O(exp(sqrt(log n * log log n))) for (Z/pZ)*
  - Sub-exponential! Much faster for large n
  - Doesn't work for elliptic curves

  These better algorithms are covered in the following tutorials.
`);

// Example 8: Summary visualization
console.log("Example 8: Algorithm Comparison");
console.log("-".repeat(50));

console.log(`
  Steps needed for different algorithms (256-bit group):

  Algorithm          | Time Complexity | Steps (n = 2^256)
  -------------------|-----------------|------------------
  Brute Force        | O(n)            | 2^256
  Baby-Giant Step    | O(sqrt(n))      | 2^128
  Pollard's rho      | O(sqrt(n))      | 2^128
  Pohlig-Hellman*    | O(sqrt(q))      | 2^128 (if q ~ n)
  Index Calculus**   | sub-exp         | ~2^60 (varies)

  * For group order with largest prime factor q
  ** Only for (Z/pZ)*, not elliptic curves

  The jump from O(n) to O(sqrt(n)) is huge!
  2^256 -> 2^128 operations is the difference between
  "impossible" and "expensive but potentially feasible."
`);

console.log("\n=== Summary ===");
console.log(`
  Brute Force Discrete Log:

  Algorithm: Try x = 0, 1, 2, ... until g^x = h

  Complexity:
  - Time: O(n) where n is the group order
  - Space: O(1)

  Practical limits:
  - n ~ 10^6: milliseconds
  - n ~ 10^9: seconds
  - n ~ 10^12: hours
  - n ~ 2^256: impossible

  When to use:
  - Educational examples
  - Very small groups
  - Bounded discrete log (small exponents)
  - Testing/verification

  Next: Baby-step Giant-step achieves O(sqrt(n)) using a table.
`);
