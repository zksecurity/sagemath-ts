/**
 * # Baby-step Giant-step Algorithm
 *
 * The baby-step giant-step (BSGS) algorithm is a time-space tradeoff
 * for solving the discrete logarithm problem. It achieves O(sqrt(n))
 * time using O(sqrt(n)) space.
 *
 * ## Algorithm Idea
 *
 * Let m = ceil(sqrt(n)). We can write any x in [0, n) as:
 *   x = i * m + j  where 0 <= i, j < m
 *
 * So g^x = g^(i*m + j) = g^(i*m) * g^j
 *
 * Rearranging: h * g^(-j) = g^(i*m)
 *
 * Algorithm:
 * 1. Baby steps: Compute g^0, g^1, ..., g^(m-1) and store in table
 * 2. Giant steps: Compute h, h*g^(-m), h*g^(-2m), ... until match
 *
 * When we find h * g^(-j) = g^(i*m), we have x = i*m + j.
 *
 * ## Complexity
 *
 * Time: O(sqrt(n)) group operations
 * Space: O(sqrt(n)) group elements
 */

import { power_mod, gcd, is_prime, isqrt, inverse_mod, factor, formatFactorization } from 'sagemath-ts';

console.log("=== Baby-step Giant-step Algorithm ===\n");

// Baby-step giant-step implementation
function babyGiantStep(g: bigint, h: bigint, p: bigint, order?: bigint): bigint | null {
  const n = order ?? p - 1n;
  const m = isqrt(n) + 1n;

  // Baby steps: build table of g^j for j = 0, 1, ..., m-1
  const babyTable = new Map<bigint, bigint>();
  let gj = 1n;

  for (let j = 0n; j < m; j++) {
    babyTable.set(gj, j);
    gj = (gj * g) % p;
  }

  // Giant step factor: g^(-m) mod p
  const gm = power_mod(g, m, p);
  const gmInv = inverse_mod(gm, p);

  // Giant steps: compute h * g^(-im) and look for match
  let gamma = h;

  for (let i = 0n; i < m; i++) {
    const j = babyTable.get(gamma);
    if (j !== undefined) {
      const x = i * m + j;
      // Verify
      if (power_mod(g, x, p) === h) {
        return x;
      }
    }
    gamma = (gamma * gmInv) % p;
  }

  return null;
}

// Example 1: Step-by-step walkthrough
console.log("Example 1: Baby-step Giant-step Walkthrough");
console.log("-".repeat(50));

const p = 29n;
const g = 2n;  // Generator
const h = 17n; // Target

const n = p - 1n;
const m = isqrt(n) + 1n;

console.log(`  Finding x such that ${g}^x = ${h} mod ${p}`);
console.log(`  Group order n = ${n}`);
console.log(`  m = ceil(sqrt(n)) = ${m}`);

console.log(`\n  Baby Steps: Computing g^j for j = 0, ..., ${m - 1n}`);
console.log("  " + "-".repeat(40));

const babyTable = new Map<bigint, bigint>();
let gj = 1n;

for (let j = 0n; j < m; j++) {
  console.log(`    j = ${j}: ${g}^${j} = ${gj}`);
  babyTable.set(gj, j);
  gj = (gj * g) % p;
}

console.log(`\n  Baby table size: ${babyTable.size} entries`);

console.log(`\n  Giant Steps: Computing h * g^(-im) for i = 0, 1, ...`);
console.log("  " + "-".repeat(40));

const gm = power_mod(g, m, p);
const gmInv = inverse_mod(gm, p);
console.log(`  g^m = ${gm}, g^(-m) = ${gmInv}`);

let gamma = h;
let foundI: bigint | null = null;
let foundJ: bigint | null = null;

for (let i = 0n; i < m; i++) {
  const j = babyTable.get(gamma);
  const match = j !== undefined ? `=== MATCH at j = ${j}!` : '';

  console.log(`    i = ${i}: h * g^(-${i}m) = ${gamma} ${match}`);

  if (j !== undefined) {
    foundI = i;
    foundJ = j;
    break;
  }

  gamma = (gamma * gmInv) % p;
}

if (foundI !== null && foundJ !== null) {
  const x = foundI * m + foundJ;
  console.log(`\n  Result: x = ${foundI} * ${m} + ${foundJ} = ${x}`);
  console.log(`  Verification: ${g}^${x} = ${power_mod(g, x, p)} (should be ${h})`);
}

// Example 2: Comparison with brute force
console.log("\nExample 2: Efficiency Comparison");
console.log("-".repeat(50));

function discreteLogBruteForce(g: bigint, h: bigint, p: bigint): { x: bigint | null; steps: number } {
  let current = 1n;
  let steps = 0;

  for (let x = 0n; x < p - 1n; x++) {
    steps++;
    if (current === h) return { x, steps };
    current = (current * g) % p;
  }

  return { x: null, steps };
}

function discreteLogBSGS(g: bigint, h: bigint, p: bigint): { x: bigint | null; steps: number } {
  const order = p - 1n;
  const m = isqrt(order) + 1n;
  let steps = 0;

  // Baby steps
  const table = new Map<bigint, bigint>();
  let gj = 1n;
  for (let j = 0n; j < m; j++) {
    table.set(gj, j);
    gj = (gj * g) % p;
    steps++;
  }

  // Giant steps
  const gmInv = inverse_mod(power_mod(g, m, p), p);
  let gamma = h;

  for (let i = 0n; i < m; i++) {
    steps++;
    const j = table.get(gamma);
    if (j !== undefined) {
      return { x: i * m + j, steps };
    }
    gamma = (gamma * gmInv) % p;
  }

  return { x: null, steps };
}

console.log("  Comparing brute force vs BSGS:");
console.log();

const testPrimes = [101n, 1009n, 10007n, 100003n];

console.log("  Prime p     | sqrt(p-1) | Brute Steps | BSGS Steps | Speedup");
console.log("  " + "-".repeat(65));

for (const prime of testPrimes) {
  const generator = 2n; // Assume 2 is a generator (or use actual generator)
  const target = power_mod(generator, (prime - 1n) / 2n, prime); // Middle element

  const brute = discreteLogBruteForce(generator, target, prime);
  const bsgs = discreteLogBSGS(generator, target, prime);

  const sqrtOrder = isqrt(prime - 1n);
  const speedup = (brute.steps / bsgs.steps).toFixed(1);

  console.log(`  ${prime.toString().padEnd(11)} | ${sqrtOrder.toString().padEnd(9)} | ${brute.steps.toString().padEnd(11)} | ${bsgs.steps.toString().padEnd(10)} | ${speedup}x`);
}

// Example 3: Why the algorithm works
console.log("\nExample 3: Why Baby-step Giant-step Works");
console.log("-".repeat(50));

console.log(`
  Key insight: Any x in [0, n) can be written as x = im + j
  where m = ceil(sqrt(n)) and 0 <= i, j < m.

  We want: g^x = h
  Substitute: g^(im + j) = h
  Rearrange: g^j = h * g^(-im)

  Baby step: Precompute all g^j for j = 0, 1, ..., m-1
  Giant step: Compute h * g^(-im) for i = 0, 1, ...

  When they match: g^j = h * g^(-im)
  Therefore: g^(im + j) = h, so x = im + j

  Why O(sqrt(n))?
  - We need at most m baby steps
  - We need at most m giant steps
  - m ~ sqrt(n), so total O(sqrt(n))
`);

// Example 4: Space-time tradeoff
console.log("Example 4: Space-Time Tradeoff");
console.log("-".repeat(50));

console.log(`
  BSGS with parameter m:

  m larger  -> More baby steps (time), less giant steps (time)
             -> Larger table (space)

  m smaller -> Fewer baby steps (time), more giant steps (time)
             -> Smaller table (space)

  Optimal: m = sqrt(n) balances both phases

  Variants:
  - m = n^(1/3): O(n^(2/3)) time, O(n^(1/3)) space
  - m = n^(2/3): O(n^(1/3)) time, O(n^(2/3)) space
`);

function bsgsWithParam(g: bigint, h: bigint, p: bigint, m: bigint): { x: bigint | null; babySteps: bigint; giantSteps: bigint } {
  const order = p - 1n;

  // Baby steps
  const table = new Map<bigint, bigint>();
  let gj = 1n;
  let babySteps = 0n;

  for (let j = 0n; j < m; j++) {
    table.set(gj, j);
    gj = (gj * g) % p;
    babySteps++;
  }

  // Giant steps
  const gmInv = inverse_mod(power_mod(g, m, p), p);
  let gamma = h;
  let giantSteps = 0n;

  const maxGiant = (order + m - 1n) / m + 1n;

  for (let i = 0n; i < maxGiant; i++) {
    giantSteps++;
    const j = table.get(gamma);
    if (j !== undefined) {
      return { x: i * m + j, babySteps, giantSteps };
    }
    gamma = (gamma * gmInv) % p;
  }

  return { x: null, babySteps, giantSteps };
}

console.log(`\n  Varying m for p = 10007 (order = 10006, sqrt ~ 100):`);
console.log();

const testP = 10007n;
const testG = 5n;
const testH = power_mod(testG, 5000n, testP);

for (const m of [10n, 50n, 100n, 200n, 500n]) {
  const result = bsgsWithParam(testG, testH, testP, m);
  const total = result.babySteps + result.giantSteps;
  console.log(`    m = ${m.toString().padEnd(4)}: baby = ${result.babySteps.toString().padEnd(4)}, giant = ${result.giantSteps.toString().padEnd(4)}, total = ${total}`);
}

// Example 5: Handling unknown group order
console.log("\nExample 5: When Group Order is Unknown");
console.log("-".repeat(50));

console.log(`
  Standard BSGS needs to know the group order n.

  What if n is unknown?
  - Use an upper bound N >= n
  - Set m = ceil(sqrt(N))
  - Algorithm still works, just uses more steps

  In (Z/pZ)*, the order is exactly p-1 when p is prime.

  For subgroups or elliptic curves, we might need to
  compute the order first (which can also be hard).
`);

// Example 6: Memory-efficient variant (Pollard's rho preview)
console.log("Example 6: Memory Considerations");
console.log("-".repeat(50));

console.log(`
  BSGS requires O(sqrt(n)) storage.

  For a 256-bit group:
  - sqrt(2^256) = 2^128 entries
  - Each entry: ~32 bytes (point) + ~32 bytes (exponent) = 64 bytes
  - Total: 64 * 2^128 bytes = 2^134 bytes

  This is way more than all storage on Earth!

  Solutions:
  1. Pollard's rho: O(sqrt(n)) time, O(1) space
  2. Distributed BSGS: Split table across many machines
  3. Time-space tradeoffs with smaller m
`);

// Example 7: Practical usage
console.log("Example 7: Practical BSGS Implementation Tips");
console.log("-".repeat(50));

console.log(`
  Optimizations for real implementations:

  1. Hash table: Use a good hash for O(1) lookups
     - JavaScript Map works well for moderate sizes

  2. Precompute g^(-m) once:
     - Compute gmInv = g^(-m) mod p
     - Multiply by gmInv instead of dividing

  3. Early termination:
     - Check for match after each giant step
     - Don't compute all baby steps if not needed

  4. Batch processing:
     - Compute multiple DLPs with same baby table
     - Amortize the baby step cost

  5. Parallelization:
     - Giant steps are independent (once table is built)
     - Can distribute across processors
`);

// Demonstrate batch processing
function bsgsBatch(g: bigint, targets: bigint[], p: bigint): Map<bigint, bigint | null> {
  const order = p - 1n;
  const m = isqrt(order) + 1n;

  // Build baby table once
  const babyTable = new Map<bigint, bigint>();
  let gj = 1n;
  for (let j = 0n; j < m; j++) {
    babyTable.set(gj, j);
    gj = (gj * g) % p;
  }

  const gmInv = inverse_mod(power_mod(g, m, p), p);
  const results = new Map<bigint, bigint | null>();

  // Solve for each target
  for (const h of targets) {
    let gamma = h;
    let found: bigint | null = null;

    for (let i = 0n; i < m && found === null; i++) {
      const j = babyTable.get(gamma);
      if (j !== undefined) {
        found = i * m + j;
      }
      gamma = (gamma * gmInv) % p;
    }

    results.set(h, found);
  }

  return results;
}

console.log("\n  Batch BSGS for multiple targets (p = 101):");
const batchP = 101n;
const batchG = 2n;
const batchTargets = [5n, 17n, 50n, 73n, 99n];

const batchResults = bsgsBatch(batchG, batchTargets, batchP);
for (const [h, x] of batchResults) {
  const verify = x !== null ? power_mod(batchG, x, batchP) : 'N/A';
  console.log(`    log_${batchG}(${h}) = ${x} (verify: ${verify})`);
}

// Example 8: Complexity summary
console.log("\nExample 8: Complexity Summary");
console.log("-".repeat(50));

console.log(`
  Baby-step Giant-step:

  Time:  O(sqrt(n)) group operations
         - O(sqrt(n)) baby steps
         - O(sqrt(n)) giant steps
         - Plus O(sqrt(n)) hash table operations

  Space: O(sqrt(n)) group elements
         - Store sqrt(n) pairs (g^j, j)
         - Each ~64-128 bytes typically

  For n = 2^256:
  - sqrt(n) = 2^128 operations
  - At 10^9 ops/sec: 10^29 seconds ~ 3*10^21 years
  - Still too slow for direct attack on modern crypto!

  But for smaller groups (n ~ 2^80 or so):
  - sqrt(n) = 2^40 operations
  - At 10^9 ops/sec: ~10^3 seconds ~ 17 minutes
  - This is why 80-bit security is considered weak now
`);

console.log("\n=== Summary ===");
console.log(`
  Baby-step Giant-step Algorithm:

  Key idea: Write x = im + j and match precomputed values

  Complexity:
  - Time: O(sqrt(n))
  - Space: O(sqrt(n))

  Advantages:
  - Deterministic (always finds answer if it exists)
  - Faster than brute force by factor of sqrt(n)
  - Can batch multiple targets

  Disadvantages:
  - Requires O(sqrt(n)) storage
  - For crypto groups, storage is prohibitive
  - Needs to know (or bound) the group order

  In practice:
  - Good for medium-sized groups (order up to ~2^50)
  - Pollard's rho preferred for larger groups (O(1) space)
  - Neither threatens 256-bit crypto groups

  Next: Pohlig-Hellman algorithm for groups with smooth order
`);
