/**
 * # Pollard's Rho Algorithm
 *
 * Pollard's rho is a probabilistic factorization algorithm that uses
 * the birthday paradox to find factors. It's much faster than trial
 * division for numbers with medium-sized factors.
 *
 * ## The Birthday Paradox
 *
 * In a group of 23 people, there's a >50% chance two share a birthday.
 * With sqrt(N) random samples from N items, expect a collision.
 *
 * ## Algorithm Idea
 *
 * 1. We want to factor n = p * q
 * 2. Generate a "pseudo-random" sequence: x_{i+1} = f(x_i) mod n
 * 3. After O(sqrt(p)) steps, we expect x_i = x_j mod p for some i != j
 * 4. This means p | (x_i - x_j) but likely n does not divide (x_i - x_j)
 * 5. So gcd(x_i - x_j, n) gives us a factor!
 *
 * ## Complexity
 *
 * Expected time: O(n^(1/4)) to find a factor of size sqrt(n)
 * This is MUCH better than trial division's O(n^(1/2))!
 *
 * ## The "Rho" Shape
 *
 * The sequence eventually cycles (finite group), looking like ρ:
 * x_0 -> x_1 -> ... -> x_t -> x_{t+1} -> ... -> x_t (cycle)
 */

import { gcd, factor, formatFactorization, is_prime, isqrt } from 'sagemath-ts';

console.log("=== Pollard's Rho Factorization Algorithm ===\n");

// Standard polynomial for Pollard's rho: f(x) = x^2 + c mod n
function pollardRhoFunction(x: bigint, c: bigint, n: bigint): bigint {
  return (x * x + c) % n;
}

// Basic Pollard's rho implementation
function pollardRhoBasic(n: bigint, c: bigint = 1n): bigint {
  if (n === 1n) return 1n;
  if (n % 2n === 0n) return 2n;

  let x = 2n;
  let y = 2n;
  let d = 1n;

  while (d === 1n) {
    x = pollardRhoFunction(x, c, n);           // Tortoise: one step
    y = pollardRhoFunction(y, c, n);
    y = pollardRhoFunction(y, c, n);           // Hare: two steps

    d = gcd(x > y ? x - y : y - x, n);
  }

  return d;
}

// Example 1: The birthday paradox
console.log("Example 1: Birthday Paradox Illustration");
console.log("-".repeat(50));

console.log(`
  The birthday paradox: In a room of n people, the probability
  that at least two share a birthday is higher than intuition suggests.

  People | P(collision)
  -------|-------------
  10     | 11.7%
  20     | 41.1%
  23     | 50.7%
  30     | 70.6%
  50     | 97.0%

  General principle: With sqrt(N) random samples from N items,
  expect a collision with constant probability.
`);

// Simulate birthday paradox
function birthdaySimulation(people: number, trials: number = 1000): number {
  let collisions = 0;

  for (let t = 0; t < trials; t++) {
    const birthdays = new Set<number>();
    let hasCollision = false;

    for (let i = 0; i < people; i++) {
      const bday = Math.floor(Math.random() * 365);
      if (birthdays.has(bday)) {
        hasCollision = true;
        break;
      }
      birthdays.add(bday);
    }

    if (hasCollision) collisions++;
  }

  return collisions / trials * 100;
}

console.log("  Simulated collision probabilities:");
for (const people of [10, 20, 23, 30, 50]) {
  const prob = birthdaySimulation(people, 1000);
  console.log(`    ${people} people: ${prob.toFixed(1)}% collision rate`);
}

// Example 2: Floyd's cycle detection
console.log("\nExample 2: Floyd's Cycle Detection (Tortoise and Hare)");
console.log("-".repeat(50));

console.log(`
  The sequence x_{i+1} = f(x_i) mod n must eventually cycle.

  Floyd's algorithm finds a cycle using O(1) space:
  - Tortoise moves one step per iteration
  - Hare moves two steps per iteration
  - When they meet, we've detected the cycle
`);

function demonstrateCycle(n: bigint, c: bigint, start: bigint, maxSteps: number = 20): void {
  const sequence: bigint[] = [start];
  let x = start;

  for (let i = 0; i < maxSteps; i++) {
    x = pollardRhoFunction(x, c, n);
    sequence.push(x);
  }

  console.log(`  Sequence mod ${n} with f(x) = x^2 + ${c}:`);
  console.log(`    ${sequence.slice(0, 15).join(" -> ")}...`);
}

demonstrateCycle(17n, 1n, 2n);

// Example 3: Basic Pollard rho in action
console.log("\nExample 3: Pollard's Rho in Action");
console.log("-".repeat(50));

function pollardRhoVerbose(n: bigint, c: bigint = 1n): bigint {
  console.log(`  Factoring n = ${n} with c = ${c}`);

  let x = 2n;
  let y = 2n;
  let d = 1n;
  let steps = 0;

  while (d === 1n) {
    x = pollardRhoFunction(x, c, n);
    y = pollardRhoFunction(y, c, n);
    y = pollardRhoFunction(y, c, n);
    steps++;

    d = gcd(x > y ? x - y : y - x, n);

    if (steps <= 10 || d !== 1n) {
      console.log(`    Step ${steps}: x=${x}, y=${y}, gcd(|x-y|, n)=${d}`);
    } else if (steps === 11) {
      console.log(`    ...`);
    }
  }

  console.log(`  Found factor ${d} in ${steps} steps`);
  return d;
}

pollardRhoVerbose(91n); // 7 * 13

// Example 4: Factoring larger numbers
console.log("\nExample 4: Factoring Larger Numbers");
console.log("-".repeat(50));

const testNumbers = [
  { n: 1001n, label: "7 * 11 * 13" },
  { n: 10403n, label: "101 * 103" },
  { n: 122987n, label: "prime check" },
  { n: 1000003n * 1000033n, label: "two 7-digit primes" },
];

for (const { n, label } of testNumbers) {
  if (is_prime(n)) {
    console.log(`  n = ${n} (${label}): PRIME`);
    continue;
  }

  const start = performance.now();
  const factor1 = pollardRhoBasic(n);
  const time = performance.now() - start;

  console.log(`  n = ${n} (${label})`);
  console.log(`    Factor found: ${factor1} (${time.toFixed(2)}ms)`);
}

// Example 5: Brent's improvement
console.log("\nExample 5: Brent's Improvement");
console.log("-".repeat(50));

console.log(`
  Brent's algorithm is an optimization of Pollard's rho:
  - Uses a different cycle detection method
  - Only computes GCD periodically (batched)
  - About 24% faster in practice

  Key ideas:
  1. Track powers of 2 in the cycle length
  2. Accumulate products before computing GCD
  3. If GCD = n, backtrack and try individual steps
`);

function brentRho(n: bigint, c: bigint = 1n): bigint {
  if (n === 1n) return 1n;
  if (n % 2n === 0n) return 2n;

  let y = 2n;
  let r = 1n;
  let q = 1n;
  let g = 1n;
  let ys = 0n;
  let x = 0n;

  while (g === 1n) {
    x = y;

    for (let i = 0n; i < r; i++) {
      y = pollardRhoFunction(y, c, n);
    }

    let k = 0n;
    while (k < r && g === 1n) {
      ys = y;
      const bound = r - k < 128n ? r - k : 128n;

      for (let i = 0n; i < bound; i++) {
        y = pollardRhoFunction(y, c, n);
        const diff = x > y ? x - y : y - x;
        q = (q * diff) % n;
      }

      g = gcd(q, n);
      k += 128n;
    }

    r *= 2n;
  }

  if (g === n) {
    // Backtrack
    while (true) {
      ys = pollardRhoFunction(ys, c, n);
      g = gcd(x > ys ? x - ys : ys - x, n);
      if (g > 1n) break;
    }
  }

  return g;
}

// Compare basic and Brent's
console.log("\n  Comparing Basic vs Brent's algorithm:");

const testN = 1000000007n * 1000000009n;
console.log(`  n = ${testN}`);

const startBasic = performance.now();
const factorBasic = pollardRhoBasic(testN);
const timeBasic = performance.now() - startBasic;

const startBrent = performance.now();
const factorBrent = brentRho(testN);
const timeBrent = performance.now() - startBrent;

console.log(`    Basic:  Found ${factorBasic} in ${timeBasic.toFixed(2)}ms`);
console.log(`    Brent:  Found ${factorBrent} in ${timeBrent.toFixed(2)}ms`);

// Example 6: Choosing the polynomial
console.log("\nExample 6: Choosing the Polynomial");
console.log("-".repeat(50));

console.log(`
  Standard choice: f(x) = x^2 + c mod n

  The constant c should NOT be 0 or -2:
  - c = 0: x^2 has a fixed point at 0
  - c = -2: x^2 - 2 satisfies x_n = x_0^(2^n) mod n (bad mixing)

  If Pollard's rho fails with one c, try another.
  Usually c = 1 works well.
`);

function factorWithRetry(n: bigint, maxAttempts: number = 5): bigint | null {
  for (let c = 1n; c <= BigInt(maxAttempts); c++) {
    if (c === 0n || c === n - 2n) continue;

    const factor1 = brentRho(n, c);

    if (factor1 !== n && factor1 !== 1n) {
      console.log(`    Found factor ${factor1} with c = ${c}`);
      return factor1;
    }
  }
  return null;
}

const hardN = 323n; // = 17 * 19
console.log(`  Factoring n = ${hardN}:`);
factorWithRetry(hardN);

// Example 7: Complexity analysis
console.log("\nExample 7: Complexity Analysis");
console.log("-".repeat(50));

console.log(`
  Pollard's rho expected complexity:

  For n = p * q where p <= q:
  - Need O(sqrt(p)) iterations to find p | gcd(x_i - x_j, n)
  - Each iteration: O(log n)^2 for multiplication
  - Overall: O(n^(1/4) * (log n)^2) to find the smallest factor

  Comparison with trial division for n = p * q, p ~ q ~ sqrt(n):

  Algorithm       | Time Complexity
  ----------------|------------------
  Trial division  | O(sqrt(n)) = O(n^(1/2))
  Pollard's rho   | O(n^(1/4))

  For a 40-digit semiprime:
  - Trial division: ~10^20 operations
  - Pollard's rho:  ~10^10 operations

  10 billion times faster!
`);

// Example 8: Repeated squaring issue
console.log("Example 8: Full Factorization Using Pollard's Rho");
console.log("-".repeat(50));

function fullFactor(n: bigint): Array<[bigint, bigint]> {
  const factors: Array<[bigint, bigint]> = [];

  if (n <= 1n) return factors;

  // Handle negative
  if (n < 0n) {
    factors.push([-1n, 1n]);
    n = -n;
  }

  // Trial division for small factors
  for (const p of [2n, 3n, 5n, 7n, 11n, 13n]) {
    let exp = 0n;
    while (n % p === 0n) {
      exp++;
      n /= p;
    }
    if (exp > 0n) factors.push([p, exp]);
  }

  // Use Pollard's rho for remaining factors
  const stack: bigint[] = [n];

  while (stack.length > 0) {
    const current = stack.pop()!;

    if (current === 1n) continue;

    if (is_prime(current)) {
      // Add to factors
      const existing = factors.find(([p]) => p === current);
      if (existing) {
        existing[1]++;
      } else {
        factors.push([current, 1n]);
      }
      continue;
    }

    // Factor using Pollard's rho
    let divisor = brentRho(current);

    // Handle failure (shouldn't happen often)
    if (divisor === current || divisor === 1n) {
      // Try different c values
      for (let c = 2n; c <= 10n; c++) {
        divisor = brentRho(current, c);
        if (divisor !== current && divisor !== 1n) break;
      }
    }

    if (divisor !== current && divisor !== 1n) {
      stack.push(divisor);
      stack.push(current / divisor);
    } else {
      // Fallback: assume prime (shouldn't happen)
      factors.push([current, 1n]);
    }
  }

  // Sort factors by prime
  factors.sort((a, b) => (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0));

  // Merge duplicates
  const merged: Array<[bigint, bigint]> = [];
  for (const [p, e] of factors) {
    const last = merged[merged.length - 1];
    if (last && last[0] === p) {
      last[1] += e;
    } else {
      merged.push([p, e]);
    }
  }

  return merged;
}

const composites = [
  360n,
  1001n,
  123456789n,
  2n ** 32n + 1n, // Fermat number F5
  1000000007n * 1000000009n * 1000000021n,
];

for (const n of composites) {
  const factors1 = fullFactor(n);
  const factors2 = factor(n);
  const str = factors1.map(([p, e]) => (e === 1n ? `${p}` : `${p}^${e}`)).join(" * ");
  console.log(`  ${n} = ${str}`);
}

console.log("\n=== Summary ===");
console.log(`
  Pollard's Rho Algorithm:

  Key ideas:
  1. Birthday paradox: Expect collision in O(sqrt(p)) samples mod p
  2. Pseudo-random walk: f(x) = x^2 + c mod n
  3. Cycle detection: Floyd's or Brent's algorithm
  4. GCD reveals factor: gcd(x_i - x_j, n) when x_i = x_j mod p

  Complexity: O(n^(1/4)) vs trial division's O(n^(1/2))

  Advantages:
  - Much faster than trial division for medium factors
  - O(1) space with Floyd's cycle detection
  - Simple to implement

  Limitations:
  - Doesn't work for prime n (no factor to find)
  - Best for factors up to ~25 digits
  - For larger factors, need quadratic sieve or NFS
`);
