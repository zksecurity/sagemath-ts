/**
 * # Computing Euler's Totient Function
 *
 * This tutorial covers efficient algorithms for computing phi(n).
 *
 * ## Methods:
 *
 * 1. Direct counting (slow, O(n) or O(n log n))
 * 2. From prime factorization (efficient for single values)
 * 3. Sieve method (efficient for computing phi for all n up to N)
 *
 * ## Key Formulas:
 *
 * phi(p^k) = p^(k-1) * (p - 1)  for prime p
 * phi(mn) = phi(m) * phi(n)     when gcd(m, n) = 1
 *
 * @module tutorial/part2/computing-phi
 */

import {
  euler_phi,
  gcd,
  factor,
  is_prime,
  is_prime_power,
  formatFactorization,
  prime_factors,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("COMPUTING EULER'S TOTIENT FUNCTION");
console.log("=".repeat(70));

// ============================================================================
// Method 1: Direct counting (naive)
// ============================================================================
console.log("\n=== Method 1: Direct Counting ===\n");

console.log("The simplest approach: count integers coprime to n.");
console.log("Time complexity: O(n log n) using GCD.\n");

function phiNaive(n: bigint): bigint {
  if (n <= 0n) return 0n;
  let count = 0n;
  for (let k = 1n; k <= n; k++) {
    if (gcd(k, n) === 1n) {
      count++;
    }
  }
  return count;
}

// Demonstrate on small values
console.log("Computing phi(n) by counting coprimes:");
for (let n = 1n; n <= 12n; n++) {
  const naive = phiNaive(n);
  const library = euler_phi(n);
  console.log(`  phi(${n.toString().padStart(2)}) = ${naive.toString().padStart(2)} [library: ${library}]`);
}

console.log(`
This method is too slow for large n!
For n = 10^9, we'd need a billion GCD computations.
`);

// ============================================================================
// Method 2: Using prime factorization
// ============================================================================
console.log("\n=== Method 2: Using Prime Factorization ===\n");

console.log("Formula: phi(n) = n * product of (1 - 1/p) for each prime p | n");
console.log("Equivalently: phi(n) = product of p^(e-1) * (p-1) for n = product of p^e\n");

function phiFromFactors(n: bigint): { phi: bigint; steps: string[] } {
  if (n <= 0n) return { phi: 0n, steps: [] };
  if (n === 1n) return { phi: 1n, steps: ["phi(1) = 1"] };

  const factors = factor(n);
  const steps: string[] = [];
  let phi = 1n;

  steps.push(`n = ${n} = ${formatFactorization(factors)}`);

  for (const [p, e] of factors) {
    if (p === -1n) continue;
    const contribution = p ** (e - 1n) * (p - 1n);
    phi *= contribution;
    steps.push(`  phi(${p}^${e}) = ${p}^${e - 1n} * (${p} - 1) = ${contribution}`);
  }

  steps.push(`  Total: phi(${n}) = ${phi}`);

  return { phi, steps };
}

const examples = [12n, 30n, 100n, 360n, 1000n];

for (const n of examples) {
  const { phi, steps } = phiFromFactors(n);
  console.log(steps.join("\n"));
  console.log(`  Verification: ${euler_phi(n)}\n`);
}

// ============================================================================
// Method 3: Product formula implementation
// ============================================================================
console.log("\n=== Method 3: Product Formula ===\n");

console.log("An alternative formula using only prime factors (not their powers):");
console.log("  phi(n) = n * product of (1 - 1/p) for each distinct prime p | n");
console.log("  phi(n) = n * product of (p - 1) / p\n");

function phiProductFormula(n: bigint): bigint {
  if (n <= 0n) return 0n;
  if (n === 1n) return 1n;

  let result = n;
  let temp = n;

  // Find all prime factors
  let p = 2n;
  while (p * p <= temp) {
    if (temp % p === 0n) {
      // Remove all factors of p
      while (temp % p === 0n) {
        temp /= p;
      }
      // Apply (1 - 1/p) = (p-1)/p
      result -= result / p;
    }
    p++;
  }

  // If temp > 1, it's a prime factor
  if (temp > 1n) {
    result -= result / temp;
  }

  return result;
}

console.log("Demonstrating product formula:");
for (const n of [12n, 30n, 100n]) {
  const primes = prime_factors(n);
  const phi = phiProductFormula(n);

  console.log(`n = ${n}:`);
  console.log(`  Prime factors: {${primes.join(", ")}}`);
  console.log(`  phi(${n}) = ${n} * ${primes.map(p => `(1 - 1/${p})`).join(" * ")}`);
  console.log(`         = ${n} * ${primes.map(p => `(${p - 1n}/${p})`).join(" * ")}`);
  console.log(`         = ${phi}`);
  console.log(`  Verification: ${euler_phi(n)}\n`);
}

// ============================================================================
// Method 4: Sieve for multiple values
// ============================================================================
console.log("\n=== Method 4: Euler's Totient Sieve ===\n");

console.log("To compute phi(1), phi(2), ..., phi(N) efficiently,");
console.log("use a sieve similar to the Sieve of Eratosthenes.\n");

function eulerPhiSieve(N: number): bigint[] {
  // Initialize phi[i] = i
  const phi: bigint[] = new Array(N + 1);
  for (let i = 0; i <= N; i++) {
    phi[i] = BigInt(i);
  }

  // Sieve
  for (let p = 2; p <= N; p++) {
    if (phi[p] === BigInt(p)) {
      // p is prime
      for (let m = p; m <= N; m += p) {
        // Apply (1 - 1/p) factor
        phi[m] = phi[m]! - phi[m]! / BigInt(p);
      }
    }
  }

  return phi;
}

const N = 30;
const sievedPhi = eulerPhiSieve(N);

console.log(`Sieved phi values for n = 1 to ${N}:`);
console.log("-".repeat(60));

// Display in rows
for (let row = 0; row < 3; row++) {
  const start = row * 10 + 1;
  const end = Math.min(start + 9, N);
  let line1 = "n:      ";
  let line2 = "phi(n): ";
  for (let n = start; n <= end; n++) {
    line1 += n.toString().padStart(4);
    line2 += sievedPhi[n]!.toString().padStart(4);
  }
  console.log(line1);
  console.log(line2);
  console.log();
}

console.log("Time complexity: O(N log log N) - same as Sieve of Eratosthenes");

// ============================================================================
// Method 5: Special cases
// ============================================================================
console.log("\n=== Method 5: Special Cases ===\n");

console.log("Some values of phi have closed-form expressions:\n");

console.log("1. phi(p) = p - 1 for prime p");
for (const p of [2n, 3n, 5n, 7n, 11n]) {
  console.log(`   phi(${p.toString().padStart(2)}) = ${euler_phi(p)} = ${p} - 1`);
}

console.log("\n2. phi(p^k) = p^(k-1) * (p - 1) for prime power p^k");
for (const [p, k] of [[2n, 4n], [3n, 3n], [5n, 2n]] as [bigint, bigint][]) {
  const pk = p ** k;
  console.log(`   phi(${p}^${k}) = phi(${pk.toString().padStart(3)}) = ${euler_phi(pk)} = ${p}^${k - 1n} * ${p - 1n}`);
}

console.log("\n3. phi(2n) = phi(n) when n is odd (and > 1)");
for (const n of [3n, 5n, 7n, 9n, 15n]) {
  const phi2n = euler_phi(2n * n);
  const phin = euler_phi(n);
  console.log(`   phi(${2n * n}) = phi(2 * ${n}) = ${phi2n}, phi(${n}) = ${phin}`);
}

console.log("\n4. phi(2n) = 2 * phi(n) when n is even");
for (const n of [4n, 6n, 8n, 10n]) {
  const phi2n = euler_phi(2n * n);
  const phin = euler_phi(n);
  console.log(`   phi(${2n * n}) = phi(2 * ${n}) = ${phi2n}, 2 * phi(${n}) = ${2n * phin}`);
}

// ============================================================================
// Method 6: Recursive formula
// ============================================================================
console.log("\n=== Method 6: Recursive Formulas ===\n");

console.log("Using the multiplicative property recursively:\n");

function phiRecursive(n: bigint, memo: Map<bigint, bigint> = new Map()): bigint {
  if (n <= 0n) return 0n;
  if (n === 1n) return 1n;

  if (memo.has(n)) return memo.get(n)!;

  // Find smallest prime factor
  let p = 2n;
  while (p * p <= n && n % p !== 0n) {
    p++;
  }

  if (p * p > n) {
    // n is prime
    memo.set(n, n - 1n);
    return n - 1n;
  }

  // Find p^k exactly dividing n
  let pk = p;
  while (n % (pk * p) === 0n) {
    pk *= p;
  }

  const m = n / pk;

  // phi(n) = phi(pk) * phi(m) since gcd(pk, m) = 1
  const phiPk = pk - pk / p; // pk - p^(k-1) = p^(k-1) * (p-1)
  const phiM = phiRecursive(m, memo);

  const result = phiPk * phiM;
  memo.set(n, result);
  return result;
}

console.log("Computing phi using recursive decomposition:");
for (const n of [60n, 120n, 360n]) {
  const factors = factor(n);
  console.log(`phi(${n}) = phi(${formatFactorization(factors)})`);
  console.log(`        = ${phiRecursive(n)}`);
  console.log(`        Verification: ${euler_phi(n)}\n`);
}

// ============================================================================
// Performance comparison
// ============================================================================
console.log("\n=== Performance Comparison ===\n");

console.log("Comparing methods for various inputs:\n");

const testValues = [100n, 1000n, 10000n, 100000n];

console.log("     n     | Naive | Factorization | Library");
console.log("-----------+-------+---------------+--------");

for (const n of testValues) {
  // Naive only for small n
  const naive = n <= 1000n ? phiNaive(n) : "-";
  const fromFactors = phiFromFactors(n).phi;
  const library = euler_phi(n);

  console.log(`${n.toString().padStart(10)} | ${naive.toString().padStart(5)} | ${fromFactors.toString().padStart(13)} | ${library.toString().padStart(6)}`);
}

console.log(`
Performance summary:
- Naive: O(n log n), only practical for n < 10^6
- Factorization: O(sqrt(n)) for trial division, faster with advanced factoring
- Sieve: O(N log log N) for all phi values up to N
- Library: Uses factorization with optimized primality tests
`);

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
Methods for Computing phi(n):

1. Direct Counting: O(n log n)
   - Count k in [1,n] with gcd(k,n) = 1
   - Simple but slow

2. From Factorization: O(sqrt(n)) with trial division
   - phi(n) = product of p^(e-1) * (p-1)
   - Most practical for single values

3. Product Formula: O(sqrt(n))
   - phi(n) = n * product of (1 - 1/p)
   - Slightly simpler implementation

4. Sieve: O(N log log N)
   - Computes phi(1), ..., phi(N)
   - Best for computing many values

5. Recursive: O(factorization time)
   - Uses multiplicative property
   - Natural with memoization

Special Cases:
- phi(p) = p - 1
- phi(p^k) = p^(k-1) * (p - 1)
- phi(2n) = phi(n) for odd n > 1
`);
