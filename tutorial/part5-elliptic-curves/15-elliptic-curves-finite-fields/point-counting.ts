/**
 * # Point Counting on Elliptic Curves
 *
 * Determining the number of points on an elliptic curve is fundamental
 * to elliptic curve cryptography. The curve order affects:
 * - Security (larger order = harder ECDLP)
 * - Protocol correctness (need to know group structure)
 * - Cofactor handling (order = h * n where n is prime)
 *
 * ## Hasse's Theorem
 *
 * For an elliptic curve E over F_q:
 *
 *   |#E(F_q) - (q + 1)| <= 2*sqrt(q)
 *
 * This means the order is roughly q, within a "square root error".
 *
 * ## The Trace of Frobenius
 *
 * Define t = q + 1 - #E(F_q). This "trace" satisfies |t| <= 2*sqrt(q).
 * The trace encodes all information about the curve order.
 *
 * @module tutorial/elliptic-curves/point-counting
 */

import { EllipticCurve, is_supersingular, FiniteFieldPrime } from 'sagemath-ts';

// ============================================================================
// The Hasse Bound
// ============================================================================

console.log("=== Hasse's Theorem ===\n");

/**
 * Hasse's Theorem (1933):
 *
 * For E/F_q, let N = #E(F_q). Then:
 *   (sqrt(q) - 1)^2 <= N <= (sqrt(q) + 1)^2
 *
 * Or equivalently: |N - (q + 1)| <= 2*sqrt(q)
 *
 * The ratio t/sqrt(q) (where t is the trace) is bounded between -2 and 2.
 */

const F101 = new FiniteFieldPrime(101n);
const q = 101n;

// Compute bounds
const sqrtQ = Math.sqrt(Number(q));
const lowerBound = Math.ceil((sqrtQ - 1) ** 2);
const upperBound = Math.floor((sqrtQ + 1) ** 2);

console.log(`For curves over F_${q}:`);
console.log(`  sqrt(q) = ${sqrtQ.toFixed(2)}`);
console.log(`  Hasse bound: ${lowerBound} <= #E <= ${upperBound}`);
console.log(`  Range size: ${upperBound - lowerBound + 1} possible orders`);
console.log("");

// Verify with several curves
console.log("Verifying Hasse bound for various curves:");

const curveData: { a: bigint; b: bigint; order: bigint; trace: bigint }[] = [];

for (let a = 1n; a <= 5n; a++) {
  for (let b = 1n; b <= 5n; b++) {
    try {
      const E = EllipticCurve(F101, [a, b]);
      const order = E.cardinality();
      const trace = q + 1n - order;
      curveData.push({ a, b, order, trace });
    } catch {
      // Singular
    }
  }
}

// Show a sample
for (const { a, b, order, trace } of curveData.slice(0, 8)) {
  const inBound = order >= BigInt(lowerBound) && order <= BigInt(upperBound);
  console.log(`  y^2 = x^3 + ${a}x + ${b}: #E = ${order}, t = ${trace}, valid: ${inBound}`);
}
console.log("");

// ============================================================================
// Trace Distribution
// ============================================================================

console.log("=== Trace Distribution ===\n");

/**
 * The Sato-Tate conjecture (now theorem) describes how traces
 * are distributed as we vary over curves:
 *
 * For random curves over F_p, the normalized trace t/(2*sqrt(p))
 * follows the semicircle distribution: f(x) = (2/pi)*sqrt(1 - x^2)
 *
 * This has important implications for:
 * - Expected security of random curves
 * - Statistical tests for curve generation
 */

const traceDistribution = new Map<bigint, number>();

// Count traces over F_101
for (let a = 0n; a < 101n; a++) {
  for (let b = 0n; b < 101n; b++) {
    try {
      const E = EllipticCurve(F101, [a, b]);
      const trace = E.trace_of_frobenius();
      traceDistribution.set(trace, (traceDistribution.get(trace) || 0) + 1);
    } catch {
      // Singular
    }
  }
}

// Show distribution
const maxTrace = 2n * BigInt(Math.ceil(sqrtQ));
console.log(`Trace distribution for curves over F_101 (|t| <= ${maxTrace}):`);

const sortedTraces = Array.from(traceDistribution.entries())
  .sort((a, b) => (a[0] < b[0] ? -1 : 1));

for (const [trace, count] of sortedTraces) {
  const bar = '#'.repeat(Math.min(Math.floor(count / 20), 40));
  console.log(`  t = ${String(trace).padStart(4)}: ${String(count).padStart(4)} ${bar}`);
}
console.log("");

// ============================================================================
// Schoof's Algorithm
// ============================================================================

console.log("=== Schoof's Algorithm ===\n");

/**
 * Schoof's algorithm (1985) computes #E(F_q) in polynomial time.
 *
 * Key insight: Work modulo small primes l and use the Chinese Remainder Theorem.
 *
 * For each prime l:
 * 1. Find the l-torsion polynomial (division polynomial)
 * 2. Compute Frobenius action on E[l]
 * 3. Determine t mod l
 *
 * Continue until we have enough primes to determine t uniquely
 * (need product > 4*sqrt(q)).
 *
 * Improvements (SEA algorithm):
 * - Schoof-Elkies-Atkin uses isogenies for "Elkies primes"
 * - Much faster in practice: quasi-quadratic time
 *
 * The library uses PARI's implementation internally.
 */

console.log("Schoof's algorithm computes #E(F_q) in polynomial time.");
console.log("");
console.log("Algorithm outline:");
console.log("  1. For small primes l = 2, 3, 5, 7, ...");
console.log("  2. Compute trace t (mod l) using l-torsion");
console.log("  3. Use CRT to combine: t (mod L) where L = product of l's");
console.log("  4. Stop when L > 4*sqrt(q)");
console.log("");

// Demonstrate the concept
console.log("Example: Computing t for E: y^2 = x^3 + 2x + 3 over F_101");
const E_example = EllipticCurve(F101, [2n, 3n]);
const t_actual = E_example.trace_of_frobenius();

console.log(`  Actual trace: t = ${t_actual}`);
console.log(`  Hasse bound: |t| <= ${Math.ceil(2 * sqrtQ)}`);
console.log(`  Need primes with product > ${Math.ceil(4 * sqrtQ)}`);
console.log(`  Using l = 2, 3, 5, 7: product = 210 > ${Math.ceil(4 * sqrtQ)}`);
console.log("");

// ============================================================================
// Computing Cardinality in Practice
// ============================================================================

console.log("=== Computing Cardinality ===\n");

// Create a curve and compute its order
const E = EllipticCurve(F101, [1n, 1n]);

console.log(`Curve: ${E.toString()}`);

// Using the library's cardinality function
const cardinality = E.cardinality();
console.log(`  Cardinality #E(F_101) = ${cardinality}`);

// Verify by enumerating points (only feasible for small fields!)
console.log("\nVerifying by enumeration (slow, for small fields only):");

let count = 1n; // Start with 1 for the point at infinity
for (let x = 0n; x < 101n; x++) {
  const xElem = F101.__call__(x);
  const rhs = xElem.pow(3).add(F101.__call__(1n).mul(xElem)).add(F101.__call__(1n));

  if (rhs.isZero()) {
    count += 1n; // One point with y = 0
  } else if (rhs.is_square()) {
    count += 2n; // Two points with y = +/- sqrt(rhs)
  }
}

console.log(`  Enumerated count: ${count}`);
console.log(`  Matches cardinality? ${count === cardinality}`);
console.log("");

// ============================================================================
// Special Cases
// ============================================================================

console.log("=== Special Cases ===\n");

/**
 * Some curves have special structure that makes counting easier:
 *
 * 1. Supersingular curves: t = 0 (mod p), so #E = p + 1 for p > 3
 * 2. Anomalous curves: #E = p (trace t = 1)
 * 3. Curves with CM: Order determined by CM discriminant
 */

console.log("Supersingular curves (t = 0 mod p):");

// Find a supersingular curve over F_101
for (let a = 0n; a < 101n; a++) {
  for (let b = 0n; b < 101n; b++) {
    try {
      const E_ss = EllipticCurve(F101, [a, b]);
      if (is_supersingular(E_ss)) {
        const order = E_ss.cardinality();
        console.log(`  y^2 = x^3 + ${a}x + ${b}: #E = ${order} = p + 1 = 102`);
        break;
      }
    } catch {
      // Singular
    }
  }
  if (is_supersingular(EllipticCurve(F101, [0n, 1n]))) break;
}

console.log("");
console.log("Anomalous curves (#E = p, INSECURE):");

// Search for anomalous curve
for (let a = 0n; a < 101n; a++) {
  for (let b = 0n; b < 101n; b++) {
    try {
      const E_anom = EllipticCurve(F101, [a, b]);
      if (E_anom.cardinality() === 101n) {
        console.log(`  y^2 = x^3 + ${a}x + ${b}: #E = ${E_anom.cardinality()} = p`);
        console.log(`  WARNING: Smart-ASS attack breaks ECDLP!`);
        break;
      }
    } catch {
      // Singular
    }
  }
}
console.log("");

// ============================================================================
// Counting Over Extension Fields
// ============================================================================

console.log("=== Counting Over Extension Fields ===\n");

/**
 * Given #E(F_q), we can compute #E(F_{q^n}) using the trace.
 *
 * The recurrence:
 *   #E(F_{q^n}) = q^n + 1 - alpha^n - beta^n
 *
 * where alpha, beta are roots of X^2 - tX + q = 0.
 *
 * Equivalently, there's a linear recurrence:
 *   a_n = t * a_{n-1} - q * a_{n-2}
 * where a_n = alpha^n + beta^n = q^n + 1 - #E(F_{q^n})
 */

const trace = E.trace_of_frobenius();

console.log(`For E: y^2 = x^3 + x + 1 over F_101:`);
console.log(`  Trace t = ${trace}`);
console.log("");

// Compute counts over extensions using the recurrence
// a_n = alpha^n + beta^n, where a_0 = 2, a_1 = t
let a_prev2 = 2n;  // a_0
let a_prev1 = trace;  // a_1

console.log("Points over extension fields:");
console.log(`  #E(F_{101^1}) = ${q + 1n - a_prev1}`);

for (let n = 2; n <= 5; n++) {
  const a_n = trace * a_prev1 - q * a_prev2;
  const q_n = q ** BigInt(n);
  const count_n = q_n + 1n - a_n;

  console.log(`  #E(F_{101^${n}}) = ${count_n}`);

  a_prev2 = a_prev1;
  a_prev1 = a_n;
}
console.log("");

// ============================================================================
// Importance for Cryptography
// ============================================================================

console.log("=== Cryptographic Importance ===\n");

console.log("Point counting is essential for:");
console.log("");
console.log("1. CURVE VALIDATION");
console.log("   - Verify claimed order is correct");
console.log("   - Check for weak curves (anomalous, low embedding degree)");
console.log("");
console.log("2. SUBGROUP SELECTION");
console.log("   - Factor #E = h * n where n is large prime");
console.log("   - h is the cofactor (should be small: 1, 2, 4, 8)");
console.log("");
console.log("3. PROTOCOL CORRECTNESS");
console.log("   - Need nP = O for all points P");
console.log("   - Need order for Pohlig-Hellman to assess security");
console.log("");

// Show factorization of curve order
console.log(`Example: E over F_101, #E = ${cardinality}`);
console.log(`  Factoring ${cardinality}...`);

// Simple trial division
let n = cardinality;
const factors: bigint[] = [];
let d = 2n;

while (d * d <= n) {
  while (n % d === 0n) {
    factors.push(d);
    n /= d;
  }
  d++;
}

if (n > 1n) {
  factors.push(n);
}

console.log(`  Factorization: ${factors.join(' * ')}`);

const largestFactor = factors[factors.length - 1] || 1n;
const cofactor = cardinality / largestFactor;

console.log(`  Largest prime factor (n): ${largestFactor}`);
console.log(`  Cofactor (h): ${cofactor}`);
console.log("");

// ============================================================================
// Exercises
// ============================================================================

console.log("=== Exercises ===\n");

console.log("1. Verify the Hasse bound for all curves over F_17");
console.log("   What is the range of possible orders?");
console.log("");

console.log("2. Find all curves over F_23 with prime order");
console.log("   (These are ideal for cryptographic use)");
console.log("");

console.log("3. Compute #E(F_{7^2}) and #E(F_{7^3}) for E: y^2 = x^3 + x + 1 over F_7");
console.log("   using the trace recurrence formula");
console.log("");
