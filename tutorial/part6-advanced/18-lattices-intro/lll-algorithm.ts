/**
 * # The LLL Algorithm
 *
 * The Lenstra-Lenstra-Lovasz (LLL) algorithm is one of the most important
 * algorithms in computational mathematics. Published in 1982, it efficiently
 * finds "short" vectors in lattices and has applications far beyond cryptography.
 *
 * ## What LLL Does
 *
 * Given a lattice basis, LLL produces a "reduced" basis where:
 * - Vectors are relatively short
 * - Vectors are nearly orthogonal
 * - First vector approximates the shortest vector (within 2^{n/2})
 *
 * ## Why LLL Matters
 *
 * - BREAKING CRYPTO: Coppersmith's attack on RSA uses LLL
 * - BUILDING CRYPTO: LLL is used in lattice scheme implementations
 * - DIOPHANTINE: Solves approximate integer relation problems
 * - POLYNOMIAL FACTORING: Original motivation for LLL
 *
 * @module tutorial/part6-advanced/18-lattices-intro/lll-algorithm
 */

import { IntegerLattice, IntegerMatrixFromEntries, lllReduce, gramSchmidt, isLLLReduced } from 'sagemath-ts';

console.log("=== The LLL Algorithm ===\n");

// ============================================================================
// Section 1: What is LLL Reduction?
// ============================================================================

console.log("--- Section 1: LLL Reduction ---\n");

/**
 * LLL-REDUCED BASIS
 *
 * A basis B = (b_1, ..., b_n) is LLL-reduced with parameter delta if:
 *
 * 1. SIZE REDUCTION:
 *    |mu_{i,j}| <= 1/2 for all i > j
 *
 *    where mu_{i,j} = <b_i, b_j*> / <b_j*, b_j*>
 *
 *    Meaning: each vector is "almost orthogonal" to previous GS vectors
 *
 * 2. LOVASZ CONDITION:
 *    For all i = 1, ..., n-1:
 *    delta * ||b_i*||^2 <= ||b_{i+1}* + mu_{i+1,i} * b_i*||^2
 *
 *    Simplifies to:
 *    delta * ||b_i*||^2 <= ||b_{i+1}*||^2 + mu_{i+1,i}^2 * ||b_i*||^2
 *
 *    Or:
 *    (delta - mu_{i+1,i}^2) * ||b_i*||^2 <= ||b_{i+1}*||^2
 *
 * PARAMETER DELTA:
 * - Must satisfy: 1/4 < delta <= 1
 * - Standard choice: delta = 3/4 or delta = 0.99
 * - Larger delta = better basis but slower algorithm
 *
 * GUARANTEE:
 * For an LLL-reduced basis with delta = 3/4:
 *   ||b_1|| <= 2^{(n-1)/2} * lambda_1(L)
 *
 * The first vector is within 2^{n/2} of the shortest!
 */

console.log("LLL-Reduced Basis Definition:\n");
console.log("1. Size Reduction: |mu_{i,j}| <= 1/2 for all i > j");
console.log("   (Each vector nearly orthogonal to previous)\n");
console.log("2. Lovasz Condition: delta * ||b_i*||^2 <= ||b_{i+1}*||^2 + mu^2 * ||b_i*||^2");
console.log("   (GS vectors don't decrease too quickly)\n");

console.log("Approximation Guarantee:");
console.log("  ||b_1|| <= 2^{(n-1)/2} * lambda_1(L)");
console.log("  First vector within factor 2^{n/2} of shortest\n");

// ============================================================================
// Section 2: The Algorithm
// ============================================================================

console.log("--- Section 2: LLL Algorithm ---\n");

/**
 * LLL ALGORITHM
 *
 * Input:  Basis B = (b_1, ..., b_n), parameter delta in (1/4, 1]
 * Output: LLL-reduced basis B' for the same lattice
 *
 * Algorithm:
 *
 * 1. Compute Gram-Schmidt orthogonalization B*, mu
 *
 * 2. k = 2
 *
 * 3. While k <= n:
 *    a. SIZE REDUCE b_k:
 *       For j = k-1, k-2, ..., 1:
 *         If |mu_{k,j}| > 1/2:
 *           b_k = b_k - round(mu_{k,j}) * b_j
 *           Update mu values
 *
 *    b. CHECK LOVASZ:
 *       If delta * ||b_{k-1}*||^2 > ||b_k*||^2 + mu_{k,k-1}^2 * ||b_{k-1}*||^2:
 *         SWAP b_k and b_{k-1}
 *         Update Gram-Schmidt
 *         k = max(k-1, 2)
 *       Else:
 *         k = k + 1
 *
 * 4. Return B
 *
 * COMPLEXITY:
 * - Number of swaps: O(n^2 log B) where B = max ||b_i||
 * - Each iteration: O(n^2) arithmetic operations
 * - Total: O(n^5 log^3 B) bit operations
 *
 * In practice, LLL is very fast - usually runs in near-linear time.
 */

console.log("LLL Algorithm (simplified):\n");
console.log("  1. Compute Gram-Schmidt");
console.log("  2. For k = 2 to n:");
console.log("       a. Size reduce b_k (make |mu| <= 1/2)");
console.log("       b. If Lovasz fails:");
console.log("            Swap b_k and b_{k-1}");
console.log("            Go back (k = k-1)");
console.log("          Else:");
console.log("            Move forward (k = k+1)");
console.log("  3. Return reduced basis\n");

console.log("Complexity:");
console.log("  O(n^5 log^3 B) bit operations");
console.log("  In practice: very fast, often near-linear\n");

// ============================================================================
// Section 3: Running LLL
// ============================================================================

console.log("--- Section 3: LLL in Action ---\n");

// Create a "bad" basis
const badBasis = IntegerMatrixFromEntries([
  [1n, 1n, 1n],
  [-1n, 0n, 2n],
  [3n, 5n, 6n]
]);

console.log("Original (bad) basis:");
for (let i = 0; i < badBasis.nrows; i++) {
  const row = [];
  let normSq = 0n;
  for (let j = 0; j < badBasis.ncols; j++) {
    const v = badBasis.get(i, j).value;
    row.push(v);
    normSq += v * v;
  }
  console.log(`  b${i+1} = (${row.join(', ')})  ||b${i+1}|| = ${Math.sqrt(Number(normSq)).toFixed(3)}`);
}

// Check GS before reduction
const gsBefore = gramSchmidt(badBasis);
console.log("\nGram-Schmidt before reduction:");
console.log(`  ||b1*|| = ${Math.sqrt(gsBefore.B[0]!).toFixed(3)}`);
console.log(`  ||b2*|| = ${Math.sqrt(gsBefore.B[1]!).toFixed(3)}`);
console.log(`  ||b3*|| = ${Math.sqrt(gsBefore.B[2]!).toFixed(3)}`);
console.log(`  mu_{2,1} = ${gsBefore.mu[1]![0]!.toFixed(3)}`);
console.log(`  mu_{3,1} = ${gsBefore.mu[2]![0]!.toFixed(3)}`);
console.log(`  mu_{3,2} = ${gsBefore.mu[2]![1]!.toFixed(3)}`);

// Apply LLL
console.log("\n--- Running LLL ---\n");
const reducedBasis = lllReduce(badBasis, { delta: 0.75, eta: 0.501 });

console.log("Reduced (good) basis:");
for (let i = 0; i < reducedBasis.nrows; i++) {
  const row = [];
  let normSq = 0n;
  for (let j = 0; j < reducedBasis.ncols; j++) {
    const v = reducedBasis.get(i, j).value;
    row.push(v);
    normSq += v * v;
  }
  console.log(`  b${i+1}' = (${row.join(', ')})  ||b${i+1}'|| = ${Math.sqrt(Number(normSq)).toFixed(3)}`);
}

// Check GS after reduction
const gsAfter = gramSchmidt(reducedBasis);
console.log("\nGram-Schmidt after reduction:");
console.log(`  ||b1*'|| = ${Math.sqrt(gsAfter.B[0]!).toFixed(3)}`);
console.log(`  ||b2*'|| = ${Math.sqrt(gsAfter.B[1]!).toFixed(3)}`);
console.log(`  ||b3*'|| = ${Math.sqrt(gsAfter.B[2]!).toFixed(3)}`);
console.log(`  mu_{2,1} = ${gsAfter.mu[1]![0]!.toFixed(3)}`);
console.log(`  mu_{3,1} = ${gsAfter.mu[2]![0]!.toFixed(3)}`);
console.log(`  mu_{3,2} = ${gsAfter.mu[2]![1]!.toFixed(3)}`);

// Verify LLL conditions
console.log("\nVerifying LLL conditions:");
console.log(`  |mu_{2,1}| <= 0.5? ${Math.abs(gsAfter.mu[1]![0]!) <= 0.501}`);
console.log(`  |mu_{3,1}| <= 0.5? ${Math.abs(gsAfter.mu[2]![0]!) <= 0.501}`);
console.log(`  |mu_{3,2}| <= 0.5? ${Math.abs(gsAfter.mu[2]![1]!) <= 0.501}`);

const delta = 0.75;
const lovasz1 = delta * gsAfter.B[0]! <= gsAfter.B[1]! + gsAfter.mu[1]![0]! ** 2 * gsAfter.B[0]!;
const lovasz2 = delta * gsAfter.B[1]! <= gsAfter.B[2]! + gsAfter.mu[2]![1]! ** 2 * gsAfter.B[1]!;
console.log(`  Lovasz(1,2) satisfied? ${lovasz1}`);
console.log(`  Lovasz(2,3) satisfied? ${lovasz2}\n`);

// ============================================================================
// Section 4: Applications of LLL
// ============================================================================

console.log("--- Section 4: Applications of LLL ---\n");

/**
 * APPLICATIONS OF LLL
 *
 * 1. POLYNOMIAL FACTORING (original application):
 *    - Factor polynomials over Q in polynomial time
 *    - Groundbreaking when published in 1982
 *
 * 2. INTEGER RELATIONS:
 *    - Given reals x_1, ..., x_n, find integers a_i with
 *      a_1*x_1 + ... + a_n*x_n = 0
 *    - Used to discover identities in computer algebra
 *
 * 3. DIOPHANTINE APPROXIMATION:
 *    - Find good rational approximations
 *    - Solve systems of linear Diophantine equations
 *
 * 4. CRYPTANALYSIS:
 *    - Coppersmith's attack: Find small roots of polynomials mod N
 *    - Break weak RSA keys
 *    - Attack knapsack cryptosystems
 *    - Solve hidden number problems
 *
 * 5. LATTICE CRYPTOGRAPHY:
 *    - Basis reduction in key generation
 *    - Decryption in some schemes
 *    - Security analysis
 *
 * 6. SIGNAL PROCESSING:
 *    - MIMO detection
 *    - GPS signal processing
 *
 * 7. LEARNING THEORY:
 *    - Solve subset-sum problems
 *    - Learning parity with noise
 */

console.log("LLL Applications:\n");
console.log("  1. Polynomial factoring (original purpose)");
console.log("  2. Integer relation finding");
console.log("  3. Diophantine approximation");
console.log("  4. Cryptanalysis (Coppersmith, knapsack attacks)");
console.log("  5. Lattice-based cryptography");
console.log("  6. Signal processing (MIMO detection)");
console.log("  7. Machine learning\n");

// ============================================================================
// Section 5: BKZ: Beyond LLL
// ============================================================================

console.log("--- Section 5: BKZ Algorithm ---\n");

/**
 * BKZ (Block Korkine-Zolotareff) ALGORITHM
 *
 * BKZ is a generalization of LLL that achieves better reduction
 * at the cost of more computation time.
 *
 * IDEA:
 * Instead of just comparing adjacent pairs (b_i, b_{i+1}),
 * BKZ considers blocks of size beta and finds the shortest
 * vector within each block using enumeration.
 *
 * PARAMETERS:
 * - Block size beta: larger = better reduction, slower
 * - BKZ-2 = LLL (special case)
 * - BKZ-n = HKZ (Hermite-Korkine-Zolotareff, optimal but slow)
 *
 * APPROXIMATION:
 * BKZ-beta achieves approximation factor roughly:
 *   gamma = (beta)^{n/(2*beta)}
 *
 * For beta = 2: gamma = 2^{n/2} (same as LLL)
 * For beta = n: gamma = 1 (exact shortest vector)
 *
 * COMPLEXITY:
 * Each block requires SVP in dimension beta:
 * - Using enumeration: 2^{O(beta^2)} per block
 * - Using sieving: 2^{O(beta)} per block
 *
 * Total time: roughly 2^{O(beta)} * poly(n)
 *
 * BKZ-2.0 (improved version):
 * - Uses pruned enumeration
 * - Incorporates preprocessing
 * - Faster in practice
 */

console.log("BKZ (Block Korkine-Zolotareff):\n");
console.log("  Generalization of LLL with block parameter beta");
console.log("  Finds SVP in blocks of size beta\n");

console.log("Trade-off:");
console.log("  beta = 2:   Same as LLL, gamma = 2^{n/2}");
console.log("  beta = 20:  gamma ~ 1.01^n, time ~2^20");
console.log("  beta = 60:  gamma ~ 1.003^n, time ~2^60");
console.log("  beta = n:   gamma = 1 (exact), time 2^{O(n)}\n");

console.log("Used in cryptanalysis:");
console.log("  Attackers run BKZ with increasing beta");
console.log("  Security parameter chosen so required beta is infeasible\n");

// ============================================================================
// Section 6: LLL Variants
// ============================================================================

console.log("--- Section 6: LLL Variants ---\n");

/**
 * LLL VARIANTS AND IMPROVEMENTS
 *
 * DEEP LLL (L^3 algorithm):
 * - Check Lovasz condition with all previous vectors, not just adjacent
 * - Better output quality
 * - Still polynomial time
 *
 * FLOATING-POINT LLL:
 * - Use floating-point arithmetic for speed
 * - Carefully track precision to ensure correctness
 * - Much faster in practice
 *
 * PROGRESSIVE BKZ:
 * - Start with small beta, increase gradually
 * - Often finds good vectors faster than fixed beta
 *
 * SLIDE REDUCTION:
 * - Variant that gives comparable results to BKZ
 * - Different block structure
 *
 * LLL FOR SPECIFIC LATTICES:
 * - Ideal lattices: Use ring structure for speedup
 * - q-ary lattices: Exploit modular structure
 *
 * IMPLEMENTATIONS:
 * - fpLLL: Reference implementation, very fast
 * - NTL: Number Theory Library, good for large integers
 * - FLINT: Fast Library for Number Theory
 * - Our sagemath-ts: Educational implementation
 */

console.log("LLL Variants:\n");
console.log("  Deep LLL:        Better quality, still polynomial");
console.log("  Floating LLL:    Faster using FP arithmetic");
console.log("  Progressive BKZ: Gradually increase block size\n");

console.log("Implementations:");
console.log("  fpLLL:      Reference C++ implementation");
console.log("  NTL:        Number Theory Library");
console.log("  FLINT:      Fast Library for Number Theory");
console.log("  sagemath-ts: Educational TypeScript implementation\n");

// ============================================================================
// Section 7: Practical Example - Coppersmith
// ============================================================================

console.log("--- Section 7: Example - Small Root Finding ---\n");

/**
 * COPPERSMITH'S METHOD (simplified)
 *
 * Problem: Given polynomial f(x) and modulus N,
 *          find small x with f(x) = 0 mod N
 *
 * This is useful for:
 * - Attacking RSA with small private exponent
 * - Factoring N when partial information is known
 * - Breaking many cryptographic schemes
 *
 * IDEA:
 * 1. Construct a lattice from f(x) and its shifts
 * 2. Apply LLL to find short vectors
 * 3. Short vectors correspond to polynomials with small roots
 *
 * EXAMPLE CONSTRUCTION:
 * For f(x) = x^2 + ax + b mod N and bound X:
 *
 * Consider polynomials: f(x), x*f(x), N, N*x
 * Write as vectors of coefficients scaled by powers of X
 *
 *           x^0    x^1    x^2    x^3
 * f(x):    [ b,     a*X,   X^2,   0   ]
 * x*f(x):  [ 0,     b*X,   a*X^2, X^3 ]
 * N:       [ N,     0,     0,     0   ]
 * N*x:     [ 0,     N*X,   0,     0   ]
 *
 * LLL finds a short vector, corresponding to a polynomial g(x)
 * with g(x_0) = 0 over the integers (not just mod N!)
 *
 * Then solve g(x) = 0 to find the small root x_0.
 */

console.log("Coppersmith's Method (simplified):\n");
console.log("  Goal: Find small x with f(x) = 0 mod N\n");
console.log("  Steps:");
console.log("    1. Build lattice from polynomial and shifts");
console.log("    2. Apply LLL to find short vector");
console.log("    3. Short vector = polynomial with integer root");
console.log("    4. Solve polynomial to find x\n");

// Simple example
console.log("Example lattice for Coppersmith:");
console.log("  Let f(x) = x - 42 mod 100, bound X = 50");
console.log("  Lattice basis:");
console.log("    [ 1,   50  ]   (coefficients of f scaled by X)");
console.log("    [ 0,   100 ]   (modular relation)");

const coppBasis = IntegerMatrixFromEntries([
  [1n, 50n],
  [0n, 100n]
]);

const coppReduced = lllReduce(coppBasis);
console.log("  After LLL:");
console.log(`    [ ${coppReduced.get(0, 0).value}, ${coppReduced.get(0, 1).value} ]`);
console.log(`    [ ${coppReduced.get(1, 0).value}, ${coppReduced.get(1, 1).value} ]`);
console.log("  Short vector might reveal the root!\n");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");

console.log("LLL Algorithm:");
console.log("  Produces a 'reduced' basis in polynomial time");
console.log("  Guarantees: vectors short, nearly orthogonal\n");

console.log("Key properties:");
console.log("  First vector: within 2^{n/2} of shortest");
console.log("  Complexity: O(n^5 log^3 B) but fast in practice\n");

console.log("BKZ extends LLL:");
console.log("  Block size beta trades time for quality");
console.log("  Used in attacks to get better approximations\n");

console.log("Applications:");
console.log("  - Cryptanalysis (Coppersmith, RSA attacks)");
console.log("  - Lattice-based cryptography");
console.log("  - Integer relation finding");
console.log("  - Polynomial factoring\n");

console.log("Next: Lattice-based cryptography (NTRU, Kyber concepts)");
