/**
 * # SVP and CVP: The Core Lattice Problems
 *
 * The Shortest Vector Problem (SVP) and Closest Vector Problem (CVP) are
 * the foundational hard problems in lattice cryptography. Understanding
 * these problems is essential for understanding why lattice crypto works.
 *
 * ## Why These Problems Matter
 *
 * 1. SECURITY FOUNDATION: Lattice crypto assumes these are hard
 * 2. QUANTUM RESISTANCE: No known quantum speedup for SVP/CVP
 * 3. PRACTICAL ATTACKS: Solving SVP/CVP breaks lattice schemes
 *
 * @module tutorial/part6-advanced/18-lattices-intro/svp-cvp
 */

import { IntegerLattice, IntegerMatrixFromEntries, lllReduce, gramSchmidt } from 'sagemath-ts';

console.log("=== SVP and CVP: Core Lattice Problems ===\n");

// ============================================================================
// Section 1: The Shortest Vector Problem (SVP)
// ============================================================================

console.log("--- Section 1: Shortest Vector Problem (SVP) ---\n");

/**
 * SHORTEST VECTOR PROBLEM (SVP)
 *
 * EXACT SVP:
 * Given: A lattice L (specified by basis B)
 * Find:  A shortest nonzero vector v in L
 *
 * That is, find v in L with ||v|| = lambda_1(L)
 *
 * APPROXIMATE SVP (gamma-SVP):
 * Given: A lattice L
 * Find:  A nonzero vector v in L with ||v|| <= gamma * lambda_1(L)
 *
 * gamma = 1: exact SVP
 * gamma > 1: allows approximation factor gamma
 *
 * HARDNESS:
 *
 * Exact SVP:
 * - NP-hard under randomized reductions (Ajtai 1998)
 * - Best known algorithms: 2^O(n) time
 *
 * gamma-SVP:
 * - gamma = poly(n): NP-hard (assuming NP != coNP)
 * - gamma = 2^sqrt(n): can be solved in poly time (LLL-based)
 * - gamma = 2^n/log(n): best quantum algorithms
 *
 * The hardness of gamma-SVP for "medium" gamma (like 2^sqrt(n))
 * is the security assumption for most lattice crypto.
 *
 * ALGORITHMS:
 *
 * 1. Enumeration: Deterministic 2^{O(n^2)} time, poly space
 * 2. Sieving: Randomized 2^{O(n)} time AND space
 * 3. LLL: Polynomial time but only gamma = 2^{n/2} approximation
 * 4. BKZ: Trade-off between time and approximation factor
 */

console.log("SVP Definition:");
console.log("  Input: Lattice L (given by basis B)");
console.log("  Output: Shortest nonzero vector in L\n");

console.log("Approximate SVP (gamma-SVP):");
console.log("  Find v with ||v|| <= gamma * lambda_1(L)");
console.log("  gamma = approximation factor\n");

console.log("Hardness:");
console.log("  Exact SVP: NP-hard");
console.log("  gamma-SVP with small gamma: also hard");
console.log("  Best algorithms: 2^O(n) time (sieving)\n");

// Create a lattice and discuss SVP
const svpBasis = IntegerMatrixFromEntries([
  [17n, 0n, 0n],
  [0n, 17n, 0n],
  [3n, 5n, 17n]
]);

const svpLattice = IntegerLattice(svpBasis);
const reducedBasis = svpLattice.reducedBasis;

console.log("Example lattice:");
console.log("  Original basis:");
for (let i = 0; i < svpBasis.nrows; i++) {
  const row = [];
  for (let j = 0; j < svpBasis.ncols; j++) {
    row.push(svpBasis.get(i, j).value);
  }
  console.log(`    (${row.join(', ')})`);
}

console.log("\n  After LLL reduction:");
for (let i = 0; i < reducedBasis.nrows; i++) {
  const row = [];
  for (let j = 0; j < reducedBasis.ncols; j++) {
    row.push(reducedBasis.get(i, j).value);
  }
  console.log(`    (${row.join(', ')})`);
}

// First vector of LLL-reduced basis approximates shortest vector
function norm(matrix: any, row: number): number {
  let sum = 0;
  for (let j = 0; j < matrix.ncols; j++) {
    const v = Number(matrix.get(row, j).value);
    sum += v * v;
  }
  return Math.sqrt(sum);
}

console.log("\n  Approximate shortest vector (LLL output):");
console.log(`    Length: ${norm(reducedBasis, 0).toFixed(3)}`);
console.log("  (May not be exact shortest, but within factor 2^{n/2})\n");

// ============================================================================
// Section 2: The Closest Vector Problem (CVP)
// ============================================================================

console.log("--- Section 2: Closest Vector Problem (CVP) ---\n");

/**
 * CLOSEST VECTOR PROBLEM (CVP)
 *
 * EXACT CVP:
 * Given: A lattice L and a target vector t in R^n
 * Find:  A lattice point v in L closest to t
 *
 * That is, find v in L minimizing ||t - v||
 *
 * APPROXIMATE CVP (gamma-CVP):
 * Given: A lattice L and target t
 * Find:  v in L with ||t - v|| <= gamma * dist(t, L)
 *
 * where dist(t, L) = min_{u in L} ||t - u||
 *
 * RELATION TO SVP:
 * - CVP is at least as hard as SVP
 * - Can reduce SVP to CVP (but not other direction in general)
 * - Both are NP-hard under randomized reductions
 *
 * ALGORITHMS:
 *
 * 1. Babai's Nearest Plane:
 *    - Use Gram-Schmidt, round to nearest lattice point
 *    - Polynomial time
 *    - Approximation: 2^{n/2} factor
 *
 * 2. Babai's Rounding:
 *    - Simpler version, slightly worse approximation
 *    - Also called "round-off" algorithm
 *
 * 3. Enumeration:
 *    - Exact solution
 *    - Exponential time: 2^{O(n)}
 *
 * 4. Embedding technique:
 *    - Reduce CVP to SVP in dimension n+1
 *    - Not always practical but theoretically nice
 */

console.log("CVP Definition:");
console.log("  Input: Lattice L, target vector t");
console.log("  Output: Lattice point v closest to t\n");

console.log("Relation to SVP:");
console.log("  CVP is at least as hard as SVP");
console.log("  Both are NP-hard\n");

console.log("Babai's Algorithm (polynomial time, approximate):");
console.log("  1. Compute Gram-Schmidt of reduced basis");
console.log("  2. Project target onto orthogonal complement");
console.log("  3. Round coefficients to nearest integer");
console.log("  4. Approximation factor: 2^{n/2}\n");

// Demonstrate CVP conceptually
console.log("CVP Example:");
console.log("  Lattice: Z^2 (integer grid)");
console.log("  Target: (2.7, 3.2)");
console.log("  Closest point: (3, 3)");
console.log("  Distance: sqrt((3-2.7)^2 + (3-3.2)^2) = 0.36\n");

// ============================================================================
// Section 3: Babai's Algorithm
// ============================================================================

console.log("--- Section 3: Babai's Algorithm ---\n");

/**
 * BABAI'S NEAREST PLANE ALGORITHM
 *
 * Given: LLL-reduced basis B = (b_1, ..., b_n), target t
 * Output: Lattice point close to t
 *
 * Algorithm:
 * 1. Compute Gram-Schmidt basis B* = (b_1*, ..., b_n*)
 * 2. Set b = t
 * 3. For i = n, n-1, ..., 1:
 *    a. c_i = round(<b, b_i*> / ||b_i*||^2)
 *    b. b = b - c_i * b_i
 * 4. Return t - b = sum_i c_i * b_i
 *
 * Why it works:
 * - Projects t onto each orthogonal direction
 * - Rounds to nearest integer in each direction
 * - LLL reduction ensures basis is "nice" enough
 *
 * APPROXIMATION GUARANTEE:
 *
 * For an LLL-reduced basis with parameter delta:
 *   ||Babai(t) - v_closest|| <= ||t - v_closest|| * (2/sqrt(4*delta - 1))^n
 *
 * With standard delta = 3/4:
 *   Approximation factor: 2^{n/2}
 *
 * Better approximation requires better basis reduction (BKZ).
 */

console.log("Babai's Nearest Plane Algorithm:\n");
console.log("  Input: LLL-reduced basis B, target t");
console.log("  Output: Lattice point close to t\n");

console.log("Algorithm:");
console.log("  1. Compute Gram-Schmidt B*");
console.log("  2. For i = n down to 1:");
console.log("       c_i = round(<t, b_i*> / ||b_i*||^2)");
console.log("       t = t - c_i * b_i");
console.log("  3. Return sum of c_i * b_i\n");

// Conceptual implementation
function babaiNearestPlane(
  basis: number[][],
  target: number[]
): { coeffs: number[], point: number[] } {
  const n = basis.length;
  const m = basis[0]!.length;

  // Gram-Schmidt
  const bStar: number[][] = [];
  const mu: number[][] = [];
  const B: number[] = [];

  for (let i = 0; i < n; i++) {
    mu.push(new Array(n).fill(0));
    mu[i]![i] = 1;

    const bi = [...basis[i]!];
    const biStar = [...bi];

    for (let j = 0; j < i; j++) {
      let dot = 0;
      for (let k = 0; k < m; k++) {
        dot += bi[k]! * bStar[j]![k]!;
      }
      mu[i]![j] = dot / B[j]!;

      for (let k = 0; k < m; k++) {
        biStar[k] = biStar[k]! - mu[i]![j]! * bStar[j]![k]!;
      }
    }

    bStar.push(biStar);
    let normSq = 0;
    for (let k = 0; k < m; k++) {
      normSq += biStar[k]! * biStar[k]!;
    }
    B.push(normSq);
  }

  // Babai's algorithm
  let b = [...target];
  const coeffs: number[] = new Array(n).fill(0);

  for (let i = n - 1; i >= 0; i--) {
    let dot = 0;
    for (let k = 0; k < m; k++) {
      dot += b[k]! * bStar[i]![k]!;
    }
    coeffs[i] = Math.round(dot / B[i]!);

    for (let k = 0; k < m; k++) {
      b[k] = b[k]! - coeffs[i]! * basis[i]![k]!;
    }
  }

  // Compute the lattice point
  const point = new Array(m).fill(0);
  for (let i = 0; i < n; i++) {
    for (let k = 0; k < m; k++) {
      point[k] = point[k]! + coeffs[i]! * basis[i]![k]!;
    }
  }

  return { coeffs, point };
}

// Demo
const demoBasis = [
  [1, 0],
  [0, 1]
];
const demoTarget = [2.7, 3.2];

const babaiResult = babaiNearestPlane(demoBasis, demoTarget);
console.log("Demo: Babai on Z^2");
console.log(`  Target: (${demoTarget.join(', ')})`);
console.log(`  Coefficients: (${babaiResult.coeffs.join(', ')})`);
console.log(`  Nearest point: (${babaiResult.point.join(', ')})\n`);

// ============================================================================
// Section 4: Computational Hardness
// ============================================================================

console.log("--- Section 4: Computational Hardness ---\n");

/**
 * COMPLEXITY OF SVP AND CVP
 *
 * EXACT PROBLEMS:
 *
 * Exact SVP:
 * - NP-hard under randomized reductions [Ajtai 1998]
 * - In NP intersect coNP (so unlikely to be NP-complete)
 * - Best deterministic: 2^{2n + o(n)} time [AKS]
 * - Best randomized: 2^{n + o(n)} time and space [ADRS, Laarhoven]
 *
 * Exact CVP:
 * - NP-hard [GMSS]
 * - At least as hard as SVP
 * - Similar algorithmic complexity
 *
 * APPROXIMATE PROBLEMS:
 *
 * gamma-SVP for gamma = 2^{n^eps}:
 * - Can be solved in poly time (LLL for eps = 0.5)
 *
 * gamma-SVP for gamma = poly(n):
 * - NP-hard assuming NP != coNP [Khot 2005]
 *
 * gamma-SVP for gamma = O(sqrt(n)):
 * - Believed hard (basis of cryptography)
 * - Not known to be NP-hard!
 *
 * QUANTUM COMPLEXITY:
 *
 * No known quantum speedup for SVP or CVP!
 * - Shor's algorithm doesn't help
 * - Grover gives at most sqrt speedup on search
 * - Best quantum: still 2^{O(n)} for exact
 *
 * This is why lattice crypto is post-quantum secure.
 */

console.log("Exact Problem Complexity:\n");
console.log("  SVP: NP-hard (randomized reduction)");
console.log("       Best algorithm: 2^n time and space\n");
console.log("  CVP: NP-hard, at least as hard as SVP");
console.log("       Similar complexity\n");

console.log("Approximate Problem Complexity:\n");
console.log("  gamma = 2^{n/2}: Polynomial (LLL)");
console.log("  gamma = poly(n): NP-hard");
console.log("  gamma = O(sqrt(n)): Believed hard (crypto assumption)\n");

console.log("Quantum Complexity:\n");
console.log("  NO known quantum speedup!");
console.log("  Still 2^{O(n)} for quantum algorithms");
console.log("  This is why lattices are POST-QUANTUM SECURE\n");

// ============================================================================
// Section 5: SVP/CVP Algorithms
// ============================================================================

console.log("--- Section 5: Algorithms Overview ---\n");

/**
 * ALGORITHM COMPARISON
 *
 * ENUMERATION (Schnorr-Euchner 1994):
 * - Deterministic
 * - Time: 2^{O(n^2)} (or 2^{O(n log n)} with preprocessing)
 * - Space: Polynomial
 * - Finds exact shortest vector
 * - Practical for n ~ 60-80
 *
 * SIEVING (AKS 2001, Laarhoven 2015):
 * - Randomized
 * - Time: 2^{O(n)}
 * - Space: 2^{O(n)} (main drawback!)
 * - Finds exact shortest vector
 * - Practical for n ~ 70-90 with huge memory
 *
 * LLL (Lenstra-Lenstra-Lovasz 1982):
 * - Polynomial time: O(n^5 log^3 B) for B-bit entries
 * - Finds 2^{n/2} approximation
 * - Practical for any n
 * - The workhorse of lattice algorithms
 *
 * BKZ (Block Korkine-Zolotareff):
 * - Combines LLL and enumeration
 * - Parameter: block size beta
 * - Time: 2^{O(beta)}
 * - Approximation: roughly beta^{n/beta}
 * - Trade-off between time and quality
 *
 * ALGORITHM SELECTION:
 *
 * For crypto attacks:
 * - Low dimension (n < 60): Enumeration
 * - Medium dimension (60 < n < 90): Sieving (if memory available)
 * - High dimension (n > 90): BKZ with increasing block size
 *
 * For crypto construction:
 * - Always use LLL for efficiency
 * - Security parameter chosen so BKZ attacks infeasible
 */

console.log("Algorithm Comparison:\n");
console.log("| Algorithm    | Time        | Space   | Approximation |");
console.log("|--------------|-------------|---------|---------------|");
console.log("| Enumeration  | 2^{n^2}     | poly    | Exact         |");
console.log("| Sieving      | 2^{0.29n}   | 2^{0.2n}| Exact         |");
console.log("| LLL          | poly        | poly    | 2^{n/2}       |");
console.log("| BKZ-beta     | 2^{O(beta)} | poly    | ~2^{n/beta}   |\n");

console.log("When to use what:\n");
console.log("  Attacking crypto: BKZ with increasing beta");
console.log("  Building crypto: LLL (fast, good enough)");
console.log("  Research/exact: Enumeration or Sieving\n");

// ============================================================================
// Section 6: Security Implications
// ============================================================================

console.log("--- Section 6: Security Implications ---\n");

/**
 * SECURITY OF LATTICE CRYPTOGRAPHY
 *
 * The security of lattice schemes depends on the hardness of
 * gamma-SVP (and related problems) for specific gamma.
 *
 * TYPICAL REDUCTION:
 *
 * Breaking Kyber/Dilithium =>
 *   Solving LWE with specific parameters =>
 *   Solving gamma-SVP with gamma ~ sqrt(n)
 *
 * CONCRETE SECURITY:
 *
 * Security level is measured by the "bit security" - the log of
 * the number of operations needed to break the scheme.
 *
 * Current estimates for SVP in dimension n:
 * - Classical: ~0.29n bits (sieving lower bound)
 * - Quantum: ~0.26n bits (estimated)
 *
 * NIST parameter choices:
 *
 * Kyber-512:  n~768,  ~128-bit security (classical)
 * Kyber-768:  n~1024, ~192-bit security
 * Kyber-1024: n~1280, ~256-bit security
 *
 * ATTACKS IN PRACTICE:
 *
 * Current records (as of 2024):
 * - Exact SVP solved up to n ~ 155
 * - approx-SVP attacks limited by memory
 *
 * Cryptographic parameters are chosen with n ~ 500-1000,
 * providing large security margins.
 *
 * OPEN QUESTIONS:
 *
 * - Exact complexity of SVP (are current algorithms optimal?)
 * - Quantum algorithms for lattice problems
 * - Structured lattice problems (Ring-LWE, Module-LWE)
 */

console.log("Security depends on SVP hardness:\n");
console.log("  Breaking Kyber => Solving LWE => Solving gamma-SVP\n");

console.log("Bit security estimates:");
console.log("  SVP in dimension n:");
console.log("    Classical: ~0.29n bits (sieving)");
console.log("    Quantum:   ~0.26n bits (estimated)\n");

console.log("NIST Kyber parameters:");
console.log("  Kyber-512:  ~128-bit security");
console.log("  Kyber-768:  ~192-bit security");
console.log("  Kyber-1024: ~256-bit security\n");

console.log("Current attack records:");
console.log("  Exact SVP solved for n ~ 155");
console.log("  Crypto uses n ~ 500-1000 (huge margin)\n");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");

console.log("SVP - Shortest Vector Problem:");
console.log("  Find shortest nonzero vector in lattice");
console.log("  NP-hard, best algorithms 2^O(n)\n");

console.log("CVP - Closest Vector Problem:");
console.log("  Find lattice point closest to target");
console.log("  At least as hard as SVP\n");

console.log("Algorithms:");
console.log("  LLL: Fast but only 2^{n/2} approximation");
console.log("  BKZ: Better approximation, slower");
console.log("  Sieving: Optimal time but exponential space\n");

console.log("Security:");
console.log("  No quantum speedup known!");
console.log("  Lattice crypto: parameters chosen for 128+ bit security");
console.log("  Large margin between attacks and crypto parameters\n");

console.log("Next: The LLL Algorithm in detail");
