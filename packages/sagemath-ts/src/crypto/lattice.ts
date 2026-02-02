/**
 * @module sage/crypto/lattice
 * @description Hard Lattice Generator
 *
 * This module contains lattice related functions relevant in cryptography.
 *
 * @see Reference: sage/crypto/lattice.py
 */

import { euler_phi } from '../arith/misc.js';
import { ValueError } from '../errors.js';
import { type RandState, current_randstate, set_random_seed } from '../misc/randstate.js';
import { type IntegerLike, toBigInt, toSafeNumber } from '../types/coercion.js';

/**
 * Type of lattice to generate.
 */
export type LatticeType = 'modular' | 'random' | 'ideal' | 'cyclotomic';

/**
 * Options for gen_lattice function.
 */
export interface GenLatticeOptions {
  /**
   * Type of lattice to generate:
   * - 'modular' (default): a class of lattices for which asymptotic worst-case
   *   to average-case connections hold
   * - 'random': special case of modular (n=1); a dense class of lattice used
   *   for testing basis reduction algorithms
   * - 'ideal': special case of modular; allows for a more compact representation
   * - 'cyclotomic': special case of ideal; allows for efficient processing
   */
  type?: LatticeType;

  /**
   * Determinant size, primal: det(L) = q^n, dual: det(L) = q^{m-n}.
   * For ideal lattices this is also the degree of the quotient polynomial.
   * @default 4
   */
  n?: IntegerLike;

  /**
   * Lattice dimension, L is a subset of Z^m.
   * @default 8
   */
  m?: IntegerLike;

  /**
   * Coefficient size, q*Z^m is a subset of L.
   * @default 11
   */
  q?: IntegerLike;

  /**
   * Randomness seed for reproducible lattice generation.
   */
  seed?: number;

  /**
   * For the type 'ideal', this determines the quotient polynomial.
   * Ignored for all other types.
   */
  quotient?: unknown;

  /**
   * Set this flag if you want a basis for q-dual(L), for example
   * for Regev's LWE bases.
   * @default false
   */
  dual?: boolean;

  /**
   * Set this flag if you want the lattice basis in NTL readable format.
   * @default false
   */
  ntl?: boolean;

  /**
   * Set this flag if you want a FreeModule_submodule_with_basis_integer
   * object instead of an integer matrix representing the basis.
   * @default false
   */
  lattice?: boolean;
}

/**
 * Generate different types of integral lattice bases of row vectors
 * relevant in cryptography.
 *
 * This function generates different types of integral lattice bases
 * of row vectors relevant in cryptography. Randomness can be set either
 * with `seed`, or by using a random seed.
 *
 * @param options - Configuration options for lattice generation
 * @returns A unique size-reduced triangular basis of row vectors for the
 *          lattice in question. Primal: lower_left, dual: lower_right.
 *          Returns a 2D array of bigints representing the matrix.
 *
 * @example
 * // Modular basis
 * const B = gen_lattice({ m: 10, seed: 42 });
 *
 * @example
 * // Random basis
 * const B = gen_lattice({ type: 'random', n: 1, m: 10, q: 11n**4n, seed: 42 });
 *
 * @example
 * // Cyclotomic bases with n=2^k are SWIFFT bases
 * const B = gen_lattice({ type: 'cyclotomic', seed: 42 });
 *
 * @example
 * // Dual modular bases are related to Regev's public-key encryption
 * const B = gen_lattice({ type: 'modular', m: 10, seed: 42, dual: true });
 *
 * @see Reference: sage/crypto/lattice.py:gen_lattice
 */
export function gen_lattice(options: GenLatticeOptions = {}): bigint[][] {
  const {
    type = 'modular',
    n: nOpt = 4n,
    m: mOpt = 8n,
    q: qOpt = 11n,
    seed,
    dual = false,
    ntl = false,
    lattice = false,
  } = options;

  const n = toSafeNumber(toBigInt(nOpt));
  const m = toSafeNumber(toBigInt(mOpt));
  const q = toBigInt(qOpt);

  // Set random seed if provided
  let randState: RandState;
  if (seed !== undefined) {
    set_random_seed(seed);
    randState = current_randstate();
  } else {
    randState = current_randstate();
  }

  // Validate parameters
  if (type === 'random' && n !== 1) {
    throw new ValueError('random bases require n = 1');
  }

  if (ntl && lattice) {
    throw new ValueError('Cannot specify ntl=True and lattice=True at the same time');
  }

  // Generate the random part of the matrix A
  // A starts as an n x n identity matrix over Z_q
  // Then we stack random rows below it

  let A: bigint[][] = [];

  // Initialize A with identity matrix (n x n)
  for (let i = 0; i < n; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < n; j++) {
      row.push(i === j ? 1n : 0n);
    }
    A.push(row);
  }

  if (type === 'random' || type === 'modular') {
    // Add (m - n) random rows
    for (let i = 0; i < m - n; i++) {
      const row: bigint[] = [];
      for (let j = 0; j < n; j++) {
        row.push(mod(randState.randint(0n, q - 1n), q));
      }
      A.push(row);
    }
  } else if (type === 'cyclotomic') {
    // For cyclotomic, find k such that euler_phi(k) = n
    // We assume n+1 <= min(euler_phi^{-1}(n)) <= 2*n
    let foundK = -1;
    for (let k = 2 * n; k > n; k--) {
      if (euler_phi(BigInt(k)) === BigInt(n)) {
        foundK = k;
        break;
      }
    }
    if (foundK === -1) {
      throw new ValueError(
        "cyclotomic bases require that n is an image of Euler's totient function"
      );
    }

    // Generate m/n random elements in the quotient ring Z_q[x]/(x^n + 1)
    // Each element is represented as a polynomial, which gives a circulant-like matrix
    for (let i = 0; i < Math.floor(m / n); i++) {
      // Generate a random polynomial of degree < n
      const poly: bigint[] = [];
      for (let j = 0; j < n; j++) {
        poly.push(mod(randState.randint(0n, q - 1n), q));
      }

      // Add the multiplication matrix of this polynomial in the quotient ring
      // For cyclotomic: x^n = -1, so multiplication by x rotates and negates
      for (let row = 0; row < n; row++) {
        const matRow: bigint[] = [];
        for (let col = 0; col < n; col++) {
          // Coefficient of x^col in poly * x^row mod (x^n + 1)
          // This is poly[(col - row) mod n] with sign adjustment for wrap-around
          const idx = (col - row + n) % n;
          let coeff = poly[idx]!;
          if (col < row) {
            // Wrapped around, negate due to x^n = -1
            coeff = mod(-coeff, q);
          }
          matRow.push(coeff);
        }
        if (A.length < m) {
          A.push(matRow);
        }
      }
    }

    // Trim to exactly m rows if needed
    A = A.slice(0, m);
  } else if (type === 'ideal') {
    throw new ValueError('ideal bases require a quotient polynomial - use cyclotomic for x^n+1');
  }

  // Convert representatives from [0, q-1] to [-(q-1)/2, (q-1)/2]
  function minrep(a: bigint): bigint {
    const halfQ = q / 2n;
    if (a > halfQ) {
      return a - q;
    }
    return a;
  }

  // Apply minrep to the bottom (m - n) rows
  const A_prime: bigint[][] = [];
  for (let i = n; i < A.length; i++) {
    A_prime.push(A[i]!.map(minrep));
  }

  // Build the final basis matrix B
  let B: bigint[][];

  if (!dual) {
    // Primal basis: [[q*I_n, 0], [A', I_{m-n}]]
    // This is an m x m lower-triangular matrix
    B = [];

    // First n rows: q on diagonal, zeros elsewhere
    for (let i = 0; i < n; i++) {
      const row: bigint[] = [];
      for (let j = 0; j < m; j++) {
        row.push(i === j ? q : 0n);
      }
      B.push(row);
    }

    // Last (m - n) rows: A' in first n columns, identity in last (m - n) columns
    for (let i = 0; i < m - n; i++) {
      const row: bigint[] = [];
      // First n columns: A'[i]
      for (let j = 0; j < n; j++) {
        row.push(A_prime[i]?.[j] ?? 0n);
      }
      // Last (m - n) columns: identity
      for (let j = 0; j < m - n; j++) {
        row.push(i === j ? 1n : 0n);
      }
      B.push(row);
    }
  } else {
    // Dual basis: [[I_{m-n}, -A'^T], [0, q*I_n]]
    // Then reverse the row order
    B = [];

    // First (m - n) rows: identity, then -A'^T
    for (let i = 0; i < m - n; i++) {
      const row: bigint[] = [];
      // Identity in first (m - n) columns
      for (let j = 0; j < m - n; j++) {
        row.push(i === j ? 1n : 0n);
      }
      // -A'^T in last n columns: -A'[j][i] for each j
      for (let j = 0; j < n; j++) {
        row.push(-(A_prime[j]?.[i] ?? 0n));
      }
      B.push(row);
    }

    // Last n rows: zeros, then q*I_n
    for (let i = 0; i < n; i++) {
      const row: bigint[] = [];
      for (let j = 0; j < m - n; j++) {
        row.push(0n);
      }
      for (let j = 0; j < n; j++) {
        row.push(i === j ? q : 0n);
      }
      B.push(row);
    }

    // Reverse row order
    for (let i = 0; i < Math.floor(m / 2); i++) {
      [B[i], B[m - i - 1]] = [B[m - i - 1]!, B[i]!];
    }
  }

  // Handle ntl format (just return as string representation)
  if (ntl) {
    // For NTL format, we'd return a string, but since TypeScript
    // prefers typed returns, we return the matrix
    // A proper implementation would format this as NTL expects
    return B;
  }

  // Handle lattice format (return as module with basis)
  if (lattice) {
    // A full implementation would return an IntegerLattice object
    // For now, return the basis matrix
    return B;
  }

  return B;
}

/**
 * Modular reduction to [0, m)
 */
function mod(a: bigint, m: bigint): bigint {
  const r = a % m;
  return r < 0n ? r + m : r;
}
