/**
 * @module sage/crypto/lattice
 * @description Hard Lattice Generator
 *
 * This module contains lattice related functions relevant in cryptography.
 *
 * @see Reference: sage/crypto/lattice.py
 */

import { euler_phi } from '../arith/misc.js';
import { TypeError as SageTypeError, ValueError } from '../errors.js';
import { type RandState, current_randstate, set_random_seed } from '../misc/randstate.js';
import { type FreeModuleIntegerLattice, IntegerLattice } from '../modules/free_module_integer.js';
import { cyclotomic_polynomial } from '../rings/finite_rings/roots_of_unity.js';
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
   * For the type 'ideal', this determines the quotient polynomial, given as
   * its coefficient list in ascending order (constant term first).  It must
   * be monic of degree `n`; e.g. `x^4 - 1` is `[-1n, 0n, 0n, 0n, 1n]`.
   * Ignored for all other types.
   *
   * @see Deviation: SageMath takes a polynomial (or symbolic expression)
   * object here; we take its coefficient list.
   */
  quotient?: IntegerLike[];

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
export function gen_lattice(
  options?: GenLatticeOptions & { ntl?: false; lattice?: false }
): bigint[][];
export function gen_lattice(options: GenLatticeOptions & { ntl: true }): string;
export function gen_lattice(
  options: GenLatticeOptions & { lattice: true }
): FreeModuleIntegerLattice;
export function gen_lattice(
  options?: GenLatticeOptions
): bigint[][] | string | FreeModuleIntegerLattice;
export function gen_lattice(
  options: GenLatticeOptions = {}
): bigint[][] | string | FreeModuleIntegerLattice {
  const {
    type = 'modular',
    n: nOpt = 4n,
    m: mOpt = 8n,
    q: qOpt = 11n,
    seed,
    quotient,
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
  } else if (type === 'cyclotomic' || type === 'ideal') {
    // Both types stack the multiplication matrices of m//n random elements of
    // R = Z_q[x]/(f), where f is the quotient polynomial (Sage:
    // `A = A.stack(R.random_element().matrix())`).
    let f: bigint[];

    if (type === 'cyclotomic') {
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
      // f = Phi_k(x), of degree euler_phi(k) = n
      f = (cyclotomic_polynomial(foundK) as number[]).map((c) => mod(BigInt(c), q));
    } else {
      if (quotient === undefined) {
        throw new ValueError('ideal bases require a quotient polynomial');
      }
      f = quotient.map((c) => mod(toBigInt(c), q));
      // Strip leading zeros to determine the degree.
      while (f.length > 0 && f[f.length - 1] === 0n) {
        f.pop();
      }
      if (f.length - 1 !== n) {
        throw new ValueError('ideal basis requires n = quotient.degree()');
      }
      if (f[n] !== 1n) {
        throw new SageTypeError('quotient must be monic');
      }
    }

    for (let i = 0; i < Math.floor(m / n); i++) {
      // Generate a random element of R (a polynomial of degree < n over Z_q)
      const poly: bigint[] = [];
      for (let j = 0; j < n; j++) {
        poly.push(mod(randState.randint(0n, q - 1n), q));
      }

      // Multiplication matrix of `poly` on the power basis 1, x, ..., x^{n-1}:
      // row j holds the coefficients of x^j * poly mod f.
      let current = [...poly];
      for (let row = 0; row < n; row++) {
        A.push([...current]);
        current = mulByXModF(current, f, n, q);
      }
    }

    // Sage slices A[n:m]; drop anything beyond m rows.
    A = A.slice(0, m);
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
    // Dual basis: block_matrix([[1, -A'^T], [0, q]]).
    // A' is (m-n) x n, so -A'^T is n x (m-n): the identity block is I_n and
    // the q block is q*I_{m-n}.  Then reverse the row order.
    B = [];

    // First n rows: I_n, then -A'^T
    for (let i = 0; i < n; i++) {
      const row: bigint[] = [];
      // Identity in first n columns
      for (let j = 0; j < n; j++) {
        row.push(i === j ? 1n : 0n);
      }
      // -A'^T in last (m - n) columns: -A'[j][i] for each j
      for (let j = 0; j < m - n; j++) {
        row.push(-A_prime[j]![i]!);
      }
      B.push(row);
    }

    // Last (m - n) rows: zeros, then q*I_{m-n}
    for (let i = 0; i < m - n; i++) {
      const row: bigint[] = [];
      for (let j = 0; j < n; j++) {
        row.push(0n);
      }
      for (let j = 0; j < m - n; j++) {
        row.push(i === j ? q : 0n);
      }
      B.push(row);
    }

    // Reverse row order
    for (let i = 0; i < Math.floor(m / 2); i++) {
      [B[i], B[m - i - 1]] = [B[m - i - 1]!, B[i]!];
    }
  }

  // Sage: `return B._ntl_()`, whose string form is the bracketed NTL matrix.
  if (ntl) {
    return `[\n${B.map((row) => `[${row.join(' ')}]`).join('\n')}\n]`;
  }

  // Sage: `return IntegerLattice(B)`
  if (lattice) {
    return IntegerLattice(B);
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

/**
 * Multiply a polynomial (coefficient list, ascending, length n) by x modulo a
 * monic polynomial `f` of degree n over Z_q.
 *
 * @param c - Coefficients of the polynomial, ascending, length n
 * @param f - Coefficients of the monic modulus, ascending, length n+1
 * @param n - Degree of the modulus
 * @param q - Modulus of the coefficient ring
 */
function mulByXModF(c: bigint[], f: bigint[], n: number, q: bigint): bigint[] {
  // Shift by one, then reduce the x^n term using x^n = -(f_0 + ... + f_{n-1} x^{n-1})
  const lead = c[n - 1] ?? 0n;
  const out = new Array<bigint>(n);
  out[0] = mod(-lead * f[0]!, q);
  for (let i = 1; i < n; i++) {
    out[i] = mod(c[i - 1]! - lead * f[i]!, q);
  }
  return out;
}
