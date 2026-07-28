/**
 * @module sage/matrix/matrix_operations
 * @description Matrix operations (determinant, rank, kernel, etc.)
 *
 * Port of: sage/matrix/matrix2.pyx
 */

import {
  ArithmeticError,
  NotImplementedError,
  TypeError,
  ValueError,
  ZeroDivisionError,
} from '../errors.js';
import { FreeModule, type FreeModuleGeneric, VectorSpace } from '../modules/free_module.js';
import type { FreeModuleElement } from '../modules/free_module_element.js';
import type { CoefficientRing, RingElement } from '../rings/polynomial/polynomial_element.js';
import { Polynomial } from '../rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';
import { Rational } from '../rings/rational.js';
import {
  echelon_form,
  hessenberg_form,
  jordan_form,
  pivots,
  rref,
} from './matrix_decompositions.js';
import { Matrix, identity_matrix, zero_matrix } from './matrix_generic.js';
import { prod_of_row_sums } from './matrix_special.js';

/**
 * Interface for field elements that support division/inverse.
 */
interface FieldElement extends RingElement {
  inverse?(): this;
  inv?(): this;
  div?(other: this): this;
}

/**
 * Get the multiplicative inverse of a field element.
 */
function getInverse<R extends FieldElement>(elem: R): R {
  if (typeof elem.inverse === 'function') {
    return elem.inverse() as R;
  }
  if (typeof elem.inv === 'function') {
    return elem.inv() as R;
  }
  throw new ArithmeticError('element does not support inverse operation');
}

/**
 * Check if a ring is a field.
 */
function isFieldRing(ring: unknown): boolean {
  if (
    ring !== null &&
    typeof ring === 'object' &&
    'is_field' in ring &&
    typeof (ring as { is_field: () => boolean }).is_field === 'function'
  ) {
    return (ring as { is_field: () => boolean }).is_field();
  }
  return false;
}

// ============================================================================
// Determinant and related operations
// ============================================================================

/**
 * Return the determinant of the matrix.
 *
 * ALGORITHM:
 *
 * For small matrices (n <= 3), uses the naive formula (`matrix2.pyx:2368-2380`).
 * For larger matrices the determinant is read off from the characteristic
 * polynomial (`matrix2.pyx:2409-2443`): `det(A) = (-1)^n * charpoly(A)[0]`.
 * The default charpoly algorithm is the division-free one of [Sei2002], so
 * this works over any commutative ring (in particular over `Z/nZ` with `n`
 * composite, where Gaussian elimination would divide by a zero divisor).
 *
 * @param matrix - The square matrix
 * @param algorithm - One of 'df' (division-free) or 'hessenberg'
 * @returns The determinant as an element of the base ring
 * @see Reference: sage/matrix/matrix2.pyx:determinant
 */
export function determinant<R extends RingElement>(
  matrix: Matrix<R>,
  algorithm?: 'df' | 'hessenberg'
): R {
  if (!matrix.is_square()) {
    throw new ArithmeticError('determinant is only defined for square matrices');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return ring.one();
  }

  if (n === 1) {
    return matrix.get(0, 0);
  }

  if (n === 2) {
    // det = a*d - b*c
    const a = matrix.get(0, 0);
    const b = matrix.get(0, 1);
    const c = matrix.get(1, 0);
    const d = matrix.get(1, 1);
    return a.mul(d).sub(b.mul(c)) as R;
  }

  if (n === 3) {
    // Sarrus rule
    const a = matrix.get(0, 0);
    const b = matrix.get(0, 1);
    const c = matrix.get(0, 2);
    const d = matrix.get(1, 0);
    const e = matrix.get(1, 1);
    const f = matrix.get(1, 2);
    const g = matrix.get(2, 0);
    const h = matrix.get(2, 1);
    const i = matrix.get(2, 2);

    // a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h
    const pos1 = a.mul(e).mul(i) as R;
    const pos2 = b.mul(f).mul(g) as R;
    const pos3 = c.mul(d).mul(h) as R;
    const neg1 = c.mul(e).mul(g) as R;
    const neg2 = b.mul(d).mul(i) as R;
    const neg3 = a.mul(f).mul(h) as R;

    return pos1.add(pos2).add(pos3).sub(neg1).sub(neg2).sub(neg3) as R;
  }

  // Generic algorithm (matrix2.pyx:2409-2443): read the determinant off the
  // characteristic polynomial.  charp[0] = det(-A) = (-1)^n det(A).
  const charp = charpoly(matrix, 'x', algorithm ?? 'df');
  const c = charp.getCoeff(0);
  return (n % 2 === 1 ? c.neg() : c) as R;
}

/**
 * Alias for determinant.
 * @see Reference: sage/matrix/matrix2.pyx:det
 */
export const det = determinant;

/**
 * Return the quantum determinant of the matrix.
 *
 * The quantum determinant of a matrix M = (m_{ij})_{i,j=1}^n is defined by
 * det_q(M) = sum_{sigma in S_n} (-q)^{l(sigma)} prod_{i} M_{sigma(i),i}
 * where l(sigma) is the number of inversions in the permutation sigma.
 *
 * When q = 1, this equals the ordinary determinant.
 * When q = -1, this equals the permanent.
 *
 * @param matrix - The square matrix
 * @param q - The parameter q
 * @returns The quantum determinant
 * @see Reference: sage/matrix/matrix2.pyx:quantum_determinant
 */
export function quantum_determinant<R extends RingElement>(matrix: Matrix<R>, q: R): R {
  if (!matrix.is_square()) {
    throw new ArithmeticError('quantum_determinant is only defined for square matrices');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return ring.one();
  }

  // Generate all permutations and compute the quantum determinant
  let result = ring.zero();
  const perms = _generatePermutations(n);

  for (const perm of perms) {
    // Compute the number of inversions
    let inversions = 0;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        if (perm[i]! > perm[j]!) {
          inversions++;
        }
      }
    }

    // Compute (-q)^{inversions}
    let sign = ring.one();
    for (let i = 0; i < inversions; i++) {
      sign = sign.mul(q.neg()) as R;
    }

    // Compute product M_{sigma(i), i}
    let prod = ring.one();
    for (let i = 0; i < n; i++) {
      prod = prod.mul(matrix.get(perm[i]!, i)) as R;
    }

    result = result.add(sign.mul(prod)) as R;
  }

  return result;
}

/**
 * Generate all permutations of [0, 1, ..., n-1].
 */
function _generatePermutations(n: number): number[][] {
  if (n === 0) return [[]];
  if (n === 1) return [[0]];

  const result: number[][] = [];

  function backtrack(current: number[], used: Set<number>) {
    if (current.length === n) {
      result.push([...current]);
      return;
    }
    for (let i = 0; i < n; i++) {
      if (!used.has(i)) {
        current.push(i);
        used.add(i);
        backtrack(current, used);
        current.pop();
        used.delete(i);
      }
    }
  }

  backtrack([], new Set());
  return result;
}

/**
 * Return the Pfaffian of a skew-symmetric matrix.
 *
 * The Pfaffian is defined for skew-symmetric matrices of even order.
 * For a 2n x 2n skew-symmetric matrix A, pf(A)^2 = det(A).
 *
 * @param matrix - A skew-symmetric square matrix of even dimension
 * @param algorithm - Algorithm to use ('definition' or 'bfl')
 * @param check - Whether to check that the matrix is skew-symmetric (default: true)
 * @returns The Pfaffian
 * @see Reference: sage/matrix/matrix2.pyx:pfaffian
 */
export function pfaffian<R extends RingElement>(
  matrix: Matrix<R>,
  algorithm: string = 'definition',
  check: boolean = true
): R {
  if (!matrix.is_square()) {
    throw new ArithmeticError('pfaffian is only defined for square matrices');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n % 2 !== 0) {
    // Pfaffian of odd-dimensional skew-symmetric matrix is 0
    return ring.zero();
  }

  if (n === 0) {
    return ring.one();
  }

  // Check skew-symmetry: A[i,j] = -A[j,i]
  if (check) {
    for (let i = 0; i < n; i++) {
      for (let j = i; j < n; j++) {
        if (i === j) {
          if (!matrix.get(i, j).isZero()) {
            throw new ValueError('matrix is not skew-symmetric: diagonal must be zero');
          }
        } else {
          const aij = matrix.get(i, j);
          const aji = matrix.get(j, i);
          if (!aij.add(aji).isZero()) {
            throw new ValueError('matrix is not skew-symmetric: A[i,j] + A[j,i] != 0');
          }
        }
      }
    }
  }

  // Use recursive definition for Pfaffian
  // pf(A) = sum_{j=2}^{n} (-1)^j * A[1,j] * pf(A_{1j})
  // where A_{1j} is the matrix with rows/cols 1 and j removed

  if (n === 2) {
    // pf([[0, a], [-a, 0]]) = a
    return matrix.get(0, 1);
  }

  // Recursive computation
  let result = ring.zero();
  for (let j = 1; j < n; j++) {
    const entry = matrix.get(0, j);
    if (!entry.isZero()) {
      // Create submatrix by removing rows and columns 0 and j
      const subRows: number[] = [];
      for (let i = 1; i < n; i++) {
        if (i !== j) subRows.push(i);
      }
      const subCols = subRows; // Same for skew-symmetric

      const subEntries: R[][] = [];
      for (const i of subRows) {
        const row: R[] = [];
        for (const k of subCols) {
          row.push(matrix.get(i, k));
        }
        subEntries.push(row);
      }

      const subMatrix = new Matrix(ring, subRows.length, subCols.length, subEntries);
      const subPf = pfaffian(subMatrix, algorithm, false);

      // (-1)^j = (-1)^(j-1+1) for 1-indexed, here j is 0-indexed so (-1)^j
      const sign = j % 2 === 0 ? ring.one().neg() : ring.one();
      result = result.add(sign.mul(entry).mul(subPf) as R) as R;
    }
  }

  return result;
}

/**
 * Return the permanent of the matrix.
 *
 * The permanent is the same as the determinant but without the sign changes.
 * For an m x n matrix with m <= n, the permanent is:
 *   per(A) = sum over all injections pi from {1,...,m} to {1,...,n}
 *            of a_{1,pi(1)} * a_{2,pi(2)} * ... * a_{m,pi(m)}
 *
 * @param matrix - The matrix
 * @param algorithm - 'Ryser' (default) or 'definition'
 * @returns The permanent
 * @see Reference: sage/matrix/matrix2.pyx:permanent
 */
export function permanent<R extends RingElement>(
  matrix: Matrix<R>,
  algorithm: 'Ryser' | 'definition' = 'Ryser'
): R {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  if (m === 0) {
    return ring.one();
  }

  if (m > n) {
    throw new ValueError(`permanent requires nrows <= ncols, but got ${m} x ${n}`);
  }

  if (algorithm === 'definition') {
    return _permanent_definition(matrix);
  }

  // Ryser's algorithm
  return _permanent_ryser(matrix);
}

/**
 * Compute permanent using Ryser's algorithm.
 * This is O(2^n * n) and is faster than the naive O(n! * n) algorithm.
 */
function _permanent_ryser<R extends RingElement>(matrix: Matrix<R>): R {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  let perm = ring.zero();

  // Iterate over all non-empty subsets of columns
  for (let r = 1; r <= m; r++) {
    const combinations = _choose_indices(n, r);
    let s = ring.zero();

    for (const cols of combinations) {
      // Compute product of row sums for selected columns
      s = s.add(prod_of_row_sums(matrix, cols)) as R;
    }

    // Apply sign and binomial coefficient
    // perm += (-1)^(m-r) * C(n-r, m-r) * s
    const sign = (m - r) % 2 === 0 ? 1 : -1;
    const binom = _binomial(n - r, m - r);

    if (sign > 0) {
      for (let i = 0; i < binom; i++) {
        perm = perm.add(s) as R;
      }
    } else {
      for (let i = 0; i < binom; i++) {
        perm = perm.sub(s) as R;
      }
    }
  }

  return perm;
}

/**
 * Compute permanent using the definition (sum over permutations).
 * This is O(n! * n) and is much slower than Ryser's algorithm.
 */
function _permanent_definition<R extends RingElement>(matrix: Matrix<R>): R {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  if (m === 0) {
    return ring.one();
  }

  let perm = ring.zero();

  // Generate all injections from {0,...,m-1} to {0,...,n-1}
  const injections = _generate_injections(m, n);

  for (const pi of injections) {
    let prod = ring.one();
    for (let i = 0; i < m; i++) {
      prod = prod.mul(matrix.get(i, pi[i]!)) as R;
    }
    perm = perm.add(prod) as R;
  }

  return perm;
}

/**
 * Generate all injections from {0,...,m-1} to {0,...,n-1}.
 */
function _generate_injections(m: number, n: number): number[][] {
  if (m === 0) {
    return [[]];
  }

  const result: number[][] = [];

  function backtrack(current: number[], used: Set<number>) {
    if (current.length === m) {
      result.push([...current]);
      return;
    }

    for (let i = 0; i < n; i++) {
      if (!used.has(i)) {
        current.push(i);
        used.add(i);
        backtrack(current, used);
        current.pop();
        used.delete(i);
      }
    }
  }

  backtrack([], new Set());
  return result;
}

/**
 * Compute binomial coefficient C(n, k).
 */
function _binomial(n: number, k: number): number {
  if (k < 0 || k > n) {
    return 0;
  }
  if (k === 0 || k === n) {
    return 1;
  }

  // Use symmetry: C(n, k) = C(n, n-k)
  if (k > n - k) {
    k = n - k;
  }

  let result = 1;
  for (let i = 0; i < k; i++) {
    result = (result * (n - i)) / (i + 1);
  }
  return Math.round(result);
}

/**
 * Generate all k-combinations of indices from {0, ..., n-1}.
 */
function _choose_indices(n: number, k: number): number[][] {
  const result: number[][] = [];

  function backtrack(start: number, current: number[]) {
    if (current.length === k) {
      result.push([...current]);
      return;
    }

    for (let i = start; i < n; i++) {
      current.push(i);
      backtrack(i + 1, current);
      current.pop();
    }
  }

  backtrack(0, []);
  return result;
}

/**
 * Return the k-th permanental minor of the matrix.
 *
 * The permanental k-minor is the sum of the permanents of all k x k submatrices.
 *
 * @param matrix - The matrix
 * @param k - The size of submatrices
 * @param algorithm - 'Ryser' (default) or 'definition'
 * @returns The sum of permanents of k x k submatrices
 * @see Reference: sage/matrix/matrix2.pyx:permanental_minor
 */
export function permanental_minor<R extends RingElement>(
  matrix: Matrix<R>,
  k: number,
  algorithm: 'Ryser' | 'definition' = 'Ryser'
): R {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  if (k === 0) {
    return ring.one();
  }

  if (k > m || k > n) {
    return ring.zero();
  }

  let pm = ring.zero();

  // Iterate over all k x k submatrices
  for (const rowIndices of _choose_indices(m, k)) {
    for (const colIndices of _choose_indices(n, k)) {
      const submatrix = _submatrix(matrix, rowIndices, colIndices);
      pm = pm.add(permanent(submatrix, algorithm)) as R;
    }
  }

  return pm;
}

/**
 * Extract a submatrix from the given rows and columns.
 */
function _submatrix<R extends RingElement>(
  matrix: Matrix<R>,
  rowIndices: number[],
  colIndices: number[]
): Matrix<R> {
  const m = rowIndices.length;
  const n = colIndices.length;
  const ring = matrix.base_ring;

  const entries: R[][] = [];
  for (let i = 0; i < m; i++) {
    entries.push([]);
    for (let j = 0; j < n; j++) {
      entries[i]!.push(matrix.get(rowIndices[i]!, colIndices[j]!));
    }
  }

  return new Matrix(ring, m, n, entries);
}

/**
 * Return all k x k minors of the matrix.
 *
 * A minor is the determinant of a k x k submatrix.
 *
 * @param matrix - The matrix
 * @param k - The size of submatrices
 * @returns List of all k x k minors
 * @see Reference: sage/matrix/matrix2.pyx:minors
 */
export function minors<R extends RingElement>(matrix: Matrix<R>, k: number): R[] {
  const m = matrix.nrows;
  const n = matrix.ncols;

  if (k < 0 || k > m || k > n) {
    return [];
  }

  const result: R[] = [];

  // Iterate over all k x k submatrices
  for (const rowIndices of _choose_indices(m, k)) {
    for (const colIndices of _choose_indices(n, k)) {
      const submatrix = _submatrix(matrix, rowIndices, colIndices);
      result.push(determinant(submatrix));
    }
  }

  return result;
}

// ============================================================================
// Inverse and pseudo-inverse
// ============================================================================

/**
 * Return the inverse of the matrix.
 *
 * Uses Gaussian elimination (row reduction) to compute the inverse.
 * Augments the matrix with the identity and reduces to get [I | A^{-1}].
 *
 * @param matrix - The invertible square matrix
 * @returns The inverse matrix
 * @throws {ArithmeticError} If the matrix is not square or not invertible
 * @see Reference: sage/matrix/matrix2.pyx:inverse
 */
export function inverse<R extends FieldElement>(matrix: Matrix<R>): Matrix<R> {
  if (!matrix.is_square()) {
    throw new ArithmeticError('inverse is only defined for square matrices');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return zero_matrix(ring, 0);
  }

  // Create augmented matrix [A | I]
  const aug: R[][] = [];
  for (let i = 0; i < n; i++) {
    aug.push([]);
    for (let j = 0; j < n; j++) {
      aug[i]!.push(matrix.get(i, j));
    }
    for (let j = 0; j < n; j++) {
      aug[i]!.push(i === j ? ring.one() : ring.zero());
    }
  }

  // Forward elimination with partial pivoting
  for (let col = 0; col < n; col++) {
    // Find pivot
    let pivotRow = -1;
    for (let i = col; i < n; i++) {
      if (!aug[i]![col]!.isZero()) {
        pivotRow = i;
        break;
      }
    }

    if (pivotRow === -1) {
      throw new ArithmeticError('matrix is singular and cannot be inverted');
    }

    // Swap rows
    if (pivotRow !== col) {
      [aug[col], aug[pivotRow]] = [aug[pivotRow]!, aug[col]!];
    }

    // Scale pivot row to make pivot = 1
    const pivot = aug[col]![col]!;
    const pivotInv = getInverse(pivot);

    for (let j = col; j < 2 * n; j++) {
      aug[col]![j] = aug[col]![j]!.mul(pivotInv) as R;
    }

    // Eliminate column entries
    for (let i = 0; i < n; i++) {
      if (i !== col && !aug[i]![col]!.isZero()) {
        const factor = aug[i]![col]!;
        for (let j = col; j < 2 * n; j++) {
          aug[i]![j] = aug[i]![j]!.sub(factor.mul(aug[col]![j]!) as R) as R;
        }
      }
    }
  }

  // Extract inverse from right half
  const result = zero_matrix(ring, n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      result.set(i, j, aug[i]![n + j]!);
    }
  }

  return result;
}

/**
 * Return the Moore-Penrose pseudoinverse of the matrix.
 *
 * The pseudoinverse A^+ satisfies the Moore-Penrose conditions:
 * - A * A^+ * A = A
 * - A^+ * A * A^+ = A^+
 * - (A * A^+)^* = A * A^+
 * - (A^+ * A)^* = A^+ * A
 *
 * For full-rank matrices, the pseudoinverse equals the left or right inverse.
 * For general matrices, this requires SVD which needs sqrt.
 *
 * @param matrix - The matrix
 * @param algorithm - Algorithm to use
 * @returns The pseudoinverse matrix
 * @see Reference: sage/matrix/matrix2.pyx:pseudoinverse
 */
export function pseudoinverse<R extends RingElement>(
  matrix: Matrix<R>,
  algorithm?: string
): Matrix<R> {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  if (m === 0 || n === 0) {
    return new Matrix(ring, n, m);
  }

  // For full column rank: A^+ = (A^T * A)^{-1} * A^T
  // For full row rank: A^+ = A^T * (A * A^T)^{-1}

  const r = rank(matrix);

  if (r === n && n <= m) {
    // Full column rank: A^+ = (A^T A)^{-1} A^T
    const AT = matrix.transpose();
    const ATA = AT.mul(matrix);
    try {
      const ATAinv = inverse(ATA);
      return ATAinv.mul(AT);
    } catch {
      // Not invertible, fall through
    }
  }

  if (r === m && m <= n) {
    // Full row rank: A^+ = A^T (A A^T)^{-1}
    const AT = matrix.transpose();
    const AAT = matrix.mul(AT);
    try {
      const AATinv = inverse(AAT);
      return AT.mul(AATinv);
    } catch {
      // Not invertible, fall through
    }
  }

  // General case requires SVD which needs sqrt
  throw new NotImplementedError(
    'pseudoinverse for rank-deficient matrices requires SVD, which needs sqrt'
  );
}

/**
 * Return the adjugate (classical adjoint) of the matrix.
 *
 * The adjugate is the transpose of the matrix of cofactors:
 * `adj(M)[i,j] = (-1)^(i+j) det(M_{j,i})`, and `adj(M) M = M adj(M) = det(M) I`.
 *
 * ALGORITHM (`matrix2.pyx:_adjugate`, Algorithm 3.1 of [Sei2002]): the
 * division-free formula `A = charpoly().shift(-1)(self)`, negated when the
 * size is even.  Unlike cofactor expansion this never divides, so it works
 * over rings such as `Z/8Z`.
 *
 * @param matrix - The square matrix
 * @returns The adjugate matrix
 * @throws {ValueError} If the matrix is not square
 * @see Reference: sage/matrix/matrix2.pyx:adjugate
 */
export function adjugate<R extends RingElement>(matrix: Matrix<R>): Matrix<R> {
  if (!matrix.is_square()) {
    throw new ValueError('must be a square matrix');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return zero_matrix(ring, 0);
  }

  const A = _evaluate_at_matrix(charpoly(matrix).shift(-1), matrix);
  return n % 2 === 1 ? A : A.neg();
}

/**
 * Alias for adjugate.
 * @see Reference: sage/matrix/matrix2.pyx:adjoint_classical
 */
export const adjoint_classical = adjugate;

// ============================================================================
// Rank and nullity
// ============================================================================

/**
 * Return the rank of the matrix.
 *
 * The rank is the number of pivot columns in the echelon form,
 * which equals the dimension of the row space (or column space).
 *
 * @param matrix - The matrix
 * @returns The rank
 * @see Reference: sage/matrix/matrix2.pyx:rank (via matrix_integer_dense.pyx)
 */
export function rank<R extends FieldElement>(matrix: Matrix<R>): number {
  // Rank = number of pivots in echelon form
  return pivots(matrix).length;
}

/**
 * Return the (left) nullity of the matrix.
 *
 * The left nullity is the dimension of the left kernel.
 * For an m x n matrix, left_nullity = m - rank.
 *
 * @param matrix - The matrix
 * @returns The left nullity
 * @see Reference: sage/matrix/matrix2.pyx:left_nullity
 */
export function left_nullity<R extends FieldElement>(matrix: Matrix<R>): number {
  // Left nullity = nrows - rank (since matrices act on row vectors from the right)
  return matrix.nrows - rank(matrix);
}

/**
 * Alias for left_nullity.
 * @see Reference: sage/matrix/matrix2.pyx:nullity
 */
export const nullity = left_nullity;

/**
 * Return the right nullity of the matrix.
 *
 * The right nullity is the dimension of the right kernel.
 * For an m x n matrix, right_nullity = n - rank.
 *
 * @param matrix - The matrix
 * @returns The right nullity
 * @see Reference: sage/matrix/matrix2.pyx:right_nullity
 */
export function right_nullity<R extends FieldElement>(matrix: Matrix<R>): number {
  // Right nullity = ncols - rank
  return matrix.ncols - rank(matrix);
}

// ============================================================================
// Kernel and Image
// ============================================================================

/**
 * Return a matrix whose rows are a basis for the right kernel of the matrix.
 *
 * The right kernel of A is the set of column vectors x such that Ax = 0.
 * Equivalently, it is the null space of A.
 *
 * ALGORITHM (`matrix2.pyx:_right_kernel_matrix_over_field`): compute the
 * echelon form and its pivots; for every non-pivot column `i` emit the vector
 * with a 1 in position `i` and `-E[r,i]` in the `r`-th pivot position.  That
 * gives the `'pivot'` basis; Sage's *default* basis over a field is
 * `'echelon'` (`matrix2.pyx:4865`), so the result is echelonized before being
 * returned.
 *
 * @param matrix - The matrix
 * @param options - Optional parameters
 * @param options.basis - Basis format: 'default'/'echelon' (default), 'pivot'
 *   or 'computed'
 * @returns A matrix whose rows span the right kernel
 * @see Reference: sage/matrix/matrix2.pyx:right_kernel_matrix
 */
export function right_kernel_matrix<R extends FieldElement>(
  matrix: Matrix<R>,
  options?: { basis?: 'default' | 'computed' | 'echelon' | 'pivot' }
): Matrix<R> {
  const basis = options?.basis ?? 'default';
  if (!['default', 'computed', 'echelon', 'pivot'].includes(basis)) {
    throw new ValueError(`matrix kernel basis format '${basis}' not recognized`);
  }
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  // Get RREF and pivots. These must be pivot COLUMNS: the loop below tests
  // column indices against pivotSet and uses the position within this list as
  // the RREF row index. `pivot_rows` returns row indices (it is
  // `pivots(transpose)`), so using it here silently produced vectors that were
  // not in the kernel.
  const E = rref(matrix);
  const pivotCols = pivots(matrix);
  const pivotSet = new Set(pivotCols);

  // Find non-pivot columns
  const nonPivotCols: number[] = [];
  for (let j = 0; j < n; j++) {
    if (!pivotSet.has(j)) {
      nonPivotCols.push(j);
    }
  }

  const kernelDim = nonPivotCols.length;

  if (kernelDim === 0) {
    // Trivial kernel - return 0 x n matrix
    return zero_matrix(ring, 0, n);
  }

  // Build kernel basis: one row for each non-pivot column
  const basisRows: R[][] = [];

  for (let k = 0; k < kernelDim; k++) {
    const col = nonPivotCols[k]!;
    const row: R[] = [];

    for (let j = 0; j < n; j++) {
      if (j === col) {
        // Put 1 in the non-pivot position
        row.push(ring.one());
      } else if (pivotSet.has(j)) {
        // Put negative of E[pivotRow][col] where pivotRow corresponds to column j
        const pivotIdx = pivotCols.indexOf(j);
        row.push(E.get(pivotIdx, col).neg() as R);
      } else {
        // Other non-pivot columns get 0
        row.push(ring.zero());
      }
    }

    basisRows.push(row);
  }

  // Create result matrix.  The vectors above form the 'pivot' basis
  // ('pivot-generic' in Sage); the default basis over a field is 'echelon'.
  const result = new Matrix(ring, kernelDim, n, basisRows);
  if (basis === 'pivot' || basis === 'computed') {
    return result;
  }
  return rref(result);
}

/**
 * Return a matrix whose rows are a basis for the left kernel of the matrix.
 *
 * The left kernel of A is the set of row vectors x such that xA = 0.
 * This is the same as the right kernel of A^T.
 *
 * @param matrix - The matrix
 * @returns A matrix whose rows span the left kernel
 * @see Reference: sage/matrix/matrix2.pyx:left_kernel_matrix
 */
export function left_kernel_matrix<R extends FieldElement>(
  matrix: Matrix<R>,
  options?: { basis?: 'default' | 'computed' | 'echelon' | 'pivot' }
): Matrix<R> {
  // Left kernel of A = right kernel of A^T
  return right_kernel_matrix(matrix.transpose(), options);
}

/**
 * Return the right kernel of the matrix as a free module.
 *
 * The right kernel of a matrix A is the set of column vectors x such that Ax = 0.
 * This returns a vector space (for fields) or free module (for PIDs) whose
 * elements are all such vectors.
 *
 * @param matrix - The matrix
 * @param options - Optional parameters
 * @param options.basis - Basis format: 'echelon' (default) or 'pivot'
 * @returns The right kernel as a free module/vector space
 * @see Reference: sage/matrix/matrix2.pyx:right_kernel
 *
 * @example
 * ```typescript
 * const F7 = GF(7n);
 * const A = MatrixSpace(F7, 2, 3).__call__([[1, 2, 3], [4, 5, 6]]);
 * const K = right_kernel(A);
 * // K is a vector space of degree 3 containing vectors x where A*x = 0
 * ```
 */
export function right_kernel<R extends FieldElement>(
  matrix: Matrix<R>,
  options?: { basis?: 'default' | 'computed' | 'echelon' | 'pivot' }
): FreeModuleGeneric {
  const ring = matrix.base_ring;
  const n = matrix.ncols;

  // Get the kernel matrix (rows are basis vectors of the kernel)
  const kernelMatrix = right_kernel_matrix(matrix, options);

  // Create the ambient space R^n
  const ambient = isFieldRing(ring) ? VectorSpace(ring, n) : FreeModule(ring, n);

  // If kernel is trivial (no rows), return zero submodule
  if (kernelMatrix.nrows === 0) {
    return ambient.span([]);
  }

  // Convert kernel matrix rows to free module elements
  const basisVectors: FreeModuleElement[] = [];
  for (let i = 0; i < kernelMatrix.nrows; i++) {
    const rowData = kernelMatrix.row(i);
    basisVectors.push(ambient.createElement(rowData));
  }

  // The default basis from right_kernel_matrix is echelonized; a 'pivot' or
  // 'computed' basis is not, so it must be echelonized by the span routine.
  const alreadyEchelonized =
    options?.basis === undefined || options.basis === 'default' || options.basis === 'echelon';
  return ambient.span(basisVectors, ring, { alreadyEchelonized });
}

/**
 * Return the left kernel of the matrix as a free module.
 *
 * The left kernel of a matrix A is the set of row vectors x such that xA = 0.
 * This is equivalent to the right kernel of A^T.
 *
 * @param matrix - The matrix
 * @param options - Optional parameters
 * @param options.basis - Basis format: 'echelon' (default) or 'pivot'
 * @returns The left kernel as a free module/vector space
 * @see Reference: sage/matrix/matrix2.pyx:left_kernel
 *
 * @example
 * ```typescript
 * const F7 = GF(7n);
 * const A = MatrixSpace(F7, 3, 2).__call__([[1, 2], [3, 4], [5, 6]]);
 * const K = left_kernel(A);
 * // K is a vector space of degree 3 containing vectors x where x*A = 0
 * ```
 */
export function left_kernel<R extends FieldElement>(
  matrix: Matrix<R>,
  options?: { basis?: 'default' | 'computed' | 'echelon' | 'pivot' }
): FreeModuleGeneric {
  // Left kernel of A = right kernel of A^T
  return right_kernel(matrix.transpose(), options);
}

/**
 * Alias for left_kernel.
 * @see Reference: sage/matrix/matrix2.pyx:kernel
 */
export const kernel = left_kernel;

/**
 * Return the image of the matrix (row space).
 *
 * The image of the homomorphism on rows defined by right multiplication
 * by this matrix is exactly the row space of the matrix.
 *
 * @param matrix - The matrix
 * @returns The image (row space) as a module
 * @see Reference: sage/matrix/matrix2.pyx:image
 *
 * @example
 * ```typescript
 * const F7 = GF(7n);
 * const A = MatrixSpace(F7, 3, 3).__call__([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
 * const img = image(A);
 * // img equals the row space of A
 * ```
 */
export function image<R extends FieldElement>(matrix: Matrix<R>): FreeModuleGeneric {
  return row_space(matrix);
}

/**
 * Return the row space of the matrix.
 *
 * The row space is the free module spanned by the rows of the matrix.
 * The basis returned is in echelon form (from the RREF of the matrix).
 *
 * @param matrix - The matrix
 * @param base_ring - Optional base ring (defaults to matrix's base ring)
 * @returns The row space as a free module/vector space
 * @see Reference: sage/matrix/matrix2.pyx:row_space
 *
 * @example
 * ```typescript
 * const F7 = GF(7n);
 * const A = MatrixSpace(F7, 3, 3).__call__([[1, 2, 3], [4, 5, 6], [5, 0, 2]]);
 * const V = row_space(A);
 * // V is a vector space of degree 3 spanned by the rows of A
 * ```
 */
export function row_space<R extends FieldElement>(
  matrix: Matrix<R>,
  base_ring?: CoefficientRing<R>
): FreeModuleGeneric {
  const ring = base_ring ?? matrix.base_ring;
  const n = matrix.ncols;

  // Create the ambient space R^n (degree is number of columns)
  const ambient = isFieldRing(ring) ? VectorSpace(ring, n) : FreeModule(ring, n);

  // Handle zero rows case
  if (matrix.nrows === 0) {
    return ambient.span([]);
  }

  // Compute the RREF to get an echelonized basis
  const E = rref(matrix);

  // The rank equals the number of non-zero rows in RREF
  const r = rank(matrix);

  // Handle zero matrix case
  if (r === 0) {
    return ambient.span([]);
  }

  // Extract the non-zero rows from RREF as the echelonized basis
  const basisVectors: FreeModuleElement[] = [];
  for (let i = 0; i < r; i++) {
    const rowData = E.row(i);
    basisVectors.push(ambient.createElement(rowData));
  }

  // Return the span with the echelonized basis
  return ambient.span(basisVectors, ring, { alreadyEchelonized: true });
}

/**
 * Return the column space of the matrix.
 *
 * The column space is the free module spanned by the columns of the matrix.
 * This is equivalent to the row space of the transpose.
 *
 * @param matrix - The matrix
 * @returns The column space as a free module/vector space
 * @see Reference: sage/matrix/matrix2.pyx:column_space
 *
 * @example
 * ```typescript
 * const F7 = GF(7n);
 * const A = MatrixSpace(F7, 3, 3).__call__([[1, 2, 3], [4, 5, 6], [5, 0, 2]]);
 * const V = column_space(A);
 * // V is a vector space of degree 3 spanned by the columns of A
 * ```
 */
export function column_space<R extends FieldElement>(matrix: Matrix<R>): FreeModuleGeneric {
  // Column space of A = row space of A^T
  return row_space(matrix.transpose());
}

// ============================================================================
// Linear System Solving
// ============================================================================

/**
 * Solve the system X * A = B for X.
 *
 * This is equivalent to solving A^T * X^T = B^T, i.e.,
 * X = (solve_right(A^T, B^T))^T.
 *
 * @param A - The matrix A
 * @param B - The right-hand side (matrix or vector)
 * @param check - Whether to check the solution
 * @param extend - Whether to extend the base ring if necessary
 * @returns The solution X
 * @see Reference: sage/matrix/matrix2.pyx:solve_left
 */
export function solve_left<R extends FieldElement>(
  A: Matrix<R>,
  B: Matrix<R>,
  check?: boolean,
  extend?: boolean
): Matrix<R> {
  // X * A = B is equivalent to A^T * X^T = B^T
  // So X^T = solve_right(A^T, B^T), hence X = solve_right(A^T, B^T)^T
  return solve_right(A.transpose(), B.transpose(), check, extend).transpose();
}

/**
 * Solve the system A * X = B for X.
 *
 * Uses Gaussian elimination (forward and back substitution) to solve the system.
 * If A is m x n and B is m x k, then X is n x k.
 *
 * @param A - The matrix A (m x n)
 * @param B - The right-hand side (m x k matrix)
 * @param check - Whether to check the solution (default: true)
 * @param extend - Whether to extend the base ring if necessary
 * @returns The solution X (n x k matrix)
 * @throws {ValueError} If the system has no solution
 * @see Reference: sage/matrix/matrix2.pyx:solve_right
 */
export function solve_right<R extends FieldElement>(
  A: Matrix<R>,
  B: Matrix<R>,
  check: boolean = true,
  extend?: boolean
): Matrix<R> {
  const m = A.nrows;
  const n = A.ncols;
  const k = B.ncols;
  const ring = A.base_ring;

  if (B.nrows !== m) {
    throw new ValueError(`incompatible dimensions: A is ${m}x${n}, B has ${B.nrows} rows`);
  }

  // Create augmented matrix [A | B]
  const aug: R[][] = [];
  for (let i = 0; i < m; i++) {
    aug.push([]);
    for (let j = 0; j < n; j++) {
      aug[i]!.push(A.get(i, j));
    }
    for (let j = 0; j < k; j++) {
      aug[i]!.push(B.get(i, j));
    }
  }

  // Forward elimination with partial pivoting
  const pivotCols: number[] = [];
  let pivotRow = 0;

  for (let col = 0; col < n && pivotRow < m; col++) {
    // Find pivot
    let found = -1;
    for (let i = pivotRow; i < m; i++) {
      if (!aug[i]![col]!.isZero()) {
        found = i;
        break;
      }
    }

    if (found === -1) {
      continue;
    }

    // Swap rows
    if (found !== pivotRow) {
      [aug[pivotRow], aug[found]] = [aug[found]!, aug[pivotRow]!];
    }

    pivotCols.push(col);

    // Scale pivot row
    const pivot = aug[pivotRow]![col]!;
    const pivotInv = getInverse(pivot);

    for (let j = col; j < n + k; j++) {
      aug[pivotRow]![j] = aug[pivotRow]![j]!.mul(pivotInv) as R;
    }

    // Eliminate below
    for (let i = pivotRow + 1; i < m; i++) {
      if (!aug[i]![col]!.isZero()) {
        const factor = aug[i]![col]!;
        for (let j = col; j < n + k; j++) {
          aug[i]![j] = aug[i]![j]!.sub(factor.mul(aug[pivotRow]![j]!) as R) as R;
        }
      }
    }

    pivotRow++;
  }

  // Check for inconsistency: if there's a row with zeros in A part but nonzero in B part
  for (let i = pivotRow; i < m; i++) {
    let allZeroA = true;
    for (let j = 0; j < n; j++) {
      if (!aug[i]![j]!.isZero()) {
        allZeroA = false;
        break;
      }
    }
    if (allZeroA) {
      for (let j = n; j < n + k; j++) {
        if (!aug[i]![j]!.isZero()) {
          throw new ValueError('matrix equation has no solutions');
        }
      }
    }
  }

  // Back substitution
  for (let i = pivotCols.length - 1; i >= 0; i--) {
    const col = pivotCols[i]!;
    // Eliminate above
    for (let row = 0; row < i; row++) {
      if (!aug[row]![col]!.isZero()) {
        const factor = aug[row]![col]!;
        for (let j = col; j < n + k; j++) {
          aug[row]![j] = aug[row]![j]!.sub(factor.mul(aug[i]![j]!) as R) as R;
        }
      }
    }
  }

  // Extract solution
  // The solution X has n rows and k columns
  // For each variable (row of X), check if it's a pivot variable or free variable
  const pivotSet = new Set(pivotCols);
  const X = zero_matrix(ring, n, k);

  for (let idx = 0; idx < pivotCols.length; idx++) {
    const col = pivotCols[idx]!;
    for (let j = 0; j < k; j++) {
      X.set(col, j, aug[idx]![n + j]!);
    }
  }

  // Verify solution if requested
  if (check) {
    const AX = A.mul(X);
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < k; j++) {
        if (!AX.get(i, j).eq(B.get(i, j))) {
          throw new ValueError('matrix equation has no solutions');
        }
      }
    }
  }

  return X;
}

// ============================================================================
// Characteristic and Minimal Polynomials
// ============================================================================

/**
 * Return the characteristic polynomial of the matrix.
 *
 * Uses the division-free algorithm from Algorithm 3.1 in [Sei2002].
 * This algorithm works over any commutative ring with identity.
 *
 * @param matrix - The square matrix
 * @param variable - Variable name (default: 'x')
 * @param algorithm - 'df' (division-free) or 'hessenberg' (default: 'df')
 * @returns The characteristic polynomial det(xI - A)
 * @see Reference: sage/matrix/matrix2.pyx:charpoly
 */
export function charpoly<R extends RingElement>(
  matrix: Matrix<R>,
  variable: string = 'x',
  algorithm: 'df' | 'hessenberg' = 'df'
): Polynomial<R> {
  if (!matrix.is_square()) {
    throw new ValueError('characteristic polynomial is only defined for square matrices');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;
  const polyRing = new PolynomialRing(ring, variable);

  // Corner case: empty matrix has characteristic polynomial 1
  if (n === 0) {
    return polyRing.one();
  }

  if (algorithm === 'hessenberg') {
    return _charpoly_hessenberg(matrix, variable);
  }

  // Use division-free algorithm (works over any ring)
  // Based on Algorithm 3.1 from [Sei2002]
  //
  // The algorithm computes coefficients F[0], ..., F[n-1] such that
  // charpoly = x^n + F[n-1]*x^{n-1} + ... + F[1]*x + F[0]

  // F[p] is the coefficient of x^{n-1-p} in the intermediate polynomial F_t
  const F: R[] = new Array(n);
  for (let i = 0; i < n; i++) {
    F[i] = ring.zero();
  }

  // a[p][i] represents vectors used in the computation
  // We only need one row of a at a time (the previous one)
  const a: R[][] = [];
  for (let p = 0; p < n; p++) {
    a.push(new Array(n).fill(ring.zero()));
  }

  // A[p] stores intermediate values
  const A: R[] = new Array(n);
  for (let i = 0; i < n; i++) {
    A[i] = ring.zero();
  }

  // F_0 = -M[0,0]
  F[0] = matrix.get(0, 0).neg() as R;

  for (let t = 1; t < n; t++) {
    // Set a[0] to be column t of M restricted to rows 0..t
    for (let i = 0; i <= t; i++) {
      a[0]![i] = matrix.get(i, t);
    }

    // Set A[0] to be the (t,t) entry
    A[0] = matrix.get(t, t);

    for (let p = 1; p < t; p++) {
      // Set a[p] to the product of M[0..t, 0..t] * a[p-1]
      for (let i = 0; i <= t; i++) {
        let s = ring.zero();
        for (let j = 0; j <= t; j++) {
          s = s.add(matrix.get(i, j).mul(a[p - 1]![j]!)) as R;
        }
        a[p]![i] = s;
      }

      // Set A[p] to be the t-th entry in a[p]
      A[p] = a[p]![t]!;
    }

    // Set A[t] to be M[t, 0..t] * a[t-1]
    let s = ring.zero();
    for (let j = 0; j <= t; j++) {
      s = s.add(matrix.get(t, j).mul(a[t - 1]![j]!)) as R;
    }
    A[t] = s;

    // Update F coefficients
    for (let p = 0; p <= t; p++) {
      s = F[p]!;
      for (let k = 0; k < p; k++) {
        s = s.sub(A[k]!.mul(F[p - k - 1]!)) as R;
      }
      F[p] = s.sub(A[p]!) as R;
    }
  }

  // Build polynomial: x^n + F[n-1]*x^{n-1} + ... + F[1]*x + F[0]
  // In the algorithm, F[p] is the coefficient of x^{n-1-p}
  // So F[0] is coeff of x^{n-1}, F[1] is coeff of x^{n-2}, ..., F[n-1] is coeff of x^0
  // We need coeffs[i] to be coefficient of x^i
  // coeffs[0] = F[n-1], coeffs[1] = F[n-2], ..., coeffs[n-1] = F[0], coeffs[n] = 1
  const coeffs: R[] = [];
  for (let i = n - 1; i >= 0; i--) {
    coeffs.push(F[i]!);
  }
  coeffs.push(ring.one()); // coefficient of x^n

  return new Polynomial<R>(coeffs, polyRing);
}

/**
 * Compute the characteristic polynomial via the Hessenberg form.
 *
 * Transforms the matrix to Hessenberg form (which preserves the characteristic
 * polynomial) and then applies the recursion of Cohen, *A Course in
 * Computational Algebraic Number Theory*, Algorithm 2.2.9.  Requires a field,
 * since the Hessenberg reduction divides.
 *
 * @see Reference: sage/matrix/matrix2.pyx:_charpoly_hessenberg
 */
function _charpoly_hessenberg<R extends RingElement>(
  matrix: Matrix<R>,
  variable: string
): Polynomial<R> {
  const n = matrix.nrows;
  const ring = matrix.base_ring;
  const polyRing = new PolynomialRing(ring, variable);

  const H = hessenberg_form(matrix as unknown as Matrix<R & FieldElement>) as unknown as Matrix<R>;

  // c is an (n+1) x (n+1) array whose rows hold the intermediate polynomials.
  const c: R[][] = [];
  for (let i = 0; i <= n; i++) {
    c.push(new Array<R>(n + 1).fill(ring.zero()));
  }
  c[0]![0] = ring.one();

  const addMultipleOfRow = (target: number, source: number, s: R): void => {
    for (let i = 0; i <= n; i++) {
      c[target]![i] = c[target]![i]!.add(s.mul(c[source]![i]!) as R) as R;
    }
  };

  for (let m = 1; m <= n; m++) {
    // Row m gets x * c[m-1] (i.e. c[m-1] shifted right) minus H[m-1,m-1]*c[m-1]
    for (let i = 1; i <= n; i++) {
      c[m]![i] = c[m - 1]![i - 1]!;
    }
    addMultipleOfRow(m, m - 1, H.get(m - 1, m - 1).neg() as R);

    let t = ring.one();
    for (let i = 1; i < m; i++) {
      t = t.mul(H.get(m - i, m - i - 1)) as R;
      addMultipleOfRow(m, m - i - 1, t.mul(H.get(m - i - 1, m - 1)).neg() as R);
    }
  }

  return new Polynomial<R>([...c[n]!], polyRing);
}

/**
 * Alias for charpoly.
 * @see Reference: sage/matrix/matrix2.pyx:characteristic_polynomial
 */
export const characteristic_polynomial = charpoly;

/**
 * Evaluate a univariate polynomial at a square matrix (Horner's rule).
 *
 * This is the matrix analogue of Sage's `f(A)` for a polynomial `f` over the
 * base ring of `A`; the constant term is multiplied by the identity matrix.
 */
function _evaluate_at_matrix<R extends RingElement>(
  poly: Polynomial<R>,
  matrix: Matrix<R>
): Matrix<R> {
  const ring = matrix.base_ring;
  const n = matrix.nrows;
  const deg = poly.degree();

  let result = zero_matrix(ring, n);
  if (deg < 0) {
    return result;
  }

  // Horner: result = ((c_d * A + c_{d-1}) * A + ...) * A + c_0
  result = identity_matrix(ring, n).scalar_mul(poly.getCoeff(deg));
  for (let i = deg - 1; i >= 0; i--) {
    result = result.mul(matrix).add(identity_matrix(ring, n).scalar_mul(poly.getCoeff(i)));
  }
  return result;
}

/**
 * Return the characteristic of a coefficient ring, or `undefined` when the ring
 * does not report one.
 */
function _ringCharacteristic<R extends RingElement>(ring: CoefficientRing<R>): bigint | undefined {
  // Some rings expose `characteristic` as a method (`QQ`), others as a plain
  // property (`GF(p)`, `Zmod(n)`).
  const anyRing = ring as unknown as {
    characteristic?: (() => bigint | number) | bigint | number;
  };
  const c = anyRing.characteristic;
  if (typeof c === 'function') {
    const v = c.call(ring);
    return typeof v === 'bigint' ? v : BigInt(v);
  }
  if (typeof c === 'bigint') {
    return c;
  }
  if (typeof c === 'number') {
    return BigInt(c);
  }
  return undefined;
}

/**
 * Test whether a polynomial over a field is squarefree.
 *
 * Mirrors `polynomial_element.pyx:is_squarefree`: a separable polynomial is
 * squarefree; in characteristic zero the converse holds, and in positive
 * characteristic we fall back on the factorization.
 */
function _is_squarefree<R extends RingElement>(f: Polynomial<R>): boolean {
  if (f.derivative().gcd(f).isConstant()) {
    return true;
  }
  const char = _ringCharacteristic(f.parent.base_ring);
  if (char === 0n) {
    return false;
  }
  return f.factor().every(([, e]) => e <= 1);
}

/**
 * Factor a polynomial, translating "no factorization here" into the message
 * that Sage's `minpoly` reports (`matrix2.pyx:3113`).
 */
function _factor_for_minpoly<R extends RingElement>(
  f: Polynomial<R>
): Array<[Polynomial<R>, number]> {
  try {
    return f.factor();
  } catch (e) {
    if (e instanceof NotImplementedError) {
      throw new NotImplementedError('minimal polynomial not implemented');
    }
    throw e;
  }
}

/**
 * Return the minimal polynomial of the matrix.
 *
 * The minimal polynomial is the monic polynomial m(x) of least degree such that
 * m(A) = 0. It divides the characteristic polynomial.
 *
 * ALGORITHM (`matrix2.pyx:3096-3128`): start from the characteristic
 * polynomial `f`.  If `f` is squarefree it *is* the minimal polynomial.
 * Otherwise the minimal polynomial is the radical of `f` multiplied by
 * `h^(n-1)` for each repeated irreducible factor `h` of `f`, where `n` is the
 * smallest exponent with `dim ker(h(A)^n) = e * deg(h)` (`e` the multiplicity
 * of `h` in `f`).
 *
 * @param matrix - The square matrix
 * @param variable - Variable name (default: 'x')
 * @returns The minimal polynomial
 * @see Reference: sage/matrix/matrix2.pyx:minpoly
 */
export function minpoly<R extends RingElement>(
  matrix: Matrix<R>,
  variable: string = 'x'
): Polynomial<R> {
  if (!matrix.is_square()) {
    throw new ValueError('minimal polynomial is only defined for square matrices');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;
  const polyRing = new PolynomialRing(ring, variable);

  if (n === 0) {
    return polyRing.one();
  }

  const f = charpoly(matrix, variable);

  // `f.is_squarefree()` is much cheaper than factoring, and when it holds f
  // *is* the minimal polynomial (matrix2.pyx:3099-3107).
  try {
    if (_is_squarefree(f)) {
      return f;
    }
  } catch (e) {
    if (!(e instanceof NotImplementedError)) {
      throw e;
    }
  }

  // Now we have to work harder: find the power of each irreducible factor of f
  // that divides the minimal polynomial.
  const factors = _factor_for_minpoly(f);

  // mp starts out as f.radical(), the product of the distinct irreducible
  // factors of f.
  let mp = polyRing.one();
  for (const [h] of factors) {
    if (h.degree() > 0) {
      mp = mp.mul(h);
    }
  }

  for (const [h, e] of factors) {
    if (e > 1) {
      // Find the power of B = h(A) whose kernel has dimension e*deg(h).
      const B = _evaluate_at_matrix(h, matrix);
      let C = B;
      let k = 1;
      while (_kernel_dimension(C) < e * h.degree()) {
        if (k === e - 1) {
          k += 1;
          break;
        }
        C = C.mul(B);
        k += 1;
      }
      mp = mp.mul(h.pow(k - 1));
    }
  }

  return mp;
}

/**
 * Dimension of the kernel of a square matrix (`n - rank`).
 */
function _kernel_dimension<R extends RingElement>(matrix: Matrix<R>): number {
  return matrix.ncols - rank(matrix as unknown as Matrix<R & FieldElement>);
}

/**
 * Alias for minpoly.
 * @see Reference: sage/matrix/matrix2.pyx:minimal_polynomial
 */
export const minimal_polynomial = minpoly;

// ============================================================================
// Eigenvalues and Eigenvectors
// ============================================================================

/**
 * Return the eigenvalues of the matrix, with multiplicity.
 *
 * Follows `_eigenvalues_sage` (`matrix2.pyx:7252`): the eigenvalues are the
 * roots of the characteristic polynomial.  With `extend = false` only the
 * roots lying in the base ring are returned.  With `extend = true` (Sage's
 * default) the roots are taken in the algebraic closure; we can only do that
 * when the characteristic polynomial already splits over the base ring, and
 * raise NotImplementedError otherwise, since algebraic closures are not
 * implemented in this port.
 *
 * @param matrix - The square matrix
 * @param extend - Whether to extend to the algebraic closure (default: true)
 * @param algorithm - Algorithm to use (currently only 'sage' supported)
 * @returns List of eigenvalues with multiplicities
 * @see Reference: sage/matrix/matrix2.pyx:eigenvalues
 * @see Reference: sage/matrix/matrix2.pyx:_eigenvalues_sage
 * @see Deviation: `extend = true` requires the charpoly to split over the
 *   base ring, because algebraic closures are not implemented.
 */
export function eigenvalues<R extends RingElement>(
  matrix: Matrix<R>,
  extend: boolean = true,
  algorithm?: string
): R[] {
  if (!matrix.is_square()) {
    throw new ValueError('eigenvalues is only defined for square matrices');
  }

  const n = matrix.nrows;
  if (n === 0) {
    return [];
  }

  // Compute characteristic polynomial
  const cp = charpoly(matrix);

  // Find roots of the characteristic polynomial using polynomial factorization.
  // The roots() method uses Berlekamp/Cantor-Zassenhaus for finite fields,
  // which is O(degree^3) instead of brute-force O(field_size * degree).
  // This mirrors SageMath's _eigenvalues_sage method (line 7279 in matrix2.pyx):
  //   return Sequence(r for r, m in self.charpoly().roots() for _ in range(m))
  let rootsWithMult: Array<[R, number]>;
  try {
    rootsWithMult = cp.roots();
  } catch (e) {
    // If roots() fails (e.g., unsupported ring type), provide informative error
    if (e instanceof NotImplementedError) {
      throw new NotImplementedError(
        `eigenvalues not implemented for matrices over ${matrix.base_ring}: ${e.message}`
      );
    }
    throw e;
  }

  // Expand roots by their multiplicities (SageMath returns eigenvalues with multiplicity)
  const eigenvalueList: R[] = [];
  for (const [root, mult] of rootsWithMult) {
    for (let i = 0; i < mult; i++) {
      eigenvalueList.push(root);
    }
  }

  if (extend && eigenvalueList.length < n) {
    throw new NotImplementedError(
      `algebraic closure is not implemented for ${matrix.base_ring}; ` +
        'the characteristic polynomial does not split over the base ring ' +
        '(pass extend = false to get only the eigenvalues in the base ring)'
    );
  }

  return eigenvalueList;
}

/**
 * Return the singular values of the matrix.
 *
 * The singular values are the square roots of the eigenvalues of A^T * A.
 * For finite fields, this concept doesn't apply directly as it requires
 * real/complex numbers and square roots.
 *
 * @param matrix - The matrix
 * @returns List of singular values
 * @see Reference: sage/matrix/matrix2.pyx:singular_values
 */
export function singular_values<R extends RingElement>(matrix: Matrix<R>): R[] {
  // Singular values require computing sqrt of eigenvalues of A^T * A
  // This requires a ring with square roots (like RDF or CDF)
  throw new NotImplementedError(
    'singular_values requires rings with sqrt (e.g., RDF, CDF). ' +
      'For finite fields, consider eigenvalues of A^T * A instead.'
  );
}

/**
 * Return the left eigenvectors of the matrix.
 *
 * Computes left eigenvectors by finding the null space of (A - lambda*I)^T
 * for each eigenvalue lambda.
 *
 * @param matrix - The square matrix
 * @param other - Optional matrix for generalized eigenvalue problem (not yet supported)
 * @param extend - Whether to extend to algebraic closure
 * @param algorithm - Algorithm to use
 * @returns List of (eigenvalue, eigenvectors, multiplicity) tuples
 * @see Reference: sage/matrix/matrix2.pyx:eigenvectors_left
 */
export function eigenvectors_left<R extends FieldElement>(
  matrix: Matrix<R>,
  other?: Matrix<R>,
  extend: boolean = true,
  algorithm?: string
): Array<[R, R[][], number]> {
  if (!matrix.is_square()) {
    throw new ValueError('eigenvectors_left is only defined for square matrices');
  }

  if (other !== undefined) {
    throw new NotImplementedError('generalized eigenvector decomposition is not yet implemented');
  }

  const n = matrix.nrows;
  if (n === 0) {
    return [];
  }

  // Get eigenvalues
  const eigenvalueList = eigenvalues(matrix, extend);

  // For each distinct eigenvalue, compute the left eigenspace
  // Left eigenvectors v satisfy v * A = lambda * v
  // which is equivalent to A^T * v^T = lambda * v^T
  // So we compute the right kernel of (A^T - lambda*I)

  const result: Array<[R, R[][], number]> = [];
  const ring = matrix.base_ring;
  const seenEigenvalues = new Set<string>();

  for (const lambda of eigenvalueList) {
    // Use string representation to check for duplicates
    const lambdaStr = String(lambda);
    if (seenEigenvalues.has(lambdaStr)) {
      continue;
    }
    seenEigenvalues.add(lambdaStr);

    // Compute A^T - lambda * I
    const AT = matrix.transpose();
    const I = identity_matrix(ring, n);
    const M = AT.sub(I.scalar_mul(lambda));

    // Compute right kernel of M (which is the left eigenspace of A for lambda)
    const kernelMatrix = right_kernel_matrix(M);

    // Extract basis vectors from kernel matrix rows
    const eigenvectors: R[][] = [];
    for (let i = 0; i < kernelMatrix.nrows; i++) {
      eigenvectors.push(kernelMatrix.row(i));
    }

    // Algebraic multiplicity: count occurrences of this eigenvalue
    let multiplicity = 0;
    for (const ev of eigenvalueList) {
      if (ev.eq(lambda)) {
        multiplicity++;
      }
    }

    if (eigenvectors.length > 0) {
      result.push([lambda, eigenvectors, multiplicity]);
    }
  }

  return result;
}

/**
 * Return the right eigenvectors of the matrix.
 *
 * Computes right eigenvectors by finding the null space of (A - lambda*I)
 * for each eigenvalue lambda.
 *
 * @param matrix - The square matrix
 * @param other - Optional matrix for generalized eigenvalue problem (not yet supported)
 * @param extend - Whether to extend to algebraic closure
 * @returns List of (eigenvalue, eigenvectors, multiplicity) tuples
 * @see Reference: sage/matrix/matrix2.pyx:eigenvectors_right
 */
export function eigenvectors_right<R extends FieldElement>(
  matrix: Matrix<R>,
  other?: Matrix<R>,
  extend: boolean = true
): Array<[R, R[][], number]> {
  if (!matrix.is_square()) {
    throw new ValueError('eigenvectors_right is only defined for square matrices');
  }

  if (other !== undefined) {
    throw new NotImplementedError('generalized eigenvector decomposition is not yet implemented');
  }

  const n = matrix.nrows;
  if (n === 0) {
    return [];
  }

  // Get eigenvalues
  const eigenvalueList = eigenvalues(matrix, extend);

  // For each distinct eigenvalue, compute the right eigenspace
  // Right eigenvectors v satisfy A * v = lambda * v
  // So we compute the right kernel of (A - lambda*I)

  const result: Array<[R, R[][], number]> = [];
  const ring = matrix.base_ring;
  const seenEigenvalues = new Set<string>();

  for (const lambda of eigenvalueList) {
    // Use string representation to check for duplicates
    const lambdaStr = String(lambda);
    if (seenEigenvalues.has(lambdaStr)) {
      continue;
    }
    seenEigenvalues.add(lambdaStr);

    // Compute A - lambda * I
    const I = identity_matrix(ring, n);
    const M = matrix.sub(I.scalar_mul(lambda));

    // Compute right kernel of M (which is the right eigenspace of A for lambda)
    const kernelMatrix = right_kernel_matrix(M);

    // Extract basis vectors from kernel matrix rows
    const eigenvectors: R[][] = [];
    for (let i = 0; i < kernelMatrix.nrows; i++) {
      eigenvectors.push(kernelMatrix.row(i));
    }

    // Algebraic multiplicity: count occurrences of this eigenvalue
    let multiplicity = 0;
    for (const ev of eigenvalueList) {
      if (ev.eq(lambda)) {
        multiplicity++;
      }
    }

    if (eigenvectors.length > 0) {
      result.push([lambda, eigenvectors, multiplicity]);
    }
  }

  return result;
}

/**
 * Eigenspace representation containing eigenvalue and basis vectors.
 */
export interface Eigenspace<R extends RingElement> {
  eigenvalue: R;
  basis: R[][];
  algebraicMultiplicity: number;
  geometricMultiplicity: number;
}

/**
 * Return the left eigenspaces of the matrix.
 *
 * The left eigenspace for an eigenvalue lambda consists of all vectors v
 * such that v * A = lambda * v.
 *
 * @param matrix - The square matrix
 * @param format - Output format ('all' or 'galois', currently only 'all' supported)
 * @param variable - Variable name for algebraic eigenvalues (not yet used)
 * @param algebraic_multiplicity - Whether to include algebraic multiplicities (default: true)
 * @returns Array of eigenspace objects containing eigenvalue, basis, and multiplicities
 * @see Reference: sage/matrix/matrix2.pyx:eigenspaces_left
 */
export function eigenspaces_left<R extends FieldElement>(
  matrix: Matrix<R>,
  format: 'all' | 'galois' = 'all',
  variable: string = 'a',
  algebraic_multiplicity: boolean = true
): Eigenspace<R>[] {
  if (!matrix.is_square()) {
    throw new ValueError('eigenspaces_left is only defined for square matrices');
  }

  // Get eigenvectors which already compute eigenspaces
  const eigenvectorData = eigenvectors_left(matrix);

  const result: Eigenspace<R>[] = [];
  for (const [eigenvalue, basis, algMult] of eigenvectorData) {
    result.push({
      eigenvalue,
      basis,
      algebraicMultiplicity: algMult,
      geometricMultiplicity: basis.length,
    });
  }

  return result;
}

/**
 * Return the right eigenspaces of the matrix.
 *
 * The right eigenspace for an eigenvalue lambda consists of all vectors v
 * such that A * v = lambda * v.
 *
 * @param matrix - The square matrix
 * @param format - Output format ('all' or 'galois', currently only 'all' supported)
 * @param variable - Variable name for algebraic eigenvalues (not yet used)
 * @param algebraic_multiplicity - Whether to include algebraic multiplicities (default: true)
 * @returns Array of eigenspace objects containing eigenvalue, basis, and multiplicities
 * @see Reference: sage/matrix/matrix2.pyx:eigenspaces_right
 */
export function eigenspaces_right<R extends FieldElement>(
  matrix: Matrix<R>,
  format: 'all' | 'galois' = 'all',
  variable: string = 'a',
  algebraic_multiplicity: boolean = true
): Eigenspace<R>[] {
  if (!matrix.is_square()) {
    throw new ValueError('eigenspaces_right is only defined for square matrices');
  }

  // Get eigenvectors which already compute eigenspaces
  const eigenvectorData = eigenvectors_right(matrix);

  const result: Eigenspace<R>[] = [];
  for (const [eigenvalue, basis, algMult] of eigenvectorData) {
    result.push({
      eigenvalue,
      basis,
      algebraicMultiplicity: algMult,
      geometricMultiplicity: basis.length,
    });
  }

  return result;
}

/**
 * Return the eigenmatrix on the left.
 *
 * Returns matrices D and P where D is a diagonal matrix of eigenvalues
 * and the rows of P are corresponding eigenvectors (or zero vectors if
 * the geometric multiplicity is less than the algebraic multiplicity).
 *
 * The matrices satisfy P * A = D * P.
 *
 * @param matrix - The square matrix A
 * @param other - Optional matrix for generalized problem (not yet supported)
 * @returns Pair (D, P) where D is diagonal and rows of P are eigenvectors
 * @see Reference: sage/matrix/matrix2.pyx:eigenmatrix_left
 */
export function eigenmatrix_left<R extends FieldElement>(
  matrix: Matrix<R>,
  other?: Matrix<R>
): [Matrix<R>, Matrix<R>] {
  if (!matrix.is_square()) {
    throw new ValueError('eigenmatrix_left is only defined for square matrices');
  }

  if (other !== undefined) {
    throw new NotImplementedError('generalized eigenmatrix decomposition is not yet implemented');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return [zero_matrix(ring, 0), zero_matrix(ring, 0)];
  }

  // Get eigenvectors
  const eigenvectorData = eigenvectors_left(matrix);

  // Build D (diagonal matrix of eigenvalues) and P (matrix of eigenvectors)
  const D = zero_matrix(ring, n);
  const P = zero_matrix(ring, n);

  let row = 0;
  for (const [eigenvalue, eigenvectors, algebraicMult] of eigenvectorData) {
    // Fill in eigenvalues and eigenvectors
    for (let i = 0; i < eigenvectors.length && row < n; i++) {
      D.set(row, row, eigenvalue);
      const vec = eigenvectors[i]!;
      for (let j = 0; j < n && j < vec.length; j++) {
        P.set(row, j, vec[j]!);
      }
      row++;
    }
    // If algebraic mult > geometric mult, fill with zeros (already done)
    // and repeat the eigenvalue
    for (let i = eigenvectors.length; i < algebraicMult && row < n; i++) {
      D.set(row, row, eigenvalue);
      // P row is already zero
      row++;
    }
  }

  return [D, P];
}

/**
 * Return the eigenmatrix on the right.
 *
 * Returns matrices D and P where D is a diagonal matrix of eigenvalues
 * and the columns of P are corresponding eigenvectors (or zero vectors if
 * the geometric multiplicity is less than the algebraic multiplicity).
 *
 * The matrices satisfy A * P = P * D.
 *
 * @param matrix - The square matrix A
 * @param other - Optional matrix for generalized problem (not yet supported)
 * @returns Pair (D, P) where D is diagonal and columns of P are eigenvectors
 * @see Reference: sage/matrix/matrix2.pyx:eigenmatrix_right
 */
export function eigenmatrix_right<R extends FieldElement>(
  matrix: Matrix<R>,
  other?: Matrix<R>
): [Matrix<R>, Matrix<R>] {
  if (!matrix.is_square()) {
    throw new ValueError('eigenmatrix_right is only defined for square matrices');
  }

  if (other !== undefined) {
    throw new NotImplementedError('generalized eigenmatrix decomposition is not yet implemented');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return [zero_matrix(ring, 0), zero_matrix(ring, 0)];
  }

  // Get eigenvectors
  const eigenvectorData = eigenvectors_right(matrix);

  // Build D (diagonal matrix of eigenvalues) and P (matrix of eigenvectors as columns)
  const D = zero_matrix(ring, n);
  const P = zero_matrix(ring, n);

  let col = 0;
  for (const [eigenvalue, eigenvectors, algebraicMult] of eigenvectorData) {
    // Fill in eigenvalues and eigenvectors
    for (let i = 0; i < eigenvectors.length && col < n; i++) {
      D.set(col, col, eigenvalue);
      const vec = eigenvectors[i]!;
      for (let j = 0; j < n && j < vec.length; j++) {
        P.set(j, col, vec[j]!);
      }
      col++;
    }
    // If algebraic mult > geometric mult, fill with zeros (already done)
    // and repeat the eigenvalue
    for (let i = eigenvectors.length; i < algebraicMult && col < n; i++) {
      D.set(col, col, eigenvalue);
      // P column is already zero
      col++;
    }
  }

  return [D, P];
}

// ============================================================================
// Matrix Properties
// ============================================================================

/**
 * Check if the matrix is the identity matrix.
 *
 * @param matrix - The matrix
 * @returns True if matrix is identity
 * @see Reference: sage/matrix/matrix2.pyx:is_one
 */
export function is_one<R extends RingElement>(matrix: Matrix<R>): boolean {
  if (!matrix.is_square()) {
    return false;
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;
  const one = ring.one();
  const zero = ring.zero();

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const entry = matrix.get(i, j);
      if (i === j) {
        if (!entry.eq(one)) return false;
      } else {
        if (!entry.eq(zero)) return false;
      }
    }
  }

  return true;
}

/**
 * Check if the matrix is a scalar matrix.
 *
 * A scalar matrix is a diagonal matrix with all diagonal entries equal.
 * If a is provided, checks if the matrix equals a*I.
 *
 * @param matrix - The matrix
 * @param a - Optional scalar to check against
 * @returns True if matrix is scalar (optionally equals a*I)
 * @see Reference: sage/matrix/matrix2.pyx:is_scalar
 */
export function is_scalar<R extends RingElement>(matrix: Matrix<R>, a?: R): boolean {
  if (!matrix.is_square()) {
    return false;
  }

  const n = matrix.nrows;
  if (n === 0) {
    return true;
  }

  const ring = matrix.base_ring;
  const zero = ring.zero();

  // Get the diagonal value (either from argument or from matrix)
  const diagonalValue = a !== undefined ? a : matrix.get(0, 0);

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const entry = matrix.get(i, j);
      if (i === j) {
        if (!entry.eq(diagonalValue)) return false;
      } else {
        if (!entry.eq(zero)) return false;
      }
    }
  }

  return true;
}

/**
 * Check if the matrix is symmetric.
 *
 * A matrix A is symmetric if A = A^T, i.e., A[i,j] = A[j,i] for all i, j.
 *
 * @param matrix - The matrix
 * @returns True if matrix is symmetric
 * @see Reference: sage/matrix/matrix2.pyx:is_symmetric
 */
export function is_symmetric<R extends RingElement>(matrix: Matrix<R>): boolean {
  if (!matrix.is_square()) {
    return false;
  }

  const n = matrix.nrows;

  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      if (!matrix.get(i, j).sub(matrix.get(j, i)).isZero()) {
        return false;
      }
    }
  }

  return true;
}

/**
 * Check if the matrix is diagonal.
 *
 * A diagonal matrix has all off-diagonal entries equal to zero.
 *
 * @param matrix - The matrix
 * @returns True if matrix is diagonal
 * @see Reference: sage/matrix/matrix2.pyx:is_diagonal
 */
export function is_diagonal<R extends RingElement>(matrix: Matrix<R>): boolean {
  const m = matrix.nrows;
  const n = matrix.ncols;

  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      if (i !== j && !matrix.get(i, j).isZero()) {
        return false;
      }
    }
  }

  return true;
}

/**
 * Check if the matrix is triangular.
 *
 * If side is 'upper', checks if all entries below the diagonal are zero.
 * If side is 'lower', checks if all entries above the diagonal are zero.
 * If side is not specified, checks if the matrix is either upper or lower triangular.
 *
 * @param matrix - The matrix
 * @param side - 'lower' or 'upper'
 * @returns True if matrix is triangular
 * @see Reference: sage/matrix/matrix2.pyx:is_triangular
 */
export function is_triangular<R extends RingElement>(
  matrix: Matrix<R>,
  side?: 'lower' | 'upper'
): boolean {
  const m = matrix.nrows;
  const n = matrix.ncols;

  function isUpperTriangular(): boolean {
    for (let i = 1; i < m; i++) {
      for (let j = 0; j < Math.min(i, n); j++) {
        if (!matrix.get(i, j).isZero()) {
          return false;
        }
      }
    }
    return true;
  }

  function isLowerTriangular(): boolean {
    for (let i = 0; i < m; i++) {
      for (let j = i + 1; j < n; j++) {
        if (!matrix.get(i, j).isZero()) {
          return false;
        }
      }
    }
    return true;
  }

  if (side === 'upper') {
    return isUpperTriangular();
  } else if (side === 'lower') {
    return isLowerTriangular();
  } else {
    return isUpperTriangular() || isLowerTriangular();
  }
}

/**
 * Check if the matrix is unitary.
 *
 * A matrix A is unitary if A * A^H = I, where A^H is the conjugate transpose.
 * For real matrices, this is the same as orthogonal (A * A^T = I).
 *
 * Note: This requires the ring to support conjugation. For rings without
 * conjugation (e.g., finite fields), this checks for orthogonality.
 *
 * @param matrix - The matrix
 * @returns True if matrix is unitary
 * @see Reference: sage/matrix/matrix2.pyx:is_unitary
 */
export function is_unitary<R extends FieldElement>(matrix: Matrix<R>): boolean {
  if (!matrix.is_square()) {
    return false;
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  // Check if A * A^T = I (for real matrices / no conjugation available)
  // A is unitary if A * A^H = I
  // For fields without conjugation, we just check A * A^T = I
  const At = matrix.transpose();

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let sum = ring.zero();
      for (let k = 0; k < n; k++) {
        sum = sum.add(matrix.get(i, k).mul(At.get(k, j))) as R;
      }
      // Check if diagonal or off-diagonal
      if (i === j) {
        if (!sum.sub(ring.one()).isZero()) {
          return false;
        }
      } else {
        if (!sum.isZero()) {
          return false;
        }
      }
    }
  }

  return true;
}

/**
 * Check if the matrix is normal (AA* = A*A).
 *
 * A matrix is normal if it commutes with its conjugate transpose.
 * For real matrices, this means AA^T = A^TA.
 *
 * @param matrix - The matrix
 * @returns True if matrix is normal
 * @see Reference: sage/matrix/matrix2.pyx:is_normal
 */
export function is_normal<R extends RingElement>(matrix: Matrix<R>): boolean {
  if (!matrix.is_square()) {
    return false;
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  // For real matrices / fields without conjugation: check AA^T = A^TA
  const At = matrix.transpose();

  // Compute A * A^T
  const AAt: R[][] = [];
  for (let i = 0; i < n; i++) {
    AAt.push([]);
    for (let j = 0; j < n; j++) {
      let sum = ring.zero();
      for (let k = 0; k < n; k++) {
        sum = sum.add(matrix.get(i, k).mul(At.get(k, j))) as R;
      }
      AAt[i]!.push(sum);
    }
  }

  // Compute A^T * A
  const AtA: R[][] = [];
  for (let i = 0; i < n; i++) {
    AtA.push([]);
    for (let j = 0; j < n; j++) {
      let sum = ring.zero();
      for (let k = 0; k < n; k++) {
        sum = sum.add(At.get(i, k).mul(matrix.get(k, j))) as R;
      }
      AtA[i]!.push(sum);
    }
  }

  // Check equality
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (!AAt[i]![j]!.sub(AtA[i]![j]!).isZero()) {
        return false;
      }
    }
  }

  return true;
}

/**
 * Check if the matrix is nilpotent.
 *
 * A matrix is nilpotent if there exists a positive integer k such that A^k = 0.
 * For an n x n matrix, this is equivalent to A^n = 0, which in turn is
 * equivalent to the characteristic polynomial being x^n.
 *
 * @param matrix - The square matrix
 * @returns True if matrix is nilpotent
 * @see Reference: sage/matrix/matrix2.pyx:is_nilpotent
 */
export function is_nilpotent<R extends RingElement>(matrix: Matrix<R>): boolean {
  if (!matrix.is_square()) {
    return false;
  }

  const n = matrix.nrows;
  if (n === 0) {
    return true; // Empty matrix is nilpotent
  }

  // Compute characteristic polynomial
  const cp = charpoly(matrix);

  // Matrix is nilpotent iff charpoly = x^n
  // i.e., all coefficients except the leading one (x^n) are zero
  for (let i = 0; i < n; i++) {
    if (!cp.coeffs[i]!.isZero()) {
      return false;
    }
  }

  return true;
}

// ============================================================================
// Base change
// ============================================================================

/**
 * Map a single element of `source` into `target`.
 *
 * Sage builds the new matrix with `M(self.list(), coerce=True)`
 * (`matrix0.pyx:1710`), i.e. it lets the coercion framework find the canonical
 * ring morphism.  This port has no coercion framework, so we first ask the
 * target ring to convert the element (`target.__call__`), and if that fails we
 * construct the canonical morphism ourselves in the only two situations where
 * one exists unconditionally:
 *
 * - `x` is a `Rational`: the (unique) morphism `QQ -> target` sends `n/d` to
 *   `target(n) * target(d)^-1`; it exists exactly when `target(d)` is a unit,
 *   and raises otherwise, just as Sage's `GF(7)(1/7)` does.
 * - `x` is an element of a quotient of `ZZ` (it carries an integral `value`):
 *   `Z/mZ -> target` is a well-defined ring map exactly when the
 *   characteristic of `target` divides `m`; `m = 0` (i.e. `ZZ`) always maps.
 *
 * Anything else raises, rather than guessing a map that may not exist.
 *
 * @throws {TypeError} when no canonical morphism is available
 */
function _coerce_entry<S extends RingElement, T extends RingElement>(
  x: S,
  source: CoefficientRing<S>,
  target: CoefficientRing<T>
): T {
  try {
    return target.__call__(x);
  } catch {
    // fall through to the canonical morphisms below
  }

  if (x instanceof Rational) {
    const [num, den] = x.asIntegerRatio();
    const n = target.__call__(num);
    const d = target.__call__(den);
    if (d.isZero()) {
      throw new ZeroDivisionError(`inverse of Mod(${den}, ...) does not exist`);
    }
    return n.mul(getInverse(d as T & FieldElement)) as T;
  }

  const value = (x as unknown as { value?: unknown }).value;
  if (typeof value === 'bigint') {
    const cs = _ringCharacteristic(source);
    const ct = _ringCharacteristic(target);
    // `Z/csZ -> target` is a ring map exactly when `ct` divides `cs`; `cs = 0`
    // is `ZZ`, which maps into everything.  An unknown characteristic on either
    // side means we cannot certify the map, so we refuse.
    if (cs !== undefined && ct !== undefined && (cs === 0n || (ct !== 0n && cs % ct === 0n))) {
      return target.__call__(value);
    }
  }

  throw new TypeError(`unable to coerce ${String(x)} into ${String(target)}`);
}

/**
 * Return the matrix obtained by coercing the entries of `matrix` into `ring`.
 *
 * Always returns a copy, as Sage does.
 *
 * @param matrix - The matrix
 * @param ring - The new base ring
 * @returns A new matrix over `ring`
 * @throws {TypeError} If some entry cannot be mapped into `ring`
 * @see Reference: sage/matrix/matrix0.pyx:change_ring
 * @see Deviation: Sage's coercion framework decides whether a morphism exists;
 *   this port has none, so `change_ring` only uses the target ring's own
 *   conversion plus the two canonical morphisms `QQ -> R` and `Z/mZ -> R`.
 */
export function change_ring<R extends RingElement, S extends RingElement>(
  matrix: Matrix<R>,
  ring: CoefficientRing<S>
): Matrix<S> {
  const source = matrix.base_ring;
  if ((ring as unknown) === (source as unknown)) {
    return matrix.copy() as unknown as Matrix<S>;
  }
  return new Matrix<S>(ring, matrix.nrows, matrix.ncols, (i, j) =>
    _coerce_entry(matrix.get(i, j), source, ring)
  );
}

/**
 * Check if the matrix is diagonalizable.
 *
 * ALGORITHM (`matrix2.pyx:12693-12722`): base-change to `base_field` when one
 * is given, then the algebraic multiplicities of the roots of the
 * characteristic polynomial must sum to `n` (i.e. the charpoly splits over the
 * base field), and for every eigenvalue the algebraic multiplicity must equal
 * the geometric multiplicity `dim ker(A - e)`.
 *
 * @param matrix - The matrix
 * @param base_field - Optional field to base-change to before testing
 * @returns True if matrix is diagonalizable
 * @throws {TypeError} If the matrix is not square (Sage: `not a square matrix`)
 * @throws {ValueError} If the base ring is not a field
 * @see Reference: sage/matrix/matrix2.pyx:is_diagonalizable
 */
export function is_diagonalizable<R extends FieldElement>(
  matrix: Matrix<R>,
  base_field?: CoefficientRing<R>
): boolean {
  if (!matrix.is_square()) {
    throw new TypeError('not a square matrix');
  }

  // `A = self.change_ring(base_field)` (matrix2.pyx:12695-12698)
  const A: Matrix<R> =
    base_field === undefined ? matrix : (change_ring(matrix, base_field) as Matrix<R>);

  const ring = A.base_ring;
  if (!isFieldRing(ring)) {
    throw new ValueError('matrix entries must be from a field');
  }

  const n = matrix.nrows;
  if (n === 0) {
    return true;
  }

  // Check that the sum of the algebraic multiplicities equals the number of rows
  const evals = charpoly(A).roots();
  let total = 0;
  for (const [, mult] of evals) {
    total += mult;
  }
  if (total < n) {
    return false;
  }

  // Check equality of algebraic and geometric multiplicity
  const I = identity_matrix(ring, n);
  for (const [e, am] of evals) {
    const gm = right_nullity(A.sub(I.scalar_mul(e)));
    if (am !== gm) {
      return false;
    }
  }

  return true;
}

/**
 * Check if the matrix is semisimple.
 *
 * A square matrix is semisimple if its minimal polynomial is squarefree,
 * equivalently if it is diagonalizable over the algebraic closure.
 *
 * @param matrix - The matrix
 * @returns True if matrix is semisimple
 * @see Reference: sage/matrix/matrix2.pyx:is_semisimple
 */
export function is_semisimple<R extends FieldElement>(matrix: Matrix<R>): boolean {
  // `self.minpoly().is_squarefree()` (matrix2.pyx:10635)
  return _is_squarefree(minpoly(matrix));
}

/**
 * Return the sign of a ring element that lies in the real numbers.
 *
 * Sage rejects base rings that cannot be seen as a subring of the real or
 * complex numbers (`matrix2.pyx:15735-15746`); we detect that by requiring the
 * elements to expose an ordering (`sign`, as `Rational`, `Integer` and
 * `RealNumber` all do).
 */
function _realSign(x: RingElement): number {
  const anyx = x as unknown as { sign?: unknown; real?: () => unknown; imag?: () => unknown };
  let s: unknown = anyx.sign;
  if (typeof s === 'function') {
    s = (s as () => unknown).call(x);
  }
  if (typeof s === 'bigint') {
    return s < 0n ? -1 : s > 0n ? 1 : 0;
  }
  if (typeof s === 'number') {
    return s < 0 ? -1 : s > 0 ? 1 : 0;
  }
  throw new ValueError(
    `Could not see ${String((x as { parent?: unknown }).parent ?? x)} as a subring of the real or complex numbers`
  );
}

/**
 * Check that the base ring can be seen as a subring of the reals/complexes.
 *
 * @throws {ValueError} mirroring `matrix2.pyx:15746`
 */
function _check_real_or_complex_subring<R extends RingElement>(ring: CoefficientRing<R>): void {
  try {
    _realSign(ring.one());
  } catch {
    throw new ValueError(
      `Could not see ${String(ring)} as a subring of the real or complex numbers`
    );
  }
}

/**
 * Return whether the matrix is Hermitian (equal to its conjugate transpose).
 *
 * @see Reference: sage/matrix/matrix2.pyx:is_hermitian
 */
export function is_hermitian<R extends RingElement>(matrix: Matrix<R>): boolean {
  if (!matrix.is_square()) {
    return false;
  }
  return matrix.eq(conjugate_transpose(matrix));
}

/**
 * Shared implementation of `is_positive_definite` and
 * `is_positive_semidefinite` (`matrix2.pyx:_is_positive_definite_or_semidefinite`).
 *
 * Sage reads the signs of the eigenvalues off the diagonal blocks of the
 * Bunch-Kaufman `block_ldlt` factorization.  We instead use the equivalent
 * exact criterion on the characteristic polynomial: writing
 * `charpoly(A) = sum_i c_i x^i`, the elementary symmetric functions of the
 * eigenvalues are `e_k = (-1)^k c_{n-k}`.  A Hermitian matrix has real
 * eigenvalues, and they are all `> 0` (resp. `>= 0`) exactly when every
 * `e_k` is `> 0` (resp. `>= 0`): if some eigenvalue were negative,
 * `p(-t) = (-1)^n sum_k e_k t^{n-k}` could not vanish for `t > 0`.
 *
 * @see Deviation: Sage uses `block_ldlt`; our `block_ldlt` is not available.
 */
function _is_positive_definite_or_semidefinite<R extends RingElement>(
  matrix: Matrix<R>,
  semi: boolean
): boolean {
  _check_real_or_complex_subring(matrix.base_ring);

  if (!is_hermitian(matrix)) {
    return false;
  }

  const n = matrix.nrows;
  if (n === 0) {
    return true; // vacuously
  }

  const p = charpoly(matrix);
  for (let k = 1; k <= n; k++) {
    // e_k = (-1)^k * c_{n-k}
    const c = p.getCoeff(n - k);
    const s = k % 2 === 0 ? _realSign(c) : -_realSign(c);
    if (semi ? s < 0 : s <= 0) {
      return false;
    }
  }

  return true;
}

/**
 * Check if the matrix is positive definite.
 *
 * By SageMath convention a positive-definite matrix must be real symmetric or
 * complex Hermitian, and the base ring must be a subring of the real or
 * complex numbers.
 *
 * @param matrix - The matrix
 * @returns True if matrix is positive definite
 * @throws {ValueError} If the base ring is not a subring of the reals/complexes
 * @see Reference: sage/matrix/matrix2.pyx:is_positive_definite
 */
export function is_positive_definite<R extends RingElement>(matrix: Matrix<R>): boolean {
  return _is_positive_definite_or_semidefinite(matrix, false);
}

/**
 * Check if the matrix is positive semidefinite.
 *
 * @param matrix - The matrix
 * @returns True if matrix is positive semidefinite
 * @throws {ValueError} If the base ring is not a subring of the reals/complexes
 * @see Reference: sage/matrix/matrix2.pyx:is_positive_semidefinite
 */
export function is_positive_semidefinite<R extends RingElement>(matrix: Matrix<R>): boolean {
  return _is_positive_definite_or_semidefinite(matrix, true);
}

/**
 * Return the similarity invariants of a square matrix over a field.
 *
 * For every monic irreducible factor `h` of the characteristic polynomial the
 * multiset of exponents of the elementary divisors `h^e` is recovered from the
 * dimensions of the kernels of `h(A)^k`: the number of elementary divisors
 * `h^e` with `e >= k` is `(dim ker h(A)^k - dim ker h(A)^(k-1)) / deg h`.
 *
 * These data determine (and are determined by) the rational canonical form, so
 * two matrices are similar exactly when their similarity invariants agree.
 *
 * @see Deviation: Sage compares `rational_form()`, which is not available here.
 */
function _similarity_invariants<R extends RingElement>(A: Matrix<R>): Map<string, number[]> {
  const invariants = new Map<string, number[]>();
  const f = charpoly(A);
  const factors = _factor_for_minpoly(f);

  for (const [h, e] of factors) {
    if (h.degree() <= 0) {
      continue;
    }
    const B = _evaluate_at_matrix(h, A);
    // ranks of B^k give the kernel dimensions
    const dims: number[] = [0];
    let C = identity_matrix(A.base_ring, A.nrows);
    for (let k = 1; k <= e; k++) {
      C = C.mul(B);
      const d = _kernel_dimension(C);
      dims.push(d);
      if (d === dims[k - 1]) {
        break;
      }
    }
    // countAtLeast[k] = number of elementary divisors h^j with j >= k
    const countAtLeast: number[] = [];
    for (let k = 1; k < dims.length; k++) {
      countAtLeast.push((dims[k]! - dims[k - 1]!) / h.degree());
    }
    // Turn into the multiset of exponents, sorted
    const exps: number[] = [];
    for (let k = 0; k < countAtLeast.length; k++) {
      const atLeastNext = k + 1 < countAtLeast.length ? countAtLeast[k + 1]! : 0;
      const exactly = countAtLeast[k]! - atLeastNext;
      for (let i = 0; i < exactly; i++) {
        exps.push(k + 1);
      }
    }
    exps.sort((a, b) => a - b);
    invariants.set(h.toString(), exps);
  }

  return invariants;
}

/**
 * Check if two matrices are similar.
 *
 * Two matrices A and B are similar if there exists an invertible P
 * such that `P^{-1} A P = B`.
 *
 * ALGORITHM (`matrix2.pyx:13047`): Sage compares the rational canonical forms
 * of the two matrices.  We compare the equivalent data — the elementary
 * divisors of the two matrices — because `rational_form` is not implemented.
 * Comparing only the characteristic and minimal polynomials is *not* enough
 * (the smallest counterexample is 6x6 with charpoly `x^6`, minpoly `x^3`
 * and block structures {3,3} vs {3,2,1}).
 *
 * @param A - First matrix
 * @param B - Second matrix
 * @param transformation - Whether to return the transformation matrix
 * @returns True if similar, or `[similar, P]` when `transformation` is set
 * @throws {TypeError} If `B` is not a matrix
 * @throws {ValueError} If the matrices are not square or have different sizes
 * @see Reference: sage/matrix/matrix2.pyx:is_similar
 */
export function is_similar<R extends FieldElement>(
  A: Matrix<R>,
  B: Matrix<R>,
  transformation?: boolean
): boolean | [boolean, Matrix<R> | null] {
  if (!(B instanceof Matrix)) {
    throw new TypeError(`similarity requires a matrix as an argument, not ${String(B)}`);
  }
  if (!A.is_square() || !B.is_square()) {
    throw new ValueError('similarity only makes sense for square matrices');
  }
  if (A.nrows !== B.nrows) {
    throw new ValueError('matrices do not have the same size');
  }

  const invA = _similarity_invariants(A);
  const invB = _similarity_invariants(B);

  let similar = invA.size === invB.size;
  if (similar) {
    for (const [h, expsA] of invA) {
      const expsB = invB.get(h);
      if (expsB === undefined || expsB.length !== expsA.length) {
        similar = false;
        break;
      }
      for (let i = 0; i < expsA.length; i++) {
        if (expsA[i] !== expsB[i]) {
          similar = false;
          break;
        }
      }
      if (!similar) {
        break;
      }
    }
  }

  if (!transformation) {
    return similar;
  }
  if (!similar) {
    return [false, null];
  }

  // Sage: "rational form routine does not provide transformation so if
  // possible, get transformations to Jordan form" (matrix2.pyx:13052-13057):
  //
  //     _, SA = A.jordan_form(transformation=True)
  //     _, SB = B.jordan_form(transformation=True)
  //     return (True, SB * SA.inverse())
  //
  // The returned `T` satisfies `A == T.inverse() * B * T` (the doctest at
  // matrix2.pyx:12831).
  try {
    const [, SA] = jordan_form(A, undefined, undefined, undefined, true) as [Matrix<R>, Matrix<R>];
    const [, SB] = jordan_form(B, undefined, undefined, undefined, true) as [Matrix<R>, Matrix<R>];
    const T = SB.mul(inverse(SA));
    if (_intertwines(T, A, B)) {
      return [true, T];
    }
  } catch {
    // Sage catches (ValueError, RuntimeError) here and moves to the algebraic
    // closure; we have none, so we fall through to the linear-algebra route.
  }

  const T = _intertwiner(A, B);
  if (T !== null) {
    return [true, T];
  }
  throw new ArithmeticError('unable to compute transformation for similar matrices');
}

/**
 * Test `T^-1 B T == A`, i.e. `B T == T A` with `T` invertible.
 */
function _intertwines<R extends FieldElement>(T: Matrix<R>, A: Matrix<R>, B: Matrix<R>): boolean {
  if (T.nrows !== A.nrows || T.ncols !== A.nrows) {
    return false;
  }
  if (!B.mul(T).sub(T.mul(A)).is_zero()) {
    return false;
  }
  return rank(T) === T.nrows;
}

/**
 * Return an invertible `T` with `T^-1 B T == A`, or `null`.
 *
 * Sage obtains the change of basis from the two Jordan forms, which requires
 * the eigenvalues to live in the base field (and raises `RuntimeError` when
 * they do not, see the `FiniteField(7^2)` example at `matrix2.pyx:12920`).
 * When that route is unavailable we solve the intertwining equation
 * `B X = X A` directly: it is a homogeneous linear system in the `n^2` entries
 * of `X`, whose solution space is a coset `X_0 * C` of the centralizer `C` of
 * `A` as soon as one invertible solution `X_0` exists — which is exactly the
 * case when `A` and `B` are similar.  A generic element of the solution space
 * is therefore invertible, so we search combinations of a kernel basis.  Every
 * candidate is verified (`B T == T A` and `rank T == n`) before it is returned,
 * so this can never produce a wrong transformation.
 *
 * @see Deviation: Sage raises `RuntimeError` when the Jordan form is
 *   unavailable; we return a (verified) transformation instead.
 */
function _intertwiner<R extends FieldElement>(A: Matrix<R>, B: Matrix<R>): Matrix<R> | null {
  const ring = A.base_ring;
  const n = A.nrows;
  if (n === 0) {
    return A.copy();
  }

  // Row (k, l), column (i, j) of the system B X - X A = 0, with X flattened
  // row-major:  coeff of x_{ij} is  [j == l] * B[k][i] - [i == k] * A[j][l].
  const N = n * n;
  const L = new Matrix<R>(ring, N, N, (row, col) => {
    const k = Math.floor(row / n);
    const l = row % n;
    const i = Math.floor(col / n);
    const j = col % n;
    let v = ring.zero();
    if (j === l) {
      v = v.add(B.get(k, i));
    }
    if (i === k) {
      v = v.sub(A.get(j, l));
    }
    return v;
  });

  const K = right_kernel_matrix(L);
  const d = K.nrows;
  if (d === 0) {
    return null;
  }

  const toMatrix = (coeffs: R[]): Matrix<R> =>
    new Matrix<R>(ring, n, n, (i, j) => {
      let v = ring.zero();
      for (let r = 0; r < d; r++) {
        v = v.add(coeffs[r]!.mul(K.get(r, i * n + j)));
      }
      return v;
    });

  // First try each basis vector on its own, then deterministic pseudo-random
  // combinations.  Over a field with q elements the proportion of invertible
  // elements of the centralizer algebra is at least prod_{i>=1}(1 - q^-i)
  // >= 0.288, so a handful of tries suffices with overwhelming probability.
  for (let r = 0; r < d; r++) {
    const coeffs = Array.from({ length: d }, (_, s) => (s === r ? ring.one() : ring.zero()));
    const T = toMatrix(coeffs);
    if (_intertwines(T, A, B)) {
      return T;
    }
  }

  let state = 0x9e3779b9;
  const nextInt = (): number => {
    state ^= state << 13;
    state >>>= 0;
    state ^= state >>> 17;
    state ^= state << 5;
    state >>>= 0;
    return state;
  };
  for (let attempt = 0; attempt < 200; attempt++) {
    const bound = 2 + Math.min(attempt, 60);
    const coeffs: R[] = [];
    for (let r = 0; r < d; r++) {
      let c = ring.zero();
      const k = nextInt() % bound;
      for (let t = 0; t < k; t++) {
        c = c.add(ring.one());
      }
      coeffs.push(c);
    }
    const T = toMatrix(coeffs);
    if (_intertwines(T, A, B)) {
      return T;
    }
  }

  return null;
}

// ============================================================================
// Matrix Norms
// ============================================================================

/**
 * Map a matrix entry to a double, mirroring Sage's `apply_map(abs, R=RDF)`.
 *
 * @throws {TypeError} If the entry has no absolute value / real embedding
 */
function _absToDouble(x: RingElement): number {
  const anyx = x as unknown as {
    abs?: () => unknown;
    toNumber?: () => number;
    valueOf?: () => unknown;
  };
  if (typeof anyx.abs !== 'function') {
    throw new TypeError(`bad operand type for abs(): '${String(x)}'`);
  }
  const a = anyx.abs() as { toNumber?: () => number; valueOf?: () => unknown };
  if (typeof a.toNumber === 'function') {
    return a.toNumber();
  }
  if (typeof a.valueOf === 'function') {
    const v = a.valueOf();
    if (typeof v === 'bigint') {
      return Number(v);
    }
    if (typeof v === 'number') {
      return v;
    }
  }
  throw new TypeError(`cannot convert ${String(x)} to a real number`);
}

/**
 * Convert a rational to the nearest double without overflowing `Number()`.
 *
 * `Rational.toNumber()` is `Number(num) / Number(den)`, which returns `NaN`
 * once both are past 2^1024; the iterates below routinely have denominators of
 * 2^120 and more, so we normalise the exponent first.
 */
function _rationalToDouble(r: Rational): number {
  const [num, den] = r.asIntegerRatio();
  if (num === 0n) {
    return 0;
  }
  const sign = num < 0n ? -1 : 1;
  const a = num < 0n ? -num : num;
  const bits = (x: bigint): number => x.toString(2).length;
  // Scale so that the integer quotient carries about 64 significant bits.
  const shift = 64 - (bits(a) - bits(den));
  const scaledNum = shift > 0 ? a << BigInt(shift) : a;
  const scaledDen = shift > 0 ? den : den << BigInt(-shift);
  const q = scaledNum / scaledDen;
  return sign * Number(q) * 2 ** -shift;
}

/**
 * Convert a matrix entry to an exact rational.
 *
 * @throws {TypeError} for base rings with no embedding into the complex
 *   numbers, mirroring the failure of Sage's `change_ring(CDF)`
 * @throws {NotImplementedError} for rings that do embed but are not exactly
 *   rational (`RR`, `CC`, number fields, ...), naming what is missing
 */
function _entryToRational(x: RingElement, ring: CoefficientRing<RingElement>): Rational {
  if (x instanceof Rational) {
    return x;
  }
  const char = _ringCharacteristic(ring);
  if (char !== undefined && char !== 0n) {
    // A finite field has no ring map to CDF, and this is Sage's own message.
    throw new TypeError(`no canonical coercion from ${String(ring)} to Complex Double Field`);
  }
  const anyx = x as unknown as { value?: unknown };
  if (typeof anyx.value === 'bigint') {
    return new Rational(anyx.value);
  }
  throw new NotImplementedError(
    `the spectral norm over ${String(ring)} is not implemented: it is computed by ` +
      'exactly isolating the largest root of the characteristic polynomial of A^H*A, ' +
      'which needs entries that are exact rationals (Sage instead runs a numerical SVD)'
  );
}

/** Evaluate a rational polynomial given by ascending coefficients. */
function _ratPolyEval(c: Rational[], x: Rational): Rational {
  let acc = new Rational(0n);
  for (let i = c.length - 1; i >= 0; i--) {
    acc = acc.mul(x).add(c[i]!);
  }
  return acc;
}

/** Derivative of a rational polynomial given by ascending coefficients. */
function _ratPolyDeriv(c: Rational[]): Rational[] {
  const out: Rational[] = [];
  for (let i = 1; i < c.length; i++) {
    out.push(c[i]!.mul(new Rational(BigInt(i))));
  }
  return _ratPolyTrim(out);
}

function _ratPolyTrim(c: Rational[]): Rational[] {
  const out = c.slice();
  while (out.length > 0 && out[out.length - 1]!.isZero()) {
    out.pop();
  }
  return out;
}

/** Remainder of `a` modulo `b` over QQ. */
function _ratPolyRem(a: Rational[], b: Rational[]): Rational[] {
  let r = _ratPolyTrim(a);
  const db = b.length - 1;
  const lb = b[db]!;
  while (r.length - 1 >= db && r.length > 0) {
    const shift = r.length - 1 - db;
    const factor = r[r.length - 1]!.div(lb);
    for (let i = 0; i <= db; i++) {
      r[i + shift] = r[i + shift]!.sub(factor.mul(b[i]!));
    }
    r = _ratPolyTrim(r);
  }
  return r;
}

/** Monic gcd over QQ. */
function _ratPolyGcd(a: Rational[], b: Rational[]): Rational[] {
  let x = _ratPolyTrim(a);
  let y = _ratPolyTrim(b);
  while (y.length > 0) {
    const r = _ratPolyRem(x, y);
    x = y;
    y = r;
  }
  if (x.length === 0) {
    return x;
  }
  const lc = x[x.length - 1]!;
  return x.map((c) => c.div(lc));
}

/** Exact quotient `a / b` over QQ (`b` divides `a`). */
function _ratPolyQuo(a: Rational[], b: Rational[]): Rational[] {
  const r = _ratPolyTrim(a).slice();
  const db = b.length - 1;
  const lb = b[db]!;
  const q: Rational[] = Array.from({ length: Math.max(0, r.length - db) }, () => new Rational(0n));
  for (let deg = r.length - 1; deg >= db; deg--) {
    const factor = r[deg]!.div(lb);
    if (factor.isZero()) {
      continue;
    }
    q[deg - db] = factor;
    for (let i = 0; i <= db; i++) {
      r[i + deg - db] = r[i + deg - db]!.sub(factor.mul(b[i]!));
    }
  }
  return _ratPolyTrim(q);
}

/**
 * Largest root of a monic rational polynomial all of whose roots are real and
 * non-negative (the characteristic polynomial of a positive semidefinite
 * Hermitian matrix).
 *
 * The distinct roots are the roots of the squarefree part `g = p / gcd(p,p')`,
 * which are simple, so Newton's method started to the right of the largest root
 * decreases monotonically to it and converges quadratically.  Every iterate is
 * rounded *up* to a multiple of `eps = trace / 2^160`, which keeps it above the
 * root (so monotonicity and the sign of `g'` are preserved) while bounding the
 * denominators.  Since `trace/n <= lambda_max <= trace`, that is a relative
 * accuracy of at most `n * 2^-160`, far below double precision, and `trace` is
 * itself the starting point of the iteration.
 */
function _largestRootPSD(p: Rational[]): Rational {
  const deg = p.length - 1;
  const zero = new Rational(0n);
  if (deg <= 0) {
    return zero;
  }
  // trace = -a_{n-1}; all eigenvalues are >= 0, so 0 <= lambda_max <= trace.
  const trace = p[deg - 1]!.neg();
  if (trace.cmp(zero) <= 0) {
    return zero; // every eigenvalue is 0
  }

  const g = _ratPolyQuo(p, _ratPolyGcd(p, _ratPolyDeriv(p)));
  const gp = _ratPolyDeriv(g);

  const eps = trace.div(new Rational(1n << 160n));
  const snapUp = (x: Rational): Rational => {
    const k = x.div(eps).ceil();
    return eps.mul(new Rational(k));
  };

  let x = trace;
  for (let iter = 0; iter < 10000; iter++) {
    const gx = _ratPolyEval(g, x);
    if (gx.isZero()) {
      return x;
    }
    const gpx = _ratPolyEval(gp, x);
    if (gpx.isZero()) {
      throw new ArithmeticError('Newton iteration for the spectral norm hit a critical point');
    }
    const next = snapUp(x.sub(gx.div(gpx)));
    if (next.cmp(x) >= 0) {
      return x;
    }
    x = next;
  }
  throw new ArithmeticError('Newton iteration for the spectral norm did not converge');
}

/**
 * Return the p-norm of the matrix, as a double (Sage returns an `RDF` number).
 *
 * - `1` -- the largest column-sum of the absolute values
 * - `2` -- (default) the Euclidean / spectral norm, i.e. the largest singular
 *   value
 * - `Infinity` -- the largest row-sum of the absolute values
 * - `'frob'` -- the Frobenius norm, `sqrt(sum of squares)`
 *
 * @param matrix - The matrix
 * @param p - The norm type: 1, 2, Infinity or 'frob' (default: 2)
 * @returns The norm as a JavaScript double
 * @throws {TypeError} If the entries have no absolute value (e.g. finite fields)
 * @see Reference: sage/matrix/matrix2.pyx:norm
 * @see Deviation: for `p = 2` Sage runs a numerical SVD of `A^H A` over `CDF`
 *   and returns `sqrt(max singular value)`.  There is no SVD in this port, so
 *   we compute the same number exactly: the singular values of the Hermitian
 *   positive semidefinite matrix `A^H A` are its eigenvalues, and we isolate
 *   the largest root of its characteristic polynomial in exact rational
 *   arithmetic before rounding to a double.  This requires exact rational
 *   entries; other base rings raise `TypeError`, as Sage's `change_ring(CDF)`
 *   does for e.g. finite fields.
 */
export function norm<R extends RingElement>(matrix: Matrix<R>, p: number | 'frob' = 2): number {
  const m = matrix.nrows;
  const n = matrix.ncols;

  if (m === 0 || n === 0) {
    return 0;
  }

  if (p === 2) {
    // matrix2.pyx:16466-16471:
    //     A = self.dense_matrix().change_ring(CDF)
    //     A = A.conjugate_transpose() * A
    //     S = A.SVD()[1]
    //     return max(S.list()).real().sqrt()
    const ring = matrix.base_ring as CoefficientRing<RingElement>;
    // Fail fast (and with Sage's error) for base rings with no embedding.
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < n; j++) {
        _entryToRational(matrix.get(i, j), ring);
      }
    }
    const M = conjugate_transpose(matrix).mul(matrix);
    const cp = charpoly(M);
    const coeffs: Rational[] = [];
    for (let i = 0; i <= cp.degree(); i++) {
      coeffs.push(_entryToRational(cp.getCoeff(i), ring));
    }
    return Math.sqrt(_rationalToDouble(_largestRootPSD(coeffs)));
  }

  // A = self.apply_map(abs, R=RDF)
  const A: number[][] = [];
  for (let i = 0; i < m; i++) {
    A.push([]);
    for (let j = 0; j < n; j++) {
      A[i]!.push(_absToDouble(matrix.get(i, j)));
    }
  }

  if (p === 1) {
    // largest column-sum
    let best = Number.NEGATIVE_INFINITY;
    for (let j = 0; j < n; j++) {
      let colSum = 0;
      for (let i = 0; i < m; i++) {
        colSum += A[i]![j]!;
      }
      if (colSum > best) {
        best = colSum;
      }
    }
    return best;
  }

  if (p === Number.POSITIVE_INFINITY) {
    // largest row-sum
    let best = Number.NEGATIVE_INFINITY;
    for (let i = 0; i < m; i++) {
      let rowSum = 0;
      for (let j = 0; j < n; j++) {
        rowSum += A[i]![j]!;
      }
      if (rowSum > best) {
        best = rowSum;
      }
    }
    return best;
  }

  if (p === 'frob') {
    let sumSquares = 0;
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < n; j++) {
        sumSquares += A[i]![j]! * A[i]![j]!;
      }
    }
    return Math.sqrt(sumSquares);
  }

  throw new NotImplementedError(`norm with p=${p} is not implemented`);
}

/**
 * Return the density of non-zero entries in the matrix.
 *
 * The density is the exact ratio of the number of nonzero positions to
 * `nrows * ncols`.  Sage returns a rational number (`matrix2.pyx:10772`), and
 * so do we.
 *
 * @param matrix - The matrix
 * @returns The density as an exact rational between 0 and 1
 * @see Reference: sage/matrix/matrix2.pyx:density
 */
export function density<R extends RingElement>(matrix: Matrix<R>): Rational {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const total = m * n;

  if (total === 0) {
    return new Rational(0n, 1n);
  }

  let nonzeroCount = 0;
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      if (!matrix.get(i, j).isZero()) {
        nonzeroCount++;
      }
    }
  }

  return new Rational(BigInt(nonzeroCount), BigInt(total));
}

// ============================================================================
// Conjugate and Hermitian
// ============================================================================

/**
 * Return the conjugate of the matrix.
 *
 * Applies complex conjugation to each entry. For real rings or fields
 * without conjugation (like finite fields), this returns a copy of the
 * original matrix.
 *
 * @param matrix - The matrix
 * @returns The conjugate matrix
 * @see Reference: sage/matrix/matrix2.pyx:conjugate
 */
export function conjugate<R extends RingElement>(matrix: Matrix<R>): Matrix<R> {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  const result: R[][] = [];
  for (let i = 0; i < m; i++) {
    result.push([]);
    for (let j = 0; j < n; j++) {
      const entry = matrix.get(i, j);
      // Try to call conjugate method if it exists
      if (typeof (entry as unknown as { conjugate?: () => R }).conjugate === 'function') {
        result[i]!.push((entry as unknown as { conjugate: () => R }).conjugate());
      } else {
        // For real fields and finite fields, conjugate is identity
        result[i]!.push(entry);
      }
    }
  }

  return new Matrix(ring, m, n, result);
}

/**
 * Return the conjugate transpose (Hermitian transpose) of the matrix.
 *
 * This is the transpose of the matrix after applying complex conjugation
 * to each entry. Also known as the Hermitian adjoint or dagger.
 *
 * For real matrices, this is the same as the regular transpose.
 *
 * @param matrix - The matrix
 * @returns The conjugate transpose
 * @see Reference: sage/matrix/matrix2.pyx:conjugate_transpose
 */
export function conjugate_transpose<R extends RingElement>(matrix: Matrix<R>): Matrix<R> {
  return conjugate(matrix).transpose();
}

/**
 * Alias for conjugate_transpose.
 * @see Reference: sage/matrix/matrix2.pyx:H
 */
export const H = conjugate_transpose;

/**
 * Return the transpose of the matrix.
 *
 * @param matrix - The matrix
 * @returns The transpose
 * @see Reference: sage/matrix/matrix2.pyx:T
 */
export function T<R extends RingElement>(matrix: Matrix<R>): Matrix<R> {
  return matrix.transpose();
}

/**
 * Return the conjugate of the matrix.
 *
 * This is an alias for the conjugate function.
 *
 * @param matrix - The matrix
 * @returns The conjugate
 * @see Reference: sage/matrix/matrix2.pyx:C
 */
export function C<R extends RingElement>(matrix: Matrix<R>): Matrix<R> {
  return conjugate(matrix);
}

// ============================================================================
// Matrix Construction Operations
// ============================================================================

/**
 * Return the direct sum of two matrices.
 *
 * The direct sum of A (m x n) and B (p x q) is the block diagonal matrix:
 *   [ A  0 ]
 *   [ 0  B ]
 * of size (m + p) x (n + q).
 *
 * @param A - First matrix (m x n)
 * @param B - Second matrix (p x q)
 * @param subdivide - Whether to add subdivisions (default: false)
 * @returns The direct sum (m + p) x (n + q)
 * @see Reference: sage/matrix/constructor.py:block_diagonal_matrix
 */
export function direct_sum<R extends RingElement>(
  A: Matrix<R>,
  B: Matrix<R>,
  subdivide: boolean = false
): Matrix<R> {
  const m = A.nrows;
  const n = A.ncols;
  const p = B.nrows;
  const q = B.ncols;
  const ring = A.base_ring;

  // Create result matrix of size (m + p) x (n + q)
  const result = zero_matrix(ring, m + p, n + q);

  // Copy A into upper-left block
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      result.set(i, j, A.get(i, j));
    }
  }

  // Copy B into lower-right block
  for (let i = 0; i < p; i++) {
    for (let j = 0; j < q; j++) {
      result.set(m + i, n + j, B.get(i, j));
    }
  }

  return result;
}

/**
 * Stack matrices vertically.
 *
 * Returns a new matrix formed by appending the matrix (or matrices) bottom
 * below self.
 *
 * @param top - The top matrix
 * @param bottom - The bottom matrix (or matrices)
 * @param subdivide - Whether to add subdivisions (default: false)
 * @returns The stacked matrix
 * @see Reference: sage/matrix/matrix1.pyx:stack
 */
export function stack<R extends RingElement>(
  top: Matrix<R>,
  bottom: Matrix<R> | Matrix<R>[],
  subdivide: boolean = false
): Matrix<R> {
  const matrices = Array.isArray(bottom) ? [top, ...bottom] : [top, bottom];
  const ring = top.base_ring;

  // Check column compatibility
  const ncols = matrices[0]!.ncols;
  for (let i = 1; i < matrices.length; i++) {
    if (matrices[i]!.ncols !== ncols) {
      throw new ArithmeticError(
        `number of columns must be the same, not ${ncols} and ${matrices[i]!.ncols}`
      );
    }
  }

  // Calculate total rows
  let totalRows = 0;
  for (const M of matrices) {
    totalRows += M.nrows;
  }

  // Create result matrix
  const result = zero_matrix(ring, totalRows, ncols);

  // Copy matrices
  let currentRow = 0;
  for (const M of matrices) {
    for (let i = 0; i < M.nrows; i++) {
      for (let j = 0; j < ncols; j++) {
        result.set(currentRow + i, j, M.get(i, j));
      }
    }
    currentRow += M.nrows;
  }

  return result;
}

/**
 * Augment matrices horizontally.
 *
 * Returns a new matrix formed by appending the matrix (or matrices) right
 * on the right side of self.
 *
 * @param left - The left matrix
 * @param right - The right matrix (or matrices)
 * @param subdivide - Whether to add subdivisions (default: false)
 * @returns The augmented matrix
 * @see Reference: sage/matrix/matrix1.pyx:augment
 */
export function augment<R extends RingElement>(
  left: Matrix<R>,
  right: Matrix<R> | Matrix<R>[],
  subdivide: boolean = false
): Matrix<R> {
  const matrices = Array.isArray(right) ? [left, ...right] : [left, right];
  const ring = left.base_ring;

  // Check row compatibility
  const nrows = matrices[0]!.nrows;
  for (let i = 1; i < matrices.length; i++) {
    if (matrices[i]!.nrows !== nrows) {
      throw new ArithmeticError(
        `number of rows must be the same, ${nrows} != ${matrices[i]!.nrows}`
      );
    }
  }

  // Calculate total columns
  let totalCols = 0;
  for (const M of matrices) {
    totalCols += M.ncols;
  }

  // Create result matrix
  const result = zero_matrix(ring, nrows, totalCols);

  // Copy matrices
  let currentCol = 0;
  for (const M of matrices) {
    for (let i = 0; i < nrows; i++) {
      for (let j = 0; j < M.ncols; j++) {
        result.set(i, currentCol + j, M.get(i, j));
      }
    }
    currentCol += M.ncols;
  }

  return result;
}
