/**
 * @module sage/coding/goppa_code
 * @description Goppa codes - classical algebraic codes used in McEliece cryptosystem
 *
 * Port of: sage/coding/goppa_code.py
 *
 * Goppa codes are a generalization of narrow-sense BCH codes.
 * They are defined by a generating polynomial g(x) over a finite field GF(q^m),
 * and a defining set L of elements from GF(q^m) that are not roots of g.
 *
 * For binary Goppa codes (used in McEliece cryptosystem):
 * - The minimum distance is >= 2t + 1 when g is squarefree, where t = deg(g)
 *   (note that `distance_bound()` reports Sage's bound, 1 + deg(g))
 * - They can correct up to t errors
 *
 * Goppa codes are used in the McEliece post-quantum cryptosystem because
 * the structure of the code (the Goppa polynomial) is hidden by a random
 * permutation, making it hard to distinguish from a random code.
 *
 * EXAMPLES:
 * ```typescript
 * const F = GFExtended(64);  // GF(2^6)
 * const R = new PolynomialRing(F, 'x');
 * const x = R.gen();
 * const g = x.pow(9).add(R.one());  // g(x) = x^9 + 1
 * const L = [...F].filter(a => !g.evaluate(a).isZero());  // support
 * const C = new GoppaCode(g, L);
 * // C is a [55, 16] Goppa code over GF(2)
 * ```
 */

import { NotImplementedError, ValueError, ZeroDivisionError } from '../errors.js';
import { Matrix, identity_matrix, zero_matrix } from '../matrix/matrix_generic.js';
import type { CoefficientRing, RingElement } from '../rings/polynomial/polynomial_element.js';
import type { Polynomial } from '../rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';

/**
 * Interface for field elements that support division/inverse.
 */
export interface FieldElement extends RingElement {
  div(other: this): this;
  inv(): this;
  pow(n: number | bigint): this;
  isZero(): boolean;
  isOne?(): boolean;
}

/**
 * Interface for a finite field.
 */
export interface FiniteField<E extends FieldElement> extends CoefficientRing<E> {
  characteristic: bigint;
  order: bigint;
  degree?: number;
  primitiveElement?(): E;
  random_element?(): E;
  [Symbol.iterator]?(): Iterator<E>;
}

/**
 * Custom error for decoding failures.
 */
export class DecodingError extends Error {
  override name = 'DecodingError';

  constructor(message: string = 'decoding failed') {
    super(message);
  }
}

/**
 * Convert a finite field element to a column vector over the prime field.
 *
 * For an element in GF(p^m), returns a column vector of m elements in GF(p)
 * representing the element in terms of the polynomial basis.
 *
 * @param element - Element of GF(p^m)
 * @param degree - Extension degree m
 * @param primeField - The prime field GF(p)
 * @returns Column vector representation
 */
function columnize<E extends FieldElement, P extends FieldElement>(
  element: E,
  degree: number,
  primeField: CoefficientRing<P>
): P[] {
  // Try to get coefficients if the element has a lift/coefficients method
  if (
    'coefficients' in element &&
    typeof (element as { coefficients: () => unknown[] }).coefficients === 'function'
  ) {
    const coeffs = (element as { coefficients: () => unknown[] }).coefficients();
    const result: P[] = [];
    for (let i = 0; i < degree; i++) {
      if (i < coeffs.length && coeffs[i] !== undefined) {
        // Convert coefficient to prime field element
        result.push(primeField.__call__(coeffs[i]) as P);
      } else {
        result.push(primeField.zero() as P);
      }
    }
    return result;
  }

  // Try lift.coeffs for Polynomial-based representation
  if ('lift' in element) {
    const lift = (element as { lift: { coeffs: readonly unknown[] } }).lift;
    if (lift && 'coeffs' in lift) {
      const coeffs = lift.coeffs;
      const result: P[] = [];
      for (let i = 0; i < degree; i++) {
        if (i < coeffs.length && coeffs[i] !== undefined) {
          result.push(primeField.__call__(coeffs[i]) as P);
        } else {
          result.push(primeField.zero() as P);
        }
      }
      return result;
    }
  }

  // Fallback: try to interpret as a single coefficient (prime field element)
  if ('value' in element) {
    const result: P[] = [primeField.__call__((element as { value: unknown }).value) as P];
    for (let i = 1; i < degree; i++) {
      result.push(primeField.zero() as P);
    }
    return result;
  }

  // Last resort: treat as constant
  const result: P[] = [primeField.__call__(element) as P];
  for (let i = 1; i < degree; i++) {
    result.push(primeField.zero() as P);
  }
  return result;
}

/**
 * Implementation of classical Goppa codes.
 *
 * A Goppa code is defined by:
 * - A generating (Goppa) polynomial g(x) over GF(q^m) of degree t
 * - A support set L = {L_0, L_1, ..., L_{n-1}} of elements from GF(q^m)
 *   such that g(L_i) != 0 for all i
 *
 * The parity check matrix H has entries:
 *   H[i][j] = coefficient of x^i in (g(L[j]))^{-1} mod g(x)
 *
 * For binary Goppa codes with squarefree g:
 * - Minimum distance >= 2t + 1
 * - Can correct up to t errors
 * - Used in McEliece cryptosystem
 *
 * @example
 * ```typescript
 * // Create a binary Goppa code over GF(2^6)
 * const F = GFExtended(64);  // GF(2^6)
 * const R = new PolynomialRing(F, 'x');
 * const x = R.gen();
 * const g = x.pow(9).add(R.one());  // g(x) = x^9 + 1, degree 9
 *
 * // Support: all elements where g(a) != 0
 * const L = [...F].filter(a => !g.evaluate(a).isZero());
 *
 * const C = new GoppaCode(g, L);
 * console.log(C.length());      // 55
 * console.log(C.dimension());   // 16 (approximately n - m*t = 55 - 6*9 = 1, but actual dimension may vary)
 * console.log(C.distance_bound()); // 10 (= 1 + deg g, as in SageMath)
 * ```
 */
export class GoppaCode<E extends FieldElement, P extends FieldElement = FieldElement> {
  private readonly _generating_pol: Polynomial<E>;
  private readonly _defining_set: E[];
  private readonly _extension_field: FiniteField<E>;
  private readonly _base_field: CoefficientRing<P>;
  private readonly _length: number;
  private readonly _extension_degree: number;

  private _parity_check_matrix: Matrix<P> | null = null;
  private _generator_matrix: Matrix<P> | null = null;
  private _dimension: number | null = null;

  /**
   * Create a Goppa code.
   *
   * @param generating_pol - Monic Goppa polynomial g(x) over GF(q^m)
   * @param defining_set - Support set L of elements from GF(q^m) where g(a) != 0
   * @throws {ValueError} If polynomial is not monic
   * @throws {ValueError} If any element of defining_set is a root of g
   */
  constructor(generating_pol: Polynomial<E>, defining_set: E[]) {
    // Verify polynomial is monic
    if (!generating_pol.is_monic()) {
      throw new ValueError('generating polynomial must be monic');
    }

    // Get the extension field from the polynomial's base ring
    this._extension_field = generating_pol.parent.base_ring as FiniteField<E>;
    this._generating_pol = generating_pol;
    this._defining_set = [...defining_set];
    this._length = defining_set.length;

    // Determine the prime subfield
    // The code is defined over GF(p), the prime subfield
    const p = this._extension_field.characteristic;
    this._extension_degree = this._extension_field.degree ?? this._computeExtensionDegree();

    // Create base field (prime field)
    // We'll use the same operations but restrict to the prime subfield
    this._base_field = this._createPrimeField(p) as CoefficientRing<P>;

    // Verify that g(a) != 0 for all a in defining_set
    for (const a of defining_set) {
      if (generating_pol.evaluate(a).isZero()) {
        throw new ValueError('defining elements cannot be roots of generating polynomial');
      }
    }
  }

  /**
   * Compute the extension degree from the field order.
   */
  private _computeExtensionDegree(): number {
    const q = this._extension_field.order;
    const p = this._extension_field.characteristic;

    // q = p^m, solve for m
    let m = 1;
    let power = p;
    while (power < q) {
      power *= p;
      m++;
    }
    return m;
  }

  /**
   * Create a prime field wrapper.
   */
  private _createPrimeField(p: bigint): CoefficientRing<P> {
    // Return a simple wrapper that uses the extension field but restricts to prime elements
    return {
      zero: () => this._extension_field.__call__(0) as unknown as P,
      one: () => this._extension_field.__call__(1) as unknown as P,
      __call__: (x: unknown) => {
        if (typeof x === 'number' || typeof x === 'bigint') {
          return this._extension_field.__call__(x) as unknown as P;
        }
        // Try to extract value from the element
        if (x && typeof x === 'object' && 'value' in x) {
          return this._extension_field.__call__((x as { value: unknown }).value) as unknown as P;
        }
        return this._extension_field.__call__(x) as unknown as P;
      },
    };
  }

  /**
   * Return the length of the code (number of codeword symbols).
   *
   * The length equals |L|, the size of the support set.
   */
  length(): number {
    return this._length;
  }

  /**
   * Return the Goppa polynomial.
   */
  goppa_polynomial(): Polynomial<E> {
    return this._generating_pol;
  }

  /**
   * Return the support (defining set) L.
   */
  support(): readonly E[] {
    return this._defining_set;
  }

  /**
   * Return the base field (prime field GF(p)) over which the code is defined.
   */
  base_field(): CoefficientRing<P> {
    return this._base_field;
  }

  /**
   * Return the extension field GF(q^m).
   */
  extension_field(): FiniteField<E> {
    return this._extension_field;
  }

  /**
   * Return the degree of the Goppa polynomial.
   */
  degree(): number {
    return this._generating_pol.degree();
  }

  /**
   * Compute the parity check matrix H.
   *
   * The parity check matrix is constructed as follows:
   * - For each element L[j] in the support, compute h[j] = 1/g(L[j])
   * - H[i][j] = h[j] * L[j]^i, expanded to the prime field
   *
   * The resulting matrix has dimensions (d * m) x n, where:
   * - d = degree of g
   * - m = extension degree
   * - n = length of code
   *
   * @returns The parity check matrix over the prime field
   */
  parity_check_matrix(): Matrix<P> {
    if (this._parity_check_matrix !== null) {
      return this._parity_check_matrix;
    }

    const g = this._generating_pol;
    const n = this._length;
    const d = g.degree();
    const m = this._extension_degree;
    const D = this._defining_set;

    // Compute h[i] = 1/g(D[i]) for each support element
    const h: E[] = [];
    for (let i = 0; i < n; i++) {
      h.push(g.evaluate(D[i]!).inv());
    }

    // Build the parity check matrix
    // It has d*m rows and n columns
    const rows: P[][] = [];

    for (let t = 0; t < d; t++) {
      // For row t, we compute h[i] * D[i]^t for each i
      // Then expand each result to m prime field elements
      for (let rowIdx = 0; rowIdx < m; rowIdx++) {
        const row: P[] = [];
        for (let i = 0; i < n; i++) {
          // Compute h[i] * D[i]^t
          let val: E;
          if (t === 0) {
            val = h[i]!;
          } else {
            val = h[i]!.mul(D[i]!.pow(t)) as E;
          }

          // Extract the rowIdx-th coefficient when expanded to prime field
          const coeffs = columnize<E, P>(val, m, this._base_field);
          row.push(coeffs[rowIdx]!);
        }
        rows.push(row);
      }
    }

    this._parity_check_matrix = new Matrix<P>(this._base_field, d * m, n, rows);
    return this._parity_check_matrix;
  }

  /**
   * Compute the dimension of the code.
   *
   * The dimension is n - rank(H), where H is the parity check matrix.
   * For a well-designed Goppa code, this is approximately n - m*t.
   */
  dimension(): number {
    if (this._dimension !== null) {
      return this._dimension;
    }

    const H = this.parity_check_matrix();
    const r = this._computeRank(H);
    this._dimension = this._length - r;
    return this._dimension;
  }

  /**
   * Compute the rank of a matrix using Gaussian elimination.
   */
  private _computeRank(matrix: Matrix<P>): number {
    const m = matrix.nrows;
    const n = matrix.ncols;

    // Make a copy of the matrix
    const M: P[][] = [];
    for (let i = 0; i < m; i++) {
      const row: P[] = [];
      for (let j = 0; j < n; j++) {
        row.push(matrix.get(i, j));
      }
      M.push(row);
    }

    let rank = 0;
    let col = 0;

    for (let row = 0; row < m && col < n; col++) {
      // Find pivot
      let pivotRow = -1;
      for (let i = row; i < m; i++) {
        if (!M[i]![col]!.isZero()) {
          pivotRow = i;
          break;
        }
      }

      if (pivotRow === -1) {
        continue;
      }

      // Swap rows
      [M[row], M[pivotRow]] = [M[pivotRow]!, M[row]!];

      const pivot = M[row]![col]!;
      const pivotInv = (pivot as FieldElement).inv() as P;

      // Scale pivot row
      for (let j = col; j < n; j++) {
        M[row]![j] = M[row]![j]!.mul(pivotInv) as P;
      }

      // Eliminate below and above
      for (let i = 0; i < m; i++) {
        if (i !== row && !M[i]![col]!.isZero()) {
          const factor = M[i]![col]!;
          for (let j = col; j < n; j++) {
            M[i]![j] = M[i]![j]!.sub(factor.mul(M[row]![j]!) as P) as P;
          }
        }
      }

      rank++;
      row++;
    }

    return rank;
  }

  /**
   * Compute the generator matrix from the parity check matrix.
   *
   * Uses Gaussian elimination to find the null space of H.
   *
   * @returns Generator matrix G such that G * H^T = 0
   */
  generator_matrix(): Matrix<P> {
    if (this._generator_matrix !== null) {
      return this._generator_matrix;
    }

    const H = this.parity_check_matrix();
    const n = this._length;
    const k = this.dimension();

    if (k === 0) {
      // Degenerate case: trivial code
      this._generator_matrix = zero_matrix(this._base_field, 0, n);
      return this._generator_matrix;
    }

    // Compute the right kernel of H (null space of H)
    // This gives us the generator matrix
    const G = this._computeKernel(H, k);
    this._generator_matrix = G;
    return this._generator_matrix;
  }

  /**
   * Compute a basis for the kernel of matrix H.
   */
  private _computeKernel(H: Matrix<P>, expectedDim: number): Matrix<P> {
    const m = H.nrows;
    const n = H.ncols;

    // Compute RREF and find pivot columns
    const M: P[][] = [];
    for (let i = 0; i < m; i++) {
      const row: P[] = [];
      for (let j = 0; j < n; j++) {
        row.push(H.get(i, j));
      }
      M.push(row);
    }

    const pivotCols: number[] = [];
    let row = 0;

    for (let col = 0; col < n && row < m; col++) {
      // Find pivot
      let pivotRow = -1;
      for (let i = row; i < m; i++) {
        if (!M[i]![col]!.isZero()) {
          pivotRow = i;
          break;
        }
      }

      if (pivotRow === -1) {
        continue;
      }

      pivotCols.push(col);

      // Swap rows
      [M[row], M[pivotRow]] = [M[pivotRow]!, M[row]!];

      const pivot = M[row]![col]!;
      const pivotInv = (pivot as FieldElement).inv() as P;

      // Scale pivot row
      for (let j = 0; j < n; j++) {
        M[row]![j] = M[row]![j]!.mul(pivotInv) as P;
      }

      // Eliminate above and below
      for (let i = 0; i < m; i++) {
        if (i !== row && !M[i]![col]!.isZero()) {
          const factor = M[i]![col]!;
          for (let j = 0; j < n; j++) {
            M[i]![j] = M[i]![j]!.sub(factor.mul(M[row]![j]!) as P) as P;
          }
        }
      }

      row++;
    }

    // Find non-pivot columns (free variables)
    const pivotSet = new Set(pivotCols);
    const freeCols: number[] = [];
    for (let j = 0; j < n; j++) {
      if (!pivotSet.has(j)) {
        freeCols.push(j);
      }
    }

    // Build kernel basis
    const kernelRows: P[][] = [];
    for (const freeCol of freeCols) {
      const vec: P[] = [];
      for (let j = 0; j < n; j++) {
        if (j === freeCol) {
          vec.push(this._base_field.one() as P);
        } else if (pivotSet.has(j)) {
          // Find which pivot row corresponds to this column
          const pivotIdx = pivotCols.indexOf(j);
          if (pivotIdx !== -1 && pivotIdx < M.length) {
            vec.push(M[pivotIdx]![freeCol]!.neg() as P);
          } else {
            vec.push(this._base_field.zero() as P);
          }
        } else {
          vec.push(this._base_field.zero() as P);
        }
      }
      kernelRows.push(vec);
    }

    return new Matrix<P>(this._base_field, kernelRows.length, n, kernelRows);
  }

  /**
   * Return a lower bound for the minimum distance of the code.
   *
   * Computed using the degree of the generating polynomial of ``self``.
   * The minimum distance is guaranteed to be bigger than or equal to this
   * bound.
   *
   * Port of `goppa_code.py:293-314`: SageMath always returns
   * `1 + deg(g)`, in every characteristic.  The `2t + 1` bound for binary
   * Goppa codes is *not* what Sage reports (its own `[8, 2]` doctest over
   * `GF(2^3)` with `g = x^2 + x + 1` gives `3`), and it is not even a valid
   * bound when `g` is not squarefree.
   *
   * @returns 1 + deg(g)
   */
  distance_bound(): number {
    return 1 + this._generating_pol.degree();
  }

  /**
   * Return the number of errors {@link decode} is guaranteed to correct.
   *
   * This has no SageMath counterpart (`sage.coding.goppa_code.GoppaCode`
   * registers no decoder at all); it reports the radius of the decoder
   * implemented here:
   *
   * - characteristic 2 with squarefree `g`: Patterson's algorithm corrects
   *   `deg(g)` errors, because `d >= 2 deg(g) + 1` for a squarefree binary
   *   Goppa polynomial;
   * - otherwise the key-equation decoder corrects
   *   `floor((distance_bound() - 1) / 2) = floor(deg(g) / 2)` errors.
   *
   * @see Deviation: port-only API, not a SageMath method
   */
  error_correction_capability(): number {
    const t = this._generating_pol.degree();
    if (this._extension_field.characteristic === 2n && this._isSquarefree(this._generating_pol)) {
      return t;
    }
    return Math.floor((this.distance_bound() - 1) / 2);
  }

  /**
   * Test whether a polynomial over the extension field is squarefree.
   */
  private _isSquarefree(f: Polynomial<E>): boolean {
    const fp = f.derivative();
    if (fp.isZero()) {
      return f.degree() === 0;
    }
    return f.gcd(fp).degree() === 0;
  }

  /**
   * Encode a message vector.
   *
   * The message is multiplied by the generator matrix to produce a codeword.
   *
   * @param message - Message vector of length k (dimension of code)
   * @returns Codeword vector of length n
   * @throws {ValueError} If message length is wrong
   */
  encode(message: P[]): P[] {
    const k = this.dimension();
    if (message.length !== k) {
      throw new ValueError(`message length must be ${k}, got ${message.length}`);
    }

    const G = this.generator_matrix();
    const n = this._length;

    // Compute message * G
    const codeword: P[] = [];
    for (let j = 0; j < n; j++) {
      let sum = this._base_field.zero() as P;
      for (let i = 0; i < k; i++) {
        sum = sum.add(message[i]!.mul(G.get(i, j)) as P) as P;
      }
      codeword.push(sum);
    }

    return codeword;
  }

  /**
   * Check if a vector is a valid codeword.
   *
   * A vector c is a codeword if H * c^T = 0.
   *
   * @param word - Vector to check
   * @returns True if word is a valid codeword
   */
  is_codeword(word: P[]): boolean {
    if (word.length !== this._length) {
      return false;
    }

    const H = this.parity_check_matrix();

    // Check H * word^T = 0
    for (let i = 0; i < H.nrows; i++) {
      let sum = this._base_field.zero() as P;
      for (let j = 0; j < this._length; j++) {
        sum = sum.add(H.get(i, j).mul(word[j]!) as P) as P;
      }
      if (!sum.isZero()) {
        return false;
      }
    }

    return true;
  }

  /**
   * Compute the syndrome of a received word.
   *
   * The syndrome S(x) = sum_i (r_i / (x - L_i)) mod g(x)
   * where r is the received word and L is the support.
   *
   * @param received - Received word
   * @returns Syndrome polynomial
   */
  syndrome(received: P[]): Polynomial<E> {
    if (received.length !== this._length) {
      throw new ValueError(`received length must be ${this._length}, got ${received.length}`);
    }

    const g = this._generating_pol;
    const R = g.parent;
    const L = this._defining_set;

    // Compute syndrome S(x) = sum_i r_i * h_i(x)
    // where h_i(x) = (product_{j != i} (x - L_j)) / g(L_i)
    // Equivalently: S(x) = sum_i r_i / (x - L_i) mod g(x)

    let S = R.zero();

    for (let i = 0; i < this._length; i++) {
      if (received[i]!.isZero()) {
        continue;
      }

      // Compute 1 / (x - L_i) mod g(x)
      const xMinusLi = R.gen().sub(R.__call__(L[i]!));

      // Extended Euclidean algorithm to find inverse mod g
      const inv = this._inverseMod(xMinusLi, g);

      // Scale by received[i] (need to convert from prime field to extension)
      const ri = this._extension_field.__call__(received[i]) as E;
      const term = inv.scalar_mul(ri);

      S = S.add(term);
    }

    // Reduce mod g
    if (S.degree() >= g.degree()) {
      const [_q, r] = S.quo_rem(g);
      S = r;
    }

    return S;
  }

  /**
   * Compute polynomial inverse modulo another polynomial.
   */
  private _inverseMod(a: Polynomial<E>, mod: Polynomial<E>): Polynomial<E> {
    const [g, s, _t] = a.xgcd(mod);

    // g should be constant if a and mod are coprime
    if (g.degree() !== 0) {
      throw new ValueError('polynomial is not invertible modulo the given modulus');
    }

    // s * a + t * mod = g (constant)
    // So s * a = g (mod mod)
    // Inverse is s / g
    const gInv = g.getCoeff(0).inv();
    return s.scalar_mul(gInv);
  }

  /**
   * Decode a received word using Patterson's algorithm (for binary Goppa codes).
   *
   * Patterson's algorithm:
   * 1. Compute syndrome S(x)
   * 2. If S(x) = 0, return received word (no errors)
   * 3. Compute T(x) = sqrt(x + S(x)^{-1}) mod g(x)
   * 4. Use extended Euclidean algorithm to find a(x), b(x) with deg(a) <= t/2, deg(b) <= (t-1)/2
   * 5. Error locator polynomial sigma(x) = a(x)^2 + x*b(x)^2
   * 6. Find roots of sigma(x) in GF(2^m) - these give error positions
   *
   * @param received - Received word (possibly with errors)
   * @returns Decoded codeword
   * @throws {DecodingError} If decoding fails
   */
  decode(received: P[]): P[] {
    if (received.length !== this._length) {
      throw new ValueError(`received length must be ${this._length}, got ${received.length}`);
    }

    // Check if already a codeword
    if (this.is_codeword(received)) {
      return [...received];
    }

    const p = this._extension_field.characteristic;

    if (p === 2n) {
      return this._patterson_decode(received);
    } else {
      // For non-binary codes, use a more general algorithm
      return this._general_decode(received);
    }
  }

  /**
   * Patterson's algorithm for binary Goppa codes.
   */
  private _patterson_decode(received: P[]): P[] {
    const g = this._generating_pol;
    const R = g.parent;
    const t = g.degree();
    const L = this._defining_set;

    // Step 1: Compute syndrome
    const S = this.syndrome(received);

    // If syndrome is zero, no errors
    if (S.isZero()) {
      return [...received];
    }

    // Step 2: Compute S^{-1} mod g
    const SInv = this._inverseMod(S, g);

    // Step 3: Compute T(x) = sqrt(x + S^{-1}(x)) mod g(x)
    // For binary fields, need to find square root
    const xPlusSInv = R.gen().add(SInv);
    const T = this._squareRootMod(xPlusSInv, g);

    // Step 4: Extended Euclidean algorithm to find a(x), b(x)
    // We want: a(x) + b(x) * T(x) = 0 mod g(x)
    // with deg(a) <= t/2 and deg(b) <= (t-1)/2

    const [a, b] = this._partialXGCD(g, T, Math.floor(t / 2));

    // Step 5: Error locator polynomial
    // sigma(x) = a(x)^2 + x * b(x)^2
    const aSq = a.mul(a);
    const bSq = b.mul(b);
    const x = R.gen();
    const sigma = aSq.add(x.mul(bSq));

    // Step 6: Find roots of sigma
    const errorPositions = this._findRoots(sigma, L);

    // Step 7: Flip the error positions
    const corrected = [...received];
    for (const pos of errorPositions) {
      // Flip bit: 0 -> 1, 1 -> 0
      const one = this._base_field.one() as P;
      corrected[pos] = corrected[pos]!.eq(one) ? (this._base_field.zero() as P) : one;
    }

    // Verify correction
    if (!this.is_codeword(corrected)) {
      throw new DecodingError('decoding failed: corrected word is not a codeword');
    }

    return corrected;
  }

  /**
   * Compute square root of a polynomial modulo g(x) in characteristic 2.
   *
   * In GF(2^m), sqrt(a) = a^{2^{m-1}}.
   * For a polynomial f(x) mod g(x), we compute f(x)^{2^{degree(g)*m - 1}}.
   */
  private _squareRootMod(f: Polynomial<E>, g: Polynomial<E>): Polynomial<E> {
    const m = this._extension_degree;
    const d = g.degree();

    // In GF(2^m)[x]/<g(x)>, the square root of f is f^{2^{dm-1}}
    // This is because x^2 = x^{2^{dm}} in this ring
    const exponent = 2n ** BigInt(d * m - 1);

    return this._powerMod(f, exponent, g);
  }

  /**
   * Compute polynomial power modulo another polynomial.
   */
  private _powerMod(base: Polynomial<E>, exp: bigint, mod: Polynomial<E>): Polynomial<E> {
    if (exp === 0n) {
      return mod.parent.one();
    }

    let result = mod.parent.one();
    let b = base.mod(mod);

    while (exp > 0n) {
      if ((exp & 1n) === 1n) {
        result = result.mul(b).mod(mod);
      }
      b = b.mul(b).mod(mod);
      exp >>= 1n;
    }

    return result;
  }

  /**
   * Partial extended GCD to find a, b with specific degree bounds.
   *
   * Returns (a, b) such that a + b*T = 0 mod g
   * with deg(a) <= bound and deg(b) <= bound - 1.
   */
  private _partialXGCD(
    g: Polynomial<E>,
    T: Polynomial<E>,
    bound: number
  ): [Polynomial<E>, Polynomial<E>] {
    const R = g.parent;

    let r0 = g;
    let r1 = T;
    let s0 = R.one();
    let s1 = R.zero();
    let t0 = R.zero();
    let t1 = R.one();

    // Stop as soon as `deg(r1) <= bound`: Patterson's algorithm needs
    // `deg(a) <= floor(t/2)`, and the loop must not run the extra step that
    // would push the degree strictly below the bound (that step throws away
    // the balanced (a, b) pair and yields a wrong error locator).
    while (r1.degree() > bound) {
      const [q, remainder] = r0.quo_rem(r1);

      const tempR = r1;
      r1 = remainder;
      r0 = tempR;

      const tempS = s1;
      s1 = s0.sub(q.mul(s1));
      s0 = tempS;

      const tempT = t1;
      t1 = t0.sub(q.mul(t1));
      t0 = tempT;
    }

    // a = r1, b = t1
    return [r1, t1];
  }

  /**
   * Find roots of a polynomial among the support elements.
   */
  private _findRoots(poly: Polynomial<E>, support: E[]): number[] {
    const positions: number[] = [];

    for (let i = 0; i < support.length; i++) {
      if (poly.evaluate(support[i]!).isZero()) {
        positions.push(i);
      }
    }

    return positions;
  }

  /**
   * Key-equation ("Sugiyama") decoding for Goppa codes in any characteristic.
   *
   * Writing `E` for the set of error positions, `sigma(x) = prod_{i in E}
   * (x - L_i)` and `omega(x) = sum_{i in E} e_i prod_{j in E, j != i}
   * (x - L_j)`, the syndrome `S(x) = sum_i r_i / (x - L_i) mod g(x)` satisfies
   * the key equation
   *
   *     sigma(x) * S(x) = omega(x)   (mod g(x)),
   *
   * with `deg(omega) < deg(sigma) <= floor(deg(g) / 2)`.  Evaluating it at
   * `L_i` gives Forney's formula `e_i = omega(L_i) / sigma'(L_i)`.
   *
   * The pair `(omega, sigma)` is recovered from the balanced partial extended
   * Euclidean algorithm on `(g, S)`; this is exactly the `_partial_xgcd` used
   * by SageMath's key-equation GRS decoder (`grs_code.py:2145-2157`), whose
   * loop runs `while r.degree() >= t.degree()`.
   */
  private _general_decode(received: P[]): P[] {
    const g = this._generating_pol;
    const L = this._defining_set;

    // Compute syndrome
    const S = this.syndrome(received);

    if (S.isZero()) {
      return [...received];
    }

    const [omega, sigma] = this._partialXGCDBalanced(g, S);

    if (sigma.isZero()) {
      throw new DecodingError('decoding failed: error locator polynomial is zero');
    }

    // Find roots of sigma among the support
    const errorPositions = this._findRoots(sigma, L);

    if (errorPositions.length !== sigma.degree()) {
      throw new DecodingError(
        `decoding failed: found ${errorPositions.length} roots but expected ${sigma.degree()}`
      );
    }

    const sigmaDerivative = sigma.derivative();
    const corrected = [...received];

    for (const pos of errorPositions) {
      const Li = L[pos]!;
      const denominator = sigmaDerivative.evaluate(Li);

      if (denominator.isZero()) {
        throw new DecodingError("decoding failed: sigma' vanishes at an error location");
      }

      // Forney: e_i = omega(L_i) / sigma'(L_i)
      const e = omega.evaluate(Li).mul(denominator.inv()) as E;
      corrected[pos] = corrected[pos]!.sub(e as unknown as P) as P;
    }

    if (!this.is_codeword(corrected)) {
      throw new DecodingError('decoding failed: corrected word is not a codeword');
    }

    return corrected;
  }

  /**
   * Balanced partial extended Euclidean algorithm.
   *
   * Port of `sage/coding/grs_code.py:2145-2157`
   * (`GRSKeyEquationSyndromeDecoder._partial_xgcd`): run the Euclidean
   * algorithm on `(a, b)` while `deg(r) >= deg(t)`, and return the pair
   * `(r, t)` reached, so that `r = a*s + b*t` for some `s`.
   */
  private _partialXGCDBalanced(a: Polynomial<E>, b: Polynomial<E>): [Polynomial<E>, Polynomial<E>] {
    const R = a.parent;

    let prevT = R.zero();
    let t = R.one();
    let prevR = a;
    let r = b;

    while (!r.isZero() && r.degree() >= t.degree()) {
      const [q] = prevR.quo_rem(r);
      const nextR = prevR.sub(q.mul(r));
      prevR = r;
      r = nextR;
      const nextT = prevT.sub(q.mul(t));
      prevT = t;
      t = nextT;
    }

    return [r, t];
  }

  /**
   * Return a string representation of this code.
   */
  toString(): string {
    const p = this._extension_field.characteristic;
    return `[${this.length()}, ${this.dimension()}] Goppa code over GF(${p})`;
  }
}

/**
 * Binary Goppa code - specialized for GF(2).
 *
 * Binary Goppa codes with a squarefree Goppa polynomial have minimum distance
 * >= 2t + 1, where t is the degree of that polynomial. This makes them
 * particularly useful for the McEliece cryptosystem.  (`distance_bound()`
 * still reports SageMath's bound, `1 + deg(g)`.)
 *
 * @example
 * ```typescript
 * const F = GFExtended(64);  // GF(2^6)
 * const R = new PolynomialRing(F, 'x');
 * const x = R.gen();
 * const g = x.pow(9).add(R.one());
 * const L = [...F].filter(a => !g.evaluate(a).isZero());
 * const C = new BinaryGoppaCode(g, L);
 * ```
 */
export class BinaryGoppaCode<
  E extends FieldElement,
  P extends FieldElement = FieldElement,
> extends GoppaCode<E, P> {
  constructor(generating_pol: Polynomial<E>, defining_set: E[]) {
    // Verify characteristic is 2
    const F = generating_pol.parent.base_ring as FiniteField<E>;
    if (F.characteristic !== 2n) {
      throw new ValueError('BinaryGoppaCode requires characteristic 2');
    }

    super(generating_pol, defining_set);
  }
}

/**
 * Create a Goppa code with automatic support generation.
 *
 * Generates the support as all field elements where g(a) != 0.
 *
 * @param g - Goppa polynomial over GF(q^m)
 * @param field - The finite field GF(q^m) (must be iterable)
 * @returns A Goppa code with maximal support
 */
export function createGoppaCode<E extends FieldElement>(
  g: Polynomial<E>,
  field: FiniteField<E>
): GoppaCode<E> {
  if (!field[Symbol.iterator]) {
    throw new ValueError('field must be iterable to generate support automatically');
  }

  const support: E[] = [];
  for (const a of field) {
    if (!g.evaluate(a).isZero()) {
      support.push(a);
    }
  }

  return new GoppaCode(g, support);
}

/**
 * Create a binary Goppa code with automatic support generation.
 *
 * @param g - Goppa polynomial over GF(2^m)
 * @param field - The finite field GF(2^m) (must be iterable)
 * @returns A binary Goppa code with maximal support
 */
export function createBinaryGoppaCode<E extends FieldElement>(
  g: Polynomial<E>,
  field: FiniteField<E>
): BinaryGoppaCode<E> {
  if (field.characteristic !== 2n) {
    throw new ValueError('binary Goppa code requires field of characteristic 2');
  }

  if (!field[Symbol.iterator]) {
    throw new ValueError('field must be iterable to generate support automatically');
  }

  const support: E[] = [];
  for (const a of field) {
    if (!g.evaluate(a).isZero()) {
      support.push(a);
    }
  }

  return new BinaryGoppaCode(g, support);
}
