/**
 * @module sage/coding/reed_muller_code
 * @description Reed-Muller codes over GF(2)
 *
 * Port of: sage/coding/reed_muller_code.py
 *
 * Reed-Muller codes are a family of linear error-correcting codes that have
 * nice recursive structure (Plotkin construction) and are used in various
 * applications including ZK systems.
 *
 * The binary Reed-Muller code RM(r, m) of order r and m variables is defined as:
 *   RM(r, m) = { (f(a) | a in GF(2)^m) | f in GF(2)[x_1,...,x_m], deg(f) <= r }
 *
 * Properties:
 * - Length: n = 2^m
 * - Dimension: k = sum_{i=0}^{r} C(m, i)
 * - Minimum distance: d = 2^{m-r}
 *
 * Special cases:
 * - RM(0, m) = repetition code of length 2^m
 * - RM(1, m) = first-order RM code (used in many applications)
 * - RM(m-1, m) = single parity check code
 * - RM(m, m) = full space GF(2)^{2^m}
 *
 * Recursive structure (Plotkin construction):
 *   RM(r, m) = {(u, u+v) : u in RM(r, m-1), v in RM(r-1, m-1)}
 */

import { binomial } from '../arith/misc.js';
import { ArithmeticError, ValueError } from '../errors.js';
import { Matrix, identity_matrix, zero_matrix } from '../matrix/matrix_generic.js';
import { GF2, GF2Element, type GF2Field } from '../rings/finite_rings/gf2.js';
import { type IntegerLike, toBigInt, toSafeNumber } from '../types/coercion.js';

/**
 * Sum of binomial coefficients: sum_{i=0}^{k} C(n, i)
 *
 * This is the dimension of RM(k, n).
 *
 * @param n - The number of variables (m)
 * @param k - The order (r), where k <= n
 * @returns Sum of binomial(n, i) for i from 0 to k
 */
function binomialSum(n: number, k: number): number {
  let sum = 1n; // binomial(n, 0) = 1
  let nCi = 1n;

  for (let i = 0; i < k; i++) {
    nCi = (nCi * BigInt(n - i)) / BigInt(i + 1);
    sum += nCi;
  }

  return Number(sum);
}

/**
 * Generate all binary vectors of length m (points in GF(2)^m).
 *
 * The ordering is lexicographic on the binary representation:
 * For m=3: (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1), ...
 *
 * @param m - Number of variables
 * @returns Array of 2^m binary vectors
 */
function generatePoints(m: number): number[][] {
  const n = 2 ** m; // 2^m
  const points: number[][] = [];

  for (let i = 0; i < n; i++) {
    const point: number[] = [];
    let rest = i;
    for (let j = 0; j < m; j++) {
      point.push(rest % 2);
      rest = Math.floor(rest / 2);
    }
    points.push(point);
  }

  return points;
}

/**
 * Generate all monomials of degree at most r in m variables.
 *
 * A monomial is represented by its exponent vector (e_1, ..., e_m)
 * where each e_i is in {0, 1} (since we're in GF(2), x_i^2 = x_i).
 *
 * The degree of a monomial is sum(e_i).
 *
 * The ordering matches SageMath's `ReedMullerVectorEncoder.generator_matrix`
 * (`reed_muller_code.py:666-680`), which loops `for degree in range(order+1)`
 * and, inside each degree, enumerates
 * `Subsets(list(range(num_of_var)) * max_individual_degree, degree,
 * submultiset=True)`.  Over GF(2) `max_individual_degree = min(r, q-1) = 1`,
 * so this is the size-`degree` subsets of `{0, ..., m-1}` in increasing
 * lexicographic order of their (ascending) index tuples: for `m = 3` and
 * `degree = 1` that is `x0, x1, x2`, and for `degree = 2` it is
 * `x0*x1, x0*x2, x1*x2`.  Sage's `GF(3)`, `r = 2`, `m = 2` doctest confirms
 * the convention (rows `1, x0, x1, x0^2, x0*x1, x1^2`).
 *
 * @param m - Number of variables
 * @param r - Maximum degree
 * @returns Array of exponent vectors for monomials of degree <= r
 */
function generateMonomials(m: number, r: number): number[][] {
  const monomials: number[][] = [];

  // Generate all 2^m possible exponent vectors
  const n = 2 ** m;
  for (let i = 0; i < n; i++) {
    const exponent: number[] = [];
    let degree = 0;
    let rest = i;

    for (let j = 0; j < m; j++) {
      const bit = rest % 2;
      rest = Math.floor(rest / 2);
      exponent.push(bit);
      degree += bit;
    }

    if (degree <= r) {
      monomials.push(exponent);
    }
  }

  // Sort by degree, then by increasing lexicographic order of the ascending
  // index tuple of the monomial's support -- i.e. at the first index where the
  // exponents differ, the monomial that *uses* that variable comes first.
  monomials.sort((a, b) => {
    const degA = a.reduce((s, x) => s + x, 0);
    const degB = b.reduce((s, x) => s + x, 0);
    if (degA !== degB) return degA - degB;
    for (let i = 0; i < a.length; i++) {
      if (a[i]! !== b[i]!) return b[i]! - a[i]!;
    }
    return 0;
  });

  return monomials;
}

/**
 * Evaluate a monomial at a point.
 *
 * The monomial is given by its exponent vector, the point by its coordinate vector.
 * In GF(2), the monomial x_1^{e_1} * ... * x_m^{e_m} evaluated at (a_1, ..., a_m)
 * equals prod_i a_i^{e_i} = prod_{e_i=1} a_i.
 *
 * @param exponent - Exponent vector of the monomial
 * @param point - Point in GF(2)^m
 * @returns 0 or 1
 */
function evaluateMonomial(exponent: number[], point: number[]): number {
  for (let i = 0; i < exponent.length; i++) {
    if (exponent[i] === 1 && point[i] === 0) {
      return 0;
    }
  }
  return 1;
}

/**
 * Binary Reed-Muller Code RM(r, m).
 *
 * A Reed-Muller code of order r and m variables over GF(2).
 *
 * The code RM(r, m) has:
 * - Length: 2^m
 * - Dimension: sum_{i=0}^{r} C(m, i)
 * - Minimum distance: 2^{m-r}
 *
 * @example
 * ```typescript
 * // Create RM(1, 3) - the first-order Reed-Muller code with 3 variables
 * const rm = new ReedMullerCode(1, 3);
 * console.log(rm.length());     // 8
 * console.log(rm.dimension());  // 4
 * console.log(rm.minimum_distance()); // 4
 *
 * // Encode a message
 * const message = [1, 0, 1, 1].map(x => GF2.__call__(x));
 * const codeword = rm.encode(message);
 *
 * // Add some errors and decode
 * codeword[0] = GF2.__call__(1 - codeword[0].value);  // flip bit
 * const decoded = rm.decode(codeword);
 * ```
 */
export class ReedMullerCode {
  private readonly _order: number;
  private readonly _num_variables: number;
  private readonly _length: number;
  private readonly _dimension: number;
  private readonly _min_distance: number;
  private _generator_matrix: Matrix<GF2Element> | null = null;
  private _monomials: number[][] | null = null;
  private _points: number[][] | null = null;

  /**
   * Create a binary Reed-Muller code RM(r, m).
   *
   * @param r - Order of the code (0 <= r <= m)
   * @param m - Number of variables
   * @throws {ValueError} If r < 0 or r > m
   * @throws {ValueError} If m < 0
   */
  constructor(r: IntegerLike, m: IntegerLike) {
    const rNum = toSafeNumber(toBigInt(r));
    const mNum = toSafeNumber(toBigInt(m));

    if (!Number.isInteger(rNum)) {
      throw new ValueError('The order of the code must be an integer');
    }
    if (!Number.isInteger(mNum)) {
      throw new ValueError('The number of variables must be an integer');
    }
    if (mNum < 0) {
      throw new ValueError('The number of variables must be non-negative');
    }
    if (rNum < 0) {
      throw new ValueError('The order must be non-negative');
    }
    if (rNum > mNum) {
      throw new ValueError(`The order must be at most ${mNum}`);
    }

    this._order = rNum;
    this._num_variables = mNum;
    // `1 << m` is a 32-bit shift in JavaScript: it returns 1 for m = 32 and a
    // negative number for m = 31.  Sage's length is q^m computed over ZZ.
    this._length = 2 ** mNum; // 2^m
    this._dimension = binomialSum(mNum, rNum);
    this._min_distance = 2 ** (mNum - rNum); // 2^(m-r)
  }

  /**
   * Return the order (r) of this Reed-Muller code.
   */
  order(): number {
    return this._order;
  }

  /**
   * Return the number of variables (m) of this Reed-Muller code.
   */
  num_variables(): number {
    return this._num_variables;
  }

  /**
   * Alias for num_variables() for API consistency.
   */
  number_of_variables(): number {
    return this._num_variables;
  }

  /**
   * Return the length (n = 2^m) of this Reed-Muller code.
   */
  length(): number {
    return this._length;
  }

  /**
   * Return the dimension (k) of this Reed-Muller code.
   *
   * The dimension is sum_{i=0}^{r} C(m, i).
   */
  dimension(): number {
    return this._dimension;
  }

  /**
   * Return the minimum distance (d = 2^{m-r}) of this Reed-Muller code.
   */
  minimum_distance(): number {
    return this._min_distance;
  }

  /**
   * Return the base field (GF(2)).
   */
  base_field(): GF2Field {
    return GF2;
  }

  /**
   * Return the code parameters as [n, k, d].
   */
  parameters(): [number, number, number] {
    return [this._length, this._dimension, this._min_distance];
  }

  /**
   * Return the decoding radius (floor((d-1)/2)).
   */
  decoding_radius(): number {
    return Math.floor((this._min_distance - 1) / 2);
  }

  /**
   * Return the code rate k/n.
   */
  rate(): number {
    return this._dimension / this._length;
  }

  /**
   * Get the monomials (cached).
   */
  private getMonomials(): number[][] {
    if (!this._monomials) {
      this._monomials = generateMonomials(this._num_variables, this._order);
    }
    return this._monomials;
  }

  /**
   * Get the evaluation points (cached).
   */
  private getPoints(): number[][] {
    if (!this._points) {
      this._points = generatePoints(this._num_variables);
    }
    return this._points;
  }

  /**
   * Return the evaluation points as arrays of GF2 elements.
   */
  evaluation_points(): GF2Element[][] {
    return this.getPoints().map((point) => point.map((x) => GF2.__call__(x)));
  }

  /**
   * Return the generator matrix of this Reed-Muller code.
   *
   * The generator matrix has dimension (k x n) where:
   * - Rows correspond to monomials of degree <= r
   * - Columns correspond to points in GF(2)^m
   * - Entry (i, j) = value of monomial i at point j
   */
  generator_matrix(): Matrix<GF2Element> {
    if (this._generator_matrix) {
      return this._generator_matrix;
    }

    const monomials = this.getMonomials();
    const points = this.getPoints();
    const k = this._dimension;
    const n = this._length;

    // Create the generator matrix
    const entries: GF2Element[][] = [];
    for (let i = 0; i < k; i++) {
      const row: GF2Element[] = [];
      for (let j = 0; j < n; j++) {
        const val = evaluateMonomial(monomials[i]!, points[j]!);
        row.push(GF2.__call__(val));
      }
      entries.push(row);
    }

    this._generator_matrix = new Matrix<GF2Element>(GF2, k, n, entries);
    return this._generator_matrix;
  }

  /**
   * Encode a message to a codeword.
   *
   * The message is interpreted as coefficients of a polynomial in the
   * monomial basis. The codeword is the evaluation of this polynomial
   * at all points in GF(2)^m.
   *
   * @param message - Array of k GF2 elements (or numbers 0/1)
   * @returns Array of n GF2 elements (the codeword)
   * @throws {ValueError} If message length is not k
   */
  encode(message: (GF2Element | number)[]): GF2Element[] {
    const k = this._dimension;
    const n = this._length;

    if (message.length !== k) {
      throw new ValueError(`message length must be ${k}, got ${message.length}`);
    }

    // Convert message to GF2 elements
    const msg: GF2Element[] = message.map((x) => (x instanceof GF2Element ? x : GF2.__call__(x)));

    // Multiply message by generator matrix
    const G = this.generator_matrix();
    const codeword: GF2Element[] = [];

    for (let j = 0; j < n; j++) {
      let sum = GF2.zero();
      for (let i = 0; i < k; i++) {
        sum = sum.add(msg[i]!.mul(G.get(i, j)));
      }
      codeword.push(sum);
    }

    return codeword;
  }

  /**
   * Encode a polynomial (given by its coefficient vector) to a codeword.
   *
   * The polynomial coefficients are given in order of the monomials:
   * constant, x_1, x_2, ..., x_1*x_2, ...
   *
   * This is the same as encode() but makes the polynomial interpretation explicit.
   *
   * @param coefficients - Polynomial coefficients
   * @returns The codeword
   */
  encode_from_polynomial(coefficients: (GF2Element | number)[]): GF2Element[] {
    return this.encode(coefficients);
  }

  /**
   * Decode a received word using Reed's majority logic decoding algorithm.
   *
   * For first-order Reed-Muller codes RM(1, m), this uses the characteristic
   * function approach with majority voting on pairs.
   *
   * For higher-order codes, this uses recursive Plotkin decoding.
   *
   * The algorithm can correct up to floor((d-1)/2) errors where d = 2^{m-r}.
   *
   * @param received - Array of n GF2 elements (the received word)
   * @returns Array of k GF2 elements (the decoded message)
   * @throws {ValueError} If received length is not n
   */
  decode(received: (GF2Element | number)[]): GF2Element[] {
    const n = this._length;
    const m = this._num_variables;
    const r = this._order;

    if (received.length !== n) {
      throw new ValueError(`received length must be ${n}, got ${received.length}`);
    }

    // Convert to array of 0/1 for easier manipulation
    const word: number[] = received.map((x) =>
      x instanceof GF2Element ? x.value : ((x % 2) + 2) % 2
    );

    // Use recursive Plotkin decoding to get the corrected codeword
    const decodedWord = this.recursiveDecode(word, r, m);

    // Extract the message from the decoded codeword
    return this.extractMessage(decodedWord);
  }

  /**
   * Recursively decode using Plotkin structure.
   *
   * For RM(r, m), a codeword has the form (u, u+v) where u in RM(r, m-1)
   * and v in RM(r-1, m-1).
   *
   * The key insight for error correction:
   * - v can be estimated from a XOR b (where a = first half, b = second half);
   *   the errors there are e_a + e_b <= t, and RM(r-1, m-1) has minimum
   *   distance 2^{m-r}, so it corrects exactly t errors and v is recovered
   * - u must then be estimated from BOTH halves: `a` carries e_a errors and
   *   `b + v` carries e_b errors, and since e_a + e_b <= t = 2^{m-r-1} - 1 at
   *   least one of them is within RM(r, m-1)'s radius 2^{m-r-2} - 1.  Both
   *   candidates are re-assembled into codewords of RM(r, m) and the one
   *   closest to the received word wins (the true codeword is the unique one
   *   within t of the received word).
   *
   * @param word - Received word as array of 0/1
   * @param r - Order parameter
   * @param m - Number of variables
   * @returns Decoded codeword as array of 0/1
   */
  private recursiveDecode(word: number[], r: number, m: number): number[] {
    const n = 2 ** m;

    // Base case: RM(0, m) - repetition code
    if (r === 0) {
      // Majority vote for the single bit
      let ones = 0;
      for (const bit of word) {
        ones += bit;
      }
      // Use > for strict majority
      const decoded = ones > n / 2 ? 1 : 0;
      return new Array(n).fill(decoded);
    }

    // Base case: RM(m, m) - full space, no correction possible
    if (r >= m) {
      return [...word];
    }

    // Base case: m = 0
    if (m === 0) {
      return [...word];
    }

    // For first-order codes, use efficient pair-based decoding
    if (r === 1) {
      const decoded = this.decodeFirstOrderRecursive(word, m);
      return decoded;
    }

    const half = n / 2;

    // Split received word: (a, b) where we received codeword (u, u+v) with errors
    const a = word.slice(0, half);
    const b = word.slice(half);

    // Step 1: Estimate v from a XOR b
    const vEstimate: number[] = [];
    for (let i = 0; i < half; i++) {
      vEstimate.push((a[i]! + b[i]!) % 2);
    }

    // Step 2: Decode v using RM(r-1, m-1)
    const vDecoded = this.recursiveDecode(vEstimate, r - 1, m - 1);

    // Step 3: two independent estimates of u -- the first half as received,
    // and the second half corrected by the decoded v.
    const uFromA = a;
    const uFromB: number[] = [];
    for (let i = 0; i < half; i++) {
      uFromB.push((b[i]! + vDecoded[i]!) % 2);
    }

    // Step 4: Decode each estimate using RM(r, m-1) and re-assemble (u, u+v)
    const assemble = (u: number[]): number[] => {
      const result: number[] = [...u];
      for (let i = 0; i < half; i++) {
        result.push((u[i]! + vDecoded[i]!) % 2);
      }
      return result;
    };

    const candidateA = assemble(this.recursiveDecode(uFromA, r, m - 1));
    const candidateB = assemble(this.recursiveDecode(uFromB, r, m - 1));

    // Step 5: keep the candidate closest to the received word
    let distA = 0;
    let distB = 0;
    for (let i = 0; i < n; i++) {
      if (candidateA[i] !== word[i]) distA++;
      if (candidateB[i] !== word[i]) distB++;
    }

    return distB < distA ? candidateB : candidateA;
  }

  /**
   * Decode first-order RM code and return codeword (for recursive use).
   */
  private decodeFirstOrderRecursive(word: number[], m: number): number[] {
    const n = 2 ** m;
    const coeffs: number[] = new Array(m + 1).fill(0);

    // Decode each linear coefficient using pairs
    for (let i = 0; i < m; i++) {
      let votes0 = 0;
      let votes1 = 0;

      for (let j = 0; j < n; j++) {
        if ((j & (1 << i)) === 0) {
          const partner = j | (1 << i);
          const estimate = (word[j]! + word[partner]!) % 2;
          if (estimate === 0) {
            votes0++;
          } else {
            votes1++;
          }
        }
      }
      coeffs[i + 1] = votes1 > votes0 ? 1 : 0;
    }

    // Decode constant term
    let votes0 = 0;
    let votes1 = 0;
    for (let j = 0; j < n; j++) {
      let linearPart = 0;
      for (let i = 0; i < m; i++) {
        if (coeffs[i + 1] === 1 && ((j >> i) & 1) === 1) {
          linearPart ^= 1;
        }
      }
      const estimate = (word[j]! + linearPart) % 2;
      if (estimate === 0) {
        votes0++;
      } else {
        votes1++;
      }
    }
    coeffs[0] = votes1 > votes0 ? 1 : 0;

    // Reconstruct codeword from coefficients
    const codeword: number[] = [];
    for (let j = 0; j < n; j++) {
      let val = coeffs[0]!;
      for (let i = 0; i < m; i++) {
        if (coeffs[i + 1] === 1 && ((j >> i) & 1) === 1) {
          val ^= 1;
        }
      }
      codeword.push(val);
    }

    return codeword;
  }

  /**
   * Extract the message from a codeword using the generator matrix.
   *
   * This finds coefficients c such that G * c = codeword.
   *
   * @param codeword - The decoded codeword as array of 0/1
   * @returns The message as array of GF2 elements
   */
  private extractMessage(codeword: number[]): GF2Element[] {
    const monomials = this.getMonomials();
    const points = this.getPoints();
    const k = this._dimension;

    // For Reed-Muller codes, we can extract coefficients directly
    // by evaluating on characteristic subsets
    const coeffs: number[] = new Array(k).fill(0);

    // Process monomials from lowest to highest degree
    const word = [...codeword];

    for (let idx = 0; idx < k; idx++) {
      const monomial = monomials[idx]!;

      // Find the unique point where this monomial equals 1 and
      // all lower-indexed monomials' variables are 0
      // For the constant (first monomial), this is the all-zeros point
      // which is point 0

      // Actually, for proper extraction, we use the structure:
      // The coefficient of monomial M is obtained by summing over
      // points where M = 1 and all variables not in M are 0

      // Find the "characteristic point" for this monomial
      let charIdx = 0;
      for (let i = 0; i < monomial.length; i++) {
        charIdx += monomial[i]! << i;
      }

      coeffs[idx] = word[charIdx]!;

      // Subtract contribution from this coefficient
      if (coeffs[idx] === 1) {
        for (let j = 0; j < this._length; j++) {
          const val = evaluateMonomial(monomial, points[j]!);
          word[j] = (word[j]! + val) % 2;
        }
      }
    }

    return coeffs.map((c) => GF2.__call__(c));
  }

  /**
   * Perform Plotkin decomposition of a codeword.
   *
   * RM(r, m) = {(u, u+v) : u in RM(r, m-1), v in RM(r-1, m-1)}
   *
   * Given a codeword c = (c_0, c_1, ..., c_{n-1}), decompose it into:
   * - u = first half = c[0:n/2]
   * - v = u + second half = c[0:n/2] + c[n/2:n]
   *
   * @param codeword - A codeword of RM(r, m)
   * @returns Tuple [u, v] where u is from RM(r, m-1) and v is from RM(r-1, m-1)
   * @throws {ValueError} If codeword length is wrong
   * @throws {ValueError} If m <= 0 (can't decompose)
   */
  plotkin_decomposition(codeword: (GF2Element | number)[]): [GF2Element[], GF2Element[]] {
    const n = this._length;
    const m = this._num_variables;

    if (m <= 0) {
      throw new ValueError('Cannot decompose RM(r, 0)');
    }

    if (codeword.length !== n) {
      throw new ValueError(`codeword length must be ${n}, got ${codeword.length}`);
    }

    const word = codeword.map((x) => (x instanceof GF2Element ? x : GF2.__call__(x)));

    const half = n / 2;

    // u = first half
    const u = word.slice(0, half);

    // v = first half XOR second half
    const v: GF2Element[] = [];
    for (let i = 0; i < half; i++) {
      v.push(u[i]!.add(word[i + half]!));
    }

    return [u, v];
  }

  /**
   * Construct a codeword from its Plotkin decomposition.
   *
   * Given u from RM(r, m-1) and v from RM(r-1, m-1),
   * return the codeword (u, u+v) in RM(r, m).
   *
   * @param u - First component from RM(r, m-1)
   * @param v - Second component from RM(r-1, m-1)
   * @returns The codeword (u, u+v)
   */
  static plotkin_compose(u: (GF2Element | number)[], v: (GF2Element | number)[]): GF2Element[] {
    if (u.length !== v.length) {
      throw new ValueError('u and v must have the same length');
    }

    const uVec = u.map((x) => (x instanceof GF2Element ? x : GF2.__call__(x)));
    const vVec = v.map((x) => (x instanceof GF2Element ? x : GF2.__call__(x)));

    // c = (u, u+v)
    const result: GF2Element[] = [...uVec];
    for (let i = 0; i < u.length; i++) {
      result.push(uVec[i]!.add(vVec[i]!));
    }

    return result;
  }

  /**
   * Return a punctured version of this Reed-Muller code.
   *
   * Puncturing removes the specified coordinate positions from all codewords.
   *
   * @param positions - Array of position indices to remove (0-indexed)
   * @returns A new PuncturedReedMullerCode
   */
  puncture(positions: number[]): PuncturedReedMullerCode {
    return new PuncturedReedMullerCode(this, positions);
  }

  /**
   * Check if this is the repetition code RM(0, m).
   */
  is_repetition_code(): boolean {
    return this._order === 0;
  }

  /**
   * Check if this is the first-order Reed-Muller code RM(1, m).
   */
  is_first_order(): boolean {
    return this._order === 1;
  }

  /**
   * Check if this is the single parity check code RM(m-1, m).
   */
  is_single_parity_check(): boolean {
    return this._order === this._num_variables - 1;
  }

  /**
   * Check if this is the full space RM(m, m) = GF(2)^{2^m}.
   */
  is_full_space(): boolean {
    return this._order === this._num_variables;
  }

  /**
   * Return the dual code.
   *
   * The dual of RM(r, m) is RM(m-r-1, m).
   *
   * @throws {ValueError} If r = m (dual would have negative order)
   */
  dual(): ReedMullerCode {
    if (this._order === this._num_variables) {
      throw new ValueError('The dual of RM(m, m) would have negative order (not defined)');
    }
    return new ReedMullerCode(
      BigInt(this._num_variables - this._order - 1),
      BigInt(this._num_variables)
    );
  }

  /**
   * Check if a word is a valid codeword.
   *
   * @param word - Array of n GF2 elements
   * @returns True if word is a codeword
   */
  contains(word: (GF2Element | number)[]): boolean {
    if (word.length !== this._length) {
      return false;
    }

    try {
      const decoded = this.decode(word);
      const reencoded = this.encode(decoded);

      const wordVec = word.map((x) => (x instanceof GF2Element ? x : GF2.__call__(x)));

      for (let i = 0; i < this._length; i++) {
        if (wordVec[i]!.value !== reencoded[i]!.value) {
          return false;
        }
      }
      return true;
    } catch {
      return false;
    }
  }

  /**
   * Compute the Hamming weight of a vector.
   */
  static hamming_weight(v: (GF2Element | number)[]): number {
    let weight = 0;
    for (const x of v) {
      const val = x instanceof GF2Element ? x.value : ((x % 2) + 2) % 2;
      weight += val;
    }
    return weight;
  }

  /**
   * Compute the Hamming distance between two vectors.
   */
  static hamming_distance(a: (GF2Element | number)[], b: (GF2Element | number)[]): number {
    if (a.length !== b.length) {
      throw new ValueError('vectors must have the same length');
    }

    let distance = 0;
    for (let i = 0; i < a.length; i++) {
      const aVal =
        a[i] instanceof GF2Element ? (a[i] as GF2Element).value : (((a[i] as number) % 2) + 2) % 2;
      const bVal =
        b[i] instanceof GF2Element ? (b[i] as GF2Element).value : (((b[i] as number) % 2) + 2) % 2;
      if (aVal !== bVal) {
        distance++;
      }
    }
    return distance;
  }

  /**
   * Return a string representation of this code.
   */
  toString(): string {
    return (
      `Binary Reed-Muller Code RM(${this._order}, ${this._num_variables}) ` +
      `with parameters [${this._length}, ${this._dimension}, ${this._min_distance}]`
    );
  }

  /**
   * Check equality with another Reed-Muller code.
   */
  eq(other: ReedMullerCode): boolean {
    return this._order === other._order && this._num_variables === other._num_variables;
  }
}

/**
 * Punctured Reed-Muller code.
 *
 * A punctured code is obtained by removing certain coordinate positions
 * from all codewords.
 */
export class PuncturedReedMullerCode {
  private readonly _base_code: ReedMullerCode;
  private readonly _punctured_positions: Set<number>;
  private readonly _remaining_positions: number[];
  private readonly _length: number;

  /**
   * Create a punctured Reed-Muller code.
   *
   * @param base_code - The original Reed-Muller code
   * @param positions - Positions to puncture (remove)
   */
  constructor(base_code: ReedMullerCode, positions: number[]) {
    this._base_code = base_code;
    this._punctured_positions = new Set(positions);

    const n = base_code.length();
    for (const pos of positions) {
      if (pos < 0 || pos >= n) {
        throw new ValueError(`Position ${pos} out of bounds [0, ${n})`);
      }
    }

    this._remaining_positions = [];
    for (let i = 0; i < n; i++) {
      if (!this._punctured_positions.has(i)) {
        this._remaining_positions.push(i);
      }
    }

    this._length = this._remaining_positions.length;
  }

  /**
   * Return the length of the punctured code.
   */
  length(): number {
    return this._length;
  }

  /**
   * Return the dimension (same as base code).
   */
  dimension(): number {
    return this._base_code.dimension();
  }

  /**
   * Encode a message using the punctured code.
   *
   * @param message - The message to encode
   * @returns Punctured codeword
   */
  encode(message: (GF2Element | number)[]): GF2Element[] {
    const full_codeword = this._base_code.encode(message);
    return this._remaining_positions.map((i) => full_codeword[i]!);
  }

  /**
   * Return the base (unpunctured) Reed-Muller code.
   */
  base_code(): ReedMullerCode {
    return this._base_code;
  }

  /**
   * Return the punctured positions.
   */
  punctured_positions(): number[] {
    return [...this._punctured_positions];
  }

  toString(): string {
    return (
      `Punctured Reed-Muller Code from RM(${this._base_code.order()}, ` +
      `${this._base_code.num_variables()}) with ${this._punctured_positions.size} positions removed`
    );
  }
}

// Factory functions for convenience

/**
 * Create a binary Reed-Muller code RM(r, m).
 *
 * @param r - Order of the code (0 <= r <= m)
 * @param m - Number of variables
 * @returns A ReedMullerCode instance
 */
export function BinaryReedMullerCode(r: IntegerLike, m: IntegerLike): ReedMullerCode {
  return new ReedMullerCode(r, m);
}

/**
 * Create a first-order Reed-Muller code RM(1, m).
 *
 * First-order Reed-Muller codes are particularly important in applications
 * as they have good error-correction properties and efficient decoding.
 *
 * @param m - Number of variables
 * @returns RM(1, m)
 */
export function FirstOrderReedMullerCode(m: IntegerLike): ReedMullerCode {
  return new ReedMullerCode(1n, m);
}

/**
 * Create a repetition code (as RM(0, m)).
 *
 * The repetition code of length 2^m has all codewords being either
 * all zeros or all ones.
 *
 * @param m - Determines the length 2^m
 * @returns RM(0, m)
 */
export function RepetitionCode(m: IntegerLike): ReedMullerCode {
  return new ReedMullerCode(0n, m);
}
