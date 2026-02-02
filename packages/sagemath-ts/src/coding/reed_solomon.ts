/**
 * @module sage/coding/reed_solomon
 * @description Reed-Solomon codes for ZK applications (FRI, STARKs)
 *
 * Port of: sage/coding/grs_code.py
 *
 * Reed-Solomon codes are maximum distance separable (MDS) codes that provide
 * optimal error correction for their rate. They are fundamental to:
 * - FRI (Fast Reed-Solomon Interactive Oracle Proof of Proximity)
 * - STARKs (Scalable Transparent Arguments of Knowledge)
 * - Traditional error-correcting codes
 *
 * A Reed-Solomon code RS(n, k) over a field F encodes a message of k symbols
 * into a codeword of n symbols by:
 * 1. Interpreting the k message symbols as coefficients of a polynomial of degree < k
 * 2. Evaluating this polynomial at n distinct points (evaluation points)
 *
 * The minimum distance is n - k + 1, allowing correction of up to floor((n-k)/2) errors.
 */

import { ValueError, ZeroDivisionError } from '../errors.js';
import type { CoefficientRing, RingElement } from '../rings/polynomial/polynomial_element.js';
import { Polynomial } from '../rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';
import { type IntegerLike, toBigInt, toSafeNumber } from '../types/coercion.js';

/**
 * A field element that supports all necessary operations for Reed-Solomon codes.
 */
export interface FieldElement extends RingElement {
  div(other: this): this;
  inv(): this;
  pow(n: number | bigint): this;
  isZero(): boolean;
  isOne?(): boolean;
}

/**
 * A finite field that supports all necessary operations for Reed-Solomon codes.
 */
export interface FiniteField<E extends FieldElement> extends CoefficientRing<E> {
  characteristic: bigint;
  order: bigint;
  multiplicative_generator?(): E;
  random_element?(): E;
  list?(): E[];
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
 * Reed-Solomon code over a finite field.
 *
 * Given n evaluation points alpha_1, ..., alpha_n in a field F, the RS code of
 * dimension k is the set of all codewords:
 *   { (f(alpha_1), ..., f(alpha_n)) | f in F[x], deg(f) < k }
 *
 * Properties:
 * - Length: n (number of evaluation points)
 * - Dimension: k (degree bound of message polynomial)
 * - Minimum distance: n - k + 1 (MDS property)
 * - Rate: k/n
 *
 * @example
 * ```typescript
 * // Create a Reed-Solomon code over GF(17) with length 8 and dimension 4
 * const F = GF(17n);
 * const rs = new ReedSolomonCode(F, 8, 4);
 *
 * // Encode a message
 * const message = [F(1n), F(2n), F(3n), F(4n)];
 * const codeword = rs.encode(message);
 *
 * // Introduce errors and decode
 * codeword[0] = F(99n);  // corrupt one position
 * const decoded = rs.decode(codeword);
 * ```
 */
export class ReedSolomonCode<E extends FieldElement> {
  readonly field: FiniteField<E>;
  readonly n: number;
  readonly k: number;
  readonly _evaluation_points: E[];
  readonly _polynomial_ring: PolynomialRing<E>;

  /**
   * Create a Reed-Solomon code.
   *
   * @param field - The finite field over which the code is defined
   * @param n - The length of the code (number of codeword symbols)
   * @param k - The dimension of the code (number of message symbols)
   * @param evaluation_points - Optional array of n distinct evaluation points.
   *   If not provided, uses consecutive powers of a primitive root: 1, g, g^2, ..., g^(n-1)
   * @throws {ValueError} If k > n or k < 1
   * @throws {ValueError} If evaluation points are not distinct
   * @throws {ValueError} If the field is too small for the code length
   */
  constructor(field: FiniteField<E>, n: IntegerLike, k: IntegerLike, evaluation_points?: E[]) {
    const nNum = toSafeNumber(toBigInt(n));
    const kNum = toSafeNumber(toBigInt(k));

    if (kNum < 1 || kNum > nNum) {
      throw new ValueError(`dimension k must satisfy 1 <= k <= n, got k=${kNum}, n=${nNum}`);
    }

    if (nNum > Number(field.order)) {
      throw new ValueError(`code length n=${nNum} exceeds field order ${field.order}`);
    }

    this.field = field;
    this.n = nNum;
    this.k = kNum;
    this._polynomial_ring = new PolynomialRing(field, 'x');

    if (evaluation_points) {
      if (evaluation_points.length !== nNum) {
        throw new ValueError(`expected ${nNum} evaluation points, got ${evaluation_points.length}`);
      }

      // Check that all points are distinct
      for (let i = 0; i < nNum; i++) {
        for (let j = i + 1; j < nNum; j++) {
          if (evaluation_points[i]!.eq(evaluation_points[j]!)) {
            throw new ValueError(
              `evaluation points must be distinct; points at indices ${i} and ${j} are equal`
            );
          }
        }
      }

      this._evaluation_points = [...evaluation_points];
    } else {
      // Generate default evaluation points as powers of a primitive root
      this._evaluation_points = this._generate_default_evaluation_points();
    }
  }

  /**
   * Generate default evaluation points as powers of a primitive root.
   * Uses 1, g, g^2, ..., g^(n-1) where g is a multiplicative generator.
   */
  private _generate_default_evaluation_points(): E[] {
    const points: E[] = [];

    if (this.field.multiplicative_generator) {
      const g = this.field.multiplicative_generator();
      let current = this.field.one() as E;

      for (let i = 0; i < this.n; i++) {
        points.push(current);
        current = current.mul(g) as E;
      }
    } else {
      // Fallback: use 0, 1, 2, ..., n-1 if possible
      for (let i = 0; i < this.n; i++) {
        points.push(this.field.__call__(i) as E);
      }
    }

    return points;
  }

  /**
   * Get the evaluation points of this code.
   */
  evaluation_points(): readonly E[] {
    return this._evaluation_points;
  }

  /**
   * Get the minimum distance of this code.
   *
   * For Reed-Solomon codes (which are MDS), d = n - k + 1.
   */
  minimum_distance(): number {
    return this.n - this.k + 1;
  }

  /**
   * Get the maximum number of errors that can be corrected.
   *
   * A Reed-Solomon code can correct up to floor((d-1)/2) = floor((n-k)/2) errors.
   */
  decoding_radius(): number {
    return Math.floor((this.n - this.k) / 2);
  }

  /**
   * Get the code rate k/n.
   */
  rate(): number {
    return this.k / this.n;
  }

  /**
   * Get the redundancy (number of parity symbols) n - k.
   */
  redundancy(): number {
    return this.n - this.k;
  }

  /**
   * Encode a message to a codeword.
   *
   * The message is interpreted as coefficients of a polynomial f(x) of degree < k:
   *   f(x) = message[0] + message[1]*x + ... + message[k-1]*x^(k-1)
   *
   * The codeword is the evaluation of f at all evaluation points:
   *   (f(alpha_1), f(alpha_2), ..., f(alpha_n))
   *
   * @param message - Array of k field elements (message polynomial coefficients)
   * @returns Array of n field elements (codeword)
   * @throws {ValueError} If message length is not k
   */
  encode(message: E[]): E[] {
    if (message.length !== this.k) {
      throw new ValueError(`message length must be ${this.k}, got ${message.length}`);
    }

    // Create the message polynomial
    const msgPoly = new Polynomial<E>(message, this._polynomial_ring);

    // Evaluate at each evaluation point
    const codeword: E[] = [];
    for (const alpha of this._evaluation_points) {
      codeword.push(msgPoly.evaluate(alpha));
    }

    return codeword;
  }

  /**
   * Encode a polynomial to a codeword by evaluating it at the evaluation points.
   *
   * @param poly - A polynomial of degree < k
   * @returns Array of n field elements (codeword)
   * @throws {ValueError} If polynomial degree >= k
   */
  encodePolynomial(poly: Polynomial<E>): E[] {
    if (poly.degree() >= this.k) {
      throw new ValueError(`polynomial degree must be < ${this.k}, got ${poly.degree()}`);
    }

    const codeword: E[] = [];
    for (const alpha of this._evaluation_points) {
      codeword.push(poly.evaluate(alpha));
    }

    return codeword;
  }

  /**
   * Decode a received word to recover the message.
   *
   * Uses the Berlekamp-Welch algorithm to correct up to floor((n-k)/2) errors.
   *
   * @param received - Array of n field elements (possibly corrupted codeword)
   * @returns Array of k field elements (decoded message)
   * @throws {DecodingError} If too many errors to correct
   * @throws {ValueError} If received length is not n
   */
  decode(received: E[]): E[] {
    if (received.length !== this.n) {
      throw new ValueError(`received length must be ${this.n}, got ${received.length}`);
    }

    const [_codeword, polynomial] = this._decode_to_code_and_message(received);

    // Extract coefficients from the polynomial
    const message: E[] = [];
    for (let i = 0; i < this.k; i++) {
      message.push(polynomial.getCoeff(i));
    }

    return message;
  }

  /**
   * Decode a received word and return both the corrected codeword and message polynomial.
   *
   * Uses a simple exhaustive search for small decoding radius, otherwise
   * attempts Gao's algorithm.
   *
   * @param received - Array of n field elements (possibly corrupted codeword)
   * @returns Tuple of [corrected_codeword, message_polynomial]
   * @throws {DecodingError} If too many errors to correct
   */
  private _decode_to_code_and_message(received: E[]): [E[], Polynomial<E>] {
    // Check if received is already a valid codeword by interpolating
    // using only the first k points and verifying
    const R = this._interpolateFirstK(received);

    // Check if this polynomial matches all received values
    let isValid = true;
    for (let i = 0; i < this.n; i++) {
      const expected = R.evaluate(this._evaluation_points[i]!);
      if (!expected.eq(received[i]!)) {
        isValid = false;
        break;
      }
    }

    if (isValid && R.degree() < this.k) {
      return [received, R];
    }

    const t = this.decoding_radius();
    if (t === 0) {
      throw new DecodingError('code has no error correction capability (decoding radius = 0)');
    }

    // Try Gao's algorithm
    const f = this._gaoDecoder(received);

    // Compute the corrected codeword
    const codeword = this.encodePolynomial(f);

    // Verify we didn't exceed the decoding radius
    let errorCount = 0;
    for (let i = 0; i < this.n; i++) {
      if (!codeword[i]!.eq(received[i]!)) {
        errorCount++;
      }
    }

    if (errorCount > t) {
      throw new DecodingError(`decoding failed: ${errorCount} errors exceeds decoding radius ${t}`);
    }

    return [codeword, f];
  }

  /**
   * Interpolate a polynomial of degree < k through the first k received values.
   */
  private _interpolateFirstK(values: E[]): Polynomial<E> {
    const points = this._evaluation_points.slice(0, this.k);
    const vals = values.slice(0, this.k);

    // Lagrange interpolation
    let result = this._polynomial_ring.zero();

    for (let i = 0; i < this.k; i++) {
      if (vals[i]!.isZero()) {
        continue;
      }

      // Build Lagrange basis polynomial L_i(x)
      let numerator = this._polynomial_ring.one();
      let denominator = this.field.one() as E;

      for (let j = 0; j < this.k; j++) {
        if (i === j) continue;

        const xMinusAlphaJ = new Polynomial<E>(
          [points[j]!.neg() as E, this.field.one() as E],
          this._polynomial_ring
        );
        numerator = numerator.mul(xMinusAlphaJ);
        denominator = denominator.mul(points[i]!.sub(points[j]!) as E) as E;
      }

      const scalar = vals[i]!.mul(denominator.inv()) as E;
      const term = numerator.scalar_mul(scalar);
      result = result.add(term);
    }

    return result;
  }

  /**
   * Gao's decoder using the extended Euclidean algorithm.
   * This works for RS codes where we have n evaluation points and want
   * to find a polynomial of degree < k that matches the received values
   * at most n - k positions with errors.
   */
  private _gaoDecoder(received: E[]): Polynomial<E> {
    // Step 1: Compute the vanishing polynomial G(x) = prod_i (x - alpha_i)
    const G = this._vanishingPolynomial();

    // Step 2: Compute interpolation polynomial R(x) through all received points
    const R = this._interpolate(received);

    // Step 3: Run partial extended Euclidean algorithm
    // We seek r(x) and t(x) such that r(x) = a(x) * G(x) + t(x) * R(x)
    // with deg(r) < (n + k) / 2 and deg(t) <= (n - k) / 2

    const stopDegree = Math.floor((this.n + this.k) / 2);

    // Extended Euclidean algorithm
    // We maintain: r_i = s_i * G + t_i * R
    let r0 = G;
    let r1 = R;
    let t0 = this._polynomial_ring.zero();
    let t1 = this._polynomial_ring.one();

    while (r1.degree() >= stopDegree) {
      if (r1.isZero()) {
        break;
      }

      const [q, remainder] = r0.quo_rem(r1);

      // Update r: r0 = r1, r1 = remainder
      r0 = r1;
      r1 = remainder;

      // Update t: t0' = t1, t1' = t0 - q * t1
      const newT = t0.sub(q.mul(t1));
      t0 = t1;
      t1 = newT;
    }

    // Now r1 = t1 * R (mod G), and we want f = r1 / t1
    if (t1.isZero()) {
      throw new DecodingError('decoding failed: error locator polynomial is zero');
    }

    const [f, rem] = r1.quo_rem(t1);
    if (!rem.isZero()) {
      throw new DecodingError('decoding failed: division has non-zero remainder (too many errors)');
    }

    if (f.degree() >= this.k) {
      throw new DecodingError('decoding failed: recovered polynomial has degree >= k');
    }

    return f;
  }

  /**
   * Compute the vanishing polynomial G(x) = prod_i (x - alpha_i).
   * This polynomial evaluates to zero at all evaluation points.
   */
  private _vanishingPolynomial(): Polynomial<E> {
    let G = this._polynomial_ring.one();
    const x = this._polynomial_ring.gen();

    for (const alpha of this._evaluation_points) {
      // Compute (x - alpha)
      const xMinusAlpha = x.sub(this._polynomial_ring.__call__(alpha));
      G = G.mul(xMinusAlpha);
    }

    return G;
  }

  /**
   * Interpolate a polynomial through the given evaluations.
   * Uses Lagrange interpolation.
   *
   * @param values - Array of n field elements (evaluations at the n evaluation points)
   * @returns The unique polynomial of degree < n interpolating the values
   */
  private _interpolate(values: E[]): Polynomial<E> {
    const points = this._evaluation_points;
    const n = points.length;

    // Lagrange interpolation
    let result = this._polynomial_ring.zero();

    for (let i = 0; i < n; i++) {
      if (values[i]!.isZero()) {
        continue;
      }

      // Build Lagrange basis polynomial L_i(x) = prod_{j != i} (x - alpha_j) / (alpha_i - alpha_j)
      let numerator = this._polynomial_ring.one();
      let denominator = this.field.one() as E;

      for (let j = 0; j < n; j++) {
        if (i === j) continue;

        // numerator *= (x - alpha_j)
        const xMinusAlphaJ = new Polynomial<E>(
          [points[j]!.neg() as E, this.field.one() as E],
          this._polynomial_ring
        );
        numerator = numerator.mul(xMinusAlphaJ);

        // denominator *= (alpha_i - alpha_j)
        denominator = denominator.mul(points[i]!.sub(points[j]!) as E) as E;
      }

      // L_i(x) = numerator / denominator
      const scalar = values[i]!.mul(denominator.inv()) as E;
      const term = numerator.scalar_mul(scalar);

      result = result.add(term);
    }

    return result;
  }

  /**
   * Compute the syndrome for a received word.
   *
   * For Reed-Solomon codes, the syndrome indicates whether a received word
   * is a valid codeword. A zero syndrome means the received word lies on
   * a polynomial of degree < k.
   *
   * The syndrome is computed by:
   * 1. Interpolating the received values to get polynomial R(x) of degree < n
   * 2. Checking the coefficients of degree k, k+1, ..., n-1
   *
   * These high-degree coefficients are zero if and only if the received word
   * is a valid codeword (i.e., corresponds to a polynomial of degree < k).
   *
   * @param received - The received word
   * @returns Array of n-k syndrome values (coefficients of x^k through x^(n-1))
   */
  syndrome(received: E[]): E[] {
    if (received.length !== this.n) {
      throw new ValueError(`received length must be ${this.n}, got ${received.length}`);
    }

    // Interpolate the received values
    const R = this._interpolate(received);

    // The syndrome is the high-degree coefficients (those that would be
    // non-zero if the polynomial has degree >= k)
    const syndromes: E[] = [];
    for (let j = this.k; j < this.n; j++) {
      syndromes.push(R.getCoeff(j));
    }

    return syndromes;
  }

  /**
   * Compute the error locator polynomial using the Berlekamp-Massey algorithm.
   *
   * Given the syndromes, finds the minimal LFSR that generates them.
   *
   * @param syndromes - Array of syndrome values
   * @returns The error locator polynomial Lambda(x)
   */
  error_locator(syndromes: E[]): Polynomial<E> {
    const t = Math.floor(syndromes.length / 2);

    // Berlekamp-Massey algorithm
    let Lambda = this._polynomial_ring.one();
    let B = this._polynomial_ring.one();
    let L = 0;
    let b = this.field.one() as E;

    for (let r = 0; r < syndromes.length; r++) {
      // Compute discrepancy
      let delta = syndromes[r]!;
      for (let j = 1; j <= L; j++) {
        const coeff = Lambda.getCoeff(j);
        if (!coeff.isZero()) {
          delta = delta.add(coeff.mul(syndromes[r - j]!) as E) as E;
        }
      }

      if (delta.isZero()) {
        // No update needed, shift B
        B = B.shift(1);
      } else {
        const T = Lambda;

        // Lambda = Lambda - (delta/b) * x * B
        const factor = delta.mul(b.inv()) as E;
        const shiftedB = B.shift(1);
        Lambda = Lambda.sub(shiftedB.scalar_mul(factor));

        if (2 * L <= r) {
          L = r + 1 - L;
          B = T;
          b = delta;
        } else {
          B = B.shift(1);
        }
      }
    }

    return Lambda;
  }

  /**
   * Compute the error evaluator polynomial.
   *
   * Omega(x) = S(x) * Lambda(x) mod x^(n-k)
   *
   * where S(x) = sum_j syndrome_j * x^j is the syndrome polynomial.
   *
   * @param syndromes - Array of syndrome values
   * @param locator - The error locator polynomial Lambda(x)
   * @returns The error evaluator polynomial Omega(x)
   */
  error_evaluator(syndromes: E[], locator: Polynomial<E>): Polynomial<E> {
    // Build syndrome polynomial S(x) = sum_j s_j * x^j
    const S = new Polynomial(syndromes, this._polynomial_ring);

    // Compute Omega = S * Lambda mod x^(n-k)
    const product = S.mul(locator);

    // Truncate to degree < n-k
    return product.truncate(this.n - this.k);
  }

  /**
   * Apply the Forney algorithm to find error values.
   *
   * Given error locator Lambda(x) and error evaluator Omega(x),
   * finds the error values at the error locations.
   *
   * For each root alpha_i^{-1} of Lambda(x) (where alpha_i is an error location):
   *   e_i = -Omega(alpha_i^{-1}) / Lambda'(alpha_i^{-1})
   *
   * @param locator - The error locator polynomial Lambda(x)
   * @param evaluator - The error evaluator polynomial Omega(x)
   * @returns Map from error position indices to error values
   */
  forney_algorithm(locator: Polynomial<E>, evaluator: Polynomial<E>): Map<number, E> {
    const errors = new Map<number, E>();
    const locatorDerivative = locator.derivative();

    for (let i = 0; i < this.n; i++) {
      const alpha = this._evaluation_points[i]!;
      const alphaInv = alpha.inv();

      // Check if alphaInv is a root of Lambda
      if (locator.evaluate(alphaInv).isZero()) {
        // Compute error value using Forney formula
        const omega = evaluator.evaluate(alphaInv);
        const lambdaPrime = locatorDerivative.evaluate(alphaInv);

        if (lambdaPrime.isZero()) {
          throw new DecodingError('Forney algorithm failed: derivative is zero at error location');
        }

        const errorValue = omega.neg().mul(lambdaPrime.inv()) as E;
        errors.set(i, errorValue);
      }
    }

    return errors;
  }

  // ============================================================
  // FRI-specific operations
  // ============================================================

  /**
   * Perform the FRI folding operation.
   *
   * In FRI (Fast Reed-Solomon Interactive Oracle Proof of Proximity),
   * folding combines pairs of evaluations using a random challenge
   * to reduce the code length by half while preserving proximity to
   * a lower-degree polynomial.
   *
   * For a codeword c of length n with evaluation points omega^0, ..., omega^(n-1):
   * The folded codeword c' of length n/2 has values:
   *   c'[i] = (c[i] + c[i + n/2]) / 2 + challenge * (c[i] - c[i + n/2]) / (2 * omega^i)
   *
   * This corresponds to:
   *   f'(x^2) = (f(x) + f(-x))/2 + challenge * (f(x) - f(-x))/(2x)
   *
   * where f is the original polynomial.
   *
   * @param codeword - Array of n field elements (must be even length)
   * @param challenge - Random field element (verifier challenge in FRI)
   * @returns Array of n/2 field elements (folded codeword)
   * @throws {ValueError} If codeword length is not even
   */
  fold(codeword: E[], challenge: E): E[] {
    if (codeword.length % 2 !== 0) {
      throw new ValueError('codeword length must be even for folding');
    }

    const halfN = codeword.length / 2;
    const folded: E[] = [];
    const two = this.field.__call__(2) as E;
    const twoInv = two.inv();

    for (let i = 0; i < halfN; i++) {
      const c0 = codeword[i]!;
      const c1 = codeword[i + halfN]!;
      const omega_i = this._evaluation_points[i]!;

      // even = (c0 + c1) / 2
      const even = c0.add(c1).mul(twoInv) as E;

      // odd = (c0 - c1) / (2 * omega_i)
      const oddNumerator = c0.sub(c1) as E;
      const oddDenominator = two.mul(omega_i) as E;
      const odd = oddNumerator.mul(oddDenominator.inv()) as E;

      // folded[i] = even + challenge * odd
      folded.push(even.add(challenge.mul(odd) as E) as E);
    }

    return folded;
  }

  /**
   * Query the codeword at specified indices.
   *
   * In FRI, the verifier queries the prover's claimed codeword at random
   * positions to verify proximity to a low-degree polynomial.
   *
   * @param codeword - The codeword to query
   * @param indices - Array of indices to query
   * @returns Array of [index, evaluation, evaluation_point] tuples
   */
  query(codeword: E[], indices: number[]): Array<[number, E, E]> {
    const results: Array<[number, E, E]> = [];

    for (const idx of indices) {
      if (idx < 0 || idx >= codeword.length) {
        throw new ValueError(`index ${idx} out of bounds [0, ${codeword.length})`);
      }

      results.push([
        idx,
        codeword[idx]!,
        this._evaluation_points[idx % this._evaluation_points.length]!,
      ]);
    }

    return results;
  }

  /**
   * Verify that a word is close to a valid codeword.
   *
   * Used in FRI to verify proximity - checks that the word differs
   * from some valid codeword in at most a bounded number of positions.
   *
   * @param word - The word to check
   * @param maxDistance - Maximum allowed Hamming distance from a codeword
   * @returns true if word is within maxDistance of a valid codeword
   */
  isClose(word: E[], maxDistance: number): boolean {
    if (word.length !== this.n) {
      return false;
    }

    try {
      const decoded = this.decode(word);
      const codeword = this.encode(decoded);

      let distance = 0;
      for (let i = 0; i < this.n; i++) {
        if (!word[i]!.eq(codeword[i]!)) {
          distance++;
          if (distance > maxDistance) {
            return false;
          }
        }
      }

      return true;
    } catch (e) {
      // If decoding fails, word is too far from any codeword
      return false;
    }
  }

  /**
   * Compute the Hamming distance between two words.
   *
   * @param a - First word
   * @param b - Second word
   * @returns Number of positions where a and b differ
   */
  hammingDistance(a: E[], b: E[]): number {
    if (a.length !== b.length) {
      throw new ValueError('words must have the same length');
    }

    let distance = 0;
    for (let i = 0; i < a.length; i++) {
      if (!a[i]!.eq(b[i]!)) {
        distance++;
      }
    }

    return distance;
  }

  /**
   * String representation.
   */
  toString(): string {
    return `[${this.n}, ${this.k}, ${this.minimum_distance()}] Reed-Solomon Code over ${this.field}`;
  }
}

/**
 * Construct a classical Reed-Solomon code.
 *
 * A classical Reed-Solomon code uses evaluation points that are consecutive
 * powers of a primitive n-th root of unity: 1, omega, omega^2, ..., omega^(n-1).
 *
 * @param field - The finite field
 * @param length - The code length n (must divide field_order - 1)
 * @param dimension - The message dimension k
 * @param primitiveRoot - Optional primitive n-th root of unity
 * @returns The Reed-Solomon code
 */
export function createClassicalReedSolomonCode<E extends FieldElement>(
  field: FiniteField<E>,
  length: IntegerLike,
  dimension: IntegerLike,
  primitiveRoot?: E
): ReedSolomonCode<E> {
  const lengthNum = toSafeNumber(toBigInt(length));
  const dimensionNum = toSafeNumber(toBigInt(dimension));
  const q = field.order;

  if ((q - 1n) % BigInt(lengthNum) !== 0n) {
    throw new ValueError(`length ${lengthNum} must divide field order - 1 = ${q - 1n}`);
  }

  let omega: E;

  if (primitiveRoot) {
    // Verify it's a primitive n-th root of unity
    if (!primitiveRoot.pow(lengthNum).eq(field.one())) {
      throw new ValueError('primitiveRoot is not an n-th root of unity');
    }
    for (let i = 1; i < lengthNum; i++) {
      if (primitiveRoot.pow(i).eq(field.one())) {
        throw new ValueError('primitiveRoot is not a primitive n-th root of unity');
      }
    }
    omega = primitiveRoot;
  } else {
    // Compute primitive n-th root of unity
    if (!field.multiplicative_generator) {
      throw new ValueError('field does not support multiplicative_generator');
    }
    const g = field.multiplicative_generator();
    const exponent = (q - 1n) / BigInt(lengthNum);
    omega = g.pow(exponent) as E;
  }

  // Generate evaluation points: 1, omega, omega^2, ..., omega^(n-1)
  const evalPoints: E[] = [];
  let current = field.one() as E;
  for (let i = 0; i < lengthNum; i++) {
    evalPoints.push(current);
    current = current.mul(omega) as E;
  }

  return new ReedSolomonCode(field, BigInt(lengthNum), BigInt(dimensionNum), evalPoints);
}

/**
 * Create a Reed-Solomon code suitable for FRI/STARKs.
 *
 * For FRI and STARKs, we typically use a code where:
 * - The evaluation domain is a multiplicative subgroup of the field
 * - The rate (k/n) is a power of 2 fraction (commonly 1/2, 1/4, 1/8)
 * - The length is a power of 2
 *
 * @param field - The finite field (order should be p where p-1 is divisible by n)
 * @param length - The code length (should be a power of 2)
 * @param dimension - The message dimension
 * @returns The Reed-Solomon code configured for FRI
 */
export function createFRIReedSolomonCode<E extends FieldElement>(
  field: FiniteField<E>,
  length: IntegerLike,
  dimension: IntegerLike
): ReedSolomonCode<E> {
  const lengthNum = toSafeNumber(toBigInt(length));
  const dimensionNum = toSafeNumber(toBigInt(dimension));

  // Verify length is a power of 2
  if ((lengthNum & (lengthNum - 1)) !== 0 || lengthNum === 0) {
    throw new ValueError(`length ${lengthNum} must be a power of 2`);
  }

  return createClassicalReedSolomonCode(field, BigInt(lengthNum), BigInt(dimensionNum));
}
