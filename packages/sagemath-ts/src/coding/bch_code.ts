/**
 * @module sage/coding/bch_code
 * @description BCH (Bose-Chaudhuri-Hocquenghem) error-correcting codes
 *
 * Port of: sage/coding/bch_code.py
 *
 * BCH codes are a class of cyclic error-correcting codes that are constructed using
 * polynomials over finite fields. They are named after Raj Bose, D.K. Ray-Chaudhuri,
 * and Alexis Hocquenghem.
 *
 * A BCH code over GF(q) with designed distance delta is a cyclic code whose
 * codewords c(x) in F[x] satisfy c(alpha^a) = 0, for all integers a in the
 * arithmetic sequence b, b + l, b + 2*l, ..., b + (delta - 2)*l.
 *
 * Key properties:
 * - The minimum distance d >= delta (the BCH bound)
 * - For narrow-sense BCH codes (b=1), tight error correction guarantees
 * - Generator polynomial is the LCM of minimal polynomials of roots
 */

import { NotImplementedError, ValueError } from '../errors.js';
import {
  type FiniteFieldElement,
  FiniteFieldExtension,
  PrimeField,
  type PrimeFieldElement,
} from '../rings/finite_rings/finite_field_extension.js';
import type { CoefficientRing, RingElement } from '../rings/polynomial/polynomial_element.js';
import { Polynomial } from '../rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';
import { type IntegerLike, toBigInt, toSafeNumber } from '../types/coercion.js';

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
 * A field element that supports all necessary operations.
 */
export interface FieldElement extends RingElement {
  div(other: this): this;
  inv(): this;
  pow(n: number | bigint): this;
  isZero(): boolean;
  isOne?(): boolean;
}

/**
 * A finite field that supports all necessary operations.
 */
export interface FiniteField<E extends FieldElement> extends CoefficientRing<E> {
  characteristic: bigint;
  order: bigint;
  degree?: number;
  cardinality?(): bigint;
  gen?(): E;
  random_element?(): E;
  [Symbol.iterator]?(): Iterator<E>;
}

/**
 * BCH code over a finite field.
 *
 * Given a finite field F = GF(q), a BCH code of length n and designed distance delta
 * is constructed by:
 * 1. Finding the splitting field GF(q^m) where n | (q^m - 1)
 * 2. Finding a primitive n-th root of unity alpha in GF(q^m)
 * 3. Computing the defining set D = {b, b+l, b+2l, ..., b+(delta-2)l}
 * 4. Computing the generator polynomial g(x) = LCM of minimal polynomials of alpha^i for i in D
 *
 * @example
 * ```typescript
 * // Create a BCH(15, 7) code over GF(2)
 * const F2 = new PrimeField(2);
 * const bch = new BCHCode(15, 7, F2);
 *
 * // Encode a message
 * const message = [F2(1), F2(0), F2(1), F2(1), F2(0)];
 * const codeword = bch.encode(message);
 *
 * // Introduce errors and decode
 * codeword[0] = codeword[0].add(F2(1));
 * const decoded = bch.decode(codeword);
 * ```
 */
export class BCHCode {
  readonly n: number; // code length
  readonly delta: number; // designed distance
  readonly baseField: PrimeField | FiniteFieldExtension; // base field GF(q)
  readonly splittingField: FiniteFieldExtension; // extension field GF(q^m)
  readonly offset: number; // b (first element in defining set), 1 for narrow-sense
  readonly jumpSize: number; // l (step size in defining set)

  private _primitiveRoot: FiniteFieldElement | null = null;
  private _definingSet: number[] | null = null;
  private _generatorPolynomial: Polynomial<PrimeFieldElement | FiniteFieldElement> | null = null;
  private _dimension: number | null = null;
  private _polynomialRing: PolynomialRing<PrimeFieldElement | FiniteFieldElement>;

  /**
   * Create a BCH code.
   *
   * @param n - Code length (must divide q^m - 1 for some m)
   * @param delta - Designed distance (must satisfy 1 <= delta <= n)
   * @param baseField - Base field GF(q) - can be a prime field or extension field
   * @param offset - First element in defining set (default: 1 for narrow-sense)
   * @param jumpSize - Jump size between elements of defining set (default: 1)
   *
   * @throws {ValueError} If delta is not in [1, n]
   * @throws {ValueError} If n does not divide q^m - 1 for any reasonable m
   */
  constructor(
    n: IntegerLike,
    delta: IntegerLike,
    baseField: PrimeField | FiniteFieldExtension,
    offset: IntegerLike = 1n,
    jumpSize: IntegerLike = 1n
  ) {
    const nNum = toSafeNumber(toBigInt(n));
    const deltaNum = toSafeNumber(toBigInt(delta));
    const offsetNum = toSafeNumber(toBigInt(offset));
    const jumpSizeNum = toSafeNumber(toBigInt(jumpSize));

    if (deltaNum < 1 || deltaNum > nNum) {
      throw new ValueError('designed_distance must belong to [1, n]');
    }

    this.n = nNum;
    this.delta = deltaNum;
    this.baseField = baseField;
    this.offset = offsetNum;
    this.jumpSize = jumpSizeNum;

    // Get the order of the base field
    const q = baseField.order;
    const p = baseField.characteristic;
    const baseDegree = baseField instanceof FiniteFieldExtension ? baseField.degree : 1;

    // Find the extension degree m such that n | (q^m - 1)
    const m = this.findExtensionDegree(nNum, q);

    if (m === -1) {
      throw new ValueError(`length ${nNum} does not divide q^m - 1 for any m <= 100 (q = ${q})`);
    }

    // Create the splitting field
    // If baseField is an extension field GF(p^k), the splitting field is GF(p^(k*m))
    // If baseField is a prime field GF(p), the splitting field is GF(p^m)
    const splittingDegree = baseDegree * m;
    this.splittingField = new FiniteFieldExtension(p, splittingDegree);

    // Create polynomial ring over base field
    const baseRing = baseField as CoefficientRing<PrimeFieldElement | FiniteFieldElement>;
    this._polynomialRing = new PolynomialRing(baseRing, 'x');
  }

  /**
   * Find the smallest m such that n | (q^m - 1).
   */
  private findExtensionDegree(n: number, q: bigint): number {
    for (let m = 1; m <= 100; m++) {
      const qm = q ** BigInt(m);
      if ((qm - 1n) % BigInt(n) === 0n) {
        return m;
      }
    }
    return -1;
  }

  /**
   * Get the effective base field order.
   */
  private getBaseFieldOrder(): bigint {
    return this.baseField.order;
  }

  /**
   * Get the degree of the splitting field relative to the base field.
   */
  splittingDegree(): number {
    const baseDegree = this.baseField instanceof FiniteFieldExtension ? this.baseField.degree : 1;
    return this.splittingField.degree / baseDegree;
  }

  /**
   * Return the code length.
   */
  length(): number {
    return this.n;
  }

  /**
   * Return the designed distance.
   */
  designed_distance(): number {
    return this.delta;
  }

  /**
   * Return the offset (b parameter).
   */
  get b(): number {
    return this.offset;
  }

  /**
   * Return the jump size (l parameter).
   */
  get l(): number {
    return this.jumpSize;
  }

  /**
   * Return the defining set D = {b, b+l, b+2l, ..., b+(delta-2)l} mod n.
   *
   * The defining set determines the roots of the generator polynomial.
   */
  defining_set(): number[] {
    if (this._definingSet !== null) {
      return this._definingSet;
    }

    const D: number[] = [];
    for (let i = 0; i < this.delta - 1; i++) {
      D.push((this.offset + this.jumpSize * i) % this.n);
    }

    this._definingSet = D;
    return D;
  }

  /**
   * Return a primitive n-th root of unity in the splitting field.
   */
  primitiveRoot(): FiniteFieldElement {
    if (this._primitiveRoot !== null) {
      return this._primitiveRoot;
    }

    // Find a primitive element of the multiplicative group
    const qm = this.splittingField.order;
    const groupOrder = qm - 1n;
    const exponent = groupOrder / BigInt(this.n);

    // Get a generator of the multiplicative group
    const g = this.splittingField.primitiveElement();

    // alpha = g^((q^m - 1) / n) is a primitive n-th root of unity
    this._primitiveRoot = g.pow(exponent);

    return this._primitiveRoot;
  }

  /**
   * Compute the minimal polynomial of an element over the base field.
   *
   * The minimal polynomial of alpha^i is the product (x - alpha^i)(x - alpha^{qi})(x - alpha^{q^2 i})...
   * taken over the cyclotomic coset of i modulo n, where q is the base field size.
   */
  private minimalPolynomial(i: number): Polynomial<PrimeFieldElement | FiniteFieldElement> {
    const alpha = this.primitiveRoot();
    // Use the base field ORDER (not characteristic) for cyclotomic cosets
    const q = Number(this.baseField.order);

    // Find the cyclotomic coset of i modulo n
    const coset: number[] = [];
    let j = i % this.n;
    const seen = new Set<number>();

    while (!seen.has(j)) {
      seen.add(j);
      coset.push(j);
      j = (j * q) % this.n;
    }

    // Build the minimal polynomial over the splitting field polynomial ring
    const extPolyRing = new PolynomialRing(this.splittingField, 'x');

    let minPoly = extPolyRing.one();
    const x = extPolyRing.gen();

    for (const k of coset) {
      // (x - alpha^k)
      const root = alpha.pow(k);
      const factor = x.sub(extPolyRing.__call__(root));
      minPoly = minPoly.mul(factor);
    }

    // The minimal polynomial has coefficients in the base field
    // Convert coefficients from extension field to base field
    if (this.baseField instanceof PrimeField) {
      const baseCoeffs: PrimeFieldElement[] = [];
      for (let deg = 0; deg <= minPoly.degree(); deg++) {
        const coeff = minPoly.getCoeff(deg);
        // Extract the constant term (coefficient of degree 0 in the polynomial representation)
        const baseCoeff = coeff.lift.getCoeff(0);
        baseCoeffs.push(baseCoeff);
      }
      return new Polynomial(baseCoeffs, this._polynomialRing as PolynomialRing<PrimeFieldElement>);
    } else {
      // Base field is an extension field
      const baseField = this.baseField as FiniteFieldExtension;
      const baseDegree = baseField.degree;
      const baseCoeffs: FiniteFieldElement[] = [];

      for (let deg = 0; deg <= minPoly.degree(); deg++) {
        const coeff = minPoly.getCoeff(deg);
        // Extract coefficients up to the base field degree
        const coeffValues: number[] = [];
        for (let k = 0; k < baseDegree; k++) {
          coeffValues.push(Number(coeff.lift.getCoeff(k).value));
        }
        baseCoeffs.push(baseField.__call__(coeffValues));
      }
      return new Polynomial(baseCoeffs, this._polynomialRing as PolynomialRing<FiniteFieldElement>);
    }
  }

  /**
   * Return the generator polynomial g(x).
   *
   * The generator polynomial is the LCM of the minimal polynomials of alpha^i
   * for all i in the defining set D.
   *
   * @returns The generator polynomial over the base field
   */
  generator_polynomial(): Polynomial<PrimeFieldElement | FiniteFieldElement> {
    if (this._generatorPolynomial !== null) {
      return this._generatorPolynomial;
    }

    const D = this.defining_set();

    if (D.length === 0) {
      // delta = 1: trivial code
      this._generatorPolynomial = this._polynomialRing.one();
      return this._generatorPolynomial;
    }

    // Find all cyclotomic cosets represented in D
    // Use base field ORDER (not characteristic)
    const q = Number(this.baseField.order);
    const representatives = new Set<number>();
    const processed = new Set<number>();

    for (const i of D) {
      if (processed.has(i)) continue;

      // Add this element as representative of its coset
      let j = i % this.n;
      if (j < 0) j += this.n;
      while (!processed.has(j)) {
        processed.add(j);
        j = (j * q) % this.n;
      }

      // Find the canonical representative (smallest in coset)
      let rep = i;
      j = (i * q) % this.n;
      while (j !== i) {
        if (j < rep) rep = j;
        j = (j * q) % this.n;
      }
      representatives.add(rep);
    }

    // Compute LCM of minimal polynomials
    let g = this._polynomialRing.one();

    for (const rep of representatives) {
      const mp = this.minimalPolynomial(rep);
      g = this.polynomialLCM(g, mp);
    }

    this._generatorPolynomial = g;
    return g;
  }

  /**
   * Compute LCM of two polynomials: LCM(a, b) = a * b / GCD(a, b)
   */
  private polynomialLCM<T extends PrimeFieldElement | FiniteFieldElement>(
    a: Polynomial<T>,
    b: Polynomial<T>
  ): Polynomial<T> {
    if (a.isZero()) return b;
    if (b.isZero()) return a;

    const gcd = a.gcd(b);
    const [quotient, _remainder] = a.mul(b).quo_rem(gcd);
    return quotient.monic();
  }

  /**
   * Return the dimension k = n - deg(g) of the code.
   */
  dimension(): number {
    if (this._dimension !== null) {
      return this._dimension;
    }

    const g = this.generator_polynomial();
    this._dimension = this.n - g.degree();
    return this._dimension;
  }

  /**
   * Return the minimum distance of the code.
   *
   * By the BCH bound, d >= delta. For many BCH codes, d = delta.
   */
  minimum_distance(): number {
    // The BCH bound guarantees d >= delta
    // Computing the exact minimum distance requires more work
    return this.delta;
  }

  /**
   * Return the decoding radius (number of errors that can be corrected).
   *
   * t = floor((delta - 1) / 2)
   */
  decoding_radius(): number {
    return Math.floor((this.delta - 1) / 2);
  }

  /**
   * Encode a message to a codeword using systematic encoding.
   *
   * For systematic encoding:
   * c(x) = m(x) * x^{n-k} - (m(x) * x^{n-k} mod g(x))
   *
   * This places the message symbols in the high-order positions.
   *
   * @param message - Array of k field elements (the message)
   * @returns Array of n field elements (the codeword)
   * @throws {ValueError} If message length is not k
   */
  encode(message: PrimeFieldElement[]): PrimeFieldElement[] {
    const k = this.dimension();

    if (message.length !== k) {
      throw new ValueError(`message length must be ${k}, got ${message.length}`);
    }

    const g = this.generator_polynomial();

    // Create message polynomial m(x)
    const m = new Polynomial(message, this._polynomialRing);

    // Compute m(x) * x^{n-k}
    const shifted = m.shift(this.n - k);

    // Compute remainder: r(x) = m(x) * x^{n-k} mod g(x)
    const [_quotient, remainder] = shifted.quo_rem(g);

    // Codeword: c(x) = m(x) * x^{n-k} - r(x)
    const codewordPoly = shifted.sub(remainder);

    // Extract coefficients
    const codeword: PrimeFieldElement[] = [];
    for (let i = 0; i < this.n; i++) {
      codeword.push(codewordPoly.getCoeff(i));
    }

    return codeword;
  }

  /**
   * Compute the syndromes of a received word.
   *
   * S_i = r(alpha^i) for i = b, b+l, b+2l, ..., b+(delta-2)*l
   *
   * The syndromes are in the splitting field GF(q^m).
   *
   * @param received - Array of n field elements
   * @returns Array of delta-1 syndrome values in the splitting field
   */
  syndrome(received: PrimeFieldElement[]): FiniteFieldElement[] {
    if (received.length !== this.n) {
      throw new ValueError(`received length must be ${this.n}, got ${received.length}`);
    }

    const alpha = this.primitiveRoot();
    const syndromes: FiniteFieldElement[] = [];

    // Create received polynomial in the splitting field
    const extPolyRing = new PolynomialRing(this.splittingField, 'x');
    const rCoeffs = received.map((c) => this.splittingField.__call__(c.value));
    const r = new Polynomial(rCoeffs, extPolyRing);

    // Compute syndromes: S_j = r(alpha^{b + l*j}) for j = 0, 1, ..., delta-2
    for (let j = 0; j < this.delta - 1; j++) {
      const exp = this.offset + this.jumpSize * j;
      const alphaExp = alpha.pow(exp);
      syndromes.push(r.evaluate(alphaExp));
    }

    return syndromes;
  }

  /**
   * Check if a word is a valid codeword (all syndromes are zero).
   */
  is_codeword(word: PrimeFieldElement[]): boolean {
    if (word.length !== this.n) {
      return false;
    }

    const syndromes = this.syndrome(word);
    return syndromes.every((s) => s.isZero());
  }

  /**
   * Decode a received word to recover the original codeword.
   *
   * Uses the Peterson-Gorenstein-Zierler (PGZ) algorithm:
   * 1. Compute syndromes
   * 2. Build syndrome matrix and solve for error locator polynomial
   * 3. Find roots using Chien search
   * 4. Correct errors at those positions
   *
   * @param received - Array of n field elements (possibly corrupted codeword)
   * @returns Array of n field elements (corrected codeword)
   * @throws {DecodingError} If too many errors to correct
   */
  decode(received: PrimeFieldElement[]): PrimeFieldElement[] {
    if (received.length !== this.n) {
      throw new ValueError(`received length must be ${this.n}, got ${received.length}`);
    }

    // Compute syndromes
    const syndromes = this.syndrome(received);

    // Check if already a valid codeword
    if (syndromes.every((s) => s.isZero())) {
      return [...received];
    }

    const t = this.decoding_radius();

    if (t === 0) {
      throw new DecodingError('code has no error correction capability');
    }

    // Find error locator polynomial using PGZ algorithm
    const Lambda = this.pgzErrorLocator(syndromes, t);

    if (Lambda === null) {
      throw new DecodingError('decoding failed: could not find error locator polynomial');
    }

    // Find error positions using Chien search
    const errorPositions = this.chienSearch(Lambda);

    if (errorPositions.length !== Lambda.degree()) {
      throw new DecodingError(
        `decoding failed: found ${errorPositions.length} roots but expected ${Lambda.degree()}`
      );
    }

    // For binary codes, error values are always 1
    // For non-binary codes, we would use Forney algorithm
    if (this.baseField.characteristic === 2n) {
      // Binary code: flip bits at error positions
      const corrected = [...received];
      for (const pos of errorPositions) {
        corrected[pos] = corrected[pos]!.add(this.baseField.one());
      }
      return corrected;
    }

    // Non-binary: use Forney's algorithm to find error values
    const errorValues = this.forneyAlgorithm(syndromes, Lambda, errorPositions);

    const corrected = [...received];
    for (let i = 0; i < errorPositions.length; i++) {
      const pos = errorPositions[i]!;
      // Convert error value from extension field to base field
      const errorVal = this.baseField.__call__(errorValues[i]!.lift.getCoeff(0).value);
      corrected[pos] = corrected[pos]!.sub(errorVal);
    }

    return corrected;
  }

  /**
   * Peterson-Gorenstein-Zierler algorithm to find the error locator polynomial.
   *
   * Builds the syndrome matrix and solves the linear system to find Lambda.
   */
  private pgzErrorLocator(
    syndromes: FiniteFieldElement[],
    maxErrors: number
  ): Polynomial<FiniteFieldElement> | null {
    const extPolyRing = new PolynomialRing(this.splittingField, 'y');

    // Try to find Lambda for v errors, v = maxErrors, maxErrors-1, ..., 1
    for (let v = maxErrors; v >= 1; v--) {
      // Build the syndrome matrix S (v x v)
      // S[i][j] = S_{i+j} for i,j = 0,...,v-1
      const S: FiniteFieldElement[][] = [];

      for (let i = 0; i < v; i++) {
        const row: FiniteFieldElement[] = [];
        for (let j = 0; j < v; j++) {
          if (i + j < syndromes.length) {
            row.push(syndromes[i + j]!);
          } else {
            row.push(this.splittingField.zero());
          }
        }
        S.push(row);
      }

      // Solve S * c = -s where s = [S_v, S_{v+1}, ..., S_{2v-1}]^T
      // and c = [Lambda_v, Lambda_{v-1}, ..., Lambda_1]^T

      const s: FiniteFieldElement[] = [];
      for (let i = 0; i < v; i++) {
        if (v + i < syndromes.length) {
          s.push(syndromes[v + i]!.neg());
        } else {
          s.push(this.splittingField.zero());
        }
      }

      const lambdaCoeffs = this.solveLinearSystem(S, s);

      if (lambdaCoeffs !== null) {
        // Build Lambda(x) = 1 + Lambda_1 * x + ... + Lambda_v * x^v
        const coeffs: FiniteFieldElement[] = [this.splittingField.one()];
        for (let i = v - 1; i >= 0; i--) {
          coeffs.push(lambdaCoeffs[i]!);
        }

        const Lambda = new Polynomial(coeffs, extPolyRing);

        // Verify: check that this Lambda is consistent
        if (Lambda.degree() === v) {
          return Lambda;
        }
      }
    }

    return null;
  }

  /**
   * Solve a linear system Ax = b using Gaussian elimination.
   * Returns null if the system is singular.
   */
  private solveLinearSystem(
    A: FiniteFieldElement[][],
    b: FiniteFieldElement[]
  ): FiniteFieldElement[] | null {
    const n = A.length;
    if (n === 0) return [];

    // Create augmented matrix
    const M: FiniteFieldElement[][] = [];
    for (let i = 0; i < n; i++) {
      M.push([...A[i]!, b[i]!]);
    }

    // Forward elimination with partial pivoting
    for (let col = 0; col < n; col++) {
      // Find pivot
      let pivotRow = -1;
      for (let row = col; row < n; row++) {
        if (!M[row]![col]!.isZero()) {
          pivotRow = row;
          break;
        }
      }

      if (pivotRow === -1) {
        return null; // Singular matrix
      }

      // Swap rows
      [M[col], M[pivotRow]] = [M[pivotRow]!, M[col]!];

      // Scale pivot row
      const pivot = M[col]![col]!;
      const pivotInv = pivot.inv();
      for (let j = col; j <= n; j++) {
        M[col]![j] = M[col]![j]!.mul(pivotInv);
      }

      // Eliminate column
      for (let row = 0; row < n; row++) {
        if (row === col) continue;
        const factor = M[row]![col]!;
        if (!factor.isZero()) {
          for (let j = col; j <= n; j++) {
            M[row]![j] = M[row]![j]!.sub(factor.mul(M[col]![j]!));
          }
        }
      }
    }

    // Extract solution
    const x: FiniteFieldElement[] = [];
    for (let i = 0; i < n; i++) {
      x.push(M[i]![n]!);
    }

    return x;
  }

  /**
   * Chien search to find the roots of the error locator polynomial.
   *
   * Tests each element alpha^{-i} for i = 0, 1, ..., n-1 to find roots.
   * The error positions are the indices i where Lambda(alpha^{-i}) = 0.
   */
  private chienSearch(Lambda: Polynomial<FiniteFieldElement>): number[] {
    const alpha = this.primitiveRoot();
    const errorPositions: number[] = [];

    for (let i = 0; i < this.n; i++) {
      // Test if alpha^{-i} is a root
      const alphaInvI = alpha.pow(-i);
      const val = Lambda.evaluate(alphaInvI);

      if (val.isZero()) {
        errorPositions.push(i);
      }
    }

    return errorPositions;
  }

  /**
   * Forney's algorithm to compute error values.
   *
   * e_i = - Omega(X_i^{-1}) / Lambda'(X_i^{-1})
   *
   * where X_i = alpha^{position_i} is the error locator.
   */
  private forneyAlgorithm(
    syndromes: FiniteFieldElement[],
    Lambda: Polynomial<FiniteFieldElement>,
    positions: number[]
  ): FiniteFieldElement[] {
    const alpha = this.primitiveRoot();
    const extPolyRing = Lambda.parent;

    // Build syndrome polynomial S(x) = S_0 + S_1*x + S_2*x^2 + ...
    const S = new Polynomial(syndromes, extPolyRing);

    // Compute error evaluator: Omega(x) = S(x) * Lambda(x) mod x^{2t}
    const product = S.mul(Lambda);
    const Omega = product.truncate(this.delta - 1);

    // Compute derivative of Lambda
    const LambdaDeriv = Lambda.derivative();

    // Compute error values
    const errorValues: FiniteFieldElement[] = [];

    for (const pos of positions) {
      const Xi = alpha.pow(pos); // error locator
      const XiInv = Xi.inv();

      const omega = Omega.evaluate(XiInv);
      const lambdaPrime = LambdaDeriv.evaluate(XiInv);

      if (lambdaPrime.isZero()) {
        throw new DecodingError('Forney algorithm failed: derivative is zero');
      }

      // In general: e = -X_i * Omega(X_i^{-1}) / Lambda'(X_i^{-1})
      // For narrow-sense codes with b=1, simplified formula
      const errorValue = omega.neg().mul(lambdaPrime.inv());
      errorValues.push(errorValue);
    }

    return errorValues;
  }

  /**
   * Decode and return the original message.
   *
   * @param received - Array of n field elements (possibly corrupted codeword)
   * @returns Array of k field elements (decoded message)
   */
  decode_to_message(received: PrimeFieldElement[]): PrimeFieldElement[] {
    const codeword = this.decode(received);

    // Extract message from systematic encoding
    // Message is in positions n-k to n-1
    const k = this.dimension();
    return codeword.slice(this.n - k);
  }

  /**
   * String representation.
   */
  toString(): string {
    return `[${this.n}, ${this.dimension()}] BCH Code over GF(${this.baseField.characteristic}) with designed distance ${this.delta}`;
  }
}

/**
 * Create a BCH code with the given parameters.
 *
 * @param n - Code length
 * @param delta - Designed distance
 * @param q - Base field size (must be prime)
 * @param narrowSense - If true (default), use narrow-sense BCH (offset=1)
 * @returns The BCH code
 */
export function createBCHCode(
  n: IntegerLike,
  delta: IntegerLike,
  q: number | bigint,
  narrowSense: boolean = true
): BCHCode {
  const nBig = toBigInt(n);
  const deltaBig = toBigInt(delta);
  const baseField = new PrimeField(q);
  const offset = narrowSense ? 1n : 0n;
  return new BCHCode(nBig, deltaBig, baseField, offset);
}

/**
 * Create a primitive BCH code.
 *
 * A primitive BCH code has length n = q^m - 1.
 *
 * @param m - Extension degree (length = q^m - 1)
 * @param delta - Designed distance
 * @param q - Base field size (must be prime)
 * @returns The primitive BCH code
 */
export function createPrimitiveBCHCode(
  m: IntegerLike,
  delta: IntegerLike,
  q: number | bigint
): BCHCode {
  const mNum = toSafeNumber(toBigInt(m));
  const deltaNum = toSafeNumber(toBigInt(delta));
  const qBig = typeof q === 'number' ? BigInt(q) : q;
  const nBig = qBig ** BigInt(mNum) - 1n;
  const baseField = new PrimeField(q);
  return new BCHCode(nBig, BigInt(deltaNum), baseField, 1n);
}

/**
 * Check if a BCH code is actually a Reed-Solomon code.
 *
 * A BCH code is a Reed-Solomon code if:
 * - n = q - 1 (where q is the base field size)
 * - The splitting field equals the base field (splitting degree = 1)
 *
 * @param bch - The BCH code to check
 * @returns true if the BCH code is a Reed-Solomon code
 */
export function isReedSolomonCode(bch: BCHCode): boolean {
  const q = bch.baseField.order;
  const n = bch.n;

  // RS code has n = q - 1 and the splitting field is the base field
  // i.e., the relative extension degree is 1
  return n === Number(q - 1n) && bch.splittingDegree() === 1;
}
