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
 * Greatest common divisor of a non-negative number and a bigint.
 */
function gcdNumber(a: number, b: bigint): bigint {
  let x = BigInt(a);
  let y = b;
  while (y !== 0n) {
    [x, y] = [y, x % y];
  }
  return x < 0n ? -x : x;
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
  private _minimumDistance: number | null = null;
  private _polynomialRing: PolynomialRing<PrimeFieldElement | FiniteFieldElement>;
  private _embeddingImage: FiniteFieldElement | null = null;
  private _sectionBasis: PrimeFieldElement[][] | null = null;

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

    // `bch_code.py:126-129`:
    //     s = Zmod(length)(q).multiplicative_order()
    //     if gcd(jump_size, q ** s - 1) != 1:
    //         raise ValueError("jump_size must be coprime with the order of "
    //                          "the multiplicative group of the splitting field")
    if (gcdNumber(Math.abs(jumpSizeNum), this.splittingField.order - 1n) !== 1n) {
      throw new ValueError(
        'jump_size must be coprime with the order of ' +
          'the multiplicative group of the splitting field'
      );
    }

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
   * Return the arithmetic progression {b, b+l, ..., b+(delta-2)l} mod n that
   * is handed to the underlying cyclic code (`bch_code.py:130-131`).
   */
  private rawDefiningSet(): number[] {
    const D: number[] = [];
    for (let i = 0; i < this.delta - 1; i++) {
      D.push((((this.offset + this.jumpSize * i) % this.n) + this.n) % this.n);
    }
    return D;
  }

  /**
   * Return the defining set of the code: the sorted union of the
   * `q`-cyclotomic cosets modulo `n` of `{b, b+l, ..., b+(delta-2)l}`.
   *
   * Port of `cyclic_code.py:440-452`:
   *
   *     cosets = Zmod(n).cyclotomic_cosets(q, D)
   *     pows = [item for l in cosets for item in l]
   *     self._defining_set = sorted(pows)
   *
   * The defining set is the full set of exponents `i` with `g(alpha^i) = 0`,
   * not just the arithmetic progression: for `n = 15`, `delta = 5`, `q = 2`
   * Sage returns `[1, 2, 3, 4, 6, 8, 9, 12]`, not `[1, 2, 3, 4]`.
   */
  defining_set(): number[] {
    if (this._definingSet !== null) {
      return this._definingSet;
    }

    const q = Number(this.baseField.order);
    const seen = new Set<number>();

    for (const i of this.rawDefiningSet()) {
      let j = i;
      while (!seen.has(j)) {
        seen.add(j);
        j = (j * q) % this.n;
      }
    }

    this._definingSet = [...seen].sort((a, b) => a - b);
    return this._definingSet;
  }

  /**
   * Return the image of the base field generator under the canonical
   * embedding of the base field into the splitting field.
   *
   * SageMath builds the splitting field as `Fsplit, FE = F.extension(s,
   * map=True)` (`cyclic_code.py:434-437`) and uses the honest ring morphism
   * `FE` and its `section()` to move between `F` and `Fsplit`.
   * `finite_field_base.pyx:1505-1508` defines that morphism by
   *
   *     alpha = E.gen()**((E.order() - 1) // (F.order() - 1))   # both Conway
   *     alpha = F.modulus().any_root(E)                         # otherwise
   *
   * and then `F.hom([alpha], codomain=E)`.  We reproduce exactly that: try
   * the Conway power first (our default moduli are Conway polynomials
   * whenever a Conway polynomial is tabulated), and fall back to a search for
   * a root of the base modulus.
   */
  private embeddingImage(): FiniteFieldElement {
    if (this._embeddingImage !== null) {
      return this._embeddingImage;
    }

    const E = this.splittingField;

    if (this.baseField instanceof PrimeField) {
      // GF(p) -> GF(p^d) sends 1 to 1; the "generator" is 1.
      this._embeddingImage = E.one();
      return this._embeddingImage;
    }

    const F = this.baseField as FiniteFieldExtension;
    const modulus = F.modulus;

    const evaluateModulus = (y: FiniteFieldElement): FiniteFieldElement => {
      let acc = E.zero();
      for (let i = modulus.degree(); i >= 0; i--) {
        acc = acc.mul(y).add(E.__call__(modulus.getCoeff(i).value));
      }
      return acc;
    };

    // Conway-compatible power, as in `finite_field_base.pyx:1506`
    const candidate = E.gen().pow((E.order - 1n) / (F.order - 1n));
    if (evaluateModulus(candidate).isZero()) {
      this._embeddingImage = candidate;
      return candidate;
    }

    // Fallback: any root of the base modulus in the splitting field
    // (Sage's `self.modulus().any_root(E)`).
    if (E.order > 1n << 22n) {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: root finding over a splitting field of order ' +
          `${E.order} (no Conway-compatible embedding available)`
      );
    }
    for (const y of E) {
      if (evaluateModulus(y).isZero()) {
        this._embeddingImage = y;
        return y;
      }
    }

    throw new ValueError('base field does not embed into the splitting field');
  }

  /**
   * Embed an element of the base field into the splitting field.
   *
   * This is SageMath's `field_embedding()` of the surrounding cyclic code
   * (`cyclic_code.py:451`).
   */
  fieldEmbedding(c: PrimeFieldElement | FiniteFieldElement): FiniteFieldElement {
    const E = this.splittingField;

    if (this.baseField instanceof PrimeField) {
      return E.__call__((c as PrimeFieldElement).value);
    }

    const alpha = this.embeddingImage();
    const F = this.baseField as FiniteFieldExtension;
    const coeffs = (c as FiniteFieldElement).lift;

    let acc = E.zero();
    for (let i = F.degree - 1; i >= 0; i--) {
      acc = acc.mul(alpha).add(E.__call__(coeffs.getCoeff(i).value));
    }
    return acc;
  }

  /**
   * The section of {@link fieldEmbedding}: map an element of the image
   * subfield back to the base field.
   *
   * This is SageMath's `C.field_embedding().section()`
   * (`cyclic_code.py:445`, `bch_code.py:417`).  Concretely: express the
   * element in the `GF(p)`-basis `1, alpha, ..., alpha^(k-1)` of the image,
   * where `alpha = ` {@link embeddingImage} and `k = [F : GF(p)]`.
   *
   * @throws {ValueError} If the element does not lie in the image subfield
   */
  fieldEmbeddingSection(y: FiniteFieldElement): PrimeFieldElement | FiniteFieldElement {
    if (this.baseField instanceof PrimeField) {
      const coeffs = y.coefficients();
      for (let i = 1; i < coeffs.length; i++) {
        if (!coeffs[i]!.isZero()) {
          throw new ValueError('element is not in the image of the field embedding');
        }
      }
      return this.baseField.__call__(coeffs[0]!.value);
    }

    const F = this.baseField as FiniteFieldExtension;
    const basis = this.sectionBasis();
    const p = this.splittingField.characteristic;
    const d = this.splittingField.degree;
    const k = F.degree;

    // Solve  sum_j c_j * basis[j] = y  over GF(p)  (basis[j] = coords of alpha^j)
    const M: bigint[][] = [];
    const target = y.coefficients();
    for (let row = 0; row < d; row++) {
      const r: bigint[] = [];
      for (let j = 0; j < k; j++) {
        r.push(basis[j]![row]!.value);
      }
      r.push(target[row]!.value);
      M.push(r);
    }

    const modInv = (a: bigint): bigint => {
      let [old_r, r] = [((a % p) + p) % p, p];
      let [old_s, s] = [1n, 0n];
      while (r !== 0n) {
        const q = old_r / r;
        [old_r, r] = [r, old_r - q * r];
        [old_s, s] = [s, old_s - q * s];
      }
      return ((old_s % p) + p) % p;
    };

    const pivotOf: number[] = [];
    let row = 0;
    for (let col = 0; col < k && row < d; col++) {
      let pivot = -1;
      for (let i = row; i < d; i++) {
        if (M[i]![col] !== 0n) {
          pivot = i;
          break;
        }
      }
      if (pivot === -1) continue;
      [M[row], M[pivot]] = [M[pivot]!, M[row]!];
      const inv = modInv(M[row]![col]!);
      for (let j = col; j <= k; j++) {
        M[row]![j] = (M[row]![j]! * inv) % p;
      }
      for (let i = 0; i < d; i++) {
        if (i === row) continue;
        const f = M[i]![col]!;
        if (f === 0n) continue;
        for (let j = col; j <= k; j++) {
          M[i]![j] = (((M[i]![j]! - f * M[row]![j]!) % p) + p) % p;
        }
      }
      pivotOf[col] = row;
      row++;
    }

    // Any row with all-zero coefficients but non-zero right-hand side means
    // y is outside the image.
    for (let i = 0; i < d; i++) {
      let allZero = true;
      for (let j = 0; j < k; j++) {
        if (M[i]![j] !== 0n) {
          allZero = false;
          break;
        }
      }
      if (allZero && M[i]![k] !== 0n) {
        throw new ValueError('element is not in the image of the field embedding');
      }
    }

    const coeffs: number[] = [];
    for (let j = 0; j < k; j++) {
      const r = pivotOf[j];
      coeffs.push(r === undefined ? 0 : Number(M[r]![k]!));
    }

    return F.__call__(coeffs);
  }

  /**
   * Coordinate vectors of `1, alpha, ..., alpha^(k-1)` over `GF(p)`.
   */
  private sectionBasis(): PrimeFieldElement[][] {
    if (this._sectionBasis !== null) {
      return this._sectionBasis;
    }

    const alpha = this.embeddingImage();
    const k = (this.baseField as FiniteFieldExtension).degree;
    const basis: PrimeFieldElement[][] = [];
    let power = this.splittingField.one();
    for (let i = 0; i < k; i++) {
      basis.push(power.coefficients());
      power = power.mul(alpha);
    }

    this._sectionBasis = basis;
    return basis;
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

    // The minimal polynomial has coefficients in the image of the base field
    // under the canonical embedding; pull them back with its section, exactly
    // as `cyclic_code.py:444-449` does (`g *= R([sec(coeff) for coeff in pol])`).
    // Truncating the splitting-field coordinate vector to `[F : GF(p)]`
    // entries is *not* the section unless the base field is prime, and
    // produces a generator polynomial that does not divide x^n - 1.
    const baseCoeffs: (PrimeFieldElement | FiniteFieldElement)[] = [];
    for (let deg = 0; deg <= minPoly.degree(); deg++) {
      baseCoeffs.push(this.fieldEmbeddingSection(minPoly.getCoeff(deg)));
    }

    return new Polynomial(baseCoeffs, this._polynomialRing);
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
   * Return the true minimum distance of the code.
   *
   * SageMath's `BCHCode` inherits `AbstractLinearCode.minimum_distance`,
   * which returns the *actual* minimum distance, not the designed one: the
   * `[23, 12]` binary Golay code `BCHCode(GF(2), 23, 5)` has designed
   * distance 5 but minimum distance 7.  Use {@link designed_distance} for the
   * BCH bound.
   *
   * Sage delegates to GAP/Guava's Brouwer-Zimmermann implementation, which we
   * have no port of; here the weights of all `q^k` codewords are enumerated,
   * which is exact but only feasible for small codes.
   *
   * @throws {NotImplementedError} If `q^k` exceeds the enumeration budget
   * @see Deviation: minimum distance by enumeration instead of Brouwer-Zimmermann
   */
  minimum_distance(): number {
    if (this._minimumDistance !== null) {
      return this._minimumDistance;
    }

    const k = this.dimension();
    const q = this.baseField.order;
    const total = q ** BigInt(k);

    if (k === 0) {
      throw new ValueError('the zero code has no minimum distance');
    }

    const BUDGET = 1n << 17n;
    if (total > BUDGET) {
      throw new NotImplementedError(
        `SAGE_NOT_IMPLEMENTED: minimum_distance for a [${this.n}, ${k}] code over ` +
          `GF(${q}) (${total} codewords exceeds the enumeration budget ${BUDGET}); ` +
          'use designed_distance() for the BCH bound'
      );
    }

    const elements = this.baseFieldElements();
    const zero = this.baseField.zero();
    let best = this.n + 1;

    const digits = new Array<number>(k).fill(0);
    for (let counter = 1n; counter < total; counter++) {
      // increment the base-q counter
      for (let i = 0; i < k; i++) {
        digits[i] = digits[i]! + 1;
        if (digits[i]! < elements.length) break;
        digits[i] = 0;
      }

      const message = digits.map((d) => elements[d]!);
      const codeword = this.encode(message);
      let weight = 0;
      for (const c of codeword) {
        if (!c.eq(zero as never)) weight++;
      }
      if (weight > 0 && weight < best) best = weight;
    }

    this._minimumDistance = best;
    return best;
  }

  /**
   * The elements of the base field, `0` first.
   */
  private baseFieldElements(): (PrimeFieldElement | FiniteFieldElement)[] {
    const q = Number(this.baseField.order);
    const out: (PrimeFieldElement | FiniteFieldElement)[] = [];
    if (this.baseField instanceof PrimeField) {
      for (let i = 0; i < q; i++) out.push(this.baseField.__call__(i));
    } else {
      for (let i = 0n; i < this.baseField.order; i++) {
        out.push((this.baseField as FiniteFieldExtension).fromInteger(i));
      }
    }
    return out;
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
  encode<T extends PrimeFieldElement | FiniteFieldElement>(message: T[]): T[] {
    const k = this.dimension();

    if (message.length !== k) {
      throw new ValueError(`message length must be ${k}, got ${message.length}`);
    }

    const g = this.generator_polynomial();

    // Create message polynomial m(x)
    const m = new Polynomial(
      message as (PrimeFieldElement | FiniteFieldElement)[],
      this._polynomialRing
    );

    // Compute m(x) * x^{n-k}
    const shifted = m.shift(this.n - k);

    // Compute remainder: r(x) = m(x) * x^{n-k} mod g(x)
    const [_quotient, remainder] = shifted.quo_rem(g);

    // Codeword: c(x) = m(x) * x^{n-k} - r(x)
    const codewordPoly = shifted.sub(remainder);

    // Extract coefficients
    const codeword: (PrimeFieldElement | FiniteFieldElement)[] = [];
    for (let i = 0; i < this.n; i++) {
      codeword.push(codewordPoly.getCoeff(i));
    }

    return codeword as T[];
  }

  /**
   * Compute the syndromes of a received word.
   *
   * S_j = r(alpha^{b + l*j}) for j = 0, 1, ..., delta-2
   *
   * The received symbols are moved into the splitting field through the
   * canonical embedding (`bch_code.py:382-401`, `bch_word_to_grs`), which is
   * the identity on a prime base field but a genuine morphism otherwise.
   *
   * @param received - Array of n field elements
   * @returns Array of delta-1 syndrome values in the splitting field
   */
  syndrome(received: (PrimeFieldElement | FiniteFieldElement)[]): FiniteFieldElement[] {
    if (received.length !== this.n) {
      throw new ValueError(`received length must be ${this.n}, got ${received.length}`);
    }

    const alpha = this.primitiveRoot();
    const syndromes: FiniteFieldElement[] = [];

    // Create received polynomial in the splitting field
    const extPolyRing = new PolynomialRing(this.splittingField, 'x');
    const rCoeffs = received.map((c) => this.fieldEmbedding(c));
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
  is_codeword(word: (PrimeFieldElement | FiniteFieldElement)[]): boolean {
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
   * 4. Correct errors at those positions with Forney's formula
   *
   * The result is checked to be a codeword, mirroring
   * `BCHUnderlyingGRSDecoder.decode_to_code`, which discards GRS decodings
   * that fall outside the BCH code (`bch_code.py:437-470`).
   *
   * @param received - Array of n field elements (possibly corrupted codeword)
   * @returns Array of n field elements (corrected codeword)
   * @throws {DecodingError} If too many errors to correct
   */
  decode<T extends PrimeFieldElement | FiniteFieldElement>(received: T[]): T[] {
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

    const corrected: (PrimeFieldElement | FiniteFieldElement)[] = [...received];

    if (this.baseField.order === 2n) {
      // Over GF(2) the only non-zero error value is 1.
      for (const pos of errorPositions) {
        corrected[pos] = corrected[pos]!.add(this.baseField.one() as never);
      }
    } else {
      // Forney's algorithm gives the error values in the splitting field;
      // they lie in the image of the base field, so pull them back with the
      // section of the field embedding.
      const errorValues = this.forneyAlgorithm(syndromes, Lambda, errorPositions);

      for (let i = 0; i < errorPositions.length; i++) {
        const pos = errorPositions[i]!;
        let errorVal: PrimeFieldElement | FiniteFieldElement;
        try {
          errorVal = this.fieldEmbeddingSection(errorValues[i]!);
        } catch {
          throw new DecodingError('decoding failed: error value is not in the base field');
        }
        corrected[pos] = corrected[pos]!.sub(errorVal as never);
      }
    }

    if (!this.is_codeword(corrected)) {
      throw new DecodingError('decoding failed: corrected word is not a codeword');
    }

    return corrected as T[];
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
   * Chien search to find the error positions.
   *
   * With jump size `l` the syndromes are `S_j = sum_i Y_i Z_i^j` where
   * `Z_i = alpha^{l * p_i}`, so the error locator vanishes at
   * `alpha^{-l * p}` for each error position `p` -- not at `alpha^{-p}`.
   * Testing `alpha^{-i}` and reporting `i` gives `l * p mod n` instead of `p`.
   */
  private chienSearch(Lambda: Polynomial<FiniteFieldElement>): number[] {
    const alpha = this.primitiveRoot();
    const errorPositions: number[] = [];

    for (let p = 0; p < this.n; p++) {
      const zInv = alpha.pow(-((this.jumpSize * p) % this.n));
      if (Lambda.evaluate(zInv).isZero()) {
        errorPositions.push(p);
      }
    }

    return errorPositions;
  }

  /**
   * Forney's algorithm to compute error values.
   *
   * With `X_i = alpha^{p_i}` and `Z_i = X_i^l`, the syndrome sequence is
   * `S_j = sum_i Y_i Z_i^j` with `Y_i = e_i X_i^b`.  Forney's formula gives
   * `Y_i = -Z_i * Omega(Z_i^{-1}) / Lambda'(Z_i^{-1})`, hence
   *
   *     e_i = Y_i * X_i^{-b} = -X_i^{l-b} * Omega(Z_i^{-1}) / Lambda'(Z_i^{-1}).
   *
   * Omitting the `X_i^{l-b}` factor (i.e. assuming the narrow-sense case
   * `b = l = 1`) scales every error value by the wrong power of `X_i`.
   * SageMath sidesteps this by decoding through the underlying GRS code
   * (`bch_code.py:239-254`, `bch_to_grs`), whose column multipliers carry
   * exactly the `alpha^{b*i}` factors.
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

    // Compute error evaluator: Omega(x) = S(x) * Lambda(x) mod x^{delta-1}
    const product = S.mul(Lambda);
    const Omega = product.truncate(this.delta - 1);

    // Compute derivative of Lambda
    const LambdaDeriv = Lambda.derivative();

    // Compute error values
    const errorValues: FiniteFieldElement[] = [];

    for (const pos of positions) {
      const Zi = alpha.pow((this.jumpSize * pos) % this.n);
      const ZiInv = Zi.inv();

      const omega = Omega.evaluate(ZiInv);
      const lambdaPrime = LambdaDeriv.evaluate(ZiInv);

      if (lambdaPrime.isZero()) {
        throw new DecodingError('Forney algorithm failed: derivative is zero');
      }

      // e = -X_i^{l - b} * Omega(Z_i^-1) / Lambda'(Z_i^-1)
      const factor = alpha.pow(((this.jumpSize - this.offset) * pos) % this.n);
      const errorValue = factor.mul(omega).mul(lambdaPrime.inv()).neg();
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
  decode_to_message<T extends PrimeFieldElement | FiniteFieldElement>(received: T[]): T[] {
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
