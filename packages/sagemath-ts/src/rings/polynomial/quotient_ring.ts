/**
 * @module sage/rings/polynomial/quotient_ring
 * @description Quotient rings R[x]/<f(x)>
 *
 * When f(x) is irreducible over a field, this constructs an extension field.
 */

import { NotImplementedError, ValueError, ZeroDivisionError } from '../../errors.js';
import { type CoefficientRing, Polynomial, type RingElement } from './polynomial_element.js';
import type { PolynomialRing } from './polynomial_ring.js';

/**
 * Element of a quotient ring R[x]/<f(x)>.
 *
 * Represented as a polynomial of degree < deg(f).
 */
export class QuotientRingElement<C extends RingElement> implements RingElement {
  readonly lift: Polynomial<C>;
  readonly parent: QuotientRing<C>;

  constructor(poly: Polynomial<C>, parent: QuotientRing<C>) {
    this.parent = parent;

    // Reduce modulo the modulus
    if (poly.degree() >= parent.modulus.degree()) {
      const [_q, r] = poly.quo_rem(parent.modulus);
      this.lift = r;
    } else {
      this.lift = poly;
    }
  }

  add(other: QuotientRingElement<C>): QuotientRingElement<C> {
    return new QuotientRingElement(this.lift.add(other.lift), this.parent);
  }

  sub(other: QuotientRingElement<C>): QuotientRingElement<C> {
    return new QuotientRingElement(this.lift.sub(other.lift), this.parent);
  }

  mul(other: QuotientRingElement<C>): QuotientRingElement<C> {
    const prod = this.lift.mul(other.lift);
    return new QuotientRingElement(prod, this.parent);
  }

  neg(): QuotientRingElement<C> {
    return new QuotientRingElement(this.lift.neg(), this.parent);
  }

  /**
   * Compute the multiplicative inverse using extended Euclidean algorithm.
   */
  inv(): QuotientRingElement<C> {
    if (this.isZero()) {
      throw new ZeroDivisionError('division by zero');
    }

    // Extended Euclidean algorithm for polynomials
    const [g, s, _t] = polyXgcd(this.lift, this.parent.modulus);

    // g should be a unit (constant) if modulus is irreducible
    if (g.degree() !== 0) {
      // Sage raises a ZeroDivisionError with this message
      // (`polynomial_quotient_ring_element.py:__invert__`).
      throw new ZeroDivisionError(
        `element ${this.lift} of quotient polynomial ring not invertible`
      );
    }

    // s * this ≡ g (mod modulus)
    // We need to divide s by g's constant term
    const gInv = g.getCoeff(0);
    const sNormalized = s.scalar_mul(divideCoeffs(s.parent.base_ring.one() as C, gInv));

    return new QuotientRingElement(sNormalized, this.parent);
  }

  div(other: QuotientRingElement<C>): QuotientRingElement<C> {
    return this.mul(other.inv());
  }

  pow(n: number | bigint): QuotientRingElement<C> {
    let exp = typeof n === 'bigint' ? n : BigInt(n);

    if (exp < 0n) {
      return this.inv().pow(-exp);
    }

    if (exp === 0n) {
      return this.parent.one();
    }

    // Binary exponentiation
    let result = this.parent.one();
    let base: QuotientRingElement<C> = this;

    while (exp > 0n) {
      if ((exp & 1n) === 1n) {
        result = result.mul(base);
      }
      base = base.mul(base);
      exp >>= 1n;
    }

    return result;
  }

  eq(other: QuotientRingElement<C> | number): boolean {
    if (typeof other === 'number') {
      const otherElem = this.parent.__call__(other);
      return this.lift.eq(otherElem.lift);
    }
    return this.lift.eq(other.lift);
  }

  isZero(): boolean {
    return this.lift.isZero();
  }

  isOne(): boolean {
    return this.lift.coeffs.length === 1 && this.lift.coeffs[0]!.eq(1);
  }

  toString(): string {
    return this.lift.toString();
  }

  repr(): string {
    return this.lift.toString();
  }
}

/**
 * Quotient ring R[x]/<f(x)>.
 */
export class QuotientRing<C extends RingElement>
  implements CoefficientRing<QuotientRingElement<C>>
{
  readonly polynomial_ring: PolynomialRing<C>;
  readonly modulus: Polynomial<C>;
  readonly degree: number;

  // Cache for iteration
  private _elements: QuotientRingElement<C>[] | null = null;
  // Cache for is_field()
  private _isField: boolean | null = null;

  constructor(polynomial_ring: PolynomialRing<C>, modulus: Polynomial<C>) {
    if (modulus.isZero()) {
      throw new ValueError('modulus cannot be zero');
    }
    // Degree-0 moduli are legal in Sage: ``R.quo(1)`` is the zero ring
    // (`polynomial_quotient_ring.py:776`, cardinality 1).

    this.polynomial_ring = polynomial_ring;
    this.modulus = modulus;
    this.degree = modulus.degree();
  }

  /**
   * Create an element from a polynomial or coefficient.
   */
  __call__(
    x: Polynomial<C> | QuotientRingElement<C> | C | number | unknown
  ): QuotientRingElement<C> {
    if (x instanceof QuotientRingElement) {
      return new QuotientRingElement(x.lift, this);
    }
    if (x instanceof Polynomial) {
      return new QuotientRingElement(x, this);
    }
    if (typeof x === 'number') {
      const coeff = this.polynomial_ring.base_ring.__call__(x) as C;
      return new QuotientRingElement(this.polynomial_ring.__call__(coeff), this);
    }
    // Assume it's a coefficient
    return new QuotientRingElement(this.polynomial_ring.__call__(x as C), this);
  }

  /**
   * Return zero element.
   */
  zero(): QuotientRingElement<C> {
    return new QuotientRingElement(this.polynomial_ring.zero(), this);
  }

  /**
   * Return one element.
   */
  one(): QuotientRingElement<C> {
    return new QuotientRingElement(this.polynomial_ring.one(), this);
  }

  /**
   * Return the generator (image of x in the quotient).
   */
  gen(): QuotientRingElement<C> {
    return new QuotientRingElement(this.polynomial_ring.gen(), this);
  }

  /**
   * Return the number of elements of this quotient ring.
   *
   * Returns `Infinity` when the quotient is infinite, mirroring Sage's
   * `+Infinity` (`polynomial_quotient_ring.py:767`).
   *
   * @see Reference: sage/rings/polynomial/polynomial_quotient_ring.py:767 (cardinality)
   */
  cardinality(): bigint | number {
    // R[x]/(c) with c a nonzero constant is the zero ring: one element.
    if (this.degree === 0) {
      return 1n;
    }
    const base = this.polynomial_ring.base_ring as { cardinality?: () => bigint };
    const baseCard = base.cardinality?.();
    if (baseCard === undefined || baseCard === 0n) {
      return Number.POSITIVE_INFINITY;
    }
    return baseCard ** BigInt(this.degree);
  }

  /**
   * Alias for {@link cardinality} (Sage: `order = cardinality`).
   */
  order(): bigint | number {
    return this.cardinality();
  }

  /**
   * Return whether this quotient ring is finite
   * (`polynomial_quotient_ring.py:819`).
   */
  is_finite(): boolean {
    if (this.degree === 0) {
      return true;
    }
    const base = this.polynomial_ring.base_ring as {
      is_finite?: () => boolean;
      cardinality?: () => bigint;
    };
    if (typeof base.is_finite === 'function') {
      return base.is_finite();
    }
    return base.cardinality?.() !== undefined;
  }

  /**
   * Iterate over all elements (only works for finite quotient rings).
   */
  *[Symbol.iterator](): Iterator<QuotientRingElement<C>> {
    if (this._elements !== null) {
      yield* this._elements;
      return;
    }

    // Generate all elements
    this._elements = [];

    // Get base ring elements
    const baseRing = this.polynomial_ring.base_ring;
    if (!baseRing[Symbol.iterator]) {
      throw new ValueError('base ring is not iterable');
    }

    const baseElements = [...(baseRing as Iterable<C>)];

    // Generate all polynomials of degree < modulus degree
    const generate = (coeffs: C[]): void => {
      if (coeffs.length === this.degree) {
        const poly = new Polynomial(coeffs, this.polynomial_ring);
        const elem = new QuotientRingElement(poly, this);
        this._elements!.push(elem);
        return;
      }

      // Prepend, so that the CONSTANT coefficient varies fastest.  Upstream
      // `polynomial_quotient_ring.py:854-880` yields `zero()` and then, for
      // i = 0 .. deg-1, `R.polynomials(of_degree=i)`, which enumerates the
      // leading coefficient slowest and the constant term fastest; the overall
      // sequence is therefore plain little-endian numeric order
      // (`list(GF(3)[x].quo(x^3-x^2-x-1))` == [0, 1, 2, xbar, xbar+1, ...]).
      for (const c of baseElements) {
        generate([c, ...coeffs]);
      }
    };

    generate([]);

    yield* this._elements;
  }

  /**
   * Number of elements.
   */
  get length(): number {
    // Force computation if needed
    if (this._elements === null) {
      [...this];
    }
    return this._elements!.length;
  }

  /**
   * Check if this is a field: the base ring must be a field and the modulus
   * must be irreducible.
   *
   * Sage: `base_ring().is_field(proof) and modulus().is_irreducible()`, with
   * `proof=False` turning a `NotImplementedError` from the irreducibility
   * test into `False`.
   *
   * @see Reference: sage/rings/polynomial/polynomial_quotient_ring.py:1006 (is_field)
   */
  is_field(proof: boolean = true): boolean {
    if (this._isField !== null) {
      return this._isField;
    }

    const base = this.polynomial_ring.base_ring;
    let ret = typeof base.is_field === 'function' ? base.is_field() : false;

    if (ret) {
      try {
        ret = this.modulus.is_irreducible();
      } catch (e) {
        if (e instanceof NotImplementedError && !proof) {
          ret = false;
        } else {
          throw e;
        }
      }
    }

    this._isField = ret;
    return ret;
  }

  toString(): string {
    return `Univariate Quotient Polynomial Ring in ${this.polynomial_ring.variable_name} over ${this.polynomial_ring.base_ring} with modulus ${this.modulus}`;
  }
}

/**
 * Extended Euclidean algorithm for polynomials.
 * Returns (g, s, t) such that s*a + t*b = g.
 */
function polyXgcd<C extends RingElement>(
  a: Polynomial<C>,
  b: Polynomial<C>
): [Polynomial<C>, Polynomial<C>, Polynomial<C>] {
  const ring = a.parent;

  let oldR = a;
  let r = b;
  let oldS = ring.one();
  let s = ring.zero();
  let oldT = ring.zero();
  let t = ring.one();

  while (!r.isZero()) {
    const [quotient, remainder] = oldR.quo_rem(r);

    const tempR = r;
    r = remainder;
    oldR = tempR;

    const tempS = s;
    s = oldS.sub(quotient.mul(s));
    oldS = tempS;

    const tempT = t;
    t = oldT.sub(quotient.mul(t));
    oldT = tempT;
  }

  return [oldR, oldS, oldT];
}

/**
 * Helper to divide coefficients.
 */
function divideCoeffs<C extends RingElement>(a: C, b: C): C {
  if ('div' in a && typeof (a as unknown as { div: (b: C) => C }).div === 'function') {
    return (a as unknown as { div: (b: C) => C }).div(b);
  }
  if ('inv' in b && typeof (b as unknown as { inv: () => C }).inv === 'function') {
    const bInv = (b as unknown as { inv: () => C }).inv();
    return a.mul(bInv) as C;
  }
  throw new ValueError('coefficient ring does not support division');
}
