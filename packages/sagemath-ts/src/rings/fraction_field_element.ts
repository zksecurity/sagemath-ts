/** Polynomial fractions; source: sage/rings/fraction_field_element.pyx.
 * @see Deviation: Native Modular Polynomial Products and Fraction Fields
 */
import {
  ArithmeticError,
  AttributeError,
  NotImplementedError,
  ValueError,
  ZeroDivisionError,
} from '../errors.js';
import type { FractionField_generic } from './fraction_field.js';
import type { Polynomial, RingElement } from './polynomial/polynomial_element.js';

import type { IntegerRing } from './integer_ring.js';
import type { Rational } from './rational.js';

export type FractionOptions = { coerce?: boolean; reduce?: boolean };
export interface FractionElement<C extends RingElement> {
  readonly _is_fraction_field_element: true;
  readonly parent: FractionField_generic<C>;
  numerator(): Polynomial<C>;
  denominator(): Polynomial<C>;
  add(b: unknown): FractionElement<C>;
  sub(b: unknown): FractionElement<C>;
  mul(b: unknown): FractionElement<C>;
  div(b: unknown): FractionElement<C>;
  neg(): FractionElement<C>;
  inv(): FractionElement<C>;
  pow(n: bigint | number): FractionElement<C>;
  eq(b: unknown): boolean;
  isZero(): boolean;
  is_square(): boolean;
  sqrt(extend?: boolean, all?: boolean): FractionElement<C> | FractionElement<C>[];
  toString(): string;
}
export type FractionElementConstructor<C extends RingElement> = new (
  parent: FractionField_generic<any>,
  n: unknown,
  d?: unknown,
  options?: FractionOptions
) => FractionElement<C>;
export function isFractionElement(x: unknown): x is FractionElement<RingElement> {
  return (
    typeof x === 'object' &&
    x !== null &&
    (x as FractionElement<RingElement>)._is_fraction_field_element === true
  );
}
export class FractionFieldElement<C extends RingElement> implements FractionElement<C> {
  readonly _is_fraction_field_element = true as const;
  protected num: Polynomial<C>;
  protected den: Polynomial<C>;
  protected _is_reduced = false;
  constructor(
    readonly parent: FractionField_generic<C>,
    numerator: unknown,
    denominator: unknown = 1n,
    options: FractionOptions = {}
  ) {
    const R = parent.ring();
    this.num = options.coerce === false ? (numerator as Polynomial<C>) : R.__call__(numerator);
    this.den = options.coerce === false ? (denominator as Polynomial<C>) : R.__call__(denominator);
    if (options.reduce !== false && parent.is_exact()) {
      try {
        this.reduce();
      } catch (e) {
        if (!(e instanceof ArithmeticError || e instanceof ZeroDivisionError)) throw e;
      }
    }
    if (this.den.isZero()) throw new ZeroDivisionError('fraction field element division by zero');
  }
  /** Sage tests numerator times denominator, including unreduced fractions. */
  is_square(root?: false): boolean;
  is_square(root: true): [boolean, FractionElement<C> | null];
  is_square(root: boolean): boolean | [boolean, FractionElement<C> | null];
  is_square(root = false): boolean | [boolean, FractionElement<C> | null] {
    const a = this.numerator(), b = this.denominator(), product = a.mul(b);
    if (!root) return product.is_square();
    const [ok, sq] = product.is_square(true);
    return ok ? [true, this.parent.__call__(sq!, b)] : [false, null];
  }
  /** Inherited commutative-ring root protocol.
   * @see Deviation: Fraction field square roots
   */
  sqrt(extend = true, all = false, name: string | null = null): FractionElement<C> | FractionElement<C>[] {
    const [ok, sq] = this.is_square(true);
    if (ok) return all ? (this.parent.characteristic() === 2n || sq!.isZero() ? [sq!] : [sq!, sq!.neg()]) : sq!;
    if (!extend) {
      if (all) return [];
      throw new ValueError(`trying to take square root of non-square ${this} with extend = False`);
    }
    if (name === null) throw new TypeError('Polynomial is not a square. You must specify the name of the square root when using the default extend = True');
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: fraction field square-root extension');
  }
  reduce(): void {
    if (this._is_reduced) return;
    try {
      const g = this.num.gcd(this.den);
      if (!unit(g)) {
        // Native floor division allocates even for a zero numerator; quo_rem retains it.
        const nativeZero =
          this.num.isZero() &&
          ['Integer Ring', 'Rational Field'].includes(String(this.parent.base_ring()));
        this.num = nativeZero ? this.parent.ring().__call__([]) : this.num.quo_rem(g)[0];
        this.den = this.den.quo_rem(g)[0];
      }
      this._is_reduced = true;
    } catch (e) {
      if (e instanceof AttributeError)
        throw new ArithmeticError('unable to reduce because lack of gcd or quo_rem algorithm');
      if (e instanceof TypeError)
        throw new ArithmeticError("unable to reduce because gcd algorithm doesn't work on input");
      if (e instanceof NotImplementedError)
        throw new ArithmeticError(
          'unable to reduce because gcd algorithm not implemented on input'
        );
      throw e;
    }
    if (!this.den.eq(1n) && unit(this.den)) {
      try {
        this.num = this.num.mul(inverseUnit(this.den));
        this.den = this.parent.ring().one();
      } catch {}
    }
  }
  protected _new(n: Polynomial<C>, d: Polynomial<C>, coerce = false): FractionFieldElement<C> {
    const Class = this.constructor as FractionElementConstructor<C>;
    return new Class(this.parent, n, d, { coerce, reduce: false }) as FractionFieldElement<C>;
  }
  add(other: unknown): FractionElement<C> {
    const b = this.parent.__call__(other),
      sn = b.numerator(),
      sd = b.denominator();
    if (this.num.isZero()) return b;
    if (sn.isZero()) return this;
    if (this.parent.is_exact())
      try {
        const d = this.den.gcd(sd);
        if (unit(d)) return this._new(this.num.mul(sd).add(this.den.mul(sn)), this.den.mul(sd));
        const rd = this.den.quo_rem(d)[0],
          ssd = sd.quo_rem(d)[0];
        let n = this.num.mul(ssd).add(rd.mul(sn));
        if (n.isZero()) return this._new(n, this.parent.ring().one());
        let den = this.den.mul(ssd);
        const e = n.gcd(d);
        if (!unit(e)) {
          n = n.quo_rem(e)[0];
          den = den.quo_rem(e)[0];
        }
        [n, den] = normalizeUnit(n, den);
        return this._new(n, den);
      } catch (e) {
        if (
          !(
            e instanceof AttributeError ||
            e instanceof NotImplementedError ||
            e instanceof TypeError
          )
        )
          throw e;
      }
    return this._new(this.num.mul(sd).add(this.den.mul(sn)), this.den.mul(sd));
  }
  sub(other: unknown): FractionElement<C> {
    return this.add(this.parent.__call__(other).neg());
  }
  mul(other: unknown): FractionElement<C> {
    const b = this.parent.__call__(other);
    let rn = this.num,
      rd = this.den,
      sn = b.numerator(),
      sd = b.denominator();
    if (rn.isZero() || sn.isZero()) return this.parent.zero();
    if (this.parent.is_exact())
      try {
        const d1 = rn.gcd(sd),
          d2 = sn.gcd(rd);
        if (!unit(d1)) {
          rn = rn.quo_rem(d1)[0];
          sd = sd.quo_rem(d1)[0];
        }
        if (!unit(d2)) {
          rd = rd.quo_rem(d2)[0];
          sn = sn.quo_rem(d2)[0];
        }
        const [n, d] = normalizeUnit(rn.mul(sn), rd.mul(sd));
        return this._new(n, d);
      } catch (e) {
        if (
          !(
            e instanceof AttributeError ||
            e instanceof NotImplementedError ||
            e instanceof TypeError
          )
        )
          throw e;
      }
    return this._new(this.num.mul(b.numerator()), this.den.mul(b.denominator()));
  }
  div(other: unknown): FractionElement<C> {
    const b = this.parent.__call__(other);
    if (b.numerator().isZero())
      throw new ZeroDivisionError('fraction field element division by zero');
    return this.mul(this._new(b.denominator(), b.numerator(), true));
  }
  neg(): FractionFieldElement<C> {
    return this._new(this.num.neg(), this.den);
  }
  inv(): FractionFieldElement<C> {
    if (this.isZero()) throw new ZeroDivisionError('Cannot invert 0');
    return this._new(this.den, this.num);
  }
  pow(exponent: bigint | number): FractionFieldElement<C> {
    const n = BigInt(exponent),
      R = this.parent.ring();
    if (n === 0n) return this._new(R.one(), R.one());
    return n > 0n
      ? this._new(this.num.pow(n) as Polynomial<C>, this.den.pow(n) as Polynomial<C>)
      : this._new(this.den.pow(-n) as Polynomial<C>, this.num.pow(-n) as Polynomial<C>);
  }
  numerator(): Polynomial<C> {
    return this.num;
  }
  denominator(): Polynomial<C> {
    return this.den;
  }
  /** Scalar hooks from fraction_field_element.pyx:845.
   * @see Deviation: Native Modular Polynomial Products and Fraction Fields
   */
  _conversion<T>(R: { __call__(x: unknown): T }): T {
    if (this.den.eq(1n)) return R.__call__(this.num);
    this.reduce();
    const num = R.__call__(this.num),
      den = R.__call__(this.den);
    if (typeof num === 'bigint' && typeof den === 'bigint') {
      if (den !== 1n && den !== -1n) throw new ArithmeticError('inverse does not exist');
      return (num * den) as T;
    }
    const d = den as unknown as {
      inverse_of_unit?: () => unknown;
      inv: () => unknown;
      modulus?: bigint;
      isUnit?: () => boolean;
    };
    // Rings.ElementMethods.inverse_of_unit checks units before attempting inversion.
    if (typeof d.modulus === 'bigint' && d.isUnit && !d.isUnit())
      throw new ArithmeticError('element is not a unit');
    const inverse = d.inverse_of_unit ? d.inverse_of_unit() : d.inv();
    return (num as unknown as { mul(x: unknown): T }).mul(inverse);
  }
  _integer_(R: IntegerRing): bigint {
    return this._conversion(R);
  }
  _rational_(): Rational {
    // Keep this module independent of runtime QQ/Polynomial imports.
    const convert = (p: Polynomial<C>): Rational =>
      String(p.parent.base_ring) === 'Rational Field'
        ? (p.parent.base_ring.__call__(p) as unknown as Rational)
        : p._rational_();
    if (this.den.eq(1n)) return convert(this.num);
    this.reduce();
    const num = convert(this.num),
      inverse = convert(this.den).inv();
    return num.mul(inverse);
  }
  isZero(): boolean {
    return this.num.isZero();
  }
  eq(other: unknown): boolean {
    try {
      const b = this.parent.__call__(other);
      return this.num.mul(b.denominator()).eq(this.den.mul(b.numerator()));
    } catch (e) {
      if (e instanceof TypeError || e instanceof ValueError) return false;
      throw e;
    }
  }
  toString(): string {
    if (this.isZero()) return '0';
    const n = String(this.num),
      d = String(this.den);
    if (this.den.eq(1n)) return n;
    return `${atomic(this.num) ? n : `(${n})`}/${atomic(this.den) && !d.includes('*') && !d.includes('/') ? d : `(${d})`}`;
  }
}
export class FractionFieldElement_1poly_field<
  C extends RingElement,
> extends FractionFieldElement<C> {
  constructor(
    parent: FractionField_generic<C>,
    n: unknown,
    d: unknown = 1n,
    options: FractionOptions = {}
  ) {
    super(parent, n, d, options);
    if (options.reduce === false) this.normalize_leading_coefficients();
  }
  override reduce(): void {
    if (this._is_reduced) return;
    super.reduce();
    this.normalize_leading_coefficients();
  }
  protected normalize_leading_coefficients(): void {
    const c = inverseCoefficient(this.den.leading_coefficient());
    this.den = this.den.monic();
    this.num = this.num.scalar_mul(c as C);
  }
}
function unit(f: Polynomial<RingElement>): boolean {
  if (f.degree() !== 0) return false;
  if (f.parent.base_ring.is_field?.()) return !f.isZero();
  const c = f.getCoeff(0);
  if (isPolynomial(c)) return unit(c);
  return c.eq(f.parent.base_ring.one()) || c.eq(f.parent.base_ring.one().neg());
}
function inverseCoefficient(c: RingElement): RingElement {
  const inv = (c as unknown as { inv?: () => RingElement }).inv;
  if (inv) return inv.call(c);
  if (isPolynomial(c)) return inverseUnit(c);
  throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: fraction coefficient inverse');
}
function inverseUnit<C extends RingElement>(f: Polynomial<C>): Polynomial<C> {
  const c = f.getCoeff(0),
    one = f.parent.base_ring.one();
  if (c.eq(one)) return f.parent.one();
  if (c.eq(one.neg())) return f.parent.one().neg();
  return f.parent.__call__([inverseCoefficient(c) as C]);
}
function normalizeUnit<C extends RingElement>(
  n: Polynomial<C>,
  d: Polynomial<C>
): [Polynomial<C>, Polynomial<C>] {
  if (!d.eq(1n) && unit(d))
    try {
      return [n.mul(inverseUnit(d)), d.parent.one()];
    } catch (e) {
      if (!(e instanceof AttributeError || e instanceof NotImplementedError)) throw e;
    }
  return [n, d];
}
function atomic(x: RingElement): boolean {
  if (isPolynomial(x)) {
    const coeffs = x.coeffs.filter((c) => !c.isZero());
    return coeffs.length === 1 && atomic(coeffs[0]!);
  }
  if ('value' in x || ('numerator' in x && typeof x.numerator === 'bigint')) return true;
  const s = String(x);
  return !s.includes('+') && !s.includes('-') && !s.includes(' ');
}

/** @internal Polynomial constructor and embedding-section inverse of unit denominators. */
export function _fraction_polynomial<C extends RingElement>(
  x: FractionElement<C>,
  section: boolean
): Polynomial<C> {
  const n = x.numerator(),
    d = x.denominator();
  if (section && x.parent.is_exact() && d.eq(1n)) return n;
  if (!unit(d))
    throw new TypeError(
      section ? 'fraction must have unit denominator' : 'denominator must be a unit'
    );
  return n.mul(inverseUnit(d));
}

function isPolynomial(x: unknown): x is Polynomial<RingElement> {
  return (
    typeof x === 'object' &&
    x !== null &&
    'coeffs' in x &&
    Array.isArray(x.coeffs) &&
    'getCoeff' in x &&
    typeof x.getCoeff === 'function'
  );
}

/** @internal Native FpT partial section; leaves generic fractions to their hooks. */
export function _fraction_native_integer(x: unknown): bigint | undefined {
  if (!isFractionElement(x)) return undefined;
  const ctor = x.constructor as unknown as {
    _section?: (x: FractionElement<RingElement>, constant: boolean) => bigint;
  };
  return ctor._section?.(x, true);
}
