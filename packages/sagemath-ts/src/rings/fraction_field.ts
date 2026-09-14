/** Fraction-field parents; source: sage/rings/fraction_field.py.
 * @see Deviation: Native Modular Polynomial Products and Fraction Fields
 */
import { TypeError, ValueError } from '../errors.js';
import type { PolynomialRing } from './polynomial/polynomial_ring.js';
import type { CoefficientRing, Polynomial, RingElement } from './polynomial/polynomial_element.js';
import type { Rational } from './rational.js';
import {
  FractionFieldElement,
  FractionFieldElement_1poly_field,
  isFractionElement,
  type FractionElement,
  type FractionElementConstructor,
} from './fraction_field_element.js';

export class FractionField_generic<C extends RingElement> {
  private _zero?: FractionElement<C>;
  private _one?: FractionElement<C>;
  constructor(
    readonly polynomial_ring: PolynomialRing<C>,
    readonly _element_class: FractionElementConstructor<C> = FractionFieldElement
  ) {}
  ring(): PolynomialRing<C> {
    return this.polynomial_ring;
  }
  is_field(): boolean {
    return true;
  }
  is_exact(): boolean {
    const exact = (r: unknown): boolean => {
      if (typeof r === 'object' && r !== null && 'variable_name' in r)
        return exact((r as unknown as { base_ring: unknown }).base_ring);
      if (String(r) === 'Real Double Field') return false;
      const f = (r as { is_exact?: () => boolean }).is_exact;
      return typeof f === 'function' ? f.call(r) : true;
    };
    return exact(this.polynomial_ring);
  }
  characteristic(): bigint {
    return _fraction_characteristic(this.polynomial_ring.base_ring);
  }
  base_ring(): CoefficientRing<C> {
    return this.polynomial_ring.base_ring;
  }
  gen(i: unknown = 0): FractionElement<C> {
    return new this._element_class(this, this.polynomial_ring.gen(i), this.polynomial_ring.one(), {
      coerce: false,
      reduce: false,
    });
  }
  ngens(): number {
    return 1;
  }
  __call__(numerator: unknown, denominator?: unknown): FractionElement<C> {
    const hasDen = denominator !== undefined && denominator !== null;
    if (Array.isArray(numerator) && numerator.length === 1) numerator = numerator[0];
    const R = this.ring();
    if (!hasDen && polynomialValue(numerator) && numerator.parent.base_ring === this) {
      if (numerator.degree() > 0) throw new TypeError(`${numerator} is not a constant polynomial`);
      return this.__call__(numerator.getCoeff(0));
    }
    if (!hasDen && isFractionElement(numerator)) {
      if (numerator.parent === this) return numerator as FractionElement<C>;
      const origin = numerator.parent.ring();
      if (R.has_coerce_map_from(origin)) {
        if (
          origin.has_coerce_map_from(R) &&
          _fraction_field_domain(R) &&
          this === R.fraction_field() &&
          numerator.parent === origin.fraction_field()
        )
          return numerator as FractionElement<C>;
        return new this._element_class(this, numerator.numerator(), numerator.denominator());
      }
    }
    if (!hasDen && isRational(numerator))
      return new this._element_class(this, numerator.numerator, numerator.denominator);
    if (hasDen && isFractionElement(numerator) && numerator.parent === this) {
      const y = this.__call__(denominator);
      return new this._element_class(
        this,
        numerator.numerator().mul(y.denominator()),
        y.numerator().mul(numerator.denominator())
      );
    }
    const d = hasDen ? denominator : R.one();
    try {
      return new this._element_class(this, numerator, d);
    } catch (e) {
      if (!(e instanceof TypeError || e instanceof ValueError)) throw e;
      if (
        !(
          typeof numerator === 'string' ||
          typeof d === 'string' ||
          isRational(numerator) ||
          isFractionElement(numerator) ||
          isRational(d) ||
          isFractionElement(d)
        )
      )
        throw e;
    }
    // After polynomial coercion fails, Sage evaluates strings in the fraction
    // generators. Only an unknown name is rewritten to the parent diagnostic.
    const evaluate = (x: unknown): unknown => {
      if (typeof x !== 'string') return x;
      try { return R._parse_fraction_string(x); }
      catch (e) {
        if ((e as Error).name === 'NameError')
          throw new TypeError(`unable to evaluate ${reprString(x)} in ${this}`);
        throw e;
      }
    };
    // fraction_field.py:resolve_fractions keeps arithmetic in the operand parents
    // before constructing in the target. Failed conversions can expose another
    // layer of coefficient denominators.
    let x: unknown = evaluate(numerator),
      y: unknown = hasDen ? evaluate(d) : this.one();
    const parentOf = fractionParent;
    while (true) {
      const x0 = x,
        y0 = y;
      const resolve = (): [unknown, unknown] => {
        const [xn, xd] = fractionParts(x0, R),
          [yn, yd] = fractionParts(y0, R);
        const retry = (f: () => [unknown, unknown]): [unknown, unknown] | null => {
          try {
            return f();
          } catch (e) {
            if (!(e instanceof TypeError || e instanceof ValueError)) throw e;
            return null;
          }
        };
        let pair = retry(() => [fractionMultiply(xn, yd), fractionMultiply(yn, xd)]);
        if (pair) return pair;
        pair = retry(() => {
          const P = parentOf(yd) as CoefficientRing<RingElement>;
          return [fractionMultiply(P.__call__(xn), yd), fractionMultiply(yn, P.__call__(xd))];
        });
        if (pair) return pair;
        pair = retry(() => {
          const P = parentOf(xd) as CoefficientRing<RingElement>;
          return [fractionMultiply(xn, P.__call__(yd)), fractionMultiply(P.__call__(yn), xd)];
        });
        if (pair) return pair;
        throw new TypeError();
      };
      try {
        [x, y] = resolve();
      } catch (e) {
        if (!(e instanceof TypeError)) throw e;
        throw new TypeError(`cannot convert ${x0}/${y0} to an element of ${this}`);
      }
      try {
        return new this._element_class(this, x, y);
      } catch (e) {
        if (!(e instanceof TypeError) || parentOf(x) === parentOf(x0)) throw e;
      }
    }
  }

  zero(): FractionElement<C> {
    return (this._zero ??= this.__call__(0n));
  }
  one(): FractionElement<C> {
    return (this._one ??= this.__call__(1n));
  }
  toString(): string {
    return `Fraction Field of ${this.polynomial_ring}`;
  }
}
export class FractionField_1poly_field<C extends RingElement> extends FractionField_generic<C> {
  constructor(
    R: PolynomialRing<C>,
    element_class: FractionElementConstructor<C> = FractionFieldElement_1poly_field
  ) {
    super(R, element_class);
  }
}
export function FractionField<C extends RingElement>(
  R: PolynomialRing<C>
): FractionField_generic<C> {
  if (R === null || R === undefined || typeof R.__call__ !== 'function')
    throw new TypeError('R must be a ring');
  if (!_fraction_field_domain(R)) throw new TypeError('R must be an integral domain');
  return R.fraction_field();
}
/** @internal Parent/category predicates for the supported polynomial coefficient tower. */
export function _fraction_field_domain(r: unknown): boolean {
  if (typeof r === 'object' && r !== null && 'variable_name' in r)
    return _fraction_field_domain((r as unknown as { base_ring: unknown }).base_ring);
  return (r as { is_field?: () => boolean }).is_field?.() === true || String(r) === 'Integer Ring';
}
export function _fraction_characteristic(r: unknown): bigint {
  if (typeof r === 'object' && r !== null && 'variable_name' in r)
    return _fraction_characteristic((r as unknown as { base_ring: unknown }).base_ring);
  const c = (r as { characteristic?: bigint | number | (() => bigint | number) }).characteristic;
  return c === undefined ? 0n : BigInt(typeof c === 'function' ? c.call(r) : c);
}

function isRational(x: unknown): x is Rational {
  return (
    typeof x === 'object' &&
    x !== null &&
    'numerator' in x &&
    typeof x.numerator === 'bigint' &&
    'denominator' in x &&
    typeof x.denominator === 'bigint'
  );
}

function polynomialValue(x: unknown): x is Polynomial<RingElement> {
  return (
    typeof x === 'object' &&
    x !== null &&
    'coeffs' in x &&
    Array.isArray(x.coeffs) &&
    'getCoeff' in x &&
    typeof x.getCoeff === 'function'
  );
}
function fractionParts<C extends RingElement>(
  x: unknown,
  R: PolynomialRing<C>
): [unknown, unknown] {
  if (isFractionElement(x)) return [x.numerator(), x.denominator()];
  if (isRational(x)) return [x.numerator, x.denominator];
  if (polynomialValue(x)) {
    if (String(R.base_ring) === 'Integer Ring' && String(x.parent.base_ring) === 'Rational Field') {
      // The rational-polynomial numerator/denominator pair is its common
      // integer coefficient vector and positive denominator. Construct the
      // integer vector directly in the eventual target polynomial parent.
      let den = 1n;
      for (const c of x.coeffs as unknown as Rational[]) {
        let a = den,
          b = c.denominator;
        while (b) [a, b] = [b, a % b];
        den = (den / a) * c.denominator;
      }
      return [
        R.__call__(
          (x.coeffs as unknown as Rational[]).map((c) => c.numerator * (den / c.denominator))
        ),
        den,
      ];
    }
    return [x, x.parent.base_ring.one()];
  }
  if (typeof x === 'boolean') return [x ? 1n : 0n, 1n];
  if (typeof x === 'bigint' || typeof x === 'number') return [x, 1n];
  throw new TypeError();
}
function fractionMultiply(a: unknown, b: unknown): unknown {
  if (typeof a === 'bigint' && typeof b === 'bigint') return a * b;
  if (polynomialValue(a)) return a.mul(polynomialValue(b) ? b : a.parent.__call__(b));
  if (polynomialValue(b)) return b.parent.__call__(a).mul(b);
  if (typeof a === 'object' && a !== null && 'mul' in a && typeof a.mul === 'function')
    return a.mul(b);
  if (typeof b === 'object' && b !== null && 'mul' in b && typeof b.mul === 'function')
    return b.mul(a);
  throw new TypeError();
}
function fractionParent(x: unknown): unknown {
  return typeof x === 'object' && x !== null && 'parent' in x ? x.parent : typeof x;
}

function reprString(s: string): string {
  const quote = s.includes("'") && !s.includes('"') ? '"' : "'";
  return quote + s.replace(/\\/g, "\\\\").replaceAll(quote, "\\" + quote)
    .replace(/\n/g, "\\n").replace(/\r/g, "\\r").replace(/\t/g, "\\t") + quote;
}
