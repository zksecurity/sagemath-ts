/**
 * Numeric coefficient subset of sage/rings/real_double.pyx.
 * @see Deviation: Polynomial Resultant Delegation and Real Double Coefficients
 */
import { NotImplementedError, ValueError, ZeroDivisionError } from '../errors.js';
import { gsl_pow_int, gsl_sf_exp, gsl_sf_log } from '@sagemath-ts/gsl-ts';
import { Integer } from './integer_ring.js';
import { type IntegerLike, toBigInt } from '../types/coercion.js';
import { Rational } from './rational.js';
import type { RingElement } from './polynomial/polynomial_element.js';

type RealDoubleOperand = RealDoubleElement | IntegerLike | Rational | number | boolean;

/** Sage coerces scalar operands to RDF before applying double arithmetic. */
function doubleOperand(parent: RealDoubleField_class, value: RealDoubleOperand): number {
  return value instanceof RealDoubleElement ? value.value : parent.__call__(value).value;
}

export class RealDoubleField_class {
  private _zero?: RealDoubleElement;
  private _one?: RealDoubleElement;
  __call__(x: unknown): RealDoubleElement {
    if (x instanceof RealDoubleElement && x.parent === this) return x;
    let value: number;
    if (x instanceof Rational) value = x.toNumber();
    else if (x instanceof Integer || typeof x === 'bigint') {
      value = Number(x instanceof Integer ? x.value : x);
    } else if (typeof x === 'number' || typeof x === 'boolean') value = Number(x);
    else throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: RDF conversion from this input type');
    return new RealDoubleElement(value, this);
  }
  zero(): RealDoubleElement {
    return (this._zero ??= new RealDoubleElement(0, this));
  }
  one(): RealDoubleElement {
    return (this._one ??= new RealDoubleElement(1, this));
  }
  is_field(): boolean {
    return true;
  }
  characteristic(): bigint {
    return 0n;
  }
  toString(): string {
    return 'Real Double Field';
  }
}
export const RDF = new RealDoubleField_class();
export class RealDoubleElement implements RingElement {
  readonly value: number;
  readonly parent: RealDoubleField_class;
  constructor(value: number, parent = RDF) {
    this.value = value;
    this.parent = parent;
  }
  add(x: RealDoubleOperand): this;
  add(x: this): this;
  add(x: RealDoubleOperand): this {
    // ModuleElement._add_long returns self for Python integer/bool zero.
    if (x === 0n || x === false) return this;
    return new RealDoubleElement(this.value + doubleOperand(this.parent, x), this.parent) as this;
  }
  sub(x: RealDoubleOperand): this;
  sub(x: this): this;
  sub(x: RealDoubleOperand): this {
    return new RealDoubleElement(this.value - doubleOperand(this.parent, x), this.parent) as this;
  }
  mul(x: RealDoubleOperand): this;
  mul(x: this): this;
  mul(x: RealDoubleOperand): this {
    // ModuleElement._mul_long returns self for Python integer/bool one.
    if (x === 1n || x === true) return this;
    return new RealDoubleElement(this.value * doubleOperand(this.parent, x), this.parent) as this;
  }
  div(x: RealDoubleOperand): this;
  div(x: this): this;
  div(x: RealDoubleOperand): this {
    return new RealDoubleElement(this.value / doubleOperand(this.parent, x), this.parent) as this;
  }
  /** Integer exponent operator: Sage real_double_element_gsl._pow_long/_pow_int.
   * @see Deviation: Polynomial Modular Powers
   */
  pow(exponentInput: IntegerLike): this {
    const n = toBigInt(exponentInput);
    if (n >= -2048n && n <= 2048n)
      return new RealDoubleElement(gsl_pow_int(this.value, Number(n)), this.parent) as this;
    let v = this.value,
      sign = n & 1n ? -1 : 1;
    const exponent = Number(n);
    if (v >= 0) {
      if (v === 1) return this;
      if (exponent === 0) return new RealDoubleElement(1, this.parent) as this;
      if (v === 0) {
        if (exponent < 0) throw new ZeroDivisionError('0.0 cannot be raised to a negative power');
        return this;
      }
      sign = 1;
    } else {
      const parity = exponent % 2;
      if (parity === 0) {
      } else if (parity === 1) sign = -1;
      else throw new ValueError('negative number cannot be raised to a fractional power');
      v = -v;
    }
    return new RealDoubleElement(sign * gsl_sf_exp(gsl_sf_log(v) * exponent), this.parent) as this;
  }
  /** Sage real_double.pyx __invert__, including signed infinities. */
  inv(): this {
    return new RealDoubleElement(1 / this.value, this.parent) as this;
  }
  neg(): this {
    return new RealDoubleElement(-this.value, this.parent) as this;
  }
  eq(x: RealDoubleOperand): boolean;
  eq(x: this | number): boolean;
  eq(x: RealDoubleOperand): boolean {
    return this.value === doubleOperand(this.parent, x);
  }
  isZero(): boolean {
    return this.value === 0;
  }
  toString(): string {
    const x = this.value;
    if (Number.isNaN(x)) return 'NaN';
    if (!Number.isFinite(x)) return x < 0 ? '-infinity' : '+infinity';
    if (x === 0) return Object.is(x, -0) ? '-0.0' : '0.0';
    if (Math.abs(x) >= 1e16 || Math.abs(x) < 1e-4) {
      const [m, e] = x.toExponential().split('e');
      const n = Number(e);
      return `${m}e${n < 0 ? '-' : '+'}${String(Math.abs(n)).padStart(2, '0')}`;
    }
    const text = String(x);
    return text.includes('.') ? text : text + '.0';
  }
}
