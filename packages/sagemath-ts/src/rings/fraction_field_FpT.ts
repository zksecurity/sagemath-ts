/** Native modular polynomial fractions; source: sage/rings/fraction_field_FpT.pyx.
 * @see Deviation: Native Modular Polynomial Products and Fraction Fields
 */
import {
  _nmod_poly_sqrt,
  n_jacobi,
  _nmod_poly_add,
  _nmod_poly_sub,
  _nmod_poly_mul,
  _nmod_poly_pow,
  _nmod_poly_gcd,
  _nmod_poly_divrem,
  _nmod_poly_make_monic,
} from '@sagemath-ts/flint-ts';
import {
  NotImplementedError,
  AssertionError,
  TypeError,
  ValueError,
  ZeroDivisionError,
  OverflowError,
} from '../errors.js';
import { FractionField_1poly_field, type FractionField_generic } from './fraction_field.js';
import type { FractionElement, FractionOptions } from './fraction_field_element.js';
import { Polynomial, type RingElement } from './polynomial/polynomial_element.js';
import { PrimeFieldElement } from './finite_rings/finite_field_extension.js';
import { FiniteFieldElement as LegacyPrimeElement } from './finite_rings/finite_field_prime.js';
import { IntegerMod } from './finite_rings/integer_mod.js';
import type { PolynomialRing } from './polynomial/polynomial_ring.js';
export class FpT<C extends RingElement> extends FractionField_1poly_field<C> {
  readonly p: bigint;
  readonly poly_ring: PolynomialRing<C>;
  constructor(R: PolynomialRing<C>) {
    super(R, FpTElement);
    this.p = this.characteristic();
    this.poly_ring = R;
    if (this.p >= 1n << 63n) throw new OverflowError('Python int too large to convert to C long');
    if (this.p <= 2n || this.p >= 46341n) throw new AssertionError('');
    const c = R.base_ring.zero();
    if (
      !(
        c instanceof PrimeFieldElement ||
        c instanceof LegacyPrimeElement ||
        c instanceof IntegerMod
      )
    )
      throw new TypeError('unsupported polynomial ring');
  }
}
export class FpTElement<C extends RingElement> implements FractionElement<C> {
  /** @internal Shared native polynomial/constant embedding-section kernel. */
  static _section<C extends RingElement>(
    x: FpTElement<C>,
    constant: boolean
  ): Polynomial<C> | bigint {
    if (x.d.length !== 1 || (constant && x.n.length > 1)) {
      x.normalize();
      if (x.d.length !== 1) throw new ValueError('not integral');
      if (constant && x.n.length > 1) throw new ValueError('not constant');
    }
    if (x.d[0] !== 1n) x.normalize();
    return constant ? (x.n[0] ?? 0n) : x.numerator();
  }

  readonly _is_fraction_field_element = true as const;
  readonly parent: FpT<C>;
  private n: bigint[];
  private d: bigint[];
  constructor(
    parent: FractionField_generic<C>,
    n: unknown,
    d: unknown = 1n,
    options: FractionOptions = {}
  ) {
    this.parent = parent as FpT<C>;
    const R = parent.ring();
    const a = options.coerce === false ? (n as Polynomial<C>) : R.__call__(n),
      b = options.coerce === false ? (d as Polynomial<C>) : R.__call__(d);
    this.n = a.coeffs.map((c) => BigInt(String(c)));
    this.d = b.coeffs.map((c) => BigInt(String(c)));
    if (!this.d.length) throw new ZeroDivisionError('fraction has denominator 0');
    if (options.reduce !== false) this.normalize();
  }
  private create(n: bigint[], d: bigint[], reduce = false): FpTElement<C> {
    const x = Object.create(FpTElement.prototype) as FpTElement<C>;
    Object.defineProperty(x, 'parent', { value: this.parent });
    Object.defineProperty(x, '_is_fraction_field_element', { value: true });
    x.n = n;
    x.d = d;
    if (reduce) x.normalize();
    return x;
  }
  private scale(): void {
    if (this.d.at(-1) === 1n) return;
    const c = this.parent.base_ring().__call__(this.d.at(-1)!) as C & { inv(): C };
    const inv = BigInt(String(c.inv()));
    this.n = _nmod_poly_mul(this.n, [inv], this.parent.p);
    this.d = _nmod_poly_mul(this.d, [inv], this.parent.p);
  }
  private normalize(): void {
    if (!this.n.length) {
      this.d = [1n];
      return;
    }
    const p = this.parent.p;
    if (this.n.length > 1 && this.d.length > 1) {
      const g = _nmod_poly_make_monic(_nmod_poly_gcd(this.n, this.d, p), p);
      if (g.length !== 1) {
        this.n = _nmod_poly_divrem(this.n, g, p)[0];
        this.d = _nmod_poly_divrem(this.d, g, p)[0];
      }
    }
    this.scale();
  }
  add(other: unknown): FpTElement<C> {
    const b = this.parent.__call__(other) as FpTElement<C>,
      p = this.parent.p;
    return this.create(
      _nmod_poly_add(_nmod_poly_mul(this.n, b.d, p), _nmod_poly_mul(b.n, this.d, p), p),
      _nmod_poly_mul(this.d, b.d, p),
      true
    );
  }
  sub(other: unknown): FpTElement<C> {
    const b = this.parent.__call__(other) as FpTElement<C>,
      p = this.parent.p;
    return this.create(
      _nmod_poly_sub(_nmod_poly_mul(this.n, b.d, p), _nmod_poly_mul(b.n, this.d, p), p),
      _nmod_poly_mul(this.d, b.d, p),
      true
    );
  }
  mul(other: unknown): FpTElement<C> {
    const b = this.parent.__call__(other) as FpTElement<C>,
      p = this.parent.p;
    return this.create(_nmod_poly_mul(this.n, b.n, p), _nmod_poly_mul(this.d, b.d, p), true);
  }
  div(other: unknown): FpTElement<C> {
    const b = this.parent.__call__(other) as FpTElement<C>,
      p = this.parent.p;
    if (!b.n.length) throw new ZeroDivisionError('');
    return this.create(_nmod_poly_mul(this.n, b.d, p), _nmod_poly_mul(this.d, b.n, p), true);
  }
  neg(): FpTElement<C> {
    return this.create(_nmod_poly_sub([], this.n, this.parent.p), this.d.slice());
  }
  inv(): FpTElement<C> {
    if (!this.n.length) throw new ZeroDivisionError('');
    return this.create(this.d.slice(), this.n.slice());
  }
  pow(exponent: bigint | number): FpTElement<C> {
    const e = BigInt(exponent);
    if (e < -(1n << 63n) || e >= 1n << 63n)
      throw new OverflowError('Python int too large to convert to C ssize_t');
    if (e < 0n && !this.n.length) throw new ZeroDivisionError('');
    const a = e < 0n ? -e : e,
      p = this.parent.p;
    const n = _nmod_poly_pow(e < 0n ? this.d : this.n, a, p),
      d = _nmod_poly_pow(e < 0n ? this.n : this.d, a, p),
      x = this.create(n, d);
    if (e < 0n) x.scale();
    return x;
  }
  /** Native FpT square root, preserving raw zero and the smaller numerator root. */
  _sqrt_or_None(): FpTElement<C> | null {
    if (!this.n.length) return this;
    const p = this.parent.p;
    const check = (f: bigint[]) => (f.length & 1) === 1 &&
      n_jacobi(f.at(-1)!, p) === 1 && n_jacobi(f[0]!, p) !== -1;
    if (!check(this.n) || !check(this.d)) return null;
    const n = _nmod_poly_sqrt(this.n, p), d = _nmod_poly_sqrt(this.d, p);
    if (n === null || d === null) return null;
    const s = this.create(n, d);
    s.scale();
    const a = s.n.at(-1)!;
    if (a > p - a) s.n = _nmod_poly_sub([], s.n, p);
    return s;
  }
  is_square(): boolean { return this._sqrt_or_None() !== null; }
  sqrt(extend = true, all = false): FpTElement<C> | FpTElement<C>[] {
    const s = this._sqrt_or_None();
    if (s === null) {
      if (extend) throw new NotImplementedError('function fields not yet implemented');
      throw new ValueError('not a perfect square');
    }
    return all ? (s.isZero() ? [s] : [s, s.neg()]) : s;
  }
  numerator(): Polynomial<C> {
    return new Polynomial(
      this.n.map((c) => this.parent.base_ring().__call__(c)),
      this.parent.ring()
    );
  }
  denominator(): Polynomial<C> {
    return new Polynomial(
      this.d.map((c) => this.parent.base_ring().__call__(c)),
      this.parent.ring()
    );
  }
  numer(): Polynomial<C> {
    return this.numerator();
  }
  denom(): Polynomial<C> {
    return this.denominator();
  }
  eq(other: unknown): boolean {
    try {
      const b = this.parent.__call__(other) as FpTElement<C>;
      return (
        this.n.length === b.n.length &&
        this.d.length === b.d.length &&
        this.n.every((v, i) => v === b.n[i]) &&
        this.d.every((v, i) => v === b.d[i])
      );
    } catch (e) {
      if (e instanceof TypeError || e instanceof ValueError) return false;
      throw e;
    }
  }
  isZero(): boolean {
    return this.n.length === 0 && this.d.length === 1 && this.d[0] === 1n;
  }
  toString(): string {
    const n = String(this.numerator());
    if (this.d.length === 1 && this.d[0] === 1n) return n;
    const d = String(this.denominator()),
      wrap = (s: string) => (/[+-]/.test(s) ? `(${s})` : s);
    return `${wrap(n)}/${wrap(d)}`;
  }
}
