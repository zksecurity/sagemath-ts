/**
 * @module sage/schemes/hyperelliptic_curves/hyperelliptic_generic
 * @description Hyperelliptic curves over a general ring
 *
 * Port of: `sage/schemes/hyperelliptic_curves/hyperelliptic_generic.py`
 *
 * A hyperelliptic curve is given by `y^2 + h(x) y = f(x)`.  The plane
 * projective model used by Sage is
 *
 * ```
 * F = y^2 z^(d-2) + F0(x,z) y z^(d-dh-1) - F1(x,z) z^(d-df)
 * ```
 *
 * with `F1 = sum f[i] x^i z^(df-i)`, `F0 = sum h[i] x^i z^(dh-i)`,
 * `d = max(df, dh+1)` (`hyperelliptic_generic.py:88-111`).
 */

import { NotImplementedError, ValueError } from '../../errors.js';
import type {
  Polynomial,
  PolynomialRingBase,
  RingElement,
} from '../../rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
import {
  type HyperellipticBaseRing,
  cardinality_of,
  characteristic_of,
  compare_elements,
  div_elements,
  field_embedding,
  is_square_of,
  sort_roots_like_sage,
  sqrt_all_of,
} from './field_ops.js';

/* -------------------------------------------------------------------------- */
/* Printing helpers                                                           */
/* -------------------------------------------------------------------------- */

/**
 * Reproduce Sage's univariate polynomial `_repr` from a list of coefficient
 * strings (`null` marks a zero coefficient).
 *
 * Port of `sage/rings/polynomial/polynomial_element.pyx:3096-3172`.
 *
 * `atomic_repr` is taken to be `False` throughout: the parenthesisation test
 * is applied to the coefficient string with a leading `-` removed, so for the
 * base rings this module supports (`GF(p)`, `QQ`, `ZZ`) the test never fires
 * and the output is identical to Sage's `atomic_repr = True` branch.
 */
export function sage_poly_repr(coeffStrings: (string | null)[], name: string): string {
  let s = ' ';
  const m = coeffStrings.length;
  for (let n = m - 1; n >= 0; n--) {
    const c = coeffStrings[n];
    if (c === null || c === undefined) {
      continue;
    }
    if (n !== m - 1) {
      s += ' + ';
    }
    let x = c;
    let y = c;
    if (y.startsWith('-')) {
      y = y.slice(1);
    }
    if (n > 0 && (y.includes('+') || y.includes('-'))) {
      x = `(${c})`;
    }
    const varPart = n > 1 ? `*${name}^${n}` : n === 1 ? `*${name}` : '';
    s += x + varPart;
  }
  s = s.split(' + -').join(' - ');
  s = s.split(' 1*').join(' ');
  s = s.split(' -1*').join(' -');
  if (s === ' ') {
    return '0';
  }
  return s.slice(1);
}

/** Render a {@link Polynomial} the way Sage prints it, in the variable `name`. */
export function poly_repr<C extends RingElement>(p: Polynomial<C>, name: string): string {
  const strs: (string | null)[] = [];
  for (let i = 0; i <= p.degree(); i++) {
    const c = p.getCoeff(i);
    strs.push(c.isZero() ? null : c.toString());
  }
  return sage_poly_repr(strs, name);
}

/* -------------------------------------------------------------------------- */
/* Points                                                                     */
/* -------------------------------------------------------------------------- */

/**
 * A point `(x : y : z)` on the plane projective model of a hyperelliptic curve.
 *
 * Sage returns `SchemeMorphism_point_projective_field` instances; those are
 * normalised by dividing through by the **last nonzero** coordinate, which is
 * what {@link HyperellipticPoint.normalize} does.
 */
export class HyperellipticPoint<C extends RingElement> {
  readonly curve: HyperellipticCurve_generic<C>;
  readonly coords: readonly [C, C, C];

  constructor(curve: HyperellipticCurve_generic<C>, coords: [C, C, C]) {
    this.curve = curve;
    this.coords = HyperellipticPoint.normalize(coords);
  }

  /** Divide through by the last nonzero coordinate. */
  static normalize<C extends RingElement>(coords: [C, C, C]): [C, C, C] {
    for (let i = 2; i >= 0; i--) {
      const c = coords[i]!;
      if (!c.isZero()) {
        if (c.eq(1)) {
          return coords;
        }
        return [
          div_elements(coords[0]!, c),
          div_elements(coords[1]!, c),
          div_elements(coords[2]!, c),
        ];
      }
    }
    throw new ValueError('[0, 0, 0] does not define a valid projective point');
  }

  /** `P[n]` in Sage. */
  get(n: number): C {
    const i = n < 0 ? n + 3 : n;
    const c = this.coords[i];
    if (c === undefined) {
      throw new ValueError(`index ${n} out of range for a projective plane point`);
    }
    return c;
  }

  eq(other: HyperellipticPoint<C>): boolean {
    return (
      this.coords[0].eq(other.coords[0]) &&
      this.coords[1].eq(other.coords[1]) &&
      this.coords[2].eq(other.coords[2])
    );
  }

  toString(): string {
    return `(${this.coords[0]} : ${this.coords[1]} : ${this.coords[2]})`;
  }
}

/* -------------------------------------------------------------------------- */
/* The curve                                                                  */
/* -------------------------------------------------------------------------- */

/**
 * Hyperelliptic curve `y^2 + h(x) y = f(x)` over a general ring.
 *
 * Instances are produced by
 * {@link import('./constructor.js').HyperellipticCurve}; the constructor of
 * this class performs no validation, exactly like Sage's
 * `HyperellipticCurve_generic.__init__`.
 */
export class HyperellipticCurve_generic<C extends RingElement> {
  readonly _f: Polynomial<C>;
  readonly _h: Polynomial<C>;
  readonly _genus: number;
  readonly _names: [string, string];

  constructor(
    f: Polynomial<C>,
    h: Polynomial<C>,
    names: [string, string] = ['x', 'y'],
    genus?: number
  ) {
    this._f = f;
    this._h = h;
    this._names = names;
    this._genus = genus ?? genus_of(f, h);
  }

  /** The polynomial ring `K[x]` the defining polynomials live in. */
  get polynomial_ring(): PolynomialRingBase<C> {
    return this._f.parent;
  }

  /** `self.base_ring()` */
  base_ring(): HyperellipticBaseRing<C> {
    return this._f.parent.base_ring as HyperellipticBaseRing<C>;
  }

  /**
   * `self.hyperelliptic_polynomials()`
   *
   * Port of `hyperelliptic_generic.py:196-210`.
   */
  hyperelliptic_polynomials(): [Polynomial<C>, Polynomial<C>] {
    return [this._f, this._h];
  }

  /** `self.genus()` (`hyperelliptic_generic.py:526-527`) */
  genus(): number {
    return this._genus;
  }

  /**
   * Return `False`: hyperelliptic curves are smooth, as checked on
   * construction (`hyperelliptic_generic.py:212-238`).
   */
  is_singular(): boolean {
    return false;
  }

  /** Return `True` (`hyperelliptic_generic.py:240-267`). */
  is_smooth(): boolean {
    return true;
  }

  /**
   * Total degree of the plane projective defining polynomial.
   *
   * `d = max(df, dh + 1)`; when `h = 0` Sage takes the `h is None` branch and
   * the degree is `df` (`hyperelliptic_generic.py:89-98`).
   */
  defining_polynomial_degree(): number {
    const df = this._f.degree();
    if (this._h.isZero()) {
      return df;
    }
    return Math.max(df, this._h.degree() + 1);
  }

  /**
   * Evaluate the plane projective defining polynomial `F(x, y, z)`.
   *
   * Port of `hyperelliptic_generic.py:89-99`.
   */
  defining_polynomial_eval(x: C, y: C, z: C): C {
    const K = this.base_ring();
    const f = this._f;
    const h = this._h;
    const df = f.degree();

    const powC = (base: C, e: number): C => {
      let r = K.one() as C;
      for (let i = 0; i < e; i++) {
        r = r.mul(base) as C;
      }
      return r;
    };

    // F1 = sum f[i] x^i z^(df - i)
    let F1 = K.zero() as C;
    for (let i = 0; i <= df; i++) {
      const c = f.getCoeff(i);
      if (c.isZero()) {
        continue;
      }
      F1 = F1.add(c.mul(powC(x, i)).mul(powC(z, df - i)) as C) as C;
    }

    if (h.isZero()) {
      // F = y^2 z^(df-2) - F1
      return y
        .mul(y)
        .mul(powC(z, df - 2))
        .sub(F1) as C;
    }

    const dh = h.degree();
    const deg = Math.max(df, dh + 1);
    let F0 = K.zero() as C;
    for (let i = 0; i <= dh; i++) {
      const c = h.getCoeff(i);
      if (c.isZero()) {
        continue;
      }
      F0 = F0.add(c.mul(powC(x, i)).mul(powC(z, dh - i)) as C) as C;
    }

    // F = y^2 z^(deg-2) + F0 y z^(deg-dh-1) - F1 z^(deg-df)
    const t1 = y.mul(y).mul(powC(z, deg - 2)) as C;
    const t2 = F0.mul(y).mul(powC(z, deg - dh - 1)) as C;
    const t3 = F1.mul(powC(z, deg - df)) as C;
    return t1.add(t2).sub(t3) as C;
  }

  /**
   * `self.point(coords, check=...)`
   *
   * When `check` is true the coordinates must satisfy the defining polynomial.
   */
  point(coords: [C, C, C], options?: { check?: boolean }): HyperellipticPoint<C> {
    const check = options?.check ?? true;
    if (check) {
      const val = this.defining_polynomial_eval(coords[0], coords[1], coords[2]);
      if (!val.isZero()) {
        throw new ValueError(`Coordinates ${coords} do not define a point on ${this}`);
      }
    }
    return new HyperellipticPoint(this, coords);
  }

  /**
   * Return `True` if `x` is the `x`-coordinate of a point on this curve.
   *
   * Port of `hyperelliptic_generic.py:269-372`.
   */
  is_x_coord(x: C | bigint | number): boolean {
    const K = this.base_ring();
    const [f, h] = this.hyperelliptic_polynomials();

    let xx: C;
    try {
      xx = K.__call__(x) as C;
    } catch {
      throw new TypeError('x must be coercible into the base ring of the curve');
    }

    // When h is zero then x is a valid coordinate if y2 is square
    if (h.isZero()) {
      const y2 = f.evaluate(xx);
      return is_square_of(K, y2);
    }

    const a = f.evaluate(xx);
    const b = h.evaluate(xx);

    if (characteristic_of(K) === 2n) {
      // F = y^2 + b y - a; x is a coordinate iff F has a root
      const R = this.polynomial_ring;
      const F = R.__call__([a.neg() as C, b, K.one() as C]);
      return F.roots().length > 0;
    }

    // D = b^2 + 4a
    const D = b.mul(b).add(int_mul(K, 4n, a)) as C;
    return is_square_of(K, D);
  }

  /**
   * Return one or all points with the given `x`-coordinate.
   *
   * Port of `hyperelliptic_generic.py:374-524`.  The `y`-coordinates are
   * sorted, which makes the choice deterministic exactly as in Sage.
   */
  lift_x(x: C | bigint | number, options?: { all?: false }): HyperellipticPoint<C>;
  lift_x(x: C | bigint | number, options: { all: true }): HyperellipticPoint<C>[];
  lift_x(
    x: C | bigint | number,
    options?: { all?: boolean }
  ): HyperellipticPoint<C> | HyperellipticPoint<C>[] {
    const all = options?.all ?? false;
    const K = this.base_ring();
    const [f, h] = this.hyperelliptic_polynomials();

    let xx: C;
    try {
      xx = K.__call__(x) as C;
    } catch {
      throw new ValueError('x must have a common parent with the base ring');
    }

    let ys: C[];
    const one = K.one() as C;

    if (h.isZero()) {
      const y2 = f.evaluate(xx);
      ys = sqrt_all_of(K, y2);
    } else {
      const a = f.evaluate(xx);
      const b = h.evaluate(xx);
      if (characteristic_of(K) === 2n) {
        const R = this.polynomial_ring;
        const F = R.__call__([a.neg() as C, b, one]);
        ys = F.roots().map(([r]) => r);
        ys.sort((u, v) => compare_elements(u, v));
      } else {
        const D = b.mul(b).add(int_mul(K, 4n, a)) as C;
        ys = sqrt_all_of(K, D).map((d) => div_elements(b.neg().add(d) as C, K.__call__(2n) as C));
        ys.sort((u, v) => compare_elements(u, v));
      }
    }

    if (ys.length > 0) {
      if (all) {
        return ys.map((y) => this.point([xx, y, one], { check: false }));
      }
      return this.point([xx, ys[0]!, one], { check: false });
    }

    if (all) {
      return [];
    }
    throw new ValueError(`No point with x-coordinate ${xx} on ${this}`);
  }

  /**
   * `self.jacobian()` (`hyperelliptic_generic.py:529-531`).
   */
  // biome-ignore lint/suspicious/noExplicitAny: broken by the constructor <-> jacobian import cycle
  jacobian(): any {
    // Imported lazily: `jacobian_generic` imports the curve classes.
    // eslint-disable-next-line @typescript-eslint/no-var-requires
    const mod = getJacobianModule();
    return new mod.HyperellipticJacobian_generic(this);
  }

  /**
   * Return an odd degree model of this curve.
   *
   * Port of `hyperelliptic_generic.py:533-613`.
   */
  odd_degree_model(): HyperellipticCurve_generic<C> {
    const [f, h] = this.hyperelliptic_polynomials();
    if (!h.isZero()) {
      throw new NotImplementedError(
        'odd_degree_model only implemented for curves in Weierstrass form'
      );
    }
    if (f.degree() % 2 !== 0) {
      return this;
    }

    const rts = f.roots().map(([r]) => r);
    if (rts.length === 0) {
      throw new ValueError('No odd degree model exists over field of definition');
    }
    // Sage takes `f.roots(multiplicities=False)[0]`; `sort_roots_like_sage`
    // reproduces the order in which `Polynomial.roots()` lists them.
    const rt = sort_roots_like_sage(rts)[0]!;

    // fnew = numerator of f((x*rt + 1)/x), i.e. sum_i f[i] (rt x + 1)^i x^(d-i)
    const R = this.polynomial_ring;
    const K = this.base_ring();
    const d = f.degree();
    const linear = R.__call__([K.one() as C, rt]); // 1 + rt*x
    let fnew = R.zero();
    for (let i = 0; i <= d; i++) {
      const c = f.getCoeff(i);
      if (c.isZero()) {
        continue;
      }
      fnew = fnew.add(
        linear
          .pow(i)
          .shift(d - i)
          .scalar_mul(c)
      );
    }

    const { HyperellipticCurve } = getConstructorModule();
    return HyperellipticCurve(fnew, R.zero(), { names: this._names });
  }

  /** `hyperelliptic_generic.py:615-635` */
  has_odd_degree_model(): boolean {
    try {
      this.odd_degree_model();
      return true;
    } catch (e) {
      if (e instanceof ValueError) {
        return false;
      }
      throw e;
    }
  }

  /**
   * `self.change_ring(R)` / `self.base_extend(R)`
   *
   * Port of `hyperelliptic_generic.py:113-143`.
   */
  change_ring<D extends RingElement>(R: HyperellipticBaseRing<D>): HyperellipticCurve_generic<D> {
    const { HyperellipticCurve } = getConstructorModule();
    const P = new PolynomialRing<D>(R, this._names[0]);
    // Sage uses the coercion model here; `field_embedding` is the canonical
    // map from the current base ring into `R`.
    const emb = field_embedding<C>(this.base_ring(), R as HyperellipticBaseRing<RingElement>);
    const conv = (p: Polynomial<C>): Polynomial<D> => P.__call__(p.coeffs.map((c) => emb(c) as D));
    return HyperellipticCurve(conv(this._f), conv(this._h), {
      names: this._names,
    }) as HyperellipticCurve_generic<D>;
  }

  base_extend<D extends RingElement>(R: HyperellipticBaseRing<D>): HyperellipticCurve_generic<D> {
    return this.change_ring(R);
  }

  /**
   * `self._magma_init_(magma)` (`hyperelliptic_generic.py:637-660`).
   */
  _magma_init_(): string {
    return `HyperellipticCurve(${poly_repr(this._f, this._names[0])}, ${poly_repr(
      this._h,
      this._names[0]
    )})`;
  }

  /**
   * String representation (`hyperelliptic_generic.py:145-167`).
   */
  toString(): string {
    const [f, h] = this.hyperelliptic_polynomials();
    const R = this.base_ring();
    const [xn, yn] = this._names;
    const y2 = `${yn}^2`;
    if (h.isZero()) {
      return `Hyperelliptic Curve over ${R} defined by ${y2} = ${poly_repr(f, xn)}`;
    }
    // h(x)*y printed as an element of K[x][y]
    const hy = sage_poly_repr([null, poly_repr(h, xn)], yn);
    return `Hyperelliptic Curve over ${R} defined by ${y2} + ${hy} = ${poly_repr(f, xn)}`;
  }

  /** Structural equality, matching Sage's comparison of the defining data. */
  eq(other: HyperellipticCurve_generic<C>): boolean {
    return (
      String(this.base_ring()) === String(other.base_ring()) &&
      this._f.eq(other._f) &&
      this._h.eq(other._h) &&
      this._names[0] === other._names[0] &&
      this._names[1] === other._names[1]
    );
  }

  /* ---------------------------------------------------------------------- */
  /* Not ported                                                             */
  /* ---------------------------------------------------------------------- */

  /** `hyperelliptic_generic.py:662-665` — needs Monsky-Washnitzer cohomology. */
  monsky_washnitzer_gens(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: monsky_washnitzer_gens (sage.schemes.hyperelliptic_curves.monsky_washnitzer is not ported)'
    );
  }

  /** `hyperelliptic_generic.py:667-682` — needs Monsky-Washnitzer cohomology. */
  invariant_differential(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: invariant_differential (sage.schemes.hyperelliptic_curves.monsky_washnitzer is not ported)'
    );
  }

  /** `hyperelliptic_generic.py:684-734` — needs power series rings in this port. */
  local_coordinates_at_nonweierstrass(): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: local_coordinates_at_nonweierstrass');
  }

  /** `hyperelliptic_generic.py:736-787` */
  local_coordinates_at_weierstrass(): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: local_coordinates_at_weierstrass');
  }

  /** `hyperelliptic_generic.py:789-842` */
  local_coordinates_at_infinity(): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: local_coordinates_at_infinity');
  }

  /** `hyperelliptic_generic.py:844-880` */
  local_coord(): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: local_coord');
  }

  /** `hyperelliptic_generic.py:882-935` — needs the generic rational-points machinery. */
  rational_points(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: rational_points (sage.schemes.curves.constructor is not ported)'
    );
  }
}

/* -------------------------------------------------------------------------- */
/* Helpers                                                                    */
/* -------------------------------------------------------------------------- */

/** `n * a`, with `n` a plain integer coerced into the base ring. */
function int_mul<C extends RingElement>(K: HyperellipticBaseRing<C>, n: bigint, a: C): C {
  return (K.__call__(n) as C).mul(a) as C;
}

/**
 * The genus attached to `(f, h)` by Sage's constructor
 * (`constructor.py:287-292`).
 */
export function genus_of<C extends RingElement>(f: Polynomial<C>, h: Polynomial<C>): number {
  const df = f.degree();
  const dh2 = 2 * h.degree();
  return dh2 < df ? floorDiv(df - 1, 2) : floorDiv(dh2 - 1, 2);
}

/** Python's `//` for integers (floor division). */
export function floorDiv(a: number, b: number): number {
  return Math.floor(a / b);
}

/** Whether `K` is a finite field of odd characteristic. */
export function base_is_finite<C extends RingElement>(K: HyperellipticBaseRing<C>): boolean {
  return cardinality_of(K) !== null;
}

/* -------------------------------------------------------------------------- */
/* Lazy module access (breaks the constructor <-> curve import cycle)         */
/* -------------------------------------------------------------------------- */

// biome-ignore lint/suspicious/noExplicitAny: late-bound modules
let _constructorModule: any = null;
// biome-ignore lint/suspicious/noExplicitAny: late-bound modules
let _jacobianModule: any = null;

/** @internal Wire up the constructor module; called by `constructor.ts`. */
// biome-ignore lint/suspicious/noExplicitAny: late-bound modules
export function _register_constructor_module(mod: any): void {
  _constructorModule = mod;
}

/** @internal Wire up the jacobian module; called by `jacobian_generic.ts`. */
// biome-ignore lint/suspicious/noExplicitAny: late-bound modules
export function _register_jacobian_module(mod: any): void {
  _jacobianModule = mod;
}

// biome-ignore lint/suspicious/noExplicitAny: late-bound modules
function getConstructorModule(): any {
  if (_constructorModule === null) {
    throw new NotImplementedError(
      'the hyperelliptic curve constructor module has not been loaded; import it from ./constructor.js'
    );
  }
  return _constructorModule;
}

// biome-ignore lint/suspicious/noExplicitAny: late-bound modules
function getJacobianModule(): any {
  if (_jacobianModule === null) {
    throw new NotImplementedError(
      'the hyperelliptic jacobian module has not been loaded; import it from ./jacobian_generic.js'
    );
  }
  return _jacobianModule;
}
