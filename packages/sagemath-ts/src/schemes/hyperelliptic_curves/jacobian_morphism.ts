/**
 * @module sage/schemes/hyperelliptic_curves/jacobian_morphism
 * @description Jacobian 'morphism' as a class in the Picard group
 *
 * Port of: `sage/schemes/hyperelliptic_curves/jacobian_morphism.py`
 *
 * This module implements the group operation in the Picard group of a
 * hyperelliptic curve, represented as divisors in Mumford representation,
 * using Cantor's algorithm.
 *
 * A divisor on the hyperelliptic curve `y^2 + y h(x) = f(x)` is stored in
 * Mumford representation, that is, as two polynomials `u(x)` and `v(x)` such
 * that
 *
 * - `u(x)` is monic,
 * - `u(x)` divides `f(x) - h(x) v(x) - v(x)^2`,
 * - `deg(v(x)) < deg(u(x)) <= g`.
 *
 * REFERENCES:
 *
 * - J. Scholten, F. Vercauteren. *An Introduction to Elliptic and Hyperelliptic
 *   Curve Cryptography and the NTRU Cryptosystem.*
 * - R. Avanzi, H. Cohen, C. Doche, G. Frey, T. Lange, K. Nguyen, and
 *   F. Vercauteren. *Handbook of Elliptic and Hyperelliptic Curve
 *   Cryptography.* CRC Press, 2005.
 */

import { ValueError } from '../../errors.js';
import type { Polynomial, RingElement } from '../../rings/polynomial/polynomial_element.js';
import { type IntegerLike, toBigInt } from '../../types/coercion.js';
import { sort_roots_like_sage } from './field_ops.js';
import { poly_repr, sage_poly_repr } from './hyperelliptic_generic.js';
import type { JacobianHomset_divisor_classes } from './jacobian_homset.js';

/** A Mumford divisor `(a(x), b(x))`. */
export type MumfordPair<C extends RingElement> = [Polynomial<C>, Polynomial<C>];

/**
 * Return the unique reduced divisor linearly equivalent to `(a, b)` on the
 * curve `y^2 = f(x)`.
 *
 * Port of `jacobian_morphism.py:122-158`.
 */
export function cantor_reduction_simple<C extends RingElement>(
  a: Polynomial<C>,
  b: Polynomial<C>,
  f: Polynomial<C>,
  genus: number
): MumfordPair<C> {
  let a2 = f.sub(b.mul(b)).quo_rem(a)[0];
  a2 = a2.monic();
  const b2 = b.neg().mod(a2);
  if (a2.degree() === a.degree()) {
    // Upstream prints "Returning ambiguous form of degree genus+1." and
    // returns (a2, b2) (`jacobian_morphism.py:151-155`).
    if (a2.degree() !== genus + 1) {
      throw new ValueError(
        `cantor_reduction_simple: expected degree ${genus + 1}, got ${a2.degree()}`
      );
    }
    return [a2, b2];
  }
  if (a2.degree() > genus) {
    return cantor_reduction_simple(a2, b2, f, genus);
  }
  return [a2, b2];
}

/**
 * Return the unique reduced divisor linearly equivalent to `(a, b)` on the
 * curve `y^2 + y h(x) = f(x)`.
 *
 * Port of `jacobian_morphism.py:161-211`.
 */
export function cantor_reduction<C extends RingElement>(
  a0: Polynomial<C>,
  b0: Polynomial<C>,
  f: Polynomial<C>,
  h: Polynomial<C>,
  genus: number
): MumfordPair<C> {
  let a = a0;
  let b = b0;
  if (a.degree() >= 2 * genus + 1) {
    throw new ValueError('cantor_reduction: deg(a) must be < 2*genus + 1');
  }
  if (b.degree() >= a.degree()) {
    throw new ValueError('cantor_reduction: deg(b) must be < deg(a)');
  }
  let k = f.sub(h.mul(b)).sub(b.mul(b));
  if (2 * a.degree() === k.degree()) {
    // must adjust b to include the point at infinity
    const g1 = a.degree();
    const R = a.parent;
    const K = R.base_ring;
    const x = R.gen();
    // r is a root of x^2 + h[g1] x - f[2 g1]
    const quad = R.__call__([f.getCoeff(2 * g1).neg() as C, h.getCoeff(g1), K.one() as C]);
    const roots = quad.roots().map(([r]) => r);
    if (roots.length === 0) {
      throw new ValueError(
        `cantor_reduction: x^2 + ${h.getCoeff(g1)}*x - ${f.getCoeff(
          2 * g1
        )} has no root in the base field`
      );
    }
    // Sage takes `.roots()[0][0]`; `sort_roots_like_sage` reproduces the order
    // in which `Polynomial.roots()` lists them.
    const r = sort_roots_like_sage(roots)[0]!;
    const xg1 = x.pow(g1);
    b = b.add(xg1.sub(xg1.mod(a)).scalar_mul(r));
    k = f.sub(h.mul(b)).sub(b.mul(b));
  }
  const [q, rem] = k.quo_rem(a);
  if (!rem.isZero()) {
    throw new ValueError('cantor_reduction: a does not divide f - h*b - b^2');
  }
  a = q.monic();
  b = b.add(h).neg().mod(a);
  if (a.degree() > genus) {
    return cantor_reduction(a, b, f, h, genus);
  }
  return [a, b];
}

/**
 * Given two reduced Mumford divisors on the Jacobian of `y^2 = f(x)`, compute
 * a (not necessarily reduced) representative of `D1 + D2`.
 *
 * Port of `jacobian_morphism.py:214-270`.
 */
export function cantor_composition_simple<C extends RingElement>(
  D1: MumfordPair<C>,
  D2: MumfordPair<C>,
  f: Polynomial<C>,
  _genus: number
): MumfordPair<C> {
  const [a1, b1] = D1;
  const [a2, b2] = D2;
  let a: Polynomial<C>;
  let b: Polynomial<C>;

  if (a1.eq(a2) && b1.eq(b2)) {
    // Duplication law
    const [d, _h1, h3] = a1.xgcd(b1.add(b1));
    const q = a1.quo_rem(d)[0];
    a = q.mul(q);
    b = b1.add(h3.mul(f.sub(b1.mul(b1)).quo_rem(d)[0])).mod(a);
  } else {
    const [d0, _u, h2] = a1.xgcd(a2);
    if (d0.degree() === 0 && d0.getCoeff(0).eq(1)) {
      a = a1.mul(a2);
      b = b2.add(h2.mul(a2).mul(b1.sub(b2))).mod(a);
    } else {
      const [d, l, h3] = d0.xgcd(b1.add(b2));
      a = a1.mul(a2).quo_rem(d.mul(d))[0];
      b = b2
        .add(l.mul(h2).mul(b1.sub(b2)).mul(a2.quo_rem(d)[0]))
        .add(h3.mul(f.sub(b2.mul(b2)).quo_rem(d)[0]))
        .mod(a);
    }
  }
  a = a.monic();
  return [a, b];
}

/**
 * Cantor composition on `y^2 + y h(x) = f(x)`.
 *
 * Port of `jacobian_morphism.py:273-352`.
 */
export function cantor_composition<C extends RingElement>(
  D1: MumfordPair<C>,
  D2: MumfordPair<C>,
  f: Polynomial<C>,
  h: Polynomial<C>,
  _genus: number
): MumfordPair<C> {
  const [a1, b1] = D1;
  const [a2, b2] = D2;
  let a: Polynomial<C>;
  let b: Polynomial<C>;

  if (a1.eq(a2) && b1.eq(b2)) {
    // Duplication law
    const [d, _h1, h3] = a1.xgcd(b1.add(b1).add(h));
    const q = a1.quo_rem(d)[0];
    a = q.mul(q);
    b = b1.add(h3.mul(f.sub(h.mul(b1)).sub(b1.mul(b1)).quo_rem(d)[0])).mod(a);
  } else {
    const [d0, _u, h2] = a1.xgcd(a2);
    if (d0.degree() === 0 && d0.getCoeff(0).eq(1)) {
      a = a1.mul(a2);
      b = b2.add(h2.mul(a2).mul(b1.sub(b2))).mod(a);
    } else {
      const e0 = b1.add(b2).add(h);
      if (e0.isZero()) {
        a = a1.mul(a2).quo_rem(d0.mul(d0))[0];
        b = b2.add(h2.mul(b1.sub(b2)).mul(a2.quo_rem(d0)[0])).mod(a);
      } else {
        const [d, l, h3] = d0.xgcd(e0);
        a = a1.mul(a2).quo_rem(d.mul(d))[0];
        b = b2
          .add(l.mul(h2).mul(b1.sub(b2)).mul(a2.quo_rem(d)[0]))
          .add(h3.mul(f.sub(h.mul(b2)).sub(b2.mul(b2)).quo_rem(d)[0]))
          .mod(a);
      }
    }
  }
  a = a.monic();
  return [a, b];
}

/**
 * An element of a Jacobian defined over a field, i.e. in
 * `J(K) = Pic^0_K(C)`.
 *
 * Port of `jacobian_morphism.py:355-877`.
 */
export class JacobianMorphism_divisor_class_field<C extends RingElement> {
  private readonly __polys: MumfordPair<C>;
  readonly parent: JacobianHomset_divisor_classes<C>;

  /**
   * Create a new Jacobian element in Mumford representation.
   *
   * .. warning:: Not for external use!  Use `J(K)([u, v])` instead.
   *
   * Port of `jacobian_morphism.py:360-406`.
   */
  constructor(
    parent: JacobianHomset_divisor_classes<C>,
    polys: MumfordPair<C>,
    options?: { check?: boolean }
  ) {
    const check = options?.check ?? true;
    this.parent = parent;
    let p = polys;
    if (check) {
      const C_ = parent.curve();
      const [f, h] = C_.hyperelliptic_polynomials();
      const [a, b] = polys;
      if (!b.mul(b).add(h.mul(b)).sub(f).mod(a).isZero()) {
        throw new ValueError(
          `Argument polys (= ${format_mumford_pair(polys, parent)}) must be divisor on curve ${C_}.`
        );
      }
      const genus = C_.genus();
      if (a.degree() > genus) {
        p = cantor_reduction(a, b, f, h, genus);
      }
    }
    this.__polys = p;
  }

  /** The pair `(a(x), b(x))` (`jacobian_morphism.py:590-612`). */
  list(): MumfordPair<C> {
    return [this.__polys[0], this.__polys[1]];
  }

  /** `P[n]` (`jacobian_morphism.py:638-667`). */
  get(n: number): Polynomial<C> {
    const i = n < 0 ? n + 2 : n;
    const p = this.__polys[i as 0 | 1];
    if (p === undefined) {
      throw new ValueError(`index ${n} out of range for a Mumford divisor`);
    }
    return p;
  }

  /** `self.scheme()` (`jacobian_morphism.py:489-516`). */
  scheme(): unknown {
    return this.parent.codomain();
  }

  /**
   * Return `true` if this divisor is *not* the additive identity
   * (`jacobian_morphism.py:730-749`: `self.__polys[0] != 1`).
   */
  isNonzero(): boolean {
    const a = this.__polys[0];
    return !(a.degree() === 0 && a.getCoeff(0).eq(1));
  }

  is_zero(): boolean {
    return !this.isNonzero();
  }

  /** `jacobian_morphism.py:669-725` */
  eq(other: JacobianMorphism_divisor_class_field<C>): boolean {
    // Sage reaches `_richcmp_` (which compares `self.scheme()`) only after the
    // coercion model has found a common parent; two divisors living in
    // *different* homset objects therefore compare unequal even when the two
    // homsets describe the same Jacobian.  Verified against SageMath:
    //
    //     sage: J1 = H.jacobian()(GF(37)); J2 = H.jacobian()(GF(37))
    //     sage: J1(H.lift_x(2)) == J2(H.lift_x(2))
    //     False
    if (this.parent.codomain() !== other.parent.codomain()) {
      return false;
    }
    return this.__polys[0].eq(other.__polys[0]) && this.__polys[1].eq(other.__polys[1]);
  }

  /**
   * The additive inverse of this divisor (`jacobian_morphism.py:751-812`).
   */
  neg(): JacobianMorphism_divisor_class_field<C> {
    if (this.is_zero()) {
      return this;
    }
    const polys = this.__polys;
    const X = this.parent;
    const [_f, h] = X.curve().hyperelliptic_polynomials();
    let D: MumfordPair<C>;
    if (h.isZero()) {
      D = [polys[0], polys[1].neg()];
    } else {
      D = [polys[0], polys[1].neg().sub(h.add(polys[0]).mod(polys[0]))];
    }
    return new JacobianMorphism_divisor_class_field(X, D, { check: false });
  }

  /**
   * A Mumford representative of `self + other` (`jacobian_morphism.py:814-843`).
   */
  add(other: JacobianMorphism_divisor_class_field<C>): JacobianMorphism_divisor_class_field<C> {
    const X = this.parent;
    const Cv = X.curve();
    const [f, h] = Cv.hyperelliptic_polynomials();
    const genus = Cv.genus();
    let D: MumfordPair<C>;
    if (h.isZero()) {
      D = cantor_composition_simple(this.__polys, other.__polys, f, genus);
      if (D[0].degree() > genus) {
        D = cantor_reduction_simple(D[0], D[1], f, genus);
      }
    } else {
      D = cantor_composition(this.__polys, other.__polys, f, h, genus);
      if (D[0].degree() > genus) {
        D = cantor_reduction(D[0], D[1], f, h, genus);
      }
    }
    return new JacobianMorphism_divisor_class_field(X, D, { check: false });
  }

  /** `jacobian_morphism.py:845-877` */
  sub(other: JacobianMorphism_divisor_class_field<C>): JacobianMorphism_divisor_class_field<C> {
    return this.add(other.neg());
  }

  /**
   * `n * self`, by double-and-add.
   *
   * Sage inherits this from `AdditiveGroupElement`; the group law itself is
   * `_add_` above.
   */
  mul(n: IntegerLike): JacobianMorphism_divisor_class_field<C> {
    let m = toBigInt(n);
    let negate = false;
    if (m < 0n) {
      negate = true;
      m = -m;
    }
    let result = this.parent.zero();
    let base: JacobianMorphism_divisor_class_field<C> = this;
    while (m > 0n) {
      if (m & 1n) {
        result = result.add(base);
      }
      base = base.add(base);
      m >>= 1n;
    }
    return negate ? result.neg() : result;
  }

  /** Format the Mumford polynomials for printing (`jacobian_morphism.py:408-430`). */
  _printing_polys(): [string, string] {
    const [xn, yn] = this.parent._printing_names();
    const [a, b] = this.__polys;
    const aStr = poly_repr(a, xn);
    // y - b(x), as an element of K[x][y]
    const minusB = b.neg();
    const bStr = minusB.isZero() ? null : poly_repr(minusB, xn);
    const yStr = sage_poly_repr([bStr, '1'], yn);
    return [aStr, yStr];
  }

  /** `jacobian_morphism.py:432-457` */
  toString(): string {
    if (this.is_zero()) {
      return '(1)';
    }
    const [a, b] = this._printing_polys();
    return `(${a}, ${b})`;
  }
}

function format_mumford_pair<C extends RingElement>(
  polys: MumfordPair<C>,
  parent: JacobianHomset_divisor_classes<C>
): string {
  const [xn] = parent._printing_names();
  return `(${poly_repr(polys[0], xn)}, ${poly_repr(polys[1], xn)})`;
}
