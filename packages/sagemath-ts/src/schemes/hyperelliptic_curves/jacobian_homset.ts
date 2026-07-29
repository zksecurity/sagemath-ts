/**
 * @module sage/schemes/hyperelliptic_curves/jacobian_homset
 * @description Rational point sets on a Jacobian
 *
 * Port of: `sage/schemes/hyperelliptic_curves/jacobian_homset.py`
 */

import { ValueError } from '../../errors.js';
import type { Polynomial, RingElement } from '../../rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
import type { HyperellipticBaseRing } from './field_ops.js';
import type { HyperellipticCurve_generic, HyperellipticPoint } from './hyperelliptic_generic.js';
import type { HyperellipticJacobian_generic } from './jacobian_generic.js';
import { JacobianMorphism_divisor_class_field, type MumfordPair } from './jacobian_morphism.js';

/** Anything `JacobianHomset.__call__` accepts. */
export type JacobianPointInput<C extends RingElement> =
  | 0
  | 0n
  | bigint
  | [Polynomial<C>, Polynomial<C>]
  | [bigint, Polynomial<C>]
  | [Polynomial<C>, bigint]
  | [bigint, bigint]
  | [0]
  | HyperellipticPoint<C>
  | [HyperellipticPoint<C>, HyperellipticPoint<C>]
  | JacobianMorphism_divisor_class_field<C>;

/**
 * The set of `S`-rational points of the Jacobian of a hyperelliptic curve.
 *
 * Port of `jacobian_homset.py:56-185`.
 */
export class JacobianHomset_divisor_classes<C extends RingElement> {
  private readonly _jacobian: HyperellipticJacobian_generic<C>;
  private readonly _value_ring: HyperellipticBaseRing<C>;
  private readonly _poly_ring: PolynomialRing<C>;

  constructor(jacobian: HyperellipticJacobian_generic<C>, S: HyperellipticBaseRing<C>) {
    this._jacobian = jacobian;
    this._value_ring = S;
    this._poly_ring = new PolynomialRing<C>(S, jacobian.curve()._names[0]);
  }

  /** `jacobian_homset.py:168-169` */
  curve(): HyperellipticCurve_generic<C> {
    return this._jacobian.curve();
  }

  /** The Jacobian itself, i.e. the codomain of these morphisms. */
  codomain(): HyperellipticJacobian_generic<C> {
    return this._jacobian;
  }

  /** `jacobian_homset.py:171-180` — return `S` for a homset `X(T)`, `T = Spec(S)`. */
  value_ring(): HyperellipticBaseRing<C> {
    return this._value_ring;
  }

  /** The polynomial ring `S[x]` Mumford polynomials live in. */
  polynomial_ring(): PolynomialRing<C> {
    return this._poly_ring;
  }

  /** @internal The `(x, y)` variable names used when printing divisors. */
  _printing_names(): [string, string] {
    return this.curve()._names;
  }

  /** The identity element `(1)` of the group. */
  zero(): JacobianMorphism_divisor_class_field<C> {
    const R = this._poly_ring;
    return new JacobianMorphism_divisor_class_field<C>(this, [R.one(), R.zero()], {
      check: false,
    });
  }

  /** `jacobian_homset.py:182-185` */
  base_extend(_R: unknown): never {
    throw new ValueError(
      'Jacobian point sets viewed as modules over rings other than ZZ not implemented'
    );
  }

  /**
   * Return a rational point in `J(K)`, given
   *
   * 0. a point `P` in `J = Jac(C)`, returning `P`;
   * 1. a point `P` on the curve `C`, returning `[P - oo]`;
   * 2. a pair of points `(P, Q)` on `C`, returning `[P - Q]`;
   * 3. a list of polynomials `(a, b)` with `b^2 + h b - f == 0 mod a`,
   *    returning `[(a(x), y - b(x))]`.
   *
   * Port of `jacobian_homset.py:69-163`.
   */
  __call__(P: JacobianPointInput<C>): JacobianMorphism_divisor_class_field<C> {
    const R = this._poly_ring;

    if (typeof P === 'bigint' || typeof P === 'number') {
      if (BigInt(P) === 0n) {
        return new JacobianMorphism_divisor_class_field<C>(this, [R.one(), R.zero()]);
      }
      throw new TypeError(`argument P (= ${P}) does not determine a divisor class`);
    }

    if (P instanceof JacobianMorphism_divisor_class_field) {
      // Sage: `isinstance(P, JacobianMorphism...) and self == P.parent()`
      if (P.parent === this) {
        return P;
      }
      throw new TypeError(`argument P (= ${P}) does not determine a divisor class`);
    }

    if (Array.isArray(P)) {
      if (P.length === 1) {
        if (BigInt(P[0] as bigint | number) === 0n) {
          return new JacobianMorphism_divisor_class_field<C>(this, [R.one(), R.zero()]);
        }
        throw new TypeError(`argument P (= ${P}) must have length 2`);
      }
      if (P.length === 2) {
        const [P1, P2] = P as [unknown, unknown];
        const isPoly = (v: unknown): v is Polynomial<C> =>
          typeof v === 'object' && v !== null && 'coeffs' in v && 'parent' in v;
        const isInt = (v: unknown): v is bigint | number =>
          typeof v === 'bigint' || typeof v === 'number';
        const isCurvePoint = (v: unknown): v is HyperellipticPoint<C> =>
          typeof v === 'object' && v !== null && 'coords' in v && 'curve' in v;

        if (isInt(P1) && isInt(P2)) {
          return new JacobianMorphism_divisor_class_field<C>(this, [
            R.__call__(this._value_ring.__call__(BigInt(P1)) as C),
            R.__call__(this._value_ring.__call__(BigInt(P2)) as C),
          ]);
        }
        if (isInt(P1) && isPoly(P2)) {
          return new JacobianMorphism_divisor_class_field<C>(this, [
            R.__call__(this._value_ring.__call__(BigInt(P1)) as C),
            P2,
          ]);
        }
        if (isPoly(P1) && isInt(P2)) {
          return new JacobianMorphism_divisor_class_field<C>(this, [
            P1,
            R.__call__(this._value_ring.__call__(BigInt(P2)) as C),
          ]);
        }
        if (isPoly(P1) && isPoly(P2)) {
          return new JacobianMorphism_divisor_class_field<C>(this, [P1, P2] as MumfordPair<C>);
        }
        if (isCurvePoint(P1) && isCurvePoint(P2)) {
          return this.__call__(P1).sub(this.__call__(P2));
        }
      }
      throw new TypeError(`argument P (= ${P}) must have length 2`);
    }

    // A point on the curve: [P - oo]
    const x0 = P.get(0);
    const y0 = P.get(1);
    const x = R.gen();
    return this.__call__([x.sub(R.__call__(x0)), R.__call__(y0)]);
  }
}
