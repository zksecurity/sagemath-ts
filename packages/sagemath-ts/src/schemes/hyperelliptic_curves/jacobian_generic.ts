/**
 * @module sage/schemes/hyperelliptic_curves/jacobian_generic
 * @description Jacobian of a general hyperelliptic curve
 *
 * Port of: `sage/schemes/hyperelliptic_curves/jacobian_generic.py`
 */

import { NotImplementedError } from '../../errors.js';
import type { RingElement } from '../../rings/polynomial/polynomial_element.js';
import type { HyperellipticBaseRing } from './field_ops.js';
import type { HyperellipticCurve_finite_field } from './hyperelliptic_finite_field.js';
import {
  type HyperellipticCurve_generic,
  _register_jacobian_module,
} from './hyperelliptic_generic.js';
import { JacobianHomset_divisor_classes } from './jacobian_homset.js';
import type { JacobianPointInput } from './jacobian_homset.js';
import type { JacobianMorphism_divisor_class_field } from './jacobian_morphism.js';

/**
 * The Jacobian of a hyperelliptic curve.
 *
 * Port of `jacobian_generic.py:22-436`.
 */
export class HyperellipticJacobian_generic<C extends RingElement> {
  private readonly _curve: HyperellipticCurve_generic<C>;
  private readonly _homsets = new Map<unknown, JacobianHomset_divisor_classes<C>>();

  constructor(curve: HyperellipticCurve_generic<C>) {
    this._curve = curve;
  }

  /** The curve this is the Jacobian of. */
  curve(): HyperellipticCurve_generic<C> {
    return this._curve;
  }

  /** The base ring of the curve. */
  base_ring(): HyperellipticBaseRing<C> {
    return this._curve.base_ring();
  }

  /**
   * Return the dimension of this Jacobian (`jacobian_generic.py:145-162`).
   */
  dimension(): bigint {
    return BigInt(this._curve.genus());
  }

  /**
   * `J(S)`: the set of `S`-rational points (`jacobian_generic.py:170-171`).
   */
  point_homset(S?: HyperellipticBaseRing<C>): JacobianHomset_divisor_classes<C> {
    const ring = S ?? this.base_ring();
    const cached = this._homsets.get(ring);
    if (cached !== undefined) {
      return cached;
    }
    const homset = new JacobianHomset_divisor_classes<C>(this, ring);
    this._homsets.set(ring, homset);
    return homset;
  }

  /**
   * Return a point of the Jacobian from a Mumford divisor
   * (`jacobian_generic.py:164-168`).
   */
  point(
    mumford: JacobianPointInput<C>,
    _options?: { check?: boolean }
  ): JacobianMorphism_divisor_class_field<C> {
    return this.point_homset().__call__(mumford);
  }

  /**
   * Sage's `J(x)`: a ring gives the point set `J(x)`, anything else a point.
   */
  __call__(
    x: HyperellipticBaseRing<C> | JacobianPointInput<C>
  ): JacobianHomset_divisor_classes<C> | JacobianMorphism_divisor_class_field<C> {
    if (is_ring(x)) {
      return this.point_homset(x as HyperellipticBaseRing<C>);
    }
    return this.point(x as JacobianPointInput<C>);
  }

  /**
   * Return the cardinality of the Jacobian
   * (`jacobian_generic.py:420-436`): `frobenius_polynomial()(1)`.
   *
   * Currently only implemented over finite fields.
   */
  cardinality(): bigint {
    const curve = this._curve as unknown as HyperellipticCurve_finite_field<C>;
    if (typeof curve.frobenius_polynomial !== 'function') {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: cardinality of a Jacobian over a non-finite base ring'
      );
    }
    const P = curve.frobenius_polynomial();
    let s = 0n;
    for (const c of P) {
      s += c;
    }
    return s;
  }

  /** `jacobian_generic.py:238-309` — needs `genus2reduction` and [Lom2019]. */
  geometric_endomorphism_algebra_is_field(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: geometric_endomorphism_algebra_is_field ' +
        '(sage.interfaces.genus2reduction and jacobian_endomorphism_utils are not ported)'
    );
  }

  /** `jacobian_generic.py:311-418` — needs `genus2reduction` and [Lom2019]. */
  geometric_endomorphism_ring_is_ZZ(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: geometric_endomorphism_ring_is_ZZ ' +
        '(sage.interfaces.genus2reduction and jacobian_endomorphism_utils are not ported)'
    );
  }

  eq(other: HyperellipticJacobian_generic<C>): boolean {
    return this._curve.eq(other._curve);
  }

  toString(): string {
    return `Jacobian of ${this._curve}`;
  }
}

/** Heuristic test for "this argument is a ring, not a point". */
function is_ring(x: unknown): boolean {
  if (typeof x !== 'object' || x === null) {
    return false;
  }
  const probe = x as { zero?: unknown; one?: unknown; __call__?: unknown };
  return (
    typeof probe.zero === 'function' &&
    typeof probe.one === 'function' &&
    typeof probe.__call__ === 'function'
  );
}

_register_jacobian_module({ HyperellipticJacobian_generic });
