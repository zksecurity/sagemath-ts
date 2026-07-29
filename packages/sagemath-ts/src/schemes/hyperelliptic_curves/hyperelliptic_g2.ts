/**
 * @module sage/schemes/hyperelliptic_curves/hyperelliptic_g2
 * @description Hyperelliptic curves of genus 2 over a general ring
 *
 * Port of: `sage/schemes/hyperelliptic_curves/hyperelliptic_g2.py`
 *
 * Upstream builds the concrete class of a curve dynamically, mixing
 * `HyperellipticCurve_g2` with the base-ring specialisation
 * (`constructor.py:335-368`).  TypeScript has no multiple inheritance, so the
 * genus-2 methods live in free functions here and are re-exposed by the three
 * concrete classes `HyperellipticCurve_g2`,
 * `HyperellipticCurve_g2_FiniteField` and `HyperellipticCurve_g2_RationalField`.
 */

import { NotImplementedError } from '../../errors.js';
import type { Polynomial, RingElement } from '../../rings/polynomial/polynomial_element.js';
import type { HyperellipticBaseRing } from './field_ops.js';
import { HyperellipticCurve_finite_field } from './hyperelliptic_finite_field.js';
import { HyperellipticCurve_generic } from './hyperelliptic_generic.js';
import { HyperellipticCurve_rational_field } from './hyperelliptic_rational_field.js';
import {
  absolute_igusa_invariants_kohel,
  absolute_igusa_invariants_wamelen,
  clebsch_invariants,
  igusa_clebsch_invariants,
} from './invariants.js';
import { HyperellipticJacobian_g2 } from './jacobian_g2.js';

/**
 * `4 f + h^2`, the sextic form whose invariants are the curve's invariants.
 */
function sextic<C extends RingElement>(curve: HyperellipticCurve_generic<C>): Polynomial<C> {
  const [f, h] = curve.hyperelliptic_polynomials();
  const K = curve.base_ring() as HyperellipticBaseRing<C>;
  const four = f.parent.__call__(K.__call__(4n) as C);
  return four.mul(f).add(h.mul(h));
}

/**
 * Return `true` if the curve is an odd degree model
 * (`hyperelliptic_g2.py:16-36`).
 */
export function g2_is_odd_degree<C extends RingElement>(
  curve: HyperellipticCurve_generic<C>
): boolean {
  const [f, h] = curve.hyperelliptic_polynomials();
  const df = f.degree();
  if (h.degree() < 3) {
    return df % 2 === 1;
  }
  if (df < 6) {
    return false;
  }
  const a0 = f.leading_coefficient();
  const c0 = h.leading_coefficient();
  const K = curve.base_ring() as HyperellipticBaseRing<C>;
  return c0
    .mul(c0)
    .add((K.__call__(4n) as C).mul(a0))
    .isZero();
}

/** `hyperelliptic_g2.py:83-118` */
export function g2_clebsch_invariants<C extends RingElement>(
  curve: HyperellipticCurve_generic<C>
): [C, C, C, C] {
  return clebsch_invariants(sextic(curve));
}

/** `hyperelliptic_g2.py:120-155` */
export function g2_igusa_clebsch_invariants<C extends RingElement>(
  curve: HyperellipticCurve_generic<C>
): [C, C, C, C] {
  return igusa_clebsch_invariants(sextic(curve));
}

/** `hyperelliptic_g2.py:157-170` */
export function g2_absolute_igusa_invariants_wamelen<C extends RingElement>(
  curve: HyperellipticCurve_generic<C>
): [C, C, C] {
  return absolute_igusa_invariants_wamelen(sextic(curve));
}

/** `hyperelliptic_g2.py:172-191` */
export function g2_absolute_igusa_invariants_kohel<C extends RingElement>(
  curve: HyperellipticCurve_generic<C>
): [C, C, C] {
  return absolute_igusa_invariants_kohel(sextic(curve));
}

/** The genus-2 methods, shared by the three concrete classes below. */
export interface HyperellipticCurve_g2_methods<C extends RingElement> {
  is_odd_degree(): boolean;
  clebsch_invariants(): [C, C, C, C];
  igusa_clebsch_invariants(): [C, C, C, C];
  absolute_igusa_invariants_wamelen(): [C, C, C];
  absolute_igusa_invariants_kohel(): [C, C, C];
  kummer_morphism(): never;
}

/** Genus 2 curve over a general ring (`hyperelliptic_g2.py:15`). */
export class HyperellipticCurve_g2<C extends RingElement>
  extends HyperellipticCurve_generic<C>
  implements HyperellipticCurve_g2_methods<C>
{
  /** `hyperelliptic_g2.py:38-58` */
  jacobian(): HyperellipticJacobian_g2<C> {
    return new HyperellipticJacobian_g2<C>(this);
  }
  is_odd_degree(): boolean {
    return g2_is_odd_degree(this);
  }
  clebsch_invariants(): [C, C, C, C] {
    return g2_clebsch_invariants(this);
  }
  igusa_clebsch_invariants(): [C, C, C, C] {
    return g2_igusa_clebsch_invariants(this);
  }
  absolute_igusa_invariants_wamelen(): [C, C, C] {
    return g2_absolute_igusa_invariants_wamelen(this);
  }
  absolute_igusa_invariants_kohel(): [C, C, C] {
    return g2_absolute_igusa_invariants_kohel(this);
  }
  /** `hyperelliptic_g2.py:60-81` — needs the Kummer surface. */
  kummer_morphism(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: kummer_morphism (sage.schemes.hyperelliptic_curves.kummer_surface is not ported)'
    );
  }
}

/** Genus 2 curve over a finite field. */
export class HyperellipticCurve_g2_FiniteField<C extends RingElement>
  extends HyperellipticCurve_finite_field<C>
  implements HyperellipticCurve_g2_methods<C>
{
  /** `hyperelliptic_g2.py:38-58` */
  jacobian(): HyperellipticJacobian_g2<C> {
    return new HyperellipticJacobian_g2<C>(this);
  }
  is_odd_degree(): boolean {
    return g2_is_odd_degree(this);
  }
  clebsch_invariants(): [C, C, C, C] {
    return g2_clebsch_invariants(this);
  }
  igusa_clebsch_invariants(): [C, C, C, C] {
    return g2_igusa_clebsch_invariants(this);
  }
  absolute_igusa_invariants_wamelen(): [C, C, C] {
    return g2_absolute_igusa_invariants_wamelen(this);
  }
  absolute_igusa_invariants_kohel(): [C, C, C] {
    return g2_absolute_igusa_invariants_kohel(this);
  }
  kummer_morphism(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: kummer_morphism (sage.schemes.hyperelliptic_curves.kummer_surface is not ported)'
    );
  }
}

/** Genus 2 curve over the rationals. */
export class HyperellipticCurve_g2_RationalField<C extends RingElement>
  extends HyperellipticCurve_rational_field<C>
  implements HyperellipticCurve_g2_methods<C>
{
  /** `hyperelliptic_g2.py:38-58` */
  jacobian(): HyperellipticJacobian_g2<C> {
    return new HyperellipticJacobian_g2<C>(this);
  }
  is_odd_degree(): boolean {
    return g2_is_odd_degree(this);
  }
  clebsch_invariants(): [C, C, C, C] {
    return g2_clebsch_invariants(this);
  }
  igusa_clebsch_invariants(): [C, C, C, C] {
    return g2_igusa_clebsch_invariants(this);
  }
  absolute_igusa_invariants_wamelen(): [C, C, C] {
    return g2_absolute_igusa_invariants_wamelen(this);
  }
  absolute_igusa_invariants_kohel(): [C, C, C] {
    return g2_absolute_igusa_invariants_kohel(this);
  }
  kummer_morphism(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: kummer_morphism (sage.schemes.hyperelliptic_curves.kummer_surface is not ported)'
    );
  }
}
