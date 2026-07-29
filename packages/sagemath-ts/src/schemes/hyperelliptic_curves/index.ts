/**
 * @module sage/schemes/hyperelliptic_curves
 * @description Hyperelliptic curves, their Jacobians and Mumford divisors
 *
 * Port of `sage/schemes/hyperelliptic_curves/`.
 *
 * @example
 * ```typescript
 * import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
 * import { PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
 * import { HyperellipticCurve } from './constructor.js';
 *
 * const K = GF(37n);
 * const R = new PolynomialRing(K, 'x');
 * const x = R.gen();
 * const H = HyperellipticCurve(x.pow(5).add(x.pow(4).scalar_mul(K.__call__(12n))));
 * const J = H.jacobian();
 * const P = J.point(H.lift_x(2n));
 * console.log(String(P.add(P)));
 * ```
 */

export {
  HyperellipticCurve,
  _parse_multivariate_defining_equation,
  type HyperellipticCurveOptions,
} from './constructor.js';

export {
  type HyperellipticBaseRing,
  absolute_trace_is_zero,
  cardinality_of,
  characteristic_of,
  compare_elements,
  degree_of,
  div_elements,
  element_to_bigint,
  field_embedding,
  is_finite_field,
  is_square_of,
  iterate_field,
  pow_element,
  sort_roots_like_sage,
  sqrt_all_of,
} from './field_ops.js';

export {
  HyperellipticCurve_generic,
  HyperellipticPoint,
  genus_of,
  poly_repr,
  sage_poly_repr,
} from './hyperelliptic_generic.js';

export {
  type CartierData,
  HyperellipticCurve_finite_field,
  type ZZPoly,
  ZZRationalFunction,
  zz_poly_repr,
} from './hyperelliptic_finite_field.js';

export { HyperellipticCurve_rational_field } from './hyperelliptic_rational_field.js';

export {
  HyperellipticCurve_g2,
  HyperellipticCurve_g2_FiniteField,
  HyperellipticCurve_g2_RationalField,
  type HyperellipticCurve_g2_methods,
  g2_absolute_igusa_invariants_kohel,
  g2_absolute_igusa_invariants_wamelen,
  g2_clebsch_invariants,
  g2_igusa_clebsch_invariants,
  g2_is_odd_degree,
} from './hyperelliptic_g2.js';

export {
  type BinaryForm,
  type UbsResult,
  Ueberschiebung,
  absolute_igusa_invariants_kohel,
  absolute_igusa_invariants_wamelen,
  clebsch_invariants,
  clebsch_to_igusa,
  diffxy,
  igusa_clebsch_invariants,
  igusa_to_clebsch,
  sextic_form,
  ubs,
} from './invariants.js';

export { HyperellipticJacobian_generic } from './jacobian_generic.js';
export { HyperellipticJacobian_g2 } from './jacobian_g2.js';
export {
  JacobianHomset_divisor_classes,
  type JacobianPointInput,
} from './jacobian_homset.js';
export {
  JacobianMorphism_divisor_class_field,
  type MumfordPair,
  cantor_composition,
  cantor_composition_simple,
  cantor_reduction,
  cantor_reduction_simple,
} from './jacobian_morphism.js';
