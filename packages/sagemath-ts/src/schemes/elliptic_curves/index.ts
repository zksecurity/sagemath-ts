/**
 * @module sage/schemes/elliptic_curves
 * @description Elliptic curves over finite fields
 *
 * This module provides functionality for:
 * - Creating elliptic curves in short and general Weierstrass form
 * - Computing curve invariants (b2, b4, b6, b8, c4, c6, discriminant, j-invariant)
 * - Point arithmetic on elliptic curves
 * - Computing curve cardinality (order)
 * - Finding generators and point orders
 *
 * There are two implementations:
 *
 * 1. **EllipticCurveFiniteField** (ell_finite_field.ts):
 *    - Optimized for finite fields (prime fields)
 *    - Supports cardinality computation, random points, generators
 *    - Short Weierstrass form: y^2 = x^3 + ax + b
 *
 * 2. **EllipticCurveGeneric** (ell_generic.ts):
 *    - Supports general Weierstrass form over any field
 *    - Computes all invariants (b2, b4, b6, b8, c4, c6)
 *    - More general but fewer finite-field-specific optimizations
 *
 * @example
 * ```typescript
 * // Using the finite field optimized implementation
 * import { EllipticCurve as EllipticCurveFF, EllipticCurveFiniteField } from './ell_finite_field.js';
 * import { FiniteFieldPrime } from '../rings/finite_rings/finite_field_prime.js';
 *
 * const F = new FiniteFieldPrime(23n);
 * const E = EllipticCurveFF(F, [1n, 1n]);  // y^2 = x^3 + x + 1
 *
 * console.log('Curve order:', E.cardinality());
 * const P = E.random_point();
 * console.log('Point order:', P.order());
 * ```
 */

// Finite field optimized implementation
export {
  EllipticCurve as EllipticCurveFF,
  EllipticCurveFiniteField,
  EllipticCurvePoint as EllipticCurvePointFF,
  type FieldElement as FieldElementFF,
  type BaseField,
  // Group structure functions
  abelian_group,
  discrete_log,
  embedding_degree,
  is_supersingular,
  is_ordinary,
  is_isogenous,
  // Torsion subgroup functions
  torsion_basis,
  torsion_subgroup as ff_torsion_subgroup,
  division_points as ff_division_points,
  _find_point_of_order,
  _is_independent,
  // Isogeny functions
  multiplication_by_m_isogeny,
  ScalarMultiplicationIsogeny,
  // Stub functions
  count_points,
  frobenius_polynomial,
  frobenius_order,
  frobenius_endomorphism,
  frobenius_discriminant,
  height_above_floor,
  endomorphism_discriminant_from_class_number,
  endomorphism_order,
  twists,
  multiplication_by_p_isogeny,
} from './ell_finite_field.js';

// Re-export EllipticCurve from ell_finite_field as the default
// (users can import the generic one explicitly if needed)
export { EllipticCurve } from './ell_finite_field.js';
export { EllipticCurvePoint } from './ell_finite_field.js';

// Generic elliptic curve implementation
export { EllipticCurveGeneric } from './ell_generic.js';

// Point functions (stubs)
export {
  weil_pairing,
  tate_pairing,
  ate_pairing,
  division_points,
  has_finite_order,
  has_infinite_order,
  height,
  archimedean_local_height,
  non_archimedean_local_height,
  elliptic_logarithm,
  padic_elliptic_logarithm,
  point_has_good_reduction,
  reduction,
  is_divisible_by,
  point_log,
} from './ell_point.js';

// Isogeny class module
export {
  IsogenyClass,
  IsogenyClassNumberField,
  IsogenyClassRational,
  isogeny_degrees_cm,
  possible_isogeny_degrees,
} from './isogeny_class.js';

// Complex multiplication module
export {
  hilbert_class_polynomial,
  is_HCP,
  cm_j_invariants,
  cm_j_invariants_and_orders,
  cm_orders,
  largest_fundamental_disc_with_class_number,
  largest_disc_with_class_number,
  discriminants_with_bounded_class_number,
  is_cm_j_invariant,
} from './cm.js';

// p-adic L-series module
export {
  pAdicLseries,
  pAdicLseriesOrdinary,
  pAdicLseriesSupersingular,
} from './padic_lseries.js';

// Formal group module
export { EllipticCurveFormalGroup } from './formal_group.js';

// Isogeny module
export {
  EllipticCurveIsogeny,
  _isogeny_determine_algorithm,
  isogeny_codomain_from_kernel,
  compute_codomain_formula,
  compute_vw_kohel_even_deg1,
  compute_vw_kohel_even_deg3,
  compute_vw_kohel_odd,
  compute_codomain_kohel,
  two_torsion_part,
  compute_isogeny_bmss,
  compute_isogeny_stark,
  compute_isogeny_kernel_polynomial,
  compute_intermediate_curves,
  compute_sequence_of_maps,
  fill_isogeny_matrix,
  unfill_isogeny_matrix,
  EllipticCurveIsogeny_from_kernel_gens,
  EllipticCurveIsogeny_from_kernel_polynomial,
} from './ell_curve_isogeny.js';

// Weierstrass morphism module
export {
  baseWI,
  WeierstrassIsomorphism,
  _isomorphisms,
  identity_morphism,
  negation_morphism,
} from './weierstrass_morphism.js';

// Torsion subgroup module
export {
  EllipticCurveTorsionSubgroup,
  torsion_bound,
  torsion_order,
  torsion_points,
  torsion_subgroup,
  _p_primary_torsion_basis,
  has_finite_order as torsion_has_finite_order,
  has_infinite_order as torsion_has_infinite_order,
  point_order,
  order_from_multiple,
  discrete_log as torsion_discrete_log,
} from './ell_torsion.js';
