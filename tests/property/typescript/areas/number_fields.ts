import { nf_ideal_binary } from '../number_field_ideal_binary.js';
import { nf_ideal_construction } from '../number_field_ideal_construction.js';
import { nf_ideal_basis } from '../number_field_ideal_basis.js';
import { nf_ideal_audit } from '../number_field_ideal.js';
import { nf_scaled_consumer } from '../nf_scaled_consumer.js';
import { nf_scaled_unit } from '../nf_scaled_unit.js';
import { nf_element } from '../number_field_element.js';
import { power_case as nf_generic_power } from '../number_field_power.js';
import { nf_index } from '../number_field_index.js';
import { nf_coercion, nf_zero, nf_float } from '../number_field_coercion.js';
/**
 * sagemath-ts side of the `number_fields` property-test area.
 *
 * Cases: tests/property/cases/number_fields.cases.json
 * SageMath counterpart: tests/property/python/areas/number_fields.py
 */

import { fundamental_discriminant } from '../../../../packages/sagemath-ts/src/arith/misc.js';
import {
  NumberField,
  QuadraticField,
  RationalPolynomial,
} from '../../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../../packages/sagemath-ts/src/rings/rational.js';

export const functions = {
  nf_ideal_binary,
  nf_ideal_construction,
  nf_ideal_basis,
  nf_ideal_audit,
  nf_scaled_consumer,
  nf_scaled_unit,
  nf_element,
  nf_generic_power,
  nf_index,
  nf_coercion,
  nf_zero,
  nf_float,
  quadratic_discriminant: (d: bigint) => {
    const K = QuadraticField.create(d);
    return K.discriminant();
  },
  nf_degree: (coeffs: bigint[]) => {
    const poly = new RationalPolynomial(coeffs.map((c) => new Rational(c)));
    const K = new NumberField(poly, 'a');
    return BigInt(K.degree());
  },
  quadratic_signature: (d: bigint) => {
    const K = QuadraticField.create(d);
    const sig = K.signature();
    return sig;
  },
  is_totally_real: (d: bigint) => {
    const K = QuadraticField.create(d);
    return K.is_totally_real();
  },
  polynomial_discriminant: (coeffs: bigint[]) => {
    const poly = new RationalPolynomial(coeffs.map((c) => new Rational(c)));
    const disc = poly.discriminant();
    // Return as integer (numerator/denominator)
    return disc.numerator / disc.denominator;
  },
  fundamental_discriminant_test: (d: bigint) => fundamental_discriminant(d),
};

import { nf_ideal_inputs, nf_cross_ideal } from '../number_field_ideal_inputs.js';
Object.assign(functions, { nf_ideal_inputs, nf_cross_ideal });

import { nf_ideal_factory, nf_ideal_intersection_rng } from '../number_field_ideal_inputs.js';
import { pari_ideal_intersection } from '../pari_ideal_intersection.js';
Object.assign(functions, { nf_ideal_factory, nf_ideal_intersection_rng, pari_ideal_intersection });

import { nf_ideal_method_input, nf_field_vector } from '../number_field_ideal_methods.js';
import {
  pari_hnfcenter,
  pari_diviiround,
  pari_centered_basis,
  pari_idealmul,
} from '../pari_ideal_adapters.js';
Object.assign(functions, {
  nf_ideal_method_input,
  nf_field_vector,
  pari_hnfcenter,
  pari_diviiround,
  pari_centered_basis,
  pari_idealmul,
});

import { nf_basis_class_number } from '../number_field_ideal_methods.js';
Object.assign(functions, { nf_basis_class_number });

import {nf_ideal_class_method,nf_prime_ideal_class}from '../number_field_ideal_classes.js';
Object.assign(functions,{nf_ideal_class_method,nf_prime_ideal_class});

import { pari_prime_valuation_data, pari_ideal_valuation } from "../pari_ideal_valuation.js";
Object.assign(functions,{pari_prime_valuation_data,pari_ideal_valuation});

import {nf_ideal_valuation,nf_ideal_valuation_input} from "../number_field_valuation.js";
Object.assign(functions,{nf_ideal_valuation,nf_ideal_valuation_input});

import {pari_prime_element_valuation}from "../pari_ideal_valuation.js";
Object.assign(functions,{pari_prime_element_valuation});

import {nf_prime_below}from "../number_field_valuation.js";
Object.assign(functions,{nf_prime_below});

import {pari_full_primedec}from "../pari_full_ideals.js";
Object.assign(functions,{pari_full_primedec});

import {pari_full_ideal}from "../pari_full_ideals.js";
Object.assign(functions,{pari_full_ideal});

import {nf_ideal_factor}from "../number_field_factor.js";
Object.assign(functions,{nf_ideal_factor});

import {nf_factor_properties}from '../number_field_factor.js';
Object.assign(functions,{nf_factor_properties});

import {pari_full_primedec_no_index}from "../pari_full_ideals.js";
Object.assign(functions,{pari_full_primedec_no_index});

import {pari_rational_polynomial}from '../pari_rational_polynomials.js';
import {nf_trace_norm}from '../number_field_trace_norm.js';
Object.assign(functions,{pari_rational_polynomial,nf_trace_norm});

import {nf_base_trace_norm}from '../number_field_trace_base.js';
Object.assign(functions,{nf_base_trace_norm});

import {nf_pari_factor_constructor,nf_irreducible_cache}from "../nf_pari_factorization.js";
Object.assign(functions,{nf_pari_factor_constructor,nf_irreducible_cache});

import { nf_legacy_polynomial_factor } from '../number_field_legacy_factor.js';
Object.assign(functions, { nf_legacy_polynomial_factor });

import { nf_automorphism_order } from '../number_field_automorphism_order.js';
Object.assign(functions, { nf_automorphism_order });

import { nf_conjugate_bound } from '../number_field_conjugate_bound.js';
Object.assign(functions, { nf_conjugate_bound });
