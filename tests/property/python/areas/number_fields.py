"""SageMath side of the ``number_fields`` property-test area.

Cases: tests/property/cases/number_fields.cases.json
"""

from sage.all import *


def nf_quadratic_discriminant(d):
    """Compute the discriminant of Q(sqrt(d))."""
    K = QuadraticField(d)
    return K.discriminant()


def nf_degree(coeffs):
    """Compute the degree of a number field defined by polynomial with given coefficients."""
    R = PolynomialRing(QQ, 'x')
    f = R(list(map(QQ, coeffs)))
    K = NumberField(f, 'a')
    return K.degree()


def nf_quadratic_signature(d):
    """Compute the signature (r1, r2) of Q(sqrt(d))."""
    K = QuadraticField(d)
    return K.signature()


def nf_is_totally_real(d):
    """Check if Q(sqrt(d)) is totally real."""
    K = QuadraticField(d)
    return K.is_totally_real()


def nf_polynomial_discriminant(coeffs):
    """Compute the discriminant of a polynomial."""
    R = PolynomialRing(QQ, 'x')
    f = R(list(map(QQ, coeffs)))
    return f.discriminant()


def nf_fundamental_discriminant(d):
    """Compute the fundamental discriminant associated to d."""
    return fundamental_discriminant(d)


FUNCTIONS = {
    'quadratic_discriminant': nf_quadratic_discriminant,
    'nf_degree': nf_degree,
    'quadratic_signature': nf_quadratic_signature,
    'is_totally_real': nf_is_totally_real,
    'polynomial_discriminant': nf_polynomial_discriminant,
    'fundamental_discriminant_test': nf_fundamental_discriminant,
}

from number_field_coercion import nf_coercion, nf_zero, nf_float
FUNCTIONS.update(nf_coercion=nf_coercion, nf_zero=nf_zero, nf_float=nf_float)

from number_field_index import nf_index
FUNCTIONS.update(nf_index=nf_index)

from number_field_element import nf_element
from number_field_power import power_case
FUNCTIONS.update(nf_element=nf_element, nf_generic_power=power_case)

from nf_scaled_unit import nf_scaled_unit
FUNCTIONS.update(nf_scaled_unit=nf_scaled_unit)

from nf_scaled_consumer import nf_scaled_consumer
FUNCTIONS.update(nf_scaled_consumer=nf_scaled_consumer)

from number_field_ideal import nf_ideal_audit
FUNCTIONS.update(nf_ideal_audit=nf_ideal_audit)

from number_field_ideal_basis import nf_ideal_basis
FUNCTIONS.update(nf_ideal_basis=nf_ideal_basis)

from number_field_ideal_binary import nf_ideal_binary
from number_field_ideal_construction import nf_ideal_construction
FUNCTIONS.update(nf_ideal_binary=nf_ideal_binary,nf_ideal_construction=nf_ideal_construction)

from number_field_ideal_inputs import nf_ideal_inputs, nf_cross_ideal
FUNCTIONS.update(nf_ideal_inputs=nf_ideal_inputs, nf_cross_ideal=nf_cross_ideal)

from number_field_ideal_inputs import nf_ideal_factory, nf_ideal_intersection_rng
from pari_ideal_intersection import native_pari_ideal_intersection
FUNCTIONS.update(nf_ideal_factory=nf_ideal_factory, nf_ideal_intersection_rng=nf_ideal_intersection_rng, pari_ideal_intersection=native_pari_ideal_intersection)

from number_field_ideal_methods import nf_ideal_method_input, nf_field_vector
from pari_ideal_adapters import pari_hnfcenter, pari_diviiround, pari_centered_basis, pari_idealmul
FUNCTIONS.update(nf_ideal_method_input=nf_ideal_method_input, nf_field_vector=nf_field_vector,
                 pari_hnfcenter=pari_hnfcenter, pari_diviiround=pari_diviiround,
                 pari_centered_basis=pari_centered_basis, pari_idealmul=pari_idealmul)

from number_field_ideal_methods import nf_basis_class_number
FUNCTIONS.update(nf_basis_class_number=nf_basis_class_number)

from number_field_ideal_classes import nf_ideal_class_method, nf_prime_ideal_class
FUNCTIONS.update(nf_ideal_class_method=nf_ideal_class_method,nf_prime_ideal_class=nf_prime_ideal_class)

from pari_ideal_valuation import pari_prime_valuation_data, pari_ideal_valuation
FUNCTIONS.update({"pari_prime_valuation_data":pari_prime_valuation_data,"pari_ideal_valuation":pari_ideal_valuation})

from number_field_valuation import nf_ideal_valuation,nf_ideal_valuation_input
FUNCTIONS.update({"nf_ideal_valuation":nf_ideal_valuation,"nf_ideal_valuation_input":nf_ideal_valuation_input})

from pari_ideal_valuation import pari_prime_element_valuation
FUNCTIONS["pari_prime_element_valuation"]=pari_prime_element_valuation

from number_field_valuation import nf_prime_below
FUNCTIONS["nf_prime_below"]=nf_prime_below

from pari_full_ideals import pari_full_primedec
FUNCTIONS["pari_full_primedec"]=pari_full_primedec

from pari_full_ideals import pari_full_ideal
FUNCTIONS["pari_full_ideal"]=pari_full_ideal

from number_field_factor import nf_ideal_factor
FUNCTIONS["nf_ideal_factor"]=nf_ideal_factor

from number_field_factor import nf_factor_properties
FUNCTIONS['nf_factor_properties']=nf_factor_properties

from pari_full_ideals import pari_full_primedec_no_index
FUNCTIONS["pari_full_primedec_no_index"]=pari_full_primedec_no_index

from pari_rational_polynomials import pari_rational_polynomial
from number_field_trace_norm import nf_trace_norm
FUNCTIONS.update(pari_rational_polynomial=pari_rational_polynomial,nf_trace_norm=nf_trace_norm)

from number_field_trace_base import nf_base_trace_norm
FUNCTIONS['nf_base_trace_norm']=nf_base_trace_norm

from nf_pari_factorization import nf_pari_factor_constructor,nf_irreducible_cache
FUNCTIONS.update(nf_pari_factor_constructor=nf_pari_factor_constructor,nf_irreducible_cache=nf_irreducible_cache)

from number_field_legacy_factor import nf_legacy_polynomial_factor
FUNCTIONS['nf_legacy_polynomial_factor'] = nf_legacy_polynomial_factor

from number_field_automorphism_order import nf_automorphism_order
FUNCTIONS['nf_automorphism_order'] = nf_automorphism_order

from number_field_conjugate_bound import nf_conjugate_bound
FUNCTIONS['nf_conjugate_bound'] = nf_conjugate_bound
