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
