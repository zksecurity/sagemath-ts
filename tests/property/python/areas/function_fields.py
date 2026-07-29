"""Live-SageMath oracle for rational function fields over prime fields."""

from sage.all import *


def _ints(xs):
    return ','.join(str(ZZ(x)) for x in xs)


def _element(K, coeffs):
    return K(K._ring(coeffs))


def _parts(f):
    return '%s/%s' % (_ints(f.numerator().list()), _ints(f.denominator().list()))


def ff_arithmetic(p, numerator, denominator, exponent):
    K = FunctionField(GF(ZZ(p)), 'x')
    f = _element(K, numerator) / _element(K, denominator)
    g = f ** ZZ(exponent)
    return 'f=%s pow=%s degree=%s vx=%s square=%s' % (
        _parts(f),
        _parts(g),
        f.degree(),
        f.valuation(K.gen()),
        f.is_square(),
    )


def ff_factor(p, numerator, denominator):
    K = FunctionField(GF(ZZ(p)), 'x')
    f = _element(K, numerator) / _element(K, denominator)
    F = f.factor()
    factors = ';'.join('%s:%s' % (_parts(q), e) for q, e in F)
    return 'unit=%s factors=%s' % (_parts(K(F.unit())), factors)


def ff_places(p, degree):
    K = FunctionField(GF(ZZ(p)), 'x')
    places = K.places(int(degree))
    return ';'.join(str(P) for P in places)


FUNCTIONS = {
    'ff_arithmetic': ff_arithmetic,
    'ff_factor': ff_factor,
    'ff_places': ff_places,
}
