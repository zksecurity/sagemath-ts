"""Live-SageMath oracle for quaternion algebras over QQ."""

from sage.all import *


def _vals(xs):
    return ','.join(str(x) for x in xs)


def quat_element(a, b, coeffs):
    A = QuaternionAlgebra(QQ, ZZ(a), ZZ(b), names=('i', 'j', 'k'))
    x = A([QQ(c) for c in coeffs])
    inv = '-' if x == 0 else _vals(list(~x))
    charpoly = _vals(x.reduced_characteristic_polynomial().list())
    return 'x=%s conj=%s tr=%s norm=%s inv=%s cp=%s' % (
        _vals(list(x)),
        _vals(list(x.conjugate())),
        x.reduced_trace(),
        x.reduced_norm(),
        inv,
        charpoly,
    )


def quat_product(a, b, left, right):
    A = QuaternionAlgebra(QQ, ZZ(a), ZZ(b), names=('i', 'j', 'k'))
    x = A([QQ(c) for c in left])
    y = A([QQ(c) for c in right])
    z = x * y
    return 'xy=%s pair=%s norm=%s anti=%s' % (
        _vals(list(z)),
        x.pair(y),
        z.reduced_norm(),
        (x * y).conjugate() == y.conjugate() * x.conjugate(),
    )


def quat_algebra(a, b):
    A = QuaternionAlgebra(QQ, ZZ(a), ZZ(b), names=('i', 'j', 'k'))
    return 'inv=%s disc=%s ram=%s definite=%s' % (
        _vals(A.invariants()),
        A.discriminant(),
        _vals(A.ramified_primes()),
        A.invariants()[0] < 0 and A.invariants()[1] < 0,
    )


FUNCTIONS = {
    'quat_element': quat_element,
    'quat_product': quat_product,
    'quat_algebra': quat_algebra,
}
