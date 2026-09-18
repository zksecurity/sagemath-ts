"""Bundled coordinate callers over distinct source/target parents."""
import ast
import copy
import json
from sage.all import GF, QQ, ZZ, PolynomialRing, EllipticCurve
from sage.structure.coerce import py_scalar_to_element
from curve_coordinate_roots import _methods, _path
from finite_polynomial_roots import _roots_univariate_polynomial


def _field(kind):
    if kind == 0: return QQ
    primes = {1: 2, 2: 3, 3: 5, 4: 7, 12: 3, 13: 2}
    if kind in primes: return GF(primes[kind])
    p, d, T, name = {5: (3, 2, [1, 0, 1], 'a'), 6: (3, 2, [1, 0, 1], 'b'),
                    7: (3, 2, [2, 1, 1], 'a'), 8: (2, 2, [1, 1, 1], 'a'),
                    9: (2, 2, [1, 1, 1], 'b'), 10: (2, 3, [1, 1, 0, 1], 'a'),
                    11: (3, 3, [1, 2, 0, 1], 'a')}[kind]
    return GF(p**d, name, modulus=PolynomialRing(GF(p), 't')(T), impl='pari_ffelt')


class _Polynomial:
    def __init__(self, f): self.f = f
    def roots(self, ring=None, multiplicities=True):
        K = self.f.base_ring()
        if K is QQ:
            from polynomial_roots import roots
            return roots(self.f, multiplicities=multiplicities)
        return _roots_univariate_polynomial(K, self.f, ring, multiplicities)


_namespace = dict(py_scalar_to_element=py_scalar_to_element,
                  PolynomialRing=lambda K, name: lambda cs: _Polynomial(PolynomialRing(K, name)(cs)))
for name in ('is_x_coord', 'lift_x'):
    node = copy.deepcopy(_methods[name]); node.decorator_list = []
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(_path), 'exec'), _namespace)


def ec_coordinate_coercion(target, source, numerator, denominator, operation):
    try:
        K = _field(target)
        E = EllipticCurve(K, [1, 0, 0, 0, 1] if K.characteristic() == 2 else [0, 0, 0, 1, 0])
        if source == 14: x = ZZ(numerator)
        elif source == 15: x = int(numerator)
        elif source == 16: x = float(numerator) / int(denominator)
        elif source == 17: x = str(numerator)
        elif source == 18: x = None
        elif source == 19: x = object()
        elif source == 20: x = bool(numerator)
        elif source == 21: x = float('nan')
        elif source == 22: x = float('inf')
        elif source == 23: x = -float('inf')
        elif source == 24: x = 'not an integer'
        elif source == 25: x = '1/2'
        elif source == 26: x = None  # undefined uses the constructor's zero default
        else:
            L = _field(source)
            if L is QQ: x = QQ(numerator) / denominator
            elif L.degree() == 1: x = L(numerator)
            else:
                n = int(numerator) % int(L.order()); ds = []
                for _ in range(L.degree()): ds.append(n % L.characteristic()); n //= L.characteristic()
                x = L(ds)
        if operation == 4: value = str(ZZ(x))
        elif operation == 3: value = str(K(x))
        elif operation == 0: value = bool(_namespace['is_x_coord'](E, x))
        else:
            pts = _namespace['lift_x'](E, x, all=operation == 1)
            if operation == 2: pts = [pts]
            value = [[str(P[0]), str(P[1]), str(P.curve().base_ring()),
                      list(map(str, P.curve().ainvs())), P.curve() is E] for P in pts]
        return json.dumps({'value': value}, separators=(',', ':'))
    except Exception as error:
        return json.dumps({'error': type(error).__name__, 'message': str(error)}, separators=(',', ':'))
