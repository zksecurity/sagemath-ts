"""Bundled coordinate callers over distinct source/target parents."""
import ast
import copy
import json
from sage.all import GF, QQ, ZZ, Zmod, PolynomialRing, EllipticCurve
from sage.structure.coerce import py_scalar_to_element
from curve_coordinate_roots import _methods, _path
from finite_polynomial_roots import _roots_univariate_polynomial


def _field(kind):
    if kind == 0: return QQ
    primes = {1: 2, 2: 3, 3: 5, 4: 7, 12: 3, 13: 2, 14: 257, 15: 65537}
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
        elif source == 27: x = Zmod(9)(numerator)
        elif source == 28: x = Zmod(6)(numerator)
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
            namespace = _extension_namespace if operation in (7, 8) else _namespace
            pts = namespace['lift_x'](E, x, all=operation in (1, 7), extend=operation in (7, 8))
            if operation in (2, 8): pts = [pts]
            value = [[str(P[0]), str(P[1]), str(P.curve().base_ring()),
                      list(map(str, P.curve().ainvs())), P.curve() is E] for P in pts]
        return json.dumps({'value': value}, separators=(',', ':'))
    except Exception as error:
        return json.dumps({'error': type(error).__name__, 'message': str(error)}, separators=(',', ':'))


# Execute the bundled base-change body, including construction and parent handling.
from sage.schemes.elliptic_curves import constructor
from sage.rings.finite_rings.finite_field_base import FiniteField
_base_namespace = dict(constructor=constructor, FiniteField=FiniteField)
_node = copy.deepcopy(_methods['base_extend']); _node.decorator_list = []
exec(compile(ast.Module(body=[_node], type_ignores=[]), str(_path), 'exec'), _base_namespace)


def ec_curve_base_change(source, target, numerator, denominator, operation):
    calls = []
    try:
        K, L = _field(source), _field(target)
        if K is QQ: c = QQ(numerator) / denominator
        elif K.degree() == 1: c = K(numerator)
        else:
            n = int(numerator) % int(K.order()); ds = []
            for _ in range(K.degree()): ds.append(n % K.characteristic()); n //= K.characteristic()
            c = K(ds)
        E = EllipticCurve(K, [1, 0, 0, c, 1] if K.characteristic() == 2 else [0, 0, 0, c, 1])
        # change_ring delegates directly to this same base_extend body.
        def convert(a):
            calls.append(str(a.parent()))
            return L(a)
        changed = _base_namespace['base_extend'](E, convert if operation >= 2 else L)
        result = {'value': [str(changed.base_ring()), list(map(str, changed.ainvs())), changed is E]}
    except Exception as error:
        result = {'error': type(error).__name__, 'message': str(error)}
    if operation >= 2: result['calls'] = calls
    return json.dumps(result, separators=(',', ':'))


_extension_namespace = dict(py_scalar_to_element=py_scalar_to_element, PolynomialRing=PolynomialRing)
_node = copy.deepcopy(_methods['lift_x']); _node.decorator_list = []
exec(compile(ast.Module(body=[_node], type_ignores=[]), str(_path), 'exec'), _extension_namespace)


def ec_lift_extension(kind, coefficients, numerator, denominator, all_points, extend):
    try:
        K = _field(kind)
        def decode(n):
            if K is QQ or K.degree() == 1: return K(n)
            n = int(n) % int(K.order()); ds = []
            for _ in range(K.degree()): ds.append(n % K.characteristic()); n //= K.characteristic()
            return K(ds)
        E = EllipticCurve(K, list(map(decode, coefficients)))
        x = decode(numerator) / denominator if K is QQ else decode(numerator)
        points = _extension_namespace['lift_x'](E, x, all=bool(all_points), extend=bool(extend))
        if not all_points: points = [points]
        value = []
        for P in points:
            M = P.curve().base_ring()
            value.append([str(P[0]), str(P[1]), str(M), list(map(str, P.curve().ainvs())),
                          P.curve() is E, bool(P in P.curve()), all(c.parent() is M for c in P),
                          list(map(str, 2*P)), list(map(str, -P))])
        return json.dumps({'value': value}, separators=(',', ':'))
    except Exception as error:
        return json.dumps({'error': type(error).__name__, 'message': str(error)}, separators=(',', ':'))
