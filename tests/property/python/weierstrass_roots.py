"""Bundled isomorphism algorithm with bundled finite-field root dispatch.

Polynomials use native Sage arithmetic; only the roots method is adapted to
execute the bundled version. Extension parents explicitly select pari_ffelt.
"""
from pathlib import Path
import json
from sage.all import GF, QQ, PolynomialRing, EllipticCurve
from finite_polynomial_roots import _roots_univariate_polynomial
from sage.schemes.elliptic_curves.weierstrass_morphism import baseWI

class _Polynomial:
    def __init__(self, f): self.f = f
    def __pow__(self, n): return _Polynomial(self.f**n)
    def __add__(self, x): return _Polynomial(self.f + (x.f if isinstance(x, _Polynomial) else x))
    __radd__ = __add__
    def __sub__(self, x): return _Polynomial(self.f - (x.f if isinstance(x, _Polynomial) else x))
    def __mul__(self, x): return _Polynomial(self.f * (x.f if isinstance(x, _Polynomial) else x))
    __rmul__ = __mul__
    def roots(self, multiplicities=False):
        if self.f.base_ring() is QQ:
            from polynomial_roots import roots as bundled_roots
            return bundled_roots(self.f, multiplicities=multiplicities)
        return _roots_univariate_polynomial(self.f.base_ring(), self.f, None, multiplicities)

def _polygen(K, name): return _Polynomial(PolynomialRing(K, name).gen())
_path = Path(__file__).resolve().parents[3] / 'reference/sage/src/sage/schemes/elliptic_curves/weierstrass_morphism.py'
_source = _path.read_text()
_source = _source[_source.index('def _isomorphisms('):_source.index('\n\nclass WeierstrassIsomorphism(')]
# Resolve the relative import while executing the original function outside its module.
_source = _source.replace('from .ell_generic import', 'from sage.schemes.elliptic_curves.ell_generic import')
_source = _source.replace('from sage.rings.polynomial.polynomial_ring import polygen', 'polygen = _polygen')
exec(compile(_source, str(_path), 'exec'))

def wm_polynomial_root_isomorphisms(p, degree, modulus, left, right, transform):
    try:
        if degree > 1:
            K = GF(p**degree, 'a', modulus=PolynomialRing(GF(p), 't')(modulus), impl='pari_ffelt')
            def decode(v):
                digits=[]
                for _ in range(degree): digits.append(v % p);v //= p
                return K(digits)
            encode = lambda v: str(v.integer_representation())
        else:
            K = GF(p) if p else QQ
            decode = K
            encode = lambda v: str(v)
        E = EllipticCurve(K, list(map(decode, left)))
        if transform:
            F = EllipticCurve(K, baseWI(*map(decode, transform))(E.ainvs()))
        else:
            F = EllipticCurve(K, list(map(decode, right)))
        rows = list(_isomorphisms(E, F))
        value = [[[encode(v) for v in row], all(v.parent() is K for v in row),
                  list(baseWI(*row)(E.ainvs())) == list(F.ainvs())] for row in rows]
        import sage.schemes.elliptic_curves.weierstrass_morphism as module
        original = module._isomorphisms
        try:
            module._isomorphisms = _isomorphisms
            try:
                iso = module.WeierstrassIsomorphism(E, None, F)
                first = {'value': [encode(v) for v in iso.tuple()]}
            except Exception as e: first = {'error': type(e).__name__, 'message': str(e)}
        finally: module._isomorphisms = original
        return json.dumps({'value': [value, first] },separators=(',', ':'))
    except Exception as e:
        return json.dumps({'error':type(e).__name__,'message':str(e)},separators=(',', ':'))

def wm_isomorphism_argument_errors(left, right):
    try:
        E = EllipticCurve(GF(5), [0, 1])
        values = [E, None, 7, [], object()]
        return json.dumps({'value': len(list(_isomorphisms(values[left], values[right])))}, separators=(',', ':'))
    except Exception as e:
        return json.dumps({'error': type(e).__name__, 'message': str(e)}, separators=(',', ':'))
