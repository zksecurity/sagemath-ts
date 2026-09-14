"""Execute the bundled Sage roots methods using native polynomial objects.

Installed Sage 10.3 propagates ArithmeticError on QQ(0). The bundled generic
method translates it to NotImplementedError. Load the actual bundled bodies;
expose Cython's private parent field through its public parent() accessor and
remove only declarations (the sparse integer algorithm is otherwise unchanged).
"""
from pathlib import Path
import re
import textwrap
import sage
from sage.all import QQ, ZZ, PolynomialRing
from sage.rings.integer_ring import IntegerRing_class
from sage.rings.rational_field import RationalField

_root = Path(__file__).resolve().parents[3] / 'reference/sage/src/sage/rings'
_source = (_root / 'polynomial/polynomial_element.pyx').read_text()
_start = _source.index('    def roots(self, ring=None, multiplicities=True, algorithm=None, **kwds):')
_end = _source.index('    def _roots_from_factorization', _start)
exec(compile(textwrap.dedent(_source[_start:_end]).replace('self._parent', 'self.parent()'), str(_root / 'polynomial/polynomial_element.pyx'), 'exec'))
_source = (_root / 'integer_ring.pyx').read_text()
_start = _source.index('    def _roots_univariate_polynomial(')
_end = _source.index('    #################################', _start)
def _declaration(match):
    declaration = match[0].strip()[5:]
    if '=' not in declaration:
        return ''
    left, right = declaration.split('=', 1)
    name = left.split(',')[-1].split()[-1]
    return '        ' + name + ' = ' + right.strip()
_body = re.sub(r'^        cdef .*$', _declaration, _source[_start:_end], flags=re.M)
exec(compile(textwrap.dedent(_body), str(_root / 'integer_ring.pyx'), 'exec'))


def polynomial_roots_fidelity(kind, coefficients, denominator):
    R = PolynomialRing(ZZ if kind == 0 else QQ, 'x')
    f = R([QQ(c)/denominator for c in coefficients])
    result = _roots_univariate_polynomial(ZZ, f) if kind == 0 else roots(f)
    return [[str(r), str(m), str(r.parent())] for r, m in result]


def polynomial_distinct_roots_fidelity(kind, coefficients, denominator):
    R = PolynomialRing(ZZ if kind == 0 else QQ, 'x')
    f = R([QQ(c)/denominator for c in coefficients])
    result = _roots_univariate_polynomial(ZZ, f, multiplicities=False) if kind == 0 else roots(f, multiplicities=False)
    return [[str(r), str(r.parent())] for r in result]
