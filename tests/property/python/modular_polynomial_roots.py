"""Bundled residue-ring roots, including native CRT/Hensel iteration order.

Only method routing is adapted: recursive polynomials keep using the bundled
bodies, and finite-field calls use the existing bundled finite-root oracle.
"""
from pathlib import Path
import textwrap
from sage.all import Zmod,PolynomialRing
from sage.arith.misc import CRT_basis
from sage.misc.mrange import cartesian_product_iterator
from sage.rings.finite_rings.finite_field_base import FiniteField
from polynomial_roots import roots as _generic_roots
from finite_polynomial_roots import _roots_univariate_polynomial as _finite_roots
from polynomial_factor_dispatch import _bundled_word_factor

_path=Path(__file__).resolve().parents[3]/'reference/sage/src/sage/rings/finite_rings/integer_mod_ring.py'
_source=_path.read_text()
_start=_source.index('    def _lift_residue_field_root(')
_end=_source.index('    def _roots_univariate_polynomial(',_start)
exec(compile(textwrap.dedent(_source[_start:_end]),str(_path),'exec'))
_start=_end;_end=_source.index('    #######################################################',_start)
exec(compile(textwrap.dedent(_source[_start:_end]),str(_path),'exec'))
_bundled_modular_roots=_roots_univariate_polynomial

class _RingProxy:
    def __init__(self,ring):self.ring=ring
    def __getattr__(self,name):return getattr(self.ring,name)
    def __call__(self,x):return self.ring(x)
    def _roots_univariate_polynomial(self,f,ring=None,multiplicities=True,algorithm=None):
        return _bundled_modular_roots(self,f,ring,multiplicities,algorithm)
    _lift_residue_field_root=staticmethod(_lift_residue_field_root)

class _ParentProxy:
    def __init__(self,parent):self.parent=parent;self.base=_RingProxy(parent.base_ring())
    def base_ring(self):return self.base
    def __getattr__(self,name):return getattr(self.parent,name)

class _PolynomialProxy:
    def __init__(self,f):self.f=f;self._parent=_ParentProxy(f.parent())
    def parent(self):return self._parent
    def __getattr__(self,name):return getattr(self.f,name)
    def __call__(self,x):return self.f(x)
    def change_ring(self,ring):return _PolynomialProxy(self.f.change_ring(ring))
    def factor(self):
        K=self.f.base_ring();p=K.characteristic()
        return _bundled_word_factor(self.f) if 2<p<2**63 else self.f.factor()
    def roots(self,ring=None,multiplicities=True,algorithm=None):
        if ring is not None and ring is not self.f.base_ring():
            return self.change_ring(ring).roots(multiplicities=multiplicities,algorithm=algorithm)
        K=self.f.base_ring()
        if not isinstance(K,FiniteField):
            return _generic_roots(self,multiplicities=multiplicities,algorithm=algorithm)
        return _finite_roots(K,self if multiplicities else self.f,None,multiplicities)

def mi_polynomial_roots(n,coefficients,multiplicities,hook=-1):
    K=Zmod(n);f=_PolynomialProxy(PolynomialRing(K,'x')(coefficients))
    if hook < 0:
        values=f.roots(multiplicities=bool(multiplicities))
    else:
        base=f.parent().base_ring()
        ring=None if hook==0 else base if hook==1 else Zmod(n+1)
        values=_bundled_modular_roots(base,f,ring,bool(multiplicities),'ignored')
    if multiplicities:
        return [[str(r),str(m),str(r.parent()),r.parent() is K.field()] for r,m in values]
    return [[str(r),str(r.parent()),r.parent() is K] for r in values]

def mi_field(n,direct,values):
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
    K=Zmod(n)
    if direct:
        # A raw Sage generic parent has no factory metadata; on prime orders its
        # first category refinement fails. Initialize the fresh comparison
        # parent with the factory's category, matching the port's category-free
        # construction profile rather than reproducing that metadata failure.
        K.is_field()
        K=IntegerModRing_generic(n,category=K.category())
    F=K.field()
    return [str(F),F.order(),F.characteristic(),F.degree(),F is K.field(),
            [[str(F(v)),F(v).parent() is F] for v in values]]

def mi_factored_order(n,direct):
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
    K=IntegerModRing_generic(n) if direct else Zmod(n)
    factors=K.factored_order()
    return [list(factors),factors is K.factored_order()]

def mi_residue_root_lift(p,e,coefficients,root):
    K=Zmod(p**e);f=PolynomialRing(K,'x')(coefficients);r=Zmod(p)(root)
    values=_lift_residue_field_root(p,e,f,f.derivative(),r)
    return [[str(v),str(v.parent()),v.parent() is K] for v in values]
