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


# Sage 10.3 lacks reduced_basis; execute the bundled 10.9 method verbatim.
# Compile only the two original methods, using installed Sage arithmetic objects.
import ast
from pathlib import Path
_source=Path(__file__).resolve().parents[4]/'reference/sage/src/sage/algebras/quatalg/quaternion_algebra.py'
_class=next(n for n in ast.parse(_source.read_text()).body if isinstance(n,ast.ClassDef) and n.name=='QuaternionFractionalIdeal_rational')
_original_reductions={}
for _method in ('reduced_basis','minimal_element'):
    _node=next(n for n in _class.body if isinstance(n,ast.FunctionDef) and n.name==_method)
    _namespace=dict(globals())
    exec(compile(ast.Module(body=[_node],type_ignores=[]),str(_source),'exec'),_namespace)
    _original_reductions[_method]=_namespace[_method]

def _ideal_proxy(native, reverse_basis=False):
    if reverse_basis:
        # Bundled 10.9 product/scale uses basis_for_quaternion_lattice(reverse=True).
        # Installed 10.3 uses the forward HNF; apply the bundled coordinate order.
        A=native.quaternion_algebra()
        G=matrix(QQ,[list(x)[::-1]for x in native.basis()][::-1])
        Z,d=G._clear_denom();H=Z._hnf_pari(0,include_zero_rows=False)
        native=A.ideal([A([QQ(c)/d for c in row][::-1])for row in H.rows()[::-1]],check=False)
    a,b=native.quaternion_algebra().invariants()
    from pari_lll_dependents import pari_lll_dependents
    import json
    class GramProxy:
        def LLL_gram(self):
            G,_=native.gram_matrix()._clear_denom()
            cols=json.loads(pari_lll_dependents(0,4,0,list(G.list())))
            U=matrix(ZZ,cols).transpose()
            if U.det()==-1:U.rescale_col(3,-1)
            return U
    class PariFormProxy:
        def qfminim(self,*args):
            G=native.quadratic_form().matrix()
            norm,vector=json.loads(pari_lll_dependents(6,4,0,[ZZ(c)for c in G.list()]))
            return ZZ(norm),[ZZ(c)for c in vector]
    class FormProxy:
        def __pari__(self):return PariFormProxy()
    class IdealProxy:
        def gram_matrix(self):return GramProxy()
        def quadratic_form(self):return FormProxy()
        def __getattr__(self,name):return getattr(native,name)
        def quaternion_algebra(self):
            class AlgebraProxy:
                def is_definite(self):return a<0 and b<0
            return AlgebraProxy()
    return IdealProxy()


def quat_lattice_reduction(a,b,den,flat,op):
    A=QuaternionAlgebra(QQ,ZZ(a),ZZ(b),names=('i','j','k'))
    native=A.ideal([A([QQ(x)/den for x in flat[i*4:i*4+4]])for i in range(4)])
    I=_ideal_proxy(native)
    if op==0:return ';'.join(_vals(list(x))for x in _original_reductions['reduced_basis'](I))
    if op==1:return _vals(list(_original_reductions['minimal_element'](I)))
    if op==2:return _vals(native.quadratic_form().matrix().list())
    from pari_lll_dependents import pari_lll_dependents
    import json
    G=native.quadratic_form().matrix()
    counts=json.loads(pari_lll_dependents(7,4,1,[ZZ(c)for c in G.list()]+[11]))
    return _vals([1]+[2*ZZ(c)for c in counts])


def quat_order_isomorphism(scale):
    from sage.algebras.quatalg.quaternion_algebra import QuaternionFractionalIdeal_rational,QuaternionOrder
    from sage.rings.rational_field import RationalField
    A=QuaternionAlgebra(QQ,-1,-19,names=('i','j','k'));i,j,k=A.gens()
    O0=A.quaternion_order([A(1),i,(i+j)/2,(1+k)/2])
    O1=A.quaternion_order([A(1),667*i,A(1)/2+j/2+9*i,(222075*i/2+333*j+k/2)/667])
    if scale!=0:
        alpha=A([1,scale,scale+1,2*scale-1])
        O1=A.quaternion_order([~alpha*x*alpha for x in O1.basis()])
    order=next(n for n in ast.parse(_source.read_text()).body if isinstance(n,ast.ClassDef)and n.name=='QuaternionOrder')
    method=next(n for n in order.body if isinstance(n,ast.FunctionDef)and n.name=='isomorphism_to')
    namespace=dict(globals(),QuaternionOrder=QuaternionOrder,RationalField=RationalField)
    exec(compile(ast.Module(body=[method],type_ignores=[]),str(_source),'exec'),namespace)
    old_min=QuaternionFractionalIdeal_rational.minimal_element
    algebra_type=type(A);old_def=algebra_type.__dict__.get('is_definite')
    try:
        algebra_type.is_definite=lambda self:all(x<0 for x in self.invariants())
        QuaternionFractionalIdeal_rational.minimal_element=lambda self:_original_reductions['minimal_element'](_ideal_proxy(self,reverse_basis=True))
        gamma=namespace['isomorphism_to'](O0,O1,conjugator=True)
        return ';'.join(_vals(list(x))for x in [gamma]+[~gamma*x*gamma for x in A.gens()])
    finally:
        QuaternionFractionalIdeal_rational.minimal_element=old_min
        if old_def is None:delattr(algebra_type,'is_definite')
        else:algebra_type.is_definite=old_def


FUNCTIONS = {
    'quat_order_isomorphism': quat_order_isomorphism,
    'quat_lattice_reduction': quat_lattice_reduction,
    'quat_element': quat_element,
    'quat_product': quat_product,
    'quat_algebra': quat_algebra,
}
