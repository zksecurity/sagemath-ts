"""Bundled p-primary caller with observed division and group-relation bindings."""
import ast
import copy
import json
from pathlib import Path
from types import SimpleNamespace
from sage.all import GF, PolynomialRing, EllipticCurve, Integer, QQ, oo
from areas.groups_modn import _bundled_group_source
generic = _bundled_group_source()

_source = Path(__file__).resolve().parents[3] / 'reference/sage/src/sage/schemes/elliptic_curves/ell_generic.py'
_tree = ast.parse(_source.read_text())
_method = next(n for c in _tree.body if isinstance(c, ast.ClassDef) for n in c.body
               if isinstance(n, ast.FunctionDef) and n.name == '_p_primary_torsion_basis')

def curve(p, T, coefficients):
    if T:
        K = GF(p**(len(T)-1), 'a', modulus=PolynomialRing(GF(p),'t')(T), impl='pari_ffelt')
        def decode(n):
            ds=[]
            for _ in range(len(T)-1): ds.append(n % p); n //= p
            return K(ds)
    else:
        K=GF(p) if p else QQ
        decode=K
    return EllipticCurve(K, list(map(decode,coefficients)))

def ec_primary_torsion(p, T, coefficients, prime, bound):
    calls=[]
    original_parent=generic.parent
    # Sage 10.3 point homsets lack zero(); 10.9's generic parser expects it.
    # Supply only this parent identity binding, retaining the bundled algorithms.
    generic.parent=lambda x: SimpleNamespace(zero=lambda:x.curve()(0)) if hasattr(x,'curve') else original_parent(x)
    try:
        E=curve(p,T,coefficients)
        class Point:
            def __init__(self, point): self.point=point
            def is_zero(self): return self.point.is_zero()
            def division_points(self, n):
                calls.append(['division_points',str(self.point),str(n)])
                return list(map(Point,self.point.division_points(n)))
            def __add__(self, other): return Point(self.point+other.point)
        def relation(P,Q,operation):
            calls.append(['linear_relation',str(P.point),str(Q.point)])
            return generic.linear_relation(P.point,Q.point,operation)
        def multiples(P,n,start,operation):
            return map(Point,generic.multiples(P.point,n,start.point,operation=operation))
        namespace=dict(Integer=Integer,generic=SimpleNamespace(linear_relation=relation,multiples=multiples))
        exec(compile(ast.Module(body=[_method],type_ignores=[]),str(_source),'exec'),namespace)
        basis=namespace['_p_primary_torsion_basis'](lambda n:Point(E(n)),prime,None if bound==-99 else bound)
        result={'value':[[str(P.point),int(k)] for P,k in basis]}
    except Exception as e: result={'error':type(e).__name__,'message':str(e)}
    finally: generic.parent=original_parent
    result['calls']=calls
    return json.dumps(result,separators=(',',':'))


_point_source = _source.with_name('ell_point.py')
_point_tree = ast.parse(_point_source.read_text())
_division_method = next(n for c in _point_tree.body if isinstance(c,ast.ClassDef) for n in c.body
                        if isinstance(n,ast.FunctionDef) and n.name=='division_points')

class _DivisionScalarBindings(ast.NodeTransformer):
    def visit_Call(self,node):
        self.generic_visit(node)
        if isinstance(node.func,ast.Attribute) and node.func.attr=='sqrt':
            node.args.insert(0,node.func.value)
            node.func=ast.Name(id='_scalar_sqrt',ctx=ast.Load())
        return node


def ec_division_points(p,T,coefficients,target,m,poly_only,known_order):
    original_parent=generic.parent
    generic.parent=lambda x: SimpleNamespace(zero=lambda:x.curve()(0)) if hasattr(x,'curve') else original_parent(x)
    try:
        E=curve(p,T,coefficients)
        K=E.base_ring()
        def decode(n):
            if not T: return K(n)
            ds=[]
            for _ in range(len(T)-1): ds.append(n%p); n//=p
            return K(ds)
        P=E(list(map(decode,target))) if target else E(0)
        if known_order: P._order=P.order()
        else:
            try: del P._order
            except AttributeError: pass
        def scalar_sqrt(value):
            if not T: return value.sqrt()
            # Bind the bundled PARI version, whose representative can differ from
            # the installed Sage library's square root; preserve caller arithmetic.
            from pari_ff_square_root import pari_ff_square_root
            raw=json.loads(pari_ff_square_root(p,T,list(map(int,value.polynomial().list())),1))
            if raw['value'] is None: raise ValueError('element is not a square')
            result=K(list(map(int,raw['value'])))
            assert result*result==value
            return result
        namespace=dict(Integer=Integer,generic=generic,oo=oo,_scalar_sqrt=scalar_sqrt)
        method=_DivisionScalarBindings().visit(copy.deepcopy(_division_method))
        exec(compile(ast.fix_missing_locations(ast.Module(body=[method],type_ignores=[])),str(_point_source),'exec'),namespace)
        result=namespace['division_points'](P,m,bool(poly_only))
        value={'points':[[str(Q),str(Q._order) if hasattr(Q,'_order') else None]for Q in result]} if isinstance(result,list) else {'polynomial':list(map(str,result.list()))}
        return json.dumps({'value':value,'target_order':str(P._order) if hasattr(P,'_order') else None},separators=(',',':'))
    except Exception as e: return json.dumps({'error':type(e).__name__,'message':str(e)},separators=(',',':'))
    finally: generic.parent=original_parent


def ec_scalar_order(modulus, m, known):
    from sage.all import Zmod, ZZ, pari, Sequence
    from sage.structure.coerce_actions import IntegerMulAction
    from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_finite_field
    from cypari2.handle_error import PariError
    # Installed Sage 10.3 selects finite-field points even over composite rings.
    # Select the scalar body by the bundled constructor's prime/field dispatch.
    name = 'EllipticCurvePoint_finite_field' if Integer(modulus).is_prime() else 'EllipticCurvePoint'
    method = next(n for c in _point_tree.body if isinstance(c,ast.ClassDef) and c.name==name
                  for n in c.body if isinstance(n,ast.FunctionDef) and n.name=='_acted_upon_')
    namespace = dict(ZZ=ZZ, pari=pari, Sequence=Sequence, PariError=PariError,
                     IntegerMulAction=IntegerMulAction,
                     EllipticCurvePoint_finite_field=EllipticCurvePoint_finite_field)
    exec(compile(ast.Module(body=[method],type_ignores=[]),str(_point_source),'exec'),namespace)
    E=EllipticCurve(Zmod(modulus),[1,0])
    P=E(0,0)
    if known: P._order=Integer(2)
    Q=namespace['_acted_upon_'](P,Integer(m),False)
    return json.dumps([str(Q),str(Q._order) if hasattr(Q,'_order') else None],separators=(',',':'))
