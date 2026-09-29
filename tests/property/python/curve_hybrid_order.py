"""Execute bundled point-order bodies and observe their real dependencies."""
import ast
import copy
import json
from pathlib import Path
from types import SimpleNamespace
from sage.all import EllipticCurve, GF, Integer
from sage.rings.factorint import factor_trial_division
from areas.groups_modn import _bundled_group_source
from curve_pari_model import pari_elliptic_model

_root=Path(__file__).resolve().parents[3]
_source=_root/'reference/sage/src/sage/schemes/elliptic_curves/ell_point.py'
_tree=ast.parse(_source.read_text())
_base=next(c for c in _tree.body if isinstance(c,ast.ClassDef) and c.name=='EllipticCurvePoint_field')
_finite=next(c for c in _tree.body if isinstance(c,ast.ClassDef) and c.name=='EllipticCurvePoint_finite_field')
_factor_source=_root/'reference/sage/src/sage/structure/factorization.py'
_factor_tree=ast.parse(_factor_source.read_text())
_complete=next(m for c in _factor_tree.body if isinstance(c,ast.ClassDef) for m in c.body if isinstance(m,ast.FunctionDef) and m.name=='is_complete_factorization')
_ns={}
exec(compile(ast.Module(body=[_complete],type_ignores=[]),str(_factor_source),'exec'),_ns)
class Factor:
    is_complete_factorization=_ns['is_complete_factorization']
    def __init__(self,F): self.__dict__['__x']=list(F)


def ec_factor_limit(n,limit,entry):
    try:
        F=Integer(n).factor(limit=Integer(limit)) if entry else factor_trial_division(Integer(n),Integer(limit))
        factors=([['-1','1']] if F.unit()==-1 else [])+[[str(p),str(e)] for p,e in F]
        return json.dumps(dict(value=factors,complete=Factor(F).is_complete_factorization()),separators=(',',':'))
    except Exception as e: return json.dumps(dict(error=type(e).__name__,message=str(e)),separators=(',',':'))


def ec_hybrid_order(p,a,n,algorithm,cache,alias=0):
    E=EllipticCurve(GF(p),[a,1]); raw=n*E(0,1); calls=[]
    cs=list(map(int,E.a_invariants()))
    class Int:
        def __init__(self,n): self.value=Integer(int(n))
        def __str__(self): return str(self.value)
        def __int__(self): return int(self.value)
        def factor(self,limit):
            calls.append(['factor',str(self),str(limit)])
            return Factor(self.value.factor(limit=limit))
    class Backend:
        def ellcard(self):
            calls.append(['ellcard'])
            return int(json.loads(pari_elliptic_model(p,cs,[],[],0,4)))
        def ellorder(self,Q,N):
            calls.append(['ellorder',str(Q.raw),str(N)])
            target=[] if Q.raw.is_zero() else list(map(int,Q.raw[:2]))
            return int(json.loads(pari_elliptic_model(p,cs,target,[],int(N),5)))
    backend=Backend()
    class Curve:
        def pari_curve(self): return backend
        def order(self):
            calls.append(['curve.order'])
            if not hasattr(self,'_order'): self._order=Int(backend.ellcard())
            return self._order
    curve=Curve()
    generic=_bundled_group_source()
    previous_parent=generic.parent
    def bounds(Q,limits):
        calls.append(['bounds',None if limits is None else list(map(str,limits))])
        return generic.order_from_bounds(Q.raw,limits)
    generic.parent=lambda x: SimpleNamespace(zero=lambda:x.curve()(0)) if hasattr(x,'curve') else previous_parent(x)
    namespace=dict(Integer=Int,generic=SimpleNamespace(order_from_bounds=bounds))
    base=ast.ClassDef(name='Base',bases=[],keywords=[],body=[copy.deepcopy(m) for m in _base.body if isinstance(m,ast.FunctionDef) and m.name in ('order','_compute_order')],decorator_list=[])
    finite=ast.ClassDef(name='Finite',bases=[ast.Name(id='Base',ctx=ast.Load())],keywords=[],body=[copy.deepcopy(m) for m in _finite.body if isinstance(m,ast.FunctionDef) and m.name=='_compute_order'],decorator_list=[])
    exec(compile(ast.fix_missing_locations(ast.Module(body=[base,finite],type_ignores=[])),str(_source),'exec'),namespace)
    class Point(namespace['Finite']):
        additive_order=namespace['Base'].order
        def __init__(self,raw): self.raw=raw
        def is_zero(self): return self.raw.is_zero()
        def curve(self): return curve
    P=Point(raw)
    if cache==1: P._order=raw.order()
    if cache==2: curve._order=Int(E.cardinality())
    try:
        selected='generic_small' if algorithm else 'hybrid'
        method=P.additive_order if alias else P.order
        value=str(method(selected));repeat=str(method(selected))
        result=dict(value=value,repeat=repeat,point_order=str(P._order),curve_order=str(curve._order) if hasattr(curve,'_order') else None)
    except Exception as e: result=dict(error=type(e).__name__,message=str(e))
    finally: generic.parent=previous_parent
    result['calls']=calls
    return json.dumps(result,separators=(',',':'))


def ec_cardinality_cache(p,a,op):
    E=EllipticCurve(GF(p),[a,1])
    try: del E._order
    except AttributeError: pass
    values=[]; states=[]
    for method in ([E.cardinality_pari,E.cardinality,E.order,E.cardinality_pari] if op==0 else [E.order,E.cardinality_pari,E.cardinality]):
        values.append(str(method()));states.append(str(E._order) if hasattr(E,'_order') else None)
    return json.dumps(dict(values=values,states=states),separators=(',',':'))
