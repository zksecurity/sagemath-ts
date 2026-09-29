"""Bundled finite-curve enumeration with live Sage group arithmetic."""
import ast
import json
from pathlib import Path
from types import MethodType
from sage.all import EllipticCurve, GF

_source=Path(__file__).resolve().parents[3]/'reference/sage/src/sage/schemes/elliptic_curves/ell_finite_field.py'
_tree=ast.parse(_source.read_text())
_methods=[n for c in _tree.body if isinstance(c,ast.ClassDef) for n in c.body
          if isinstance(n,ast.FunctionDef) and n.name in ('points','_points_via_group_structure')]
_namespace={}
exec(compile(ast.Module(body=_methods,type_ignores=[]),str(_source),'exec'),_namespace)

def ec_finite_points(p,a,b):
    E=EllipticCurve(GF(p),[a,b])
    calls=[]
    class Curve:
        def __call__(self,x): return E(x)
        def abelian_group(self):
            calls.append('abelian_group')
            return E.abelian_group()
    curve=Curve()
    for name in ('points','_points_via_group_structure'):
        setattr(curve,name,MethodType(_namespace[name],curve))
    points=curve.points()
    same=points is curve.points()
    try:
        points[0]=E(0)
        mutation=None
    except Exception as e: mutation=[type(e).__name__,str(e)]
    group=E.abelian_group()
    value=[str(P) for P in points]
    assert value==[str(P) for P in E.points()]
    return json.dumps(dict(value=value,calls=calls,same=same,mutation=mutation,
                          group_cached=group is E.abelian_group(),
                          gens_updated=list(E.gens())==[g.element() for g in group.gens()]),separators=(',',':'))
