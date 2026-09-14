"""Execute bundled division_points with the bundled distinct-root algorithm."""
from pathlib import Path
import textwrap,json
from sage.all import Integer,GF,EllipticCurve,oo
from sage.groups import generic
from finite_polynomial_roots import _roots_univariate_polynomial

_path=Path(__file__).resolve().parents[3]/'reference/sage/src/sage/schemes/elliptic_curves/ell_point.py'
_source=_path.read_text();_start=_source.index('    def division_points(');_end=_source.index('    def _divide_out(',_start)
def _distinct_roots(g):
    return _roots_univariate_polynomial(g.base_ring(),g,None,False)
exec(compile(textwrap.dedent(_source[_start:_end]).replace('g.roots(multiplicities=False)','_distinct_roots(g)'),str(_path),'exec'))
_bundled_division_points=division_points

def elliptic_division_points_order(p,coefficients,point,m):
    try:
        K=GF(p);E=EllipticCurve(K,coefficients);P=E(point) if point else E(0)
        values=_bundled_division_points(P,m)
        return json.dumps({'value':[[[str(c) for c in Q],all(c.parent() is K for c in Q)] for Q in values]},separators=(',',':'))
    except Exception as e:
        return json.dumps({'error':type(e).__name__,'message':str(e)},separators=(',',':'))
