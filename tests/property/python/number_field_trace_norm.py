"""Native no-argument trace/norm and the port's existing absolute-value aliases."""
from sage.all import *
def nf_trace_norm(op,cs,flat,d):
    R=PolynomialRing(QQ,'x');K=NumberField(R(cs),'a',check=False);a=K([QQ(c)/d for c in flat])
    return str(a.trace()if op in (0,3,5)else a.norm())
