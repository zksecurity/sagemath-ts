"""Sage scalar coercion and zero-division oracle, with unembedded fields."""
import json
from sage.all import *

def _field(d):
    x=polygen(QQ)
    return NumberField(x**3-2 if d==0 else x**2-d, 'a')

def nf_coercion(d, fn, c, n, den, kind):
    K=_field(d)
    a=(c*K.gen()+2)/3
    q=QQ(n)/den
    scalar=int(n) if kind==0 else ZZ(n) if kind==1 else q if kind==2 else K(q)
    if fn==0: r=[str(v) for v in (a+scalar).list()]
    elif fn==1: r=[str(v) for v in (a-scalar).list()]
    elif fn==2: r=[str(v) for v in (a/scalar).list()]
    elif fn==3: r=bool(a==scalar)
    else: r=[str(v) for v in K(scalar).list()]
    return json.dumps(r,separators=(',',':'))

def nf_zero(d, fn, kind):
    K=_field(d)
    z=K.zero()
    scalar=int(0) if kind==0 else ZZ(0) if kind==1 else QQ.zero() if kind==2 else z
    r=~z if fn==0 else K.gen()/scalar if fn==1 else z/scalar if fn==2 else z**(-1)
    return json.dumps([str(v) for v in r.list()],separators=(',',':'))

def nf_float(d, n, den, scale):
    v=float(n)/float(den) if den else float('inf') if n>0 else -float('inf') if n<0 else float('nan')
    v*=2.0**int(scale)
    return json.dumps([str(x) for x in _field(d)(v).list()],separators=(',',':'))
