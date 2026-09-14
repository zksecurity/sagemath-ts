"""Native double arithmetic and bundled generic.py operation schedules."""
import struct,math,operator
from sage.all import *

def _from_bits(b):return struct.unpack('>d',struct.pack('>Q',int(b)))[0]
def _bits(x):
    x=float(x)
    return 'NaN' if math.isnan(x) else str(struct.unpack('>Q',struct.pack('>d',x))[0])
def rdf_scalar(fn,left,n,d,right,kind):
    a=RDF(_from_bits(left))
    b=int(n) if kind==0 else ZZ(n) if kind==1 else QQ(n)/d if kind==2 else _from_bits(right) if kind==3 else RDF(_from_bits(right)) if kind==4 else bool(n)
    r=[operator.add,operator.sub,operator.mul,operator.truediv,operator.eq][int(fn)](a,b)
    return str(r).lower() if isinstance(r,bool) else _bits(r)
def rdf_multiple(g,value,n,mode):return _bits(g.multiple(RDF(_from_bits(value)),n,operation='+' if mode==0 else '*'))

def rdf_alias(fn,left,n,kind):
    import json
    a=RDF(_from_bits(left))
    b=int(n) if kind==0 else ZZ(n) if kind==1 else QQ(n) if kind==2 else float(n) if kind==3 else RDF(n) if kind==4 else bool(n)
    r=[operator.add,operator.sub,operator.mul,operator.truediv][int(fn)](a,b)
    return json.dumps([_bits(r),r is a],separators=(',',':'))
