import json
from sage.all import ZZ
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.arith.power import generic_power

def power_case(value,exponent,kind):
    trace=[]
    class Ring(Parent):
        def one(self):trace.append(['one']);return Box(1)
    ring=Ring()
    class Box(Element):
        def __init__(self,v):Element.__init__(self,ring);self.v=int(v)%101
        def __mul__(self,b):trace.append(['mul',self.v,b.v]);return Box(self.v*b.v)
        def __invert__(self):
            trace.append(['inv',self.v])
            if not self.v:raise ZeroDivisionError('zero inverse')
            return Box(pow(self.v,-1,101))
    a=Box(value)
    try:
        r=generic_power(a,int(exponent) if kind==0 else ZZ(exponent));result=[str(r.v),r is a];error=None
    except Exception as e:result=None;error=type(e).__name__+': '+str(e)
    return json.dumps([result,error,trace],separators=(',',':'))
