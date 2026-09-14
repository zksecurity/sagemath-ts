"""Trace genuine group operations through bundled Sage generic algorithms."""
import json

def group_parent(g, fn, mode, modulus, value, n, lb, ub, d, variant=0):
    modulus=int(modulus); trace=[]; additive=int(mode)%2==0;active=False
    class Box:
        def __init__(self,v):self.v=int(v)%modulus
        @staticmethod
        def zero():
            trace.append(['zero'])
            if active and variant==1:raise ValueError('identity sentinel')
            return Box(0)
        @staticmethod
        def one():
            trace.append(['one'])
            if active and variant==1:raise ValueError('identity sentinel')
            return Box(1)
        def __add__(self,b):
            trace.append(['add',self.v,b.v]);return Box(self.v+b.v)
        def __mul__(self,b):
            if isinstance(b,Box):trace.append(['mul',self.v,b.v]);return Box(self.v*b.v)
            trace.append(['scale',self.v,int(b)]);return Box(self.v*int(b))
        def __rmul__(self,b):return self*b
        def __pow__(self,k):
            trace.append(['pow',self.v,int(k)]);return Box(pow(self.v,int(k),modulus))
        def __neg__(self):trace.append(['neg',self.v]);return Box(-self.v)
        def __invert__(self):trace.append(['inv',self.v]);return Box(pow(self.v,-1,modulus))
        def __eq__(self,b):
            trace.append(['eq',self.v,b.v])
            if active and variant==2:raise TypeError('equality sentinel')
            return isinstance(b,Box) and self.v==b.v
        def __bool__(self):trace.append(['bool',self.v]);return bool(self.v)
        def is_one(self):trace.append(['is_one',self.v]);return self.v==1
        def __hash__(self):return 0 if mode>=4 else hash(self.v)
        def __str__(self):return 'collision' if mode>=4 else str(self.v)
    a=Box(value);target=Box(int(value)*int(n) if additive else pow(int(value),int(n),modulus))
    kw={'operation':'+' if additive else '*'}
    if mode%4>=2:kw={'operation':'other','identity':Box.zero() if additive else Box.one(),'inverse':lambda x:-x if additive else ~x,'op':lambda x,y:x+y if additive else x*y}
    trace.clear();active=True
    try:
        if fn==0:r=g.multiple(a,n,**kw);r=str(r)
        elif fn==1:r=str(g.bsgs(a,target,(lb,ub),**kw))
        elif fn==2:r=bool(g.has_order(a,n,operation=kw['operation']))
        elif fn==3:r=str(g.order_from_bounds(a,(lb,ub),d,**kw))
        elif fn==4:r=str(g.order_from_multiple(a,n,**kw))
        elif fn==6:r=str(g.order_from_bounds(a,None,d,**kw))
        else:r=str(g.discrete_log(target,a,ord=(modulus if additive else modulus-1),**kw))
        error=None
    except Exception as e:r=None;error=type(e).__name__+': '+str(e)
    return json.dumps([r,error,trace],separators=(',',':'))

