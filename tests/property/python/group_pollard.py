"""Bundled Pollard algorithms with parent, mutability and hash-collision controls."""
import json
from sage.all import ZZ,set_random_seed
from sage.misc.prandom import getrandbits

def pollard(g,algorithm,seed,mode,modulus,value,n,order,lb,ub,mutable,collision,fail_parent):
    set_random_seed(ZZ(seed));modulus=int(modulus);trace=[];active=False;additive=mode%2==0
    def check(x):
        if mutable==2 and not x.frozen:raise TypeError('mutable element')
    class Box:
        def __init__(self,v):self.v=int(v)%modulus;self.frozen=False
        @staticmethod
        def zero():
            if active:
                trace.append(['zero'])
                if fail_parent:raise ValueError('identity sentinel')
            return Box(0)
        @staticmethod
        def one():
            if active:
                trace.append(['one'])
                if fail_parent:raise ValueError('identity sentinel')
            return Box(1)
        def __add__(self,b):trace.append(['add',self.v,b.v]);return Box(self.v+b.v)
        def __mul__(self,b):
            if isinstance(b,Box):trace.append(['mul',self.v,b.v]);return Box(self.v*b.v)
            trace.append(['scale',self.v,int(b)]);return Box(self.v*int(b))
        def __rmul__(self,b):return self*b
        def __pow__(self,k):trace.append(['pow',self.v,int(k)]);return Box(pow(self.v,int(k),modulus))
        def __neg__(self):trace.append(['neg',self.v]);return Box(-self.v)
        def __invert__(self):trace.append(['inv',self.v]);return Box(pow(self.v,-1,modulus))
        def __eq__(self,b):return isinstance(b,Box) and self.v==b.v
        def __hash__(self):check(self);return 0 if collision==1 else self.v%4 if collision==2 else self.v%17 if collision==3 else hash(self.v)
        def __str__(self):check(self);return 'collision' if collision==1 else str(self.v%4) if collision==2 else str(self.v%17) if collision==3 else str(self.v)
        def set_immutable(self):
            trace.append(['freeze',self.v]);self.frozen=True
            if mutable==3:raise ValueError('freeze sentinel')
    if not mutable:del Box.set_immutable
    base=Box(value);target=Box(int(value)*int(n) if additive else pow(int(value),int(n),modulus))
    kw={'operation':'+' if additive else '*'}
    if mode>=2:kw={'operation':'other','identity':Box.zero() if additive else Box.one(),'inverse':lambda x:-x if additive else ~x,'op':lambda x,y:x+y if additive else x*y}
    def hs(x):
        check(x)
        if algorithm==0:return int(x.v*x.v+1)
        h=5381
        for c in str(x):h=((h*33)^ord(c))&0xffffffff
        return h
    kw['hash_function']=hs;trace.clear();active=True
    try:
        r=str(g.discrete_log_lambda(target,base,(ZZ(lb),ZZ(ub)),**kw) if algorithm==0 else g.discrete_log_rho(target,base,ZZ(order),**kw));error=None
    except Exception as e:r=None;error=type(e).__name__+': '+str(e)
    return json.dumps([r,error,trace,base.frozen,target.frozen,str(getrandbits(64))],separators=(',',':'))

from sage.all import matrix,GF,identity_matrix,zero_matrix

def pollard_matrix(g,seed,mode,n,lb,ub):
    set_random_seed(ZZ(seed));base=matrix(GF(2),[[0,1],[1,1]]);identity=zero_matrix(GF(2),2) if mode==0 else identity_matrix(GF(2),2)
    op=lambda x,y:x+y if mode==0 else x*y
    inverse=lambda x:-x if mode==0 else ~x
    target=identity
    for _ in range(int(n)):target=op(target,base)
    def hs(x):
        v=0
        for i in range(2):
            for j in range(2):v=2*v+int(x[i,j])
        return v
    try:r=str(g.discrete_log_lambda(target,base,(ZZ(lb),ZZ(ub)),operation='other',identity=identity,inverse=inverse,op=op,hash_function=hs));error=None
    except Exception as e:r=None;error=type(e).__name__+': '+str(e)
    return json.dumps([r,error,base.is_immutable(),target.is_immutable(),str(getrandbits(64))],separators=(',',':'))
