"""Native integer/rational coercion and bundled generic-group oracle."""
import json,operator
from sage.all import *

def integer_rational(fn,n,m,den,kind):
    a=ZZ(n)
    b=QQ(m)/den if kind==0 else int(m) if kind==1 else ZZ(m)
    r=[operator.add,operator.sub,operator.mul,operator.eq,operator.lt,operator.le,operator.gt,operator.ge][int(fn)](a,b)
    return json.dumps([str(r).lower() if isinstance(r,bool) else str(r), 'bool' if isinstance(r,bool) else 'Rational' if r.parent() is QQ else 'Integer'],separators=(',',':'))

def integer_group(g,kind,fn,mode,value,n):
    a=ZZ(value) if kind==0 else QQ(value)
    operation='+' if mode==0 else '*'
    target=n*a if operation=='+' else a**n
    if fn==0:r=g.multiple(a,n,operation=operation)
    elif fn==1:r=g.bsgs(a,target,(0,64),operation=operation)
    elif fn==2:r=bool(g.has_order(a,n,operation=operation))
    elif fn==3:r=g.order_from_multiple(a,n,operation=operation)
    elif fn==4:r=g.discrete_log(target,a,ord=4,operation=operation)
    elif fn==5:r=g.order_from_bounds(a,(1,16),operation=operation)
    else:
        op,id,inv,mul=g._parse_group_def(a.parent(),operation,None,None,None)
        r=[str(id),str(inv(a)),str(g.multiple(a,n,operation=op,identity=id,inverse=inv,op=mul)),bool(a==id)]
        return json.dumps(r,separators=(',',':'))
    return json.dumps(r if isinstance(r,bool) else str(r),separators=(',',':'))
