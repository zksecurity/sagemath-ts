"""Sage trace/norm base arguments, return parents and error contracts."""
from sage.all import *
import json
def nf_base_trace_norm(op,n,cs,d,mode):
    n=int(n);R=PolynomialRing(QQ,'x');x=R.gen();K=NumberField(x**n-2,'a',check=False);B=NumberField(x-3,'b',check=False)
    a=K([QQ(c)/d for c in cs])
    args=[]if mode==0 else[None]if mode==1 else[QQ]if mode==2 else[K]if mode==3 else[B]if mode==4 else[NumberField(x**7-3,'b',check=False)]if mode==5 else[QQ,QQ]if mode==6 else[4]if mode==7 else[0]if mode==8 else['bad']if mode==9 else[True]if mode==10 else[1.5]if mode==11 else[[]]if mode==12 else[{}]if mode==13 else[Integer(4)]if mode==14 else[QQ(1)/2]
    z=(a.trace if op==0 else a.norm)(*args)
    return json.dumps({'parent':'QQ'if z.parent()is QQ else'K'if z.parent()is K else'B','coefficients':[str(z)]if z.parent()is QQ else[str(c)for c in z.list()]},separators=(',',':'))
