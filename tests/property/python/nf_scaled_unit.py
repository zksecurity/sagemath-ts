import json
from sage.all import *
def nf_scaled_unit(op,cs,num,den,exponent):
    R=PolynomialRing(QQ,'x');K=NumberField(R([QQ(c)*num/den for c in cs]),'a');a=K.gen()
    try:
        U=K.unit_group()
        if op==0:
            if K.signature()[0]:
                f=K.polynomial().monic();b,c=f[1],f[0];D=K.discriminant();u=pari(D).quadunit();uv=[QQ(u.component(2)),QQ(u.component(3))];s=(2*a+b)/sqrt((b*b-4*c)/D);w=(1+s)/2 if D%4==1 else s/2;r=[str(c)for c in (uv[0]+uv[1]*w).list()]
            else:r=sorted([','.join(map(str,z.list()))for z in K.roots_of_unity()])
        elif op==1:
            v=[ZZ(1)]+([ZZ(exponent)]if K.signature()[0] else[]);r=list(map(str,U.log(U.exp(v))))
        elif op==2:r=str(int(round(float(K.regulator())*10**8)))
        else:r=str(U.zeta_order())
        error=None
    except Exception as e:r=None;error=type(e).__name__+': '+str(e)
    return json.dumps([r,error],separators=(',',':'))
