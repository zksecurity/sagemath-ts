"""Native number-field element and defining-polynomial behavior."""
import json
from sage.all import *

def nf_element(fn,fcoeffs,fden,coeffs,den,exponent,kind):
    R=PolynomialRing(QQ,'x');K=NumberField(R([QQ(c)/fden for c in fcoeffs]),'a');a=K(R([QQ(c)/den for c in coeffs]))
    scalar=lambda x:str(x)
    poly=lambda x:[str(c)for c in x.list()]
    elem=lambda x:[str(c)for c in x.list()]
    try:
        if fn==0:r=elem(a)
        elif fn==1:r=elem(K.gen())
        elif fn==2:r=poly(a.polynomial())
        elif fn==3:r=poly(K.defining_polynomial())
        elif fn==4:r=poly(K.polynomial())
        elif fn==5:r=str(K)
        elif fn==6:r=str(a)
        elif fn==7:r=scalar(a.norm())
        elif fn==8:r=scalar(a.trace())
        elif fn==9:r=poly(a.charpoly())
        elif fn==10:r=poly(a.minpoly())
        elif fn==11:r=bool(a.is_integral())
        elif fn==12:r=bool(a.is_unit())
        elif fn==13:r=bool(a.is_integral() and abs(a.norm())==1)
        elif fn==14:r=str(a.denominator())
        elif fn==15:r=elem(a.numerator() if K.degree()==2 else a*a.denominator()) # numerator extension to all absolute fields
        elif fn==16:
            result=a**(int(exponent) if kind==0 else ZZ(exponent));r=[elem(result),result is a]
        elif fn==17:r=scalar(a.absolute_norm())
        elif fn==18:r=scalar(a.relative_norm())
        elif fn in [19,20]:r=scalar(a.trace()) # port-only trace aliases
        elif fn==21:r=a.parent() is K
        elif fn==22:r=bool(a.is_zero())
        elif fn==23:r=bool(a.is_one())
        else:r=poly(a.minpoly()) # port-only minimal_polynomial alias
        error=None
    except Exception as e:r=None;error=type(e).__name__+': '+str(e)
    return json.dumps([r,error],separators=(',',':'))
