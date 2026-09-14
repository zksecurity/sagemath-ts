"""Sage's integer-factor degree dispatch and cached irreducibility observations."""
import json
from sage.all import QQ, PolynomialRing, NumberField, pari

def nf_pari_factor_constructor(n,kind,scale):
    n,kind,scale=int(n),int(kind),int(scale)
    R=PolynomialRing(QQ,'x');x=R.gen()
    f=5*x**n-2 if kind==5 else x**n-2 if kind==0 else x**n-1 if kind==1 else (x-1)*(x**(n-1)-2)if kind==2 else x**(n-2)*(x-1)**2 if kind==3 else x**n+x+1
    f*=QQ(scale)/7;result=[bool(f.is_irreducible())]
    try:result.append(['degree',str(NumberField(f,'a').degree())])
    except Exception as e:result.append([type(e).__name__,str(e)])
    return json.dumps(result,separators=(',',':'))

def nf_irreducible_cache(P):
    f=PolynomialRing(QQ,'x')(P);saved=pari.getrand()
    try:
        pari.setrand(1);first=bool(f.is_irreducible());pari.setrand(2);before=pari.getrand();second=bool(f.is_irreducible())
        return json.dumps([first,second,bool(before==pari.getrand())],separators=(',',':'))
    finally:pari.setrand(saved)
