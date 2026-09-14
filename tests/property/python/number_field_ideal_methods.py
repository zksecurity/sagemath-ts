"""Sage coercion at ideal arithmetic and element-membership boundaries."""
import json
from sage.all import *
from number_field_ideal_inputs import _field, _lattice
from sage.rings.number_field.number_field_ideal import NumberFieldIdeal

def nf_ideal_method_input(cs,ideal_kind,op,kind,n,d,other_cs=None):
    K=_field(cs);L=_field(cs if other_cs is None else other_cs,'b');q=QQ(n)/d
    I=[K.ideal(0),K.ideal(1),K.ideal(2),K.ideal(QQ(1)/2),K.ideal(K.gen()),K.ideal(3,K.gen()+1)][int(ideal_kind)]
    values=[int(n),ZZ(n),q,float(q),K(q),K.gen()+q,[q],[],[q]*K.degree(),[q,K.gen()],K.ideal(q),K.ideal(0),L(q),L.gen()+q,L.ideal(q),L.ideal(L.gen()+q),[0],K.ideal(3,K.gen()+1),NumberFieldIdeal(K,[K(q)])]
    x=values[int(kind)]
    try:
        r=(x in I)if op==0 else I.divides(x)if op==1 else I.is_coprime(x)if op==2 else _lattice(I+x)if op==3 else _lattice(I*x)if op==4 else _lattice(I/x)
        error=None
    except Exception as e:r=None;error=type(e).__name__+': '+str(e)
    return json.dumps([r,error],separators=(',',':'))

def nf_field_vector(cs,kind,length,n,d):
    K=_field(cs);L=_field(cs,'b')
    try:
        value=int(n)if kind==0 else ZZ(n)if kind==1 else QQ(n)/d if kind==2 else float(n)/float(d)if kind==3 and d else float('inf')if kind==3 and n>0 else -float('inf')if kind==3 and n<0 else float('nan')if kind==3 else K(QQ(n)/d)if kind==4 else L(QQ(n)/d)if kind==5 else K.gen()if kind==6 else L.gen()
        r=[str(c)for c in K([value]*int(length)).list()];error=None
    except Exception as e:r=None;error=type(e).__name__+': '+str(e)
    return json.dumps([r,error],separators=(',',':'))

def nf_basis_class_number(cs,fn,fd,op):
    R=PolynomialRing(QQ,'x');K=NumberField(R([QQ(c)*fn/fd for c in cs]),'a')
    try:r=str(K.class_number()if op==0 else K.class_group().order());error=None
    except Exception as e:r=None;error=type(e).__name__+': '+str(e)
    return json.dumps([r,error],separators=(',',':'))
