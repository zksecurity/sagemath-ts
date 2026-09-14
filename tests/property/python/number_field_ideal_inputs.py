"""Ideal input/class/zero semantics, compared in canonical power-basis coordinates."""
import json
from sage.all import *
from sage.rings.number_field.number_field_ideal import NumberFieldIdeal, NumberFieldFractionalIdeal

def nf_ideal_intersection_rng(cs,fn,fd,aa,ad,ac,bb,bd,bc,seed):
    R=PolynomialRing(QQ,'x');K=NumberField(R([QQ(c)*fn/fd for c in cs]),'a');n=K.degree()
    def ideal(flat,d,count):return K.ideal([K(R([QQ(c)/d for c in flat[j*n:(j+1)*n]]))for j in range(count)])
    I=ideal(aa,ad,ac);J=ideal(bb,bd,bc)
    K.pari_nf();I.pari_hnf();J.pari_hnf()
    pari.setrand(seed)
    result=I.intersection(J);next_value=pari.random(ZZ(2)**128)
    return json.dumps([_lattice(result),str(next_value)],separators=(',',':'))

def _field(cs,name='a'):
    R=PolynomialRing(QQ,'x');return NumberField(R(cs),name)
def _lattice(I):
    if I.is_zero():return []
    rows=[el.list()for el in I.basis()];d=lcm([c.denominator()for row in rows for c in row])
    H=matrix(ZZ,[[ZZ(c*d)for c in row]for row in rows]).hermite_form()
    return [[str(QQ(c)/d)for c in row]for row in H.rows()]
def nf_ideal_inputs(op,cs,aa,bb,d):
    K=_field(cs);R=PolynomialRing(QQ,'x');a=K(R(aa))/d;b=K(R(bb))/d;I=K.ideal(a);J=K.ideal(b)
    try:
        if op==0:L=I.intersection(J)
        elif op==1:L=I.intersection(QQ(3)/d)
        elif op==2:L=I.intersection(b)
        elif op==3:L=I.intersection([b,QQ(3)/d])
        elif op==4:L=I.intersection([])
        elif op==5:L=I.intersection([0])
        elif op==6:L=K.fractional_ideal(0)
        elif op==7:L=NumberFieldFractionalIdeal(K,[K(0)])
        elif op==8:L=K.ideal([a,b])
        elif op==9:L=K.ideal(I)
        elif op==10:L=K.fractional_ideal(I)
        elif op==11:L=K.fractional_ideal([a,b])
        elif op==12:L=NumberFieldFractionalIdeal(K,[])
        elif op==13:L=K.ideal(NumberFieldIdeal(K,[K(3)]))
        elif op==14:L=I.intersection(NumberFieldIdeal(K,[K(3)]))
        elif op==15:return json.dumps([str(K.ideal(0)),None],separators=(',',':'))
        elif op==16:return json.dumps([str(NumberFieldIdeal(K,[K(3)])),None],separators=(',',':'))
        elif op==17:L=I+J
        elif op==18:L=K.ideal(I+J)
        else:L=I.intersection(I+J)
        r=[_lattice(L),L is I,isinstance(L,NumberFieldFractionalIdeal)]
        if op==17:r.append([[str(c)for c in v.list()]for v in L.gens()])
        error=None
    except Exception as e:r=None;error=type(e).__name__+': '+str(e)
    return json.dumps([r,error],separators=(',',':'))
def nf_cross_ideal(cs,mode,op):
    K=_field([2,0,1]);L=_field(cs,'b');v=[L(0),L(2),L.gen(),L.ideal(0),L.ideal(2),L.ideal(L.gen())][int(mode)]
    try:
        r=[str(c)for c in K(v).list()]if op==0 else _lattice(K.ideal(v)if op==1 else K.fractional_ideal(v));error=None
    except Exception as e:r=None;error=type(e).__name__+': '+str(e)
    return json.dumps([r,error],separators=(',',':'))

def nf_ideal_factory(cs,op,kind,layout,n,d):
    K=_field(cs)
    try:
        value=int(n)if kind==0 else ZZ(n)if kind==1 else QQ(n)/d if kind==2 else K(QQ(n)/d)if kind==3 else float(n)/float(d)if d else float('inf')if n>0 else -float('inf')if n<0 else float('nan')
        args=[value]if layout==0 else [[value]]if layout==1 else [value,value]if layout==2 else []if layout==3 else [[]]
        I=K.ideal(*args)if op==0 else K.fractional_ideal(*args)if op==1 else NumberFieldIdeal(K,args)if op==2 else NumberFieldFractionalIdeal(K,args)
        r=[_lattice(I),isinstance(I,NumberFieldFractionalIdeal),[[str(c)for c in v.list()]for v in I.gens()]];error=None
    except Exception as e:r=None;error=type(e).__name__+': '+str(e)
    return json.dumps([r,error],separators=(',',':'))
