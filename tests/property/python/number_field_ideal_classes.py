"""Method availability on direct base ideals versus factory fractional ideals."""
import json
from sage.all import *
from sage.rings.number_field.number_field_ideal import NumberFieldIdeal
from number_field_ideal_inputs import _field, _lattice


def nf_ideal_class_method(cs, factory, n, d, op):
    K=_field(cs);q=QQ(n)/d
    I=K.ideal(q) if factory else NumberFieldIdeal(K,[K(q)])
    if op==0:r=~I
    elif op==1:r=I.numerator()
    elif op==2:r=I.denominator()
    elif op==3:r=I.divides(K.ideal(6))
    elif op==4:r=I.is_coprime(K.ideal(3))
    elif op==5:r=I**int(-1)
    elif op==6:r=I**ZZ(-1)
    elif op==7:r=I/K.ideal(2)
    elif op==8:r=I/int(2)
    elif op==9:r=I**2
    elif op==10:r=[[str(P.norm()),str(e)]for P,e in I.factor()]
    elif op==11:r=I.ramification_index()
    elif op==12:r=I.residue_class_degree()
    elif op==13:
        F=I.residue_field();r=[str(F.characteristic()),str(F.order()),str(F.degree())]
    elif op==14:r=I.is_maximal()
    elif op==15:r=I.is_integral()
    elif op==16:r=I**0
    else:r=I**1
    if isinstance(r,NumberFieldIdeal):r=[type(r).__name__,_lattice(r)]
    elif not isinstance(r,bool) and isinstance(r,(int,Integer)):r=str(r)
    return json.dumps(r,separators=(',',':'))


def nf_prime_ideal_class(cs,p,op):
    K=_field(cs);out=[]
    for I in K.primes_above(p):
        if op==0:r=type(I).__name__
        elif op==1:r=[type(~I).__name__,str((~I).norm())]
        elif op==2:r=[str(I.ramification_index()),str(I.residue_class_degree())]
        elif op==3:r=I.divides(p)
        elif op==4:r=I.is_coprime(p)
        else:r=[_lattice(I.numerator()),_lattice(I.denominator())]
        out.append([_lattice(I),r])
    # Prime enumeration order is audited separately; compare each canonical ideal.
    return json.dumps(sorted(out,key=lambda row:json.dumps(row[0],separators=(',',':'))),separators=(',',':'))
