"""Sage ideal valuation: receiver meaning, fractional inputs and coercion precedence."""
import json
from sage.all import *
from sage.rings.number_field.number_field_ideal import NumberFieldIdeal
from number_field_ideal_inputs import _field,_lattice

def nf_ideal_valuation(cs,p,flat,d,count,base):
    K=_field(cs);n=K.degree();gens=[K([QQ(c)/d for c in flat[i*n:(i+1)*n]])for i in range(count)]
    I=NumberFieldIdeal(K,gens) if base else K.ideal(gens)
    out=[[_lattice(P),str(I.valuation(P))]for P in K.primes_above(p)]
    return json.dumps(sorted(out,key=lambda row:json.dumps(row[0],separators=(',',':'))),separators=(',',':'))

def nf_ideal_valuation_input(cs,receiver,kind,n,d):
    K=_field(cs);L=_field(cs,'b');q=QQ(n)/d
    I=[K.ideal(0),K.ideal(1),K.ideal(6),K.ideal(K.gen()+1),NumberFieldIdeal(K,[K(0)]),NumberFieldIdeal(K,[K(6)])][receiver]
    values=[int(n),ZZ(n),q,float(q),K(q),K.gen()+q,[],[q],[q]*K.degree(),K.ideal(q),K.ideal(0),NumberFieldIdeal(K,[K(q)]),L(q),L.gen()+q,L.ideal(q),L.prime_above(2),'2',[K.gen(),q]]
    try:return str(I.valuation(values[kind]))
    except ValueError as error:
        # Only multi-generator display is normalized: the pre-existing BNF
        # pretty-printing dependency remains open. Preserve error class/message.
        if kind in (8,17) and str(error).startswith('p (= '):
            P=K.ideal(values[kind]);raise ValueError('p (= '+json.dumps(_lattice(P),separators=(',',':'))+')'+str(error).rsplit(')',1)[1]) from None
        raise

def nf_prime_below(cs,p):
    K=_field(cs)
    return json.dumps(sorted([[_lattice(P),str(P.smallest_integer())]for P in K.primes_above(p)],key=lambda row:json.dumps(row[0],separators=(',',':'))),separators=(',',':'))
