"""Bundled indexing branches via native absolute/quadratic field classes."""
from sage.all import *

def nf_index(coeffs, denominator, value, n, den, kind):
    x=polygen(QQ)
    K=NumberField(sum(QQ(c)/denominator*x**i for i,c in enumerate(coeffs)), 'a')
    a=K.zero() if value==0 else K(QQ(2)/3) if value==1 else K.gen()+QQ(2)/3
    # Integral host-number indices adapt to Python ints, as for other array APIs.
    i=int(n) if kind==1 or kind==0 and den==1 else ZZ(n) if kind==2 else float(n)/float(den) if den else float('nan') if n==0 else float('inf') if n>0 else -float('inf')
    return str(a[i])
