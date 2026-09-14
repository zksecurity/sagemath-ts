"""Native scalar actions and generic group calls over number fields."""
import json
from sage.all import *

def _field(d):
    if d: return QuadraticField(d, 'a')
    x = polygen(QQ)
    return NumberField(x**3-2, 'a')

def nf_scalar(d, c, n, den, kind):
    K = _field(d)
    a = (K.gen()+c)/3
    scalar = QQ(n)/den if kind == 4 else ZZ(n)
    return json.dumps([str(x) for x in (a*scalar).list()], separators=(',', ':'))

def nf_group(g, d, fn, mode, n, variant):
    K = _field(d)
    op = '+' if mode == 0 else '*'
    a = K.gen() if variant == 0 else K(-1) if variant == 1 else K.zero() if op == '+' else K.one()
    target = n*a if op == '+' else a**n
    if fn == 0: r = [str(x) for x in g.multiple(a,n,operation=op).list()]
    elif fn == 1: r = str(g.bsgs(a,target,(0,64),operation=op))
    elif fn == 2: r = bool(g.has_order(a,n,operation=op))
    elif fn == 3: r = str(g.order_from_multiple(a,n,operation=op))
    elif fn in (4,5): r = str(g.discrete_log(target,a,ord=4,operation=op))
    else: r = str(g.order_from_bounds(a,(1,16),operation=op))
    return json.dumps(r, separators=(',', ':'))
