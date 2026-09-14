"""Legacy integer factor helpers versus bundled PARI, including random state."""
import json
from math import gcd
from functools import reduce
from pari_qx_factor import pari_qx_factor


def nf_legacy_polynomial_factor(op, coefficients):
    f = list(map(int, coefficients))
    while f and not f[-1]:
        f.pop()
    if int(op) == 1:
        value, state = json.loads(pari_qx_factor(14, 2, 0, 0, 0, 1, f, [], []))
        result = value == '1'
    else:
        # The convenience API returns primitive factors without content or
        # multiplicities, and an empty list for zero/constant polynomials.
        content = reduce(gcd, f, 0)
        if f and f[-1] < 0:
            content = -content
        if content:
            f = [c // content for c in f]
        value, state = json.loads(pari_qx_factor(5, 2, 0, 0, 0, 1, f, [], []))
        result = value[0] if len(f) > 1 else []
    return json.dumps([result, state], separators=(',', ':'))
