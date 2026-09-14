# Installed Sage 10.3 lacks the signal guard present in bundled nmod_poly_linkage.pxi.
# This invokes the original implementation; no polynomial algorithm is reproduced here.
from cysignals.signals cimport sig_on, sig_off

def guarded_gcd(a, b):
    sig_on()
    try:
        return a.gcd(b)
    finally:
        sig_off()

def guarded_xgcd(a, b):
    sig_on()
    try:
        return a.xgcd(b)
    finally:
        sig_off()
