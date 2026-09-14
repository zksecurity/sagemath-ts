"""Call the original PARI comparison on explicitly allocated real records."""
from cysignals.signals cimport sig_on, sig_off
from cypari2.types cimport GEN, pari_sp
from cypari2.paridecl cimport cgetr, evalexpo, evalsigne, cmprr
cdef extern from "pari/pari.h":
    pari_sp get_avma()
    void set_avma(pari_sp)

cdef GEN record(long p,long e,m,long s):
    cdef GEN x=cgetr(p//64+2)
    cdef long i
    cdef unsigned long word
    x[1]=evalsigne(s)|evalexpo(e)
    for i in range(2,p//64+2):
        word=(m>>(p-64*(i-1)))&((1<<64)-1)
        x[i]=<long>word
    return x

def native_compare(long p,long e,m,long s,long q,long f,n,long t):
    cdef pari_sp saved
    sig_on();saved=get_avma()
    try:return cmprr(record(p,e,m,s),record(q,f,n,t))
    finally:set_avma(saved);sig_off()
