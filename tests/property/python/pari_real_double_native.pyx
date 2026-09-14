"""PARI binary64 conversions, including native overflow errors, under sig_on."""
from cysignals.signals cimport sig_on, sig_off
from cypari2.types cimport GEN, pari_sp
from cypari2.paridecl cimport dbltor, rtodbl, rtor, cgetr, evalexpo, evalsigne, lg, signe, expo
cdef extern from "pari/pari.h":
    pari_sp get_avma()
    void set_avma(pari_sp)
cdef extern from *:
    "typedef union {double value; unsigned long bits;} audit_pari_double;"
    ctypedef union audit_pari_double:
        double value
        unsigned long bits

def native_double(long op,long p,long e,m,long sign):
    cdef GEN x,z
    cdef pari_sp saved
    cdef audit_pari_double bits
    cdef long i,words,s,ex
    cdef unsigned long word
    sig_on();saved=get_avma()
    try:
        if op in [0,2]:
            bits.bits=m;z=dbltor(bits.value)
            if op==2 and p:z=rtor(z,(p+63)//64+2)
            words=lg(z);s=signe(z);ex=expo(z);mantissa=0
            if s:
                for i in range(2,words):
                    word=<unsigned long>z[i];mantissa=(mantissa<<64)+word
            return [s,str(ex),str(mantissa),(words-2)*64]
        words=p//64+2;x=cgetr(words);x[1]=evalsigne(sign)|evalexpo(e)
        for i in range(2,words):
            word=(m>>(p-64*(i-1)))&((1<<64)-1)
            x[i]=<long>word
        bits.value=rtodbl(x)
        return format(bits.bits,'016x')
    finally:set_avma(saved);sig_off()
