"""Native logarithms, logarithm constants and exponentials on complete real records."""
from cysignals.signals cimport sig_on, sig_off
from cypari2.types cimport GEN, pari_sp
from cypari2.paridecl cimport cgetr, evalexpo, evalsigne, lg, signe, expo, logr_abs, mpexp, mplog2
cdef extern from "pari/pari.h":
    pari_sp get_avma()
    void set_avma(pari_sp)
    long EXPNEWTON_LIMIT
    long LOGAGM_LIMIT
    long INVNEWTON_LIMIT

def native_logexp(long op,long p,long e,m,long s,long q):
    cdef GEN x,z
    cdef pari_sp saved
    cdef long i,words,sign,ex,old_exp,old_log,old_inv
    cdef unsigned long word
    global EXPNEWTON_LIMIT,LOGAGM_LIMIT,INVNEWTON_LIMIT
    sig_on();saved=get_avma()
    old_exp=EXPNEWTON_LIMIT;old_log=LOGAGM_LIMIT;old_inv=INVNEWTON_LIMIT
    EXPNEWTON_LIMIT=4224//64+2;LOGAGM_LIMIT=384//64+2;INVNEWTON_LIMIT=4800//64+2
    try:
        x=cgetr(p//64+2);x[1]=evalsigne(s)|evalexpo(e)
        for i in range(2,p//64+2):
            word=(m>>(p-64*(i-1)))&((1<<64)-1);x[i]=<long>word
        if op<2:z=logr_abs(x)
        elif op<4:z=mpexp(x)
        else:z=mplog2((q+63)//64+2)
        words=lg(z);sign=signe(z);ex=expo(z);mantissa=0
        if sign:
            for i in range(2,words):
                word=<unsigned long>z[i];mantissa=(mantissa<<64)+word
        return [sign,str(ex),str(mantissa),(words-2)*64]
    finally:
        EXPNEWTON_LIMIT=old_exp;LOGAGM_LIMIT=old_log;INVNEWTON_LIMIT=old_inv
        set_avma(saved);sig_off()
