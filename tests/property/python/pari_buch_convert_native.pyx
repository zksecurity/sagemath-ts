"""Original constructors, allocation conversions and scalar real operations."""
from cysignals.signals cimport sig_on, sig_off
from cypari2.types cimport GEN, pari_sp
from cypari2.paridecl cimport (itor, rtor, real_0, real_1, real_0_bit, negr, absr,
    shiftr, truncr, trunc2nr, mulur, cgetr, cgeti, evalexpo, evalsigne, lg, signe, expo)
cdef extern from "pari/pari.h":
    pari_sp get_avma()
    void set_avma(pari_sp)

def native_convert(long op,long p,long e,m,long s,long q,n):
    cdef GEN x,z,a
    cdef pari_sp saved
    cdef long i,words,sign,ex
    cdef unsigned long word
    sig_on();saved=get_avma()
    try:
        x=cgetr(p//64+2);x[1]=evalsigne(s)|evalexpo(e)
        for i in range(2,p//64+2):
            word=(m>>(p-64*(i-1)))&((1<<64)-1);x[i]=<long>word
        if op==0:
            mag=abs(n);words=(mag.bit_length()+63)//64+2;a=cgeti(words)
            a[1]=evalsigne(1 if n>0 else -1 if n<0 else 0)|words
            for i in range(2,words):
                word=(mag>>(64*(i-2)))&((1<<64)-1);a[i]=<long>word
            z=itor(a,(q+63)//64+2)
        elif op==1:z=real_0((q+63)//64+2)
        elif op==2:z=real_1((q+63)//64+2)
        elif op==3:z=real_0_bit(e)
        elif op==4:z=rtor(x,(q+63)//64+2)
        elif op==5:z=negr(x)
        elif op==6:z=absr(x)
        elif op==7:z=shiftr(x,n)
        elif op in [8,12]:
            z=truncr(x) if op==8 else trunc2nr(x,n);words=z[1]&((1<<56)-1);value=0
            for i in range(words-1,1,-1):
                word=<unsigned long>z[i];value=(value<<64)+word
            return ['integer',str(-value if signe(z)<0 else value)]
        elif op==9:return ['integer',str(signe(x))]
        elif op==10:return ['integer',str(expo(x))]
        else:z=mulur(n,x)
        words=lg(z);sign=signe(z);ex=expo(z);mantissa=0
        if sign:
            for i in range(2,words):
                word=<unsigned long>z[i];mantissa=(mantissa<<64)+word
        return ['real',sign,str(ex),str(mantissa),(words-2)*64]
    finally:set_avma(saved);sig_off()
