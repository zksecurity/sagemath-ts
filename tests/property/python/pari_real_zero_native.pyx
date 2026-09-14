"""Native allocated-zero operations under cypari2's signal/error handler."""
from cysignals.signals cimport sig_on, sig_off
from cypari2.types cimport GEN, pari_sp
from cypari2.paridecl cimport (itor, rtor, stoi, gen_0,
    real_0, real_0_bit, real_1, setexpo, mulir, divir, invr, divrr,
    mulrr, sqrr, addrr, rdivii, lg, signe, expo)

cdef extern from "pari/pari.h":
    pari_sp get_avma()
    void set_avma(pari_sp)

def native_zero(long op, long p, long q, long e, long n):
    cdef GEN x, z
    cdef pari_sp saved
    cdef long sign, exponent, words, i
    cdef unsigned long word
    sig_on()
    saved=get_avma()
    try:
        # Installed headers use word lengths; bundled headers use bit precision.
        x=itor(gen_0,p//64+2) if p else real_0_bit(e)
        setexpo(x,e)
        if op==0:z=itor(gen_0,q//64+2)
        elif op==1:z=rtor(x,q//64+2)
        elif op==2:z=mulir(stoi(n),x)
        elif op==3:z=divir(stoi(n),x)
        elif op==4:z=invr(x)
        elif op==5:z=divrr(real_1(q//64+2),x)
        elif op==6:z=mulrr(x,real_1(q//64+2))
        elif op==7:z=sqrr(x)
        elif op==8:z=addrr(x,real_1(q//64+2))
        elif op==9:z=rdivii(gen_0,stoi(n),q//64+2)
        else:z=real_0(q//64+2)
        sign=signe(z);exponent=expo(z);words=lg(z)
        mantissa=0
        # Zero payload words are unspecified; allocation and accuracy are meaningful.
        if sign:
            for i in range(2,words):
                word=<unsigned long>z[i]
                mantissa=(mantissa<<64)+word
        return [sign,str(exponent),str(mantissa),(words-2)*64]
    finally:
        set_avma(saved)
        sig_off()
