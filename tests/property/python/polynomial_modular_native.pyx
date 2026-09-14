# distutils: language = c++
from cysignals.signals cimport sig_on,sig_off
from sage.rings.integer cimport Integer
from sage.libs.ntl.ntl_ZZ_pEX cimport ntl_ZZ_pEX
from sage.libs.ntl.types cimport ZZ_pEX_c, ZZ_pEX_Modulus_c, ZZ_c
from sage.libs.ntl.convert cimport mpz_to_ZZ
cdef extern from "NTL/ZZ_pEX.h":
    void build "NTL::build"(ZZ_pEX_Modulus_c&,const ZZ_pEX_c&) except +
    void PowerMod "NTL::PowerMod"(ZZ_pEX_c&,const ZZ_pEX_c&,const ZZ_c&,const ZZ_pEX_Modulus_c&) except +
    void PowerXMod "NTL::PowerXMod"(ZZ_pEX_c&,const ZZ_c&,const ZZ_pEX_Modulus_c&) except +

def native_extension_modular_power(ntl_ZZ_pEX a,exponent,ntl_ZZ_pEX m,int method):
    cdef ntl_ZZ_pEX out=a._new()
    cdef ZZ_c E
    cdef ZZ_pEX_Modulus_c M
    cdef Integer e=Integer(exponent)
    a.c.restore_c()
    mpz_to_ZZ(&E,e.value)
    sig_on()
    try:
        build(M,m.x)
        if method==0:PowerMod(out.x,a.x,E,M)
        else:PowerXMod(out.x,E,M)
    finally:sig_off()
    return out

from sage.libs.ntl.ntl_GF2X cimport ntl_GF2X
from sage.libs.ntl.types cimport GF2X_c,GF2XModulus_c
cdef extern from "NTL/GF2X.h":
    void binary_build "NTL::build"(GF2XModulus_c&,const GF2X_c&) except +
    void binary_rem "NTL::rem"(GF2X_c&,const GF2X_c&,const GF2X_c&) except +
    void binary_PowerMod "NTL::PowerMod"(GF2X_c&,const GF2X_c&,const ZZ_c&,const GF2XModulus_c&) except +
    void binary_PowerXMod "NTL::PowerXMod"(GF2X_c&,const ZZ_c&,const GF2XModulus_c&) except +
    long binary_IsX "NTL::IsX"(const GF2X_c&)

def native_binary_modular_power(ntl_GF2X a,exponent,ntl_GF2X m):
    cdef ntl_GF2X out=ntl_GF2X.__new__(ntl_GF2X)
    cdef GF2X_c reduced
    cdef GF2XModulus_c M
    cdef ZZ_c E
    cdef Integer e=Integer(exponent)
    mpz_to_ZZ(&E,e.value)
    sig_on()
    try:
        binary_rem(reduced,a.x,m.x)
        binary_build(M,m.x)
        if binary_IsX(reduced):binary_PowerXMod(out.x,E,M)
        else:binary_PowerMod(out.x,reduced,E,M)
    finally:sig_off()
    return out
