# Direct original FLINT multiplication/power oracle; no polynomial arithmetic is reproduced.
from sage.libs.flint.types cimport nmod_poly_t, nmod_poly_struct
from sage.libs.flint.nmod_poly cimport (
    nmod_poly_init, nmod_poly_clear, nmod_poly_set_coeff_ui,
    nmod_poly_get_coeff_ui, nmod_poly_length, nmod_poly_mul,
    nmod_poly_mul_classical, nmod_poly_mul_KS, nmod_poly_mul_KS2,
    nmod_poly_mul_KS4, nmod_poly_pow, nmod_poly_add, nmod_poly_sub, nmod_poly_neg,
)

def native_product(unsigned long p,a,b,int method,bint same=False):
    cdef nmod_poly_t A,B,R
    cdef nmod_poly_struct * other
    nmod_poly_init(A,p);nmod_poly_init(B,p);nmod_poly_init(R,p)
    try:
        for i,c in enumerate(a):nmod_poly_set_coeff_ui(A,i,c%p)
        for i,c in enumerate(b):nmod_poly_set_coeff_ui(B,i,c%p)
        if same:other=A
        else:other=B
        if method==0:nmod_poly_mul(R,A,other)
        elif method==1:nmod_poly_mul_classical(R,A,other)
        elif method==2:nmod_poly_mul_KS(R,A,other,0)
        elif method==3:nmod_poly_mul_KS2(R,A,other)
        elif method==4:nmod_poly_mul_KS4(R,A,other)
        elif method==5:nmod_poly_add(R,A,other)
        elif method==6:nmod_poly_sub(R,A,other)
        else:nmod_poly_neg(R,A)
        return [str(nmod_poly_get_coeff_ui(R,i)) for i in range(nmod_poly_length(R))]
    finally:
        nmod_poly_clear(A);nmod_poly_clear(B);nmod_poly_clear(R)

def native_power(unsigned long p,a,unsigned long e):
    cdef nmod_poly_t A,R
    nmod_poly_init(A,p);nmod_poly_init(R,p)
    try:
        for i,c in enumerate(a):nmod_poly_set_coeff_ui(A,i,c%p)
        nmod_poly_pow(R,A,e)
        return [str(nmod_poly_get_coeff_ui(R,i)) for i in range(nmod_poly_length(R))]
    finally:
        nmod_poly_clear(A);nmod_poly_clear(R)
