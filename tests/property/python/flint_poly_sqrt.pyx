from cysignals.signals cimport sig_on, sig_off
from sage.libs.flint.types cimport nmod_poly_t
from sage.libs.flint.nmod_poly cimport (
    nmod_poly_init, nmod_poly_clear, nmod_poly_set_coeff_ui,
    nmod_poly_get_coeff_ui, nmod_poly_length, nmod_poly_sqrt,
    nmod_poly_sqrt_series, nmod_poly_invsqrt_series, nmod_poly_inv_series,
    nmod_poly_mulhigh, nmod_poly_mulhigh_classical,
)
def native_poly(int op,a,b,unsigned long p,long n):
    cdef nmod_poly_t A,B,C
    cdef int success=1
    nmod_poly_init(A,p);nmod_poly_init(B,p);nmod_poly_init(C,p)
    try:
        for i,c in enumerate(a):nmod_poly_set_coeff_ui(A,i,c%p)
        for i,c in enumerate(b):nmod_poly_set_coeff_ui(C,i,c%p)
        sig_on()
        try:
            if op==0:success=nmod_poly_sqrt(B,A)
            elif op==1:nmod_poly_sqrt_series(B,A,n)
            elif op==2:nmod_poly_invsqrt_series(B,A,n)
            elif op==3:nmod_poly_inv_series(B,A,n)
            elif op==4:nmod_poly_mulhigh(B,A,C,n)
            elif op==5:nmod_poly_mulhigh_classical(B,A,C,n)
        finally:sig_off()
        if not success:return None
        return [nmod_poly_get_coeff_ui(B,i)for i in range(nmod_poly_length(B))]
    finally:nmod_poly_clear(A);nmod_poly_clear(B);nmod_poly_clear(C)

from sage.libs.flint.types cimport nmod_poly_factor_t
from sage.libs.flint.nmod_poly_factor cimport nmod_poly_factor_init,nmod_poly_factor_clear,nmod_poly_factor_squarefree

def native_squarefree(a,unsigned long p):
    cdef nmod_poly_t A
    cdef nmod_poly_factor_t F
    nmod_poly_init(A,p);nmod_poly_factor_init(F)
    try:
        for i,c in enumerate(a):nmod_poly_set_coeff_ui(A,i,c%p)
        sig_on()
        try:nmod_poly_factor_squarefree(F,A)
        finally:sig_off()
        return [([nmod_poly_get_coeff_ui(F.p+i,j)for j in range(nmod_poly_length(F.p+i))],F.exp[i])for i in range(F.num)]
    finally:nmod_poly_clear(A);nmod_poly_factor_clear(F)


from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint
from sage.rings.polynomial.polynomial_rational_flint cimport Polynomial_rational_flint
from sage.libs.flint.fmpz_poly cimport fmpz_poly_lcm
from sage.libs.flint.fmpq_poly cimport fmpq_poly_lcm

def native_lcm(int op,a,b):
    from sage.all import ZZ,QQ
    cdef Polynomial_integer_dense_flint A=ZZ['x'](a),B=ZZ['x'](b),C=ZZ['x'](0)
    cdef Polynomial_rational_flint D=QQ['x'](a),E=QQ['x'](b),F=QQ['x'](0)
    sig_on()
    try:
        if op==0:fmpz_poly_lcm(C._poly,A._poly,B._poly)
        else:fmpq_poly_lcm(F._poly,D._poly,E._poly)
    finally:sig_off()
    return C.list()if op==0 else (F.numerator().list(),F.denominator())
