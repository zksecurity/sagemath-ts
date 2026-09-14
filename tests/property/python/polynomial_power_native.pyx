from cysignals.signals cimport sig_on,sig_off
from libc.stdio cimport fflush,stdout

def safe_power(f,e):
    sig_on()
    try:
        return f**e
    finally:
        sig_off()

from sage.rings.integer cimport Integer
from sage.libs.flint.types cimport fmpz_poly_t, fmpz_t, fmpq_poly_t
from sage.libs.flint.flint cimport flint_free
from sage.libs.flint.fmpz_poly cimport (
    fmpz_poly_init, fmpz_poly_clear, fmpz_poly_set_str, fmpz_poly_get_str, fmpz_poly_length, fmpz_poly_get_coeff_fmpz,
    fmpz_poly_pow, fmpz_poly_pow_binomial, fmpz_poly_pow_multinomial, fmpz_poly_pow_binexp,
)
from sage.libs.flint.fmpz cimport fmpz_init, fmpz_clear, fmpz_get_mpz, fmpz_get_str
from sage.libs.flint.fmpq_poly cimport (
    fmpq_poly_init,fmpq_poly_clear,fmpq_poly_set_str,fmpq_poly_pow,
    fmpq_poly_get_numerator,fmpq_poly_get_denominator,
)

cdef poly_strings(fmpz_poly_t poly):
    # Native hexadecimal serialization avoids quadratic decimal conversion.
    cdef fmpz_t c
    cdef char * raw
    fmpz_init(c)
    try:
        result=[]
        for i in range(fmpz_poly_length(poly)):
            fmpz_poly_get_coeff_fmpz(c,poly,i)
            raw=fmpz_get_str(NULL,16,c)
            try:result.append((<bytes>raw).decode())
            finally:flint_free(raw)
        return result
    finally:fmpz_clear(c)

def native_integer_power(a,unsigned long e,int method):
    cdef fmpz_poly_t A,R
    fmpz_poly_init(A);fmpz_poly_init(R)
    try:
        raw=(str(len(a))+'  '+' '.join(map(str,a))).encode()
        if fmpz_poly_set_str(A,raw):raise ValueError('invalid polynomial fixture')
        sig_on()
        try:
            if method in (0,4):fmpz_poly_pow(R,A,e)
            elif method==1:fmpz_poly_pow_binomial(R,A,e)
            elif method==2:fmpz_poly_pow_multinomial(R,A,e)
            else:fmpz_poly_pow_binexp(R,A,e)
        finally:sig_off()
        return poly_strings(R)
    finally:fmpz_poly_clear(A);fmpz_poly_clear(R)

def native_rational_power(a,den,unsigned long e):
    cdef fmpq_poly_t A,R
    cdef fmpz_poly_t N
    cdef fmpz_t D
    cdef Integer d=Integer.__new__(Integer)
    fmpq_poly_init(A);fmpq_poly_init(R);fmpz_poly_init(N);fmpz_init(D)
    try:
        raw=(str(len(a))+'  '+' '.join(str(c)+'/'+str(den) for c in a)).strip().encode()
        if fmpq_poly_set_str(A,raw):raise ValueError('invalid polynomial fixture')
        sig_on()
        try:fmpq_poly_pow(R,A,e)
        finally:sig_off()
        fmpq_poly_get_numerator(N,R);fmpq_poly_get_denominator(D,R)
        fmpz_get_mpz(d.value,D)
        return [poly_strings(N),str(d)]
    finally:fmpq_poly_clear(A);fmpq_poly_clear(R);fmpz_poly_clear(N);fmpz_clear(D)

from sage.libs.flint.fmpz_poly cimport fmpz_poly_mullow,fmpz_poly_pow_trunc,fmpz_poly_inv_series
from sage.libs.flint.fmpq_poly cimport fmpq_poly_mullow,fmpq_poly_inv_series
from sage.libs.flint.types cimport nmod_poly_t
from sage.libs.flint.nmod_poly cimport nmod_poly_init,nmod_poly_clear,nmod_poly_set_coeff_ui,nmod_poly_get_coeff_ui,nmod_poly_length,nmod_poly_mullow,nmod_poly_pow_trunc

def native_polynomial_series(int kind,unsigned long p,a,da,b,db,unsigned long e,long n):
    cdef fmpz_poly_t A,B,C,N
    cdef fmpq_poly_t QA,QB,QC
    cdef nmod_poly_t MA,MB,MC
    cdef fmpz_t D
    cdef Integer denominator=Integer.__new__(Integer)
    if kind in (2,4):
        nmod_poly_init(MA,p);nmod_poly_init(MB,p);nmod_poly_init(MC,p)
        try:
            for i,c in enumerate(a):nmod_poly_set_coeff_ui(MA,i,int(c)%p)
            for i,c in enumerate(b):nmod_poly_set_coeff_ui(MB,i,int(c)%p)
            sig_on()
            try:
                if kind==2:nmod_poly_mullow(MC,MA,MB,n)
                else:nmod_poly_pow_trunc(MC,MA,e,n)
            finally:sig_off()
            return [[format(nmod_poly_get_coeff_ui(MC,i),'x') for i in range(nmod_poly_length(MC))],'1']
        finally:nmod_poly_clear(MA);nmod_poly_clear(MB);nmod_poly_clear(MC)
    if kind in (0,3,5):
        fmpz_poly_init(A);fmpz_poly_init(B);fmpz_poly_init(C)
        try:
            raw=(str(len(a))+'  '+' '.join(map(str,a))).strip().encode()
            if fmpz_poly_set_str(A,raw):raise ValueError('invalid polynomial fixture')
            raw=(str(len(b))+'  '+' '.join(map(str,b))).strip().encode()
            if fmpz_poly_set_str(B,raw):raise ValueError('invalid polynomial fixture')
            sig_on()
            try:
                if kind==0:fmpz_poly_mullow(C,A,B,n)
                elif kind==3:fmpz_poly_pow_trunc(C,A,e,n)
                else:fmpz_poly_inv_series(C,A,n)
            finally:sig_off()
            return [poly_strings(C),'1']
        finally:fmpz_poly_clear(A);fmpz_poly_clear(B);fmpz_poly_clear(C)
    fmpq_poly_init(QA);fmpq_poly_init(QB);fmpq_poly_init(QC);fmpz_poly_init(N);fmpz_init(D)
    try:
        raw=(str(len(a))+'  '+' '.join(str(c)+'/'+str(da) for c in a)).strip().encode()
        if fmpq_poly_set_str(QA,raw):raise ValueError('invalid polynomial fixture')
        raw=(str(len(b))+'  '+' '.join(str(c)+'/'+str(db) for c in b)).strip().encode()
        if fmpq_poly_set_str(QB,raw):raise ValueError('invalid polynomial fixture')
        sig_on()
        try:
            if kind==1:fmpq_poly_mullow(QC,QA,QB,n)
            else:fmpq_poly_inv_series(QC,QA,n)
        finally:sig_off()
        fmpq_poly_get_numerator(N,QC);fmpq_poly_get_denominator(D,QC)
        fmpz_get_mpz(denominator.value,D)
        return [poly_strings(N),str(denominator)]
    finally:fmpq_poly_clear(QA);fmpq_poly_clear(QB);fmpq_poly_clear(QC);fmpz_poly_clear(N);fmpz_clear(D)

from sage.libs.flint.fmpz_poly cimport fmpz_poly_mul
from sage.libs.flint.fmpq_poly cimport fmpq_poly_mul

def native_full_product(a,da,b,db,int rational,int same):
    cdef fmpz_poly_t A,B,R
    cdef fmpq_poly_t QA,QB,QR
    cdef fmpz_t D
    cdef Integer d=Integer.__new__(Integer)
    fmpz_poly_init(A);fmpz_poly_init(B);fmpz_poly_init(R)
    fmpq_poly_init(QA);fmpq_poly_init(QB);fmpq_poly_init(QR);fmpz_init(D)
    try:
        if rational:
            raw=(str(len(a))+'  '+' '.join(str(c)+'/'+str(da) for c in a)).strip().encode()
            if fmpq_poly_set_str(QA,raw):raise ValueError('invalid polynomial fixture')
            raw=(str(len(b))+'  '+' '.join(str(c)+'/'+str(db) for c in b)).strip().encode()
            if fmpq_poly_set_str(QB,raw):raise ValueError('invalid polynomial fixture')
            sig_on()
            try:
                if same:fmpq_poly_mul(QR,QA,QA)
                else:fmpq_poly_mul(QR,QA,QB)
            finally:sig_off()
            fmpq_poly_get_numerator(R,QR);fmpq_poly_get_denominator(D,QR)
            fmpz_get_mpz(d.value,D)
            return [poly_strings(R),str(d)]
        raw=(str(len(a))+'  '+' '.join(map(str,a))).encode()
        if fmpz_poly_set_str(A,raw):raise ValueError('invalid polynomial fixture')
        raw=(str(len(b))+'  '+' '.join(map(str,b))).encode()
        if fmpz_poly_set_str(B,raw):raise ValueError('invalid polynomial fixture')
        sig_on()
        try:
            if same:fmpz_poly_mul(R,A,A)
            else:fmpz_poly_mul(R,A,B)
        finally:sig_off()
        return [poly_strings(R),'1']
    finally:
        fmpz_poly_clear(A);fmpz_poly_clear(B);fmpz_poly_clear(R)
        fmpq_poly_clear(QA);fmpq_poly_clear(QB);fmpq_poly_clear(QR);fmpz_clear(D)

from sage.libs.gsl.math cimport gsl_pow_int
from sage.libs.gsl.log cimport gsl_sf_log_e
from sage.libs.gsl.exp cimport gsl_sf_exp_e
from sage.libs.gsl.types cimport gsl_sf_result

def native_gsl_power_value(double x,int n,int method):
    cdef gsl_sf_result r
    if method==0:return gsl_pow_int(x,n)
    if method==1:gsl_sf_log_e(x,&r)
    else:gsl_sf_exp_e(x,&r)
    return r.val

from sage.libs.flint.nmod_poly cimport (
    nmod_poly_inv_series,nmod_poly_reverse,nmod_poly_powmod_ui_binexp,
    nmod_poly_powmod_fmpz_binexp_preinv,nmod_poly_powmod_x_fmpz_preinv,
)
from sage.libs.flint.fmpz cimport fmpz_set_mpz

def native_modular_power(unsigned long p,a,m,exponent,int method,long precision):
    cdef nmod_poly_t A,M,I,R
    cdef fmpz_t E
    cdef Integer e=Integer(exponent)
    nmod_poly_init(A,p);nmod_poly_init(M,p);nmod_poly_init(I,p);nmod_poly_init(R,p);fmpz_init(E)
    try:
        for i,c in enumerate(a):nmod_poly_set_coeff_ui(A,i,int(c)%p)
        for i,c in enumerate(m):nmod_poly_set_coeff_ui(M,i,int(c)%p)
        fmpz_set_mpz(E,e.value)
        sig_on()
        try:
            if method==3:nmod_poly_inv_series(R,A,precision)
            elif method==0:nmod_poly_powmod_ui_binexp(R,A,int(exponent),M)
            else:
                nmod_poly_reverse(I,M,nmod_poly_length(M))
                nmod_poly_inv_series(I,I,nmod_poly_length(M))
                if method==1:nmod_poly_powmod_fmpz_binexp_preinv(R,A,E,M,I)
                else:nmod_poly_powmod_x_fmpz_preinv(R,E,M,I)
        finally:sig_off()
        return [format(nmod_poly_get_coeff_ui(R,i),'x') for i in range(nmod_poly_length(R))]
    finally:fflush(stdout);nmod_poly_clear(A);nmod_poly_clear(M);nmod_poly_clear(I);nmod_poly_clear(R);fmpz_clear(E)
