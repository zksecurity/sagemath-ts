# distutils: language = c++
from cysignals.signals cimport sig_on,sig_off
from sage.libs.flint.types cimport fmpz_poly_t,fmpz_t,fmpz
from sage.libs.flint.fmpz_poly cimport fmpz_poly_init,fmpz_poly_clear,fmpz_poly_set_str,_fmpz_poly_evaluate_fmpz
from sage.libs.flint.fmpz cimport fmpz_init,fmpz_clear,fmpz_set_str,fmpz_set_ui,fmpz_get_str
from sage.libs.flint.flint cimport flint_free

cdef hexadecimal(fmpz_t value):
    cdef char * raw=fmpz_get_str(NULL,16,value)
    try:return (<bytes>raw).decode()
    finally:flint_free(raw)

# BUNDLED_EVALUATION_ROUTINES

def native_evaluation(int kind,a,den,num,xden):
    cdef fmpz_poly_t A
    cdef fmpz_t D,XN,XD,RN,RD
    fmpz_poly_init(A)
    fmpz_init(D);fmpz_init(XN);fmpz_init(XD);fmpz_init(RN);fmpz_init(RD)
    try:
        raw=(str(len(a))+'  '+' '.join(map(str,a))).encode()
        if fmpz_poly_set_str(A,raw):raise ValueError('invalid polynomial fixture')
        fmpz_set_str(D,str(den).encode(),10);fmpz_set_str(XN,str(num).encode(),10);fmpz_set_str(XD,str(xden).encode(),10)
        sig_on()
        try:
            if kind==0:
                _fmpz_poly_evaluate_fmpz(RN,A[0].coeffs,A[0].length,XN);fmpz_set_ui(RD,1)
            elif kind==1:audit_fmpz_poly_evaluate_fmpq(RN,RD,A[0].coeffs,A[0].length,XN,XD)
            elif kind==2:audit_fmpq_poly_evaluate_fmpq(RN,RD,A[0].coeffs,D,A[0].length,XN,XD)
            else:audit_fmpq_poly_evaluate_fmpz(RN,RD,A[0].coeffs,D,A[0].length,XN)
        finally:sig_off()
        return [hexadecimal(RN),hexadecimal(RD)]
    finally:
        fmpz_poly_clear(A)
        fmpz_clear(D);fmpz_clear(XN);fmpz_clear(XD);fmpz_clear(RN);fmpz_clear(RD)

from sage.libs.ntl.ntl_ZZ_pEX cimport ntl_ZZ_pEX
from sage.libs.ntl.ntl_ZZ_pE cimport ntl_ZZ_pE
from sage.libs.ntl.types cimport ZZ_pEX_c,ZZ_pE_c
cdef extern from "NTL/ZZ_pEX.h":
    void ntl_eval "NTL::eval"(ZZ_pE_c&,const ZZ_pEX_c&,const ZZ_pE_c&) except +
def native_extension_evaluation(ntl_ZZ_pEX a,ntl_ZZ_pE x):
    cdef ntl_ZZ_pE out=x._new()
    a.c.restore_c()
    sig_on()
    try:ntl_eval(out.x,a.x,x.x)
    finally:sig_off()
    return out

from sage.libs.flint.types cimport fmpq_poly_t,nmod_poly_t
from sage.libs.flint.fmpz_poly cimport fmpz_poly_compose,fmpz_poly_taylor_shift,fmpz_poly_length,fmpz_poly_get_coeff_fmpz
from sage.libs.flint.fmpq_poly cimport fmpq_poly_init,fmpq_poly_clear,fmpq_poly_set_str,fmpq_poly_compose,fmpq_poly_get_numerator,fmpq_poly_get_denominator
from sage.libs.flint.nmod_poly cimport nmod_poly_init,nmod_poly_clear,nmod_poly_set_coeff_ui,nmod_poly_compose,nmod_poly_length,nmod_poly_get_coeff_ui

def native_composition(int kind,a,da,b,db,unsigned long p):
    cdef fmpz_poly_t A,B,R
    cdef fmpq_poly_t QA,QB,QR
    cdef nmod_poly_t MA,MB,MR
    cdef fmpz_t D,C
    fmpz_poly_init(A);fmpz_poly_init(B);fmpz_poly_init(R);fmpz_init(D);fmpz_init(C)
    try:
        if kind==2:
            nmod_poly_init(MA,p);nmod_poly_init(MB,p);nmod_poly_init(MR,p)
            try:
                for i,c in enumerate(a):nmod_poly_set_coeff_ui(MA,i,int(c%p))
                for i,c in enumerate(b):nmod_poly_set_coeff_ui(MB,i,int(c%p))
                sig_on()
                try:nmod_poly_compose(MR,MA,MB)
                finally:sig_off()
                return [[format(nmod_poly_get_coeff_ui(MR,i),'x') for i in range(nmod_poly_length(MR))],'1']
            finally:nmod_poly_clear(MA);nmod_poly_clear(MB);nmod_poly_clear(MR)
        elif kind==1:
            fmpq_poly_init(QA);fmpq_poly_init(QB);fmpq_poly_init(QR)
            try:
                raw=(str(len(a))+'  '+' '.join(str(c)+'/'+str(da) for c in a)).strip().encode()
                if fmpq_poly_set_str(QA,raw):raise ValueError('invalid rational polynomial fixture')
                raw=(str(len(b))+'  '+' '.join(str(c)+'/'+str(db) for c in b)).strip().encode()
                if fmpq_poly_set_str(QB,raw):raise ValueError('invalid rational polynomial fixture')
                sig_on()
                try:fmpq_poly_compose(QR,QA,QB)
                finally:sig_off()
                fmpq_poly_get_numerator(R,QR);fmpq_poly_get_denominator(D,QR)
            finally:fmpq_poly_clear(QA);fmpq_poly_clear(QB);fmpq_poly_clear(QR)
        else:
            raw=(str(len(a))+'  '+' '.join(map(str,a))).encode()
            if fmpz_poly_set_str(A,raw):raise ValueError('invalid integer polynomial fixture')
            if kind==3:fmpz_set_str(C,str(b[0]).encode(),10)
            else:
                raw=(str(len(b))+'  '+' '.join(map(str,b))).encode()
                if fmpz_poly_set_str(B,raw):raise ValueError('invalid integer polynomial fixture')
            sig_on()
            try:
                if kind==3:fmpz_poly_taylor_shift(R,A,C)
                else:fmpz_poly_compose(R,A,B)
            finally:sig_off()
            fmpz_set_ui(D,1)
        values=[]
        for i in range(fmpz_poly_length(R)):
            fmpz_poly_get_coeff_fmpz(C,R,i);values.append(hexadecimal(C))
        return [values,hexadecimal(D)]
    finally:fmpz_poly_clear(A);fmpz_poly_clear(B);fmpz_poly_clear(R);fmpz_clear(D);fmpz_clear(C)

def native_integer_double(value):
    cdef fmpz_t n
    fmpz_init(n)
    try:
        fmpz_set_str(n,str(value).encode(),10)
        return auditfmpz_get_d(n)
    finally:fmpz_clear(n)
