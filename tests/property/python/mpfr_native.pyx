"""Serialize actual MPFR setter state, binary64 rounding and decimal digit output."""
from sage.libs.mpfr cimport mpfr_t,mpfr_init2,mpfr_clear,mpfr_set_str,mpfr_set_d,mpfr_set_z,mpfr_get_d,mpfr_get_str,mpfr_free_str,mpfr_get_z_exp,mpfr_nan_p,mpfr_inf_p,mpfr_zero_p,mpfr_signbit,MPFR_RNDN
from sage.rings.integer cimport Integer

cdef state_frame(mpfr_t x,int precision,int digits,int status):
    import struct,math
    cdef Integer z=Integer(0)
    cdef long exponent=0,decimal_exponent=0
    cdef char* formatted=NULL
    cdef double number
    try:
        kind='nan' if mpfr_nan_p(x) else 'inf' if mpfr_inf_p(x) else 'zero' if mpfr_zero_p(x) else 'finite'
        if kind=='finite':exponent=mpfr_get_z_exp(z.value,x)+precision
        number=mpfr_get_d(x,MPFR_RNDN)
        formatted=mpfr_get_str(NULL,&decimal_exponent,10,digits,x,MPFR_RNDN)
        return [status,kind,-1 if mpfr_signbit(x) else 1,str(abs(z)),int(exponent),precision,'NaN' if math.isnan(number) else struct.pack('>d',number).hex(),formatted.decode('ascii'),int(decimal_exponent)]
    finally:
        if formatted!=NULL:mpfr_free_str(formatted)

def native_conversion(str text,int precision,int digits):
    cdef mpfr_t x
    cdef bytes encoded=text.encode('ascii')
    cdef int status
    mpfr_init2(x,precision)
    try:
        status=mpfr_set_str(x,encoded,10,MPFR_RNDN)
        return state_frame(x,precision,digits,status)
    finally:mpfr_clear(x)

def native_sequence(texts,int precision,int digits):
    cdef mpfr_t x
    cdef bytes encoded
    cdef int status
    mpfr_init2(x,precision)
    result=[]
    try:
        for text in texts:
            encoded=text.encode('ascii');status=mpfr_set_str(x,encoded,10,MPFR_RNDN)
            result.append(state_frame(x,precision,digits,status))
        return result
    finally:mpfr_clear(x)

def native_double(unsigned long long raw,int precision,int digits):
    import struct
    cdef mpfr_t x
    cdef double value=struct.unpack('>d',raw.to_bytes(8,'big'))[0]
    cdef int status
    mpfr_init2(x,precision)
    try:
        status=mpfr_set_d(x,value,MPFR_RNDN)
        return state_frame(x,precision,digits,status)
    finally:mpfr_clear(x)

def native_double_sequence(raws,int precision,int digits):
    import struct
    cdef mpfr_t x
    cdef double value
    cdef int status
    mpfr_init2(x,precision)
    result=[]
    try:
        for raw in raws:
            value=struct.unpack('>d',int(raw).to_bytes(8,'big'))[0]
            status=mpfr_set_d(x,value,MPFR_RNDN)
            result.append(state_frame(x,precision,digits,status))
        return result
    finally:mpfr_clear(x)


def native_integer(value,int precision,int digits):
    cdef mpfr_t x
    cdef Integer integer=Integer(value)
    cdef int status
    mpfr_init2(x,precision)
    try:
        status=mpfr_set_z(x,integer.value,MPFR_RNDN)
        return state_frame(x,precision,digits,status)
    finally:mpfr_clear(x)


from sage.libs.mpfr cimport mpfr_sgn,mpfr_number_p,mpfr_integer_p,mpfr_cmp_si,mpfr_rint,mpfr_round,mpfr_trunc,mpfr_ceil,mpfr_floor,mpfr_rnd_t,MPFR_RNDZ,MPFR_RNDU,MPFR_RNDD,MPFR_RNDA

# Sage 10.3's pxd omits roundeven and declares get_z void; the native MPFR
# header exports both with int return values, used here without adaptation.
from sage.libs.gmp.types cimport mpz_t
cdef extern from "mpfr.h":
    int mpfr_roundeven(mpfr_t,mpfr_t)
    int native_get_z "mpfr_get_z"(mpz_t,mpfr_t,mpfr_rnd_t)

def native_observer(str text,int precision,int target,int operation,int mode,long integer):
    cdef mpfr_t x,r
    cdef Integer z=Integer(0)
    cdef bytes encoded=text.encode('ascii')
    cdef int status
    cdef int rounding_code=[MPFR_RNDN,MPFR_RNDZ,MPFR_RNDU,MPFR_RNDD,MPFR_RNDA,-1][mode]
    cdef mpfr_rnd_t rounding=<mpfr_rnd_t>rounding_code
    mpfr_init2(x,precision);mpfr_init2(r,target)
    try:
        mpfr_set_str(x,encoded,10,MPFR_RNDN)
        if operation==0:
            return [bool(mpfr_nan_p(x)),bool(mpfr_inf_p(x)),bool(mpfr_number_p(x)),bool(mpfr_integer_p(x)),mpfr_sgn(x),mpfr_cmp_si(x,integer)]
        if operation==7:
            status=native_get_z(z.value,x,rounding)
            return [str(z),status]
        if operation==8:
            status=mpfr_rint(x,x,rounding)
            return state_frame(x,precision,0,status)
        if operation==1:status=mpfr_rint(r,x,rounding)
        elif operation==2:status=mpfr_roundeven(r,x)
        elif operation==3:status=mpfr_round(r,x)
        elif operation==4:status=mpfr_trunc(r,x)
        elif operation==5:status=mpfr_ceil(r,x)
        else:status=mpfr_floor(r,x)
        return state_frame(r,target,0,status)
    finally:mpfr_clear(x);mpfr_clear(r)


from sage.libs.mpfr cimport mpfr_set_exp

def native_rint_extreme(int sign,long exponent,int precision,int target,int mode):
    cdef mpfr_t x,r
    cdef Integer z=Integer(sign*3),mantissa=Integer(0)
    cdef int status,rounding_code=[MPFR_RNDN,MPFR_RNDZ,MPFR_RNDU,MPFR_RNDD,MPFR_RNDA,-1][mode]
    cdef long result_exponent=0
    mpfr_init2(x,precision);mpfr_init2(r,target)
    try:
        mpfr_set_z(x,z.value,MPFR_RNDN);mpfr_set_exp(x,exponent)
        status=mpfr_rint(r,x,<mpfr_rnd_t>rounding_code)
        kind='zero' if mpfr_zero_p(r) else 'finite'
        if kind=='finite':result_exponent=mpfr_get_z_exp(mantissa.value,r)+target
        return [status,kind,-1 if mpfr_signbit(r) else 1,str(abs(mantissa)),int(result_exponent),target]
    finally:mpfr_clear(x);mpfr_clear(r)


cdef extern from "mpfr.h":
    int native_set "mpfr_set"(mpfr_t,mpfr_t,mpfr_rnd_t)
    int native_frac "mpfr_frac"(mpfr_t,mpfr_t,mpfr_rnd_t)
    int native_cmp "mpfr_cmp"(mpfr_t,mpfr_t)

def native_copy_fraction(str text,int precision,int target,bint alias):
    cdef mpfr_t x,y,z
    cdef bytes encoded=text.encode('ascii')
    cdef int status,comparison
    mpfr_init2(x,precision);mpfr_init2(y,target);mpfr_init2(z,target)
    try:
        mpfr_set_str(x,encoded,10,MPFR_RNDN)
        if alias:
            status=native_set(x,x,MPFR_RNDN)
            copied=state_frame(x,precision,0,status);comparison=native_cmp(x,x)
            status=native_frac(x,x,MPFR_RNDN)
            return [copied,comparison,state_frame(x,precision,0,status)]
        status=native_set(y,x,MPFR_RNDN)
        copied=state_frame(y,target,0,status);comparison=native_cmp(x,y)
        status=native_frac(z,x,MPFR_RNDN)
        return [copied,comparison,state_frame(z,target,0,status)]
    finally:mpfr_clear(x);mpfr_clear(y);mpfr_clear(z)

cdef extern from "mpfr.h":
    int native_add "mpfr_add"(mpfr_t,mpfr_t,mpfr_t,mpfr_rnd_t)
    int native_mul "mpfr_mul"(mpfr_t,mpfr_t,mpfr_t,mpfr_rnd_t)
    long native_get_emin "mpfr_get_emin"()
    long native_get_emax "mpfr_get_emax"()

cdef compact_frame(mpfr_t x,int precision,int status):
    cdef Integer mantissa=Integer(0)
    cdef long exponent=0
    kind='nan' if mpfr_nan_p(x) else 'inf' if mpfr_inf_p(x) else 'zero' if mpfr_zero_p(x) else 'finite'
    if kind=='finite':exponent=mpfr_get_z_exp(mantissa.value,x)+precision
    return [status,kind,-1 if mpfr_signbit(x) else 1,str(abs(mantissa)),int(exponent),precision]

def native_binary(str left,str right,int lp,int rp,int target,int operation,int alias,long le=0,long re=0):
    cdef mpfr_t x,y,z
    cdef bytes l=left.encode('ascii'),r=right.encode('ascii')
    cdef int status
    mpfr_init2(x,lp);mpfr_init2(y,rp);mpfr_init2(z,target)
    try:
        mpfr_set_str(x,l,10,MPFR_RNDN);mpfr_set_str(y,r,10,MPFR_RNDN)
        if le:mpfr_set_exp(x,le)
        if re:mpfr_set_exp(y,re)
        mpfr_set_d(z,-0.0,MPFR_RNDN)
        if alias==1:
            status=native_add(x,x,y,MPFR_RNDN) if operation==0 else native_mul(x,x,y,MPFR_RNDN)
            return compact_frame(x,lp,status)
        if alias==2:
            status=native_add(y,x,y,MPFR_RNDN) if operation==0 else native_mul(y,x,y,MPFR_RNDN)
            return compact_frame(y,rp,status)
        status=native_add(z,x,y,MPFR_RNDN) if operation==0 else native_mul(z,x,y,MPFR_RNDN)
        return compact_frame(z,target,status)
    finally:mpfr_clear(x);mpfr_clear(y);mpfr_clear(z)

def native_exponent_limits():
    return [native_get_emin(),native_get_emax()]
