from sage.libs.m4ri cimport mzd_t,mzd_init,mzd_free,mzd_row,mzd_mul,mzd_mul_m4rm,mzd_mul_naive,mzd_add,mzd_init_window,mzd_free_window,m4ri_word
from libc.stdlib cimport malloc,free
from cysignals.signals cimport sig_on,sig_off
cdef extern from "m4ri/m4ri.h":
    void mzd_make_table(const mzd_t*,int,int,int,mzd_t*,int*)

cdef mzd_t* packed_matrix(int nr,int nc,values) except NULL:
    cdef mzd_t* a=mzd_init(nr,nc)
    cdef int i,j
    cdef m4ri_word* row
    if nc==0:return a
    for i in range(nr):
        row=mzd_row(a,i)
        for j in range(a.width):row[j]=int((values[i]>>(64*j))&((1<<64)-1))
    return a

cdef packed_rows(mzd_t* a):
    cdef int i,j
    cdef m4ri_word* row
    cdef object value
    cdef object one=1
    out=[]
    if a.ncols==0:return [str(a.nrows),str(a.ncols),['0']*a.nrows]
    for i in range(a.nrows):
        row=mzd_row(a,i);value=0
        for j in range(a.width):value|=int(row[j])<<(64*j)
        out.append(str(value&((one<<a.ncols)-1)))
    return [str(a.nrows),str(a.ncols),out]

def native_m4ri(int method,int m,int k,int n,aa,bb,int parameter,bint square):
    cdef mzd_t *a=packed_matrix(m,k,aa)
    cdef mzd_t *b=NULL
    cdef mzd_t *c=NULL
    cdef int *lookup=NULL
    cdef bint window=False
    try:
        if method==4:
            c=mzd_init_window(a,parameter,64,m,k);window=True
        elif method==5:
            c=mzd_init(1<<parameter,k)
            lookup=<int*>malloc((1<<parameter)*sizeof(int))
            if lookup==NULL:raise MemoryError()
            mzd_make_table(a,0,0,parameter,c,lookup)
            return [packed_rows(c),[int(lookup[i]) for i in range(1<<parameter)]]
        else:
            b=a if square else packed_matrix(m if method==3 else k,k if method==3 else n,bb)
            if m==0 or k==0 or (method!=3 and n==0):
                c=mzd_init(m,k if method==3 else n)
                return packed_rows(c)
            sig_on()
            try:
                if method==0:c=mzd_mul(NULL,a,b,parameter)
                elif method==1:c=mzd_mul_m4rm(NULL,a,b,parameter)
                elif method==2:c=mzd_mul_naive(NULL,a,b)
                else:c=mzd_add(NULL,a,b)
            finally:sig_off()
        return packed_rows(c)
    finally:
        if lookup!=NULL:free(lookup)
        if c!=NULL:
            if window:mzd_free_window(c)
            else:mzd_free(c)
        if b!=NULL and b!=a:mzd_free(b)
        mzd_free(a)

from sage.libs.m4ri cimport mzd_submatrix
def native_submatrix(int m,int n,rows,int r0,int c0,int r1,int c1):
    cdef mzd_t* a=packed_matrix(m,n,rows)
    cdef mzd_t* b=NULL
    try:
        # Sage skips native extraction of empty slices; unaligned zero-width
        # native destinations do not have a valid row word to dereference.
        if r0==r1 or c0==c1:b=mzd_init(r1-r0,c1-c0)
        else:b=mzd_submatrix(NULL,a,r0,c0,r1,c1)
        return packed_rows(b)
    finally:
        if b!=NULL:mzd_free(b)
        mzd_free(a)
