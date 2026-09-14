# Direct native FLINT oracle for the word kernels used by polynomial square roots.
from sage.libs.flint.ulong_extras cimport (
    n_sqrtmod, n_jacobi_unsigned, n_jacobi, n_is_square,
    n_preinvert_limb, n_powmod2_ui_preinv,
)
cdef extern from "flint/ulong_extras.h":
    int _n_jacobi_unsigned(unsigned long x, unsigned long y, unsigned int r) noexcept

def native_word(int op,a,unsigned long p,unsigned long e,unsigned int sign):
    if op==0:return str(n_sqrtmod(a,p))
    if op==1:return str(n_jacobi_unsigned(a,p))
    if op==2:return str(_n_jacobi_unsigned(a,p,sign))
    if op==3:return str(bool(n_is_square(a)))
    if op==4:return str(n_powmod2_ui_preinv(a,e,p,n_preinvert_limb(p)))
    if op==5:return str(n_preinvert_limb(p))
    if op==6:return str(n_jacobi(a,p))
    raise ValueError('unknown native word operation')
