# cython: infer_types=False
"""Original binary splitting and unsigned rational atanh kernels."""
from cysignals.signals cimport sig_on, sig_off
from cypari2.types cimport GEN, pari_sp
from cypari2.paridecl cimport lg, signe, expo
cdef extern from "pari/pari.h":
    pari_sp get_avma()
    void set_avma(pari_sp)
# Installed 2.15.4 attempts allocation for a double-rounded-equal u/v ratio.
# Compile the bundled 2.18.1 get_nmax/atanhuu source with only symbol renaming.
# Installed precision macros translate the wrapper's word length to native bits.
cdef extern from *:
    r"""
    #include <pari/pari.h>
    #include <pari/paripriv.h>
static long
audit_get_nmax(double u, double v, long prec)
{
  double d = 2 * log2(((double)v) / u); /* can be 0 due to rounding */
  long nmax = -1;
  if (d)
  {
    d = ceil(prec2nbits(prec) / d);
    if (dblexpo(d) < BITS_IN_LONG) nmax = (long)d;
  }
  return nmax;
}
/* atanh(u/v) using binary splitting, 0 < u < v */
GEN
audit_atanhuu(ulong u, ulong v, long prec)
{
  GEN u2 = sqru(u), v2 = sqru(v);
  long i, nmax = audit_get_nmax((double)u, (double)v, prec);
  struct abpq_res R;
  struct abpq A;
  if (nmax < 0) pari_err_OVERFLOW("atanhuu");
  abpq_init(&A, nmax); /* nmax satisfies (2n+1) (v/u)^2n > 2^bitprec */
  A.a[0] = A.b[0] = gen_1;
  A.p[0] = utoipos(u);
  A.q[0] = utoipos(v);
  for (i = 1; i <= nmax; i++)
  {
    A.a[i] = gen_1;
    A.b[i] = utoipos((i<<1)+1);
    A.p[i] = u2;
    A.q[i] = v2;
  }
  abpq_sum(&R, 0, nmax, &A);
  return rdivii(R.T, mulii(R.B,R.Q),prec);
}
    """
    GEN audit_atanhuu(unsigned long,unsigned long,long)
def native_split(long op,unsigned long n,unsigned long a,long b):
    cdef GEN z
    cdef pari_sp saved
    cdef long j,sign,ex,words
    cdef unsigned long word
    cdef object mantissa
    sig_on();saved=get_avma()
    try:
        z=audit_atanhuu(n,a,(b+63)//64+2)
        sign=signe(z);ex=expo(z);words=lg(z);mantissa=0
        if sign:
            for j in range(2,words):
                word=<unsigned long>z[j];mantissa=(mantissa<<64)+word
        return [sign,str(ex),str(mantissa),(words-2)*64]
    finally:set_avma(saved);sig_off()
