# cython: infer_types=False
"""Bundled exp1r_abs with the installed PARI precision ABI."""
from cysignals.signals cimport sig_on, sig_off
from cypari2.types cimport GEN, pari_sp
from cypari2.paridecl cimport lg, signe, expo
cdef extern from "pari/pari.h":
    pari_sp get_avma()
    void set_avma(pari_sp)
# Original excerpt SHA-256: 02b8f1c3012ed197886462a4063f0bc487843098a5d73bbc26ed816ff51bbf75
# PARI 2.15 lacks the bundled nonpositive-d and tiny-result branches.
# Only the symbol is renamed; installed macros translate word precision.
cdef extern from *:
    r"""
    #include <pari/pari.h>
    #include <pari/paripriv.h>
GEN
audit_exp1r_abs(GEN x)
{
  long l = realprec(x), a = expo(x), b = prec2nbits(l), L, i, n, m, B;
  GEN y, p2, X;
  pari_sp av;
  double d;

  if (b + a <= 0) return mpabs(x);

  y = cgetr(l); av = avma;
  B = b/3 + BITS_IN_LONG + (BITS_IN_LONG*BITS_IN_LONG)/ b;
  d = a/2.; m = (long)(d + sqrt(d*d + B)); /* >= 0 */
  if (m < (-a) * 0.1) m = 0; /* not worth it */
 /* Multiplication is quadratic in this range (l is small, otherwise we
  * use logAGM + Newton). Set Y = 2^(-e-a) x, compute truncated series
  * sum_{k <= n} Y^k/k!: this costs roughly
  *    m b^2 + sum_{k <= n} (k e + BITS_IN_LONG)^2
  * bit operations with n ~ b/e, |x| <  2^(1+a), |Y| < 2^(1-e) , m = e+a and
  * b bits of accuracy needed, so
  *    B := (b / 3 + BITS_IN_LONG + BITS_IN_LONG^2 / b) ~ m(m-a)
  * we want b ~ 3 m (m-a) or m~b+a hence
  *     m = min( a/2 + sqrt(a^2/4 + B),  b + a )
  * NB: e ~ (b/3)^(1/2) as b -> oo
  *
  * Truncate the sum at k = n (>= 1), the remainder is
  *   sum_{k >= n+1} Y^k / k! < Y^(n+1) / (n+1)! (1-Y) < Y^(n+1) / n!
  * We want Y^(n+1) / n! <= Y 2^-b, hence -n log_2 |Y| + log_2 n! >= b
  *   log n! ~ (n + 1/2) log(n+1) - (n+1) + log(2Pi)/2,
  * error bounded by 1/6(n+1) <= 1/12. Finally, we want
  * n (-1/log(2) -log_2 |Y| + log_2(n+1)) >= b  */
  d = m-dbllog2(x)-1/M_LN2; /* ~ -log_2 Y - 1/log(2) */
  while (d <= 0) { d++; m++; } /* d < 0 can occur from expm1 */
  L = l + nbits2extraprec(m);
  b += m;
  n = (long)(b / d); /* > 0 */
  if (n == 1)
    n = (long)(b / (d + log2((double)n+1))); /* log ~ const in small ranges */
  while (n*(d+log2((double)n+1)) < b) n++; /* expect few corrections */

  X = rtor(x,L); shiftr_inplace(X, -m); setsigne(X, 1);
  if (n == 1) p2 = X;
  else
  {
    long s = 0, l1 = nbits2prec((long)(d + n + 16));
    GEN unr = real_1(L);
    pari_sp av2;

    p2 = cgetr(L); av2 = avma;
    for (i=n; i>=2; i--, set_avma(av2))
    { /* compute X^(n-1)/n! + ... + X/2 + 1 */
      GEN p1, p3;
      setprec(X,l1); p3 = divru(X,i);
      l1 += nbits2extraprec(dvmdsBIL(s - expo(p3), &s)<<TWOPOTBITS_IN_LONG);
      if (l1>L) l1=L;
      setprec(unr,l1); p1 = addrr_sign(unr,1, i == n? p3: mulrr(p3,p2),1);
      setprec(p2,l1); affrr(p1,p2); /* p2 <- 1 + (X/i)*p2 */
    }
    setprec(X,L); p2 = mulrr(X,p2);
  }

  B = prec2nbits(L);
  for (i = 1; i <= m; i++)
  {
    if (realprec(p2) > L) setprec(p2,L);
    if (expo(p2) < -B)
      shiftr_inplace(p2, 1); /* 2 + p2 ~ 2 and may blow up accuracy */
    else
      p2 = mulrr(p2, addsr(2,p2));
  }
  affrr_fixlg(p2,y); return gc_const(av,y);
}

    """
    GEN audit_exp1r_abs(GEN)
from cypari2.paridecl cimport cgetr, evalexpo, evalsigne

def native_exp1(long p,long e,m,long s):
    cdef GEN x,z
    cdef pari_sp saved
    cdef long i,words,sign,ex
    cdef unsigned long word
    sig_on();saved=get_avma()
    try:
        x=cgetr(p//64+2);x[1]=evalsigne(s)|evalexpo(e)
        for i in range(2,p//64+2):
            word=(m>>(p-64*(i-1)))&((1<<64)-1);x[i]=<long>word
        z=audit_exp1r_abs(x)
        words=lg(z);sign=signe(z);ex=expo(z);mantissa=0
        if sign:
            for i in range(2,words):
                word=<unsigned long>z[i];mantissa=(mantissa<<64)+word
        return [sign,str(ex),str(mantissa),(words-2)*64]
    finally:set_avma(saved);sig_off()
