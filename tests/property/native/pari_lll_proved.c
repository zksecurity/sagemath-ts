/* Original arbitrary-precision PARI LLL stage, with optional Gram data. */
#include "lll.c"
#ifndef LONG_IS_64BIT
#error "Real LLL oracle requires the port 64-bit PARI ABI"
#endif
#include <stdio.h>

static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}

static GEN read_matrix(GEN v, long m, long n)
{
  GEN A = cgetg(n + 1, t_MAT);
  for (long j = 1; j <= n; j++) {
    gel(A, j) = cgetg(m + 1, t_COL);
    for (long i = 1; i <= m; i++) gcoeff(A, i, j) = gel(v, (i - 1) * n + j);
  }
  return A;
}

static void matrix_columns(GEN A)
{
  if (!A) { printf("null"); return; }
  printf("[");
  for (long j = 1; j < lg(A); j++) {
    if (j > 1) printf(",");
    printf("[");
    for (long i = 1; i < lg(gel(A, j)); i++) {
      if (i > 1) printf(",");
      pari_printf("\"%Ps\"", gcoeff(A, i, j));
    }
    printf("]");
  }
  printf("]");
}

static void norms(GEN r)
{
  if (!r) { printf("null"); return; }
  printf("[");
  for (long i = 1; i < lg(r); i++) {
    GEN x = gel(r, i);
    if (i > 1) printf(",");
    GEN mantissa = gen_0;
    if (signe(x)) for (long j=2;j<lg(x);j++) mantissa=addii(shifti(mantissa,64),utoi((ulong)x[j]));
    printf("[%ld,\"%ld\",",signe(x),expo(x));
    pari_printf("\"%Ps\",%ld]",mantissa,(lg(x)-2)*64);
  }
  printf("]");
}

int main(void)
{
  long mode, m, n, dn, dd, en, ed, keep, track, want, precision;
  static char data[4000000];
  pari_init(256000000, 500000);
  while (scanf("%ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %3999999s",
               &mode, &m, &n, &dn, &dd, &en, &ed, &keep, &track, &want, &precision, data) == 12) {
    pari_sp av = avma;
    pari_CATCH(CATCH_ALL) {
      char *error = pari_err2str(pari_err_last());
      print_error(error);
      pari_free(error);
    }
    pari_TRY {
      GEN B = read_matrix(gp_read_str(data), m, n), G = NULL;
      GEN U = track ? matid(n) : NULL, r = NULL;
      if (mode) {
        G = zeromatcopy(n, n);
        for (long j = 1; j <= n; j++)
          for (long i = 1; i <= n; i++)
            gcoeff(G, i, j) = ZV_dotproduct(gel(B, i), gel(B, j));
      }
      if (mode == 2) B = NULL;
      long z = fplll(&G, &B, &U, want ? &r : NULL,
                         (double)dn / dd, (double)en / ed, keep, precision);
      printf("[%ld,", z);
      matrix_columns(G);
      printf(",");
      matrix_columns(B);
      printf(",");
      matrix_columns(U);
      printf(",");
      norms(r);
      printf("]\n");
    }
    pari_ENDCATCH;
    set_avma(av);
    fflush(stdout);
  }
  pari_close();
  return 0;
}
