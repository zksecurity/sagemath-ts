/* Compile the unmodified original, including its private fast stage.
 * This oracle audits the AArch64 Clang contraction profile documented in DESIGN.md.
 */
#if !defined(__aarch64__) || !defined(__clang__)
#error "The fast-stage transcript requires the audited AArch64 Clang profile"
#endif
#include "lll.c"
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

static double bits_double(GEN x)
{
  uint64_t u = itou(x);
  double d;
  memcpy(&d, &u, 8);
  return d;
}

static void double_bits(double d)
{
  uint64_t u;
  memcpy(&u, &d, 8);
  if (isnan(d)) printf("\"NaN\"");
  else printf("\"%llu\"", (unsigned long long)u);
}

static void scalar(long op, GEN v)
{
  double a = bits_double(gel(v, 1));
  int e;
  printf("[");
  if (op == 1) double_bits(fma(a, bits_double(gel(v, 2)), bits_double(gel(v, 3))));
  else if (op == 2) double_bits(ldexp(a, itos(gel(v, 2))));
  else {
    double d = frexp(a, &e);
    double_bits(d);
    printf(",%d", e);
  }
  printf("]\n");
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

int main(void)
{
  long op, m, n, dn, dd, en, ed, keep, track;
  static char data[4000000];
  pari_init(256000000, 500000);
  while (scanf("%ld %ld %ld %ld %ld %ld %ld %ld %ld %3999999s",
               &op, &m, &n, &dn, &dd, &en, &ed, &keep, &track, data) == 10) {
    pari_sp av = avma;
    pari_CATCH(CATCH_ALL) { printf("ERROR\n"); }
    pari_TRY {
      if (op) { scalar(op, gp_read_str(data)); }
      else {
        GEN B = read_matrix(gp_read_str(data), m, n), U = track ? matid(n) : NULL;
        long z = fplll_fast(&B, &U, (double)dn / dd, (double)en / ed, keep);
        printf("[%ld,", z);
        matrix_columns(B);
        printf(",");
        matrix_columns(U);
        printf("]\n");
      }
    }
    pari_ENDCATCH;
    set_avma(av);
    fflush(stdout);
  }
  pari_close();
  return 0;
}
