/* Bundled PARI oracle for LLL callers; matrices travel as flat row-major vectors.
 * Output matrices retain PARI's column order. Precision uses the bundled bit ABI.
 */
#include "pari.h"
#include "paripriv.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static void error(const char *s)
{
  printf("ERROR \"");
  for (; *s; s++) {
    unsigned char c = *s;
    if (c == '"' || c == '\\') { putchar('\\'); putchar(c); }
    else if (c < 32) printf("\\u%04x", c);
    else putchar(c);
  }
  printf("\"\n");
}

static void value(GEN x)
{
  if (!x) { printf("null"); return; }
  if (typ(x) == t_INT) { pari_printf("\"%Ps\"", x); return; }
  printf("[");
  for (long i = 1; i < lg(x); i++) {
    if (i > 1) printf(",");
    value(gel(x, i));
  }
  printf("]");
}

int main(void)
{
  long op, n, flag;
  static char data[400000];
  pari_init(256000000, 500000);
  while (scanf("%ld %ld %ld %399999s", &op, &n, &flag, data) == 4) {
    pari_sp av = avma;
    pari_CATCH(CATCH_ALL) {
      char *s = pari_err2str(pari_err_last());
      error(s);
      pari_free(s);
    } pari_TRY {
      GEN v = gp_read_str(data);
      if (op == 5) {
        uint64_t bits = itou(gel(v, 1));
        double d;
        memcpy(&d, &bits, sizeof(d));
        GEN f = algdep(dbltor(d), n);
        value(typ(f) == t_INT ? mkvec(f) : RgX_to_RgV(f, lg(f) - 2));
      } else {
        long rows = n ? (lg(v) - 1 - (op == 3 ? 2 : op == 7 ? 1 : 0)) / n : 0;
        GEN M = cgetg(n + 1, t_MAT), U = NULL, H;
        for (long j = 1; j <= n; j++) {
          gel(M, j) = cgetg(rows + 1, t_COL);
          for (long i = 1; i <= rows; i++) gcoeff(M, i, j) = gel(v, (i - 1) * n + j);
        }
        if (op == 0) value(lllgramint(M));
        else if (op == 1) {
          H = ZM_hnflll(M, &U, flag);
          value(mkvec2(H, U));
        } else if (op == 2) {
          H = ZM_snf_group(M, NULL, &U);
          value(mkvec2(H, U));
        } else if (op == 3)
          value(ZM_lll(M, gtodouble(gdiv(gel(v, lg(v) - 2), gel(v, lg(v) - 1))),
                       LLL_IM | LLL_GRAM));
        else if (op == 4) value(ZM_hnflll(M, NULL, flag));
        else if (op == 6) value(qfminim0(M, NULL, NULL, 1, DEFAULTPREC));
        else if (op == 7) {
          GEN counts = qfrep0(M, gel(v, lg(v) - 1), flag);
          printf("[");
          for (long i = 1; i < lg(counts); i++) {
            if (i > 1) printf(",");
            printf("\"%ld\"", counts[i]);
          }
          printf("]");
        }
      }
      printf("\n");
    } pari_ENDCATCH;
    set_avma(av);
    fflush(stdout);
  }
  pari_close();
  return 0;
}
