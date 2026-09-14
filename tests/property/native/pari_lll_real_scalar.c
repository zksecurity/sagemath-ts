/* Native roundr_safe and abscmprr on exact real mantissa records. */
#include "pari.h"
#include "paripriv.h"
#ifndef LONG_IS_64BIT
#error "Real LLL scalar oracle requires the port 64-bit PARI ABI"
#endif
#include <stdio.h>

static GEN real_value(long p, long e, long s, char *m)
{
  GEN x = itor(gp_read_str(m), p);
  setsigne(x, s);
  setexpo(x, e);
  return x;
}

int main(void)
{
  long op, p, e, s, q, f, t;
  static char m[20000], n[20000];
  pari_init(256000000, 500000);
  while (scanf("%ld %ld %ld %19999s %ld %ld %ld %19999s %ld",
               &op, &p, &e, m, &s, &q, &f, n, &t) == 9) {
    pari_sp av = avma;
    pari_CATCH(CATCH_ALL) { printf("ERROR\n"); }
    pari_TRY {
      GEN x = real_value(p, e, s, m), y = real_value(q, f, t, n);
      if (op == 0) pari_printf("%Ps\n", roundr_safe(x));
      else printf("%d\n", abscmprr(x, y));
    }
    pari_ENDCATCH;
    set_avma(av);
    fflush(stdout);
  }
  pari_close();
  return 0;
}
