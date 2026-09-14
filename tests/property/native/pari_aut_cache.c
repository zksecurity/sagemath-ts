/* Exact compiled-PARI oracle for the public automorphism cache regressions. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>

static void polynomial(GEN p) {
  putchar('[');
  for (long i=2; i<lg(p); i++) {
    if (i>2) putchar(',');
    pari_printf("%Ps", gel(p,i));
  }
  putchar(']');
}
static void error_json(const char *s) {
  printf("ERROR \"");
  for (; *s; s++) {
    unsigned char c=*s;
    if (c=='"' || c=='\\') { putchar('\\'); putchar(c); }
    else if (c<32) printf("\\u%04x", c);
    else putchar(c);
  }
  puts("\"");
}
int main(void) {
  long op;
  static char sp[16000], sn[16000], sx[2000000], sa[2000000], st[2000000];
  pari_init(256000000, 500000);
  while (scanf("%ld %15999s %15999s %1999999s %1999999s %1999999s", &op, sp, sn, sx, sa, st)==6) {
    pari_sp av=avma;
    pari_CATCH(CATCH_ALL) {
      char *e=pari_err2str(pari_err_last());
      error_json(e); pari_free(e);
    } pari_TRY {
      GEN p=gp_read_str(sp), n=gp_read_str(sn);
      GEN x=gtopolyrev(gp_read_str(sx),0), a=gtopolyrev(gp_read_str(sa),0);
      GEN T=gtopolyrev(gp_read_str(st),0), R;
      if (op==0) R=FpXQ_autpow(x,itou(n),T,p);
      else if (op==1) R=FpXQ_autpowers(x,itos(n),T,p);
      else R=FpXQ_auttrace(mkvec2(x,a),itou(n),T,p);
      printf("OK ");
      if (op) {
        putchar('[');
        for (long i=1; i<lg(R); i++) {
          if (i>1) putchar(',');
          polynomial(gel(R,i));
        }
        putchar(']');
      } else polynomial(R);
      putchar('\n');
    } pari_ENDCATCH;
    set_avma(av); fflush(stdout);
  }
  pari_close();
  return 0;
}
