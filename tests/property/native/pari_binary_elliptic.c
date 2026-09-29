/* Direct bundled PARI binary-curve kernels, with packed coefficient transport. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static GEN packed(GEN n) {
  if (!signe(n)) return pol0_F2x(evalvarn(1));
  long d=expi(n); GEN v=cgetg(d+2,t_VEC);
  for(long i=0;i<=d;i++) gel(v,i+1)=bittest(n,i)?gen_1:gen_0;
  return ZX_to_F2x(gtopolyrev(v,1));
}
static GEN point(GEN v) { return lg(v)==1?ellinf():mkvec2(packed(gel(v,1)),packed(gel(v,2))); }
static void coefficient(GEN x) { if(!x) printf("null"); else pari_printf("%Ps",ZX_Z_eval(F2x_to_ZX(x),gen_2)); }
static void printpoint(GEN P) {
  if(ell_is_inf(P)) {printf("[]");return;}
  printf("[");coefficient(gel(P,1));printf(",");coefficient(gel(P,2));printf("]");
}
static void error(const char *s) {
  printf("ERROR \"");
  for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
  printf("\"\n");
}
int main(void) {
  char *line=NULL;size_t capacity=0;pari_init(64000000,500000);
  while(getline(&line,&capacity,stdin)>0) {
    pari_sp av=avma;
    pari_CATCH(CATCH_ALL) {char *s=pari_err2str(pari_err_last());error(s);pari_free(s);}
    pari_TRY {
      GEN v=gp_read_str(line),T=packed(gel(v,1)),a=gel(v,2),P=point(gel(v,3)),Q,n=gel(v,5),ch=gel(v,6),R;
      long op=itos(gel(v,7));int paired=op>=10;
      if(paired) {
        long mode=itos(gel(a,1));GEN a3=packed(gel(a,2)),a4=packed(gel(a,3));
        if(!mode) {
          GEN x=gel(P,1),y=gel(P,2),x2;
          if(!lgpol(x)) x=pol1_F2x(T[1]);
          P=mkvec2(x,y);x2=F2xq_sqr(x,T);
          a=F2xq_div(F2x_add(F2x_add(F2xq_sqr(y,T),F2xq_mul(x,y,T)),F2x_add(F2xq_mul(x2,x,T),pol1_F2x(T[1]))),x2,T);
        } else a=mkvec3(a3,a4,F2xq_inv(a3,T));
        Q=F2xqE_mul(P,gel(v,4),a,T);op-=10;
      } else {
        Q=point(gel(v,4));
        a=typ(a)==t_INT?packed(a):mkvec3(packed(gel(a,1)),packed(gel(a,2)),packed(gel(a,3)));
      }
      if(op==7)R=F2xq_invsafe(gel(P,1),T);
      else if(op==8)R=F2xq_inv(gel(P,1),T);
      else if(op==9)R=F2xq_div(gel(P,1),gel(P,2),T);
      else if(op==0)R=F2xqE_add(P,Q,a,T);
      else if(op==1)R=F2xqE_dbl(P,a,T);
      else if(op==2)R=F2xqE_neg(P,a,T);
      else if(op==3)R=F2xqE_sub(P,Q,a,T);
      else if(op==4)R=F2xqE_mul(P,n,a,T);
      else {GEN c=cgetg(5,t_VEC);for(long i=1;i<=4;i++)gel(c,i)=packed(gel(ch,i));R=op==5?F2xqE_changepoint(P,c,T):F2xqE_changepointinv(P,c,T);}
      printf("OK ");
      if(paired){printf("[");printpoint(Q);printf(",");printpoint(R);printf("]");}
      else if(op>=7)coefficient(R);else printpoint(R);
      printf("\n");
    } pari_ENDCATCH;
    set_avma(av);fflush(stdout);
  }
  free(line);pari_close();return 0;
}
