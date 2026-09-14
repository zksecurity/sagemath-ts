#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static GEN dense(GEN v,long var){GEN x=cgetg(lg(v)+1,t_POL);x[1]=evalsigne(1)|evalvarn(var);for(long i=1;i<lg(v);i++)gel(x,i+1)=typ(gel(v,i))==t_VEC?gtopolyrev(gel(v,i),1):gel(v,i);return FpXX_renormalize(x,lg(x));}
static GEN packed(GEN n){if(!signe(n))return pol0_F2x(evalvarn(1));long d=expi(n);GEN v=cgetg(d+2,t_VEC);for(long i=0;i<=d;i++)gel(v,i+1)=bittest(n,i)?gen_1:gen_0;return ZX_to_F2x(gtopolyrev(v,1));}
static GEN binary(GEN v){GEN x=cgetg(lg(v)+1,t_POL);x[1]=evalsigne(1)|evalvarn(0);for(long i=1;i<lg(v);i++)gel(x,i+1)=packed(gel(v,i));return F2xX_renormalize(x,lg(x));}
static void print_coefficient(GEN c){if(typ(c)==t_INT){pari_printf("%Ps",c);return;}printf("[");for(long j=2;j<lg(c);j++){if(j>2)printf(",");pari_printf("%Ps",gel(c,j));}printf("]");}
static void print_outer(GEN R,long mode){printf("[");for(long i=2;i<lg(R);i++){if(i>2)printf(",");if(mode==2)pari_printf("%Ps",ZX_Z_eval(F2x_to_ZX(gel(R,i)),gen_2));else print_coefficient(mode==1?Flx_to_ZX(gel(R,i)):gel(R,i));}printf("]");}


static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}

int main(void){long op;static char sp[16000],sc[16000],sa[2000000],sb[2000000];pari_init(256000000,500000);
 while(scanf("%ld %15999s %15999s %1999999s %1999999s",&op,sp,sc,sa,sb)==5){pari_sp av=avma;
 pari_CATCH(CATCH_ALL){char*e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
 pari_TRY{GEN p=gp_read_str(sp),c=gp_read_str(sc),a=gtopolyrev(gp_read_str(sa),0),b=gtopolyrev(gp_read_str(sb),0),R;
 R=op==0?FpX_red(a,p):op==1?FpX_add(a,b,p):op==2?FpX_sub(a,b,p):op==3?FpX_neg(a,p):FpX_Fp_mul(a,c,p);
 printf("OK ");print_coefficient(R);putchar(10);
 }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }pari_close();return 0;}
