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

int main(void){
 long mode,code;static char sp[16000],seed[16000],st[2000000],ss[2000000],sx[2000000];pari_init(256000000,500000);
 while(scanf("%ld %ld %15999s %15999s %1999999s %1999999s %1999999s",&mode,&code,sp,seed,st,ss,sx)==7){
  pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
  pari_TRY{
   GEN p=gp_read_str(sp),T=gtopolyrev(gp_read_str(st),1),S=dense(gp_read_str(ss),0),x=dense(gp_read_str(sx),0),R;
   pari_init_rand();setrand(gp_read_str(seed));
   if(mode==1){ulong q=itou(p);T=ZX_to_Flx(T,q);S=ZXX_to_FlxX(S,q,1);x=ZXX_to_FlxX(x,q,1);
    if(code>=100)T=Flx_get_red(T,q);if(code%100==10)S=mkvec2(FlxqX_invBarrett(S,T,q),S);else if(code%100==20)S=FlxqX_get_red(S,T,q);
    R=FlxqXQ_minpoly(x,S,T,q);
   }else{
    if(code>=100)T=FpX_get_red(T,p);if(code%100==10)S=mkvec2(FpXQX_invBarrett(S,T,p),S);else if(code%100==20)S=FpXQX_get_red(S,T,p);
    R=FpXQXQ_minpoly(x,S,T,p);
   }
   printf("OK ");print_outer(R,mode);pari_printf("|%Ps\n",getrand());
  }pari_ENDCATCH;
  set_avma(av);fflush(stdout);
 }
 pari_close();return 0;
}
