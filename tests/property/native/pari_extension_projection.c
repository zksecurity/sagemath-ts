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

static GEN raw_word(GEN a){GEN r=cgetg(lg(a),t_VECSMALL);r[1]=evalvarn(varn(a));for(long i=2;i<lg(a);i++)r[i]=itou(gel(a,i));return r;}
int main(void){
 long mode,code,n;static char sp[16000],seed[16000],st[2000000],sa[2000000],sb[2000000];pari_init(256000000,500000);
 while(scanf("%ld %ld %15999s %15999s %ld %1999999s %1999999s %1999999s",&mode,&code,sp,seed,&n,st,sa,sb)==8){
  pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
  pari_TRY{
   GEN p=gp_read_str(sp),T=gtopolyrev(gp_read_str(st),1),a=dense(gp_read_str(sa),0),b=dense(gp_read_str(sb),0),R;long op=code%10;
   pari_init_rand();setrand(gp_read_str(seed));
   if(mode==1){ulong q=itou(p);T=op==0?raw_word(T):ZX_to_Flx(T,q);a=ZXX_to_FlxX(a,q,1);b=ZXX_to_FlxX(b,q,1);if(code>=10)T=Flx_get_red(T,q);
    R=op==0?random_FlxqX(n,0,T,q):op==1?FlxqX_dotproduct(a,b,T,q):op==2?FlxqXn_mul(a,b,n,T,q):FlxqXn_sqr(a,n,T,q);
   }else{
    if(code>=10)T=FpX_get_red(T,p);
    R=op==0?random_FpXQX(n,0,T,p):op==1?FpXQX_dotproduct(a,b,T,p):op==2?FpXQXn_mul(a,b,n,T,p):FpXQXn_sqr(a,n,T,p);
   }
   printf("OK ");if(op==1)print_coefficient(mode==1?Flx_to_ZX(R):R);else print_outer(R,mode);pari_printf("|%Ps\n",getrand());
  }pari_ENDCATCH;
  set_avma(av);fflush(stdout);
 }
 pari_close();return 0;
}
