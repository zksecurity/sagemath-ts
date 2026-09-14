#include <math.h>
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

int main(void) {
 long mode,op;static char sp[16000],st[2000000],sP[2000000],sx[2000000],sV[2000000];
 pari_init(256000000,500000);
 while(scanf("%ld %ld %15999s %1999999s %1999999s %1999999s %1999999s",&mode,&op,sp,st,sP,sx,sV)==7){
  pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
  pari_TRY{
   GEN p=gp_read_str(sp),T,P,x,V=gp_read_str(sV),R;long cached=op/10;op%=10;
   if(mode>=2){
    T=packed(gp_read_str(st));P=mode==3?packed(gel(gp_read_str(sP),1)):binary(gp_read_str(sP));x=packed(gp_read_str(sx));
    for(long i=1;i<lg(V);i++)gel(V,i)=packed(gel(V,i));
    if(op==9)R=F2xq_powers(x,brent_kung_optpow(get_F2x_degree(T)-1,lgpol(P),1),T);
    else R=mode==3?(op?F2x_F2xqV_eval(P,V,T):F2x_F2xq_eval(P,x,T)):(op?F2xY_F2xqV_evalx(P,V,T):F2xY_F2xq_evalx(P,x,T));
   }else{
    T=gtopolyrev(gp_read_str(st),1);P=dense(gp_read_str(sP),0);x=gtopolyrev(gp_read_str(sx),1);
    for(long i=1;i<lg(V);i++)gel(V,i)=gtopolyrev(gel(V,i),1);
    if(mode==1){
     ulong q=itou(p);T=ZX_to_Flx(T,q);P=ZXX_to_FlxX(P,q,1);x=ZX_to_Flx(x,q);
     for(long i=1;i<lg(V);i++)gel(V,i)=ZX_to_Flx(gel(V,i),q);
     if(cached)T=Flx_get_red(T,q);
     R=op==9?Flxq_powers(x,brent_kung_optpow(get_Flx_degree(T)-1,lgpol(P),1),T,q):op?FlxY_FlxqV_evalx(P,V,T,q):FlxY_Flxq_evalx(P,x,T,q);
    }else{
     if(cached)T=FpX_get_red(T,p);
     R=op==9?FpXQ_powers(x,brent_kung_optpow(get_FpX_degree(T)-1,lgpol(P),1),T,p):op?FpXY_FpXQV_evalx(P,V,T,p):FpXY_FpXQ_evalx(P,x,T,p);
    }
   }
   printf("OK ");
   if(op==9){putchar('[');for(long i=1;i<lg(R);i++){if(i>1)putchar(',');if(mode==2)pari_printf("%Ps",ZX_Z_eval(F2x_to_ZX(gel(R,i)),gen_2));else print_coefficient(mode==1?Flx_to_ZX(gel(R,i)):gel(R,i));}putchar(']');}
   else if(mode==3)pari_printf("%Ps",ZX_Z_eval(F2x_to_ZX(R),gen_2));else print_outer(R,mode);putchar(10);
  }pari_ENDCATCH;
  set_avma(av);fflush(stdout);
 }
 pari_close();return 0;
}
