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
int main(void){long mode,op;static char sp[16000],st[2000000],sa[2000000],sb[2000000];pari_init(256000000,500000);
while(scanf("%ld %ld %15999s %1999999s %1999999s %1999999s",&mode,&op,sp,st,sa,sb)==6){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char *e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
pari_TRY{GEN p=gp_read_str(sp),T,A,B,R,u=NULL,v=NULL;GEN *pu=op==1?&u:NULL;
if(mode==2){T=packed(gp_read_str(st));A=binary(gp_read_str(sa));B=binary(gp_read_str(sb));R=op==0?F2xqX_gcd(A,B,T):op==2?F2xqX_halfgcd(A,B,T):F2xqX_extgcd(A,B,T,pu,&v);}
else{T=gtopolyrev(gp_read_str(st),1);A=dense(gp_read_str(sa),0);B=dense(gp_read_str(sb),0);
if(mode==1){ulong q=itou(p);T=ZX_to_Flx(T,q);A=ZXX_to_FlxX(A,q,1);B=ZXX_to_FlxX(B,q,1);R=op==0?FlxqX_gcd(A,B,T,q):op==2?FlxqX_halfgcd(A,B,T,q):FlxqX_extgcd(A,B,T,q,pu,&v);}
else {R=op==0?FpXQX_gcd(A,B,T,p):op==2?FpXQX_halfgcd(A,B,T,p):FpXQX_extgcd(A,B,T,p,pu,&v);}}
printf("OK ");if(op==0)print_outer(R,mode);else if(op==2){printf("[");for(long i=1;i<=2;i++){if(i>1)printf(",");printf("[");for(long j=1;j<=2;j++){if(j>1)printf(",");print_outer(gcoeff(R,i,j),mode);}printf("]");}printf("]");}else{printf("[");print_outer(R,mode);if(op==1){printf(",");print_outer(u,mode);}printf(",");print_outer(v,mode);printf("]");}printf("\n");
}pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
