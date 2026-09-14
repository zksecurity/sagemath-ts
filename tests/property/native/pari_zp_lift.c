#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}



static void poly(GEN f){printf("[");for(long i=2;i<lg(f);i++){if(i>2)printf(",");pari_printf("%Ps",gel(f,i));}printf("]");}
int main(void){long op,e;static char sf[2000000],sa[2000000],st[2000000],sq[2000000],sp[30000];pari_init(256000000,500000);
while(scanf("%ld %1999999s %1999999s %1999999s %1999999s %29999s %ld",&op,sf,sa,st,sq,sp,&e)==7){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{GEN f=gtopolyrev(gp_read_str(sf),0),a=gp_read_str(sa),t=gtopolyrev(gp_read_str(st),0),flat=gp_read_str(sq),p=gp_read_str(sp),z,Q;long count=0;
for(long i=1;i<lg(flat);){long n=itos(gel(flat,i));i+=n+1;count++;}Q=cgetg(count+1,t_VEC);
for(long i=1,j=1;i<lg(flat);j++){long n=itos(gel(flat,i++));GEN v=cgetg(n+1,t_VEC);for(long k=1;k<=n;k++)gel(v,k)=gel(flat,i++);gel(Q,j)=gtopolyrev(v,0);}
/* Match the initial native evaluation/inverse phases, but stop before
 * unchecked exact division when the root precondition is violated. */
if(op==5 && e>1){GEN S=gtopolyrev(a,0),TT=FpX_get_red(t,powiu(p,e)),q=sqri(p),Tq=FpXT_red(TT,q),Tq2=FpXT_red(Tq,p);
(void)FpXQ_inv(FpX_FpXQ_eval(FpX_deriv(f,p),S,Tq2,p),Tq2,p);
GEN value=FpX_FpXQ_eval(FpX_red(f,q),S,Tq,q);
for(long i=2;i<lg(value);i++)if(signe(modii(gel(value,i),p))){printf("GUARD Hensel lifting requires exact polynomial division\n");goto done;}}
switch(op){case 0:z=ZpX_liftroot(f,gel(a,1),p,e);break;case 1:z=ZpX_liftroots(f,a,p,e);break;case 2:z=ZpX_roots(f,p,e);break;case 3:z=ZpX_liftfact(f,Q,powiu(p,e),p,e);break;case 4:z=bezout_lift_fact(f,Q,p,e);break;default:z=ZpX_ZpXQ_liftroot(f,gtopolyrev(a,0),t,p,e);}
printf("OK ");if(op==0)pari_printf("%Ps",z);else if(op==5)poly(z);else{printf("[");for(long i=1;i<lg(z);i++){if(i>1)printf(",");if(op==3||op==4)poly(gel(z,i));else pari_printf("%Ps",gel(z,i));}printf("]");}printf("\n");done:;}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
