/* Bundled PARI polynomial GCD, Bezout and half-GCD oracle. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void print_poly(GEN f){printf("[");for(long i=2;i<lg(f);i++){if(i>2)printf(",");pari_printf("%Ps",gel(f,i));}printf("]");}
static void print_matrix(GEN M){printf("[");for(long i=1;i<=2;i++){if(i>1)printf(",");printf("[");for(long j=1;j<=2;j++){if(j>1)printf(",");print_poly(gcoeff(M,i,j));}printf("]");}printf("]");}
struct audit_power {GEN p;ulong squares,multiplies;};
static GEN audit_square(void *data,GEN x){struct audit_power *s=data;s->squares++;return modii(addis(mulis(x,3),1),s->p);}
static GEN audit_multiply(void *data,GEN x,GEN y){struct audit_power *s=data;s->multiplies++;return modii(addis(addii(mulis(x,5),mulis(y,7)),1),s->p);}
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}

int main(void){long op;static char sp[16000],sa[2000000],sb[2000000];pari_init(256000000,500000);
 while(scanf("%ld %15999s %1999999s %1999999s",&op,sp,sa,sb)==4){pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *error=pari_err2str(pari_err_last());print_error(error);pari_free(error);}
  pari_TRY{GEN p=gp_read_str(sp),a=gtopolyrev(gp_read_str(sa),0),b=gtopolyrev(gp_read_str(sb),0),r,u=NULL,v=NULL;
   if(op==15){struct audit_power state={p,0,0};GEN e=lgpol(b)?gel(b,2):gen_0;GEN value=gen_pow_i(lgpol(a)?gel(a,2):gen_0,e,&state,audit_square,audit_multiply);r=mkvec3(value,utoi(state.squares),utoi(state.multiplies));}
   else if(op==14){GEN e=lgpol(b)?gel(b,2):gen_0;r=FpXQ_pow(a,e,RgX_shift_shallow(b,-1),p);}
   else if(op==11)r=FpX_split_part(a,p);
   else if(op==12)r=stoi(FpX_nbroots(a,p));
   else if(op==9)r=FpX_halfgcd(a,b,p);
   else if(op==10){ulong q=itou(p);r=FlxM_to_ZXM(Flx_halfgcd(ZX_to_Flx(a,q),ZX_to_Flx(b,q),q));}
   else if(op==0)r=FpX_gcd(a,b,p);
   else if(op==1)r=FpX_extgcd(a,b,p,&u,&v);
   else if(op==2)r=FpX_halfgcd_all(a,b,p,&u,&v);
   else if(op==6)r=FpXQ_inv(a,b,p);
   else if(op==7||op==8){GEN second=FpX_div(a,b,p),factors=mkvec2(b,second);r=op==7?ZpX_liftfact(a,factors,powiu(p,3),p,3):bezout_lift_fact(a,factors,p,3);}
   else{ulong q=itou(p);a=ZX_to_Flx(a,q);b=ZX_to_Flx(b,q);
    if(op==3)r=Flx_gcd(a,b,q);
    else if(op==4)r=Flx_extgcd(a,b,q,&u,&v);
    else r=Flx_halfgcd_all(a,b,q,&u,&v);
    if(op==5)r=FlxM_to_ZXM(r);else r=Flx_to_ZX(r);
    if(u)u=Flx_to_ZX(u);if(v)v=Flx_to_ZX(v);
   }
   printf("OK ");if(u)printf("[");if(op==15){printf("[");for(long i=1;i<lg(r);i++){if(i>1)printf(",");pari_printf("%Ps",gel(r,i));}printf("]");}
   else if(op==12)pari_printf("%Ps",r);else if(op==2||op==5||op==9||op==10)print_matrix(r);else if(op==7||op==8){printf("[");for(long i=1;i<lg(r);i++){if(i>1)printf(",");print_poly(gel(r,i));}printf("]");}else print_poly(r);
   if(u){printf(",");print_poly(u);printf(",");print_poly(v);printf("]");}printf("\n");
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }pari_close();return 0;}
