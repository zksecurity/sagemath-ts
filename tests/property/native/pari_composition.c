/* Original bundled PARI modular composition and supporting schedule oracle. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void print_poly(GEN f){printf("[");for(long i=2;i<lg(f);i++){if(i>2)printf(",");pari_printf("%Ps",gel(f,i));}printf("]");}
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}

int main(void){
    long op;static char sp[16000],sn[100],sq[1000000],sx[1000000],st[1000000];
    pari_init(256000000,500000);
    while(scanf("%ld %15999s %99s %999999s %999999s %999999s",&op,sp,sn,sq,sx,st)==6){
        pari_sp av=avma;
        pari_CATCH(CATCH_ALL){char *error=pari_err2str(pari_err_last());print_error(error);pari_free(error);}
        pari_TRY{
            GEN p=gp_read_str(sp),ni=gp_read_str(sn),Q=gtopolyrev(gp_read_str(sq),0),x=gtopolyrev(gp_read_str(sx),0),T=gtopolyrev(gp_read_str(st),0),r;
            long n=op==5?0:itos(ni);
            if(op==0)r=FpX_FpXQ_eval(Q,x,T,p);
            else if(op==1)r=FpX_FpXQV_eval(Q,FpXQ_powers(x,n,T,p),T,p);
            else if(op==2||op==3){
                ulong pp=itou(p);GEN q=ZX_to_Flx(Q,pp),xx=ZX_to_Flx(x,pp),tt=ZX_to_Flx(T,pp);
                r=Flx_to_ZX(op==2?Flx_Flxq_eval(q,xx,tt,pp):Flx_FlxqV_eval(q,Flxq_powers(xx,n,tt,pp),tt,pp));
            }
            else if(op==8)r=FpX_FpXQV_eval(Q,cgetg(1,t_VEC),T,p);
            else if(op==9){ulong pp=itou(p);r=Flx_to_ZX(Flx_FlxqV_eval(ZX_to_Flx(Q,pp),cgetg(1,t_VEC),ZX_to_Flx(T,pp),pp));}
            else if(op==4)r=stoi(brent_kung_optpow(itos(p),n,itos(polcoef_i(Q,0,-1))));
            else if(op==5)r=FpXQ_auttrace(mkvec2(Q,x),itou(ni),T,p);
            else if(op==6)r=FpX_mul(Q,x,p);
            else r=FpX_Fp_mul(Q,polcoef_i(x,0,-1),p);
            printf("OK ");
            if(op==4)pari_printf("%Ps",r);
            else if(op==5){printf("[");print_poly(gel(r,1));printf(",");print_poly(gel(r,2));printf("]");}
            else print_poly(r);
            printf("\n");
        }pari_ENDCATCH;
        set_avma(av);fflush(stdout);
    }
    pari_close();return 0;
}
