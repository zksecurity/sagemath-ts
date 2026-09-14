/* Bundled PARI polynomial product/square ABI oracle. */
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
    long op;static char sp[16000],sa[2000000],sb[2000000];
    pari_init(256000000,500000);
    while(scanf("%ld %15999s %1999999s %1999999s",&op,sp,sa,sb)==4){
        pari_sp av=avma;
        pari_CATCH(CATCH_ALL){char *error=pari_err2str(pari_err_last());print_error(error);pari_free(error);}
        pari_TRY{
            GEN p=gp_read_str(sp),a=gtopolyrev(gp_read_str(sa),0),b=gtopolyrev(gp_read_str(sb),0),r;
            if(op==0)r=FpX_mul(a,b,p);
            else if(op==1)r=FpX_sqr(a,p);
            else if(op==2)r=ZX_mul(a,b);
            else if(op==3)r=ZX_sqr(a);
            else if(op==4){ulong q=itou(p);r=Flx_to_ZX(Flx_mul(ZX_to_Flx(a,q),ZX_to_Flx(b,q),q));}
            else if(op==5){ulong q=itou(p);r=Flx_to_ZX(Flx_sqr(ZX_to_Flx(a,q),q));}
            else if(op==6)r=ZX_mul(a,a);
            else r=FpX_mul(a,a,p);
            printf("OK ");print_poly(r);printf("\n");
        }pari_ENDCATCH;
        set_avma(av);fflush(stdout);
    }
    pari_close();return 0;
}
