/* Isolated original PARI quotient-polynomial helper oracle. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void print_poly(GEN f) {
    printf("[");
    for(long i=2;i<lg(f);i++){if(i>2)printf(",");pari_printf("%Ps",gel(f,i));}
    printf("]");
}
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}

int main(void) {
    long op,n;char sp[16000],sx[16000],st[16000];
    pari_init(128000000,500000);
    while(scanf("%ld %15999s %ld %15999s %15999s",&op,sp,&n,sx,st)==5) {
        pari_sp av=avma;
        pari_CATCH(CATCH_ALL) {
            char *error=pari_err2str(pari_err_last());print_error(error);pari_free(error);
        } pari_TRY {
            GEN p=gp_read_str(sp),x=gtopolyrev(gp_read_str(sx),0),T=gtopolyrev(gp_read_str(st),0),r;
            if(op==0)r=FpXQ_powers(x,n,T,p);
            else if(op==1)r=FpXQ_autpow(x,(ulong)n,T,p);
            else r=FpXQ_autpowers(x,n,T,p);
            printf("OK ");
            if(op==1)print_poly(r);
            else {printf("[");for(long i=1;i<lg(r);i++){if(i>1)printf(",");print_poly(gel(r,i));}printf("]");}
            printf("\n");
        } pari_ENDCATCH;
        set_avma(av);fflush(stdout);
    }
    pari_close();return 0;
}
