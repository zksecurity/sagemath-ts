/* ABI driver: all generator and sampling algorithms come from bundled PARI. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void print_poly(GEN x,int packed,int small) {
    if(packed) {
        GEN z=gen_0;
        for(long i=F2x_degree(x);i>=0;i--)z=addii(shifti(z,1),utoi(F2x_coeff(x,i)));
        pari_printf("%Ps",z);return;
    }
    long start=small==2?1:2;
    for(long i=start;i<lg(x);i++) {
        if(i>start)printf(",");
        if(small==2)printf("%ld",x[i]);
        else if(small)printf("%lu",(ulong)x[i]);
        else pari_printf("%Ps",gel(x,i));
    }
}
int main(void) {
    long op,count,length;char seed[4096],limit[4096];
    pari_init(64000000,500000);
    while(scanf("%ld %4095s %4095s %ld %ld",&op,seed,limit,&count,&length)==5) {
        pari_sp av=avma;
        pari_CATCH(CATCH_ALL) {
            char *error=pari_err2str(pari_err_last());printf("ERROR %s\n",error);pari_free(error);
        } pari_TRY {
            pari_init_rand();
            if(strcmp(seed,"-1"))setrand(strtoi(seed));
            GEN n=strtoi(limit);
            if(op==9)pari_init_rand();
            printf("OK ");
            for(long i=0;i<count;i++) {
                if(i)printf(";");
                switch(op) {
                    case 0:case 9:printf("%lu",pari_rand());break;
                    case 1:printf("%ld",random_bits(itos(n)));break;
                    case 2:printf("%lu",random_Fl(itou(n)));break;
                    case 3:pari_printf("%Ps",randomi(n));break;
                    case 4:print_poly(random_F2x(itos(n),0),1,0);break;
                    case 5:print_poly(random_zv(itos(n)),0,2);break;
                    case 6:print_poly(random_Flx(length,0,itou(n)),0,1);break;
                    case 7:print_poly(random_FpX(length,0,n),0,0);break;
                    case 8:{GEN old=getrand();ulong a=pari_rand();setrand(old);printf("%lu,%lu",a,pari_rand());break;}
                }
            }
            pari_printf("|%Ps\n",getrand());
        } pari_ENDCATCH;
        set_avma(av);fflush(stdout);
    }
    pari_close();return 0;
}
