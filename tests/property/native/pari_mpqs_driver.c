#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
int main(){char line[1000];pari_init(512000000,500000);while(fgets(line,sizeof(line),stdin)){pari_sp av=avma;GEN a=gp_read_str(line),N=gel(a,1);setrand(gel(a,2));GEN r=mpqs(N);printf("OK [");if(!r)printf("null");else if(typ(r)==t_INT)pari_printf("[[%Ps,1],[%Ps,1]]",r,diviiexact(N,r));else {printf("[");for(long i=1;i<lg(r);i+=3){if(i>1)printf(",");pari_printf("[%Ps,%Ps]",gel(r,i),gel(r,i+1));}printf("]");}pari_printf(",%Ps]\n",getrand());set_avma(av);fflush(stdout);}pari_close();}