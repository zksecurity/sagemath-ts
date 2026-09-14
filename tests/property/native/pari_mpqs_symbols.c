#include "pari.h"
#include <stdio.h>
int main(){char line[4000];pari_init(64000000,500000);while(fgets(line,sizeof(line),stdin)){pari_sp av=avma;GEN a=gp_read_str(line);long op=itos(gel(a,1));GEN x=gel(a,2);ulong y=itou(gel(a,3));printf("OK %ld\n",op?kroiu(x,y):krouu(itou(x),y));set_avma(av);fflush(stdout);}pari_close();}
