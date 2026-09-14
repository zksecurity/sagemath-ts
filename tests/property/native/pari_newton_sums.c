/* Native integer/rational polynomial Newton sums, including supplied prefixes. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static GEN polynomial(GEN v){long l=lg(v);GEN P=cgetg(l+1,t_POL);P[1]=evalsigne(1)|evalvarn(0);for(long i=1;i<l;i++)gel(P,i+1)=gel(v,i);return normalizepol_lg(P,l+1);}
static void value(GEN x,int rational){if(!x){printf("null");return;}if(typ(x)==t_INT||typ(x)==t_FRAC){if(rational){printf("[");pari_printf("\"%Ps\",\"%Ps\"",numer_i(x),denom_i(x));printf("]");}else pari_printf("\"%Ps\"",x);return;}printf("[");for(long i=1;i<lg(x);i++){if(i>1)printf(",");value(gel(x,i),rational);}printf("]");}
static void error(const char*s){printf("ERROR ");putchar(34);for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}putchar(34);putchar(10);}
int main(void){static char ns[4000000],ps[4000000],cs[4000000];long op,n,usecache;pari_init(256000000,500000);
while(scanf("%ld %ld %3999999s %ld %3999999s %3999999s",&op,&n,ns,&usecache,ps,cs)==6){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*s=pari_err2str(pari_err_last());error(s);pari_free(s);}
pari_TRY{GEN N=gp_read_str(ns),P=polynomial(gp_read_str(ps)),cache=gp_read_str(cs),y=NULL,r;
if(usecache){long len=lg(cache)-1;if(!signe(N)){y=cgetg(len/2+1,t_COL);for(long i=1;i<lg(y);i++)gel(y,i)=Qdivii(gel(cache,2*i-1),gel(cache,2*i));}else y=gtocol(cache);}
r=op?polsym(P,n):polsym_gen(P,y,n,NULL,signe(N)?N:NULL);value(r,op||!signe(N));printf("\n");}pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
