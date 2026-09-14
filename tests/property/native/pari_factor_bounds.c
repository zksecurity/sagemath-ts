/* Direct scalar dependencies and private QX factorization bound routines. */
#include "QX_factor.c"
#include <stdio.h>
static GEN polynomial(GEN v){long l=lg(v);GEN P=cgetg(l+1,t_POL);P[1]=evalsigne(1)|evalvarn(0);for(long i=1;i<l;i++)gel(P,i+1)=gel(v,i);return normalizepol_lg(P,l+1);}
static void value(GEN x){if(!x){printf("null");return;}if(typ(x)==t_INT){pari_printf("\"%Ps\"",x);return;}if(typ(x)==t_REAL){GEN m=gen_0;if(signe(x))for(long i=2;i<lg(x);i++)m=addii(shifti(m,64),utoi((ulong)x[i]));printf("[%ld,\"%ld\",",signe(x),expo(x));pari_printf("\"%Ps\",%ld]",m,(lg(x)-2)*64);return;}printf("[");for(long i=1;i<lg(x);i++){if(i>1)printf(",");value(gel(x,i));}printf("]");}
static void error(const char*s){printf("ERROR ");putchar(34);for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}putchar(34);putchar(10);}
int main(void){static char data[4000000];long op,p,shift;ulong n;pari_init(256000000,500000);
while(scanf("%ld %ld %ld %lu %3999999s",&op,&p,&shift,&n,data)==5){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*s=pari_err2str(pari_err_last());error(s);pari_free(s);}
pari_TRY{GEN a=gp_read_str(data),x,r;
if(op==0){x=n==0?shiftr(itor(a,p),shift):n==1?a:gdiv(gel(a,1),gel(a,2));r=ceil_safe(x);}
else if(op==1)r=powruhalf(shiftr(itor(a,p),shift),n);
else if(op==6)r=vecbinomial(n);
else if(op==10){x=Qdivii(gel(a,1),gel(a,2));r=mkvec2(numer_i(x),denom_i(x));}
else if(op==9)r=ZX_Z_eval(polynomial(vecslice(a,2,lg(a)-1)),gel(a,1));
else if(op==7||op==8){x=shiftr(itor(gel(a,2),p),shift);r=stoi(op==7?cmpir(gel(a,1),x):cmpri(x,gel(a,1)));}
else {x=polynomial(a);r=op==2?Mignotte_bound(x):op==3?Beauzamy_bound(x):op==4?factor_bound(x):root_bound(x);}
value(r);printf("\n");}pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
