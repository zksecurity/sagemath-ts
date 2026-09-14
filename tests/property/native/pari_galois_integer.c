#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}


/* Number-returning facade: display the exact integer represented by binary64. */
static void word(ulong x){printf("%.0f",(double)x);}
int main(void){long op;static char sa[30000],sb[30000],sk[30000];pari_init(256000000,500000);
while(scanf("%ld %29999s %29999s %29999s",&op,sa,sb,sk)==4){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
pari_TRY{GEN a=gp_read_str(sa),b=gp_read_str(sb),k=gp_read_str(sk),z=NULL;ulong u,v;long value=0;
if(op<4){u=itou(a);if(op<2){v=itou(b);u=op==0?ugcd(u,v):ulcm(u,v);}else if(op==2)z=factoru(u);}
else value=logint0(a,b,NULL);
printf("OK ");if(op<2)word(u);else if(op==2){printf("[");for(long i=1;i<lg(gel(z,1));i++){if(i>1)printf(",");printf("[");word(uel(gel(z,1),i));printf(",");word(uel(gel(z,2),i));printf("]");}printf("]");}
else if(op==3){forprime_t T;u_forprime_init(&T,u,ULONG_MAX);printf("[");for(long i=0;i<itos(k);i++){if(i)printf(",");word(u_forprime_next(&T));}printf("]");}
else printf("%ld",value);printf("\n");}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
