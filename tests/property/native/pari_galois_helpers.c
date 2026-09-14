#include "../../../reference/pari/src/basemath/galconj.c"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}




static void print_gen(GEN v){long t=typ(v);if(t==t_POL){printf("[");for(long i=2;i<lg(v);i++){if(i>2)printf(",");print_gen(gel(v,i));}printf("]");return;}if(t==t_INT){pari_printf("%Ps",v);return;}printf("[");for(long i=1;i<lg(v);i++){if(i>1)printf(",");if(t==t_VECSMALL)printf("%ld",v[i]);else print_gen(gel(v,i));}printf("]");}
static GEN matrix(GEN flat,long n){long k=(lg(flat)-1)/n;GEN M=cgetg(k+1,t_MAT);for(long i=1;i<=k;i++){GEN c=cgetg(n+1,t_COL);for(long j=1;j<=n;j++)gel(c,j)=gel(flat,(i-1)*n+j);gel(M,i)=c;}return M;}

int main(void){long op;static char sn[30000],sa[200000],sb[200000],sf[30000];pari_init(256000000,500000);
while(scanf("%ld %29999s %199999s %199999s %29999s",&op,sn,sa,sb,sf)==5){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{GEN N=gp_read_str(sn),A=gp_read_str(sa),B=gp_read_str(sb),z=NULL;ulong result=0;
if(op<=1){GEN a=gtopolyrev(A,0),b=gtopolyrev(B,0);if(op==0)z=ZX_sub(a,b);else result=gequal(a,b);}
else if(op==2)z=pol_xn(itos(N),0);
else if(op==3){long n=itos(N);GEN M=matrix(A,n);for(long i=1;i<lg(M);i++)gel(M,i)=ZV_to_zv(gel(M,i));z=galoisfindgroups(M,ZV_to_zv(B),itos(gp_read_str(sf)));}
else{ulong n=itou(N);result=op==4?radicalu(n):eulerphiu(n);}
printf("OK ");if(z)print_gen(z);else if(op==1||op==6)printf("%lu",result);else printf("%.0f",(double)result);printf("\n");}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
