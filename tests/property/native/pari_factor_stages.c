#include "../../../reference/pari/src/basemath/ifactor1.c"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}




static void print_gen(GEN v){long t=typ(v);if(t==t_POL){printf("[");for(long i=2;i<lg(v);i++){if(i>2)printf(",");print_gen(gel(v,i));}printf("]");return;}if(t==t_INT){pari_printf("%Ps",v);return;}printf("[");for(long i=1;i<lg(v);i++){if(i>1)printf(",");if(t==t_VECSMALL)printf("%ld",v[i]);else print_gen(gel(v,i));}printf("]");}
static GEN matrix(GEN flat,long n){long k=(lg(flat)-1)/n;GEN M=cgetg(k+1,t_MAT);for(long i=1;i<=k;i++){GEN c=cgetg(n+1,t_COL);for(long j=1;j<=n;j++)gel(c,j)=gel(flat,(i-1)*n+j);gel(M,i)=c;}return M;}

static void factors(GEN z,GEN n){if(!z){printf("0");return;}if(typ(z)==t_INT){print_gen(mkvec2(z,diviiexact(n,z)));return;}printf("[");long first=1;GEN product=gen_1;for(long i=1;i<lg(z);i+=3)for(long e=0;e<itos(gel(z,i+1));e++){if(!first)printf(",");first=0;print_gen(gel(z,i));product=mulii(product,gel(z,i));}GEN cofactor=diviiexact(n,product);if(!equali1(cofactor)){if(!first)printf(",");print_gen(cofactor);}printf("]");}
int main(void){long op,rounds,seed;static char sn[30000];pari_init(256000000,500000);
while(scanf("%ld %29999s %ld %ld",&op,sn,&rounds,&seed)==4){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{GEN n=gp_read_str(sn),z=op==0?squfof(n):op==1?pollardbrent(n):Z_pollardbrent(n,rounds,seed);printf("OK ");if(op==2){if(z)print_gen(z);else printf("0");}else factors(z,n);printf("\n");}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
