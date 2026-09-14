#include "../../../reference/pari/src/basemath/mpqs.c"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}




static void print_gen(GEN v){long t=typ(v);if(t==t_POL){printf("[");for(long i=2;i<lg(v);i++){if(i>2)printf(",");print_gen(gel(v,i));}printf("]");return;}if(t==t_INT){pari_printf("%Ps",v);return;}printf("[");for(long i=1;i<lg(v);i++){if(i>1)printf(",");if(t==t_VECSMALL)printf("%ld",v[i]);else print_gen(gel(v,i));}printf("]");}
static GEN matrix(GEN flat,long n){long k=(lg(flat)-1)/n;GEN M=cgetg(k+1,t_MAT);for(long i=1;i<=k;i++){GEN c=cgetg(n+1,t_COL);for(long j=1;j<=n;j++)gel(c,j)=gel(flat,(i-1)*n+j);gel(M,i)=c;}return M;}


int main(void){long q,mode;static char sn[1000],sy1[1000],sy2[1000],sr1[2000],sr2[2000];pari_init(256000000,500000);
while(scanf("%999s %ld %999s %1999s %999s %1999s %ld",sn,&q,sy1,sr1,sy2,sr2,&mode)==7){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{mpqs_handle_t h={0};mpqs_FB_entry_t fb[9]={0};h.N=gp_read_str(sn);h.size_of_FB=6;h.FB=fb;long ps[]={0,0,2,3,5,7,11,13};for(long i=2;i<8;i++)fb[i].fbe_p=ps[i];
GEN a=gp_read_str(sr1),b=gp_read_str(sr2);GEN u=cgetg(lg(a),t_VECSMALL),v=cgetg(lg(b),t_VECSMALL);for(long i=1;i<lg(a);i++)u[i]=itos(gel(a,i));for(long i=1;i<lg(b);i++)v[i]=itos(gel(b,i));
GEN z=combine_large_primes(&h,q,mkvec2(gp_read_str(sy1),u),mkvec2(gp_read_str(sy2),v),mode);printf("OK ");if(z)print_gen(z);else printf("0");printf("\n");}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
