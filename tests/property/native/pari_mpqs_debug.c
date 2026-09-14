#define MPQS_DEBUG
#include "../../../reference/pari/src/basemath/mpqs.c"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}




static void print_gen(GEN v){long t=typ(v);if(t==t_POL){printf("[");for(long i=2;i<lg(v);i++){if(i>2)printf(",");print_gen(gel(v,i));}printf("]");return;}if(t==t_INT){pari_printf("%Ps",v);return;}printf("[");for(long i=1;i<lg(v);i++){if(i>1)printf(",");if(t==t_VECSMALL)printf("%ld",v[i]);else print_gen(gel(v,i));}printf("]");}
static GEN matrix(GEN flat,long n){long k=(lg(flat)-1)/n;GEN M=cgetg(k+1,t_MAT);for(long i=1;i<=k;i++){GEN c=cgetg(n+1,t_COL);for(long j=1;j<=n;j++)gel(c,j)=gel(flat,(i-1)*n+j);gel(M,i)=c;}return M;}



int main(void){char *line=NULL;size_t cap=0;pari_init(128000000,500000);
while(getline(&line,&cap,stdin)>0){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{GEN a=gp_read_str(line);long op=itos(gel(a,1));mpqs_handle_t h={0};mpqs_FB_entry_t fb[9]={0};h.N=gel(a,2);h.size_of_FB=6;h.FB=fb;long ps[]={0,0,2,3,5,7,11,13};for(long i=2;i<8;i++)fb[i].fbe_p=ps[i];
GEN v=gel(a,4),w=cgetg(lg(v),t_VECSMALL);for(long i=1;i<lg(v);i++)w[i]=itos(gel(v,i));
if(!op){GEN z=mpqs_factorback(&h,w);printf("OK ");print_gen(z);putchar(10);}else {mpqs_check_rel(&h,mkvec2(gel(a,3),w),itou(gel(a,5)),itos(gel(a,6)));printf("OK 1\n");}}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
