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



int main(){char*line=NULL;size_t cap=0;pari_init(256000000,500000);while(getline(&line,&cap,stdin)>0){pari_sp av=avma;GEN a=gp_read_str(line);mpqs_handle_t h={0};long ok=mpqs_class_init(&h,gel(a,1),itos(gel(a,2))),rounds=itos(gel(a,3));GEN v=gel(a,4),missing=NULL;if(lg(v)>1){missing=cgetg(lg(v),t_VECSMALL);for(long i=1;i<lg(v);i++)missing[i]=itos(gel(v,i));}printf("OK [%ld,[",ok);hashtable frel,lprel;hash_init_GEN(&frel,h.target_rels,gequal,1);hash_init_ulong(&lprel,h.target_rels,1);if(ok){h.index_j=(mpqs_uint32_t)-1;h.candidates=(long*)stack_malloc(2016*sizeof(long));for(long i=0;i<rounds;i++){if(i)printf(",");ok=mpqs_self_init(&h,missing);if(!ok){printf("[0]");break;}mpqs_sieve(&h);long tc=mpqs_eval_sieve(&h);GEN factor=mpqs_eval_cand(&h,tc,&frel,&lprel,MPQS_MODE_CLASSGROUP);printf("[1,%ld,%lu,",tc,frel.nb);if(factor)print_gen(factor);else printf("0");printf("]");}}printf("],");print_gen(hash_keys_GEN(&frel));printf("]\n");set_avma(av);fflush(stdout);}pari_close();}
