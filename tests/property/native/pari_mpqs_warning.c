#include "pari.h"
#include <stdio.h>
#include <stdarg.h>
static long warnings=0;
static void audit_warn(int code,const char*fmt,...){char s[1000];va_list args;va_start(args,fmt);vsnprintf(s,sizeof(s),fmt,args);va_end(args);if(strcmp(s,"MPQS: wrong relation found after Gauss")==0)warnings++;pari_warn(code,"%s",s);}
#define pari_warn audit_warn
#define MPQS_DEBUG
#include "../../../reference/pari/src/basemath/mpqs.c"
#undef pari_warn
static void factor_result(GEN r){if(!r){printf("null");return;}printf("[");for(long i=1;i<lg(r);i+=3){if(i>1)printf(",");pari_printf("[%Ps,%Ps]",gel(r,i),gel(r,i+1));}printf("]");}
int main(){char*line=NULL;size_t cap=0;pari_init(128000000,500000);while(getline(&line,&cap,stdin)>0){pari_sp av=avma;GEN a=gp_read_str(line),packed=gel(a,2);mpqs_handle_t h={0};mpqs_FB_entry_t fb[9]={0};h.N=gel(a,1);h.FB=fb;h.size_of_FB=6;long ps[]={0,0,2,3,5,7,11,13};for(long i=2;i<8;i++)fb[i].fbe_p=ps[i];hashtable frel;hash_init_GEN(&frel,20,gequal,1);for(long i=1;i<lg(packed);){GEN Y=gel(packed,i++);long n=itos(gel(packed,i++));GEN r=cgetg(n+1,t_VECSMALL);for(long j=1;j<=n;j++)r[j]=itos(gel(packed,i++));frel_add(&frel,mkvec2(Y,r));}warnings=0;GEN r=mpqs_solve_linear_system(&h,&frel);printf("OK [");factor_result(r);printf(",%ld]\n",warnings);set_avma(av);fflush(stdout);}pari_close();}
