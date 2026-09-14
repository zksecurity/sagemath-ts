/* Original private bounded recombination, with fixed lifted-factor inputs. */
#include "QX_factor.c"
#include <stdio.h>
static GEN polynomial(GEN v){long l=lg(v);GEN P=cgetg(l+1,t_POL);P[1]=evalsigne(1)|evalvarn(0);for(long i=1;i<l;i++)gel(P,i+1)=gel(v,i);return normalizepol_lg(P,l+1);}
static GEN factors(GEN v){long count=0,j=1;while(j<lg(v)){j+=itos(gel(v,j))+1;count++;}GEN f=cgetg(count+1,t_VEC);for(long i=1,j=1;i<=count;i++){long n=itos(gel(v,j++));gel(f,i)=polynomial(vecslice(v,j,j+n-1));j+=n;}return f;}
static void value(GEN x){if(!x){printf("null");return;}if(typ(x)==t_INT){pari_printf("\"%Ps\"",x);return;}printf("[");long start=typ(x)==t_POL?2:1;for(long i=start;i<lg(x);i++){if(i>start)printf(",");if(typ(x)==t_VECSMALL)printf("\"%ld\"",x[i]);else value(gel(x,i));}printf("]");}
static void error(const char*s){printf("ERROR ");putchar(34);for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}putchar(34);putchar(10);}
int main(void){static char ps[4000000],bs[4000000],xs[4000000],fs[4000000];long op,a,b,k;pari_init(256000000,500000);
while(scanf("%ld %3999999s %ld %ld %ld %3999999s %3999999s %3999999s",&op,ps,&a,&b,&k,bs,xs,fs)==8){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*s=pari_err2str(pari_err_last());error(s);pari_free(s);}
pari_TRY{GEN p=gp_read_str(ps),B=gp_read_str(bs),v=gp_read_str(xs),f=gp_read_str(fs),x=polynomial(v),r=NULL;
if(op==0||op==6)r=ZX_divides_i(x,polynomial(f),op==6||signe(B)<0?NULL:B);
else if(op==1){GEN pa,pb;long aa,bb;int flag=cmbf_precs(p,gel(v,1),B,&aa,&bb,&pa,&pb);r=mkvec5(stoi(flag),stoi(aa),stoi(bb),pa,pb);}
else if(op==2){long maxK;int done;GEN L=cmbf(x,factors(f),B,p,a,b,k,&maxK,&done);r=mkvec4(gel(L,1),gel(L,2),stoi(maxK),stoi(done));}
else if(op==3)r=centermodii(gel(v,1),p,a?B:NULL);
else if(op==4)r=FpXV_prod(factors(f),p);
else if(op==5)r=stoi(cmbf_maxK(a));
else if(op==7){if(!umodiu(leading_coeff(x),itou(p))||!FpX_is_squarefree(FpX_red(x,p),p))r=NULL;else{GEN A=factor_bound(x),bb=mului(degpol(x),sqri(mulii(absi(leading_coeff(x)),root_bound(x)))),pa,pb;long aa,ab;cmbf_precs(p,A,bb,&aa,&ab,&pa,&pb);setrand(gen_1);GEN ff=gel(FpX_factor(FpX_normalize(FpX_red(x,p),p),p),1);ff=ZpX_liftfact(x,ff,pa,p,aa);r=mkvec4(stoi(aa),stoi(ab),A,ff);}}
value(r);printf("\n");}pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
