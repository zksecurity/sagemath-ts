/* Original private/public integer factorization stages and random state. */
#include "QX_factor.c"
#include <stdio.h>
static GEN polynomial(GEN v){long l=lg(v);GEN P=cgetg(l+1,t_POL);P[1]=evalsigne(1)|evalvarn(0);for(long i=1;i<l;i++)gel(P,i+1)=gel(v,i);return normalizepol_lg(P,l+1);}
static GEN nested(GEN v,int matrix){long count=0,j=1;while(j<lg(v)){j+=itos(gel(v,j))+1;count++;}GEN f=cgetg(count+1,matrix?t_MAT:t_VEC);for(long i=1,j=1;i<=count;i++){long n=itos(gel(v,j++));GEN c=vecslice(v,j,j+n-1);gel(f,i)=matrix?gtocol(c):polynomial(c);j+=n;}return f;}
static void value(GEN x){if(!x){printf("null");return;}if(typ(x)==t_INT){pari_printf("\"%Ps\"",x);return;}printf("[");long start=typ(x)==t_POL?2:1;for(long i=start;i<lg(x);i++){if(i>start)printf(",");if(typ(x)==t_VECSMALL)printf("\"%ld\"",x[i]);else value(gel(x,i));}printf("]");}
static void error(const char*s){printf("ERROR ");putchar(34);for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}putchar(34);putchar(10);}
int main(void){static char ps[4000000],bs[4000000],xs[4000000],fs[4000000],ms[4000000];long op,a,b,k;pari_init(256000000,500000);
while(scanf("%ld %3999999s %ld %ld %ld %3999999s %3999999s %3999999s %3999999s",&op,ps,&a,&b,&k,bs,xs,fs,ms)==9){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*s=pari_err2str(pari_err_last());error(s);pari_free(s);}
pari_TRY{GEN p=gp_read_str(ps),B=gp_read_str(bs),x=polynomial(gp_read_str(xs)),f=gp_read_str(fs),M=gp_read_str(ms),r=NULL;setrand(gen_1);
if(op==0){GEN R,H=ZX_gcd_all(x,polynomial(f),&R);r=mkvec2(H,R);}
else if(op==1){GEN E,F=ZX_squff(x,&E);r=mkvec2(F,E);}
else if(op==2){pari_timer timer;r=utoi(pick_prime(x,k,&timer));}
else if(op==3)r=ZX_DDF_max(x,k);
else if(op==4)r=ZX_DDF(x);
else if(op==5)r=ZX_factor(x);
else if(op==6)r=stoi(ZX_is_irred(x));
else if(op==7)r=combine_factors(x,nested(f,0),p,k);
else if(op==8)r=LLL_cmbf(x,nested(f,0),p,powiu(p,a),B,a,b);
else if(op==9)r=chk_factors(x,nested(M,1),B,nested(f,0),powiu(p,a));
else if(op==10)r=chk_factors_get(signe(B)?B:NULL,nested(f,0),gtocol(M),NULL,p);
else if(op==16)r=stoi(numberofconjugates(x,a));
else if(op==17)r=stoi(Flx_is_squarefree(ZX_to_Flx(x,itou(p)),itou(p)));
else if(op==18)r=Flx_to_ZX(Flx_deriv(ZX_to_Flx(x,itou(p)),itou(p)));
/* Sage factor/cached-irreducibility PARI effects. NTL calls do not alter PARI's
 * stream; a=0 represents that branch, whose factors are checked in Sage/NTL. */
else if(op==19){r=cgetg(k+1,t_VEC);for(long i=1;i<=k;i++){if(b)QX_factor(gdiv(x,B));else{GEN c=ZX_content(x);Z_factor(c);if(a)ZX_factor(ZX_Z_divexact(x,c));}gel(r,i)=getrand();}}
else if(op==15)r=galoisconj(x,NULL);
else if(op==14)r=stoi(polisirreducible(x));
else if(op==12)r=QX_factor(gdiv(x,B));
else if(op==13){GEN A=factor_bound(x),bound2=mului(degpol(x),sqri(mulii(absi(leading_coeff(x)),root_bound(x)))),pa,pb;long aa,bb,maxK;int done;cmbf_precs(p,A,bound2,&aa,&bb,&pa,&pb);GEN F=gel(FpX_factor(FpX_normalize(FpX_red(x,p),p),p),1);F=ZpX_liftfact(x,F,pa,p,aa);GEN L=cmbf(x,F,A,p,aa,bb,degpol(x)-1,&maxK,&done),R=gel(L,1),G=gel(L,2);long last=lg(R)-1;GEN ff=gel(G,last);if(maxK>0&&lg(ff)-1>2*maxK){GEN P=gel(R,last);if(last!=1)A=factor_bound(P);r=mkvec5(P,ff,stoi(aa),stoi(maxK),A);}}

printf("[");value(r);printf(",");value(getrand());printf("]\n");}pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
