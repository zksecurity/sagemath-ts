/* Direct bundled PARI ABI driver for packed binary arithmetic and kernel bases. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static GEN unpack(GEN x) {
    long n=lgefint(x)-2;
    GEN f=cgetg(n+2,t_VECSMALL);f[1]=0;
    for(long i=0;i<n;i++)f[i+2]=*int_W(x,i);
    return f;
}
static GEN packed(GEN x) {
    GEN z=gen_0;
    for(long i=lg(x)-1;i>=2;i--)z=addii(shifti(z,64),utoi((ulong)x[i]));
    return z;
}
static GEN pack_columns(GEN columns) {
    GEN v=cgetg(lg(columns),t_VEC);
    for(long i=1;i<lg(columns);i++)gel(v,i)=packed(gel(columns,i));
    return v;
}
static void print_vec(GEN v) {
    for(long i=1;i<lg(v);i++){if(i>1)printf(",");pari_printf("%Ps",gel(v,i));}
}
int main(void) {
    long op,d,e;char sa[16000],sb[16000],sc[16000];
    pari_init(128000000,500000);
    while(scanf("%ld %15999s %15999s %15999s %ld %ld",&op,sa,sb,sc,&d,&e)==6) {
        pari_sp av=avma;
        pari_CATCH(CATCH_ALL) {
            char *error=pari_err2str(pari_err_last());printf("ERROR %s\n",error);pari_free(error);
        } pari_TRY {
            GEN ai=strtoi(sa),bi=strtoi(sb),ci=strtoi(sc);
            GEN a=unpack(ai),b=unpack(bi),c=unpack(ci),result=NULL,matrix=NULL;
            long alias=0,kind=0;
            switch(op) {
                case 0:result=mkvec(packed(F2x_mul(a,b)));break;
                case 1:result=mkvec(packed(F2x_sqr(a)));break;
                case 2:result=mkvec(packed(F2x_sqrt(a)));break;
                case 3:result=mkvec(packed(F2x_gcd(a,b)));break;
                case 4:{GEN r,q=F2x_divrem(a,b,&r);result=mkvec2(packed(q),packed(r));break;}
                case 5:result=mkvec(packed(F2x_rem(a,b)));break;
                case 6:result=mkvec(packed(F2x_deriv(a)));break;
                case 7:result=mkvec(packed(F2x_Frobenius(a)));break;
                case 8:result=pack_columns(F2x_matFrobenius(a));break;
                case 9:{
                    long cols=itos(bi);matrix=zero_F2m_copy(d,cols);
                    for(long i=1;i<=cols;i++)for(long j=1;j<=d;j++)if(bittest(ai,(i-1)*d+j-1))F2m_set(matrix,j,i);
                    GEN r=e==2?F2m_ker(matrix):F2m_ker_sp(matrix,e);
                    alias=r==matrix;
                    if(!r)kind=2;
                    else if(typ(r)==t_VECSMALL){kind=1;result=mkvec(packed(r));}
                    else result=pack_columns(r);
                    break;
                }
                case 10:result=pack_columns(F2xq_powers(a,d,b));break;
                case 11:result=mkvec(stoi(F2x_degree(a)));break;
                case 12:{GEN r;long v=F2x_valrem(a,&r);result=mkvec2(stoi(v),packed(r));break;}
                case 13:result=mkvec(packed(F2xq_mul(a,b,c)));break;
                case 14:result=mkvec(packed(F2xq_sqr(a,b)));break;
                case 15:result=mkvec(packed(F2x_add(a,b)));break;
                case 16:case 17:case 18:{
                    setrand(bi);
                    if(op==17)result=pack_columns(F2x_factor_squarefree(a));
                    else {
                        GEN fac=op==16?F2x_factor(a):F2x_ddf(a);
                        GEN polys=gel(fac,1),exps=gel(fac,2);
                        result=cgetg(1+2*(lg(polys)-1),t_VEC);
                        for(long i=1;i<lg(polys);i++){gel(result,2*i-1)=packed(gel(polys,i));gel(result,2*i)=stoi(exps[i]);}
                    }
                    break;
                }
            }
            printf("OK ");
            if(op==9){
                printf(kind==2?"null":kind==1?"v:":"m:");
                if(result)print_vec(result);printf("|");print_vec(pack_columns(matrix));printf("|%ld",alias);
            }else print_vec(result);
            if(op>=16)pari_printf("|%Ps",getrand());
            printf("\n");
        } pari_ENDCATCH;
        set_avma(av);fflush(stdout);
    }
    pari_close();return 0;
}
