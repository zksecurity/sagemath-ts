/* Direct bundled PARI matrix multiplication and dispatch-boundary oracle. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static GEN entry(long i,long j,long bits,long seed) {
    if(!bits||!seed||(seed%4==3&&(i+3*j)%4))return gen_0;
    GEN value=addii(int2n(bits-1),stoi((i*17+j*23+seed*31)%257));
    value=modii(value,int2n(bits));
    return seed%2==0&&(i+j)%2?negi(value):value;
}
static GEN matrix(long rows,long cols,long bits,long seed) {
    GEN a=cgetg(cols+1,t_MAT);
    for(long j=1;j<=cols;j++){
        gel(a,j)=cgetg(rows+1,t_COL);
        for(long i=1;i<=rows;i++)gcoeff(a,i,j)=entry(i,j,bits,seed);
    }
    return a;
}
static void print_matrix(GEN a) {
    printf("[");
    for(long j=1;j<lg(a);j++){
        if(j>1)printf(",");printf("[");
        for(long i=1;i<lg(gel(a,j));i++){if(i>1)printf(",");pari_printf("%Ps",gcoeff(a,i,j));}
        printf("]");
    }
    printf("]");
}
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}

int main(void) {
    long op,m,n,k,ab,bb,seed;char sp[16000];
    pari_init(512000000,500000);
    while(scanf("%ld %15999s %ld %ld %ld %ld %ld %ld",&op,sp,&m,&n,&k,&ab,&bb,&seed)==8){
        pari_sp av=avma;
        pari_CATCH(CATCH_ALL){char *error=pari_err2str(pari_err_last());print_error(error);pari_free(error);}
        pari_TRY{
            GEN p=gp_read_str(sp),a=matrix(m,n,ab,seed),b=matrix(n,k,bb,seed?seed+5:0),r;
            if(op==0)r=ZM_mul(a,b);
            else if(op==1){ulong q=itou(p);r=Flm_to_ZM(Flm_mul(ZM_to_Flm(a,q),ZM_to_Flm(b,q),q));}
            else if(op==2)r=FpM_mul(a,b,p);
            else r=F2m_to_ZM(F2m_mul(ZM_to_F2m(a),ZM_to_F2m(b)));
            printf("OK ");print_matrix(r);printf("\n");
        }pari_ENDCATCH;
        set_avma(av);fflush(stdout);
    }
    pari_close();return 0;
}
