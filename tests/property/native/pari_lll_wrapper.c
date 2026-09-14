/* Bundled PARI adaptive LLL and FLATTER. Internal functions from lll.c. */
#include "lll.c"
#include <stdio.h>
#ifndef LONG_IS_64BIT
#error "Adaptive LLL oracle requires 64-bit PARI"
#endif
static void value(GEN x) {
  if (!x) { printf("null"); return; }
  if (typ(x)==t_INT) { pari_printf("\"%Ps\"",x); return; }
  if (typ(x)==t_REAL) {
    GEN m=gen_0;
    if(signe(x))for(long i=2;i<lg(x);i++)m=addii(shifti(m,64),utoi((ulong)x[i]));
    printf("[%ld,\"%ld\",",signe(x),expo(x));pari_printf("\"%Ps\",%ld]",m,(lg(x)-2)*64);return;
  }
  printf("[");for(long i=1;i<lg(x);i++){if(i>1)printf(",");value(gel(x,i));}printf("]");
}
static void error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}

int main(void){
 long op,flag,mode,m,n,p,shift;static char data[4000000];pari_init(256000000,500000);
 while(scanf("%ld %ld %ld %ld %ld %ld %ld %3999999s",&op,&flag,&mode,&m,&n,&p,&shift,data)==8){
  pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *s=pari_err2str(pari_err_last());error(s);pari_free(s);}
  pari_TRY{
   GEN v=gp_read_str(data),A=cgetg(n+1,t_MAT),r;
   for(long j=1;j<=n;j++){
    gel(A,j)=cgetg(m+1,t_COL);
    for(long i=1;i<=m;i++){
     GEN x=gel(v,(i-1)*n+j);
     gcoeff(A,i,j)=mode==1||(mode==2&&(i+j)%2)?shiftr(itor(x,p),shift):x;
    }
   }
   if(op==0)r=ZM_lll(A,.99,flag);
   else if(op==1)r=lllfp(A,.99,flag);
   else if(op==2){GEN T=NULL;long s,pot;GEN M=flat(A,flag,&T,&s,&pot);r=mkvec4(M,T?T:gen_0,stoi(s),stoi(pot));}
   else if(op==3)r=ZM_flatter(A,flag);
   else if(op==4)r=ZM_flatter_rank(A,ZM_rank(A),flag);
   else if(op==5)r=ZM_flattergram(A,flag);
   else {GEN a=gel(v,1),b=gel(v,2),c=gel(v,3),U,Q=redimagsl2(mkqfb(a,b,c,qfb_disc3(a,b,c)),&U);r=mkvec2(mkvec3(gel(Q,1),gel(Q,2),gel(Q,3)),U);}
   value(r);printf("\n");
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }
 pari_close();return 0;
}
