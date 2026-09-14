#include "../../../reference/pari/src/basemath/ifactor1.c"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}




static void print_gen(GEN v){long t=typ(v);if(t==t_POL){printf("[");for(long i=2;i<lg(v);i++){if(i>2)printf(",");print_gen(gel(v,i));}printf("]");return;}if(t==t_INT){pari_printf("%Ps",v);return;}printf("[");for(long i=1;i<lg(v);i++){if(i>1)printf(",");if(t==t_VECSMALL)printf("%ld",v[i]);else print_gen(gel(v,i));}printf("]");}
static GEN matrix(GEN flat,long n){long k=(lg(flat)-1)/n;GEN M=cgetg(k+1,t_MAT);for(long i=1;i<=k;i++){GEN c=cgetg(n+1,t_COL);for(long j=1;j<=n;j++)gel(c,j)=gel(flat,(i-1)*n+j);gel(M,i)=c;}return M;}


static GEN
bounded_ellfacteur(GEN N, int insist, long rounds)
{
  const long size = expi(N) + 1;
  pari_sp av = avma;
  struct ECM E;
  long nbc, dsn, dsnmax, rep = 0;
  if (insist)
  {
    const long DSNMAX = numberof(TB1)-1;
    dsnmax = (size >> 2) - 10;
    if (dsnmax < 0) dsnmax = 0;
    else if (dsnmax > DSNMAX) dsnmax = DSNMAX;
    E.seed = 1 + (nbcmax<<7)*(size&0xffff); /* seed for choice of curves */

    dsn = (size >> 3) - 5;
    if (dsn < 0) dsn = 0; else if (dsn > 47) dsn = 47;
    /* pick up the torch where noninsistent stage would have given up */
    nbc = dsn + (dsn >> 2) + 9; /* 8 or more curves in parallel */
    nbc &= ~3; /* 4 | nbc */
  }
  else
  {
    dsn = (size - 140) >> 3;
    if (dsn < 0)
    {
#ifndef __EMX__ /* unless DOS/EMX: MPQS's disk access is abysmally slow */
      if (DEBUGLEVEL >= 4)
        err_printf("ECM: number too small to justify this stage\n");
      return NULL; /* too small, decline the task */
#endif
      dsn = 0;
    } else if (dsn > 12) dsn = 12;
    rep = (size <= 248 ?
           (size <= 176 ? (size - 124) >> 4 : (size - 148) >> 3) :
           (size - 224) >> 1);
#ifdef __EMX__ /* DOS/EMX: extra rounds (shun MPQS) */
    rep += 20;
#endif
    dsnmax = 72;
    /* Use disjoint sets of curves for non-insist and insist phases; moreover,
     * repeated calls acting on factors of the same original number should try
     * to use fresh curves. The following achieves this */
    E.seed = 1 + (nbcmax<<3)*(size & 0xf);
    nbc = -1;
  }
  ECM_init(&E, N, nbc);
  if (DEBUGLEVEL >= 4)
  {
    timer_start(&E.T);
    err_printf("ECM: working on %ld curves at a time; initializing", E.nbc);
    if (!insist)
    {
      if (rep == 1) err_printf(" for one round");
      else          err_printf(" for up to %ld rounds", rep);
    }
    err_printf("...\n");
  }
  if (dsn > dsnmax) dsn = dsnmax;
  for(long round=0; round<rounds; round++)
  {
    ulong B1 = insist? TB1[dsn]: TB1_for_stage[dsn];
    GEN g = ECM_loop(&E, N, B1);
    if (g)
    {
      if (DEBUGLEVEL >= 4)
        err_printf("ECM: time = %6ld ms\n\tfound factor = %Ps\n",
                   timer_delay(&E.T), g);
      return gc_GEN(av, g);
    }
    if (dsn < dsnmax)
    {
      if (insist) dsn++;
      else { dsn += 2; if (dsn > dsnmax) dsn = dsnmax; }
    }
    if (!insist && !--rep)
    {
      if (DEBUGLEVEL >= 4)
        err_printf("ECM: time = %6ld ms,\tellfacteur giving up.\n",
                   timer_delay(&E.T));
      return gc_NULL(av);
    }
  }
  return gc_NULL(av);
}

int main(void){long op,nbc,seed,rounds;ulong B1;static char sn[30000];pari_init(256000000,500000);
while(scanf("%ld %29999s %ld %ld %lu %ld",&op,sn,&nbc,&seed,&B1,&rounds)==6){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{GEN N=gp_read_str(sn),g=NULL;if(op==2)g=ellfacteur(N,0);else if(op==3)g=bounded_ellfacteur(N,1,rounds);else{struct ECM E;E.seed=seed;ECM_init(&E,N,nbc);for(long i=0;i<rounds;i++){g=ECM_loop(&E,N,B1);if(g)break;}}
printf("OK ");print_gen(g?g:gen_0);printf("\n");}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
