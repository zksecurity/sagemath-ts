/* Independent pinned-M4RI cutoff oracle. No SageMath process survives a native abort. */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <m4ri/m4ri.h>
#include <m4ri/ple_russian.h>
static mzd_t *fixture(int m,int n,uint64_t seed) {
  mzd_t *a=mzd_init(m,n);
  if(n==0)return a;
  for(int i=0;i<m;i++) {
    word *row=mzd_row(a,i);
    for(int j=0;j<a->width;j++) {seed=seed*UINT64_C(6364136223846793005)+UINT64_C(1442695040888963407);row[j]=seed;}
    if(n%64)row[a->width-1]&=(UINT64_C(1)<<(n%64))-1;
  }
  return a;
}
int main(int argc,char **argv) {
  if(argc==7 && strcmp(argv[1],"basic")==0) {
    int m=atoi(argv[2]),n=atoi(argv[3]),right=atoi(argv[4]),op=atoi(argv[6]);
    uint64_t seed=strtoull(argv[5],NULL,10);
    mzd_t *a=fixture(m,n,seed),*b=fixture(m,op==1?n:right,seed+1),*c;
    c=op==0?mzd_transpose(NULL,a):op==1?mzd_add(NULL,a,b):mzd_concat(NULL,a,b);
    printf("%d %d\n",c->nrows,c->ncols);
    for(int i=0;i<c->nrows;i++) {
      if(c->ncols==0)puts("0");
      else {const word *row=mzd_row(c,i);for(int j=c->width-1;j>=0;j--)printf("%016" PRIx64,(uint64_t)row[j]);putchar('\n');}
    }
    mzd_free(a);mzd_free(b);mzd_free(c);return 0;
  }
  if(argc==9 && strcmp(argv[1],"solve")==0) {
    int m=atoi(argv[2]),n=atoi(argv[3]),bc=atoi(argv[4]),pattern=atoi(argv[6]),cutoff=atoi(argv[7]),check=atoi(argv[8]);
    uint64_t seed=strtoull(argv[5],NULL,10);
    mzd_t *a=fixture(m,n,seed),*b=fixture(m>n?m:n,bc,seed+1);
    if(pattern==1) mzd_set_ui(a,0);
    if(pattern==2) mzd_set_ui(a,1);
    if(pattern==3) {
      mzd_t *x=fixture(n,bc,seed+2),*product=mzd_mul(NULL,a,x,0);
      mzd_set_ui(b,0);mzd_copy(b,product);mzd_free(x);mzd_free(product);
    }
    int ret=mzd_solve_left(a,b,cutoff,check);printf("%d\n",ret);
    for(int pass=0;pass<2;pass++) {
      mzd_t *z=pass?b:a;
      for(int i=0;i<z->nrows;i++) {
        if(z->ncols==0)puts("0");
        else {const word *row=mzd_row(z,i);for(int j=z->width-1;j>=0;j--)printf("%016" PRIx64,(uint64_t)row[j]);putchar('\n');}
      }
    }
    mzd_free(a);mzd_free(b);return 0;
  }
  if(argc==7 && strcmp(argv[1],"kernel")==0) {
    int m=atoi(argv[2]),n=atoi(argv[3]),pattern=atoi(argv[5]),cutoff=atoi(argv[6]);
    mzd_t *a=fixture(m,n,strtoull(argv[4],NULL,10));
    if(pattern==1)mzd_set_ui(a,0);
    if(pattern==2)mzd_set_ui(a,1);
    mzd_t *b=mzd_kernel_left_pluq(a,cutoff);printf("%d\n",b?b->ncols:-1);
    for(int pass=0;pass<2;pass++) {
      mzd_t *z=pass?b:a;if(!z)continue;
      for(int i=0;i<z->nrows;i++) {
        if(z->ncols==0)puts("0");
        else {const word *row=mzd_row(z,i);for(int j=z->width-1;j>=0;j--)printf("%016" PRIx64,(uint64_t)row[j]);putchar('\n');}
      }
    }
    mzd_free(a);if(b)mzd_free(b);return 0;
  }
  if(argc==6 && strcmp(argv[1],"inverse")==0) {
    int n=atoi(argv[2]),pattern=atoi(argv[4]),k=atoi(argv[5]);
    mzd_t *a=fixture(n,n,strtoull(argv[3],NULL,10));
    if(pattern)for(int i=0;i<n;i++)for(int j=0;j<n;j++)
      if(pattern==2 || j<=i)mzd_write_bit(a,i,j,pattern==1 && i==j);
    mzd_t *b=mzd_inv_m4ri(NULL,a,k);
    for(int i=0;i<n;i++) {
      const word *row=mzd_row(b,i);for(int j=b->width-1;j>=0;j--)printf("%016" PRIx64,(uint64_t)row[j]);putchar('\n');
    }
    mzd_free(b);mzd_free(a);return 0;
  }
  if(argc==9 && strcmp(argv[1],"echelon_pattern")==0) {
    int m=atoi(argv[2]),n=atoi(argv[3]),pattern=atoi(argv[5]),algorithm=atoi(argv[6]),full=atoi(argv[7]),k=atoi(argv[8]);
    mzd_t *a=fixture(m,n,strtoull(argv[4],NULL,10));
    for(int i=0;i<m;i++)for(int j=0;j<n;j++) {
      int old=mzd_read_bit(a,i,j),v=pattern==0?old:pattern==1?0:pattern==2?j==n-1:pattern==3?(i<400?i==j:old):pattern==4?(j>=n-128?old:i==j):pattern==5?(j<128?0:old):(i>=m/2?0:old);
      mzd_write_bit(a,i,j,v);
    }
    int rank=0;
    if(m&&n)rank=algorithm==0?mzd_echelonize_m4ri(a,full,k):algorithm==1?mzd_echelonize(a,full):mzd_echelonize_pluq(a,full);
    printf("%d\n",rank);
    for(int i=0;i<m;i++) {
      if(n==0)puts("0");
      else {const word *row=mzd_row(a,i);for(int j=a->width-1;j>=0;j--)printf("%016" PRIx64,(uint64_t)row[j]);putchar('\n');}
    }
    mzd_free(a);return 0;
  }
  if(argc==8 && strcmp(argv[1],"echelon")==0) {
    int m=atoi(argv[2]),n=atoi(argv[3]),algorithm=atoi(argv[5]),full=atoi(argv[6]),k=atoi(argv[7]);
    mzd_t *a=fixture(m,n,strtoull(argv[4],NULL,10));
    int rank=0;
    if(m&&n)rank=algorithm==0?mzd_echelonize_m4ri(a,full,k):algorithm==1?mzd_echelonize(a,full):mzd_echelonize_pluq(a,full);
    printf("%d\n",rank);
    for(int i=0;i<m;i++) {
      if(n==0)puts("0");
      else {const word *row=mzd_row(a,i);for(int j=a->width-1;j>=0;j--)printf("%016" PRIx64,(uint64_t)row[j]);putchar('\n');}
    }
    mzd_free(a);return 0;
  }
  if(argc==7 && strcmp(argv[1],"ple")==0) {
    int m=atoi(argv[2]),n=atoi(argv[3]),kind=atoi(argv[5]),k=atoi(argv[6]);
    mzd_t *a=fixture(m,n,strtoull(argv[4],NULL,10));
    mzp_t *p=mzp_init(m),*q=mzp_init(n);int rank=0;
    rank=kind==0?_mzd_ple_russian(a,p,q,k):kind==1?mzd_ple(a,p,q,k):_mzd_ple_naive(a,p,q);
    printf("%d\n",rank);
    for(int i=0;i<m;i++)printf("%d ",p->values[i]);putchar('\n');
    for(int i=0;i<n;i++)printf("%d ",q->values[i]);putchar('\n');
    for(int i=0;i<m;i++) {
      if(n==0)puts("0");
      else {const word *row=mzd_row(a,i);for(int j=a->width-1;j>=0;j--)printf("%016" PRIx64,(uint64_t)row[j]);putchar('\n');}
    }
    mzp_free(p);mzp_free(q);mzd_free(a);return 0;
  }
  if(argc==7 && strcmp(argv[1],"trsm")==0) {
    int m=atoi(argv[2]),n=atoi(argv[3]),upper=atoi(argv[5]),cutoff=atoi(argv[6]);
    uint64_t seed=strtoull(argv[4],NULL,10);
    mzd_t *a=fixture(m,m,seed),*b=fixture(m,n,seed+1);
    for(int i=0;i<m;i++)for(int j=0;j<m;j++)if(i==j|| (upper?j<i:j>i))mzd_write_bit(a,i,j,i==j);
    if(m&&n){if(upper)mzd_trsm_upper_left(a,b,cutoff);else mzd_trsm_lower_left(a,b,cutoff);}
    for(int i=0;i<m;i++) {
      if(n==0)puts("0");
      else {const word *row=mzd_row(b,i);for(int j=b->width-1;j>=0;j--)printf("%016" PRIx64,(uint64_t)row[j]);putchar('\n');}
    }
    mzd_free(a);mzd_free(b);return 0;
  }
  if(argc==8 && strcmp(argv[1],"swap")==0) {
    int m=atoi(argv[2]),n=atoi(argv[3]),axis=atoi(argv[5]),i=atoi(argv[6]),j=atoi(argv[7]);
    mzd_t *a=fixture(m,n,strtoull(argv[4],NULL,10));
    if(axis)mzd_col_swap(a,i,j);else mzd_row_swap(a,i,j);
    mzd_free(a);return 0;
  }
  if(argc==7 && strcmp(argv[1],"factorization")==0) {
    int m=atoi(argv[2]),n=atoi(argv[3]),kind=atoi(argv[5]),k=atoi(argv[6]);
    mzd_t *a=fixture(m,n,strtoull(argv[4],NULL,10));
    mzp_t *p=mzp_init(m),*q=mzp_init(n);int rank;
    rank=kind==0?mzd_ple(a,p,q,k):kind==1?_mzd_ple_russian(a,p,q,k):kind==2?_mzd_ple_naive(a,p,q):kind==3?mzd_pluq(a,p,q,k):kind==4?_mzd_pluq_russian(a,p,q,k):_mzd_pluq_naive(a,p,q);
    printf("%d\n",rank);
    for(int i=0;i<m;i++)printf("%d ",p->values[i]);putchar('\n');
    for(int i=0;i<n;i++)printf("%d ",q->values[i]);putchar('\n');
    for(int i=0;i<m;i++) {
      if(n==0)puts("0");
      else {const word *row=mzd_row(a,i);for(int j=a->width-1;j>=0;j--)printf("%016" PRIx64,(uint64_t)row[j]);putchar('\n');}
    }
    mzp_free(p);mzp_free(q);mzd_free(a);return 0;
  }
  if(argc!=7)return 2;
  if(strcmp(argv[1],"density")==0) {
    int m=atoi(argv[2]),n=atoi(argv[3]),pattern=atoi(argv[4]),resolution=atoi(argv[6]);
    uint64_t seed=strtoull(argv[5],NULL,10);
    mzd_t *a=fixture(m,n,seed);
    if(pattern!=5)for(int i=0;i<m;i++)for(int j=0;j<n;j++) {
      int bit=pattern==0?0:pattern==1?1:pattern==2?j>=n/2:pattern==3?j<n/2:pattern==4?(i*n+j)%2:(uint64_t)i*n+j<seed;
      mzd_write_bit(a,i,j,bit);
    }
    printf("%a\n",mzd_density(a,resolution));mzd_free(a);return 0;
  }
  int m=atoi(argv[1]),k=atoi(argv[2]),n=atoi(argv[3]),cutoff=atoi(argv[5]),square=atoi(argv[6]);
  uint64_t seed=strtoull(argv[4],NULL,10);
  mzd_t *a=fixture(m,k,seed),*b=square?a:fixture(k,n,seed+1),*c;
  c=(m==0||k==0||n==0)?mzd_init(m,n):mzd_mul(NULL,a,b,cutoff);
  for(int i=0;i<m;i++) {
    if(n==0)puts("0");
    else {const word *row=mzd_row(c,i);for(int j=c->width-1;j>=0;j--)printf("%016" PRIx64,(uint64_t)row[j]);putchar('\n');}
  }
  mzd_free(c);if(b!=a)mzd_free(b);mzd_free(a);return 0;
}
