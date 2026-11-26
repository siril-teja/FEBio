/* Standalone MKL PARDISO benchmark on a dumped lower-CSC (== upper CSR) symmetric matrix.
   Uses the same settings FEBio's PardisoSolver uses (mtype -2, default iparm, 1 refinement step). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl.h>
#include <omp.h>
static double now(){ return omp_get_wtime(); }
int main(int argc, char** argv){
  FILE* f = fopen(argv[1], "rb"); long long n, nnz;
  fread(&n,8,1,f); fread(&nnz,8,1,f);
  int* ia = malloc(sizeof(int)*(n+1)); int* ja = malloc(sizeof(int)*nnz); double* a = malloc(sizeof(double)*nnz);
  fread(ia,4,n+1,f); fread(ja,4,nnz,f); fread(a,8,nnz,f); fclose(f);
  void* pt[64] = {0}; int iparm[64] = {0}; int mtype=-2, maxfct=1, mnum=1, nrhs=1, msglvl=0, err=0, N=(int)n, phase;
  pardisoinit(pt, &mtype, iparm);
  iparm[18]=0; iparm[34]=1; /* zero-based */
  double t=now(); phase=11;
  pardiso(pt,&maxfct,&mnum,&mtype,&phase,&N,a,ia,ja,NULL,&nrhs,iparm,&msglvl,NULL,NULL,&err);
  printf("pardiso analyze %.3fs err=%d nnz(L)=%d\n", now()-t, err, iparm[17]);
  for(int r=0;r<3;r++){ t=now(); phase=22;
    pardiso(pt,&maxfct,&mnum,&mtype,&phase,&N,a,ia,ja,NULL,&nrhs,iparm,&msglvl,NULL,NULL,&err);
    printf("pardiso factor[%d] %.3fs err=%d perturbed=%d\n", r, now()-t, err, iparm[13]); }
  double *b=malloc(8*n), *x=malloc(8*n), *xt=malloc(8*n);
  for(long i=0;i<n;i++){ xt[i]=((i*7919)%101)/101.0-0.5; b[i]=0; }
  for(long j=0;j<n;j++) for(int p=ia[j];p<ia[j+1];p++){ int i=ja[p]; b[i]+=a[p]*xt[j]; if(i!=j) b[j]+=a[p]*xt[i]; }
  iparm[7]=1; t=now(); phase=33;
  pardiso(pt,&maxfct,&mnum,&mtype,&phase,&N,a,ia,ja,NULL,&nrhs,iparm,&msglvl,b,x,&err);
  double ts=now()-t, e=0, xn=0; for(long i=0;i<n;i++){ e+=(x[i]-xt[i])*(x[i]-xt[i]); xn+=xt[i]*xt[i]; }
  printf("pardiso solve %.4fs rel error %.3e (refine steps %d)\n", ts, sqrt(e/xn), iparm[6]);
  phase=-1; pardiso(pt,&maxfct,&mnum,&mtype,&phase,&N,a,ia,ja,NULL,&nrhs,iparm,&msglvl,NULL,NULL,&err);
  return 0; }
