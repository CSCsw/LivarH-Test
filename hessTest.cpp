#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <sys/time.h>

#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>
#include <adolc/hessian/edge_main.h>

#define COMPUT_GRAPH 1
#define PRE_ACC 0

#define tag 1
//#define edge_pushing 1
#define _compare_with_full

//#define _PRINTOUT
#define def_tol (0.00001)


extern int get_num_ind();
extern void get_initial_value(double *x);
extern adouble func_eval(adouble *x);

//void printmat(char* kette, int n, int m, double** M);
void compare_matrix(int n, double** H, int nnz, unsigned int *r, unsigned int *c, double *v){
  int i,j;
  for(i=0;i<nnz;i++){
    H[r[i]][c[i]]-=v[i];
  }
for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if (fabs(H[i][j])>def_tol){
        printf("WRONG!\n");
        exit(-1);
        return;
      }
    }
  }
  printf("CORRECT!\n");
}
int main(int argc, char *argv[]) {
  int n=get_num_ind();
  int i,j;
  struct timeval tv1,tv2;
  adouble *xad;
  adouble fad;
  double f;
  double *x;
  x=new double[n];
  xad=new adouble[n];
get_initial_value(x);

  printf("evaluating the function...");
trace_on(tag);
  for(i=0;i<n;i++)
  {
    xad[i] <<= x[i];  
  }
  fad=func_eval(xad); 
  fad >>= f;
trace_off();
  printf("done!\n");
//  printf("function value  =<%10.20f>\n",f);
//  function(tag,1,n,x,&f);
//  printf("adolc func value=<%10.20f>\n",f);
//tape_doc(tag,1,n,x,&f);
#ifdef _compare_with_full
  double **H;
  H = myalloc2(n,n);
  printf("computing full hessain....");
  gettimeofday(&tv1,NULL);
  hessian(tag,n,x,H);
  printf("done\n");
  gettimeofday(&tv2,NULL);
  printf("Computing the full hessian cost %10.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);
#ifdef _PRINTOUT
    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
        printf("H[%d][%d]=<%10.10f>",i,j,H[i][j]);
      }
      printf("\n");
    }
    printf("\n");
#endif
#endif

#ifdef edge_pushing
  unsigned int    *rind  = NULL;
  unsigned int    *cind  = NULL;
  double *values = NULL;
  int nnz;
  int options[2];
  options[0]=PRE_ACC;
  options[1]=COMPUT_GRAPH;
  gettimeofday(&tv1,NULL);
//  edge_hess(tag, 1, n, x, &nnz, &rind, &cind, &values, options);
  sparse_hess(tag,n,0,x, &nnz, &rind, &cind, &values, options);
  gettimeofday(&tv2,NULL);
  printf("Sparse Hessian: edge pushing cost %10.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);

#ifdef _PRINTOUT
  for(i=0;i<nnz;i++){
    printf("<%d,%d>:<%10.10f>\n",cind[i],rind[i],values[i]);
//    printf("%d %d \n", rind[i], cind[i]);
  }
#endif
#endif

#ifdef _compare_with_full
#ifdef edge_pushing
  compare_matrix(n,H,nnz,cind,rind,values);
#endif
  myfree2(H);
#endif

#ifdef edge_pushing
  printf("nnz=%d\n", nnz);
  free(rind); rind=NULL;
  free(cind); cind=NULL;
  free(values); values=NULL;
#endif
  delete[] x;
  delete[] xad;
  return 0;
}


