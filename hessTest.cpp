#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <sys/time.h>

#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>
#include <adolc/hessian/edge_main.h>


#define tag 1


#define LIVARH
//#define LIVARHACC
//#define DIRECT
//#define INDIRECT

#define COMPARE_WITH_FULL_HESS

//#define PRINT_RESULTS

#define def_tol (0.000001)


// Basically, one can implement his own objective function by implementing
// the following 3 functions:
// Get the number of independent variables
extern int get_num_ind();
// Get the initial values of independent variables
extern void get_initial_value(double *x);
// The objective function
extern adouble func_eval(adouble *x);

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

#ifdef COMPARE_WITH_FULL_HESS
  double **H;
  H = myalloc2(n,n);
  printf("computing full hessain....");
  gettimeofday(&tv1,NULL);
  hessian(tag,n,x,H);
  printf("done\n");
  gettimeofday(&tv2,NULL);
  printf("Computing the full hessian cost %10.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);

#ifdef PRINT_RESULTS
    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
        printf("H[%d][%d]=<%10.10f>",i,j,H[i][j]);
      }
      printf("\n");
    }
    printf("\n");
#endif
#endif

  unsigned int    *rind  = NULL;
  unsigned int    *cind  = NULL;
  double *values = NULL;
  int nnz;
  int options[2];

#ifdef LIVARH
  options[0]=0;
  options[1]=1;
  gettimeofday(&tv1,NULL);
  edge_hess(tag, 1, n, x, &nnz, &rind, &cind, &values, options);
  gettimeofday(&tv2,NULL);
  printf("Sparse Hessian: LivarH cost %10.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);
#endif

#ifdef LIVARHACC
  options[0]=1;
  options[1]=1;
  gettimeofday(&tv1,NULL);
  edge_hess(tag, 1, n, x, &nnz, &rind, &cind, &values, options);
  gettimeofday(&tv2,NULL);
  printf("Sparse Hessian: LivarHACC cost %10.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);
#endif

// Sparse ADOL-C drivers report the upper matrix
#ifdef DIRECT
  options[0]=0;
  options[1]=1;
  gettimeofday(&tv1,NULL);
  sparse_hess(tag, n, 0, x, &nnz, &cind, &rind, &values, options);
  gettimeofday(&tv2,NULL);
  printf("Sparse Hessian: direct recovery cost %10.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);
#endif

#ifdef INDIRECT
  options[0]=0;
  options[1]=0;
  gettimeofday(&tv1,NULL);
  sparse_hess(tag, n, 0, x, &nnz, &cind, &rind, &values, options);
  gettimeofday(&tv2,NULL);
  printf("Sparse Hessian: indirect recovery cost %10.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);
#endif

#ifdef PRINT_RESULTS
  for(i=0;i<nnz;i++){
    printf("<%d,%d>:<%10.10f>\n",rind[i],cind[i],values[i]);
  }
#endif

#ifdef COMPARE_WITH_FULL_HESS
  compare_matrix(n,H,nnz,rind,cind,values);
  myfree2(H);
#endif

  free(rind); rind=NULL;
  free(cind); cind=NULL;
  free(values); values=NULL;

  delete[] x;
  delete[] xad;
  return 0;
}


