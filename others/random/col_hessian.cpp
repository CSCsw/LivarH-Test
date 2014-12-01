
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

//#define direct 1
//#define indirect 1
#define LIVARHNO 1
//#define LIVARHACC 1
//#define _compare_with_full

//#define _PRINTOUT
#define def_tol (0.00001)

double hessian_crc;
int sparsity_crc_row;
int sparsity_crc_col;

extern int get_num_ind();
extern void get_initial_value(double *x);
extern adouble func_eval(adouble *x);

void printmat(char* kette, int n, int m, double** M);

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
/****************************************************************************/
/*******                function evaluation                   ***************/
/****************************************************************************/
    struct timeval start, end;
    long mtime, seconds, useconds;    
    int  n,m;	int i,j;
    long double c0,c1, Mytime;
struct timeval tv1,tv2;
    n=get_num_ind();

// Variables associated to function calculation and the functions tape!
    adouble * xad;
    xad = new adouble [n];
    adouble fad;
    double f;
    double *x = new double [n]; 
/*---------------Taping the function---------------------------------------*/

get_initial_value(x);
    trace_on(tag);
    for(i=0;i<n;i++) {
      xad[i] <<= x[i];				  // active independs
    } 
    fad = func_eval(xad);
    fad >>= f;
    trace_off();
#ifdef _compare_with_full
    double **H;
    H = myalloc2(n,n);
printf("computing full hessain....");
gettimeofday(&tv1,NULL);
    hessian(tag,n,x,H);
printf("done\n");
gettimeofday(&tv2,NULL);
printf("Computing the full hessian cost %.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);
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


/****************************************************************************/
/*******       sparse Hessians, complete driver              ***************/
/****************************************************************************/
#ifdef direct
    unsigned int    *rind  = NULL;
    unsigned int    *cind  = NULL;
    double *values = NULL;
    int nnz;
    int options[2];
    options[0] = 0;          /*                               safe mode (default) */ 
    options[1] = 1;          /*                         direct recovery (1) */ 
#endif

#ifdef indirect
    unsigned int    *rind0  = NULL;
    unsigned int    *cind0  = NULL;
    double *values0 = NULL;
    int nnz0;
    int options0[2];
    options0[0] = 0;          /*                               safe mode (default) */ 
    options0[1] = 0;          /*                       indirect recovery (default) */ 
#endif

#ifdef LIVARHNO
    unsigned int    *rind1  = NULL;
    unsigned int    *cind1  = NULL;
    double *values1 = NULL;
    int nnz1;
    int options1[2];
#endif

#ifdef LIVARHACC
    unsigned int    *rind2  = NULL;
    unsigned int    *cind2  = NULL;
    double *values2 = NULL;
    int nnz2;
    int options2[2];
#endif

#ifdef direct 
gettimeofday(&tv1,NULL);
    sparse_hess(tag, n, 0, x, &nnz, &rind, &cind, &values, options);
gettimeofday(&tv2,NULL);
printf("Sparse Hessian: STAR cost %.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);
#ifdef _PRINTOUT
  for(i=0;i<nnz0;i++){
    printf("<%d,%d>:<%10.10f>\n",rind[i],cind[i],values[i]);
  }
#endif
  free(rind); rind=NULL;
  free(cind); cind=NULL;
  free(values); values=NULL;
#endif

#ifdef indirect 

gettimeofday(&tv1,NULL);
    sparse_hess(tag, n, 0, x, &nnz0, &rind0, &cind0, &values0, options0);
gettimeofday(&tv2,NULL);
printf("Sparse Hessian: ACYCLIC cost %.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);
#ifdef _PRINTOUT
  for(i=0;i<nnz0;i++){
    printf("<%d,%d>:<%10.10f>\n",rind0[i],cind0[i],values0[i]);
  }
#endif
  free(rind0); rind0=NULL;
  free(cind0); cind0=NULL;
  free(values0); values0=NULL;
#endif

#ifdef LIVARHNO

    options1[0]=0;
    options1[1]=1;
gettimeofday(&tv1,NULL);
    edge_hess(tag, 1, n, x, &nnz1, &rind1, &cind1, &values1, options1);
gettimeofday(&tv2,NULL);
printf("Sparse Hessian: LIVARHNO cost %.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);
#ifdef _PRINTOUT
  for(i=0;i<nnz1;i++){
    printf("<%d,%d>:<%10.10f>\n",rind1[i],cind1[i],values1[i]);
  }
#endif
#ifdef _compare_with_full
    compare_matrix(n, H, nnz1, rind1, cind1, values1);
#endif
    int *counter = (int*)malloc(sizeof(int)*n);
    for (i=0;i<n;i++) {
        counter[i]=0;
    }
    for(i=0; i<nnz1; i++) {
      counter[rind1[i]]++;
      if (rind1[i] != cind1[i]) {
          counter[cind1[i]]++;
      }
    }
    int max=0;
    int min=n;
    double s=0;
    for(i=0;i<n;i++) {
      max = (max>counter[i])?max:counter[i];
      min = (min<counter[i])?min:counter[i];
      s=s+counter[i];  
    }
    printf("#nnz: max=%d, min=%d, avg=%.5f\n", max, min, s/n);

  free(rind1); rind1=NULL;
  free(cind1); cind1=NULL;
  free(values1); values1=NULL;
#endif

#ifdef LIVARHACC

    options2[0]=1;
    options2[1]=1;
gettimeofday(&tv1,NULL);
    edge_hess(tag, 1, n, x, &nnz2, &rind2, &cind2, &values2, options2);
gettimeofday(&tv2,NULL);
printf("Sparse Hessian: LIVARHACC cost %.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);
#ifdef _PRINTOUT
  for(i=0;i<nnz2;i++){
    printf("<%d,%d>:<%10.10f>\n",rind2[i],cind2[i],values2[i]);
  }
#endif
#ifdef _compare_with_full
    compare_matrix(n, H, nnz2, rind2, cind2, values2);
#endif
  free(rind2); rind2=NULL;
  free(cind2); cind2=NULL;
  free(values2); values2=NULL;
#endif

  delete[] x;
  delete[] xad;
}


