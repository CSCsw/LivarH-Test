/*----------------------------------------------------------------------------
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     sparse_hessian.cpp
 Revision: $Id: sparse_hessian.cpp 203 2011-02-23 11:33:40Z kulshres $
 Contents: example for computation of sparse hessians

 Copyright (c) Andrea Walther, Andreas Griewank, Andreas Kowarz, 
               Hristo Mitev, Sebastian Schlenkrich, Jean Utke, Olaf Vogel
  
 This file is part of ADOL-C. This software is provided as open source.
 Any use, reproduction, or distribution of the software constitutes 
 recipient's acceptance of the terms of the accompanying license file.
 
---------------------------------------------------------------------------*/

#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <sys/time.h>

#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>
#include <adolc/sparse/edge_main.h>

#define COMPUT_GRAPH 0
#define PRE_ACC 0

#define tag 1
//#define _SAPARATE 1
#define XXX 1
//#define direct 1
#define indirect 1
//#define edge_pushing 1
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

int main(int argc, char *argv[]) {
//    int n=6;
//    double f, x[6];
//    adouble fad, xad[6];
//    int i, j;

/****************************************************************************/
/*******                function evaluation                   ***************/
/****************************************************************************/
    struct timeval start, end;
    long mtime, seconds, useconds;    
    int  n,m;	int i,j;
    long double c0,c1, Mytime;
struct timeval tv1,tv2;
    n=get_num_ind();
    //----------------------------------------*/
// Variables associated to function calculation and the functions tape!
/*--------------------------------------------------------------------------*/
        adouble * xad;
        xad = new adouble [n];
	adouble fad;
	double f;
	double *x = new double [n]; 
/*---------------Taping the function---------------------------------------*/

get_initial_value(x);
	
printf("evaluating the function...");
     trace_on(tag);

     for(i=0;i<n;i++) 
	  {
		xad[i] <<= x[i];				  // active independs        //  
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

#ifdef edge_pushing
    unsigned int    *rind1  = NULL;
    unsigned int    *cind1  = NULL;
    double *values1 = NULL;
//    double **J=new double* [1];
//    J[0]= new double [n];
    int nnz1;
    int options1[2];

    options1[0] = 0;          /*                               safe mode (default) */ 
    options1[1] = 2;          /*                                   edge_pushing*/
#endif

#ifdef direct 
printf("\n**********Safe Mode, Direct recovery*******************\n");

gettimeofday(&tv1,NULL);
    sparse_hess(tag, n, 0, x, &nnz, &rind, &cind, &values, options);
gettimeofday(&tv2,NULL);
printf("Sparse Hessian: direct recovery cost %10.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);
#ifdef _PRINTOUT
  for(i=0;i<nnz0;i++){
    printf("<%d,%d>:<%10.10f>\n",rind[i],cind[i],values[i]);
  }
#endif
#endif

#ifdef indirect 
printf("\n**********Safe Mode, indirect recovery*******************\n");

gettimeofday(&tv1,NULL);
    sparse_hess(tag, n, 0, x, &nnz0, &rind0, &cind0, &values0, options0);
gettimeofday(&tv2,NULL);
printf("Sparse Hessian: indirect recovery cost %10.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);
#ifdef _PRINTOUT
  for(i=0;i<nnz0;i++){
    printf("<%d,%d>:<%10.10f>\n",rind0[i],cind0[i],values0[i]);
  }
#endif
#endif

#ifdef edge_pushing
printf("\n**********Safe Mode, Edge Pushing*******************\n");

    options1[0]=PRE_ACC;
    options1[1]=COMPUT_GRAPH;
gettimeofday(&tv1,NULL);
    edge_hess(tag, 1, n, x, &nnz1, &rind1, &cind1, &values1, options1);
//    jacobian(tag,1,n,x,J);
gettimeofday(&tv2,NULL);
printf("Sparse Hessian: edge pushing cost %10.6f seconds\n",(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000);

#ifdef _PRINTOUT
  for(i=0;i<nnz1;i++){
    printf("<%d,%d>:<%10.10f>\n",rind1[i],cind1[i],values1[i]);
  }
/*
  for(i=0;i<n;i++){
    printf("Adjoints[%d]=<%10.10f>\n",i,J[0][i]);
  }
*/
#endif
#endif

#ifdef direct
    free(rind); rind=NULL;
    free(cind); cind=NULL;
    free(values); values=NULL;
#endif

#ifdef indirect
    free(rind0); rind0=NULL;
    free(cind0); cind0=NULL;
    free(values0); values0=NULL;
#endif
#ifdef edge_pushing
    free(rind1); rind1=NULL;
    free(cind1); cind1=NULL;
    free(values1); values1=NULL;
#endif

  delete[] x;
  delete[] xad;
}


