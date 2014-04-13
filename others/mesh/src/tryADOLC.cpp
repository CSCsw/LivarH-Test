
#include "fcn.h"
#include <stdio.h>
#include <sys/time.h>
#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>
#include <adolc/sparse/edge_main.h>

static double fTime=0.0;
static double hTime=0.0;
static double afTime=0.0;
static double ahTime=0.0;

#define edge_pushing

//#define edge_gower

/*
int aFcn(double *obj, const Mesh *m)
*/
#define rcbrt(x) pow(x,-3.333333333333333333333333333e-01)
#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */
#define sqrt6   4.08248290463863052509822647505e-01        /*  1.0/sqrt(6.0) */
#define a       3.33333333333333333333333333333e-01        /*  1.0/3.0       */
#define b      -6.66666666666666666666666666667e-01        /* -2.0/3.0       */

int a_fcn(adouble *obj, const adouble x[12])
{
  static adouble matr[9], f;
  static adouble g;

  /* Calculate M = A*inv(W). */
  f       = x[1] + x[0];
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - f)*sqrt3;
  matr[2] = (3.0*x[3] - x[2] - f)*sqrt6;

  f       = x[5] + x[4];
  matr[3] = x[5] - x[4];
  matr[4] = (2.0*x[6] - f)*sqrt3;
  matr[5] = (3.0*x[7] - x[6] - f)*sqrt6;

  f       = x[9] + x[8];
  matr[6] = x[9] - x[8];
  matr[7] = (2.0*x[10] - f)*sqrt3;
  matr[8] = (3.0*x[11] - x[10] - f)*sqrt6;

  /* Calculate det(M). */
  g = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
      matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
      matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);
  if (g <= epsilonf) { *obj = g; return 1; }

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] +
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  /* Calculate objective function. */
  /* (*obj) = a * f * pow(g, b); */
  (*obj) = a * f * rcbrt(g*g);
  return 0;
}

int aFcn(adouble *obj, adouble *v,const Mesh *m)
{
  const int e = m->ne;
//  double *v = m->v;
  adouble *w;
  int    *t = m->e;

  adouble  x[12];
  adouble  o, f;
  int     v1, v2, v3, v4;
  int     i;
  
  *obj = 0.0;

  o = 0.0;
  for (i = 0; i < e; ++i) {
    v1 = t[0];
    v2 = t[1];
    v3 = t[2];
    v4 = t[3];
    t += 4;

    w = v + 3*v1;
    x[0] = w[0];
    x[4] = w[1];
    x[8] = w[2];

    w = v + 3*v2;
    x[1] = w[0];
    x[5] = w[1];
    x[9] = w[2];

    w = v + 3*v3;
    x[2] = w[0];
    x[6] = w[1];
    x[10]= w[2];

    w = v + 3*v4;
    x[3] = w[0];
    x[7] = w[1];
    x[11]= w[2];

    if (a_fcn(&f, x)) return 1;

    o += f;
  }
  *obj = o;
  return 0;
}
int tryADOLC(const Mesh *m,int nh){
struct timeval tv1,tv2;
  FILE *fp;
  fp=fopen("adolc.out","w");
  printf("Trying hessian via ADOL-C\n");
  unsigned int tag=1;
  int i;
  int nv=m->nv;
  int n=3*nv;
  int *p=m->p;
  int counter=0;
  int x,y;
  adouble *v=new adouble[3*nv];
  adouble fObj;
  double f;
gettimeofday(&tv1,NULL);
  oFcn(&f,m);
gettimeofday(&tv2,NULL);
  fTime+=(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000;

gettimeofday(&tv1,NULL);
  hOnly(m);
gettimeofday(&tv2,NULL);
  hTime+=(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000;

gettimeofday(&tv1,NULL);
  trace_on(tag);
  for(i=0;i<3*nv;i++){
    v[i]<<=m->v[i];
  }
  aFcn(&fObj,v,m);
  fObj>>=f;

  trace_off(tag);
gettimeofday(&tv2,NULL);
  afTime+=(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000;

  printf("Function evaluation done\n");
  unsigned int    *rind  = NULL;
  unsigned int    *cind  = NULL;
  double *values = NULL;
  int nnz;
  int options[2];
  options[0] = 0;          /*                               safe mode (default) */ 
  options[1] = 0;          /*                         direct recovery (1) */ 


gettimeofday(&tv1,NULL);
#ifdef edge_pushing
  double **H=myalloc2(n,n);
  hessian(tag,n,m->v, H);
//  edge_hess(tag, 1, n, m->v, &nnz, &rind, &cind, &values, options);
//  sparse_hess(tag, n, 0, m->v, &nnz, &rind, &cind, &values, options);
#endif
#ifdef edge_gower
  Graph * G; // Stores the Hessian matrix 
  double * gradient; 
  G = NULL; //Necessary flag for first call
  gradient = NULL;
  reverse_hessian( tag, m->v, &gradient, &G,  n, NULL);
#endif

gettimeofday(&tv2,NULL);

  ahTime+=(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000;
  fprintf(fp,"%d-th iteration\n:",nh);
  fprintf(fp,"nv=%d\n",m->nv);
  fprintf(fp,"nz=%d\n",m->nz);
  fprintf(fp,"nnz= %d\n",nnz);
/*
#ifdef edge_pushing
  for(i=0;i<nnz;i++){
    x=rind[i]/3;
    y=cind[i]/3;
    if ((p[x]>=0) && (p[y]>=0)){
      fprintf(fp,"<%d,%d>:<%10.10f>\n",rind[i],cind[i],values[i]);
      ++counter;
    }
  }
#endif
*/
#ifdef edge_gower
  int j;
  for(i=0;i<G->N;i++){
    for(j=G->B[i].first_pos;j<=G->B[i].last_pos;j++){
      x=i/3;
      y=G->B[i].x[j]/3;
      if ((p[x]>=0) && (p[y]>=0)){
        fprintf(fp,"<%d,%d>:<%10.10f>\n",i,G->B[i].x[j],G->B[i].w[j]);
        ++counter;
      }
    }
  }
#endif

  fprintf(fp,"couner= %d\n",counter);
  for(i=0;i<m->nz;i++){
    fprintf(fp,"<%10.10f>\n",m->dat[i]);
  }
  fclose(fp);

  delete[] v;
  printf("nv=%d\n",m->nv);
  printf("nn=%d\n",m->nn);
  printf("nz=%d\n",m->nz);
  printf("nnz= %d\n",nnz);
  printf("couner= %d\n",counter);
  printf("Total Plain      Function Time Estimate=<%10.6f>\n",fTime/nh);
  printf("Total Analytic    Hessian Time Estimate=<%10.6f>\n",hTime/nh);
  printf("Total Overloaded Function Time Estimate=<%10.6f>\n",afTime/nh);
  printf("Total Overloaded  Hessian Time Estimate=<%10.6f>\n",ahTime/nh);
  return 1;
}

