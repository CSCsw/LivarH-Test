
#include "fcn.h"
#include <stdio.h>
#include <sys/time.h>
#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>
#include <adolc/hessian/edge_main.h>

static double fTime=0.0;
static double hTime=0.0;
static double afTime=0.0;
static double ahTime=0.0;

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

  g = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
      matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
      matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);
  if (g <= epsilonf) { *obj = g; return 1; }

  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] +
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

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
int main(int argc, char **argv){
    int nh=1;
    Mesh *m = NULL;

//Write result to file or not
    int outputFlag=0;
    int hessFlag=0;
    int options[2];
    options[0]=0;
    options[1]=1;

    char meshFile[20]="gear.mesh";
    char hessFile[20]="hess.out";
    char edgeFile[20]="edge.out";
    if (argc>1){
        printf("Reading mesh data=%d, analytic_hess>>%s, edge_hess>>%s (should be the same)\n",meshFile, hessFile,edgeFile);
        outputFlag=1;
        if (argc>2){
            hessFlag=atoi(argv[2]);
        }
        if (argc>4){
            options[0]=atoi(argv[3]);
            options[1]=atoi(argv[4]);
        }
    }
    if (readMesh(meshFile, &m)) {
        freeMesh(&m);
        return -1;
    }
    hMesh(m);
    printf("read mesh done..\n");
    struct timeval tv1,tv2;
    FILE *fp1;
    FILE *fp2;
    if (outputFlag==1){
        fp1=fopen(hessFile,"w");
        fp2=fopen(edgeFile,"w");
    }
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
//Active Section
    trace_on(tag);
    for(i=0;i<3*nv;i++){
        v[i]<<=m->v[i];
    }
    aFcn(&fObj,v,m);
    fObj>>=f;

    trace_off(tag);
//Done
gettimeofday(&tv2,NULL);
    afTime+=(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000;

    printf("Function evaluation done\n");
    unsigned int    *rind  = NULL;
    unsigned int    *cind  = NULL;
    double *values = NULL;
    int nnz;
gettimeofday(&tv1,NULL);
    if (hessFlag==0){
        edge_hess(tag, 1, n, m->v, &nnz, &rind, &cind, &values, options);
    }
    else{
        sparse_hess(tag, n, 0, m->v, &nnz, &rind, &cind, &values, options);
    }
gettimeofday(&tv2,NULL);
    ahTime+=(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000;
    if (outputFlag==1){
        for(i=0;i<nnz;i++){
            x=rind[i]/3;
            y=cind[i]/3;
            if ((p[x]>=0) && (p[y]>=0)){
//                fprintf(fp,"<%d,%d>:<%10.10f>\n",rind[i],cind[i],values[i]);
                fprintf(fp2,"%15.8f\n",values[i]);
                ++counter;
            }
        }
        for(i=0;i<m->nz;i++){
            fprintf(fp1,"%15.8f\n",m->dat[i]);
        }
        fclose(fp1);
        fclose(fp2);
    }
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
    return 0;
}
int tryADOLC(const Mesh *m,int nh){
}
