/*
Structured Jacobian problem
*/
#include <math.h>
#include <adolc/adolc.h>

#define NUM_IND	10000

#define max(a,b) ((a>b)?a:b)
#define min(a,b) ((a<b)?a:b)
int get_num_ind(){
  return NUM_IND;
}
void get_initial_value(double *x){
  int i;
  for(i=0;i<NUM_IND;i++){
      x[i]=-1;
  }
}
adouble func_eval(adouble *x){
  adouble fad=0;
  adouble t=0;
  int n=NUM_IND;
  int k;
  adouble l=-x[n-3]+0.5*x[n-2]-x[n-1];
  for(k=0;k<n;k++){
    t=-2*x[k]*x[k]+3*x[k]+l;
    if (k!=0) { t+=-x[k-1]; }
    if (k!=n-1) { t+=-2*x[k+1]; }
    fad+=t*t;
  }
  return fad;
}

