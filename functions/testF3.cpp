/*
Potra and Rheinboldt boundary value problem
F3 in the paper
*/
#include <math.h>
#include <adolc/adolc.h>

#define NUM_IND	200

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
  for(k=1;k<n-1;k++){
    if (k < n/2) {
        t=2*x[k]-x[k-1]-x[k+1]+(pow(x[k],2)+x[k]+0.1+x[k+n/2]-1.2);
    } else {
        t=2*x[k]-x[k-1]-x[k+1]+(0.2*pow(x[k-n/2],2)+x[k]*x[k]+2*x[k]-0.6);
    }
    fad+=t*t;
  }
  return fad;
}

