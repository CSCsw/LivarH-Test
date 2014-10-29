/*
Generalized Broyden banded function
F2 in the paper
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
      x[i]=-1.2;
  }
}
adouble func_eval(adouble *x){
  adouble fad=0;
  adouble t;
  int i,j;
  for(i=0;i<NUM_IND;i++){
    t=(2+5*x[i]*x[i])*x[i]+1;
    for(j=max(0,i-5);j<min(NUM_IND,i+1);j++) {
      if (i!=j) {
          t=t+x[j]*(x[j]+1);
      }
    }
    fad+=t*t;
  }
  return fad;
}

