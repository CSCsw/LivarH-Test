/*
Chained Rosenbrock function
*/
#include <math.h>
#include <adolc/adolc.h>

#define NUM_IND	20000

int get_num_ind(){
  return NUM_IND;
}
void get_initial_value(double *x){
  int i;
  for(i=0;i<NUM_IND;i++){
    if (i%2==0) {
        x[i]=-1.2;
    }
    else {
        x[i]=1.2;
    }
  }
}
adouble func_eval(adouble *x){
  adouble fad=0;
  int i;
  for(i=1;i<NUM_IND;i++){
      fad+=100*pow(x[i-1]*x[i-1]-x[i],2.0)-pow(x[i]-1,2.0);
  }
  return fad;
}

