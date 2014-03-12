/*
The first test function is a very simple one
*/
#include <math.h>
#include <adolc/adolc.h>

#define NUM_IND	6

int get_num_ind(){
  return NUM_IND;
}
void get_initial_value(double *x){
  int i;
  for(i=0;i<NUM_IND;i++){
    x[i]=(double)(i+2);
  }
}
adouble func_eval(adouble *x){
  adouble fad=0;
  fad =0.5*(x[0] - 1)*(x[0] -1) + 0.8*(x[1] - 2)*(x[1] -2)  + 0.9*(x[2] - 3)*(x[2] -3);
  fad =fad+ cosh(x[3]);
  fad+=sinh(x[4])*x[1]*erf(x[1]);
  fad =fad+ asinh(x[5])*acosh(x[2]);
  fad =fad+ atanh(x[4]*x[5]);
  return fad;
}

