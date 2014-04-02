/*
Test advector in ADOL-C
As R-value, have to explicitly cast?
*/
#include <math.h>
#include <adolc/adolc.h>
#include <adolc/advector.h> 

#define NUM_IND	6

int get_num_ind(){
  return NUM_IND;
}
void get_initial_value(double *x){
  int i;
  for(i=0;i<NUM_IND;i++){
    x[i]=(double)(i);
  }
}
adouble func_eval(adouble *x){
  adouble fad=0;
  x[2]=(x[1]>x[0])*x[1];
  fad=x[2]*x[1];
  return fad;
}


