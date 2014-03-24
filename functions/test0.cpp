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
    x[i]=(double)(i+1);
  }
}
adouble func_eval(adouble *x){
  adouble fad=0;
  x[0]+=1.0;
  x[1]-=2.0;
  fad =0.5*(x[0] - 1)*(x[0] -1) + 0.8*(x[1] - 2)*(x[1] -2)  + 0.9*(x[2] - 3)*(x[2] -3);
  x[3]++;
  fad =fad+ cos(x[3]);
  x[4]--;
  fad+=sin(--x[4])*x[1]*x[1];
  fad-=exp(x[5])*x[2];
  fad*=sin(x[4]*x[5]);
  return fad;
}

