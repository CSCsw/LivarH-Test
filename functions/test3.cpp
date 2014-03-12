/*
Test condassign & condassign_s in ADOL-C
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
  x[0]=x[0]-x[1];
  x[1]=x[2]-x[1];
  condassign(fad,x[0],x[3]*x[4],x[4]*x[5]);
  condassign(x[0],x[1],x[2]*x[3],x[3]*x[4]);
  fad-=x[0]*x[3];
  condassign(x[1],x[5],x[2]+x[3]);
  fad/=x[1];
  return fad;
}

