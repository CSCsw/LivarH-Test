/*
Advanced Branching
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
  adouble t1,t2;
  t1=x[3]*x[4];
  t2=x[4]*x[5];
  x[0]=(x[0]!=x[0])*t1;
  x[3]=sin((x[3]<=x[4])*t2);
  x[2]=(x[1]+(x[1]>=x[5]))*x[3]*x[4];
  x[4]=x[0]*x[1]*x[2];
  x[5]=x[3]*x[4]*x[5];
  fad=fad+x[0]+x[1]+x[2]+x[3]+x[4]+x[5];
  return fad;
}

