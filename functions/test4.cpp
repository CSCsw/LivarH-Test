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
  advector adv(NUM_IND);
  adv[x[0]]=x[0];
  adv[x[1]]=x[1];
  adv[x[2]]=x[2];
  condassign(adv[x[3]],x[0]-x[1],x[0],(adouble)adv[x[1]]);
  condassign(adv[x[4]],x[2]-x[1],(adouble)adv[x[2]],(adouble)adv[x[3]]);
  condassign(adv[x[5]],x[3]-x[2],x[5]);
  adv[x[0]]+=cosh((adouble)adv[x[0]]);
  adv[x[1]]-=0.2;
  adv[x[1]-x[0]]=exp((adouble)adv[x[5]]);
  fad+=adv[0]*adv[1]*adv[2]*adv[3]*adv[4]*adv[5];
  return fad;
}


