/*
This test contains some special cases.
Where the monotonic assumption does not hold.
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
  x[0]=x[0]/x[0]; //Then x[0]=1.0 const
  x[1]=x[0]*x[4]*x[5];
  x[3]=sin(x[1])*cosh(x[3]);
  fad=pow(x[0],x[1]);
  fad+=x[1];
  fad++;
  fad-=x[0]*(-x[2]);
  fad--;
  fad=sqrt(fad);
  return fad;
}

