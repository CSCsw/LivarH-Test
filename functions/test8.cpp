/*
Generalized Broyden banded function
*/
#include <math.h>
#include <adolc/adolc.h>

#define NUM_IND	10

#define max(a,b) ((a>b)?a:b)
#define min(a,b) ((a<b)?a:b)
int get_num_ind(){
  return NUM_IND;
}
void get_initial_value(double *x){
  int i;
  for(i=0;i<NUM_IND;i++){
      x[i]=sin(i);
  }
}
adouble func_eval(adouble *x){
  adouble fad=0;
  adouble t=0;
  int i,j,k;
  int m=5*NUM_IND;
  int b,c;
  for(k=0;k<m;k++){
    b=5-(k*4)/(m);
    c=k%5+1;
    i=k%(NUM_IND/2);
    j=i+NUM_IND/2;
    t=pow(pow(x[i],1.0)-pow(x[j],b),c);
    fad+=t*t;
  }
  return fad;
}

