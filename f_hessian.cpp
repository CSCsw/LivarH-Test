
#include <math.h>
#include <cstdlib>
#include <cstdio>

#include <adolc/adolc.h>

#define NUM_IND 100

int get_num_ind(){
  return NUM_IND;
}

void get_initial_value(double *x){
  int i;

  for(i=0;i<NUM_IND;i++){
    x[i]=i+1;
  }

}
/***************************************************************************/
adouble func_eval(adouble * x){
  int  i, j;
//  adouble fad=x[0]*x[1]*x[2];
  adouble fad=0;  
  for (i=0; i<NUM_IND-1; i++) {
    fad+= 100.*(x[i]*x[i]-x[i+1])*(x[i]*x[i]-x[i+1]) +(x[i]-1.)*(x[i] - 1.);
  }
   return(fad);
}

