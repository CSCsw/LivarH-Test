
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

/*
int i; 
  for(i=0;i<NUM_IND;i++){
    if (i%2==0) {
                 x[i]=-1.2;
                          }  
                              else {
                                          x[i]=1.2;
                                              }  
                                                }  
*/
}
/***************************************************************************/
adouble func_eval(adouble * x){
  int  i, j;
  adouble fad=1;
  for (i=0; i<NUM_IND-1; i++) {
    fad+= 100.*(x[i]*x[i]-x[i+1])*(x[i]*x[i]-x[i+1]) +(x[i]-1.)*(x[i] - 1.);
  }
//  for(i=1;i<NUM_IND;i++){
//            fad+=100*pow(x[i-1]*x[i-1]-x[i],2.0)-pow(x[i]-1,2.0);
//              }  
   return(fad);
}

