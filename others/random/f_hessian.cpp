
#include <math.h>
#include <cstdlib>
#include <cstdio>

#include <adolc/adolc.h>

#define NUM_IND 10

#define K_NNZ 1

#define EXP_1 x[i]*x[r]
#define EXP_2 EXP_1+EXP_1
#define EXP_4 EXP_2+EXP_2

#define LINE_1 fad+=x[i]*x[r];
#define LINE_2 LINE_1;LINE_1;
#define LINE_4 LINE_2;LINE_2;
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
  srand(12345);
  int  i, j, k;
  int r;
  adouble fad=0;  
  for(i = 0; i < NUM_IND; ++i) {
    for(j = 0; j < K_NNZ; j++) {
      r = rand() % NUM_IND;
      printf("r = %d\n", r);
//      fad += x[i]*x[r];
//      fad += EXP_4;
      LINE_4;
    }
  }
  return(fad);
}

