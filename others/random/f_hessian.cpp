
#include <math.h>
#include <cstdlib>
#include <cstdio>

#include <adolc/adolc.h>

//#define NUM_IND 20000

//#define K_NNZ 3

#define EXP_1 1.0+x[i]*x[r]
#define EXP_2 EXP_1+EXP_1
#define EXP_4 EXP_2+EXP_2
#define EXP_8 EXP_4+EXP_4
#define EXP_16 EXP_8+EXP_8
#define EXP_32 EXP_16+EXP_16
#define EXP_64 EXP_32+EXP_32
#define EXP_128 EXP_64+EXP_64
#define EXP_256 EXP_128+EXP_128

#define STMT_1_1 fad+=EXP_1;
#define STMT_2_1 STMT_1_1;STMT_1_1;
#define STMT_4_1 STMT_2_1;STMT_2_1;
#define STMT_8_1 STMT_4_1;STMT_4_1;
#define STMT_16_1 STMT_8_1;STMT_8_1;
#define STMT_32_1 STMT_16_1;STMT_16_1;
#define STMT_64_1 STMT_32_1;STMT_32_1;
#define STMT_128_1 STMT_64_1;STMT_64_1;
#define STMT_256_1 STMT_128_1;STMT_128_1;

#define STMT_1_2 fad+=EXP_2;
#define STMT_2_2 STMT_1_2;STMT_1_2;
#define STMT_4_2 STMT_2_2;STMT_2_2;
#define STMT_8_2 STMT_4_2;STMT_4_2;
#define STMT_16_2 STMT_8_2;STMT_8_2;
#define STMT_32_2 STMT_16_2;STMT_16_2;
#define STMT_64_2 STMT_32_2;STMT_32_2;
#define STMT_128_2 STMT_64_2;STMT_64_2;
#define STMT_256_2 STMT_128_2;STMT_128_2;

#define STMT_1_4 fad+=EXP_4;
#define STMT_2_4 STMT_1_4;STMT_1_4;
#define STMT_4_4 STMT_2_4;STMT_2_4;
#define STMT_8_4 STMT_4_4;STMT_4_4;
#define STMT_16_4 STMT_8_4;STMT_8_4;
#define STMT_32_4 STMT_16_4;STMT_16_4;
#define STMT_64_4 STMT_32_4;STMT_32_4;
#define STMT_128_4 STMT_64_4;STMT_64_4;
#define STMT_256_4 STMT_128_4;STMT_128_4;

#define STMT_1_8 fad+=EXP_8;
#define STMT_2_8 STMT_1_8;STMT_1_8;
#define STMT_4_8 STMT_2_8;STMT_2_8;
#define STMT_8_8 STMT_4_8;STMT_4_8;
#define STMT_16_8 STMT_8_8;STMT_8_8;
#define STMT_32_8 STMT_16_8;STMT_16_8;
#define STMT_64_8 STMT_32_8;STMT_32_8;
#define STMT_128_8 STMT_64_8;STMT_64_8;
#define STMT_256_8 STMT_128_8;STMT_128_8;

#define STMT_1_16 fad+=EXP_16;
#define STMT_2_16 STMT_1_16;STMT_1_16;
#define STMT_4_16 STMT_2_16;STMT_2_16;
#define STMT_8_16 STMT_4_16;STMT_4_16;
#define STMT_16_16 STMT_8_16;STMT_8_16;
#define STMT_32_16 STMT_16_16;STMT_16_16;
#define STMT_64_16 STMT_32_16;STMT_32_16;
#define STMT_128_16 STMT_64_16;STMT_64_16;
#define STMT_256_16 STMT_128_16;STMT_128_16;


#define STMT_1_32 fad+=EXP_32;
#define STMT_2_32 STMT_1_32;STMT_1_32;
#define STMT_4_32 STMT_2_32;STMT_2_32;
#define STMT_8_32 STMT_4_32;STMT_4_32;
#define STMT_16_32 STMT_8_32;STMT_8_32;

#define STMT_1_64 fad+=EXP_64;
#define STMT_2_64 STMT_1_64;STMT_1_64;
#define STMT_4_64 STMT_2_64;STMT_2_64;

#define STMT_1_128 fad+=EXP_128;
#define STMT_2_128 STMT_1_128;STMT_1_128;

#define STMT_1_256 fad+=EXP_256;


// #define LOOP_BODY STMT_4_16 // the loop body

int get_num_ind(){
  return NUM_IND;
}

void get_initial_value(double *x){
  int i;
  for(i=0;i<NUM_IND;i++){
    x[i]=i+1;
  }
}
adouble func_eval(adouble * x){
  srand(12345);
  int  i, j, k;
  int r;
  adouble fad=0;  
  for(i = 0; i < NUM_IND; ++i) {
    for(j = 0; j < K_NNZ; j++) {
      r = rand() % NUM_IND;
      LOOP_BODY;
    }
  }
  return(fad);
}

