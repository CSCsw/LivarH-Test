#!/bin/sh
. ./../../adpath.sh
TEST_CONF="-D NUM_IND=20000 -D K_NNZ=3"

echo "Testing: Performance dependence on number of basic operations per statement."
echo ""
for loop_body in STMT_128_1 STMT_64_2 STMT_32_4 STMT_16_8 STMT_8_16 STMT_4_32 STMT_2_64 STMT_1_128; do
  echo "configuration: $TEST_CONF $loop_body"
  export LD_LIBRARY_PATH=$AD_LIBS_PATH:$OLD_LD_PATH
  for testmethod in LIVARH DIRECT INDIRECT; do
    command="$CXX $TEST_CONF -D LOOP_BODY=$loop_body -D $testmethod $CXXFLAGS -I$AD_INCL_PATH -o hessTest hessTest.cpp f_hessian.cpp -L$AD_LIBS_PATH -ladolc"; 
    #echo "$command";
    $command
    sh -c "./hessTest";
  done
  testmethod="LIVARHACC"
  export LD_LIBRARY_PATH=$PREACC_LIBS_PATH:$OLD_LD_PATH
  command="$CXX $TEST_CONF -D LOOP_BODY=$loop_body -D $testmethod $CXXFLAGS -I$PREACC_INCL_PATH -o hessTest hessTest.cpp f_hessian.cpp -L$PREACC_LIBS_PATH -ladolc"; 
  #echo "$command";
  $command
  sh -c "./hessTest";
  echo ""
done
rm hessTest *.tap
