#!/bin/sh
. ./../../adpath.sh

LOOP_BODY="STMT_2_8"

echo "Testing: Performance dependence on sparsity of the final Hessian."
echo ""
for test_conf in "-D NUM_IND=30000 -D K_NNZ=2" "-D NUM_IND=20000 -D K_NNZ=3" "-D NUM_IND=15000 -D K_NNZ=4" "-D NUM_IND=10000 -D K_NNZ=6" ; do
  echo "configuration: $test_conf $LOOP_BODY"
  export LD_LIBRARY_PATH=$AD_LIBS_PATH:$OLD_LD_PATH
  for testmethod in LIVARH DIRECT INDIRECT; do
    command="$CXX $test_conf -D LOOP_BODY=$LOOP_BODY -D $testmethod $CXXFLAGS -I$AD_INCL_PATH -o hessTest hessTest.cpp f_hessian.cpp -L$AD_LIBS_PATH -ladolc"; 
    #echo "$command";
    $command
    sh -c "./hessTest";
  done
  testmethod="LIVARHACC"
  export LD_LIBRARY_PATH=$PREACC_LIBS_PATH:$OLD_LD_PATH
  command="$CXX $test_conf -D LOOP_BODY=$LOOP_BODY -D $testmethod $CXXFLAGS -I$PREACC_INCL_PATH -o hessTest hessTest.cpp f_hessian.cpp -L$PREACC_LIBS_PATH -ladolc"; 
  #echo "$command";
  $command
  sh -c "./hessTest";
  echo ""
done
rm hessTest *.tap
