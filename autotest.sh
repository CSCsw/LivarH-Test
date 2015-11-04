#!/bin/sh
. ./adpath.sh
for testfunc in ./functions/testF1.cpp ./functions/testF2.cpp ./functions/testF3.cpp ./functions/testF4.cpp ; do
  echo "testing: $testfunc"
  export LD_LIBRARY_PATH=$AD_LIBS_PATH:$OLD_LD_PATH
  for testmethod in COMPUTE_FULL_HESS DIRECT INDIRECT LIVARH; do
    command="$CXX -D $testmethod $CXXFLAGS -I$AD_INCL_PATH -o hessTest hessTest.cpp $testfunc -L$AD_LIBS_PATH -ladolc"; 
    #echo "$command";
    $command
    sh -c "./hessTest";
  done
  testmethod="LIVARHACC"
  export LD_LIBRARY_PATH=$PREACC_LIBS_PATH:$OLD_LD_PATH
  command="$CXX -D $testmethod $CXXFLAGS -I$PREACC_INCL_PATH -o hessTest hessTest.cpp $testfunc -L$PREACC_LIBS_PATH -ladolc"; 
  #echo "$command";
  $command
  sh -c "./hessTest";
  echo ""
done
export LD_LIBRARY_PATH=$OLD_LD_PATH
rm hessTest *.tap
