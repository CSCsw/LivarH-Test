#!/bin/sh
CXX="/usr/local/bin/g++"
CXXFLAGS="-O3"
AD_LIBS="-L${HOME}/packages/adolc_edge/lib64 -ladolc"
AD_INCL="-I${HOME}/packages/adolc_edge/include/ -I/usr/local/include/"

for testfunc in ./functions/* ; do
  command="$CXX $CXXFLAGS $AD_INCL -o hessTest hessTest.cpp $testfunc $AD_LIBS"; 
  echo "$command";
  $command
  sh -c "./hessTest > temp.out";
  result=$(grep -i "CORRECT!" temp.out);
  echo "$result";
  if [ "$result" = "CORRECT!" ]
    then echo "$testfunc           CORRECT!";
    else echo "$testfunc             WRONG!";
  fi 
done
rm hessTest temp.out
