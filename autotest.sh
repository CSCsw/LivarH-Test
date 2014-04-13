#!/bin/sh
CXX="/usr/local/bin/g++"
CXXFLAGS="-O3"
AD_LIBS="-L${HOME}/adolc_edge/lib64 -ladolc"
AD_INCL="-I${HOME}/adolc_edge/include/ -I/usr/local/include/"

function trapFunc {
  echo "$testfunc   Error";
}

trap 'trapFunc $testfunc ; exit' INT TERM EXIT ERR
for testfunc in ./functions/* ; do
  command="$CXX $CXXFLAGS $AD_INCL -o hessTest hessTest.cpp $testfunc $AD_LIBS"; 
#  echo "$command";
  $command
  sh -c "./hessTest > temp.out";
  echo "${testfunc} OK!";
done
rm hessTest temp.out

trap '' EXIT
