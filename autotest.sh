#!/bin/sh
CXX="g++"
CXXFLAGS="-O3"
ADOLCLIBS="-L${HOME}/adolc_edge/lib -ladolc"
ADOLCINCL="-I${HOME}/adolc_edge/include/ -I/usr/local/include/"

for testfunc in ./functions/* ; do
  command="$CXX $CXXFLAGS $ADOLCINCL -o hessTest hessTest.cpp $testfunc $ADOLCLIBS"; 
  echo "$command";
  $command
  sh -c "./hessTest";
done
rm hessTest
