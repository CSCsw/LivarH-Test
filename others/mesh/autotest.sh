#!/bin/sh
CXX="g++"
CXXFLAGS="-O3"
AD_LIBS="-L${HOME}/adolc_edge/lib -ladolc"
AD_INCL="-I${HOME}/adolc_edge/include/ -I/usr/local/include/"
meshdata="gear.mesh"

for testMethod in 0 1
do
  for option0 in 0 1
  do
    for option1 in 0 1
    do
        make > temp.out
        command="./meshHess $meshdata $testMethod $option0 $option1";
        echo "$command";
        $command > output.txt;
        sort hess.out > hess.sorted;
        sort edge.out > edge.sorted;
        cmp -s hess.sorted edge.sorted > /dev/null;
        if [ $? -eq 1 ]; then
          echo "WRONG!"
        else
          echo "messHess($meshdata, $testMethod, $option0, $option1) CORRECT"
        fi
    done
  done
done
rm meshHess *.out *.sorted *.tap
