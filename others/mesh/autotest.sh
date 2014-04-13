#!/bin/sh
CXX="/usr/local/bin/g++"
CXXFLAGS="-O3"
AD_LIBS="-L${HOME}/adolc_edge/lib64 -ladolc"
AD_INCL="-I${HOME}/adolc_edge/include/ -I/usr/local/include/"

function trapFunc {
echo "$testfunc   Error";
}

trap 'trapFunc $testfunc ; exit' EXIT ERR
for testMethod in {0,1}; do
    for testIndex in {0,1}; do
        make > temp.out
        command="./meshHess 0 $testMethod 0 $testIndex";
#        echo "$command";
        $command > temp.out;
        sh -c "sort hess.out > hess.sorted";
        sh -c "sort edge.out > edge.sorted";
        [[ `diff hess.sorted edge.sorted` ]] &&  (echo "wrong!") || (echo "meshHess($testMethod, 0, $testIndex) OK!");
    done
done
rm meshHess *.out *.sorted *.tap
trap '' EXIT
