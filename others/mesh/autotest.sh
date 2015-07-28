#!/bin/sh
. ./../../adpath.sh

for meshdata in gear.mesh duct10.mesh duct8.mesh ; do
  echo "testing mesh : $meshdata"
  export LD_LIBRARY_PATH=$AD_LIBS_PATH:$OLD_LD_PATH
  for testmethod in LIVARH DIRECT INDIRECT; do
    command="$CXX -D $testmethod -O3 -Wno-unused-result -I$AD_INCL_PATH -I./include -o meshHess ./src/meshHess.cpp ./src/fcn3e_am.c ./src/elem3.c ./src/tet.c ./src/mesh3.c ./src/mesh.c ./src/opt3.c ./src/pre3.c -L$AD_LIBS_PATH -ladolc"
#    echo "$command";
    $command
    sh -c "./meshHess $meshdata"
  done

  testmethod="LIVARHACC"
  export LD_LIBRARY_PATH=$PREACC_LIBS_PATH:$OLD_LD_PATH
  command="$CXX -D $testmethod -O3 -Wno-unused-result -I$PREACC_INCL_PATH -I./include -o meshHess ./src/meshHess.cpp ./src/fcn3e_am.c ./src/elem3.c ./src/tet.c ./src/mesh3.c ./src/mesh.c ./src/opt3.c ./src/pre3.c -L$PREACC_LIBS_PATH -ladolc"
#  echo "$command";
  $command
  sh -c "./meshHess $meshdata"
  echo ""
done
export LD_LIBRARY_PATH=$OLD_LD_PATH
rm meshHess *.tap
