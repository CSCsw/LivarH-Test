all:	meshHess

meshHess:	FORCE
			/usr/local/bin/g++ -O3 -I/Users/wang970/adolc_edge/include/ -I/usr/local/include/ -I./include -o meshHess ./src/meshHess.cpp ./src/fcn3e_am.c ./src/elem3.c ./src/tet.c ./src/mesh3.c ./src/mesh.c ./src/opt3.c ./src/pre3.c -L/Users/wang970/adolc_edge/lib64 -ladolc
			
FORCE:
