all:	meshHess

meshHess:	FORCE
	g++ -g -O3 -I/home/zsdfe/adolc_edge/include/ -I/usr/local/include/ -I./include -o meshHess ./src/meshHess.cpp ./src/fcn3e_am.c ./src/elem3.c ./src/tet.c ./src/mesh3.c ./src/mesh.c ./src/opt3.c ./src/pre3.c -L/home/zsdfe/adolc_edge/lib -ladolc
			
FORCE:
