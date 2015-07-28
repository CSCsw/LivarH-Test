CXX=g++

#PATH TO YOUR ADOL-C
ADOLCPATH=/home/muwang/adolc_edge
ADOLCLIBS=$(ADOLCPATH)/lib64

all:	meshHess

meshHess:	FORCE
	$(CXX) -O3 -Wno-unused-result -I$(ADOLCPATH)/include/ -I/usr/local/include/ -I./include -o meshHess ./src/meshHess.cpp ./src/fcn3e_am.c ./src/elem3.c ./src/tet.c ./src/mesh3.c ./src/mesh.c ./src/opt3.c ./src/pre3.c -L$(ADOLCLIBS) -ladolc

clean: FORCE
	rm meshHess *.tap *.out *.sorted			
FORCE:
