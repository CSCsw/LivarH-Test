CXX = g++
CXXFLAGS = -O3

ADOLCPATH=/home/zsdfe/adolc_edge
ADOLCLIBS=${ADOLCPATH}/lib
ADOLCINCL=$(ADOLCPATH)/include

all: col_hessian

col_hessian: col_hessian.o f_hessian.o
	$(CXX) $(CXXFLAGS) col_hessian.o f_hessian.o -L$(ADOLCLIBS) -ladolc -o col_hessian

f_hessian.o: f_hessian.cpp
	$(CXX) $(CXXFLAGS) -I$(ADOLCINCL) -c -o f_hessian.o f_hessian.cpp

col_hessian.o: col_hessian.cpp	
	$(CXX) $(CXXFLAGS) -I$(ADOLCINCL) -c -o col_hessian.o col_hessian.cpp

clean:
	rm *.o col_hessian

