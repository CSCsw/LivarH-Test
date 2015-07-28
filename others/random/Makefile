CXX = g++
CXXFLAGS = -O3

ADOLCPATH=/home/zsdfe/adolc_edge
ADOLCLIBS=${ADOLCPATH}/lib
ADOLCINCL=$(ADOLCPATH)/include

all: hessTest

hessTest: hessTest.o f_hessian.o
	$(CXX) $(CXXFLAGS) hessTest.o f_hessian.o -L$(ADOLCLIBS) -ladolc -o hessTest

f_hessian.o: f_hessian.cpp
	$(CXX) $(CXXFLAGS) -I$(ADOLCINCL) -c -o f_hessian.o f_hessian.cpp

hessTest.o: hessTest.cpp	
	$(CXX) $(CXXFLAGS) -I$(ADOLCINCL) -c -o hessTest.o hessTest.cpp

clean:
	rm *.o *.tap hessTest

