# This file is distributed under the Common Public License.

# CHANGEME: Additional flags for compilation (e.g., include flags)

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile.   #
##########################################################################
#g++ -I/home/eldronin/adolc_edge/include/adolc -L/home/eldronin/adolc_edge/lib 
#-Wl,--rpath -Wl,/home/eldronin/adolc_edge/lib -ladolc -lColPack *.cpp  -o p
# C++ Compiler command
CXX = g++

# C++ Compiler options
#CXXFLAGS = -g -O0 -v -da -Q 
CXXFLAGS = -g

# Libraries necessary to link with 
LIBS_col = -L${HOME}/adolc_edge/lib64 -ladolc

# -lColPack

# Necessary Include dirs (we use the CYGPATH_W variables to allow
# compilation with Windows compilers)
INCL_col =  -I${HOME}/adolc_edge/include/ $(ADDINCFLAGS)

all: col_hessian

col_hessian: col_hessian.o f_hessian.o
	$(CXX) $(CXXFLAGS) col_hessian.o f_hessian.o $(LIBS_col) -o col_hessian

f_hessian.o: f_hessian.cpp
	$(CXX) $(CXXFLAGS) $(INCL_col) -c -o f_hessian.o f_hessian.cpp

col_hessian.o: col_hessian.cpp	
	$(CXX) $(CXXFLAGS) $(INCL_col) -c -o col_hessian.o col_hessian.cpp

clean:
	rm *.o col_hessian

