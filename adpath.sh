#C++ Compiler
CXX="g++"

#Optimization level
CXXFLAGS="-O3 -w"

#Current LD_LIBRARY_PATH, don't need to change
OLD_LD_PATH="${LD_LIBRARY_PATH}"

#The path to the default adolc installation, with --enable-sparse
AD_LIBS_PATH="${HOME}/adolc_edge/lib64"
AD_INCL_PATH="${HOME}/adolc_edge/include/"

#The path to the adolc installation with --enable-preacc
PREACC_LIBS_PATH="${HOME}/adolc_preacc/lib64"
PREACC_INCL_PATH="${HOME}/adolc_preacc/include/"
