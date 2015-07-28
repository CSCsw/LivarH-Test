A collection of simple tests to check the consistency of the hessian computed by:

hessian() in ADOL-C and edge_hess() (which implements the edge_pushing algorithm)

edge_hess() has almost the same parameters as sparse_hess()
Except edge_hess() has the dependents check but does not have a repeat flag, and a different meaning for options[2]:

    options[0] = 0 (Default)     \\ Disable Preaccumulation
               = 1               \\ Enable Preaccumulation 
                                 \\ (need --enable-preacc when configure ADOL-C)

    options[0] = 1               \\ Deprecated

Include the header "adolc/hessian/edge_main.h"
And then call edge_hess() in the same way as if sparse_hess() is called.
## Configuration : adpath.sh
Change the content in adpath.sh according to your installation. The default values are:

//C++ compiler
CXX="g++" 

//Optimization level
CXXFLAGS="-O3"

//Current LD_LIBRARY_PATH, don't change
OLD_LD_PATH="${LD_LIBRARY_PATH}"

//The path to the default adolc installation with --enable-sparse. (Change this according to your installation)
AD_LIBS_PATH="${HOME}/adolc_edge/lib64"
AD_INCL_PATH="${HOME}/adolc_edge/include"

//The path to the adolc installation with preaccumulation, --enable-preacc. (Change this according to your installation)
PREACC_LIBS_PATH="${HOME}/adolc_preacc/lib64"
PREACC_INCL_PATH="${HOME}/adolc_preacc/include"

The setting is shared by all following scripts.

## hessTest.cpp & functions/test*.cpp : simple & synthetic test functions 

test1-test6: simple functionality tests, some of them may require flags "--enable-atrig-erf --enable-advanced-branching" when configure ADOL-C.

testF1-testF4: the four tests used in performance measurement.

The autotest.sh script is provided to automatically run the testF1-testF4 using the setting in the paper. (Table 6 in the paper). Simply type:
./autotest.sh

Note: only overall runtime is reported, for breakdown performance for Direct/Indirect methods in ADOL-C, please instrument sparse_hess() in "ADOL-C/src/sparse/sparsedriver.cpp".

##others/mesh/ : mesh optimization problem
The three dataset used in the paper are : gear.mesh, duct12.mesh, duct8.mesh.Also a "autotest.sh" is provided to regenerate results on these cases. (Table 11 in the paper)

Note: the script only run the sparse methods in this case. Because the full hessian takes extrodinally long time to finish. 

##others/random/ : random hessian generator
autotestA.sh : Test the performance dependency on number of basic operations per statement. (Figure 10 in the paper)

autotestB.sh : Test the scalability with respect to function complexity. (Figure 12 in the paper)

autotestC.sh : Test the performance dependency on sparsity of the final Hessian. (Figure 13 in the paper)
