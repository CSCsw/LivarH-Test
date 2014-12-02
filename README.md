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

hessTest.cpp & functions/test*.cpp : simple & synthetic test functions 

test1-test6: simple functionality tests, some of them may require flags "--enable-atrig-erf --enable-advanced-branching" when configure ADOL-C.

testF1-testF4: the four tests used in performance measurement.

Please refer to autotest.sh to see how to run all the tests automatically.

Please refer to hessTest.cpp to see how to choose the algorithms.


others/mesh/ : mesh optimization problem
A Makefile is provided, user need to change the ADOLC path in that file.
Then compile with "make" command. To run the test, use
"./meshHess $testdata $testMethod $option0 $option1"
For example:
"./meshHess gear.mesh 0 0 1"
$testdata: the mesh data to be use, we reported "gear.mseh", "duct12.mesh" and "duct8.mesh".
$testMethod: 0 means edge_hess(), 1 means ADOL-C sparse_hess()
$option0/$option1: options[2] provided with the selected methods.
For sparse_hess() :
options[0] = 0 (safe mode, default) / 1 (tight mode)
options[1] = 0 (indirect recovery, default) /1 (direct recovery)
For edge_hess():
options[0] = 0 (LivarH) / 1 (LivarHAcc)
options[1] == 1 (fixed as 1, not used anymore)


others/random/ : random hessian generator
