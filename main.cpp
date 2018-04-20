#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <vector>
#include "include/Parameters.h"
#include "include/Gutzwiller.h"
#include "include/File_op.h"

using namespace std;
using namespace file_op;

typedef complex<double> cd;
typedef vector<cd> cd1D;
typedef vector<cd1D> cd2D;
typedef vector<cd2D> cd3D;
typedef vector<cd3D> cd4D;

int main(int argc , char* argv[])
{
    printf("git test\n");
    //input: cluster_size, compnr, repnr, gentype (Gell-Mann matrices or Tensor Operators), tol (10e-tol), maxit
    Parameters params(argc, argv);
    Gutzwiller gutzwiller(params);
    gutzwiller.run();
//    print_params(params,gutzwiller.get_iter());
//    print_energy(params,gutzwiller.get_eigval(2));
//    print_n(params,gutzwiller.calc_n());
    print_phi(params,gutzwiller.get_phi_out());
//    printf("testtesttest\n");
//    print_spin(params,gutzwiller.get_spin_out());

    return 0;
}
