#ifndef FILE_OP_H
#define FILE_OP_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <complex.h>
#include <vector>
#include "Eigstatval.h"
#include "Parameters.h"

using namespace std;

typedef complex<double> cd;
typedef vector<cd> cd1D;
typedef vector<cd1D > cd2D;
typedef vector<cd2D> cd3D;
typedef vector<cd3D> cd4D;

namespace file_op
{
    void print_params(Parameters params, int iter);
    void print_phi(Parameters params, cd4D phiout);
    void print_spin(Parameters params, cd4D spinout);
    void print_n(Parameters params, cd3D nout);
    void print_energy(Parameters params, cd eigval);
}

#endif // FILE_OP_H
