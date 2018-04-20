#ifndef GUTZWILLER_H
#define GUTZWILLER_H

#include <complex.h>
#include <vector>
#include "Parameters.h"
#include <math.h>
#include "../Utils.h"
#include "Eigstatval.h"

using namespace std;
using namespace utils;

typedef complex<double> cd;
typedef vector<cd> cd1D;
typedef vector<cd1D > cd2D;
typedef vector<cd2D> cd3D;
typedef vector<cd3D> cd4D;
typedef vector<int> int1D;

const cd u = 0.5;

class Gutzwiller
{
    public:
        Gutzwiller(Parameters params);
        virtual ~Gutzwiller();
        void run();
        void diag();
        cd3D calc_n();
        cd4D get_phi_out() { return phi_out; }
        cd4D get_spin_out() { return spin_out;}
        int get_iter() { return nriter; }
        cd get_eigval(int i){return gdiag.get_eigval(i);}
        cd calc_energy(Parameters params);
    protected:
    private:
        Parameters gparams;
        Eigstatval gdiag;
        cd4D phi_out;
        cd4D spin_out;
        int nriter;
        cd2D calc_hamiltonian(cd1D phitemp);
        cd1D calc_neighbours(int isite1, int isite2, int iunit);
        void calc_phi(int isite1, int isite2, int iunit);
        void calc_spin(int isite1, int isite2, int iunit);
        bool comp_phi(cd1D phi_in, cd1D phi_out);
        bool notconv = false;
        double dist=0;
};

#endif // GUTZWILLER_H
