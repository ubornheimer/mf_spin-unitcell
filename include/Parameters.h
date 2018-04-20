#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <complex.h>
#include <vector>
#include <random>

using namespace std;

typedef complex<double> cd;
typedef vector<cd> cd1D;
typedef vector<cd1D > cd2D;
typedef vector<cd2D> cd3D;
typedef vector<cd3D> cd4D;

class Parameters
{
    public:
        Parameters();
        Parameters(int argc, char* argv[]);
        virtual ~Parameters();
//        int get_cluster_size() { return cluster_size; }
        int get_repnr() { return repnr; }
        int get_compnr() { return compnr; }
        int get_tol(){return tol;}
        int get_maxit(){return maxit;}
        int get_minit(){return minit;}
        int get_gentype(){return gentype;}
        int get_linknr(){return linknr; }
        int get_unitcellnr(){return unitcellnr;}
        int get_sitenr(){return sitenr;}
        cd4D get_phi_in() { return phi_in; }
        cd3D get_coeffgen() { return coeffgen; }
        cd3D get_gens(int tt) {if(tt==0) return GellMann; else if(tt==1) return TensorOps; else printf("error returning generators!");}
        cd3D get_spin() { return SpinOps;}
    protected:
    private:
        int unitcellnr; // size of the cluster
        int compnr;
        int linknr; // number of different links to the cluster
        int repnr; // size of representation
        int gentype;
        int tol;
        int maxit;
        int minit;
        int sitenr;
        double a0;
        double a2;
        cd cI = cd(0.0,1.0);
        cd4D phi_in;
        cd3D coeffgen;
        cd4D read_phiin();
        cd3D read_coeffgen();
        cd2D ffile_coeffgen(string filename);
        cd3D GellMann = {{{1./sqrt(3.),0.,0.},{0.,1./sqrt(3.),0.},{0.,0.,1./sqrt(3.)}},
                         {{0.,1./sqrt(2.),0.},{1./sqrt(2.),0.,0.},{0.,0.,0.}},
                         {{0.,-cI/sqrt(2.),0.},{cI/sqrt(2.),0.,0.},{0.,0.,0.}},
                         {{1./sqrt(2.),0.,0.},{0.,-1./sqrt(2.),0.},{0.,0.,0.}},
                         {{0.,0.,1./sqrt(2.)},{0.,0.,0.},{1./sqrt(2.),0.,0.}},
                         {{0.,0.,-cI/sqrt(2.)},{0.,0.,0.},{cI/sqrt(2.),0.,0.}},
                         {{0.,0.,0.},{0.,0.,1./sqrt(2.)},{0.,1./sqrt(2.),0.}},
                         {{0.,0.,0.},{0.,0.,-cI/sqrt(2.)},{0.,cI/sqrt(2.),0.}},
                         {{1./sqrt(6.),0.,0.},{0.,1./sqrt(6.),0.},{0.,0.,-2./sqrt(6.)}}};
        cd3D TensorOps = {{{1./sqrt(3.),0.,0.},{0.,1./sqrt(3.),0.},{0.,0.,1./sqrt(3.)}},
                         {{0.,1./sqrt(2.),0.},{0.,0.,1./sqrt(2.)},{0.,0.,0.}},
                         {{1./sqrt(2.),0.,0.},{0.,0.,0.},{0.,0.,-1./sqrt(2.)}},
                         {{0.,0.,0.},{1./sqrt(2.),0.,0.},{0.,1./sqrt(2.),0.}},
                         {{0.,0.,1.},{0.,0.,0.},{0.,0.,0.}},
                         {{0.,-1./sqrt(2.),0.},{0.,0.,1./sqrt(2.)},{0.,0.,0.}},
                         {{1./sqrt(6.),0.,0.},{0.,-2./sqrt(6.),0.},{0.,0.,1./sqrt(6.)}},
                         {{0.,0.,0.},{1./sqrt(2.),0.,0.},{0.,-1./sqrt(2.),0.}},
                         {{0.,0.,0.},{0.,0.,0.},{1.,0.,0.}}};
        cd3D SpinOps = {{{0.,1./sqrt(2.),0.},{1./sqrt(2.),0.,1./sqrt(2.)},{0.,1./sqrt(2.),0.}},
                        {{0., cI/sqrt(2.), 0.},{cI/sqrt(2.),0.,-cI/sqrt(2.)},{0.,-cI/sqrt(2.),0.}},
                        {{1.,0.,0.},{0.,0.,0.},{0.,0.,-1.}}};
};

#endif // PARAMETERS_H
