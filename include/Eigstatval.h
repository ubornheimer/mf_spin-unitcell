#ifndef EIGSTATVAL_H
#define EIGSTATVAL_H

#include <complex.h>
#include <vector>

using namespace std;

typedef complex<double> cd;
typedef vector<cd> cd1D;
typedef vector<cd1D > cd2D;

class Eigstatval
{
    public:
        Eigstatval();
        Eigstatval(cd1D eigval_in, cd2D eigstat_in);
        virtual ~Eigstatval();
        cd get_eigval(int i){return eigval[i];}
        cd2D get_eigstat(){return eigstat;}
        cd1D get_eigstat(int n){return eigstat[n];}
        void set_eigval(cd1D in){eigval = in;}
        void set_eigstat(cd2D in){eigstat = in;}
    protected:
    private:
        cd1D eigval;
        cd2D eigstat;
};

#endif // EIGSTATVAL_H
