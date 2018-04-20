#include "../include/Eigstatval.h"

Eigstatval::Eigstatval()
{
}

Eigstatval::Eigstatval(cd1D eigval_in, cd2D eigstat_in)
{
    eigval = eigval_in;
    eigstat = eigstat_in;
}

Eigstatval::~Eigstatval()
{
    //dtor
}
