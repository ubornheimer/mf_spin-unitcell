#include "../include/Gutzwiller.h"

Gutzwiller::Gutzwiller(Parameters params)
{
    gparams = params;
    phi_out = gparams.get_phi_in();
    spin_out = gparams.get_phi_in();
}

Gutzwiller::~Gutzwiller()
{}

void Gutzwiller::run()
{
    int maxit = gparams.get_maxit();
    int minit = gparams.get_minit();
    int compnr = gparams.get_compnr();
    int lsize = gparams.get_sitenr();
    int unitcnr = gparams.get_unitcellnr();
    int iter = 0;

    cd4D phi_in = phi_out;
    cd4D phi_temp = phi_out;

    do
    {
        notconv = true;
        for (int isite1 = 0; isite1 < lsize; isite1++)
        {
            for (int isite2 = 0; isite2 < lsize; isite2++)
            {
                for (int iunit = 0; iunit < unitcnr; iunit++)
                {
                    //printf("site %i %i %i\n",isite1,isite2,iunit);
                    //for (int gen = 0; gen < 9; gen++)
                    //    printf("phi_in %f %f\n", real(phi_in[isite1][isite2][iunit][gen]), imag(phi_in[isite1][isite2][iunit][gen]));
                    gdiag = diagonalize(calc_hamiltonian(phi_in[isite1][isite2][iunit]),compnr);
                    calc_phi(isite1, isite2, iunit);
                    //printf("diagonalized and calc_phi\n");
                }
            }
        }

        printf("diagonalization and calc_phi of all sites completed\n");

        notconv = false;
        dist = 0;

        for (int isite1 = 0; isite1 < lsize; isite1++)
        {
            for (int isite2 = 0; isite2 < lsize; isite2++)
            {
                for (int iunit = 0; iunit < unitcnr; iunit++)
                {
                    //printf("site %i %i %i\n",isite1,isite2,iunit);
                    phi_temp[isite1][isite2][iunit]=calc_neighbours(isite1,isite2,iunit);
                    //printf("calc neighbours successfull\n");
                    if (notconv == false)
                    {
                        printf("site %i %i %i\n",isite1,isite2,iunit);
                        notconv = comp_phi(phi_in[isite1][isite2][iunit],phi_temp[isite1][isite2][iunit]);
                        //calc_spin(isite1,isite2,iunit);
                    }
                }
            }
        }

        printf("distance %f\n", real(dist));

        phi_out = phi_in+ u*(phi_temp-phi_in); //if u=0: nothing happens, if u=1: new mf-phi is the new input
        phi_in = phi_out;
        iter++;
        if (notconv == false) {printf("converged!\n"); }
        printf("iteration %i\n", iter);
        printf(".................\n");
    }
    while((notconv or (iter <= minit)) and (iter <= maxit));
//
//    nriter = iter;
//
//    printf("\n");
//
//    if(iter > maxit)
//        printf("did not converge! maxit reached!\n");
//    else
//        printf("nr of iterations: %i\n", iter);
}

cd2D Gutzwiller::calc_hamiltonian(cd1D phitemp)
{
    int repnr = gparams.get_repnr();
    int gentype = gparams.get_gentype();
    int compnr = gparams.get_compnr();
    cd3D gencoeff = gparams.get_coeffgen();
    cd3D gens = gparams.get_gens(gentype);
    //single-site mean-field
    //two-site mean-field
    //cd2D ham = cron(zero_matrix(compnr),zero_matrix(compnr));
    cd2D uu = unit_matrix(compnr);
    cd2D ham = uu;

//    for (int i=0; i<repnr; i++)
//    {

        for (int j=0; j<repnr; j++)
        {
            //printf("calc_ham %f\n",phitemp[j]);
            ham += 2.*phitemp[j]*gens[j]-phitemp[j]*phitemp[j]*uu;//?
        }
//    }

//    for (int i=0; i<3; i++)
//    {
//        for(int j=0; j<3; j++)
//            printf("%f %f\t",real(ham[i][j]), imag(ham[i][j]));
//        printf("\n");
//    }
//
    return ham;
}

void Gutzwiller::calc_phi(int isite1, int isite2, int iunit)
{
    int repnr = gparams.get_repnr();
    int compnr = gparams.get_compnr();
    int gentype = gparams.get_gentype();
    cd3D gens = gparams.get_gens(gentype);
    cd1D expvals;
    cd3D gencoeff = gparams.get_coeffgen();
    cd3D crontemp;

    //printf("calc_phi test");

    //calculate the expectation values of the generators
    if (notconv==true)
    {
        cd1D state = gdiag.get_eigstat(0);

//        for (int m=0; m<compnr*compnr; m++)
//            printf("eigstat %f\n",state[m]);

        for (int n=0; n<repnr; n++)
        {
            crontemp = vecunit(1,compnr);
            crontemp[0] = gens[n];
            cd test = conj(state)*cron(crontemp)*state;
            expvals.push_back(test);
            //printf("calc_phi %f %f\n", real(test), imag(test));
//            printf("calc_phi %i %i %i %f\n", isite1,isite2,iunit,real(expvals[n]));
        }

        //printf("phi %i %i %i %f\n", isite1,isite2,iunit,real(expvals[1]));
    }
    else
    {
        cd4D phi_in = gparams.get_phi_in();
        expvals = phi_in[isite1][isite2][iunit];
        printf("converged phi\n");
    }

    phi_out[isite1][isite2][iunit] = expvals;
    //printf("calc_phi %i %i %i\n", isite1, isite2, iunit);

    //calculate the mean-field
//    for (int i=0; i<repnr; i++)
//    {
//        phimftemp=0;
//        for (int j=0; j<repnr; j++)
//        {
//            for (int nn=0; nn<linknr; nn++)
//            {isite1
//                phimftemp+=expvals[j]*gencoeff[nn][i][j];
//            }
//        }
////        printf("phimftemp %f\n",real(phimftemp));
//        phi_out[0][i]=phimftemp;
//    }
}

void Gutzwiller::calc_spin(int isite1, int isite2, int iunit)
{
    int compnr = gparams.get_compnr();
    cd3D spin = gparams.get_spin();
    cd1D expvals;
    cd3D crontemp;

    //calculate the expectation values of the generators
    cd1D state = gdiag.get_eigstat(0);

    for (int n=0; n<compnr; n++)
    {
        crontemp = vecunit(1,compnr);
        crontemp[0] = spin[n];
        cd test = conj(state)*cron(crontemp)*state;
        expvals.push_back(test);
        //printf("phi[%i][%i]: %f, %f\n", i,n,real(test), imag(test));
    }


    spin_out[isite1][isite2][iunit] = expvals;

    //calculate the mean-field
//    for (int i=0; i<repnr; i++)
//    {
//        phimftemp=0;
//        for (int j=0; j<repnr; j++)
//        {
//            for (int nn=0; nn<linknr; nn++)
//            {
//                phimftemp+=expvals[j]*gencoeff[nn][i][j];
//            }
//        }
////        printf("phimftemp %f\n",real(phimftemp));
//        phi_out[0][i]=phimftemp;
//    }
}

cd1D Gutzwiller::calc_neighbours(int isite1, int isite2, int iunit)
{
    cd1D phitemp;
    int repnr = gparams.get_repnr();
    int sitenr = gparams.get_sitenr();
    cd3D gencoeff = gparams.get_coeffgen();

    int unext;

    if (iunit == 0)
        unext = -1;
    else
        unext = 1;

//    printf("neighbours %i %i %i\n", isite1, isite2, (iunit+1)%2);
//    printf("neighbours %i %i %i\n", (isite1+unext+sitenr)%sitenr, isite2, (iunit+1)%2);
//    printf("neighbours %i %i %i\n", isite1, (isite2+unext+sitenr)%sitenr, (iunit+1)%2);

    //calculate the mean-field
    for (int i=0; i<repnr; i++)
    {
        cd phimftemp = 0;
        for (int j=0; j<repnr; j++)
        {
            phimftemp+=phi_out[isite1][isite2][(iunit+1)%2][j]*gencoeff[0][i][j];
            //printf("neigh %i %i %i %i %f %f\n", isite1, isite2, (iunit+1)%2, j, abs(phi_out[isite1][isite2][(iunit+1)%2][j]), gencoeff[0][i][j]);
            phimftemp+=phi_out[(isite1+unext+sitenr)%sitenr][isite2][(iunit+1)%2][j]*gencoeff[1][i][j];
            //printf("neigh %i %i %i %i %f %f\n", (isite1+unext+sitenr)%sitenr, isite2, (iunit+1)%2, j, abs(phi_out[(isite1+unext+sitenr)%sitenr][isite2][(iunit+1)%2][j]), gencoeff[1][i][j]);
            phimftemp+=phi_out[isite1][(isite2+unext+sitenr)%sitenr][(iunit+1)%2][j]*gencoeff[2][i][j];
            //printf("neigh %i %i %i %i %f %f\n", isite1, (isite2+unext+sitenr)%sitenr, (iunit+1)%2, j, abs(phi_out[isite1][(isite2+unext+sitenr)%sitenr][(iunit+1)%2][j]), gencoeff[2][i][j]);

        }
        //printf("calc_neighbours %i %i %i %i %f %f\n",isite1, isite2, iunit, i, real(1./6.*phimftemp), imag(1./6.*phimftemp));
        phitemp.push_back(1./6.*phimftemp);
    }

    //printf("phimf %i %i %i %f\n", isite1, isite2, iunit, phitemp[1]);

    return phitemp;
}

bool Gutzwiller::comp_phi(cd1D phi_old, cd1D phi_new)
{
    bool retnotconv = false;
    int repnr = gparams.get_repnr();
    int tol = gparams.get_tol();

    double ddist=0;

    for(int j=0; j<repnr; j++)
    {
        printf("phi_in %f, phi_out %f\n",real(phi_old[j]),real(phi_new[j]));
        ddist = abs(phi_old[j]-phi_new[j]);
        if(ddist>pow(10,-tol))
            retnotconv = true;
        if(ddist>dist)
            dist = ddist;
    }
    //printf(".................\n");
    return retnotconv;
}

cd3D Gutzwiller::calc_n()
{
    cd3D result;
//    cd2D temp;
//    cd1D ttemp;
//    int lsize = gparams.get_lattice_size();
//    int focknr = gparams.get_focknr();
//    int compnr = gparams.get_compnr();
//    cd1D eigstat;
//    cd3D crontemp;
//
//    for(int xsite = 0 ; xsite < lsize; xsite++)
//    {
//        result.push_back(temp);
//        for(int ysite = 0; ysite < lsize; ysite++)
//        {
//            result[xsite].push_back(ttemp);
//            for(int comp=0; comp< compnr; comp++)
//            {
//                eigstat = gdiag[xsite][ysite].get_eigstat(0);
//                crontemp = vecunit(compnr,focknr);
//                crontemp[comp] = n_matrix(focknr);
//                result[xsite][ysite].push_back(conj(eigstat)*cron(crontemp)*eigstat);
//            }
//        }
//    }

    return result;
}

cd Gutzwiller::calc_energy(Parameters params)
{
    cd eigval;
    //cd eigval = gdiag.get_eigval(0);
//    int lsize = params.get_cluster_size();
//    cd2D J = params.get_J();
//    int compnr = params.get_compnr();
//    cd3D phi = params.get_phi_in();
//
////    for(int i=0; i<lsize; i++)
////    {
////        for(int j=0; j<lsize; j++)
////        {
//    int i=0;
//    int j=0;
//    for(int k=0; k<compnr; k++)
//    {
//        cd nn1 = conj(phi[i][j][k])*phi[i][(j+1)%lsize][k];
//        cd nn2 = conj(phi[i][j][k])*phi[i][(j+lsize-1)%lsize][k];
//        cd nn3 = conj(phi[i][j][k])*phi[(i+1)%lsize][j][k];
//        cd nn4 = conj(phi[i][j][k])*phi[(i+lsize-1)%lsize][j][k];
//        eigval += J[k][k]/8.*(nn1+conj(nn1)+nn2+conj(nn2)+nn3+conj(nn3)+nn4+conj(nn4));
//    }
////        }
////    }

    return eigval;
}
