#include "../include/File_op.h"

namespace file_op
{
    void print_params(Parameters params, int iter)
    {
        FILE * paramfile;
        paramfile = fopen("params.txt","a");

        fprintf(paramfile, "run of meanfield_spinunitcell\n");
        fprintf(paramfile, "a cluster mean-field treatment of a SW Hamiltonian in the HC lattice\n");
//        fprintf(paramfile, "cluster of total size: %i\n", params.get_cluster_size());
        fprintf(paramfile, "number of components: %i\n", params.get_compnr());
        fprintf(paramfile, "number of elements in representation: %i\n", params.get_repnr());
        fprintf(paramfile, "type of representation: %i\n", params.get_gentype());
        fprintf(paramfile, "tolerance for conversion 10e-tol: %i\n", params.get_tol());
        fprintf(paramfile, "maximum iterations allowed: %i\n",params.get_maxit());
        fprintf(paramfile, "iterations used: %i\n", iter);

        fprintf(paramfile, "------------------------------------------------------\n");
        fclose(paramfile);
    }

    void print_phi(Parameters params, cd4D phi_out)
    {
        FILE * existsfparfile = fopen("sfpar.txt", "r");
        FILE * sfparfile = fopen("sfpar.txt", "a");
        FILE * tempsfparfile;
        int cellnr = params.get_unitcellnr();
        int sitenr = params.get_sitenr();
        int totalsites = cellnr*sitenr;

//        int csize = params.get_cluster_size();
        int repnr = params.get_repnr();

        if(existsfparfile == NULL)
        {
            fprintf(sfparfile, "Final order parameter for each site in the cluster and each element of the representation\n");
            for (int i=0; i<cellnr; i++)
            {
                for (int j=0; j>sitenr; j++)
                {
                    fprintf(sfparfile, "site %i, %i \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t",j,i);
                }
            }
            fprintf(sfparfile,"\n");
            for (int i=0; i<totalsites; i++)
            {
                for (int j=0; j<repnr; j++)
                {
                    fprintf(sfparfile, "el %i\t\t",j);
                }
            }
            fprintf(sfparfile,"\n");
        }
        else
            fclose(existsfparfile);

        tempsfparfile = fopen("tempphi.txt", "w");

        for (int i=0; i<sitenr; i++)
        {
            for (int l=0; l<sitenr; l++)
            {
                for (int k=0; k<cellnr; k++)
                {
                    for (int j=0; j<repnr; j++)
                    {
                        fprintf(sfparfile, "%.6f\t", real(phi_out[i][l][k][j]));
                        fprintf(tempsfparfile, "%.10f\t", real(phi_out[i][l][k][j]));
                    }
                }
            }
        }

        fprintf(tempsfparfile, "\n");
        fprintf(sfparfile, "\n");
        fclose(sfparfile);
        fclose(tempsfparfile);
    }

    void print_spin(Parameters params, cd4D spin_out)
    {
        FILE * existspinfile = fopen("spin.txt", "r");
        FILE * spinfile = fopen("spin.txt", "a");
        int cellnr = params.get_unitcellnr();
        int sitenr = params.get_sitenr();
        int totalsites = cellnr*sitenr;
        int compnr = params.get_compnr();

        if(existspinfile == NULL)
        {
            fprintf(spinfile, "spin expectation values for each site in the lattice and each direction\n");
            for (int i=0; i<cellnr; i++)
            {
                for (int j=0; j>sitenr; j++)
                {
                    fprintf(spinfile, "site %i, %i \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t",j,i);
                }
            }
            fprintf(spinfile,"\n");
            for (int i=0; i<totalsites; i++)
            {
                for (int j=0; j<compnr; j++)
                {
                    fprintf(spinfile, "dir %i\t\t",j);
                }
            }
            fprintf(spinfile,"\n");
        }
        else
            fclose(existspinfile);

        for (int i=0; i<sitenr; i++)
        {
            for (int l=0; l<sitenr; l++)
            {
                for (int k=0; k<cellnr; k++)
                {
                    for (int j=0; j<compnr; j++)
                    {
                        fprintf(spinfile, "%.6f\t", real(spin_out[i][l][k][j]));
                    }
                }
            }
        }

        fprintf(spinfile, "\n");
        fclose(spinfile);
    }

    void print_n(Parameters params, cd3D nout)
    {
//        FILE * existdensfile = fopen("dens.txt", "r");
//        FILE * densfile = fopen("dens.txt", "a");
//
//        int lsize = params.get_lattice_size();
//        int components = params.get_compnr();
//
//        if(existdensfile == NULL)
//        {
//            fprintf(densfile, "Final density for each site and component\n");
//            for(int xsite = 0; xsite < lsize; xsite++)
//            {
//                for(int ysite = 0; ysite <lsize; ysite++)
//                    fprintf(densfile, "site %i%i\t\t\t\t", xsite, ysite);
//            }
//            fprintf(densfile, "\n");
//            for(int site = 0; site < lsize*lsize; site++)
//            {
//                for(int comp = 0; comp < components; comp++)
//                    fprintf(densfile, "comp %i\t\t", comp);
//            }
//            fprintf(densfile, "\n");
//        }
//        else
//            fclose(existdensfile);
//
//        for(int xsite = 0; xsite < lsize; xsite++)
//        {
//            for(int ysite = 0; ysite <lsize; ysite++)
//            {
//                for(int comp = 0; comp < components; comp++)
//                    fprintf(densfile, "%.6f\t", real(nout[xsite][ysite][comp]));
//
//            }
//        }
//
//        fprintf(densfile, "\n");
//        fclose(densfile);
    }

    void print_energy(Parameters params, cd eigval)
    {
        printf("%f\n", real(eigval));
        FILE * existenergyfile = fopen("energy.txt", "r");
        FILE * energyfile = fopen("energy.txt", "a");

        if(existenergyfile == NULL)
            fprintf(energyfile, "total energy with quadratic mean-field correction\n");
        else
            fclose(existenergyfile);

        fprintf(energyfile, "%.6f\t", real(eigval));

        fprintf(energyfile, "\n");
        fclose(energyfile);
    }
}
