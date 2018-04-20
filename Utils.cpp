#include "Utils.h"

namespace utils
{
    Eigstatval diagonalize(cd2D diag_array, int lastEW_in)
    {
//        printf("diagonalize matrix:\n");
        int dim = diag_array.size();

//        for (int kk = 0; kk<dim; kk++)
//        {
//            for (int ll=0; ll<dim; ll++) printf("%f %f\t",real(diag_array[kk][ll]),imag(diag_array[kk][ll]));
//            printf("\n");
//        }

        int last_EW = lastEW_in;

        //Input variables for zheevr
        char JOBZ = 'V';                        //wenn EW und EZ berechnet werden sollen
        char RANGE = 'I';                       //wenn nur gewisse Anzahl der EW werden berechnet
        char UPLO = 'U';                        //unter Benutzung der oberen Dreiecksmatrix
        int n = dim;                      //Order of matrix to diagonalize
        double *W;                              //(output)array of the Eigenvalues. The first M elements contain the
        //selected eigenvalues in ascending order.
        cd *Z;                      //(output)array of orth.norm. Eigenvectors(column i = EV to W(i))
        int LDA = dim;                            //Leading dimension of matrix to diagonalize
        cd *WORK;
        double *RWORK;
        int LWORK, LIWORK, LRWORK;
        int *IWORK;
        double VL = 0.0, VU = 0.0;              //spielt nur eine Rolle, wenn RANGE = 'V'
        int IL = 1;                             //Index des kleinsten zu berechnenden EW's
        int IU = last_EW;                       //Index des hoechsten zu berechnenden Ew's
        char safemin = 'S';                     //S for safe minimum. It means high relative accuracy in future release ;)
        double ABSTOL = dlamch_(&safemin);      //The absolute error tolerance for the eigenvalues
        int M ;                                 //The total number of eigenvalues found, if RANGE = 'I', M = IU-IL+1
        M = IU-IL+1;
        int *ISUPPZ;

        LIWORK = -1;                            //after first call of zheevr in WORK array the optimal length of the
        LRWORK = -1;                            //arrays will be stored. The arrays WORK, LWORK, IWORK must be allocated
        LWORK = -1;                             //than and the routine zheevr has to be called again for diagonalisation

        int INFO;                               //(output) = 0:  successful exit.
        //< 0:  if INFO = -i, the i-th argument had an illegal value.
        //> 0:  if INFO = i, the algorithm failed to converge; i
        //off-diagonal elements of an intermediate tridiagonal
        //form did not converge to zero.

        //dummy variables for finding optimal length of arrays
        cd dummy_work[2];
        int dummy_iwork[2];
        double dummy_rwork[2];
        cd *matrix;

        //Allocate all arrays
        ISUPPZ = (int *) calloc(2*last_EW, sizeof(int));
        Z = (cd*) calloc (dim*last_EW, sizeof(cd));
        W = (double*) calloc (dim, sizeof(double));
        matrix = (cd*) calloc(dim*dim, sizeof(cd));

        int index = 0;
        for ( int j = 0; j < dim; j++)
        {
            for (int i = 0; i < dim; i++)
            {
                matrix[index] = diag_array[i][j]; //schreibt die Eintraege aus Hmatrix spaltenweise ab
                //printf("lapack matrix %i %f\n", index, matrix[index]);
                index++;
            }
        }
        // fclose(subfile);
        // first call for finding optimal sizes for arrays
        zheevr_(&JOBZ, &RANGE, &UPLO, &n, matrix, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, W, Z, &n, ISUPPZ, dummy_work, &LWORK, dummy_rwork, &LRWORK, dummy_iwork, &LIWORK, &INFO);

        // Here we allocate what the workspace query returned.
        LWORK = (int) real(dummy_work[0]);
        LRWORK = (int) dummy_rwork[0];
        LIWORK = dummy_iwork[0];

        WORK = (cd*) calloc (LWORK, sizeof(cd));
        RWORK = (double*) calloc (LRWORK, sizeof(double));
        IWORK = (int*) calloc (LIWORK, sizeof(int));

        // diagonalisation
        zheevr_(&JOBZ, &RANGE, &UPLO, &n, matrix, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, W, Z, &n, ISUPPZ, WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK, &INFO);

        cd2D gstate;
        cd1D eigval;

        // schreib den Grundzst (erste N*N-1 Stellen im Z-array) in ein globales Array
        for (int i = 0; i < dim; i++)
        {
            cd1D prelimgstate;
            gstate.push_back(prelimgstate);
            for (int j=0; j<dim; j++)
            {
                gstate[i].push_back(Z[i*dim+j]);
//                printf("diagonalize %f %f\n", real(Z[i*dim+j]), imag(Z[i*dim+j]));
                //printf("diagonalize2 %f\n", gstate[i][j]);
            }
        }

        for(int i = 0; i < dim; i++)
        {
            eigval.push_back(W[i]);
//            printf("eigval %f %f\n",real(W[i]), imag(W[i]));
            //printf("eigval2 %f\n", eigval[i]);
        }

        Eigstatval result(eigval, gstate);

        free(Z);
        free(W);
        free(WORK);
        free(RWORK);
        free(IWORK);
        free(ISUPPZ);
        free(matrix);

        return result;
    }

    cd2D zero_matrix(int dim)
    {
        cd2D matrix;
        for (int i =0; i<dim; i++)
        {
            cd1D temp;
            for (int j=0; j<dim; j++)
            {
                temp.push_back(0);
            }
            matrix.push_back(temp);
        }
        return matrix;
    }

    cd2D unit_matrix(int dim)
    {
        cd2D umatrix=zero_matrix(dim);

        for (int i=0; i<dim; i++)
        {
            umatrix[i][i]=cd(1,0);
        }

        return umatrix;
    }

    cd2D n_matrix(int dim)
    {
        cd2D nmatrix=zero_matrix(dim);

        for (int i=0; i<dim; i++)
        {
            nmatrix[i][i]=cd(i,0);
        }

        return nmatrix;
    }

    cd2D a_matrix(int dim)
    {
        cd2D amatrix=zero_matrix(dim);

        for (int i=1; i<dim; i++)
        {
            amatrix[i][i-1]=cd(sqrt(i),0);
        }

        return amatrix;
    }

    cd2D adag_matrix(int dim)
    {
        cd2D admatrix=zero_matrix(dim);

        for (int i=1; i<dim; i++)
        {
            admatrix[i-1][i]=cd(sqrt(i),0);
        }

        return admatrix;
    }

    cd2D cron(cd2D a1, cd2D a2)
    {
        int a11size = a1.size();
        int a12size = a1[0].size();
        int a21size = a2.size();
        int a22size = a2[0].size();

        cd2D result;
        cd1D row;

        for(int ia11 = 0; ia11<a11size;ia11++)
        {
            for(int ia21 = 0; ia21<a21size; ia21++)
            {
                result.push_back(row);
                for(int ia12 = 0; ia12<a12size; ia12++)
                {
                    for(int ia22=0; ia22<a22size; ia22++)
                    {
                        result[ia11*a21size+ia21].push_back(a1[ia11][ia12]*a2[ia21][ia22]);
                    }
                }
            }
        }

        return result;
    }

    cd2D cron(cd3D a)
    {
        cd2D result = a[0];
        cd2D temp = a[0];

        for(unsigned int i=1; i < a.size(); i++)
        {
            result = cron(temp,a[1]);
            temp = result;
        }

        return result;
    }

    cd1D conj(cd1D a)
    {
        cd1D result;
        for(unsigned int i=0; i<a.size(); i++)
            result.push_back(conj(a[i]));
        return result;
    }

    cd3D vecunit(int dim1, int dim23)
    {
        cd3D result;
        for (int i=0; i<dim1; i++)
        {
            result.push_back(unit_matrix(dim23));
        }
        return result;
    }

    cd2D mult(cd a, cd2D b)
    {
        cd2D result;
        cd1D temp;

        for (unsigned int i=0; i < b.size(); i++)
        {
            result.push_back(temp);
            for (unsigned int j=0; j < b[i].size(); j++)
            {
                result[i].push_back(a*b[i][j]);
            }
        }

        return result;
    }
}
