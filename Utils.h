#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <complex.h>
#include <vector>
#include "include/Eigstatval.h"

using namespace std;

typedef complex<double> cd;
typedef vector<cd> cd1D;
typedef vector<cd1D > cd2D;
typedef vector<cd2D> cd3D;
typedef unsigned int unint;

//for diagonalization routine
extern "C" double dlamch_(char *cmach);
extern "C" void   zheevr_(char * JOBZ, char * RANGE, char * UPLO,int* n, cd* A, int* LDA, double* VL, double* VU, int* IL,
                          int* IU, double *ABSTOL, int* M, double *W, cd *Z, int *LDZ, int * ISUPPZ,cd* WORK, int *LWORK,
                          double *RWORK, int *LRWORK, int* IWORK, int *LIWORK, int *INFO );


namespace utils
{
    Eigstatval diagonalize(cd2D diag_array, int lastEW_in);

    cd2D zero_matrix(int dim);
    cd2D unit_matrix(int dim);
    cd2D n_matrix(int dim);
    cd2D a_matrix(int dim);
    cd2D adag_matrix(int dim);
    cd3D vecunit(int dim1, int dim23);
    cd2D cron(cd2D a1, cd2D a2);
    cd2D cron(cd3D a);
    cd1D conj(cd1D a);

    template <typename T>
    vector<T> operator+(const vector<T>& a, const vector<T>& b)
    {
        vector<T> result;

        if (a.size()==b.size())
        {
            for (int i=0; i < a.size(); i++)
            {
                result.push_back(a[i]+b[i]);
            }
        }
        else
        {
            printf("plus 1d 1d: size of vectors doesnt match!");
            exit(1);
        }

        return result;
    }

    template <typename T>
    vector<vector<T>> operator+(const vector<vector<T>>& a, const vector<vector<T>>& b)
    {
        vector<vector<T>> result;

        if (a.size()==b.size() and a[0].size()==b[0].size())
        {
            vector<T> temp;

            for (int i=0; i < a.size(); i++)
            {
                result.push_back(temp);
                for (int j=0; j < a[0].size(); j++)
                {
                    result[i].push_back(a[i][j]+b[i][j]);
                }
            }
        }
        else
        {
            printf("plus 2d 2d: size of cd2D matrices doesnt match!");
            exit(1);
        }

        return result;
    }

    template <typename T>
    vector<vector<vector<T>>> operator+(const vector<vector<vector<T>>>& a, const vector<vector<vector<T>>>& b)
    {
        vector<vector<vector<T>>> result;

        if (a.size()==b.size() and a[0].size()==b[0].size() and a[0][0].size()==b[0][0].size())
        {

            for (int h=0; h<a.size(); h++)
            {
                vector<vector<T>> tttemp;
                for (int i=0; i < a[h].size(); i++)
                {
                    vector<T> ttemp;
                    for (int j=0; j < a[h][i].size(); j++)
                    {
                        ttemp.push_back(a[h][i][j]+b[h][i][j]);
                    }
                    tttemp.push_back(ttemp);
                }
                result.push_back(tttemp);
            }
        }
        else
        {
            printf("plus 3d 3d: size of cd3D matrices doesnt match!");
            exit(1);
        }

        return result;
    }

    template <typename T>
    vector<vector<vector<vector<T>>>> operator+(const vector<vector<vector<vector<T>>>>& a, const vector<vector<vector<vector<T>>>>& b)
    {
        vector<vector<vector<vector<T>>>> result;

        if (a.size()==b.size() and a[0].size()==b[0].size() and a[0][0].size()==b[0][0].size() and a[0][0][0].size()==b[0][0][0].size())
        {

            for (int g=0; g<a.size();g++)
            {
                vector<vector<vector<T>>> ttttemp;
                for (int h=0; h<a[g].size(); h++)
                {
                    vector<vector<T>> tttemp;
                    for (int i=0; i < a[g][h].size(); i++)
                    {
                        vector<T> ttemp;
                        for (int j=0; j < a[g][h][i].size(); j++)
                        {
                            ttemp.push_back(a[g][h][i][j]+b[g][h][i][j]);
                        }
                        tttemp.push_back(ttemp);
                    }
                    ttttemp.push_back(tttemp);
                }
                result.push_back(ttttemp);
            }
        }
        else
        {
            printf("plus 4d 4d: size of cd4D matrices doesnt match!");
            exit(1);
        }

        return result;
    }



    template <typename T>
    vector<vector<T>> operator+=(vector<vector<T>>& a, const vector<vector<T>>& b)
    {
        if (a.size()==b.size() and a[0].size()==b[0].size())
        {
            for (unint i=0; i < a.size(); i++)
            {
                for (unint j=0; j < a[0].size(); j++)
                {
                    a[i][j] = a[i][j] + b[i][j];
                }
            }
        }
        else
        {
            printf("plus equal 2d 2d: size of cd2D matrices doesnt match!");
            exit(1);
        }

        return a;
    }

    template <typename T>
    vector<vector<T>> operator-=(vector<vector<T>>& a, const vector<vector<T>>& b)
    {
        if (a.size()==b.size() and a[0].size()==b[0].size())
        {
            for (int i=0; i < a.size(); i++)
            {
                for (int j=0; j < a[0].size(); j++)
                {
                    a[i][j] = a[i][j] - b[i][j];
                }
            }
        }
        else
        {
            printf("plus equal 2d 2d: size of cd2D matrices doesnt match!");
            exit(1);
        }

        return a;
    }

    template <typename T>
    vector<T> operator-(const vector<T>& a, const vector<T>& b)
    {
        vector<T> result;

        if (a.size()==b.size())
        {
            for (int i=0; i < a.size(); i++)
            {
                result.push_back(a[i]-b[i]);
            }
        }
        else
        {
            printf("minus 1d 1d: size of vectors doesnt match!");
            exit(1);
        }

        return result;
    }

    template <typename T>
    vector<vector<T>> operator-(const vector<vector<T>>& a, const vector<vector<T>>& b)
    {
        vector<vector<T>> result;

        if (a.size()==b.size() and a[0].size()==b[0].size())
        {
            vector<T> temp;

            for (unint i=0; i < a.size(); i++)
            {
                result.push_back(temp);
                for (unint j=0; j < a[0].size(); j++)
                {
                    result[i].push_back(a[i][j]-b[i][j]);
                }
            }
        }
        else
        {
            printf("plus 2d 2d: size of cd2D matrices doesnt match!");
            exit(1);
        }

        return result;
    }

    template <typename T>
    vector<T> operator*(const T& a, const vector<T>& b)
    {
        vector<T> result;

        for (int i=0; i < b.size(); i++)
        {
            result.push_back(a*b[i]);
        }

        return result;
    }

    template <typename T>
    vector<vector<T>> operator*(const T& a, const vector<vector<T>>& b)
    {
        vector<vector<T>> result;
        vector<T> temp;

        for (unint i=0; i < b.size(); i++)
        {
            result.push_back(temp);
            for (unint j=0; j < b[i].size(); j++)
            {
                result[i].push_back(a*b[i][j]);
            }
        }

        return result;
    }

    template <typename T>
    vector<vector<vector<T>>> operator*(const T& a, const vector<vector<vector<T>>>& b)
    {
        vector<vector<vector<T>>> result;

        for (unint k=0; k<b.size();k++)
        {
            vector<vector<T>> tttemp;
            for (unint i=0; i < b[k].size(); i++)
            {
                vector<T> ttemp;
                for (unint j=0; j < b[k][i].size(); j++)
                {
                    ttemp.push_back(a*b[k][i][j]);
                }
                tttemp.push_back(ttemp);
            }
            result.push_back(tttemp);
        }

        return result;
    }

    template <typename T>
    vector<vector<vector<vector<T>>>> operator*(const T& a, const vector<vector<vector<vector<T>>>>& b)
    {
        vector<vector<vector<vector<T>>>> result;

        for (unint h=0; h<b.size();h++)
        {
            vector<vector<vector<T>>> ttttemp;
            for (unint k=0; k<b[h].size();k++)
            {
                vector<vector<T>> tttemp;
                for (unint i=0; i < b[h][k].size(); i++)
                {
                    vector<T> ttemp;
                    for (unint j=0; j < b[h][k][i].size(); j++)
                    {
                        ttemp.push_back(a*b[h][k][i][j]);
                    }
                    tttemp.push_back(ttemp);
                }
                ttttemp.push_back(tttemp);
            }
            result.push_back(ttttemp);
        }

        return result;
    }

    template <typename T>
    vector<T> operator*(const vector<T>& a, const vector<vector<T>>& b)
    {
        vector<T> result;

        if (a.size()==b.size())
        {
            T ttemp;
            for (unint i=0; i < a.size(); i++)
            {
                ttemp = 0;
                for (unint k=0; k < a.size(); k++)
                {
                    ttemp += a[k]*b[k][i];
                }
                result.push_back(ttemp);
            }
        }
        else
        {
            printf("multiply 1d 2d: size of matrices doesnt match! sizes are %lu and %lux%lu",a.size(),b.size(),b[0].size());
            exit(1);
        }

        return result;
    }

    template <typename T>
    T operator*(const vector<T>& a, const vector<T>& b)
    {
        T result = 0;

        if (a.size()==b.size())
        {
            for (unint i=0; i < a.size(); i++)
            {
                result += a[i]*b[i];
            }
        }
        else
        {
            printf("multiply 1d 1d: size of vectors doesnt match!");
            exit(1);
        }

        return result;
    }

    template <typename T>
    vector<T> operator*(const vector<vector<T>>& a, const vector<T>& b)
    {
        vector<T> result;

        if (a.size()==b.size())
        {
            T ttemp;
            for (int i=0; i < a.size(); i++)
            {
                ttemp = 0;
                for (int k=0; k < a.size(); k++)
                {
                    ttemp += a[i][k]*b[k];
                }
                result.push_back(ttemp);
            }
        }
        else
        {
            printf("multiply 1d 2d: size of matrices doesnt match!");
            exit(1);
        }

        return result;
    }

    template <typename T>
    vector<vector<T>> operator*(const vector<vector<T>>& a, const vector<vector<T>>& b)
    {
        vector<vector<T>> result;
        vector<T> temp;

        if (a.size()==b.size() and a[0].size()==b[0].size() and a.size()==a[0].size())
        {
            T ttemp = 0;

            for (unint i=0; i < a.size(); i++)
            {
                result.push_back(temp);
                for (unint j=0; j < a[0].size(); j++)
                {
                    ttemp = 0;
                    for (unint k=0; k < a.size(); k++)
                    {
                        ttemp += a[i][k]*b[k][j];
                    }
                    result[i].push_back(ttemp);
                }
            }
        }
        else
        {
            printf("multiply 2d 2d: size of cd2D matrices doesnt match!");
            exit(1);
        }

        return result;
    }
}

#endif // UTILS_H_INCLUDED
