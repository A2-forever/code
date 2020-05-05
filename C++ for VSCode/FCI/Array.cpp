#ifndef ARRAY_CPP_
#define ARRAY_CPP_

#include "Array.h"


double **double2_new(const int row, const int column)//创建二维double数组
{
    double **M = new double *[row];
    for( int i=0; i<row; i++ )
        M[i] = new double [column];

    return M;
}

double **double2_new(const long long row, const long long column)//创建二维double数组,long long
{
    double **M = new double *[row];
    for( long long i=0; i<row; i++ )
        M[i] = new double [column];

    return M;
}

double ***double3_new(const int r1, const int r2, const int r3)//创建三维double数组
{
    double ***M = new double **[r1];
    for( int i=0; i<r1; i++ )
        M[i] = double2_new(r2,r3);

    return M;
}

double ****double4_new(const int r1, const int r2, const int r3, const int r4)//创建四维double数组
{
    double ****M = new double ***[r1];
    for( int i=0; i<r1; i++ )
        M[i] = double3_new(r2,r3,r4);

    return M;
}


void double2_delete(double **M, const int row)//释放二维double数组
{
    for(int i=0;i<row;i++)
        delete [] M[i];
    delete []M;

}

void double2_delete(double **M, const long long row)//释放二维double数组,long long
{
    for(long long i=0;i<row;i++)
        delete [] M[i];
    delete []M;

}

void double3_delete(double ***M, const int r1, const int r2)//释放三维double数组
{
    for(int i=0;i<r2;i++)
        double2_delete(M[i],r2);
    delete []M;
}

void double4_delete(double ****M, const int r1, const int r2, const int r3)//释放四维double数组
{
    for(int i=0;i<r1;i++)
        double3_delete(M[i],r2,r3);
    delete []M;
}



int **int2_new(const int row, const int column)//创建int二维数组
{
    int **M = new int *[row];
    for( int i=0; i<row; i++ )
        M[i] = new int [column];

    return M;
}

int ***int3_new(const int r1, const int r2, const int r3)//创建三维int数组
{
    int ***M = new int **[r1];
    for( int i=0; i<r1; i++ )
        M[i] = int2_new(r2,r3);

    return M;
}

int ****int4_new(const int r1, const int r2, const int r3, const int r4)//创建四维int数组
{
    int ****M = new int ***[r1];
    for( int i=0; i<r1; i++ )
        M[i] = int3_new(r2,r3,r4);

    return M;
}


void int2_delete(int **M, const int row)//释放二维int数组
{
    for(int i=0;i<row;i++)
        delete [] M[i];
    delete []M;

}

void int3_delete(int ***M, const int r1, const int r2)//释放三维int数组
{
    for(int i=0;i<r2;i++)
        int2_delete(M[i],r2);
    delete []M;
}

void int4_delete(int ****M, const int r1, const int r2, const int r3)//释放四维int数组
{
    for(int i=0;i<r1;i++)
        int3_delete(M[i],r2,r3);
    delete []M;
}


void output(double **M,const int n)//输出二维int数组
{
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            std::cout<<M[i][j]<<"\t";
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

#endif