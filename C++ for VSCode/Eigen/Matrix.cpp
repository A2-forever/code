#include <cmath>
#include "Matrix.h"
#define eps 1e-6

using namespace std;

namespace MATRIX
{
    Matrix::Matrix()//构造函数，默认行列为0
    {
        row = 0;
        column = 0;
    }


    Matrix::Matrix(const Matrix &A)//两个矩阵之间相互赋值
    {
        row = A.row;
        column = A.column;
        for (int i = 0; i < row; i++)
            for (int j = 0; j < column; j++)
                M[i][j] = A.M[i][j];
    }

    Matrix::Matrix(const double  **N, int r,int c)//将一个数组赋值给矩阵，数组N大小为rxc
    {
        row = r;
        column = c;
        for (int i = 0; i < row; i++)
            for (int j = 0; j < column; j++)
                M[i][j] = N[i][j];
    }

    void Matrix::set_position(int r,int c)//更改矩阵的行列
    {
        row = r;
        column = c;
    }

    void Matrix::set_Num(int r,int c, float k)//指定某一位置的值,M[r,c]=k
    {
        M[r - 1][c - 1] = k;
    }


    void Matrix::operator=(const Matrix &A)//矩阵赋值
    {
        row = A.row;
        column = A.column;
        for (int i = 0; i < row; i++)
            for (int j = 0; j < column; j++)
                M[i][j] = A.M[i][j];
    }

    Matrix Matrix::operator+(const Matrix &A) const//矩阵加法
    {
        Matrix *New_Matrix=new Matrix;
        &New_Matrix.set_position(row, column);
        for (int i = 0; i < row; i++)
            for (int j = 0; j < column; j++)
                &New_Matrix.M[i][j] = A.M[i][j] + M[i][j];

        return &New_Matrix;
    }

    Matrix Matrix::operator-(const Matrix &A) const//矩阵减法
    {
        Matrix *New_Matrix=new Matrix;
        &New_Matrix.set_position(row, column);
        for (int i = 0; i < row; i++)
            for (int j = 0; j < column; j++)
                &New_Matrix.M[i][j] = M[i][j] - A.M[i][j];

        return &New_Matrix;
    }

    Matrix Matrix::operator*(const Matrix &A) const//矩阵乘法
    {
        Matrix *New_Matrix=new Matrix;
        &New_Matrix.set_position(row, column);
        &New_Matrix = Matrix_multiply(M, A.M, row, column, A.row);

        return &New_Matrix;
    }

    Matrix Matrix::operator*(const double n) const//矩阵与数字相乘
    {
        Matrix *New_Matrix=new Matrix;
        &New_Matrix.set_position(row, column);
        for (int i = 0; i < row; i++)
            for (int j = 0; j < column; j++)
                &New_Matrix.M[i][j] = M[i][j] * n;

        return &New_Matrix;
    }

    void Matrix::swag(int r1,int r2,char ch='c')//交换矩阵的行列，默认为列
    {
        if(ch=='c')
        {
            for (int i = 0; i < column; i++)
                swag(M[i][r1],M[i][r2]);
        }
        else if(ch=='r')
        {
            for (int i = 0; i < row; i++)
                swag(M[r1][i],M[r2][i]);
        }
    }

    void Matrix::Primary_Matrix(int r,int c, double k,char ch='c')//矩阵的初等变换，默认为列变换，第c列（行）乘以k，加到第r列（行）上
    {
        if(ch=='c')
        {
            for (int i = 0; i < column; i++)
                M[i][r] += M[i][c] * k;
        }
        else if(ch=='r')
        {
            for (int i = 0; i < row; i++)
                M[r][i] += M[c][i] * k;
        }
    }

    Matrix Matrix_multiply(const double **M, const double **N, unsigned int r, unsigned int c, unsigned int c1)//对于两个数组相乘
    {
        Matrix New_Matrix;
        New_Matrix.set_position(r, c1);
        double temp = 0;
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c1; j++)
            {
                temp = 0;
                for (int k = 0; k < c; k++)
                    temp += M[i][k] * N[k][j];
                New_Matrix.set_Number(i, j, temp);
            }
        return New_Matrix;
    }

    void swag(double &a,double &b)//两个变量的值
    {
        double temp=a;
        a=b;
        b=temp;
    }

} // namespace MATRIX